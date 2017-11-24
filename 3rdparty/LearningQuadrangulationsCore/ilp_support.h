#ifndef ILP_SUPPORT
#define ILP_SUPPORT

#include "SCIP/scip_solver.h"
#include "patch.h"
#include "polychord_support.h"

#define TOLERANCE 1

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The ILPSupport class provides some helpful functions for building integer linear problems.
 */
template < typename PolyMeshType >
class ILPSupport {
public:
  typedef pl::Patch<PolyMeshType>                           Patch;
  typedef pl::PolychordSupport<PolyMeshType>                PolychordSupport;
  typedef typename PolychordSupport::FlowPoints             FlowPoints;
  typedef typename PolychordSupport::FlowContainer          FlowContainer;
  typedef typename PolychordSupport::FlowMap                FlowMap;
  typedef typename PolychordSupport::FlowMapConstIterator   FlowMapConstIterator;
  typedef SCIP_Solver<num_type>                             ILPSolver;
  typedef typename ILPSolver::VarType                       ILPVarType;
  typedef std::vector<num_type>                             ILPArrayType;
  typedef std::vector<ILPArrayType>                         ILPMatrixType;

  /**
   * @brief ConstraintILPbyFlows adds some constraints to an existing problem drived by flows.
   * @param patch The input Patch where to check the flows.
   * @param nRotations The number of rotations to make to the patch.
   * @param boundaryLength The boundary constraints.
   * @param flows The flow constraints.
   * @param objectiveFun The initial and final objective function.
   * @param lowerBounds The initial and final variable lower bounds.
   * @param upperBounds The initial and final variable upper bounds.
   * @param constraints The initial and final constraints.
   * @param rhs The initial and final right hand sides of the constraints.
   * @return true if at least a solution has been found, false otherwise.
   */
  static bool ConstraintILPbyFlows(const Patch &patch, const n_corners_type nRotations,
                                   const std::vector<num_type> &boundaryLength, const FlowMap &flows,
                                   ILPArrayType &objectiveFun, ILPArrayType &lowerBounds, ILPArrayType &upperBounds,
                                   ILPMatrixType &constraints, ILPArrayType &rhs, const u_small_int_type tolerance = TOLERANCE) {
    if (patch.isNull() || boundaryLength.empty() || nRotations >= boundaryLength.size() || flows.empty())
      return false;
    if (objectiveFun.empty() || lowerBounds.empty() || constraints.empty() || rhs.empty())
      return false;
    if (objectiveFun.size() != lowerBounds.size() || lowerBounds.size() != patch.getBoundaryVarVec().size())
      return false;
    if (constraints.size() != boundaryLength.size() || rhs.size() != boundaryLength.size())
      return false;
    for (size_t i = 0; i < constraints.size(); ++i)
      if (constraints[i].size() != objectiveFun.size())
        return false;

    n_corners_type nSides = boundaryLength.size();

    // first step: serialize flows and map boundary vars to flow topologies
    FlowTopologyVarsMap flowTopologyVarsMap;
    n_corners_type startSide, endSide;
    var_label_type boundaryID, currentNVariables, newDummyVar1, newDummyVar2, newMiddleDummyVar1, newMiddleDummyVar2;
    VarsArray varsArray;
    SerializedFlows serializedFlows;
    FlowData flowData;
    for (FlowMapConstIterator flowsIt = flows.begin(); flowsIt != flows.end(); ++flowsIt) {
      const FlowTopology &flowTopology = flowsIt->first;
      const FlowContainer &flowContainer = flowsIt->second;
      if (flowTopology.startSide >= nSides || flowTopology.endSide >= nSides)
        continue;
      // collect boundary vars
      varsArray.clear();
      startSide = (flowTopology.startSide + nRotations) % nSides;
      endSide = (flowTopology.endSide + nRotations) % nSides;
      for (size_t edge = 0; edge < patch.getBoundaryIDsVec()[startSide].size(); ++edge) {
        boundaryID = patch.getBoundaryIDsVec()[startSide][edge];
        if ((patch.getBoundaryIDsEnds()[boundaryID].first == startSide && patch.getBoundaryIDsEnds()[boundaryID].second == endSide) ||
            (patch.getBoundaryIDsEnds()[boundaryID].second == startSide && patch.getBoundaryIDsEnds()[boundaryID].first == endSide))
          varsArray.push_back(boundaryID);
      }
      if (varsArray.empty())  // if no polychords in these sides, no constraints may be satisfied
        return false;
      flowTopologyVarsMap.insert(std::make_pair(flowTopology, varsArray));

      // collect flow data
      for (size_t j = 0; j < flowContainer.size(); ++j) {
        const FlowPoints &flowPoints = flowContainer[j];
        // check edge indices
        if (flowPoints.startEdge < 0 || flowPoints.startEdge >= boundaryLength[flowTopology.startSide] ||
            flowPoints.endEdge < 0 || flowPoints.endEdge >= boundaryLength[flowTopology.endSide])
          continue;
        // add to serialized flows
        flowData.flowTopology = flowTopology;
        flowData.startEdge = flowPoints.startEdge;
        flowData.endEdge = flowPoints.endEdge;
        serializedFlows.push_back(flowData);
      }
    }
    // if there aren't any serialized flows (probably only loops, which cannot be handled here), return
    if (serializedFlows.empty())
      return true;

    // second step: enqueue and serialize ILP problems
    ILPSolver ilp_solver;
    ILPArrayType solution;
    num_type objectiveVal;
    ILPData ilp_data;
    ilp_data.objectiveFun = objectiveFun;
    ilp_data.lowerBounds = lowerBounds;
    ilp_data.upperBounds = upperBounds;
    ilp_data.constraints = constraints;
    ilp_data.rhs = rhs;
    size_t flowIndex = 0;
    FlowTopologyVarsMapIterator flowTopologyVarsMapIt;
    ProblemStack problemStack;
    problemStack.push(Problem(0, ilp_data));
    Problem currentProblem, newProblem;
    ILPArrayType newConstraint;
    num_type startEdge, endEdge;
    while (!problemStack.empty()) {
      currentProblem = problemStack.top();
      problemStack.pop();
      flowIndex = currentProblem.serializedFlowIndex;
      flowTopologyVarsMapIt = flowTopologyVarsMap.find(serializedFlows.at(flowIndex).flowTopology);
      assert(flowTopologyVarsMapIt != flowTopologyVarsMap.end());
      const FlowTopology &flowTopo = flowTopologyVarsMapIt->first;
      startSide = (flowTopo.startSide + nRotations) % nSides;
      endSide = (flowTopo.endSide + nRotations) % nSides;
      startEdge = serializedFlows.at(flowIndex).startEdge;
      endEdge = serializedFlows.at(flowIndex).endEdge;
      VarsArray &vars = flowTopologyVarsMapIt->second;
      for (size_t v = 0; v < vars.size(); ++v) {
        ilp_data = currentProblem.ilp_data;
        currentNVariables = ilp_data.objectiveFun.size();

        // add two new dummy variables x_b1 and x_b2
        newDummyVar1 = currentNVariables;
        ++currentNVariables;
        newDummyVar2 = currentNVariables;
        ++currentNVariables;
        // resize lower bounds for x_b1 and x_b2
        ilp_data.lowerBounds.resize(currentNVariables, 0);
        // resize upper bounds for x_b1 and x_b2 (infinity)
        ilp_data.upperBounds.resize(currentNVariables, std::numeric_limits<num_type>::max());
        // add two new dummy variables x_m1 and x_m2
        if (tolerance > 0) {
          newMiddleDummyVar1 = currentNVariables;
          ++currentNVariables;
          newMiddleDummyVar2 = currentNVariables;
          ++currentNVariables;
          // resize lower bounds for x_m1 and x_m2
          ilp_data.lowerBounds.resize(currentNVariables, -(num_type)tolerance);
          // resize upper bounds for x_m1 and x_m2
          ilp_data.upperBounds.resize(currentNVariables, +(num_type)tolerance);
        }
        // resize all the current arrays
        ilp_data.objectiveFun.resize(currentNVariables, 1);
        for (size_t c = 0; c < ilp_data.constraints.size(); ++c)
          ilp_data.constraints[c].resize(currentNVariables, 0);
        // substitution: x_b = x_b1 + 1 + x_b2
        // make new constraint
        newConstraint.assign(currentNVariables, 0);
        // variable change constraint 1
        newConstraint[vars[v]] = 1;
        newConstraint[newDummyVar1] = -1;
        newConstraint[newDummyVar2] = -1;
        // add it
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        ilp_data.constraints.push_back(std::move(newConstraint));
#else
        ilp_data.constraints.push_back(newConstraint);
#endif
        // add rhs of this constraints
        ilp_data.rhs.resize(ilp_data.rhs.size() + 1, 1);
        // make first interval constraint on start side
        newConstraint.assign(currentNVariables, 0);
        // add terms before first new dummy variable
        for (size_t edge = 0; patch.getBoundaryIDsVec()[startSide][edge] != vars[v]; ++edge)
          ++newConstraint[patch.getBoundaryIDsVec()[startSide][edge]];
        // add the first dummy variable x_b1
        newConstraint[newDummyVar1] = 1;
        // add the first middle dummy variable x_m1
        if (tolerance > 0)
          newConstraint[newMiddleDummyVar1] = -1;
        // add it
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        ilp_data.constraints.push_back(std::move(newConstraint));
#else
        ilp_data.constraints.push_back(newConstraint);
#endif
        // add rhs to this constraint
        ilp_data.rhs.resize(ilp_data.rhs.size() + 1, startEdge);
        // make last interval constraint on start side
        newConstraint.assign(currentNVariables, 0);
        // add terms after second new dummy variable (or the first new dummy one, if startSide == endSide)
        for (size_t edge = patch.getBoundaryIDsVec()[startSide].size() - 1; patch.getBoundaryIDsVec()[startSide][edge] != vars[v] ; --edge)
          ++newConstraint[patch.getBoundaryIDsVec()[startSide][edge]];
        // add the second dummy variable
        if (startSide != endSide)
          newConstraint[newDummyVar2] = 1;
        else
          newConstraint[newDummyVar1] = 1;
        // add the middle dummy variable
        if (tolerance > 0) {
          if (startSide != endSide)
            newConstraint[newMiddleDummyVar1] = 1;  // x_m1 for startSide != endSide
          else
            newConstraint[newMiddleDummyVar2] = 1;  // x_m2 for startSide == endSide
        }
        // add it
#if __cplusplus >= 201103L || _MSC_VER >= 1800
        ilp_data.constraints.push_back(std::move(newConstraint));
#else
        ilp_data.constraints.push_back(newConstraint);
#endif
        // add rhs to this constraint
        if (startSide != endSide)
          ilp_data.rhs.resize(ilp_data.rhs.size() + 1, boundaryLength[flowTopo.startSide] - startEdge - 1);
        else
          ilp_data.rhs.resize(ilp_data.rhs.size() + 1, boundaryLength[flowTopo.startSide] - endEdge - 1);

        if (startSide != endSide) {
          // make first interval constraint on end side
          newConstraint.assign(currentNVariables, 0);
          // add terms before second new dummy variable
          for (size_t edge = 0; patch.getBoundaryIDsVec()[endSide][edge] != vars[v]; ++edge)
            ++newConstraint[patch.getBoundaryIDsVec()[endSide][edge]];
          // add the second dummy variable x_b2
          newConstraint[newDummyVar2] = 1;
          // add the second middle dummy variable x_m2
          if (tolerance > 0)
            newConstraint[newMiddleDummyVar2] = -1;
          // add it
#if __cplusplus >= 201103L || _MSC_VER >= 1800
          ilp_data.constraints.push_back(std::move(newConstraint));
#else
          ilp_data.constraints.push_back(newConstraint);
#endif
          // add rhs to this constraint
          ilp_data.rhs.resize(ilp_data.rhs.size() + 1, endEdge);
          // make last interval constraint on end side
          newConstraint.assign(currentNVariables, 0);
          // add terms after second new dummy variable
          for (size_t edge = patch.getBoundaryIDsVec()[endSide].size() - 1; patch.getBoundaryIDsVec()[endSide][edge] != vars[v] ; --edge)
            ++newConstraint[patch.getBoundaryIDsVec()[endSide][edge]];
          // add the second dummy variable x_b1
          newConstraint[newDummyVar1] = 1;
          // add the second middle dummy variable x_m2
          if (tolerance > 0)
            newConstraint[newMiddleDummyVar2] = 1;
          // add it
#if __cplusplus >= 201103L || _MSC_VER >= 1800
          ilp_data.constraints.push_back(std::move(newConstraint));
#else
          ilp_data.constraints.push_back(newConstraint);
#endif
          // add rhs to this constraint
          ilp_data.rhs.resize(ilp_data.rhs.size() + 1, boundaryLength[flowTopo.endSide] - endEdge - 1);
        }

        // build and push on stack new problem
        if (flowIndex + 1 < serializedFlows.size()) {
          newProblem.serializedFlowIndex = flowIndex + 1;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
          newProblem.ilp_data = std::move(ilp_data);
#else
          newProblem.ilp_data = ilp_data;
#endif
          problemStack.push(newProblem);
        } else {
          // generate the problem
          ilp_solver.newProblem(ilp_data.objectiveFun, std::vector<ILPVarType>(1, ILPVarType::INTEGER), ilp_data.lowerBounds, ilp_data.upperBounds);
          // set the system
          ilp_solver.addConstraints(ilp_data.constraints, ilp_data.rhs, ilp_data.rhs);
          // solve!
          ilp_solver.solve(solution, objectiveVal);
          // if there's some solution, return true
          if (!solution.empty()) {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
            objectiveFun = std::move(ilp_data.objectiveFun);
            lowerBounds = std::move(ilp_data.lowerBounds);
            upperBounds = std::move(ilp_data.upperBounds);
            constraints = std::move(ilp_data.constraints);
            rhs = std::move(ilp_data.rhs);
#else
            objectiveFun = ilp_data.objectiveFun;
            lowerBounds = ilp_data.lowerBounds;
            upperBounds = ilp_data.upperBounds;
            constraints = ilp_data.constraints;
            rhs = ilp_data.rhs;
#endif
            return true;
          }
        }
      }
    }

    return false;
  }

private:
  typedef std::deque<var_label_type>        VarsArray;
  typedef std::map<FlowTopology, VarsArray> FlowTopologyVarsMap;
  typedef FlowTopologyVarsMap::iterator     FlowTopologyVarsMapIterator;

  /**
   * @brief The FlowData struct stores the only important flow data to serialize.
   */
  struct FlowData {
    FlowTopology      flowTopology;
    num_type          startEdge;
    num_type          endEdge;
  };

  typedef std::deque<FlowData>              SerializedFlows;

  /**
   * @brief The ILPData struct stores the problem of a ILP run.
   */
  struct ILPData {
    ILPArrayType          objectiveFun;
    ILPArrayType          lowerBounds;
    ILPArrayType          upperBounds;
    ILPMatrixType         constraints;
    ILPArrayType          rhs;

    /**
     * @brief ILPData Default constructor.
     */
    ILPData() { }

    /**
     * @brief ILPData Copy constructor.
     * @param d The ILPData to copy.
     */
    ILPData(const ILPData &d) : objectiveFun(d.objectiveFun),
                                lowerBounds(d.lowerBounds),
                                upperBounds(d.upperBounds),
                                constraints(d.constraints),
                                rhs(d.rhs) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
    /**
     * @brief ILPData Move constructor.
     * @param d The ILPData to move.
     */
    ILPData(ILPData&& d) : objectiveFun(std::move(d.objectiveFun)),
                           lowerBounds(std::move(d.lowerBounds)),
                           upperBounds(std::move(d.upperBounds)),
                           constraints(std::move(d.constraints)),
                           rhs(std::move(d.rhs)) { }
#endif

    /**
     * @brief operator = Copy assignment.
     * @param d The ILPData to copy.
     * @return A reference to this.
     */
    ILPData & operator = (const ILPData &d) {
      if (this != &d) {
        objectiveFun = d.objectiveFun;
        lowerBounds = d.lowerBounds;
        upperBounds = d.upperBounds;
        constraints = d.constraints;
        rhs = d.rhs;
      }
      return *this;
    }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
    /**
     * @brief operator = Move assignment.
     * @param d The ILPData to move.
     * @return A reference to this.
     */
    ILPData & operator = (ILPData&& d) {
      if (this != &d) {
        objectiveFun = std::move(d.objectiveFun);
        lowerBounds = std::move(d.lowerBounds);
        upperBounds = std::move(d.upperBounds);
        constraints = std::move(d.constraints);
        rhs = std::move(d.rhs);
      }
      return *this;
    }
#endif
  };

  /**
   * @brief The Problem struct stores the current ILP problem and the next serialized flow to add as constraint.
   */
  struct Problem {
    size_t                    serializedFlowIndex;
    ILPData                   ilp_data;

    /**
     * @brief Problem Default constructor.
     */
    Problem() : serializedFlowIndex(0) { }

    /**
     * @brief Problem COnstructor initializer.
     * @param i The next index.
     * @param d The current ILPData.
     */
    Problem(const size_t i, const ILPData &d) : serializedFlowIndex(i),
                                                ilp_data(d) { }

    /**
     * @brief Problem Copy constructor.
     * @param p The Problem to copy.
     */
    Problem(const Problem &p) : serializedFlowIndex(p.serializedFlowIndex),
                                ilp_data(p.ilp_data) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
    /**
     * @brief Problem Move constructor.
     * @param p The Problem to move.
     */
    Problem(Problem&& p) : serializedFlowIndex(p.serializedFlowIndex),
                           ilp_data(std::move(p.ilp_data)) { }
#endif

    /**
     * @brief operator = Copy assignment.
     * @param p The Problem to copy.
     * @return A reference to this.
     */
    Problem & operator = (const Problem &p) {
      if (this != &p) {
        serializedFlowIndex = p.serializedFlowIndex;
        ilp_data = p.ilp_data;
      }
      return *this;
    }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
    /**
     * @brief operator = Move assignemnt.
     * @param p The Problem to move.
     * @return A reference to this.
     */
    Problem & operator = (Problem&& p) {
      if (this != &p) {
        serializedFlowIndex = p.serializedFlowIndex;
        ilp_data = std::move(p.ilp_data);
      }
      return *this;
    }
#endif
  };

  typedef std::stack<Problem>               ProblemStack;
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // ILP_SUPPORT

