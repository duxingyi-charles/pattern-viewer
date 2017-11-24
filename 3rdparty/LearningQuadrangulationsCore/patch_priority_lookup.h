#ifndef PATCH_PRIORITY_LOOKUP_H
#define PATCH_PRIORITY_LOOKUP_H

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include <memory>
#endif
#include "patch_db_server.h"
#include "ilp_support.h"
#include "predicate.h"
#include "geometry_support.h"
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  #include "string_support.h"
#endif

namespace vcg {
namespace tri {

// namespace pl: patch learning
namespace pl {

/**
 * @brief The PatchPriorityLookup class provides methods for looking up patches with given boundary and flow constraints.
 * Make a query with PatchPriorityLookup::findPatchesFromBoundaryLength,
 * then retrieve the extracted patches with PatchPriorityLookup::at.
 */
template < typename PolyMeshType >
class PatchPriorityLookup {
public:
  typedef pl::PatchDBServer<PolyMeshType>                 PatchDBServer;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef std::shared_ptr<PatchDBServer>                  PatchDBServerPointer;
#else
  typedef PatchDBServer*                                  PatchDBServerPointer;
#endif
  typedef typename PatchDBServer::Patch                   Patch;
  typedef typename Patch::PosType                         PosType;
  typedef typename Patch::FacePointer                     FacePointer;
  typedef typename Patch::VertexPointer                   VertexPointer;
  typedef std::deque<Patch>                               PatchContainer;
  typedef pl::PatchCompare<PolyMeshType>                  PatchCompare;
  typedef typename PatchCompare::CompareCriteria          CompareCriteria;
  typedef pl::PolychordSupport<PolyMeshType>              PolychordSupport;
  typedef typename PolychordSupport::FlowPoints           FlowPoints;
  typedef typename PolychordSupport::FlowContainer        FlowContainer;
  typedef typename PolychordSupport::FlowMap              FlowMap;
  typedef typename PolychordSupport::FlowMapConstIterator FlowMapConstIterator;
  typedef pl::PatchGeometrizator<PolyMeshType>            PatchGeometrizator;
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef std::shared_ptr<PatchGeometrizator>             PatchGeometrizatorPointer;
#else
  typedef PatchGeometrizator *                            PatchGeometrizatorPointer;
#endif

  /**
   * @brief PatchPriorityLookup Default constructor.
   */
  PatchPriorityLookup() : _patchDBServer(0),
                          _minValence(2),
                          _maxValence(std::numeric_limits<valence_type>::max()),
                          _maxSingularities(std::numeric_limits<num_type>::max()),
                          _maxConcavities(std::numeric_limits<n_corners_type>::max()),
                          _minSizeThresholdILP(0),
                          _flowConstraintsILPtolerance(TOLERANCE) { }

  /**
   * @brief PatchPriorityLookup Constructor with a DB server parameter.
   * @param patchDBServer A pointer to the DB server.
   */
  PatchPriorityLookup(const PatchDBServerPointer &patchDBServer) : _patchDBServer(patchDBServer),
                                                                   _minValence(2),
                                                                   _maxValence(std::numeric_limits<valence_type>::max()),
                                                                   _maxSingularities(std::numeric_limits<num_type>::max()),
                                                                   _maxConcavities(std::numeric_limits<n_corners_type>::max()),
                                                                   _minSizeThresholdILP(0),
                                                                   _flowConstraintsILPtolerance(TOLERANCE) { }

  /**
   * @brief PatchPriorityLookup Copy constructor.
   * @param plu The PatchPriorityLookup to copy.
   */
  PatchPriorityLookup(const PatchPriorityLookup &plu) : _patchDBServer(plu._patchDBServer),
                                                        _patches(plu._patches),
                                                        _compare(plu._compare),
                                                        _minValence(plu._minValence),
                                                        _maxValence(plu._maxValence),
                                                        _maxSingularities(plu._maxSingularities),
                                                        _maxConcavities(plu._maxConcavities),
                                                        _minSizeThresholdILP(plu._minSizeThresholdILP),
                                                        _flowConstraintsILPtolerance(plu._flowConstraintsILPtolerance) { }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief PatchPriorityLookup Move constructor.
   * @param plu The PatchPriorityLookup to move.
   */
  PatchPriorityLookup(PatchPriorityLookup&& plu) : _patchDBServer(plu._patchDBServer),
                                                   _patches(std::move(plu._patches)),
                                                   _compare(std::move(plu._compare)),
                                                   _minValence(plu._minValence),
                                                   _maxValence(plu._maxValence),
                                                   _maxSingularities(plu._maxSingularities),
                                                   _maxConcavities(plu._maxConcavities),
                                                   _minSizeThresholdILP(plu._minSizeThresholdILP),
                                                   _flowConstraintsILPtolerance(plu._flowConstraintsILPtolerance) {
    plu._minValence = valence_type(2);
    plu._maxValence = std::numeric_limits<valence_type>::max();
    plu._maxSingularities = std::numeric_limits<num_type>::max();
    plu._maxConcavities = std::numeric_limits<n_corners_type>::max();
    plu._minSizeThresholdILP = std::numeric_limits<num_type>::max();
    plu._flowConstraintsILPtolerance = TOLERANCE;
  }
#endif

  /**
   * @brief ~PatchPriorityLookup Default destructor.
   */
  virtual ~PatchPriorityLookup() { }

  /**
   * @brief operator = Copy assignment operator.
   * @param plu The PatchPriorityLookup to assign.
   * @return A reference to this.
   */
  PatchPriorityLookup & operator = (const PatchPriorityLookup &plu) {
    if (&plu != this) {
      _patchDBServer = plu._patchDBServer;
      _patches = plu._patches;
      _compare = plu._compare;
      _uniquePatches = plu._uniquePatches;
      _minValence = plu._minValence;
      _maxValence = plu._maxValence;
      _maxSingularities = plu._maxSingularities;
      _maxConcavities = plu._maxConcavities;
      _minSizeThresholdILP = plu._minSizeThresholdILP;
      _flowConstraintsILPtolerance = plu._flowConstraintsILPtolerance;
    }
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param plu The PatchPriorityLookup to move.
   * @return A reference to this.
   */
  PatchPriorityLookup & operator = (PatchPriorityLookup&& plu) {
    if (&plu != this) {
      _patchDBServer = plu._patchDBServer;
      _patches = std::move(plu._patches);
      _compare = std::move(plu._compare);
      _uniquePatches = std::move(plu._uniquePatches);
      _minValence = plu._minValence;
      _maxValence = plu._maxValence;
      _maxSingularities = plu._maxSingularities;
      _maxConcavities = plu._maxConcavities;
      _minSizeThresholdILP = plu._minSizeThresholdILP;
      _flowConstraintsILPtolerance = plu._flowConstraintsILPtolerance;
      plu._minValence = valence_type(2);
      plu._maxValence = std::numeric_limits<valence_type>::max();
      plu._maxSingularities = std::numeric_limits<num_type>::max();
      plu._maxConcavities = std::numeric_limits<n_corners_type>::max();
      plu._minSizeThresholdILP = std::numeric_limits<num_type>::max();
      plu._flowConstraintsILPtolerance = TOLERANCE;
    }
    return *this;
  }
#endif

  /**
   * @brief setPatchDBServer changes the patch DB server.
   * @param patchDBServer
   */
  virtual void setPatchDBServer(const PatchDBServerPointer &patchDBServer) {
    _patches.clear();
    _patchDBServer = patchDBServer;
  }

  /**
   * @brief getSortCriteria gets the current criteria for comparing and sorting patches.
   * @return The current sequence of criteria.
   */
  virtual const std::list<CompareCriteria> & getSortCriteria() const {
    return _compare.getCriteria();
  }

  /**
   * @brief setSortCriteria sets the criteria for comparing and sorting patches.
   * @param criteria The sequence of criteria to copy.
   */
  virtual void setSortCriteria(const std::list<CompareCriteria> &criteria) {
    _compare.setCriteria(criteria);
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief setSortCriteria sets the criteria for comparing and sorting patches.
   * @param criteria The sequence of criteria to move.
   */
  virtual void setSortCriteria(std::list<CompareCriteria>&& criteria) {
    _compare.setCriteria(std::move(criteria));
  }
#endif

  /**
   * @brief setMinValence See _minValence.
   * @param minValence
   */
  virtual void setMinValence(const valence_type minValence = 2) {
    _minValence = minValence;
    if (_minValence < valence_type(2))
      _minValence = valence_type(2);
  }

  /**
   * @brief getMinValence See _minValence.
   * @return
   */
  virtual valence_type getMinValence() const {
    return _minValence;
  }

  /**
   * @brief setMaxValence See _maxValence.
   * @param maxValence
   */
  virtual void setMaxValence(const valence_type maxValence = std::numeric_limits<valence_type>::max()) {
    _maxValence = maxValence;
    if (_maxValence < valence_type(2))
      _maxValence = valence_type(2);
  }

  /**
   * @brief getMaxValence See _maxValence.
   * @return
   */
  virtual valence_type getMaxValence() const {
    return _maxValence;
  }

  /**
   * @brief setMaxSigularities See _maxSingularities.
   * @param maxSingularities
   */
  virtual void setMaxSingularities(const num_type maxSingularities = std::numeric_limits<num_type>::max()) {
    _maxSingularities = maxSingularities;
    if (_maxSingularities < 0)
      _maxSingularities = 0;
  }

  /**
   * @brief getMaxSingularities See _maxSingularities.
   * @return
   */
  virtual num_type getMaxSingularities() const {
    return _maxSingularities;
  }

  /**
   * @brief setMaxConcavities See _maxConcavities.
   * @param maxConcavities
   */
  virtual void setMaxConcavities(const n_corners_type maxConcavities = std::numeric_limits<n_corners_type>::max()) {
    _maxConcavities = maxConcavities;
  }

  /**
   * @brief getMaxConcavities See _maxConcavities.
   * @return
   */
  virtual n_corners_type getMaxConcavities() const {
    return _maxConcavities;
  }

  /**
   * @brief setMinSizeThresholdILP See _minSizeThresholdILP.
   * @param th
   */
  virtual void setMinSizeThresholdILP(const num_type th) {
    _minSizeThresholdILP = th;
    if (_minSizeThresholdILP < 0)
      _minSizeThresholdILP = 0;
  }

  /**
   * @brief getMinSizeThresholdILP See _minSizeThresholdILP.
   * @return
   */
  virtual num_type getMinSizeThresholdILP() const {
    return _minSizeThresholdILP;
  }

  /**
   * @brief setFlowConstraintsILPtolerance See _flowConstraintsILPtolerance.
   * @param tolerance
   */
  virtual void setFlowConstraintsILPtolerance(const u_small_int_type tolerance) {
    _flowConstraintsILPtolerance = tolerance;
  }

  /**
   * @brief getFlowConstraintsILPtolerance See _flowConstraintsILPtolerance.
   * @return
   */
  virtual u_small_int_type getFlowConstraintsILPtolerance() const {
    return _flowConstraintsILPtolerance;
  }

  /**
   * @brief numberOfPatches gives the number of patches satisfying the given boundary constraints.
   * @return The number of patches.
   */
  virtual count_type numberOfPatches() const {
    return _patches.size();
  }

  /**
   * @brief findPatchesFromBoundaryLength makes a new query to find Patches with nCorners and boundaryLength and a flows.
   * @param nCorners The number of corners (it may be 0).
   * @param boundaryLength The boundary's sides' length.
   * @param flows A list of flows to check (patches without such flows, if any, are discarded).
   * @param geometrizator An object responsible for assigning geometry.
   * @param predicate A predicate whose evaluation allows for the continuation of the extraction of patch.
   */
  virtual void findPatchesFromBoundaryLength(const n_corners_type nCorners, const std::vector<num_type> &boundaryLength,
                                             const FlowMap &flows = FlowMap(),
                                             const PatchGeometrizatorPointer &geometrizator = PatchGeometrizatorPointer(0),
                                             const PredicatePointer &predicate = PredicatePointer(0)) {
    // clear the container
    _patches.clear();

    if (!_patchDBServer)
      return;

    // check valence range
    if (_maxValence < _minValence)
      std::swap(_minValence, _maxValence);

    num_type perimeter = std::accumulate(boundaryLength.begin(), boundaryLength.end(), 0);

    // first filter: perimeter and singularities and concavities
    std::list<PatchIndexRotationType> patchesList, tmpPatchList;
    structureFilter(patchesList, nCorners, boundaryLength, perimeter, predicate, 0, _patchDBServer->numberOfCachedPatches(nCorners));

    // second filter: flows
    if (!flows.empty()) {
      flowFilter(tmpPatchList, nCorners, patchesList.begin(), patchesList.end(), flows, predicate);
      patchesList.clear();
      patchesList.splice(patchesList.end(), tmpPatchList);
    }

    // last filter: ILP and uniqueness, plus sorting
    ilpFilter(nCorners, boundaryLength, perimeter, patchesList.begin(), patchesList.end(), flows, geometrizator, predicate);
  }

  /**
   * @brief sortPatches sorts the previously found patches.
   * @param flows A list of flows to check (patches without such flows, if any, are discarded).
   * @param geometrizator An object responsible for assigning geometry.
   * @param predicate A predicate whose evaluation allows for the continuation of the extraction of patch.
   */
  virtual void sortPatches(const FlowMap &flows,
                           const PatchGeometrizatorPointer &geometrizator = PatchGeometrizatorPointer(0),
                           const PredicatePointer &predicate = PredicatePointer(0)) {
    // copy buffer
#if __cplusplus >= 201103L || _MSC_VER >= 1800
    PatchContainer tmpPatches(std::move(_patches));
#else
    PatchContainer tmpPatches(_patches);
#endif
    _patches.clear();
    for (size_t i = 0; i < tmpPatches.size(); ++i) {
      if (predicate && !predicate->operator()()) {
        _patches.clear();
        return;
      }
      // compute new geometry
      if (geometrizator)
        geometrizator->geometrizePatch(tmpPatches[i], flows);
      // add in new position
#if __cplusplus >= 201103L || _MSC_VER >= 1800
      addNewPatch(std::move(tmpPatches[i]));
#else
      addNewPatch(tmpPatches[i]);
#endif
    }
  }

  /**
   * @brief at gives access to the patch at index.
   * @param index The index of the patch.
   * @return A reference to the the requested Patch. A reference to a NULL Patch is returned
   * in case of a out-of-bounds index.
   */
  virtual Patch & at(const count_type index) {
    // a static object to return in case of error
    static Patch staticPatch;
    staticPatch.reset();

    // check for indices
    if (index >= _patches.size())
      return staticPatch;

    // return the requested Patch
    return _patches.at(index);
  }

  /**
   * @brief at gives access to the patch at index.
   * @note This method can be used during execution. It is useful for dynamic updating on-screen results.
   * @param index The index of the patch.
   * @param patch The output Patch. In case of out-of-bounds, it is set to NULL.
   */
  virtual void at(const count_type index, Patch &patch) const {
    patch.reset();

    if (index < _patches.size())
      patch = _patches.at(index);
  }

  /**
   * @brief at gives access to the range of patches at [index, index+offset-1].
   * @param index The starting index of the range.
   * @param offset The number of patches in the range.
   * @param patches The output container of patches
   */
  virtual void at(const count_type index, const count_type offset, PatchContainer &patches) const {
    patches.clear();

    for (count_type i = 0; i < offset && index + i < _patches.size(); ++i)
      patches.push_back(_patches.at(index + i));
  }

  /**
   * @brief at gives access to the range of patches at [index, index+offset-1].
   * @param index The starting index of the range.
   * @param offset The number of patches in the range.
   * @param patches The output container of patches
   * @param cmp A PatchCompare object useful to sort the offset patches by any other criteria.
   */
  virtual void at(const count_type index, const count_type offset, PatchContainer &patches, const PatchCompare &cmp) const {
    at(index, offset, patches);
    std::sort(patches.begin(), patches.end(), cmp);
  }

  /**
   * @brief clear empties the buffer of patches.
   */
  virtual void clear() {
    _patches.clear();
  }



protected:
#if __cplusplus >= 201103L || _MSC_VER >= 1800
  typedef Hasher::HashValue                           HashValue;
  typedef std::unordered_set<HashValue>               TopologySetType;
#else
  typedef std::set<std::string>                       TopologySetType;
#endif
  typedef std::pair<count_type, n_corners_type>       PatchIndexRotationType;
  typedef pl::ILPSupport<PolyMeshType>                ILPSupport;
  typedef typename ILPSupport::ILPSolver              ILPSolver;
  typedef typename ILPSupport::ILPVarType             ILPVarType;
  typedef typename ILPSupport::ILPArrayType           ILPArrayType;
  typedef typename ILPSupport::ILPMatrixType          ILPMatrixType;

  /**
   * @brief structureFilter lists all the patches (and their rotations) that pass the perimeter and singularities and concavities test.
   * @param outPatches The list of passed patches.
   * @param nCorners The number of corners.
   * @param boundaryLenVec The boundary sides length to test against.
   * @param perimeter The perimeter of the boundary, which must be equal to the sum of boundaryLenVec[i], for each i.
   * @param predicate A predicate whose evaluation allows for the continuation of the process.
   * @param startIndex The patches starting index.
   * @param endIndex The patches ending index (not included).
   */
  void structureFilter(std::list<PatchIndexRotationType> &outPatches,
                       const n_corners_type nCorners,
                       const std::vector<num_type> &boundaryLenVec,
                       const num_type perimeter,
                       const PredicatePointer &predicate,
                       const count_type startIndex = 0,
                       const count_type endIndex = std::numeric_limits<count_type>::max()) {
    outPatches.clear();
    if (!_patchDBServer)
      return;
    if (perimeter < 4)
      return;
    if (perimeter % 2 != 0)
      return;
    std::vector< std::pair<valence_type,num_type> > singularityVec;
    valence_type minValence, maxValence;
    // for each patch
    for (count_type index = startIndex; index < endIndex && index < _patchDBServer->numberOfCachedPatches(nCorners); ++index) {
      _patchDBServer->at(nCorners, index).getSingularityVec(singularityVec);
      if (singularityVec.empty()) {
        minValence = _minValence;
        maxValence = _maxValence;
      } else {
        minValence = singularityVec.front().first;
        maxValence = singularityVec.back().first;
      }
      if (_patchDBServer->at(nCorners, index).getNumberOfSingularities() <= _maxSingularities &&
          minValence >= _minValence && maxValence <= _maxValence &&
          _patchDBServer->at(nCorners, index).getNumberOfConcaveCorners() <= _maxConcavities &&
          _patchDBServer->at(nCorners, index).getPerimeter() <= perimeter) {
        // for each rotation, only those patches that pass the filter are added
        for (n_corners_type nRotations = 0; nRotations < (nCorners > 0 ? nCorners : 1); ++nRotations) {
          if (predicate && !predicate->operator()()) {
            outPatches.clear();
            return;
          }
          if (compareSidesLength(boundaryLenVec, _patchDBServer->at(nCorners, index).getBoundaryLenVec(), nRotations))
            outPatches.push_back(PatchIndexRotationType(index, nRotations));
        }
      }
    }
  }

  /**
   * @brief flowFilter lists all patches (and their rotations) that pass the flow test.
   * The flow test is passed if there are the polychords from/to sthe sides specified in
   * flows, or if there is any loop when flows specifies it.
   * @param outPatches The list of passed patches.
   * @param nCorners The number of corners.
   * @param firstIt The first patch to check.
   * @param lastIt The last patch to check (not included).
   * @param flows The list of flows to check against.
   * @param predicate A predicate whose evaluation allows for the continuation of the process.
   */
  void flowFilter(std::list<PatchIndexRotationType> &outPatches,
                  const n_corners_type nCorners,
                  std::list<PatchIndexRotationType>::const_iterator firstIt,
                  std::list<PatchIndexRotationType>::const_iterator lastIt,
                  const FlowMap &flows,
                  const PredicatePointer &predicate = PredicatePointer(0)) const {
    outPatches.clear();

    if (!_patchDBServer)
      return;

    FlowMapConstIterator flowsIt;
    n_corners_type startSide, endSide, nSides = (nCorners > 0 ? nCorners : 1);
    var_label_type boundaryID = 0;
    PolyMeshType mesh;
    PosType startPos;
    n_corners_type nc = 0;
    std::deque<PosType> corners;
    std::deque<PosType> polychords;
    std::vector<bool> isFlowPassed(flows.size(), false);
    size_t fInd = 0;
    bool pass = true;

    // for each patch
    for(; firstIt != lastIt; ++firstIt) {
      if (predicate && !predicate->operator()()) {
        outPatches.clear();
        return;
      }
      // get polychords
      _patchDBServer->at(nCorners, firstIt->first).getMeshAndStartCorner(mesh, startPos);
      startPos.FlipE();
      while (!startPos.IsBorder()) {
        startPos.FlipF();
        startPos.FlipE();
      }
      startPos.FlipV();
      nc = PolychordSupport::FindPolygonCorners(mesh, corners, false, false, startPos);
      assert(nc == nCorners);
      // for each flow
      fInd = 0;
      for (flowsIt = flows.begin(); flowsIt != flows.end(); ++flowsIt) {
        if (predicate && !predicate->operator()()) {
          outPatches.clear();
          return;
        }
        isFlowPassed[fInd] = false;
        // border-to-border case
        if (flowsIt->first.startSide < nSides && flowsIt->first.endSide < nSides) {
          startSide = (flowsIt->first.startSide + firstIt->second) % nSides;
          endSide = (flowsIt->first.endSide + firstIt->second) % nSides;
          for (size_t edge = 0; edge < _patchDBServer->at(nCorners, firstIt->first).getBoundaryIDsVec()[startSide].size(); ++edge) {
            boundaryID = _patchDBServer->at(nCorners, firstIt->first).getBoundaryIDsVec()[startSide][edge];
            if ((_patchDBServer->at(nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].first == startSide
                 && _patchDBServer->at(nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].second == endSide) ||
                (_patchDBServer->at(nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].second == startSide
                 && _patchDBServer->at(nCorners, firstIt->first).getBoundaryIDsEnds()[boundaryID].first == endSide)) {
              isFlowPassed[fInd] = true;
              break;
            }
          }
        } else {    // loop case
          PolychordSupport::FindPolychordsWithTopology(mesh, corners, flowsIt->first, polychords);
          if (polychords.size() >= flowsIt->second.size())
            isFlowPassed[fInd] = true;
        }
        ++fInd;
      }
      pass = true;
      for (fInd = 0; fInd < isFlowPassed.size(); ++fInd)
        pass = pass && isFlowPassed[fInd];
      // if pass, push out
      if (pass)
        outPatches.push_back(*firstIt);
    }
  }

  /**
   * @brief ilpFilter adds all the patches that pass the ILP filter and that are unique.
   * @param nCorners The number of corners.
   * @param boundaryLength The boundary's sides' length condition.
   * @param perimeter It is equal to the sum of the boundary's sides' length.
   * @param firstIt The first patch to check.
   * @param lastIt The last patch to check (not included).
   * @param flowFilters A list of flows to check (patches without such flows, if any, are discarded).
   * @param geometrizator The responsible of assigning geometry to patches.
   * @param predicate A predicate whose evaluation allows for the continuation of the process.
   */
  void ilpFilter(const n_corners_type nCorners,
                 const std::vector<num_type> &boundaryLength,
                 const num_type perimeter,
                 std::list<PatchIndexRotationType>::const_iterator firstIt,
                 std::list<PatchIndexRotationType>::const_iterator lastIt,
                 const FlowMap &flows,
                 const PatchGeometrizatorPointer &geometrizator,
                 const PredicatePointer &predicate = PredicatePointer(0)) {
    if (!_patchDBServer)
      return;

    ILPMatrixType solutions;
    Patch patch;
    std::vector<num_type> boundaryVarVec;
    std::deque<PosType> corners;
    PosType startPos;
    n_corners_type nc = 0;
    std::deque<PosType> polychords;
    FlowMapConstIterator flowsIt;
    bool pass = true;
    for (; firstIt != lastIt; ++firstIt) {
      if (predicate && !predicate->operator()()) {
        _patches.clear();
        return;
      }
      // solve
      solveILP(nCorners, boundaryLength, *firstIt, flows, solutions);
      if (!solutions.empty()) {
        // init patch
        patch = _patchDBServer->at(nCorners, firstIt->first);
        patch.setNumberOfRotations(firstIt->second);   ///dxy say: prone to be wrong for new-style pattern
          ///dxy add: save old LenVec
          auto oldLenVec = patch.getBoundaryLenVec();
          ///
        patch.setBoundaryLenVec(boundaryLength);   ///dxy say: may be wrong
        patch.setPerimeter(perimeter);
        boundaryVarVec = patch.getBoundaryVarVec();
        // for each solution
        for (size_t i = 0; i < solutions.size(); ++i) {
          if (predicate && !predicate->operator()()) {
            _patches.clear();
            return;
          }
          // set current solution : set boundaryVarVec
          for (size_t v = 0; v < boundaryVarVec.size(); ++v)
            boundaryVarVec[v] = (num_type)solutions[i][v];
          patch.setBoundaryVarVec(boundaryVarVec);
            ///dxy add: recover the old LenVec temporarily
            patch.setBoundaryLenVec(oldLenVec);
            ///
          patch.buildExpandedTopologyStr();   //dxy say: the LenVec has been expanded, this is wrong!
            ///dxy add: recover the new LenVec
            patch.setBoundaryLenVec(boundaryLength);
            ///
          // if new, add it
          pass = false;
          if (isUnique(patch.getExpandedTopologyStr())) {
            pass = true;
            startPos = patch.getStartCorner();
            startPos.FlipE();
            while (!startPos.IsBorder()) {
              startPos.FlipF();
              startPos.FlipE();
            }
            startPos.FlipV();
            nc = PolychordSupport::FindPolygonCorners(patch.getMesh(), corners, false, false, startPos);
            assert(nc == nCorners);
            for (flowsIt = flows.begin(); flowsIt != flows.end(); ++flowsIt) {
              PolychordSupport::FindPolychordsWithTopology(patch.getMesh(), corners, flowsIt->first, polychords);
              if (polychords.size() < flowsIt->second.size()) {
                pass = false;
                break;
              }
            }
            if (pass && geometrizator)
              pass = geometrizator->geometrizePatch(patch, flows);
            if (pass)
              addNewPatch(patch);
          }
#ifdef SPLIT_POLYCHORD_LOOPS
          if (pass) {
            flowsIt = flows.find(FlowTopology());
            if (flowsIt != flows.end() && flowsIt->second.size() == 1) {
              Patch splitLoopPatch = patch;
              startPos = splitLoopPatch.getStartCorner();
              startPos.FlipE();
              while (!startPos.IsBorder()) {
                startPos.FlipF();
                startPos.FlipE();
              }
              startPos.FlipV();
              nc = PolychordSupport::FindPolygonCorners(splitLoopPatch.getMesh(), corners, false, false, startPos);
              assert(nc == nCorners);
              PolychordSupport::FindPolychordsWithTopology(splitLoopPatch.getMesh(), corners, FlowTopology(), polychords);
              if (polychords.size() == 1) {
                std::vector<FacePointer *> facesToUpdate;
                std::vector<VertexPointer *> verticesToUpdate;
                facesToUpdate.push_back(&splitLoopPatch.getStartCorner().F());
                verticesToUpdate.push_back(&splitLoopPatch.getStartCorner().V());
                for (size_t i = 0; i < polychords.size(); ++i) {
                  facesToUpdate.push_back(&polychords[i].F());
                  verticesToUpdate.push_back(&polychords[i].V());
                }
                for (size_t i = 0; i < polychords.size(); ++i)
                  vcg::tri::PolychordCollapse<PolyMeshType>::SplitPolychord(splitLoopPatch.getMesh(), polychords[i],
                                                                            (splitLoopPatch.getMesh().FN() > 40 ? 3 : 2),
                                                                            facesToUpdate, verticesToUpdate);
                std::string expandedTopology;
                Patch::BuildTopologyStr(splitLoopPatch.getMesh(), splitLoopPatch.getStartCorner(), expandedTopology);
                splitLoopPatch.setExpandedTopologyStr(expandedTopology);
                if (isUnique(expandedTopology)) {
                  pass = geometrizator && geometrizator->geometrizePatch(splitLoopPatch, flows);
                  if (pass)
                    addNewPatch(splitLoopPatch);
                }
              }
            }
          }
#endif
        }
      }
    }
  }

  /**
   * @brief solveILP solves the ILP associated to a given ConvexPatch and a given boundary's sides' length.
   * @param nCorners The number of corners.
   * @param boundaryLength The boundary's sides' length constraints.
   * @param patchIndexRotation The index of the patch and the number of rotations.
   * @param flows The FlowMap to strongly constraint the problem.
   * @param solutions The output matrix of solutions, which may be empty if no solution exists.
   */
  void solveILP(const n_corners_type nCorners,
                const std::vector<num_type> &boundaryLength,
                const PatchIndexRotationType &patchIndexRotation,
                const FlowMap &flows,
                ILPMatrixType &solutions) {
    solutions.clear();

    // set up the system
    var_label_type nVariables = _patchDBServer->at(nCorners, patchIndexRotation.first).getBoundaryVarVec().size();
    // objective function
    ILPArrayType objectiveFun(nVariables, 1);
    // lower bounds
    ILPArrayType lowerBounds(nVariables, 1);
    // upper bounds
    ILPArrayType upperBounds(nVariables, std::numeric_limits<num_type>::max());
    // objectives to find
    ILPArrayType objectives;
    // the solver
    ILPSolver ilp_solver;

    // the constraints matrix
    ILPMatrixType constraints(boundaryLength.size(), ILPArrayType(nVariables, 0));
    for (var_label_type v = 0; v < nVariables; ++v) {
      ++constraints[_patchDBServer->at(nCorners, patchIndexRotation.first).getBoundaryIDsEnds()[v].first][v];
      ++constraints[_patchDBServer->at(nCorners, patchIndexRotation.first).getBoundaryIDsEnds()[v].second][v];
    }
    // the left- and right-hand-sides of constraints based on the rotation
    ILPArrayType rhs(boundaryLength.size(), 0);
    for (size_t i = 0; i < boundaryLength.size(); ++i)
      rhs[(i + patchIndexRotation.second) % boundaryLength.size()] = boundaryLength[i];

    // try to restrict the problem by adding constraints with the flows
    if (!flows.empty() && *std::min_element(boundaryLength.begin(), boundaryLength.end()) > _minSizeThresholdILP) {
      bool found = ILPSupport::ConstraintILPbyFlows(_patchDBServer->at(nCorners, patchIndexRotation.first), patchIndexRotation.second,
                                                    boundaryLength, flows, objectiveFun, lowerBounds, upperBounds, constraints, rhs,
                                                    _flowConstraintsILPtolerance);
      if (!found) {
        return;
      }
    }

    // generate the problem
    ilp_solver.newProblem(objectiveFun, std::vector<ILPVarType>(1, ILPVarType::INTEGER), lowerBounds, upperBounds);
    // set the system
    ilp_solver.addConstraints(constraints, rhs, rhs);
    // solve!
    ilp_solver.solve(solutions, objectives);
  }

  /**
   * @brief compareSidesLength compares if a boundary sides length is greater or equal to another one starting from nRotations.
   * It corresponds to:
   * for all sides, lBoundaryLenVec[side] >= rBoundaryLenVec[side + nRotations], where the sum is modulo rBoundaryLenVec.size().
   * @param lBoundaryLenVec The left boundary sides length.
   * @param rBoundaryLenVec The right boundary sides length.
   * @param nRotations The starting side for lBoundaryLenVec.
   * @return true if lBoundaryLenVec is greater or equal to rBoundaryLenVec, false otherwise.
   */
  bool compareSidesLength(const std::vector<num_type> &lBoundaryLenVec, const std::vector<num_type> &rBoundaryLenVec,
                          const n_corners_type nRotations) const {
    for (size_t side = 0; side < lBoundaryLenVec.size() && side < rBoundaryLenVec.size(); ++side)
      if (lBoundaryLenVec[side] < rBoundaryLenVec[(side + nRotations) % rBoundaryLenVec.size()])
        return false;
    return true;
  }

  /**
   * @brief isUnique checks if a given (expanded) patch is new or has been already found.
   * @param expandedTopology The topology string of the expanded patch.
   * @return true if is the first, false otherwise.
   */
  bool isUnique(const std::string &expandedTopology) {
#if __cplusplus >= 201103L || _MSC_VER >= 1800
    Hasher::HashValue hashValue = Hasher::StringToHashValue(expandedTopology);
    if (_uniquePatches.find(hashValue) == _uniquePatches.end()) {
      _uniquePatches.insert(hashValue);
#else
    if (_uniquePatches.find(expandedTopology) == _uniquePatches.end()) {
      _uniquePatches.insert(expandedTopology);
#endif
      return true;
    }
    return false;
  }

  /**
   * @brief addNewPatch pushes a new patch back into the container.
   * @param patch The new Patch to insert.
   */
  void addNewPatch(const Patch &patch) {
    // push patch
    typename PatchContainer::iterator pos = std::lower_bound(_patches.begin(), _patches.end(), patch, _compare);
    _patches.insert(pos, patch);
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief addNewPatch pushes a new patch back into the container.
   * @param patch The new ConvexPatch to insert.
   */
  void addNewPatch(Patch&& patch) {
    // push patch
    typename PatchContainer::iterator pos = std::lower_bound(_patches.begin(), _patches.end(), patch, _compare);
    _patches.insert(pos, std::move(patch));
  }
#endif

  /**
   * @brief The ILPdata struct stores the problem of a ILP run.
   */
  struct ILPdata {
    ILPArrayType          objectiveFun;
    ILPArrayType          lowerBounds;
    ILPMatrixType         constraints;
    ILPArrayType          rhs;
    ILPdata              *parent;
    std::deque<ILPdata *> children;
  };


  PatchDBServerPointer      _patchDBServer;                 ///< DB cache server.
  PatchContainer            _patches;                       ///< Container of patches.
  PatchCompare              _compare;                       ///< Compare object to select sort creteria.
  TopologySetType           _uniquePatches;                 ///< Avoids to duplicate patches.
  valence_type              _minValence;                    ///< Patches having singularities with valence under this threshold are discarded.
  valence_type              _maxValence;                    ///< Patches having singularities with valence above this threshold are discarded.
  num_type                  _maxSingularities;              ///< Patches having more singularities than this threshold are discarded.
  n_corners_type            _maxConcavities;                ///< Patches with more than this number are discarded.
  num_type                  _minSizeThresholdILP;           ///< When the shortest side is longer than this, constraint the problem with the flows.
  u_small_int_type          _flowConstraintsILPtolerance;   ///< A tolerance for flow-driven constraints on ILP.
};

} // end namespace pl

} // end namespace tri
} // end namespace vcg

#endif // PATCH_PRIORITY_LOOKUP_H
