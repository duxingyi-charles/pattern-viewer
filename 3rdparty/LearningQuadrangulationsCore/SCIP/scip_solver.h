#ifndef SCIP_SOLVER_H
#define SCIP_SOLVER_H

#include <vector>
#include <limits>
#include <iostream>
#include <sstream>
#include "SCIP/scip/scip.h"
#include "SCIP/scip/scipdefplugins.h"
#include "SCIP/scip/struct_var.h"

/////////////////////////////// code copied from SCIP/scipoptsuite-3.1.0/scip-3.1.0/src/scip/cons_countsols.c
typedef SCIP_Longint         Int;
/* constraint handler properties */
#define CONSHDLR_NAME          "countsols"
/** creates and adds a constraint which cuts off the solution from the feasibility
 *  region
 *
 *  input:
 *  - scip            : SCIP main data structure
 *  - sol             : solution to cut off
 *  - conshdlrdata    : constraint handler data
 */
#define CUTOFF_CONSTRAINT(x) SCIP_RETCODE x (SCIP* scip, SCIP_SOL* sol, SCIP_CONSHDLRDATA* conshdlrdata)
/** constraint handler data */
struct SCIP_ConshdlrData
{
   /* solution data and statistic variables */
   SCIP_SPARSESOL**      solutions;          /**< array to store all solutions */
   int                   nsolutions;         /**< number of solution stored */
   int                   ssolutions;         /**< size of the solution array */
   int                   feasST;             /**< number of non trivial feasible subtrees */
   int                   nDiscardSols;       /**< number of discard solutions */
   int                   nNonSparseSols;     /**< number of non sparse solutions */
   Int                   nsols;              /**< number of solutions */
   CUTOFF_CONSTRAINT((*cutoffSolution));     /**< method for cutting of a solution */

   /* constraint handler parameters */
   SCIP_Longint          sollimit;           /**< counting stops, if the given number of solutions were found (-1: no limit) */
   SCIP_Bool             active;             /**< constraint handler active */
   SCIP_Bool             discardsols;        /**< allow to discard solutions */
   SCIP_Bool             sparsetest;         /**< allow to check for sparse solutions */
   SCIP_Bool             collect;            /**< should the solutions be collected */

   SCIP_Bool             warning;            /**< was the warning messages already posted? */

   /* specific problem data */
   SCIP_HASHMAP*         hashmap;            /**< hashmap to store position of active transformed problem variable in our vars array */
   SCIP_VAR**            allvars;            /**< array containing a copy of all variables before presolving */
   SCIP_VAR**            vars;               /**< array containing a copy of all active variables (after presolving) */
   int                   nallvars;           /**< number of all variables in the problem */
   int                   nvars;              /**< number of all active variables in the problem */
   SCIP_Bool             continuous;         /**< are there continuous variables */
};
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief The SCIP_Solver class is a ILP solver which is able to count and get all the counted solutions.
 * It uses the SCIP library with the SoPlex LP solver. See @ref http://scip.zib.de
 */
template < typename ScalarType >
class SCIP_Solver {
public:
  /**
   * @brief The VarType enum defines the type of a variable, which can be binary, integer or continuous.
   */
  enum VarType {
    BINARY            = SCIP_VARTYPE_BINARY,      ///< Binary variable in {0, 1}.
    INTEGER           = SCIP_VARTYPE_INTEGER,     ///< Integer variable in {lb, ..., ub}.
    IMPLICIT_INTEGER  = SCIP_VARTYPE_IMPLINT,     ///< Continuous variable that is always integer.
    CONTINUOUS        = SCIP_VARTYPE_CONTINUOUS   ///< Continuous variable in [lb, ub].
  };

  /**
   * @brief SCIP_Solver is the default constructor.
   */
  SCIP_Solver() : _scip(0) { }

  /**
   * @brief SCIP_Solver Safe copy constructor.
   * @param ilp
   */
  SCIP_Solver(const SCIP_Solver &ilp) : _scip(0) {
    *this = ilp;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief SCIP_Solver Move constructor.
   * @param ilp
   */
  SCIP_Solver(SCIP_Solver&& ilp) : _scip(0) {
    *this = ilp;
  }
#endif

  /**
   * @brief operator = Safe assignment operator.
   * @return
   */
  SCIP_Solver & operator = (const SCIP_Solver &/*ilp*/) {
    release();
    return *this;
  }

#if __cplusplus >= 201103L || _MSC_VER >= 1800
  /**
   * @brief operator = Move assignment operator.
   * @param ilp
   * @return
   */
  SCIP_Solver & operator = (SCIP_Solver&& ilp) {
    if (&ilp != this) {
      release();
      _scip = ilp._scip;
      ilp._scip = 0;
      _variables = std::move(ilp._variables);
      _constraints = std::move(ilp._constraints);
      _constraintsMtx = std::move(ilp._constraintsMtx);
      _lhs = std::move(ilp._lhs);
      _rhs = std::move(ilp._rhs);
    }
    return *this;
  }
#endif

  /**
   * @brief ~SCIP_Solver is the destructor.
   */
  ~SCIP_Solver() {
    release();
  }

  /**
   * @brief newProblem constructs a new ILP with objective, whose variables have lowerBounds and upperBounds.
   *
   * A new set of objective.size() variables is created. Variables not involved in the objective must have
   * coefficient 0.
   *
   * @param objective Vector of coefficients of the objective function.
   * @param variableTypes defines the type of each variable.
   * @param lowerBounds The lower bound of each variable.
   * @param upperBounds The upper bound of each variable.
   * @param toMaximize true if the objective must be maximixed, false otherwise.
   */
  void newProblem(const std::vector<ScalarType> &objective,
                  const std::vector<VarType> &variableTypes = std::vector<VarType>(1, VarType::CONTINUOUS),
                  const std::vector<ScalarType> &lowerBounds = std::vector<ScalarType>(1, 0),
                  const std::vector<ScalarType> &upperBounds = std::vector<ScalarType>(1, std::numeric_limits<ScalarType>::max()),
                  const bool toMaximize = false) {
    // check input
    if (objective.empty()) {
      std::cout << "Objective function coefficients vector can not be empty." << std::endl;
      return;
    }
    if (lowerBounds.empty()) {
      std::cout << "Lower bounds vector can not be empty." << std::endl;
      return;
    }
    if (upperBounds.empty()) {
      std::cout << "Upper bounds vector can not be empty." << std::endl;
      return;
    }
    if (variableTypes.empty()) {
      std::cout << "Variable type vector can not be empty." << std::endl;
      return;
    }

    // initialize
    initialize();

    size_t nVariables = objective.size();
    std::vector<ScalarType> lb = lowerBounds;
    if (lb.size() != nVariables)
      lb.resize(nVariables, lowerBounds.back());
    std::vector<ScalarType> ub = upperBounds;
    if (ub.size() != nVariables)
      ub.resize(nVariables, upperBounds.back());
    std::vector<VarType> vt = variableTypes;
    if (vt.size() != nVariables)
      vt.resize(nVariables, variableTypes.back());

    // create the problem
    SCIP_CALL_ABORT(SCIPcreateProb(_scip, "Problem", 0, 0, 0, 0, 0, 0, 0));
    // set the objective sense
    SCIP_CALL_ABORT(SCIPsetObjsense(_scip, toMaximize ? SCIP_OBJSENSE_MAXIMIZE : SCIP_OBJSENSE_MINIMIZE));
    // create and add variables
    std::string str;
    char c = 'x';
    std::stringstream stream;
    _variables.resize(nVariables, 0);
    SCIP_Real lobn = 0, upbn = 0;
    for (size_t i = 0; i < nVariables; ++i) {
      // define the name of the variable
      str.clear();
      stream.str(str);
      stream << c << i;
      str = stream.str();
      // set lower bound
      lobn = lb[i] == -std::numeric_limits<ScalarType>::max() ? -SCIPinfinity(_scip) : static_cast<SCIP_Real>(lb[i]);
      // set upper bound
      upbn = ub[i] == std::numeric_limits<ScalarType>::max() ? SCIPinfinity(_scip) : static_cast<SCIP_Real>(ub[i]);
      // create a new variable
      SCIP_CALL_ABORT(SCIPcreateVar(_scip, &_variables[i], str.c_str(), lobn, upbn, objective[i],
        static_cast<SCIP_VARTYPE>(vt[i]), TRUE, FALSE, 0, 0, 0, 0, 0));
      // add it
      SCIP_CALL_ABORT(SCIPaddVar(_scip, _variables[i]));
    }
  }
  
  /**
   * @brief addConstraints adds rows to the system of constraints.
   *
   * SCIP manages constraints of type: lhs <= cX <= rhs
   * where X is a vector of variables, c a vector of row coefficients,
   * lhs and rhs are the left-hand-side and right-hand-side, respectively.
   *
   * @param coefficients Matrix of coefficients. It must have a number of columns equal to the number of variables.
   * @param lhs Left-hand-side.
   * @param rhs Right-hand-side.
   */
  void addConstraints(const std::vector< std::vector<ScalarType> > &coefficients,
                      const std::vector<ScalarType> &lhs,
                      const std::vector<ScalarType> &rhs) {
    // check input
    if (_scip == 0)
      return;
    if (coefficients.empty()) {
      std::cout << "No constraints added. These can not be empty." << std::endl;
      return;
    }
    for (size_t i = 0; i < coefficients.size(); ++i)
      if (coefficients[i].size() != _variables.size()) {
        std::cout << "No constraints added. Wrong number of variables." << std::endl;
        return;
      }
    if (lhs.size() != coefficients.size()) {
      std::cout << "No constraints added. Wrong LHS." << std::endl;
      return;
    }
    if (rhs.size() != coefficients.size()) {
      std::cout << "No constraints added. Wrong RHS." << std::endl;
      return;
    }
    // copy constraints
    size_t oldSize = _constraints.size();
    _constraints.resize(oldSize + coefficients.size(), 0);
    _constraintsMtx.insert(_constraintsMtx.end(), coefficients.begin(), coefficients.end());
    _lhs.insert(_lhs.end(), lhs.begin(), lhs.end());
    _rhs.insert(_rhs.end(), rhs.begin(), rhs.end());
    // build SCIP constraints
    std::string str;
    std::stringstream stream;
    SCIP_Real l, r, t;
    for (size_t i = 0; i < coefficients.size(); ++i) {
      // define the name of the constraint
      str.clear();
      stream.str(str);
      stream << "Constraint_" << oldSize + i;
      str = stream.str();
      l = lhs[i] == -std::numeric_limits<ScalarType>::max() ? -SCIPinfinity(_scip) : lhs[i];
      r = rhs[i] == std::numeric_limits<ScalarType>::max() ? SCIPinfinity(_scip) : rhs[i];
      if (l > r) {
        t = l;
        r = l;
        l = t;
      }
      // create the constraint
      SCIP_CALL_ABORT(SCIPcreateConsLinear(_scip, &_constraints[oldSize + i], str.c_str(), 0, 0, 0, l, r,
        TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE));
      // add the variables' coefficients to the constraint
      for (size_t j = 0; j < _variables.size(); ++j)
        SCIP_CALL_ABORT(SCIPaddCoefLinear(_scip, _constraints[oldSize + i], _variables[j], coefficients[i][j]));
      // add the constraint
      SCIP_CALL_ABORT(SCIPaddCons(_scip, _constraints[oldSize + i]));
    }
  }

  /**
   * @brief solve searches the optimal solution of the current problem, if any.
   * @param solution The vector of variable's values.
   * @param objectiveVal The optimal objective value.
   * @return 1 if an optimal solution is found, 0 otherwise.
   */
  int solve(std::vector<ScalarType> &solution, ScalarType &objectiveVal) {
    solution.clear();
    // check status
    if (_scip == 0)
      return 0;

    SCIP_CALL_ABORT(SCIPsolve(_scip));
    SCIP_SOL *sol = SCIPgetBestSol(_scip);
    if (sol != 0) {
      objectiveVal = static_cast<ScalarType>(SCIPgetSolOrigObj(_scip, sol));
      solution.resize(_variables.size());
      for (size_t v = 0; v < _variables.size(); ++v)
        solution[v] = static_cast<ScalarType>(SCIPgetSolVal(_scip, sol, _variables[v]));
      release();
      return 1;
    }
    release();
    return 0;
  }

  /**
   * @brief solve counts and collects all the feasible solutions of the current problem.
   * @param solutions Matrix with a number of columns equal to the number of variables. Each row is a solution.
   * @param objectiveVals Vector of objective values, one for each solution.
   * @return The number of counted solutions.
   */
  int solve(std::vector< std::vector<ScalarType> > &solutions, std::vector<ScalarType> &objectiveVals) {
    solutions.clear();
    objectiveVals.clear();
    // check status
    if (_scip == 0)
      return 0;

    // set default params
    SCIP_CALL_ABORT(SCIPsetParamsCountsols(_scip));
    // set all solutions collectable (meaning keep them in memory)
    SCIP_CALL_ABORT(SCIPsetBoolParam(_scip, "constraints/countsols/collect", TRUE));

    // solve and count
    SCIP_CALL_ABORT(SCIPcount(_scip));
    int nSols = getSolutions(solutions, objectiveVals);
    release();
    return nSols;
  }

private:
  /**
   * @brief initialize creates the solver and its initial structures.
   */
  void initialize() {
    // first release
    release();
    // create a SCIP object
    SCIP_CALL_ABORT(SCIPcreate(&_scip));
    // load plugins
    SCIP_CALL_ABORT(SCIPincludeDefaultPlugins(_scip));
    SCIP_CALL_ABORT(SCIPsetMessagehdlr(_scip, 0));
  }

  /**
   * @brief release frees the structures needed by the solver.
   */
  void release() {
    if (_scip == 0)
      return;
    releaseVariables();
    releaseConstraints();
    SCIP_CALL_ABORT(SCIPfree(&_scip));
    _scip = 0;
    _lhs.clear();
    _rhs.clear();
  }

  /**
   * @brief releaseVariables frees the structures for the variables.
   */
  void releaseVariables() {
    if (_scip == 0)
      return;
    for (size_t i = 0; i < _variables.size(); ++i)
      SCIP_CALL_ABORT(SCIPreleaseVar(_scip, &_variables[i]));
    _variables.clear();
  }

  /**
   * @brief releaseConstraints frees the structures for the constraints.
   */
  void releaseConstraints() {
    if (_scip == 0)
      return;
    for (size_t i = 0; i < _constraints.size(); ++i)
      SCIP_CALL_ABORT(SCIPreleaseCons(_scip, &_constraints[i]));
    _constraints.clear();
  }

  /**
   * @brief getSolutions collects all the counted solutions, if any.
   *
   * @note This code is adapted from the SCIP code in the cons_countsols.c file.
   *
   * @param solutions Matrix with a number of columns equal to the number of variables. Each row is a solution.
   * @param objectiveVals Vector of objective values, one for each solution.
   * @return The number of counted solutions.
   */
  int getSolutions(std::vector< std::vector<ScalarType> > &solutions, std::vector<ScalarType> &objectiveVals) {
    // reset input
    solutions.clear();
    objectiveVals.clear();

    SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;
    int nsparsesols;

    SCIP_Bool valid = FALSE;
    SCIP_Longint nsols = SCIPgetNCountedSols(_scip, &valid);

    // find the countsols constraint handler
    conshdlr = SCIPfindConshdlr(_scip, CONSHDLR_NAME);
    assert(conshdlr != 0);

    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != 0);

    nsparsesols = conshdlrdata->nsolutions;

    if (!valid || nsols == 0 || nsparsesols == 0) {
      return 0;
    } else {
      SCIP_SPARSESOL** sparsesols;
      SCIP_VAR** origvars;
      SCIP_VAR** allvars;
      int norigvars;
      int nactivevars;

      // get sparse solutions defined over the active variables
      nactivevars = conshdlrdata->nvars;
      sparsesols = conshdlrdata->solutions;

      // get original problem variables
      SCIP_CALL_ABORT(SCIPallocBufferArray(_scip, &origvars, SCIPgetNOrigVars(_scip)));

      norigvars = 0;
      for (int v = 0; v < SCIPgetNOrigVars(_scip); ++v)
        if (SCIPvarGetType(SCIPgetOrigVars(_scip)[v]) != SCIP_VARTYPE_CONTINUOUS) {
        origvars[norigvars] = SCIPgetOrigVars(_scip)[v];
        ++norigvars;
        }
      assert(norigvars == conshdlrdata->nallvars);
      assert((size_t)norigvars == _variables.size());

      // resize output
      // duplicate buffer for original variables
      SCIP_CALL_ABORT(SCIPduplicateBufferArray(_scip, &allvars, conshdlrdata->allvars, norigvars));

      // expand and copy solution
      SCIP_SPARSESOL* sparsesol;
      SCIP_VAR** vars;
      SCIP_Real* scalars;
      SCIP_Longint* sol;
      int solInd = 0;

      // get memory to store active solution
      SCIP_CALL_ABORT(SCIPallocBufferArray(_scip, &sol, nactivevars + 1));
      SCIP_CALL_ABORT(SCIPallocBufferArray(_scip, &vars, nactivevars + 1));
      SCIP_CALL_ABORT(SCIPallocBufferArray(_scip, &scalars, nactivevars + 1));

      // resize the output
      solutions.resize(nsols);
      for (size_t i = 0; i < solutions.size(); ++i)
        solutions[i].resize(_variables.size());
      objectiveVals.resize(nsols);

      // loop over all sparse solutions
      for (int s = 0; s < nsparsesols; ++s) {
        sparsesol = sparsesols[s];
        assert(sparsesol != 0);
        assert(SCIPsparseSolGetNVars(sparsesol) == nactivevars);

        // get first solution of the sparse solution
        SCIPsparseSolGetFirstSol(sparsesol, sol, nactivevars);

        // loop over the current sparse solution
        do {
          SCIP_Real objval = 0.0;

          // extract active variables
          for (int v = 0; v < norigvars; ++v) {
            SCIP_Real constant;
            SCIP_Real realvalue;
            int requiredsize;
            int nvars;
            int idx;

            vars[0] = allvars[v];     // this is the culprit!
            scalars[0] = 1.0;
            nvars = 1;
            constant = 0.0;

            SCIP_CALL_ABORT(SCIPgetProbvarLinearSum(_scip, vars, scalars, &nvars, norigvars, &constant, &requiredsize, TRUE));
            assert(requiredsize <= norigvars);
            assert(nvars <= nactivevars);

            realvalue = constant;

            for (int i = 0; i < nvars; ++i) {
              assert(SCIPhashmapExists(conshdlrdata->hashmap, vars[i]));
              idx = ((int)(size_t)SCIPhashmapGetImage(conshdlrdata->hashmap, vars[i])) - 1;
              assert(0 <= idx && idx < nactivevars);
              assert(conshdlrdata->vars[idx] == vars[i]);

              objval += SCIPvarGetObj(vars[i]) * sol[idx];
              realvalue += scalars[i] * sol[idx];
            }
            assert(SCIPisIntegral(_scip, realvalue));

            // store variable's value
            solutions[solInd][v] = static_cast<ScalarType>(realvalue);
          }

          // transform objective value into original problem space
          objectiveVals[solInd] = SCIPretransformObj(_scip, objval);

          // next solution
          ++solInd;
        } while (SCIPsparseSolGetNextSol(sparsesol, sol, nactivevars));
      }

      // free buffer arrays
      SCIPfreeBufferArray(_scip, &scalars);
      SCIPfreeBufferArray(_scip, &vars);
      SCIPfreeBufferArray(_scip, &sol);
      SCIPfreeBufferArray(_scip, &allvars);
      SCIPfreeBufferArray(_scip, &origvars);

      // returns the number of solutions
      return nsols;
    }
  }


  SCIP                                   *_scip;            ///< The main SCIP structure.
  std::vector<SCIP_VAR *>                 _variables;       ///< Vector of SCIP variables.
  std::vector<SCIP_CONS *>                _constraints;     ///< Vector of SCIP constraints.
  std::vector< std::vector<ScalarType> >  _constraintsMtx;  ///< Total input constraints.
  std::vector<ScalarType>                 _lhs;             ///< Total input left-hand-sides.
  std::vector<ScalarType>                 _rhs;             ///< Total input right-hand-sides.
};

#endif // SCIP_SOLVER_H
