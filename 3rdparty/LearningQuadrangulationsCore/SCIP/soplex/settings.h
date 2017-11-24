/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the class library                   */
/*       SoPlex --- the Sequential object-oriented simPlex.                  */
/*                                                                           */
/*    Copyright (C) 1996-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SoPlex is distributed under the terms of the ZIB Academic Licence.       */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SoPlex; see the file COPYING. If not email to soplex@zib.de.  */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file  settings.h
 * @brief Settings class
 */

#ifndef _SETTINGS_H_
#define _SETTINGS_H_

#include <string.h>

#include "spxdefines.h"
#include "soplex.h"

namespace soplex
{

/**@class Settings
 * @brief   Class of parameter settings
 * @ingroup Settings
 */
class Settings
{

public:

   //**@name Parameters */
   //@{

   /// boolean parameters
   typedef enum
   {
      /// should lifting be used to reduce range of nonzero matrix coefficients?
      LIFTING = 0,

      /// should LP be transformed to equality form before a rational solve?
      EQTRANS = 1,

      /// number of boolean parameters
      BOOLPARAM_COUNT = 2
   } BoolParam;

   /// integer parameters
   typedef enum
   {
      /// objective sense
      OBJSENSE = 0,

      /// type of computational form, i.e., column or row representation
      REPRESENTATION = 1,

      /// type of algorithm, i.e., enter or leave
      ALGORITHM = 2,

      /// type of LU update
      FACTOR_UPDATE_TYPE = 3,

      /// maximum number of updates without fresh factorization
      FACTOR_UPDATE_MAX = 4,

      /// iteration limit (-1 if unlimited)
      ITERLIMIT = 5,

      /// refinement limit (-1 if unlimited)
      REFLIMIT = 6,

      /// stalling refinement limit (-1 if unlimited)
      STALLREFLIMIT = 7,

      /// display frequency
      DISPLAYFREQ = 8,

      /// verbosity level
      VERBOSITY = 9,

      /// type of simplifier
      SIMPLIFIER = 10,

      /// type of scaler
      SCALER = 11,

      /// type of starter used to create crash basis
      STARTER = 12,

      /// type of pricer
      PRICER = 13,

      /// type of ratio test
      RATIOTESTER = 14,

      /// mode for synchronizing real and rational LP
      SYNCMODE = 15,

      /// mode for reading LP files
      READMODE = 16,

      /// mode for iterative refinement strategy
      SOLVEMODE = 17,

      /// mode for a posteriori feasibility checks
      CHECKMODE = 18,

      /// mode for hyper sparse pricing
      HYPER_PRICING = 19,

      /// number of integer parameters
      INTPARAM_COUNT = 20
   } IntParam;

   /// values for parameter OBJSENSE
   enum
   {
      /// minimization
      OBJSENSE_MINIMIZE = -1,

      /// maximization
      OBJSENSE_MAXIMIZE = 1
   };

   /// values for parameter REPRESENTATION
   enum
   {
      /// automatic choice according to number of rows and columns
      REPRESENTATION_AUTO = 0,

      /// column representation Ax - s = 0, lower <= x <= upper, lhs <= s <= rhs
      REPRESENTATION_COLUMN = 1,

      /// row representation (lower,lhs) <= (x,Ax) <= (upper,rhs)
      REPRESENTATION_ROW = 2
   };

   /// values for parameter ALGORITHM
   enum
   {
      /// entering algorithm, i.e., primal simplex for column and dual simplex for row representation
      ALGORITHM_ENTER = 0,

      /// leaving algorithm, i.e., dual simplex for column and primal simplex for row representation
      ALGORITHM_LEAVE = 1
   };

   /// values for parameter FACTOR_UPDATE_TYPE
   enum
   {
      /// product form update
      FACTOR_UPDATE_TYPE_ETA = 0,

      /// Forrest-Tomlin type update
      FACTOR_UPDATE_TYPE_FT = 1
   };

   /// values for parameter VERBOSITY
   enum
   {
      /// only error output
      VERBOSITY_ERROR = 0,

      /// only error and warning output
      VERBOSITY_WARNING = 1,

      /// only error, warning, and debug output
      VERBOSITY_DEBUG = 2,

      /// standard verbosity level
      VERBOSITY_NORMAL = 3,

      /// high verbosity level
      VERBOSITY_HIGH = 4,

      /// full verbosity level
      VERBOSITY_FULL = 5
   };

   /// values for parameter SIMPLIFIER
   enum
   {
      /// no simplifier
      SIMPLIFIER_OFF = 0,

      /// automatic choice
      SIMPLIFIER_AUTO = 1
   };

   /// values for parameter SCALER
   enum
   {
      /// no scaler
      SCALER_OFF = 0,

      /// equilibrium scaling on rows or columns
      SCALER_UNIEQUI = 1,

      /// equilibrium scaling on rows and columns
      SCALER_BIEQUI = 2,

      /// geometric mean scaling on rows and columns, max 1 round
      SCALER_GEO1 = 3,

      /// geometric mean scaling on rows and columns, max 8 rounds
      SCALER_GEO8 = 4
   };

   /// values for parameter STARTER
   enum
   {
      /// slack basis
      STARTER_OFF = 0,

      /// greedy crash basis weighted by objective, bounds, and sides
      STARTER_WEIGHT = 1,

      /// crash basis from a greedy solution
      STARTER_SUM = 2,

      /// generic solution-based crash basis
      STARTER_VECTOR = 3
   };

   /// values for parameter PRICER
   enum
   {
      /// automatic pricer
      PRICER_AUTO = 0,

      /// Dantzig pricer
      PRICER_DANTZIG = 1,

      /// partial multiple pricer based on Dantzig pricing
      PRICER_PARMULT = 2,

      /// devex pricer
      PRICER_DEVEX = 3,

      /// steepest edge pricer with initialization to unit norms
      PRICER_QUICKSTEEP = 4,

      /// steepest edge pricer with exact initialization of norms
      PRICER_STEEP = 5
   };

   /// values for parameter RATIOTESTER
   enum
   {
      /// textbook ratio test without stabilization
      RATIOTESTER_TEXTBOOK = 0,

      /// standard Harris ratio test
      RATIOTESTER_HARRIS = 1,

      /// modified Harris ratio test
      RATIOTESTER_FAST = 2,

      /// bound flipping ratio test for long steps in the dual simplex
      RATIOTESTER_BOUNDFLIPPING = 3
   };

   /// values for parameter SYNCMODE
   enum
   {
      /// store only real LP
      SYNCMODE_ONLYREAL = 0,

      /// automatic sync of real and rational LP
      SYNCMODE_AUTO = 1,

      /// user sync of real and rational LP
      SYNCMODE_MANUAL = 2
   };

   /// values for parameter READMODE
   enum
   {
      /// standard floating-point parsing
      READMODE_REAL = 0,

      /// rational parsing
      READMODE_RATIONAL = 1
   };

   /// values for parameter SOLVEMODE
   enum
   {
      /// apply standard floating-point algorithm
      SOLVEMODE_REAL = 0,

      /// decide depending on tolerances whether to apply iterative refinement
      SOLVEMODE_AUTO = 1,

      /// force iterative refinement
      SOLVEMODE_RATIONAL = 2
   };

   /// values for parameter CHECKMODE
   enum
   {
      /// floating-point check
      CHECKMODE_REAL = 0,

      /// decide according to READMODE
      CHECKMODE_AUTO = 1,

      /// rational check
      CHECKMODE_RATIONAL = 2
   };

   /// values for parameter HYPER_PRICING
   enum
   {
      /// never
      HYPER_PRICING_OFF = 0,

      /// decide according to problem size
      HYPER_PRICING_AUTO = 1,

      /// always
      HYPER_PRICING_ON = 2
   };

   /// real parameters
   typedef enum
   {
      /// primal feasibility tolerance
      FEASTOL = 0,

      /// dual feasibility tolerance
      OPTTOL = 1,

      /// general zero tolerance
      EPSILON_ZERO = 2,

      /// zero tolerance used in factorization
      EPSILON_FACTORIZATION = 3,

      /// zero tolerance used in update of the factorization
      EPSILON_UPDATE = 4,

      /// pivot zero tolerance used in factorization
      EPSILON_PIVOT = 5,

      /// infinity threshold
      INFTY = 6,

      /// time limit in seconds (INFTY if unlimited)
      TIMELIMIT = 7,

      /// lower limit on objective value
      OBJLIMIT_LOWER = 8,

      /// upper limit on objective value
      OBJLIMIT_UPPER = 9,

      /// working tolerance for feasibility in floating-point solver during iterative refinement
      FPFEASTOL = 10,

      /// working tolerance for optimality in floating-point solver during iterative refinement
      FPOPTTOL = 11,

      /// maximum increase of scaling factors between refinements
      MAXSCALEINCR = 12,

      /// lower threshold in lifting (nonzero matrix coefficients with smaller absolute value will be reformulated)
      LIFTMINVAL = 13,

      /// upper threshold in lifting (nonzero matrix coefficients with larger absolute value will be reformulated)
      LIFTMAXVAL = 14,

      /// sparse pricing threshold (\#violations < dimension * SPARSITY_THRESHOLD activates sparse pricing)
      SPARSITY_THRESHOLD = 15,

      /// number of real parameters
      REALPARAM_COUNT = 16
   } RealParam;

   /// returns boolean parameter value
   bool boolParam(const BoolParam param) const;

   /// returns integer parameter value
   int intParam(const IntParam param) const;

   /// returns real parameter value
   Real realParam(const RealParam param) const;

   /// returns current parameter settings
   const Settings& settings() const;

   /// sets boolean parameter value; returns true on success
   bool setBoolParam(const BoolParam param, const bool value, const bool quiet = false, const bool init = false);

   /// sets integer parameter value; returns true on success
   bool setIntParam(const IntParam param, const int value, const bool quiet = false, const bool init = false);

   /// sets real parameter value; returns true on success
   bool setRealParam(const RealParam param, const Real value, const bool quiet = false, const bool init = false);

   /// sets parameter settings; returns true on success
   bool setSettings(const Settings& newSettings, const bool quiet = false, const bool init = false);

   /// print non-default parameter values
   void printUserSettings();

   /// writes settings file; returns true on success
   bool saveSettingsFile(const char* filename, const bool onlyChanged = false) const;

   /// reads settings file; returns true on success
   bool loadSettingsFile(const char* filename);

   /// parses one setting string and returns true on success; note that string is modified
   bool parseSettingsString(char* line);

   //@}

   /// default constructor initializing default settings
   Settings(bool init = false);

   /// copy constructor
   Settings(const Settings& settings);

   /// assignment operator
   Settings& operator=(const Settings& settings);

private:

   /// array of names for boolean parameters
   static std::string _boolParamName[BOOLPARAM_COUNT];

   /// array of names for integer parameters
   static std::string _intParamName[INTPARAM_COUNT];

   /// array of names for real parameters
   static std::string _realParamName[REALPARAM_COUNT];

   /// array of descriptions for boolean parameters
   static std::string _boolParamDescription[BOOLPARAM_COUNT];

   /// array of descriptions for integer parameters
   static std::string _intParamDescription[INTPARAM_COUNT];

   /// array of descriptions for real parameters
   static std::string _realParamDescription[REALPARAM_COUNT];

   /// array of default values for boolean parameters
   bool _boolParamDefault[BOOLPARAM_COUNT];

   /// array of default values for integer parameters
   int _intParamDefault[INTPARAM_COUNT];

   /// array of default values for real parameters
   Real _realParamDefault[REALPARAM_COUNT];

   /// array of lower bounds for int parameter values
   int _intParamLower[INTPARAM_COUNT];

   /// array of upper bounds for int parameter values
   int _intParamUpper[INTPARAM_COUNT];

   /// array of lower bounds for real parameter values
   Real _realParamLower[REALPARAM_COUNT];

   /// array of upper bounds for real parameter values
   Real _realParamUpper[REALPARAM_COUNT];

   /// have static arrays been initialized?
   bool _defaultsAndBoundsInitialized;

   /// array of current boolean parameter values
   bool _boolParamValues[BOOLPARAM_COUNT];

   /// array of current integer parameter values
   int _intParamValues[INTPARAM_COUNT];

   /// array of current real parameter values
   Real _realParamValues[REALPARAM_COUNT];


};
}

#endif // _SETTINGS_H_
