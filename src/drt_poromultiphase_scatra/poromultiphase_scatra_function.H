/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_function.H

 \brief Managing and evaluating of (reaction) functions for poromultiphase_scatra
        problems
        Contains hard-coded derivatives and evaluation routines for the most commonly
        used reaction functions to speed up calculation as opposed to the VariableExprFunctions,
        which use FAD

   \level 3

   \maintainer  Johannes Kremheller
                kremheller@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "../drt_lib/drt_function.H"

#ifndef SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_FUNCTION_H_
#define SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_FUNCTION_H_

namespace POROMULTIPHASESCATRA
{
  // abstract class derived from VariableExprFunction since reaction functions are cast to
  // VariableExprFunction in scatra_reaction_coupling and fluidporo_multiphase_singlereaction
  class PoroMultiPhaseScaTraFunction : public DRT::UTILS::VariableExprFunction
  {
   public:
    PoroMultiPhaseScaTraFunction();

    /*!

    \brief add expression to an existing ExprFunction in order to extend
           it to a vector-valued spatial function.

           Every call to AddExpr adds one more component to the
           vector-valued function.

    \param buf (i) (c-string) expression to be parsed during evaluation of this component
    \param constants   (i) vector containing constants defined by a pair (name and value)

    */
    void AddExpr(std::string buf, std::vector<std::pair<std::string, double>> constants)
    {
      dserror("PoroMultiPhaseScaTraFunction does not implement AddExpr");
    }

    /*!

    \brief evaluate function at given time and position

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param t     (i) The point in time in which the
                     function will be evaluated
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated

    */
    double Evaluate(const int index, const double* x, const double t)
    {
      dserror("PoroMultiPhaseScaTraFunction does not implement Evaluate(index,x,t)");
      return 0.0;
    }

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable

    */
    double Evaluate(const int index, const std::vector<std::pair<std::string, double>>& variables)
    {
      dserror("PoroMultiPhaseScaTraFunction does not implement Evaluate(index,variables)");
      return 0.0;
    }

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    virtual double Evaluate(const int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants) = 0;

    /*!

    \brief evaluate derivatives function at given position in space

    \param index (i) For vector-valued functions, index defines the
                     function-component which should be evaluated
                     For scalar functions, index is always set to 0
    \param x     (i) The point in 3-dimensional space in which the
                     function will be evaluated
    \param t     (i) Absolute time in which the
                     function will be evaluated

    */
    std::vector<double> EvaluateSpatialDerivative(const int index, const double* x, const double t)
    {
      dserror("PoroMultiPhaseScaTraFunction does not implement EvaluateSpatialDerivative");
      std::vector<double> emptyvector;
      return emptyvector;
    }

    /*!

    \brief  at given position in space for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable

    */
    std::vector<double> EvaluateDerivative(
        const int index, const std::vector<std::pair<std::string, double>>& variables)
    {
      dserror("PoroMultiPhaseScaTraFunction does not implement EvaluateDerivative(variables)");
      std::vector<double> emptyvector;
      return emptyvector;
    }

    /*!

    \brief  at given position in space for a given set of variables and constants

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    virtual std::vector<double> EvaluateDerivative(int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants) = 0;

    /*!

    \brief  check for correct order of input parameters variables and constants
            this check is performed only once

    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    virtual void CheckOrder(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants) = 0;

    /*!
    \brief  check if a given name is a variable --> this is just a dummy method overwriting
            the one in VariableExprFunction since we do not need it here

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param varname  (i)  Variable name to be checked

    */

    bool IsVariable(const int index, const std::string& varname) const { return false; }

    /*!
    \brief  add a variable to the parser object --> this is just a dummy method overwriting
            the one in VariableExprFunction since we do not need it here

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param varname  (i)  Variable name

    \param varvalue (i)  Value of the variable

    */

    void AddVariable(const int index, const std::string& varname, double varvalue) { return; }

    /*!
    \brief  parse the underlying expression
    */
    void ParseExpressions()
    {
      dserror("PoroMultiPhaseScaTraFunction does not implement ParseExpressions");
      return;
    }

    // return the number of components --> Reaction functions always have to have only one
    // component
    int NumberComponents() { return 1; }

   protected:
    // function parameters for the specific POROMULTIPHASESCATRA_FUNCTION (have to be read in from
    // input file)
    std::vector<double> myfunct_params_;
    // indicating if the order has already been check for this function
    bool order_checked_;
  };

  // standard growth law for tumor cells <--> IF (with lysis) and pressure dependency:
  // +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  // | TUMOR_GROWTH_LAW_HEAVISIDE | |
  // (gamma_T_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2))*(1-phi2)*porosity*S2
  // - lambda*phi2*porosity*S2 | | with phi1: mass fraction of oxygen, phi2: mass fraction of
  // necrotic tumor cells | | Furthermore, we assume that phase 1: healthy cells, phase2: tumor
  // cells, phase3: IF | | | | (possible) INPUT DEFINITION: | | POROMULTIPHASESCATRA_FUNCTION
  // TUMOR_GROWTH_LAW_HEAVISIDE NUMPARAMS 5 PARAMS gamma_T_growth 9.6e-6 w_nl_crit 2.0e-6
  // w_nl_env 4.2e-6 lambda 0.0 p_t_crit 1.0e9                  |
  // +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

  class TumorGrowthLawHeaviside : public PoroMultiPhaseScaTraFunction
  {
   public:
    /*!


    \brief Default constructor creating empty object. Add the function parameters (read from the
    input file) to the function

    */
    TumorGrowthLawHeaviside(std::vector<std::pair<std::string, double>> funct_params);

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    double Evaluate(const int index, const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  at given position in space for a given set of variables and constants

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    std::vector<double> EvaluateDerivative(int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);


    /*!

    \brief  check for correct order of input parameters variables and constants
            this check is performed only once

    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    void CheckOrder(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);
  };

  // standard necrosis law for tumor growth model:
  // +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  // | NECROSIS_LAW_HEAVISIDE | |
  // (1-phi2)*porosity*S2*(-gamma_t_necr*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside(-(phi1-w_nl_crit)/(w_nl_env-w_nl_crit))+
  // delta_a_t*heaviside(p2-p_t_crit) )               | | with phi1: mass fraction of oxygen, phi2:
  // mass fraction of necrotic tumor cells, S2: volume fraction of tumor cells | | Furthermore, we
  // assume that phase 1: healthy cells, phase2: tumor cells, phase3: IF | | | | (possible) INPUT
  // DEFINITION with exactly 5 function parameters: | | POROMULTIPHASESCATRA_FUNCTION
  // NECROSIS_LAW_HEAVISIDE NUMPARAMS 5 PARAMS gamma_t_necr 9.6e-6 w_nl_crit 2.0e-6 w_nl_env 4.2e-6
  // delta_a_t 0.0 p_t_crit 1.0e9                     |
  // +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

  class NecrosisLawHeaviside : public PoroMultiPhaseScaTraFunction
  {
   public:
    /*!


    \brief Default constructor creating empty object. Add the function parameters (read from the
    input file) to the function

    */
    NecrosisLawHeaviside(std::vector<std::pair<std::string, double>> funct_params);

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    double Evaluate(const int index, const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  at given position in space for a given set of variables and constants

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    std::vector<double> EvaluateDerivative(int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  check for correct order of input parameters variables and constants
            this check is performed only once

    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    void CheckOrder(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);
  };

  // standard oxygen consumption law for tumor growth model:
  // +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  // | OXYGEN_CONSUMPTION_LAW_HEAVISIDE | |
  // porosity*(1-phi2)*S2*(gamma_nl_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2)+
  // gamma_0_nl*sin(pi/2.0*phi1/w_nl_env)) | | with phi1: mass fraction of oxygen, phi2: mass
  // fraction of necrotic tumor cells, S2: volume fraction of tumor cells | | Furthermore, we assume
  // that phase 1: healthy cells, phase2: tumor cells, phase3: IF | | | | (possible) INPUT
  // DEFINITION with exactly 5 function parameters: | | POROMULTIPHASESCATRA_FUNCTION
  // OXYGEN_CONSUMPTION_LAW_HEAVISIDE NUMPARAMS 5 PARAMS gamma_nl_growth 2.4e-7 gamma_0_nl 6e-7
  // w_nl_crit 2.0e-6 w_nl_env 4.2e-6 p_t_crit 1.0e9                   |
  // +--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

  class OxygenConsumptionLawHeaviside : public PoroMultiPhaseScaTraFunction
  {
   public:
    /*!


    \brief Default constructor creating empty object. Add the function parameters (read from the
    input file) to the function

    */
    OxygenConsumptionLawHeaviside(std::vector<std::pair<std::string, double>> funct_params);

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    double Evaluate(const int index, const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  at given position in space for a given set of variables and constants

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    std::vector<double> EvaluateDerivative(int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  check for correct order of input parameters variables and constants
            this check is performed only once

    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    void CheckOrder(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);
  };

  // standard growth law for tumor cells <--> IF (with lysis) and pressure dependency as introduced
  // into balance of mass of oxygen:
  // +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  // | TUMOR_GROWTH_LAW_HEAVISIDE_OXY | | phi1*S2*porosity*((
  // gamma_T_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2))*(1-phi2)-
  // lambda*phi2)      | | with phi1: mass fraction of oxygen, phi2: mass fraction of necrotic tumor
  // cells, S2: volume fraction of tumor cells | | Furthermore, we assume that phase 1: healthy
  // cells, phase2: tumor cells, phase3: IF | | | | (possible) INPUT DEFINITION: | |
  // POROMULTIPHASESCATRA_FUNCTION TUMOR_GROWTH_LAW_HEAVISIDE_OXY NUMPARAMS 5 PARAMS
  // gamma_T_growth 9.6e-6 w_nl_crit 2.0e-6 w_nl_env 4.2e-6 lambda 0.0 p_t_crit 1.0e9              |
  // +-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

  class TumorGrowthLawHeavisideOxy : public PoroMultiPhaseScaTraFunction
  {
   public:
    /*!


    \brief Default constructor creating empty object. Add the function parameters (read from the
    input file) to the function

    */
    TumorGrowthLawHeavisideOxy(std::vector<std::pair<std::string, double>> funct_params);

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    double Evaluate(const int index, const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  at given position in space for a given set of variables and constants

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    std::vector<double> EvaluateDerivative(int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  check for correct order of input parameters variables and constants
            this check is performed only once

    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    void CheckOrder(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);
  };

  // standard growth law for tumor cells <--> IF (with lysis) and pressure dependency as introduced
  // into balance of mass of necrotic cells:
  // +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
  // | TUMOR_GROWTH_LAW_HEAVISIDE_NECRO | | porosity*S2*(((
  // gamma_T_growth*(phi1-w_nl_crit)/(w_nl_env-w_nl_crit)*heaviside((phi1-w_nl_crit)/(w_nl_env-w_nl_crit))*heaviside(p_t_crit-p2))*(1-phi2)
  // - lambda*phi2)*phi2 + lambda*phi2)   | | with phi1: mass fraction of oxygen, phi2: mass
  // fraction of necrotic tumor cells | | Furthermore, we assume that phase 1: healthy cells,
  // phase2: tumor cells, phase3: IF | | | | (possible) INPUT DEFINITION: | |
  // POROMULTIPHASESCATRA_FUNCTION TUMOR_GROWTH_LAW_HEAVISIDE_NECRO NUMPARAMS 5 PARAMS
  // gamma_T_growth 9.6e-6 w_nl_crit 2.0e-6 w_nl_env 4.2e-6 lambda 0.0 p_t_crit 1.0e9 |
  // +---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

  class TumorGrowthLawHeavisideNecro : public PoroMultiPhaseScaTraFunction
  {
   public:
    /*!


    \brief Default constructor creating empty object. Add the function parameters (read from the
    input file) to the function

    */
    TumorGrowthLawHeavisideNecro(std::vector<std::pair<std::string, double>> funct_params);

    /*!

    \brief evaluate function for a given set of variables

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    double Evaluate(const int index, const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  at given position in space for a given set of variables and constants

    \param index    (i)  For vector-valued functions, index defines the
                         function-component which should be evaluated
    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */
    std::vector<double> EvaluateDerivative(int index,
        const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);

    /*!

    \brief  check for correct order of input parameters variables and constants
            this check is performed only once

    \param variables (i) A vector containing a pair (variablename, value) for each variable
    \param constants (i) A vector containing a pair (variablename, value) for each constant

    */

    void CheckOrder(const std::vector<std::pair<std::string, double>>& variables,
        const std::vector<std::pair<std::string, double>>& constants);
  };

}  // namespace POROMULTIPHASESCATRA

#endif /* SRC_DRT_POROMULTIPHASE_SCATRA_POROMULTIPHASE_SCATRA_FUNCTION_H_ */
