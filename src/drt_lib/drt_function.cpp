/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluating of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#include <Sacado.hpp>
#include "drt_function.H"
#include "drt_function_manager.H"
#include "drt_functionvariables.H"
#include "drt_globalproblem.H"
#include "drt_linedefinition.H"
#include "../drt_io/io.H"

namespace
{
  /// creates a vector of times from a given NUMPOINTS and TIMERANGE
  std::vector<double> CreateTimesFromTimeRange(
      const std::vector<double>& timerange, const int& numpoints)
  {
    std::vector<double> times;

    // get initial and final time
    double t_initial = timerange[0];
    double t_final = timerange[1];

    // build the vector of times
    times.push_back(t_initial);
    int n = 0;
    double dt = (t_final - t_initial) / (numpoints - 1);
    while (times[n] + dt <= t_final + 1.0e-14)
    {
      if (times[n] + 2 * dt <= t_final + 1.0e-14)
      {
        times.push_back(times[n] + dt);
      }
      else
      {
        times.push_back(t_final);
      }
      ++n;
    }

    return times;
  }

  /// returns a vector of times either from NUMPOINTS and TIMERANGE or from TIMES of a line
  /// definition
  std::vector<double> returnTimeVector(const Teuchos::RCP<DRT::INPUT::LineDefinition> timevar)
  {
    // read the number of points
    int numpoints;
    timevar->ExtractInt("NUMPOINTS", numpoints);

    // read whether times are defined by number of points or by vector
    bool bynum = timevar->HasString("BYNUM");

    // read times
    std::vector<double> times;
    if (bynum)  // times defined by number of points
    {
      // read the time range
      std::vector<double> timerange;
      timevar->ExtractDoubleVector("TIMERANGE", timerange);

      // create time vector from number of points and time range
      times = CreateTimesFromTimeRange(timerange, numpoints);
    }
    else  // times defined by vector
    {
      timevar->ExtractDoubleVector("TIMES", times);
    }

    // check if the times are in ascending order
    if (!std::is_sorted(times.begin(), times.end()))
      dserror("the TIMES must be in ascending order");

    return times;
  }

  /// set the values of the variables
  void SetVariableValuesInExpressionDeriv(
      const std::vector<std::pair<std::string, double>>& variables,
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>>> exprd,
      const int index)
  {
    // number of variables
    auto numvariables = static_cast<int>(variables.size());

    // counter for variable numbering
    int counter = 0;

    // define Fad object for evaluation
    using FAD = Sacado::Fad::DFad<double>;

    // set the values of the variables
    std::vector<std::pair<std::string, double>>::const_iterator it;
    for (it = variables.begin(); it != variables.end(); it++)
    {
      // for 1st order derivatives
      FAD varfad(numvariables, counter, it->second);
      // set the value in expression
      exprd[index]->SetValue(it->first, varfad);
      // update counter
      counter++;
    }
  }

  /// set the values of the variables in expression
  void SetVariableValuesInExpression(const std::vector<std::pair<std::string, double>>& variables,
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>> expr, const int index)
  {
    // set the values of the variables
    std::vector<std::pair<std::string, double>>::const_iterator it;
    for (it = variables.begin(); it != variables.end(); it++)
    {
      expr[index]->SetValue(it->first, it->second);
    }
  }

  /// set the values of the constants in expression derivative
  void SetConstantsValuesInExpressivDeriv(
      const std::vector<std::pair<std::string, double>>& constants,
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>>> exprd,
      const int index)
  {
    // set the values of the constants
    std::vector<std::pair<std::string, double>>::const_iterator it;
    for (it = constants.begin(); it != constants.end(); it++)
    {
      if (exprd[index]->IsVariable(it->first))
        // set the value in expression
        exprd[index]->SetValue(it->first, it->second);
    }
  }

  /// set the values of the constants in expression
  void SetConstantsValuesInExpression(const std::vector<std::pair<std::string, double>>& constants,
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>> expr, const int index)
  {
    // set the values of the constants
    std::vector<std::pair<std::string, double>>::const_iterator it;
    for (it = constants.begin(); it != constants.end(); it++)
    {
      if (expr[index]->IsVariable(it->first))
        // set the value in expression
        expr[index]->SetValue(it->first, it->second);
    }
  }

  /// check if index is in range of the dimensions of the expression, otherwise throw error
  void AssertVariableIndexInDimensionOfExpression(
      const int index, std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>> expr)
  {
    if (index > (int)expr.size() - 1 || index < 0)
      dserror("Tried to add a variable to a function in a not available dimension.");
  }

  /// needs further description: Finds a colum index depending on a time that is either contained in
  /// variables at a certain row index or not
  /// --> not fully understood yet
  unsigned int FindColIndexForVariableDefinition(
      std::vector<std::vector<Teuchos::RCP<DRT::UTILS::FunctionVariable>>>& variables,
      const int index_row, const double time)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int index_col = 0;
    bool containtime = false;
    while (!containtime)
    {
      containtime = variables[index_row][index_col]->ContainTime(time);
      if (index_col == variables[index_row].size())
      {
        dserror("the variable %d is not defined at time %f", index_row, time);
      }
      if (!containtime)
      {
        ++index_col;
      }
    }

    return index_col;
  }

  /// needs further description: Sets the values of exprdd depending on various inputs
  /// --> not fully understood yet
  void SetValuesOfExprDerivDeriv(std::size_t index_mod,
      std::vector<std::vector<Teuchos::RCP<DRT::UTILS::FunctionVariable>>> variables,
      const double* x, const double t, int dim_,
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>>>
          exprdd)
  {
    // define Fad object for evaluation
    using FAD = Sacado::Fad::DFad<Sacado::Fad::DFad<double>>;

    // define FAD variables
    // arguments are: x, y, z, and t
    const int number_of_arguments = 4;
    // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
    const int fad_size = number_of_arguments + static_cast<int>(variables.size());
    FAD xfad(fad_size, 0, x[0]);
    FAD yfad(fad_size, 1, x[1]);
    FAD zfad(fad_size, 2, x[2]);
    FAD tfad(fad_size, 3, t);

    xfad.val() = Sacado::Fad::DFad<double>(fad_size, 0, x[0]);
    yfad.val() = Sacado::Fad::DFad<double>(fad_size, 1, x[1]);
    zfad.val() = Sacado::Fad::DFad<double>(fad_size, 2, x[2]);
    tfad.val() = Sacado::Fad::DFad<double>(fad_size, 3, t);

    std::vector<FAD> fadvectvars(variables.size());
    for (int i = 0; i < static_cast<int>(variables.size()); ++i)
    {
      // find the right definition of the variable according to the hierarchy
      unsigned int n = FindColIndexForVariableDefinition(variables, i, t);

      fadvectvars[i] = FAD(fad_size, number_of_arguments + i, variables[i][n]->Value(t));
      fadvectvars[i].val() =
          Sacado::Fad::DFad<double>(fad_size, number_of_arguments + i, variables[i][n]->Value(t));
    }

    // set spatial variables
    switch (dim_)
    {
      case 3:
      {
        exprdd[index_mod]->SetValue("x", xfad);
        exprdd[index_mod]->SetValue("y", yfad);
        exprdd[index_mod]->SetValue("z", zfad);
        break;
      }
      case 2:
      {
        exprdd[index_mod]->SetValue("x", xfad);
        exprdd[index_mod]->SetValue("y", yfad);
        exprdd[index_mod]->SetValue("z", 0);
        break;
      }
      case 1:
      {
        exprdd[index_mod]->SetValue("x", xfad);
        exprdd[index_mod]->SetValue("y", 0);
        exprdd[index_mod]->SetValue("z", 0);
        break;
      }
      default:
        dserror("Problem dimension has to be 1, 2, or 3.");
        break;
    }

    // set temporal variable
    exprdd[index_mod]->SetValue("t", tfad);

    // set the values of the variables at time t
    for (unsigned int i = 0; i < variables.size(); ++i)
    {
      exprdd[index_mod]->SetValue(variables[i][0]->Name(), fadvectvars[i]);
    }
  }

  std::size_t FindModifiedIndex(
      const int index, std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>> expr)
  {
    std::size_t index_mod = index;

    if (expr.size() == 1)
    {
      index_mod = 0;
    }

    return index_mod;
  }
}  // namespace

Teuchos::RCP<DRT::UTILS::Function> DRT::UTILS::TryCreateFunctionFunction(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int index_current_funct_in_manager)
{
  if (function_lin_def->HaveNamed("VARFUNCTION"))
  {
    Teuchos::RCP<DRT::UTILS::VariableExprFunction> vecfunc =
        Teuchos::rcp(new DRT::UTILS::VariableExprFunction());

    std::string component;
    function_lin_def->ExtractString("VARFUNCTION", component);

    std::vector<std::pair<std::string, double>> constants;
    if (function_lin_def->HaveNamed("CONSTANTS"))
    {
      function_lin_def->ExtractPairOfStringAndDoubleVector("CONSTANTS", constants);
    }

    vecfunc->AddExpr(component, constants);
    return vecfunc;
  }
  else
  {
    return Teuchos::RCP<DRT::UTILS::VariableExprFunction>(NULL);
  }
}

Teuchos::RCP<DRT::UTILS::Function> DRT::UTILS::TryCreateBasicFunction(
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs)
{
  // define a new vector of functions
  Teuchos::RCP<ExprFunction> vecfunc = Teuchos::rcp(new ExprFunction());

  // evaluate the maximum component and the number of variables
  int maxcomp = 0;
  int maxvar = -1;
  for (const auto& ith_function_lin_def : functions_lin_defs)
  {
    ith_function_lin_def->ExtractInt("COMPONENT", maxcomp);
    ith_function_lin_def->ExtractInt("VARIABLE", maxvar);
  }

  // evaluate the number of rows used for the definition of the variables
  std::size_t numrowsvar = functions_lin_defs.size() - maxcomp - 1;

  // define a vector of strings
  std::vector<std::string> functstring(maxcomp + 1);

  // read each row where the components of the i-th function are defined
  for (int n = 0; n <= maxcomp; ++n)
  {
    // update the current row
    Teuchos::RCP<DRT::INPUT::LineDefinition> functcomp = functions_lin_defs[n];

    // check the validity of the n-th component
    int compid = 0;
    functcomp->ExtractInt("COMPONENT", compid);
    if (compid != n) dserror("expected COMPONENT %d but got COMPONENT %d", n, compid);


    // read the expression of the n-th component of the i-th function
    functcomp->ExtractString("FUNCTION", functstring[n]);
  }

  // define the structure functvarvector
  std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> functvarvector;

  // define the structure functvar
  std::vector<Teuchos::RCP<FunctionVariable>> functvar;

  // define the structure vardef
  Teuchos::RCP<FunctionVariable> vardef;

  int vardefinition = 1;
  int varidold = -1;

  // read each row where the variables of the i-th function are defined
  for (std::size_t j = 1; j <= numrowsvar; ++j)
  {
    // update the current row
    Teuchos::RCP<DRT::INPUT::LineDefinition> timevar = functions_lin_defs[maxcomp + j];

    // read the number of the variable
    int varid;
    timevar->ExtractInt("VARIABLE", varid);

    // evaluate the number of the definition for the variable
    if (varid == varidold)
    {
      ++vardefinition;
    }
    else
    {
      vardefinition = 1;
    }

    // update the old varid
    varidold = varid;

    // read the name of the variable
    std::string varname;
    timevar->ExtractString("NAME", varname);

    // read the type of the variable
    std::string vartype;
    timevar->ExtractString("TYPE", vartype);

    // read periodicity data
    struct periodicstruct periodicdata
    {
    };
    periodicdata.periodic = timevar->HasString("PERIODIC");
    if (periodicdata.periodic)
    {
      timevar->ExtractDouble("T1", periodicdata.t1);
      timevar->ExtractDouble("T2", periodicdata.t2);
    }
    else
    {
      periodicdata.t1 = 0;
      periodicdata.t2 = 0;
    }

    // distinguish the type of the variable
    if (vartype == "expression")
    {
      std::string description;
      timevar->ExtractString("DESCRIPTION", description);
      vardef = Teuchos::rcp(new ParsedFunctionVariable(varname, description));
    }
    else if (vartype == "linearinterpolation")
    {
      // read times
      std::vector<double> times = returnTimeVector(timevar);

      // read values
      std::vector<double> values;
      timevar->ExtractDoubleVector("VALUES", values);

      vardef = Teuchos::rcp(new LinearInterpolationVariable(varname, times, values, periodicdata));
    }
    else if (vartype == "multifunction")
    {
      // read times
      std::vector<double> times = returnTimeVector(timevar);

      // read descriptions (strings separated with spaces)
      std::vector<std::string> description_vec;
      timevar->ExtractStringVector("DESCRIPTION", description_vec);

      // check if the number of times = number of descriptions + 1
      std::size_t numtimes = times.size();
      std::size_t numdescriptions = description_vec.size();
      if (numtimes != numdescriptions + 1)
        dserror("the number of TIMES and the number of DESCRIPTIONs must be consistent");

      vardef =
          Teuchos::rcp(new MultiFunctionVariable(varname, times, description_vec, periodicdata));
    }
    else if (vartype == "fourierinterpolation")
    {
      // read times
      std::vector<double> times = returnTimeVector(timevar);

      // read values
      std::vector<double> values;
      timevar->ExtractDoubleVector("VALUES", values);

      vardef = Teuchos::rcp(new FourierInterpolationVariable(varname, times, values, periodicdata));
    }
    else
    {
      dserror("unknown variable type");
    }

    // insert the variable in the vector of the variables of the function
    if (vardefinition == 1)
    {
      if (varid != 0)
      {
        functvarvector.push_back(functvar);
        functvar.clear();
      }
    }
    functvar.push_back(vardef);

    if (j == numrowsvar)
    {
      functvarvector.push_back(functvar);
    }
  }

  // add the expressions to the function vector
  for (int n = 0; n <= maxcomp; ++n)
  {
    vecfunc->AddExpr(functstring[n], functvarvector);
  }

  return vecfunc;
}


DRT::UTILS::ExprFunction::ExprFunction()
{
  dim_ = DRT::Problem::Instance()->NDim();
  expr_.clear();
  exprdd_.clear();
  variables_.clear();
  isparsed_ = false;
}

void DRT::UTILS::ExprFunction::AddExpr(const std::string& buf,
    const std::vector<std::vector<Teuchos::RCP<FunctionVariable>>>& variables)
{
  variables_ = variables;

  Teuchos::RCP<DRT::PARSER::Parser<double>> parser =
      Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));
  parser->AddVariable("x", 0);
  parser->AddVariable("y", 0);
  parser->AddVariable("z", 0);
  parser->AddVariable("t", 0);

  Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>> parserdd =
      Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>(buf));
  parserdd->AddVariable("x", 0);
  parserdd->AddVariable("y", 0);
  parserdd->AddVariable("z", 0);
  parserdd->AddVariable("t", 0);

  // add variables from the defined VARIABLES
  for (unsigned int i = 0; i < variables.size(); ++i)
  {
    parser->AddVariable(variables_[i][0]->Name(), 0);
    parserdd->AddVariable(variables_[i][0]->Name(), 0);
  }

  // save parsers
  expr_.push_back(parser);
  exprdd_.push_back(parserdd);

  isparsed_ = false;
}

double DRT::UTILS::ExprFunction::Evaluate(const int index, const double* x, double t)
{
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  if (not isparsed_) ParseExpressions();

  std::size_t index_mod = FindModifiedIndex(index, expr_);

  // set spatial variables
  if (dim_ > 0) expr_[index_mod]->SetValue("x", x[0]);
  if (dim_ > 1) expr_[index_mod]->SetValue("y", x[1]);
  if (dim_ > 2) expr_[index_mod]->SetValue("z", x[2]);

  // set temporal variable
  expr_[index_mod]->SetValue("t", t);

  // set the values of the variables at time t
  for (unsigned int i = 0; i < variables_.size(); ++i)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int n = FindColIndexForVariableDefinition(variables_, i, t);

    expr_[index_mod]->SetValue(variables_[i][0]->Name(), variables_[i][n]->Value(t));
  }

  // evaluation of F = F ( x, y, z, t, v1, ..., vn )
  return expr_[index_mod]->Evaluate();
}

std::vector<double> DRT::UTILS::ExprFunction::EvaluateSpatialDerivative(
    const int index, const double* x, const double t)
{
  // parse expression if not already parsed
  if (not isparsed_) ParseExpressions();

  std::size_t index_mod = FindModifiedIndex(index, expr_);

  SetValuesOfExprDerivDeriv(index_mod, variables_, x, t, dim_, exprdd_);

  // evaluation of derivatives
  using FAD = Sacado::Fad::DFad<Sacado::Fad::DFad<double>>;
  FAD fdfad = exprdd_[index_mod]->Evaluate();

  // result vector
  std::vector<double> res(3, 0.0);
  for (int d = 0; d < 3; ++d)
  {
    res[d] = fdfad.dx(d).val();
  }

  // return derivatives
  return res;
}

std::vector<double> DRT::UTILS::ExprFunction::EvaluateTimeDerivative(
    const int index, const double* x, const double t, const unsigned deg)
{
  // resulting vector holding
  std::vector<double> res(deg + 1);

  // add the value at time t
  res[0] = Evaluate(index, x, t);

  // parse expression if not already parsed
  if (not isparsed_) ParseExpressions();

  std::size_t index_mod = FindModifiedIndex(index, expr_);

  SetValuesOfExprDerivDeriv(index_mod, variables_, x, t, dim_, exprdd_);

  using FAD = Sacado::Fad::DFad<Sacado::Fad::DFad<double>>;
  FAD fdfad;
  const int number_of_arguments = 4;

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // evaluation of derivatives
    fdfad = exprdd_[index_mod]->Evaluate();

    // evaluation of dF/dt applying the chain rule:
    // dF/dt = dF*/dt + sum_i(dF/dvi*dvi/dt)
    double fdfad_dt = fdfad.dx(3).val();
    for (int i = 0; i < static_cast<int>(variables_.size()); ++i)
    {
      // find the right definition of the variable according to the hierarchy
      unsigned int n = FindColIndexForVariableDefinition(variables_, i, t);

      fdfad_dt +=
          fdfad.dx(number_of_arguments + i).val() * variables_[i][n]->TimeDerivativeValue(t);
    }

    res[1] = fdfad_dt;
  }

  // add the 2nd time derivative at time t
  if (deg >= 2)
  {
    // evaluation of d^2F/dt^2 applying the chain rule:
    // d^2F/dt^2 = d(dF*/dt)/dt + sum_i{
    //                                [d(dF*/dt)/dvi + d(dF/dvi)/dt + sum_j[d(dF/dvi)/dvj * dvj/dt]]
    //                                * dvi/dt +
    //                                 + dF/dvi * d^2vi/dt^2
    //                                }
    double fdfad_dt2 = fdfad.dx(3).dx(3);
    std::vector<double> fdfad_dt2_term(variables_.size());

    for (int i = 0; i < static_cast<int>(variables_.size()); ++i)
    {
      fdfad_dt2_term[i] = 0;

      unsigned int n = FindColIndexForVariableDefinition(variables_, i, t);

      fdfad_dt2_term[i] += fdfad.dx(3).dx(number_of_arguments + i);
      fdfad_dt2_term[i] += fdfad.dx(number_of_arguments + i).dx(3);
      for (int j = 0; j < static_cast<int>(variables_.size()); ++j)
      {
        unsigned int m = FindColIndexForVariableDefinition(variables_, j, t);

        fdfad_dt2_term[i] += fdfad.dx(number_of_arguments + i).dx(number_of_arguments + j) *
                             variables_[j][m]->TimeDerivativeValue(t);
      }
      fdfad_dt2_term[i] *= variables_[i][n]->TimeDerivativeValue(t);
      fdfad_dt2_term[i] +=
          fdfad.dx(number_of_arguments + i).val() * variables_[i][n]->TimeDerivativeValue(t, 2);
      fdfad_dt2 += fdfad_dt2_term[i];
    }

    res[2] = fdfad_dt2;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  return res;
}

bool DRT::UTILS::ExprFunction::IsVariable(const int index, const std::string& varname) const
{
  AssertVariableIndexInDimensionOfExpression(index, expr_);

  return expr_[index]->IsVariable(varname);
}

void DRT::UTILS::ExprFunction::AddVariable(
    const int index, const std::string& varname, double varvalue)
{
  AssertVariableIndexInDimensionOfExpression(index, expr_);

  if (isparsed_) dserror("Function has already been parsed! Variables can no longer be added!");

  expr_[index]->AddVariable(varname, varvalue);
  exprdd_[index]->AddVariable(varname, varvalue);
}

void DRT::UTILS::ExprFunction::ParseExpressions()
{
  // define iterators
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>>::iterator it;
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>>>::
      iterator itdd;

  // loop over expressions and parse them
  for (it = expr_.begin(); it != expr_.end(); it++) (*it)->ParseFunction();

  // loop over expressions for derivatives and parse them
  for (itdd = exprdd_.begin(); itdd != exprdd_.end(); itdd++) (*itdd)->ParseFunction();

  isparsed_ = true;
}


void DRT::UTILS::VariableExprFunction::AddExpr(
    const std::string& buf, const std::vector<std::pair<std::string, double>>& constants)
{
  // do the almost same as the expression function (base class) but do not yet parse!

  // build the parser for the function evaluation
  Teuchos::RCP<DRT::PARSER::Parser<double>> parser =
      Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));

  // build the parser for the function derivative evaluation
  Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>> parserd =
      Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<double>>(buf));

  // add constants
  for (const auto& constant : constants) parser->AddVariable(constant.first, constant.second);
  for (const auto& constant : constants) parserd->AddVariable(constant.first, constant.second);

  // save the parsers
  expr_.push_back(parser);
  exprd_.push_back(parserd);

  isparsed_ = false;
}

bool DRT::UTILS::VariableExprFunction::IsVariable(int index, const std::string& varname) const
{
  AssertVariableIndexInDimensionOfExpression(index, expr_);

  return expr_[index]->IsVariable(varname);
}

void DRT::UTILS::VariableExprFunction::AddVariable(
    int index, const std::string& varname, double varvalue)
{
  AssertVariableIndexInDimensionOfExpression(index, expr_);

  if (isparsed_) dserror("Function has already been parsed! Variables can no longer be added!");

  expr_[index]->AddVariable(varname, varvalue);
  exprd_[index]->AddVariable(varname, varvalue);
}

void DRT::UTILS::VariableExprFunction::ParseExpressions()
{
  // define iterators
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>>::iterator it;
  std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>>>::iterator itd;

  // loop over expressions and parse them
  for (it = expr_.begin(); it != expr_.end(); it++) (*it)->ParseFunction();

  // loop over expressions for derivatives and parse them
  for (itd = exprd_.begin(); itd != exprd_.end(); itd++) (*itd)->ParseFunction();

  isparsed_ = true;
}

double DRT::UTILS::VariableExprFunction::Evaluate(
    const int index, const std::vector<std::pair<std::string, double>>& variables)
{
  if (not isparsed_) ParseExpressions();

  SetVariableValuesInExpression(variables, expr_, index);

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}

double DRT::UTILS::VariableExprFunction::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  if (not isparsed_) ParseExpressions();

  SetVariableValuesInExpression(variables, expr_, index);

  SetConstantsValuesInExpression(constants, expr_, index);

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}

std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateDerivative(
    const int index, const std::vector<std::pair<std::string, double>>& variables)
{
  if (not isparsed_) ParseExpressions();

  SetVariableValuesInExpressionDeriv(variables, exprd_, index);

  // number of variables
  auto numvariables = static_cast<int>(variables.size());

  // evaluate the expression
  Sacado::Fad::DFad<double> fdfad = exprd_[index]->Evaluate();

  // resulting vector
  std::vector<double> res(numvariables);

  // fill the result vector
  for (int i = 0; i < numvariables; i++) res[i] = fdfad.dx(i);

  return res;
}

std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateDerivative(int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  if (not isparsed_) ParseExpressions();

  SetVariableValuesInExpressionDeriv(variables, exprd_, index);

  SetConstantsValuesInExpressivDeriv(constants, exprd_, index);

  // number of variables
  auto numvariables = static_cast<int>(variables.size());

  // evaluate the expression
  Sacado::Fad::DFad<double> fdfad = exprd_[index]->Evaluate();

  // resulting vector
  std::vector<double> res(numvariables);

  // fill the result vector
  for (int i = 0; i < numvariables; i++) res[i] = fdfad.dx(i);

  return res;
}

double DRT::UTILS::VariableExprFunction::Evaluate(const int index, const double* x, const double t)
{
  std::vector<std::pair<std::string, double>> variables;
  variables.reserve(dim_);

  switch (dim_)
  {
    case 3:
    {
      variables.emplace_back("x", x[0]);  //-x_[index]
      variables.emplace_back("y", x[1]);  //-y_[index]
      variables.emplace_back("z", x[2]);  //-z_[index]
    }
    case 2:
    {
      variables.emplace_back("x", x[0]);  //-x_[index]
      variables.emplace_back("y", x[1]);  //-y_[index]
    }
    case 1:
    {
      variables.emplace_back("x", x[0]);  //-x_[index]
    }
  }

  variables.emplace_back("t", t);

  return Evaluate(index, variables);
}

std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateSpatialDerivative(
    int index, const double* x, const double t)
{
  // arguments are: x, y, z, and t
  const int number_of_arguments = 4;
  std::vector<std::pair<std::string, double>> variables(number_of_arguments);

  variables[0] = std::pair<std::string, double>("x", x[0]);
  variables[1] = std::pair<std::string, double>("y", x[1]);
  variables[2] = std::pair<std::string, double>("z", x[2]);
  variables[3] = std::pair<std::string, double>("t", t);

  return EvaluateDerivative(index, variables);
}
