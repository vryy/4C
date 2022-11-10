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
#include "drt_linedefinition.H"
#include "drt_parser.H"
#include "io.H"

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
  std::vector<double> ExtractTimeVector(const DRT::INPUT::LineDefinition& timevar)
  {
    // read the number of points
    int numpoints;
    timevar.ExtractInt("NUMPOINTS", numpoints);

    // read whether times are defined by number of points or by vector
    bool bynum = timevar.HasString("BYNUM");

    // read respectively create times vector
    std::vector<double> times;
    if (bynum)  // times defined by number of points
    {
      // read the time range
      std::vector<double> timerange;
      timevar.ExtractDoubleVector("TIMERANGE", timerange);

      // create time vector from number of points and time range
      times = CreateTimesFromTimeRange(timerange, numpoints);
    }
    else  // times defined by vector
    {
      timevar.ExtractDoubleVector("TIMES", times);
    }

    // check if the times are in ascending order
    if (!std::is_sorted(times.begin(), times.end()))
      dserror("the TIMES must be in ascending order");

    return times;
  }

  /// converts the values of variables from type double to FAD double and returns the modified
  /// vector of name-value-pairs
  std::vector<std::pair<std::string, Sacado::Fad::DFad<double>>> ConvertVariableValuesToFADObjects(
      const std::vector<std::pair<std::string, double>>& variables)
  {
    // prepare return vector
    std::vector<std::pair<std::string, Sacado::Fad::DFad<double>>> variables_FAD;

    // number of variables
    auto numvariables = static_cast<int>(variables.size());

    // counter for variable numbering
    int counter = 0;

    // set the values of the variables
    for (const auto& [name, value] : variables)
    {
      // FAD object for 1st order derivatives
      Sacado::Fad::DFad<double> varfad(numvariables, counter, value);

      // create name-value-pairs with values now of type FAD double and add to vector
      variables_FAD.emplace_back(name, varfad);

      // update counter
      counter++;
    }
    return variables_FAD;
  }

  /// set the values of the variables or constants in expression or in first derivative of
  /// expression
  template <typename Expression, typename ValueType>
  void SetValuesInExpressionOrExpressionFirstDeriv(std::vector<Expression>& expr_or_exprd,
      const int index, const std::vector<std::pair<std::string, ValueType>>& variables_or_constants)
  {
    // set the values of the variables or constants
    for (const auto& [name, value] : variables_or_constants)
    {
      if (expr_or_exprd[index]->IsVariable(name))
        // set the value in expression
        expr_or_exprd[index]->SetValue(name, value);
    }
  }

  /// sets the values of the variables in second derivative of expression
  template <int dim>
  void SetValuesInExpressionSecondDeriv(
      const Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>>& exprdd,
      const std::vector<Teuchos::RCP<DRT::UTILS::FunctionVariable>>& variables, const double* x,
      const double t)
  {
    static_assert(dim >= 0 && dim <= 3, "Spatial dimension has to be 1, 2 or 3.");

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
      fadvectvars[i] = FAD(fad_size, number_of_arguments + i, variables[i]->Value(t));
      fadvectvars[i].val() =
          Sacado::Fad::DFad<double>(fad_size, number_of_arguments + i, variables[i]->Value(t));
    }

    // initialize spatial variables with zero
    exprdd->SetValue("x", 0);
    exprdd->SetValue("y", 0);
    exprdd->SetValue("z", 0);

    if constexpr (dim > 0) exprdd->SetValue("x", xfad);
    if constexpr (dim > 1) exprdd->SetValue("y", yfad);
    if constexpr (dim > 2) exprdd->SetValue("z", zfad);

    // set temporal variable
    exprdd->SetValue("t", tfad);

    // set the values of the variables at time t
    for (unsigned int i = 0; i < variables.size(); ++i)
    {
      exprdd->SetValue(variables[i]->Name(), fadvectvars[i]);
    }
  }

  /// evaluate an expression and assemble to the result vector
  std::vector<double> EvaluateAndAssembleExpressionToResultVector(
      const std::vector<std::pair<std::string, double>>& variables, const int index,
      std::vector<Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>>> exprd)
  {
    // number of variables
    auto numvariables = static_cast<int>(variables.size());

    // evaluate the expression
    Sacado::Fad::DFad<double> fdfad = exprd[index]->Evaluate();

    // resulting vector
    std::vector<double> res(numvariables);

    // fill the result vector
    for (int i = 0; i < numvariables; i++) res[i] = fdfad.dx(i);

    return res;
  }

  /// check if index is in range of the dimensions of the expression, otherwise throw an error
  void AssertVariableIndexInDimensionOfExpression(const int index,
      const std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>>& expr,
      const std::string& varName)
  {
    if (index > (int)expr.size() - 1 || index < 0)
      dserror(
          "Tried to access variable <%s> for expression at index %d but only %d expressions are "
          "available.",
          varName.c_str(), index, expr.size());
  }

  /// modifies the index to zero in case the expression is of size one
  std::size_t FindModifiedIndex(
      const std::size_t index, const std::vector<Teuchos::RCP<DRT::PARSER::Parser<double>>>& expr)
  {
    std::size_t index_mod = index;

    if (expr.size() == 1)
    {
      index_mod = 0;
    }

    return index_mod;
  }

  template <int dim>
  Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> CreateVariableExprFunction(
      const std::string& component, const std::vector<std::pair<std::string, double>>& constants)
  {
    auto vecfunc = Teuchos::rcp(new DRT::UTILS::VariableExprFunction<dim>());
    vecfunc->AddExpr(component, constants);
    return vecfunc;
  }
}  // namespace

template <int dim>
Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateVariableExprFunction(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int index_current_funct_in_manager)
{
  (void)index_current_funct_in_manager;
  if (function_lin_def->HaveNamed("VARFUNCTION"))
  {
    std::string component;
    function_lin_def->ExtractString("VARFUNCTION", component);

    std::vector<std::pair<std::string, double>> constants;
    if (function_lin_def->HaveNamed("CONSTANTS"))
    {
      function_lin_def->ExtractPairOfStringAndDoubleVector("CONSTANTS", constants);
    }

    return CreateVariableExprFunction<dim>(component, constants);
  }
  else
  {
    return Teuchos::null;
  }
}


template <int dim>
Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateExprFunction(
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs)
{
  // evaluate the maximum component and the number of variables
  int maxcomp = 0;
  int maxvar = -1;
  bool found_function_of_space_time(false);
  for (const auto& ith_function_lin_def : functions_lin_defs)
  {
    ith_function_lin_def->ExtractInt("COMPONENT", maxcomp);
    ith_function_lin_def->ExtractInt("VARIABLE", maxvar);
    if (ith_function_lin_def->HaveNamed("SYMBOLIC_FUNCTION_OF_SPACE_TIME"))
      found_function_of_space_time = true;
  }

  if (!found_function_of_space_time) return Teuchos::null;

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
    functcomp->ExtractString("SYMBOLIC_FUNCTION_OF_SPACE_TIME", functstring[n]);
  }

  std::map<int, std::vector<Teuchos::RCP<FunctionVariable>>> variable_pieces;

  // read each row where the variables of the i-th function are defined
  for (std::size_t j = 1; j <= numrowsvar; ++j)
  {
    // update the current row
    Teuchos::RCP<DRT::INPUT::LineDefinition> line = functions_lin_defs[maxcomp + j];

    // read the number of the variable
    int varid;
    line->ExtractInt("VARIABLE", varid);


    const auto variable = [&line]() -> Teuchos::RCP<DRT::UTILS::FunctionVariable>
    {
      // read the name of the variable
      std::string varname;
      line->ExtractString("NAME", varname);

      // read the type of the variable
      std::string vartype;
      line->ExtractString("TYPE", vartype);

      // read periodicity data
      periodicstruct periodicdata{};

      periodicdata.periodic = line->HasString("PERIODIC");
      if (periodicdata.periodic)
      {
        line->ExtractDouble("T1", periodicdata.t1);
        line->ExtractDouble("T2", periodicdata.t2);
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
        line->ExtractString("DESCRIPTION", description);
        return Teuchos::rcp(new ParsedFunctionVariable(varname, description));
      }
      else if (vartype == "linearinterpolation")
      {
        // read times
        std::vector<double> times = ExtractTimeVector(*line);

        // read values
        std::vector<double> values;
        line->ExtractDoubleVector("VALUES", values);

        return Teuchos::rcp(new LinearInterpolationVariable(varname, times, values, periodicdata));
      }
      else if (vartype == "multifunction")
      {
        // read times
        std::vector<double> times = ExtractTimeVector(*line);

        // read descriptions (strings separated with spaces)
        std::vector<std::string> description_vec;
        line->ExtractStringVector("DESCRIPTION", description_vec);

        // check if the number of times = number of descriptions + 1
        std::size_t numtimes = times.size();
        std::size_t numdescriptions = description_vec.size();
        if (numtimes != numdescriptions + 1)
          dserror("the number of TIMES and the number of DESCRIPTIONs must be consistent");

        return Teuchos::rcp(
            new MultiFunctionVariable(varname, times, description_vec, periodicdata));
      }
      else if (vartype == "fourierinterpolation")
      {
        // read times
        std::vector<double> times = ExtractTimeVector(*line);

        // read values
        std::vector<double> values;
        line->ExtractDoubleVector("VALUES", values);

        return Teuchos::rcp(new FourierInterpolationVariable(varname, times, values, periodicdata));
      }
      else
      {
        dserror("unknown variable type");
        return Teuchos::null;
      }
    }();

    variable_pieces[varid].emplace_back(variable);
  }

  std::vector<Teuchos::RCP<FunctionVariable>> functvarvector;

  for (const auto& [id, pieces] : variable_pieces)
  {
    // exactly one variable piece -> can be added directly
    if (pieces.size() == 1) functvarvector.emplace_back(pieces[0]);
    // multiple pieces make up this variable -> join them in a PiecewiseVariable
    else
    {
      const auto& name = pieces.front()->Name();

      const bool names_of_all_pieces_equal = std::all_of(
          pieces.begin(), pieces.end(), [&name](auto& var) { return var->Name() == name; });
      if (not names_of_all_pieces_equal)
        dserror("Variable %d has a piece-wise definition with inconsistent names.", id);

      functvarvector.emplace_back(Teuchos::rcp(new PiecewiseVariable(name, pieces)));
    }
  }

  return Teuchos::rcp(new ExprFunction<dim>(functstring, functvarvector));
}

template <int dim>
DRT::UTILS::ExprFunction<dim>::ExprFunction(const std::vector<std::string>& expressions,
    std::vector<Teuchos::RCP<FunctionVariable>> variables)
    : variables_(std::move(variables))
{
  const auto addVariablesToParser = [this](auto& parser)
  {
    parser->AddVariable("x", 0);
    parser->AddVariable("y", 0);
    parser->AddVariable("z", 0);
    parser->AddVariable("t", 0);

    for (const auto& var : variables_) parser->AddVariable(var->Name(), 0);
  };

  for (const auto& expression : expressions)
  {
    {
      auto parser = Teuchos::rcp(new DRT::PARSER::Parser<ValueType>(expression));
      addVariablesToParser(parser);
      parser->ParseFunction();
      expr_.push_back(parser);
    }

    {
      auto parser = Teuchos::rcp(new DRT::PARSER::Parser<SecondDerivativeType>(expression));
      addVariablesToParser(parser);
      parser->ParseFunction();
      exprdd_.push_back(parser);
    }
  }
}

template <int dim>
double DRT::UTILS::ExprFunction<dim>::Evaluate(
    const double* x, const double t, const std::size_t component)
{
  std::size_t component_mod = FindModifiedIndex(component, expr_);

  if (component_mod < 0 || component_mod >= expr_.size())
    dserror("There are %d expressions but tried to access index %d", expr_.size(), component);

  // set spatial variables
  if constexpr (dim > 0) expr_[component_mod]->SetValue("x", x[0]);
  if constexpr (dim > 1) expr_[component_mod]->SetValue("y", x[1]);
  if constexpr (dim > 2) expr_[component_mod]->SetValue("z", x[2]);

  // set temporal variable
  expr_[component_mod]->SetValue("t", t);

  // set the values of the variables at time t
  for (const auto& variable : variables_)
  {
    expr_[component_mod]->SetValue(variable->Name(), variable->Value(t));
  }

  // evaluate F = F ( x, y, z, t, v1, ..., vn )
  return expr_[component_mod]->Evaluate();
}

template <int dim>
std::vector<double> DRT::UTILS::ExprFunction<dim>::EvaluateSpatialDerivative(
    const double* x, const double t, const std::size_t component)
{
  std::size_t component_mod = FindModifiedIndex(component, expr_);

  if (component_mod < 0 || component_mod >= expr_.size())
    dserror("There are %d expressions but tried to access component %d", expr_.size(), component);

  SetValuesInExpressionSecondDeriv<dim>(exprdd_[component_mod], variables_, x, t);

  // The expression evaluates to an FAD object for up to second derivatives
  SecondDerivativeType fdfad = exprdd_[component_mod]->Evaluate();

  // Here we return the first spatial derivatives given by FAD index 0, 1 and 2
  return {fdfad.dx(0).val(), fdfad.dx(1).val(), fdfad.dx(2).val()};
}

template <int dim>
std::vector<double> DRT::UTILS::ExprFunction<dim>::EvaluateTimeDerivative(
    const double* x, const double t, const unsigned deg, const std::size_t component)
{
  // result vector
  std::vector<double> res(deg + 1);

  std::size_t component_mod = FindModifiedIndex(component, expr_);

  SetValuesInExpressionSecondDeriv<dim>(exprdd_[component_mod], variables_, x, t);

  // FAD object for evaluation of derivatives
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> fdfad;

  const int number_of_arguments = 4;

  // add the value at time t
  res[0] = Evaluate(x, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // evaluation of derivatives
    fdfad = exprdd_[component_mod]->Evaluate();

    // evaluation of dF/dt applying the chain rule:
    // dF/dt = dF*/dt + sum_i(dF/dvi*dvi/dt)
    double fdfad_dt = fdfad.dx(3).val();                           // 1) dF*/dt
    for (int i = 0; i < static_cast<int>(variables_.size()); ++i)  // 2) sum_i{...}
    {
      fdfad_dt += fdfad.dx(number_of_arguments + i).val() * variables_[i]->TimeDerivativeValue(t);
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
    double fdfad_dt2 = fdfad.dx(3).dx(3);                   // 1) add d(dF*/dt)/dt to d^2F/dt^2
    std::vector<double> fdfad_dt2_term(variables_.size());  // prepare sum_i{...}

    for (int i = 0; i < static_cast<int>(variables_.size()); ++i)
    {
      fdfad_dt2_term[i] = 0;

      fdfad_dt2_term[i] += fdfad.dx(3).dx(number_of_arguments + i);  // ... + d(dF*/dt)/dvi
      fdfad_dt2_term[i] += fdfad.dx(number_of_arguments + i).dx(3);  // ... + d(dF/dvi)/dt

      for (int j = 0; j < static_cast<int>(variables_.size()); ++j)  // prepare + sum_j{...}
      {
        fdfad_dt2_term[i] +=
            fdfad.dx(number_of_arguments + i).dx(number_of_arguments + j) *  // d(dF/dvi)/dvj ...
            variables_[j]->TimeDerivativeValue(t);                           // ... * dvj/dt
      }

      fdfad_dt2_term[i] *= variables_[i]->TimeDerivativeValue(t);  // ... * dvi/dt

      fdfad_dt2_term[i] += fdfad.dx(number_of_arguments + i).val() *  /// ... + dF/dvi ...
                           variables_[i]->TimeDerivativeValue(t, 2);  /// ... * d^2vi/dt^2

      fdfad_dt2 += fdfad_dt2_term[i];  // 2) add sum_i{...} to d^2F/dt^2
    }

    res[2] = fdfad_dt2;
  }

  // error for higher derivatives
  if (deg >= 3)
  {
    dserror("Higher time derivatives than second not supported!");
  }

  // return derivatives
  return res;
}


template <int dim>
void DRT::UTILS::VariableExprFunction<dim>::AddExpr(
    const std::string& buf, const std::vector<std::pair<std::string, double>>& constants)
{
  // do the almost same as the expression function (base class) but do not yet parse!

  // build the parser for the function evaluation
  auto parser = Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));

  // build the parser for the function derivative evaluation
  auto parserd = Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<double>>(buf));

  // add constants
  for (const auto& constant : constants) parser->AddVariable(constant.first, constant.second);
  for (const auto& constant : constants) parserd->AddVariable(constant.first, constant.second);

  // save the parsers
  expr_.push_back(parser);
  exprd_.push_back(parserd);

  isparsed_ = false;
}

template <int dim>
bool DRT::UTILS::VariableExprFunction<dim>::IsVariable(int index, const std::string& varname) const
{
  AssertVariableIndexInDimensionOfExpression(index, expr_, varname);

  return expr_[index]->IsVariable(varname);
}

template <int dim>
void DRT::UTILS::VariableExprFunction<dim>::AddVariable(
    int index, const std::string& varname, double varvalue)
{
  AssertVariableIndexInDimensionOfExpression(index, expr_, varname);

  if (isparsed_) dserror("Function has already been parsed! Variables can no longer be added!");

  expr_[index]->AddVariable(varname, varvalue);
  exprd_[index]->AddVariable(varname, varvalue);
}

template <int dim>
void DRT::UTILS::VariableExprFunction<dim>::ParseExpressions()
{
  // loop over expressions and parse them
  for (auto& parser : expr_) parser->ParseFunction();

  // loop over expressions for derivatives and parse them
  for (auto& parser : exprd_) parser->ParseFunction();

  isparsed_ = true;
}

template <int dim>
double DRT::UTILS::VariableExprFunction<dim>::Evaluate(
    const double* x, const double t, const std::size_t component)
{
  std::vector<std::pair<std::string, double>> variables;
  variables.reserve(dim_);

  // set spatial variables
  if (dim_ < 0 or dim_ > 3) dserror("Problem dimension has to be 1, 2, or 3.");
  if (dim_ > 0) variables.emplace_back("x", x[0]);
  if (dim_ > 1) variables.emplace_back("y", x[1]);
  if (dim_ > 2) variables.emplace_back("z", x[2]);

  // set temporal variable
  variables.emplace_back("t", t);

  return Evaluate(component, variables);
}

template <int dim>
std::vector<double> DRT::UTILS::VariableExprFunction<dim>::EvaluateSpatialDerivative(
    const double* x, const double t, const std::size_t component)
{
  // arguments are: x, y, z, and t
  const int number_of_arguments = 4;
  std::vector<std::pair<std::string, double>> variables(number_of_arguments);

  variables[0] = std::pair<std::string, double>("x", x[0]);
  variables[1] = std::pair<std::string, double>("y", x[1]);
  variables[2] = std::pair<std::string, double>("z", x[2]);
  variables[3] = std::pair<std::string, double>("t", t);

  return EvaluateDerivative(component, variables);
}

template <int dim>
double DRT::UTILS::VariableExprFunction<dim>::Evaluate(
    const int index, const std::vector<std::pair<std::string, double>>& variables)
{
  if (not isparsed_) ParseExpressions();

  SetValuesInExpressionOrExpressionFirstDeriv(expr_, index, variables);

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}

template <int dim>
double DRT::UTILS::VariableExprFunction<dim>::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  if (not isparsed_) ParseExpressions();

  SetValuesInExpressionOrExpressionFirstDeriv(expr_, index, variables);

  SetValuesInExpressionOrExpressionFirstDeriv(expr_, index, constants);

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}

template <int dim>
std::vector<double> DRT::UTILS::VariableExprFunction<dim>::EvaluateDerivative(
    const int index, const std::vector<std::pair<std::string, double>>& variables)
{
  if (not isparsed_) ParseExpressions();

  auto variables_FAD = ConvertVariableValuesToFADObjects(variables);
  SetValuesInExpressionOrExpressionFirstDeriv(exprd_, index, variables_FAD);

  return EvaluateAndAssembleExpressionToResultVector(variables, index, exprd_);
}

template <int dim>
std::vector<double> DRT::UTILS::VariableExprFunction<dim>::EvaluateDerivative(int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  if (not isparsed_) ParseExpressions();

  auto variables_FAD = ConvertVariableValuesToFADObjects(variables);
  SetValuesInExpressionOrExpressionFirstDeriv(exprd_, index, variables_FAD);

  SetValuesInExpressionOrExpressionFirstDeriv(exprd_, index, constants);

  return EvaluateAndAssembleExpressionToResultVector(variables, index, exprd_);
}

// explicit instantiations

template class DRT::UTILS::ExprFunction<1>;
template class DRT::UTILS::ExprFunction<2>;
template class DRT::UTILS::ExprFunction<3>;

template class DRT::UTILS::VariableExprFunction<1>;
template class DRT::UTILS::VariableExprFunction<2>;
template class DRT::UTILS::VariableExprFunction<3>;

template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateExprFunction<1>(
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateExprFunction<2>(
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateExprFunction<3>(
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs);

template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateVariableExprFunction<1>(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int index_current_funct_in_manager);
template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateVariableExprFunction<2>(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int index_current_funct_in_manager);
template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateVariableExprFunction<3>(
    Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def, DRT::UTILS::FunctionManager& manager,
    const int index_current_funct_in_manager);
