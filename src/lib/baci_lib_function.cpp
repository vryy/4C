/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluating of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#include "baci_lib_function.H"

#include "baci_io.H"
#include "baci_lib_function_manager.H"
#include "baci_lib_functionvariables.H"
#include "baci_lib_linedefinition.H"

#include <Sacado.hpp>

#include <utility>

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



  /// sets the values of the variables in second derivative of expression
  template <int dim>
  void SetValuesInExpressionSecondDeriv(
      const std::vector<Teuchos::RCP<DRT::UTILS::FunctionVariable>>& variables, const double* x,
      const double t,
      std::map<std::string, Sacado::Fad::DFad<Sacado::Fad::DFad<double>>>& variable_values)
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

    if constexpr (dim == 1)
    {
      variable_values.emplace("x", xfad);
      variable_values.emplace("y", 0);
      variable_values.emplace("z", 0);
    }
    else if constexpr (dim == 2)
    {
      variable_values.emplace("x", xfad);
      variable_values.emplace("y", yfad);
      variable_values.emplace("z", 0);
    }
    if constexpr (dim == 3)
    {
      variable_values.emplace("x", xfad);
      variable_values.emplace("y", yfad);
      variable_values.emplace("z", zfad);
    }

    // set temporal variable
    variable_values.emplace("t", tfad);

    // set the values of the variables at time t
    for (unsigned int i = 0; i < variables.size(); ++i)
    {
      variable_values.emplace(variables[i]->Name(), fadvectvars[i]);
    }
  }


  /// evaluate an expression and assemble to the result vector
  std::vector<double> EvaluateAndAssembleExpressionToResultVector(
      const std::map<std::string, Sacado::Fad::DFad<double>>& variables,
      const std::size_t component,
      std::vector<Teuchos::RCP<DRT::UTILS::SymbolicExpression<double>>> expr,
      const std::map<std::string, double>& constant_values)
  {
    // number of variables
    auto numvariables = static_cast<int>(variables.size());

    // evaluate the expression
    Sacado::Fad::DFad<double> fdfad = expr[component]->FirstDerivative(variables, constant_values);

    // resulting vector
    std::vector<double> res(numvariables);

    // fill the result vector
    for (int i = 0; i < numvariables; i++) res[i] = fdfad.dx(i);

    return res;
  }

  /// modifies the component to zero in case the expression is of size one
  std::size_t FindModifiedComponent(const std::size_t component,
      const std::vector<Teuchos::RCP<DRT::UTILS::SymbolicExpression<double>>>& expr)
  {
    return (expr.size() == 1) ? 0 : component;
  }

  //! throw an error if a constant given in the input file is a primary variables
  template <typename T>
  void AssertValidInput(const std::map<std::string, T>& variable_values,
      const std::vector<std::pair<std::string, double>>& constants_from_input)
  {
    const bool all_constants_from_input_valid =
        std::all_of(constants_from_input.begin(), constants_from_input.end(),
            [&](const auto& var_name) { return variable_values.count(var_name.first) == 0; });

    if (!all_constants_from_input_valid)
    {
      const auto join_keys = [](const auto& map)
      {
        return std::accumulate(map.begin(), map.end(), std::string(),
            [](const std::string& acc, const auto& v)
            { return acc.empty() ? v.first : acc + ", " + v.first; });
      };

      dserror(
          "It is not allowed to set primary variables of your problem as constants in the "
          "VARFUNCTION.\n\n"
          "Variables passed to Evaluate: %s \n"
          "Constants from Input: %s",
          join_keys(variable_values).c_str(), join_keys(constants_from_input).c_str());
    }
  }
}  // namespace



template <int dim>
Teuchos::RCP<DRT::UTILS::FunctionOfAnything> DRT::UTILS::TryCreateSymbolicFunctionOfAnything(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs)
{
  if (function_line_defs.size() != 1) return Teuchos::null;

  const auto& function_lin_def = function_line_defs.front();

  if (function_lin_def->HaveNamed("VARFUNCTION"))
  {
    std::string component;
    function_lin_def->ExtractString("VARFUNCTION", component);

    std::vector<std::pair<std::string, double>> constants;
    if (function_lin_def->HaveNamed("CONSTANTS"))
    {
      function_lin_def->ExtractPairOfStringAndDoubleVector("CONSTANTS", constants);
    }

    return Teuchos::rcp(new DRT::UTILS::SymbolicFunctionOfAnything<dim>(component, constants));
  }
  else
  {
    return Teuchos::null;
  }
}



template <int dim>
Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime> DRT::UTILS::TryCreateSymbolicFunctionOfSpaceTime(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs)
{
  // evaluate the maximum component and the number of variables
  int maxcomp = 0;
  int maxvar = -1;
  bool found_function_of_space_time(false);
  for (const auto& ith_function_lin_def : function_line_defs)
  {
    ith_function_lin_def->ExtractInt("COMPONENT", maxcomp);
    ith_function_lin_def->ExtractInt("VARIABLE", maxvar);
    if (ith_function_lin_def->HaveNamed("SYMBOLIC_FUNCTION_OF_SPACE_TIME"))
      found_function_of_space_time = true;
  }

  if (!found_function_of_space_time) return Teuchos::null;

  // evaluate the number of rows used for the definition of the variables
  std::size_t numrowsvar = function_line_defs.size() - maxcomp - 1;

  // define a vector of strings
  std::vector<std::string> functstring(maxcomp + 1);

  // read each row where the components of the i-th function are defined
  for (int n = 0; n <= maxcomp; ++n)
  {
    // update the current row
    Teuchos::RCP<DRT::INPUT::LineDefinition> functcomp = function_line_defs[n];

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
    Teuchos::RCP<DRT::INPUT::LineDefinition> line = function_line_defs[maxcomp + j];

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

  return Teuchos::rcp(new SymbolicFunctionOfSpaceTime<dim>(functstring, functvarvector));
}

template <int dim>
DRT::UTILS::SymbolicFunctionOfSpaceTime<dim>::SymbolicFunctionOfSpaceTime(
    const std::vector<std::string>& expressions,
    std::vector<Teuchos::RCP<FunctionVariable>> variables)
    : variables_(std::move(variables))
{
  for (const auto& expression : expressions)
  {
    {
      auto symbolicexpression =
          Teuchos::rcp(new DRT::UTILS::SymbolicExpression<ValueType>(expression));
      expr_.push_back(symbolicexpression);
    }
  }
}

template <int dim>
double DRT::UTILS::SymbolicFunctionOfSpaceTime<dim>::Evaluate(
    const double* x, const double t, const std::size_t component)
{
  std::size_t component_mod = FindModifiedComponent(component, expr_);

  if (component_mod < 0 || component_mod >= expr_.size())
    dserror("There are %d expressions but tried to access component %d", expr_.size(), component);

  // create map for variables
  std::map<std::string, double> variable_values;

  // set spatial variables
  if constexpr (dim > 0) variable_values.emplace("x", x[0]);
  if constexpr (dim > 1) variable_values.emplace("y", x[1]);
  if constexpr (dim > 2) variable_values.emplace("z", x[2]);

  // set temporal variable
  variable_values.emplace("t", t);

  // set the values of the variables at time t
  for (const auto& variable : variables_)
  {
    variable_values.emplace(variable->Name(), variable->Value(t));
  }

  // evaluate F = F ( x, y, z, t, v1, ..., vn )
  return expr_[component_mod]->Value(variable_values);
}

template <int dim>
std::vector<double> DRT::UTILS::SymbolicFunctionOfSpaceTime<dim>::EvaluateSpatialDerivative(
    const double* x, const double t, const std::size_t component)
{
  std::size_t component_mod = FindModifiedComponent(component, expr_);

  if (component_mod < 0 || component_mod >= expr_.size())
    dserror("There are %d expressions but tried to access component %d", expr_.size(), component);


  // variables
  std::map<std::string, Sacado::Fad::DFad<Sacado::Fad::DFad<double>>> variable_values;

  SetValuesInExpressionSecondDeriv<dim>(variables_, x, t, variable_values);

  // The expression evaluates to an FAD object for up to second derivatives
  SecondDerivativeType fdfad = expr_[component_mod]->SecondDerivative(variable_values, {});

  // Here we return the first spatial derivatives given by FAD component 0, 1 and 2
  return {fdfad.dx(0).val(), fdfad.dx(1).val(), fdfad.dx(2).val()};
}

template <int dim>
std::vector<double> DRT::UTILS::SymbolicFunctionOfSpaceTime<dim>::EvaluateTimeDerivative(
    const double* x, const double t, const unsigned deg, const std::size_t component)
{
  // result vector
  std::vector<double> res(deg + 1);

  std::size_t component_mod = FindModifiedComponent(component, expr_);

  // variables
  std::map<std::string, Sacado::Fad::DFad<Sacado::Fad::DFad<double>>> variable_values;

  SetValuesInExpressionSecondDeriv<dim>(variables_, x, t, variable_values);

  // FAD object for evaluation of derivatives
  Sacado::Fad::DFad<Sacado::Fad::DFad<double>> fdfad;

  const int number_of_arguments = 4;

  // add the value at time t
  res[0] = Evaluate(x, t, component);

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // evaluation of derivatives

    fdfad = expr_[component_mod]->SecondDerivative(variable_values, {});

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
DRT::UTILS::SymbolicFunctionOfAnything<dim>::SymbolicFunctionOfAnything(
    const std::string& component, std::vector<std::pair<std::string, double>> constants)
    : constants_from_input_(std::move(constants))
{
  // build the parser for the function evaluation
  auto symbolicexpression = Teuchos::rcp(new DRT::UTILS::SymbolicExpression<double>(component));

  // save the parsers
  expr_.push_back(symbolicexpression);
}



template <int dim>
double DRT::UTILS::SymbolicFunctionOfAnything<dim>::Evaluate(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const std::size_t component)
{
  // create map for variables
  std::map<std::string, double> variable_values;

  // convert vector of pairs to map variables_values
  std::copy(
      variables.begin(), variables.end(), std::inserter(variable_values, variable_values.begin()));

  // add constants
  std::copy(
      constants.begin(), constants.end(), std::inserter(variable_values, variable_values.end()));

  if (constants_from_input_.size() != 0)
  {
    // check if constants_from_input are valid
    AssertValidInput<double>(variable_values, constants_from_input_);

    // add constants from input
    std::copy(constants_from_input_.begin(), constants_from_input_.end(),
        std::inserter(variable_values, variable_values.end()));
  }

  // evaluate the function and return the result
  return expr_[component]->Value(variable_values);
}


template <int dim>
std::vector<double> DRT::UTILS::SymbolicFunctionOfAnything<dim>::EvaluateDerivative(
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants, const std::size_t component)
{
  auto variables_FAD = ConvertVariableValuesToFADObjects(variables);

  std::map<std::string, Sacado::Fad::DFad<double>> variable_values;

  std::map<std::string, double> constant_values;

  // convert vector of pairs to map variables_values
  std::copy(variables_FAD.begin(), variables_FAD.end(),
      std::inserter(variable_values, variable_values.begin()));

  // add constants
  std::copy(
      constants.begin(), constants.end(), std::inserter(constant_values, constant_values.begin()));


  if (constants_from_input_.size() != 0)
  {
    // check if constants_from_input are valid
    AssertValidInput<Sacado::Fad::DFad<double>>(variable_values, constants_from_input_);

    // add constants from input
    std::copy(constants_from_input_.begin(), constants_from_input_.end(),
        std::inserter(constant_values, constant_values.end()));
  }
  return EvaluateAndAssembleExpressionToResultVector(
      variable_values, component, expr_, constant_values);
}


// explicit instantiations

template class DRT::UTILS::SymbolicFunctionOfSpaceTime<1>;
template class DRT::UTILS::SymbolicFunctionOfSpaceTime<2>;
template class DRT::UTILS::SymbolicFunctionOfSpaceTime<3>;

template class DRT::UTILS::SymbolicFunctionOfAnything<1>;
template class DRT::UTILS::SymbolicFunctionOfAnything<2>;
template class DRT::UTILS::SymbolicFunctionOfAnything<3>;

template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime>
DRT::UTILS::TryCreateSymbolicFunctionOfSpaceTime<1>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime>
DRT::UTILS::TryCreateSymbolicFunctionOfSpaceTime<2>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfSpaceTime>
DRT::UTILS::TryCreateSymbolicFunctionOfSpaceTime<3>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);

template Teuchos::RCP<DRT::UTILS::FunctionOfAnything>
DRT::UTILS::TryCreateSymbolicFunctionOfAnything<1>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfAnything>
DRT::UTILS::TryCreateSymbolicFunctionOfAnything<2>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
template Teuchos::RCP<DRT::UTILS::FunctionOfAnything>
DRT::UTILS::TryCreateSymbolicFunctionOfAnything<3>(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs);
