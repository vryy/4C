/*! \file
\brief Interface for functions of time
\level 0
*/

#include "lib_function_of_time.H"

#include "lib_linedefinition.H"
#include "lib_symbolic_expression.H"

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
}  // namespace

DRT::UTILS::SymbolicFunctionOfTime::SymbolicFunctionOfTime(
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

double DRT::UTILS::SymbolicFunctionOfTime::Evaluate(
    const double time, const std::size_t component) const
{
  std::map<std::string, ValueType> variable_values;

  variable_values.emplace("t", time);

  for (const auto& variable : variables_)
  {
    variable_values.emplace(variable->Name(), variable->Value(time));
  }

  return expr_[component]->Value(variable_values);
}

double DRT::UTILS::SymbolicFunctionOfTime::EvaluateDerivative(
    const double time, const std::size_t component) const
{
  std::map<std::string, FirstDerivativeType> variable_values;

  // define FAD variables
  // argument is only time
  const int number_of_arguments = 1;
  // we consider a function of the type F = F ( t, v1(t), ..., vn(t) )
  const int fad_size = number_of_arguments + static_cast<int>(variables_.size());
  FirstDerivativeType tfad(fad_size, 0, time);

  std::vector<FirstDerivativeType> fadvectvars(variables_.size());
  for (int i = 0; i < static_cast<int>(variables_.size()); ++i)
  {
    fadvectvars[i] =
        FirstDerivativeType(fad_size, number_of_arguments + i, variables_[i]->Value(time));
    fadvectvars[i].val() = variables_[i]->Value(time);
  }

  // set temporal variable
  variable_values.emplace("t", tfad);

  // set the values of the variables at time t
  for (unsigned int i = 0; i < variables_.size(); ++i)
  {
    variable_values.emplace(variables_[i]->Name(), fadvectvars[i]);
  }

  auto f_dfad = expr_[component]->FirstDerivative(variable_values, {});

  double f_dt = f_dfad.dx(0);
  for (int i = 0; i < static_cast<int>(variables_.size()); ++i)
  {
    f_dt += f_dfad.dx(number_of_arguments + i) * variables_[i]->TimeDerivativeValue(time);
  }

  return f_dt;
}

Teuchos::RCP<DRT::UTILS::FunctionOfTime> DRT::UTILS::TryCreateFunctionOfTime(
    const std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>& function_line_defs)
{
  // evaluate the maximum component and the number of variables
  int maxcomp = 0;
  int maxvar = -1;
  bool found_function_of_time(false);
  for (const auto& ith_function_lin_def : function_line_defs)
  {
    ith_function_lin_def->ExtractInt("COMPONENT", maxcomp);
    ith_function_lin_def->ExtractInt("VARIABLE", maxvar);
    if (ith_function_lin_def->HaveNamed("SYMBOLIC_FUNCTION_OF_TIME")) found_function_of_time = true;
  }

  if (!found_function_of_time) return Teuchos::null;

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
    functcomp->ExtractString("SYMBOLIC_FUNCTION_OF_TIME", functstring[n]);
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

  return Teuchos::rcp(new SymbolicFunctionOfTime(functstring, functvarvector));
}