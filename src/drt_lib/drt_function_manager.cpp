/*----------------------------------------------------------------------*/
/*! \file

\brief Managing of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#include "drt_function.H"
#include "function_of_time.H"
#include "drt_linedefinition.H"
#include "../drt_fluid/fluid_functions.H"
#include "../drt_fluid_xfluid/xfluid_functions.H"
#include "../drt_fluid_xfluid/xfluid_functions_combust.H"
#include "../drt_structure_new/str_functions.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_function.H"
#include "drt_function_library.H"
#include "../drt_io/io.H"
#include "drt_globalproblem.H"


namespace
{
  template <typename T>
  using CreateFunctionType = Teuchos::RCP<T> (*)(
      Teuchos::RCP<DRT::INPUT::LineDefinition>, DRT::UTILS::FunctionManager&, const int);

  template <typename T>
  using CreateMultiLineFunctionType = Teuchos::RCP<T> (*)(
      std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>);

  /**
   * Utility function that takes a function object returning a Teuchos::RCP<T> and erases its return
   * type via std::any. In addition, if the returned object would be Teuchos::null, discard it and
   * return an empty std::any instead.
   */
  template <typename T>
  std::function<std::any(
      Teuchos::RCP<DRT::INPUT::LineDefinition>, DRT::UTILS::FunctionManager&, const int)>
  WrapFunction(CreateFunctionType<T> fun)
  {
    return [fun](Teuchos::RCP<DRT::INPUT::LineDefinition> linedef, DRT::UTILS::FunctionManager& fm,
               const int index) -> std::any
    {
      Teuchos::RCP<T> created = fun(linedef, fm, index);
      if (created == Teuchos::null)
        return {};
      else
        return created;
    };
  }

  template <typename T>
  std::function<std::any(std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>)>
  WrapMultiLineFunction(CreateMultiLineFunctionType<T> fun)
  {
    return [fun](std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> multiline_def) -> std::any
    {
      Teuchos::RCP<T> created = fun(multiline_def);
      if (created == Teuchos::null)
        return {};
      else
        return created;
    };
  }

  template <int dim>
  void FillFunctions(DRT::INPUT::DatFileReader& reader, std::vector<std::any>& functions,
      DRT::UTILS::FunctionManager& functionManager)
  {
    Teuchos::RCP<DRT::INPUT::Lines> lines = DRT::UTILS::FunctionManager::ValidFunctionLines();

    // Read FUNCT sections starting from FUNCT1 until the first empty one is encountered.
    // This implies that the FUNCT sections must form a contiguous range in the input file.
    // Otherwise, the read fails later.
    for (int funct_suffix = 1;; ++funct_suffix)
    {
      // Read lines belonging to section FUNCT<i> in the input file
      std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> section_line_defs =
          lines->Read(reader, funct_suffix);

      // Stop reading as soon as the first FUNCT section in the input file is empty
      if (section_line_defs.empty()) break;

      Teuchos::RCP<DRT::INPUT::LineDefinition> first_line_in_section = section_line_defs.front();

      // List all known TryCreate functions in a vector, so they can be called with a unified
      // syntax below. Also, erase their exact return type, since we can only store std::any.
      std::vector<std::function<std::any(
          Teuchos::RCP<DRT::INPUT::LineDefinition>, DRT::UTILS::FunctionManager&, const int)>>
          try_create_function_vector{WrapFunction(DRT::UTILS::TryCreateVariableExprFunction<dim>),
              WrapFunction(POROMULTIPHASESCATRA::TryCreatePoroFunction<dim>),
              WrapFunction(STR::TryCreateStructureFunction),
              WrapFunction(FLD::TryCreateFluidFunction),
              WrapFunction(DRT::UTILS::TryCreateCombustFunction),
              WrapFunction(DRT::UTILS::TryCreateXfluidFunction),
              WrapFunction(DRT::UTILS::TryCreateLibraryFunction),
              WrapFunction(DRT::UTILS::TryCreateLibraryFunctionScalar)};

      bool found_function = false;

      // First, parse functions that read a single line
      for (const auto& try_create_function : try_create_function_vector)
      {
        auto special_funct =
            try_create_function(first_line_in_section, functionManager, funct_suffix);
        if (special_funct.has_value())
        {
          functions.emplace_back(special_funct);
          found_function = true;
          break;
        }
      }

      // List all known multi-line functions in a vector, so they can be called with a unified
      // syntax below. Also, erase their exact return type, since we can only store std::any.
      std::vector<std::function<std::any(std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>>)>>
          try_create_multiline_function_vector{
              WrapMultiLineFunction(DRT::UTILS::TryCreateExprFunction<dim>),
              WrapMultiLineFunction(DRT::UTILS::TryCreateFunctionOfTime)};

      // If we haven't found the function by now, try to parse a multi-line function.
      if (!found_function)
      {
        for (const auto& try_create_multiline_function : try_create_multiline_function_vector)
        {
          auto basic_funct = try_create_multiline_function(section_line_defs);
          if (basic_funct.has_value())
          {
            functions.emplace_back(basic_funct);
            found_function = true;
            break;
          }
        }
      }

      if (!found_function)
      {
        dserror("Could not create any function from the given function line definition.");
      }
    }
  }
}  // namespace

void PrintFunctionDatHeader()
{
  DRT::UTILS::FunctionManager functionmanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = functionmanager.ValidFunctionLines();

  lines->Print(std::cout);
}

void DRT::UTILS::AddValidFunctionFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines)
{
  DRT::INPUT::LineDefinition onecomponentexpr;
  onecomponentexpr.AddNamedString("SYMBOLIC_FUNCTION_OF_SPACE_TIME");

  DRT::INPUT::LineDefinition symbolic_function_of_time;
  symbolic_function_of_time.AddNamedString("SYMBOLIC_FUNCTION_OF_TIME");

  DRT::INPUT::LineDefinition componentexpr;
  componentexpr.AddNamedInt("COMPONENT").AddNamedString("SYMBOLIC_FUNCTION_OF_SPACE_TIME");

  DRT::INPUT::LineDefinition variableexpr;
  variableexpr.AddNamedInt("VARIABLE")
      .AddNamedString("NAME")
      .AddNamedString("TYPE")
      .AddOptionalNamedString("DESCRIPTION")
      .AddOptionalNamedInt("NUMPOINTS")
      .AddOptionalNamedString("BYNUM")
      .AddOptionalNamedDoubleVector("TIMERANGE", 2)
      .AddOptionalNamedDoubleVector("TIMES", "NUMPOINTS")
      .AddOptionalNamedDoubleVector("VALUES", "NUMPOINTS")
      .AddOptionalNamedString("PERIODIC")
      .AddOptionalNamedDouble("T1")
      .AddOptionalNamedDouble("T2");

  DRT::INPUT::LineDefinition variableexprmulti;
  variableexprmulti.AddNamedInt("VARIABLE")
      .AddNamedString("NAME")
      .AddNamedString("TYPE")
      .AddOptionalNamedInt("NUMPOINTS")
      .AddOptionalNamedString("BYNUM")
      .AddOptionalNamedDoubleVector("TIMERANGE", 2)
      .AddOptionalNamedDoubleVector("TIMES", "NUMPOINTS")
      .AddOptionalNamedDoubleVector("VALUES", "NUMPOINTS")
      .AddOptionalNamedStringVector("DESCRIPTION", "NUMPOINTS")  // only NUMPOINTS-1 are taken
      .AddOptionalNamedString("PERIODIC")
      .AddOptionalNamedDouble("T1")
      .AddOptionalNamedDouble("T2");

  lines->Add(onecomponentexpr);
  lines->Add(symbolic_function_of_time);
  lines->Add(componentexpr);
  lines->Add(variableexpr);
  lines->Add(variableexprmulti);

  DRT::INPUT::LineDefinition varfunct;
  varfunct.AddNamedString("VARFUNCTION")
      .AddOptionalNamedInt("NUMCONSTANTS")
      .AddOptionalNamedPairOfStringAndDoubleVector("CONSTANTS", "NUMCONSTANTS");

  lines->Add(varfunct);
}

Teuchos::RCP<DRT::INPUT::Lines> DRT::UTILS::FunctionManager::ValidFunctionLines()
{
  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("FUNCT"));

  DRT::UTILS::AddValidFunctionFunctionLines(lines);
  DRT::UTILS::AddValidLibraryFunctionLines(lines);
  STR::AddValidStructureFunctionLines(lines);
  FLD::AddValidFluidFunctionLines(lines);
  DRT::UTILS::AddValidCombustFunctionLines(lines);
  DRT::UTILS::AddValidXfluidFunctionLines(lines);
  POROMULTIPHASESCATRA::AddValidPoroFunctionLines(lines);

  return lines;
}


void DRT::UTILS::FunctionManager::ReadInput(DRT::INPUT::DatFileReader& reader)
{
  functions_.clear();

  switch (DRT::Problem::Instance()->NDim())
  {
    case 1:
      FillFunctions<1>(reader, functions_, *this);
      break;
    case 2:
      FillFunctions<2>(reader, functions_, *this);
      break;
    case 3:
      FillFunctions<3>(reader, functions_, *this);
      break;
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
  }
}
