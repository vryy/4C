/*----------------------------------------------------------------------*/
/*! \file

\brief Managing of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#include "baci_fluid_functions.H"
#include "baci_fluid_xfluid_functions.H"
#include "baci_fluid_xfluid_functions_combust.H"
#include "baci_io.H"
#include "baci_lib_function.H"
#include "baci_lib_function_library.H"
#include "baci_lib_function_of_time.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_linedefinition.H"
#include "baci_poromultiphase_scatra_function.H"
#include "baci_structure_new_functions.H"

#include <stdexcept>


namespace
{
  using LineDefinitionVector = std::vector<DRT::INPUT::LineDefinition>;

  using TypeErasedFunctionCreator = std::function<std::any(const LineDefinitionVector&)>;

  template <typename T>
  using FunctionCreator = Teuchos::RCP<T> (*)(const LineDefinitionVector&);

  /**
   * Utility function that takes a function object returning a Teuchos::RCP<T> and erases its return
   * type via std::any. In addition, if the returned object would be Teuchos::null, discard it and
   * return an empty std::any instead.
   */
  template <typename T>
  TypeErasedFunctionCreator WrapFunction(FunctionCreator<T> fun)
  {
    return [fun](const LineDefinitionVector& linedefs) -> std::any
    {
      Teuchos::RCP<T> created = fun(linedefs);
      if (created == Teuchos::null)
        return {};
      else
        return created;
    };
  }

  template <int dim>
  void FillFunctions(DRT::INPUT::DatFileReader& reader, std::vector<std::any>& functions)
  {
    DRT::INPUT::Lines lines = DRT::UTILS::FunctionManager::ValidFunctionLines();

    // Read FUNCT sections starting from FUNCT1 until the first empty one is encountered.
    // This implies that the FUNCT sections must form a contiguous range in the input file.
    // Otherwise, the read fails later.
    for (int funct_suffix = 1;; ++funct_suffix)
    {
      // Read lines belonging to section FUNCT<i> in the input file
      std::vector<DRT::INPUT::LineDefinition> section_line_defs = lines.Read(reader, funct_suffix);

      // Stop reading as soon as the first FUNCT section in the input file is empty
      if (section_line_defs.empty()) break;

      // List all known TryCreate functions in a vector, so they can be called with a unified
      // syntax below. Also, erase their exact return type, since we can only store std::any.
      std::vector<TypeErasedFunctionCreator> try_create_function_vector{
          WrapFunction(DRT::UTILS::TryCreateSymbolicFunctionOfAnything<dim>),
          WrapFunction(POROMULTIPHASESCATRA::TryCreatePoroFunction<dim>),
          WrapFunction(STR::TryCreateStructureFunction), WrapFunction(FLD::TryCreateFluidFunction),
          WrapFunction(DRT::UTILS::TryCreateCombustFunction),
          WrapFunction(DRT::UTILS::TryCreateXfluidFunction),
          WrapFunction(DRT::UTILS::TryCreateLibraryFunctionScalar),
          WrapFunction(DRT::UTILS::TryCreateSymbolicFunctionOfSpaceTime<dim>),
          WrapFunction(DRT::UTILS::TryCreateFunctionOfTime)};

      const bool found_and_inserted_function = std::invoke(
          [&]()
          {
            for (const auto& try_create_function : try_create_function_vector)
            {
              auto maybe_function = try_create_function(section_line_defs);
              if (maybe_function.has_value())
              {
                functions.emplace_back(maybe_function);
                return true;
              }
            }
            return false;
          });

      // If we didn't find the function, we don't know how to read the lines we were given.
      if (not found_and_inserted_function)
      {
        std::stringstream ss;
        for (const auto& line : section_line_defs)
        {
          ss << "\n";
          line.Print(ss);
        }

        dserror("Could not create any function from the following function line definition:\n%s",
            ss.str().c_str());
      }
    }
  }
}  // namespace

void PrintFunctionDatHeader()
{
  DRT::UTILS::FunctionManager functionmanager;
  DRT::INPUT::Lines lines = functionmanager.ValidFunctionLines();

  lines.Print(std::cout);
}

void DRT::UTILS::AddValidFunctionFunctionLines(DRT::INPUT::Lines& lines)
{
  using namespace DRT::INPUT;

  LineDefinition onecomponentexpr =
      LineDefinition::Builder().AddNamedString("SYMBOLIC_FUNCTION_OF_SPACE_TIME").Build();

  LineDefinition symbolic_function_of_time =
      LineDefinition::Builder().AddNamedString("SYMBOLIC_FUNCTION_OF_TIME").Build();

  LineDefinition componentexpr = LineDefinition::Builder()
                                     .AddNamedInt("COMPONENT")
                                     .AddNamedString("SYMBOLIC_FUNCTION_OF_SPACE_TIME")
                                     .Build();

  LineDefinition variableexprmulti =
      LineDefinition::Builder()
          .AddNamedInt("VARIABLE")
          .AddNamedString("NAME")
          .AddNamedString("TYPE")
          .AddOptionalNamedInt("NUMPOINTS")
          .AddOptionalNamedString("BYNUM")
          .AddOptionalNamedDoubleVector("TIMERANGE", 2)
          .AddOptionalNamedDoubleVector("TIMES", LengthFromIntNamed("NUMPOINTS"))
          .AddOptionalNamedDoubleVector("VALUES", LengthFromIntNamed("NUMPOINTS"))
          .AddOptionalNamedStringVector("DESCRIPTION",
              // Special case where only NUMPOINTS-1 are taken
              [](const LineDefinition& already_read_line)
              {
                try
                {
                  int length;
                  already_read_line.ExtractInt("NUMPOINTS", length);
                  return length - 1;
                }
                catch (const std::runtime_error& e)
                {
                  // When NUMPOINTS is not set, then we still allow for a single DESCRIPTION entry
                  return 1;
                }
              })
          .AddOptionalNamedString("PERIODIC")
          .AddOptionalNamedDouble("T1")
          .AddOptionalNamedDouble("T2")
          .Build();

  lines.Add(onecomponentexpr);
  lines.Add(symbolic_function_of_time);
  lines.Add(componentexpr);
  lines.Add(variableexprmulti);

  LineDefinition varfunct = LineDefinition::Builder()
                                .AddNamedString("VARFUNCTION")
                                .AddOptionalNamedInt("NUMCONSTANTS")
                                .AddOptionalNamedPairOfStringAndDoubleVector(
                                    "CONSTANTS", LengthFromIntNamed("NUMCONSTANTS"))
                                .Build();

  lines.Add(varfunct);
}

DRT::INPUT::Lines DRT::UTILS::FunctionManager::ValidFunctionLines()
{
  DRT::INPUT::Lines lines(
      "FUNCT", "Definition of functions for various cases, mainly boundary conditions");

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
      FillFunctions<1>(reader, functions_);
      break;
    case 2:
      FillFunctions<2>(reader, functions_);
      break;
    case 3:
      FillFunctions<3>(reader, functions_);
      break;
    default:
      dserror("Unsupported dimension %d.", DRT::Problem::Instance()->NDim());
  }
}
