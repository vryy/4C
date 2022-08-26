/*----------------------------------------------------------------------*/
/*! \file

\brief Managing of space- and/or time-dependent functions

\level 0

*/
/*----------------------------------------------------------------------*/

#include "drt_function.H"
#include "drt_linedefinition.H"
#include "../drt_fluid/fluid_functions.H"
#include "../drt_fluid_xfluid/xfluid_functions.H"
#include "../drt_fluid_xfluid/xfluid_functions_combust.H"
#include "../drt_structure_new/str_functions.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_function.H"
#include "drt_function_library.H"
#include "../drt_io/io.H"
#include "drt_globalproblem.H"

void PrintFunctionDatHeader()
{
  DRT::UTILS::FunctionManager functionmanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = functionmanager.ValidFunctionLines();

  lines->Print(std::cout);
}

void DRT::UTILS::AddValidFunctionFunctionLines(Teuchos::RCP<DRT::INPUT::Lines> lines)
{
  DRT::INPUT::LineDefinition onecomponentexpr;
  onecomponentexpr.AddNamedString("FUNCTION");

  DRT::INPUT::LineDefinition componentexpr;
  componentexpr.AddNamedInt("COMPONENT").AddNamedString("FUNCTION");

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

template <int dim>
void FillFunctions(DRT::INPUT::DatFileReader& reader,
    std::vector<Teuchos::RCP<DRT::UTILS::TemporaryFunctionInterface>>& functions,
    DRT::UTILS::FunctionManager& functionManager)
{
  Teuchos::RCP<DRT::INPUT::Lines> lines = DRT::UTILS::FunctionManager::ValidFunctionLines();

  // test for as many functions as there are
  for (int i = 1;; ++i)
  {
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions_lin_defs =
        lines->Read(reader, i);

    if (functions_lin_defs.empty())
    {
      break;
    }
    else
    {
      Teuchos::RCP<DRT::INPUT::LineDefinition> function_lin_def = functions_lin_defs[0];

      bool found_function = false;

      // list all known TryCreate functions in a vector so they can be called with a unified syntax
      // below
      std::vector<std::function<Teuchos::RCP<DRT::UTILS::Function>(
          Teuchos::RCP<DRT::INPUT::LineDefinition>, DRT::UTILS::FunctionManager&, const int)>>
          try_create_function_vector{DRT::UTILS::TryCreateVariableExprFunction<dim>,
              POROMULTIPHASESCATRA::TryCreatePoroFunction<dim>, STR::TryCreateStructureFunction,
              FLD::TryCreateFluidFunction, DRT::UTILS::TryCreateCombustFunction,
              DRT::UTILS::TryCreateXfluidFunction, DRT::UTILS::TryCreateLibraryFunction};

      for (const auto& try_create_function : try_create_function_vector)
      {
        auto special_funct = try_create_function(function_lin_def, functionManager, i);
        if (special_funct != Teuchos::null)
        {
          functions.emplace_back(special_funct);
          found_function = true;
          break;  // jumps out of for statement
        }
      }

      if (!found_function)
      {
        auto basic_funct = DRT::UTILS::TryCreateExprFunction<dim>(functions_lin_defs);
        if (basic_funct != Teuchos::null)
        {
          functions.emplace_back(basic_funct);
        }
        else
        {
          dserror("Could not create any function from the given function line definition.");
        }
      }
    }
  }
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

DRT::UTILS::Function& DRT::UTILS::FunctionManager::Funct(int num)
{
  // ensure that desired function is available (prevents segmentation fault)
  if (functions_.size() < (unsigned int)(num + 1) || num < 0)
    dserror("function %d not available", num + 1);

  return dynamic_cast<Function&>(*(functions_[num]));
}