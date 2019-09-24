/*----------------------------------------------------------------------*/
/*! \file

\brief Managing and evaluating of space- and/or time-dependent functions

\level 0

\maintainer Martin Kronbichler

*/
/*----------------------------------------------------------------------*/

#include <Sacado.hpp>
#include "drt_function.H"
#include "drt_functionvariables.H"
#include "drt_globalproblem.H"
#include "standardtypes_cpp.H"
#include "drt_linedefinition.H"
#include "../drt_fluid/fluid_functions.H"
#include "../drt_structure_new/str_functions.H"
#include "../drt_poromultiphase_scatra/poromultiphase_scatra_function.H"
#include "../drt_fluid_xfluid/xfluid_functions.H"
#include "../drt_fluid_xfluid/xfluid_functions_combust.H"
#include "../drt_io/io.H"
#include "drt_function_library.H"

/*----------------------------------------------------------------------*/
//! Print function
/*----------------------------------------------------------------------*/
void PrintFunctionDatHeader()
{
  DRT::UTILS::FunctionManager functionmanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = functionmanager.ValidFunctionLines();

  lines->Print(std::cout);
}


Teuchos::RCP<DRT::INPUT::Lines> DRT::UTILS::FunctionManager::ValidFunctionLines()
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

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("FUNCT"));
  lines->Add(onecomponentexpr);
  lines->Add(componentexpr);
  lines->Add(variableexpr);
  lines->Add(variableexprmulti);

  DRT::INPUT::LineDefinition varfunct;
  varfunct.AddNamedString("VARFUNCTION")
      .AddOptionalNamedInt("NUMCONSTANTS")
      .AddOptionalNamedPairOfStringAndDoubleVector("CONSTANTS", "NUMCONSTANTS");

  DRT::INPUT::LineDefinition linelin;
  linelin.AddNamedDoubleVector("LINE_LIN", 8);

  DRT::INPUT::LineDefinition radiuslin;
  radiuslin.AddNamedDoubleVector("RADIUS_LIN", 8);

  DRT::INPUT::LineDefinition radiusquad;
  radiusquad.AddNamedDoubleVector("RADIUS_QUAD", 6);

  DRT::INPUT::LineDefinition beltrami;
  beltrami.AddTag("BELTRAMI").AddNamedDouble("c1");

  DRT::INPUT::LineDefinition channelweaklycompressible;
  channelweaklycompressible.AddTag("CHANNELWEAKLYCOMPRESSIBLE");

  DRT::INPUT::LineDefinition correctiontermchannelweaklycompressible;
  correctiontermchannelweaklycompressible.AddTag("CORRECTIONTERMCHANNELWEAKLYCOMPRESSIBLE");

  DRT::INPUT::LineDefinition channelweaklycompressiblefourier3;
  channelweaklycompressiblefourier3.AddTag("CHANNELWEAKLYCOMPRESSIBLEFOURIER3");

  DRT::INPUT::LineDefinition correctiontermchannelweaklycompressiblefourier3;
  correctiontermchannelweaklycompressiblefourier3.AddTag(
      "CORRECTIONTERMCHANNELWEAKLYCOMPRESSIBLEFOURIER3");

  DRT::INPUT::LineDefinition bodyforcechannelweaklycompressiblefourier3;
  bodyforcechannelweaklycompressiblefourier3.AddTag("BODYFORCECHANNELWEAKLYCOMPRESSIBLEFOURIER3");

  DRT::INPUT::LineDefinition weaklycompressiblepoiseuille;
  weaklycompressiblepoiseuille.AddTag("WEAKLYCOMPRESSIBLE_POISEUILLE")
      .AddNamedInt("MAT")
      .AddNamedDouble("L")
      .AddNamedDouble("R")
      .AddNamedDouble("U");

  DRT::INPUT::LineDefinition weaklycompressiblepoiseuilleforce;
  weaklycompressiblepoiseuilleforce.AddTag("WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE")
      .AddNamedInt("MAT")
      .AddNamedDouble("L")
      .AddNamedDouble("R")
      .AddNamedDouble("U");

  DRT::INPUT::LineDefinition weaklycompressiblemanufacturedflow;
  weaklycompressiblemanufacturedflow.AddTag("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressiblemanufacturedflowforce;
  weaklycompressiblemanufacturedflowforce.AddTag("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW_FORCE")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennecfd;
  weaklycompressibleetiennecfd.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_CFD").AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennecfdforce;
  weaklycompressibleetiennecfdforce.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_FORCE")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennecfdviscosity;
  weaklycompressibleetiennecfdviscosity.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_VISCOSITY")
      .AddNamedInt("MAT");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsifluid;
  weaklycompressibleetiennefsifluid.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID")
      .AddNamedInt("MAT_FLUID")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsifluidforce;
  weaklycompressibleetiennefsifluidforce.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE")
      .AddNamedInt("MAT_FLUID")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsifluidviscosity;
  weaklycompressibleetiennefsifluidviscosity
      .AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY")
      .AddNamedInt("MAT_FLUID")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsistructure;
  weaklycompressibleetiennefsistructure.AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition weaklycompressibleetiennefsistructureforce;
  weaklycompressibleetiennefsistructureforce
      .AddTag("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE_FORCE")
      .AddNamedInt("MAT_STRUC");

  DRT::INPUT::LineDefinition kimmoin;
  kimmoin.AddTag("KIM-MOIN");

  DRT::INPUT::LineDefinition bochevup;
  bochevup.AddTag("BOCHEV-UP");

  DRT::INPUT::LineDefinition bochevrhs;
  bochevrhs.AddTag("BOCHEV-RHS");

  DRT::INPUT::LineDefinition beltramiup;
  beltramiup.AddTag("BELTRAMI-UP").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition beltramigradu;
  beltramigradu.AddTag("BELTRAMI-GRADU").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition beltramirhs;
  beltramirhs.AddTag("BELTRAMI-RHS")
      .AddNamedInt("MAT")
      .AddNamedInt("ISSTAT")
      .AddNamedInt("ISSTOKES");

  DRT::INPUT::LineDefinition kimmoinup;
  kimmoinup.AddTag("KIMMOIN-UP").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition kimmoingradu;
  kimmoingradu.AddTag("KIMMOIN-GRADU").AddNamedInt("MAT").AddNamedInt("ISSTAT");

  DRT::INPUT::LineDefinition kimmoinrhs;
  kimmoinrhs.AddTag("KIMMOIN-RHS").AddNamedInt("MAT").AddNamedInt("ISSTAT").AddNamedInt("ISSTOKES");

  DRT::INPUT::LineDefinition kimmoinstress;
  kimmoinstress.AddTag("KIMMOIN-STRESS")
      .AddNamedInt("MAT")
      .AddNamedInt("ISSTAT")
      .AddNamedDouble("AMPLITUDE");

  DRT::INPUT::LineDefinition turbboulayer;
  turbboulayer.AddTag("TURBBOULAYER");

  DRT::INPUT::LineDefinition turbboulayerbfs;
  turbboulayerbfs.AddTag("TURBBOULAYER-BFS");

  DRT::INPUT::LineDefinition turbboulayeroracles;
  turbboulayeroracles.AddTag("TURBBOULAYERORACLES");

  DRT::INPUT::LineDefinition jefferyhamel;
  jefferyhamel.AddTag("JEFFERY-HAMEL");

  DRT::INPUT::LineDefinition womersley;
  womersley.AddTag("WOMERSLEY")
      .AddNamedInt("Local")
      .AddNamedInt("MAT")
      .AddNamedInt("CURVE")
      .AddNamedString("FSI");

  DRT::INPUT::LineDefinition localwomersley;
  localwomersley.AddTag("WOMERSLEY")
      .AddNamedInt("Local")
      .AddNamedDouble("Radius")
      .AddNamedInt("MAT")
      .AddNamedInt("CURVE");

  DRT::INPUT::LineDefinition cylinder3d;
  cylinder3d.AddNamedDouble("CYLINDER_3D");

  DRT::INPUT::LineDefinition controlledrotation;
  controlledrotation.AddTag("CONTROLLEDROTATION")
      .AddNamedString("FILE")
      .AddNamedString("TYPE")
      .AddNamedDoubleVector("ORIGIN", 3);

  DRT::INPUT::LineDefinition accelerationprofile;
  accelerationprofile.AddTag("ACCELERATIONPROFILE").AddNamedString("FILE");

  DRT::INPUT::LineDefinition ramptovalue;
  ramptovalue.AddTag("RAMPTOVALUE")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("STARTTIME")
      .AddNamedDouble("DURATION")
      .AddNamedString("TYPE");

  DRT::INPUT::LineDefinition nodenormal;
  nodenormal.AddTag("NODENORMAL")
      .AddNamedString("GEOMETRY")
      .AddNamedDoubleVector("ORIGIN", 3)
      .AddNamedDouble("RADIUS")
      .AddNamedDouble("CYLINDERHEIGHT")
      .AddNamedDoubleVector("ORIENTATION", 3)
      .AddNamedDouble("CASSINIA");

  DRT::INPUT::LineDefinition poromultiphasescatra_funct;
  poromultiphasescatra_funct.AddNamedString("POROMULTIPHASESCATRA_FUNCTION")
      .AddOptionalNamedInt("NUMPARAMS")
      .AddOptionalNamedPairOfStringAndDoubleVector("PARAMS", "NUMPARAMS");

  DRT::INPUT::LineDefinition fastpolynomial_funct;
  fastpolynomial_funct.AddTag("FASTPOLYNOMIAL")
      .AddNamedInt("NUMCOEFF")
      .AddNamedDoubleVector("COEFF", "NUMCOEFF");

  DRT::INPUT::LineDefinition translatedfunction_funct;
  translatedfunction_funct.AddTag("TRANSLATEDFUNCTION").AddNamedInt("ORIGIN").AddNamedInt("LOCAL");

  lines->Add(translatedfunction_funct);

  lines->Add(varfunct);
  lines->Add(linelin);
  lines->Add(radiuslin);
  lines->Add(radiusquad);
  lines->Add(beltrami);
  lines->Add(channelweaklycompressible);
  lines->Add(correctiontermchannelweaklycompressible);
  lines->Add(channelweaklycompressiblefourier3);
  lines->Add(correctiontermchannelweaklycompressiblefourier3);
  lines->Add(bodyforcechannelweaklycompressiblefourier3);
  lines->Add(weaklycompressiblepoiseuille);
  lines->Add(weaklycompressiblepoiseuilleforce);
  lines->Add(weaklycompressiblemanufacturedflow);
  lines->Add(weaklycompressiblemanufacturedflowforce);
  lines->Add(weaklycompressibleetiennecfd);
  lines->Add(weaklycompressibleetiennecfdforce);
  lines->Add(weaklycompressibleetiennecfdviscosity);
  lines->Add(weaklycompressibleetiennefsifluid);
  lines->Add(weaklycompressibleetiennefsifluidforce);
  lines->Add(weaklycompressibleetiennefsifluidviscosity);
  lines->Add(weaklycompressibleetiennefsistructure);
  lines->Add(weaklycompressibleetiennefsistructureforce);
  lines->Add(kimmoin);
  lines->Add(bochevup);
  lines->Add(bochevrhs);
  lines->Add(beltramiup);
  lines->Add(beltramigradu);
  lines->Add(beltramirhs);
  lines->Add(kimmoinup);
  lines->Add(kimmoingradu);
  lines->Add(kimmoinrhs);
  lines->Add(kimmoinstress);
  lines->Add(turbboulayer);
  lines->Add(turbboulayerbfs);
  lines->Add(turbboulayeroracles);
  lines->Add(jefferyhamel);
  lines->Add(womersley);
  lines->Add(localwomersley);
  lines->Add(cylinder3d);
  lines->Add(controlledrotation);
  lines->Add(accelerationprofile);
  lines->Add(ramptovalue);
  lines->Add(nodenormal);
  lines->Add(poromultiphasescatra_funct);
  lines->Add(fastpolynomial_funct);
  // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  DRT::UTILS::CombustValidFunctionLines(lines);
  DRT::UTILS::XfluidValidFunctionLines(lines);

  return lines;
}


void DRT::UTILS::FunctionManager::ReadInput(DRT::INPUT::DatFileReader& reader)
{
  functions_.clear();

  Teuchos::RCP<DRT::INPUT::Lines> lines = ValidFunctionLines();

  // test for as many functions as there are
  for (int i = 1;; ++i)
  {
    std::vector<Teuchos::RCP<DRT::INPUT::LineDefinition>> functions = lines->Read(reader, i);

    if (functions.size() == 0)
      break;

    else
    {
      Teuchos::RCP<DRT::INPUT::LineDefinition> function = functions[0];

      // Old function +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (function->HaveNamed("BELTRAMI"))
      {
        double c1;
        function->ExtractDouble("c1", c1);

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiFunction(c1)));
      }
      else if (function->HaveNamed("CHANNELWEAKLYCOMPRESSIBLE"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::ChannelWeaklyCompressibleFunction()));
      }
      else if (function->HaveNamed("CORRECTIONTERMCHANNELWEAKLYCOMPRESSIBLE"))
      {
        functions_.push_back(
            Teuchos::rcp(new FLD::CorrectionTermChannelWeaklyCompressibleFunction()));
      }
      else if (function->HaveNamed("CHANNELWEAKLYCOMPRESSIBLEFOURIER3"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::ChannelWeaklyCompressibleFourier3Function()));
      }
      else if (function->HaveNamed("CORRECTIONTERMCHANNELWEAKLYCOMPRESSIBLEFOURIER3"))
      {
        functions_.push_back(
            Teuchos::rcp(new FLD::CorrectionTermChannelWeaklyCompressibleFourier3Function()));
      }
      else if (function->HaveNamed("BODYFORCECHANNELWEAKLYCOMPRESSIBLEFOURIER3"))
      {
        functions_.push_back(
            Teuchos::rcp(new FLD::BodyForceChannelWeaklyCompressibleFourier3Function()));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_POISEUILLE"))
      {
        // read data
        int mat_id = -1;
        double L = 0.0;
        double R = 0.0;
        double U = 0.0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractDouble("L", L);
        function->ExtractDouble("R", R);
        function->ExtractDouble("U", U);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_POISEUILLE");
        if (L <= 0) dserror("Please give a (reasonable) 'L' in WEAKLYCOMPRESSIBLE_POISEUILLE");
        if (R <= 0) dserror("Please give a (reasonable) 'R' in WEAKLYCOMPRESSIBLE_POISEUILLE");
        if (U <= 0) dserror("Please give a (reasonable) 'U' in WEAKLYCOMPRESSIBLE_POISEUILLE");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressiblePoiseuilleFunction(mat_id, L, R, U)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE"))
      {
        // read data
        int mat_id = -1;
        double L = 0.0;
        double R = 0.0;
        double U = 0.0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractDouble("L", L);
        function->ExtractDouble("R", R);
        function->ExtractDouble("U", U);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");
        if (L <= 0)
          dserror("Please give a (reasonable) 'L' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");
        if (R <= 0)
          dserror("Please give a (reasonable) 'R' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");
        if (U <= 0)
          dserror("Please give a (reasonable) 'U' in WEAKLYCOMPRESSIBLE_POISEUILLE_FORCE");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressiblePoiseuilleForceFunction(mat_id, L, R, U)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW"))
      {
        // read data
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressibleManufacturedFlowFunction(mat_id)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW_FORCE"))
      {
        // read data
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_MANUFACTUREDFLOW_FORCE");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressibleManufacturedFlowForceFunction(mat_id)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_CFD"))
      {
        // read data
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_ETIENNE_CFD");

        functions_.push_back(Teuchos::rcp(new FLD::WeaklyCompressibleEtienneCFDFunction(mat_id)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_FORCE"))
      {
        // read data
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_ETIENNE_CFD_FORCE");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressibleEtienneCFDForceFunction(mat_id)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_CFD_VISCOSITY"))
      {
        // read data
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0)
          dserror("Please give a (reasonable) 'MAT' in WEAKLYCOMPRESSIBLE_ETIENNE_CFD_VISCOSITY");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressibleEtienneCFDViscosityFunction(mat_id)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID"))
      {
        // read data
        int mat_id_fluid = -1;
        int mat_id_struc = -1;

        function->ExtractInt("MAT_FLUID", mat_id_fluid);
        function->ExtractInt("MAT_STRUC", mat_id_struc);

        if (mat_id_fluid <= 0)
          dserror("Please give a (reasonable) 'MAT_FLUID' in WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID");
        if (mat_id_struc <= 0)
          dserror("Please give a (reasonable) 'MAT_STRUC' in WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID");

        functions_.push_back(Teuchos::rcp(
            new FLD::WeaklyCompressibleEtienneFSIFluidFunction(mat_id_fluid, mat_id_struc)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE"))
      {
        // read data
        int mat_id_fluid = -1;
        int mat_id_struc = -1;

        function->ExtractInt("MAT_FLUID", mat_id_fluid);
        function->ExtractInt("MAT_STRUC", mat_id_struc);

        if (mat_id_fluid <= 0)
          dserror(
              "Please give a (reasonable) 'MAT_FLUID' in "
              "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE");
        if (mat_id_struc <= 0)
          dserror(
              "Please give a (reasonable) 'MAT_STRUC' in "
              "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_FORCE");

        functions_.push_back(Teuchos::rcp(
            new FLD::WeaklyCompressibleEtienneFSIFluidForceFunction(mat_id_fluid, mat_id_struc)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY"))
      {
        // read data
        int mat_id_fluid = -1;
        int mat_id_struc = -1;

        function->ExtractInt("MAT_FLUID", mat_id_fluid);
        function->ExtractInt("MAT_STRUC", mat_id_struc);

        if (mat_id_fluid <= 0)
          dserror(
              "Please give a (reasonable) 'MAT_FLUID' in "
              "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY");
        if (mat_id_struc <= 0)
          dserror(
              "Please give a (reasonable) 'MAT_STRUC' in "
              "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_FLUID_VISCOSITY");

        functions_.push_back(
            Teuchos::rcp(new FLD::WeaklyCompressibleEtienneFSIFluidViscosityFunction(
                mat_id_fluid, mat_id_struc)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE"))
      {
        // read data
        int mat_id_struc = -1;

        function->ExtractInt("MAT_STRUC", mat_id_struc);

        if (mat_id_struc <= 0)
          dserror(
              "Please give a (reasonable) 'MAT_STRUC' in WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE");

        functions_.push_back(
            Teuchos::rcp(new STR::WeaklyCompressibleEtienneFSIStructureFunction(mat_id_struc)));
      }
      else if (function->HaveNamed("WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE_FORCE"))
      {
        // read data
        int mat_id_struc = -1;

        function->ExtractInt("MAT_STRUC", mat_id_struc);

        if (mat_id_struc <= 0)
          dserror(
              "Please give a (reasonable) 'MAT_STRUC' in "
              "WEAKLYCOMPRESSIBLE_ETIENNE_FSI_STRUCTURE_FORCE");

        functions_.push_back(Teuchos::rcp(
            new STR::WeaklyCompressibleEtienneFSIStructureForceFunction(mat_id_struc)));
      }
      else if (function->HaveNamed("KIM-MOIN"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::KimMoinFunction()));
      }
      else if (function->HaveNamed("BOCHEV-UP"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::BochevUPFunction()));
      }
      else if (function->HaveNamed("BOCHEV-RHS"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::BochevRHSFunction()));
      }
      else if (function->HaveNamed("BELTRAMI-UP"))
      {
        // read material
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-UP");

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiUP(mat_id)));
      }
      else if (function->HaveNamed("BELTRAMI-GRADU"))
      {
        // read material
        int mat_id = -1;

        function->ExtractInt("MAT", mat_id);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiGradU(mat_id)));
      }
      else if (function->HaveNamed("BELTRAMI-RHS"))
      {
        // read material
        int mat_id = -1;
        int is_stokes = 0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractInt("ISSTOKES", is_stokes);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in BELTRAMI-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::BeltramiRHS(mat_id, (bool)is_stokes)));
      }
      else if (function->HaveNamed("KIMMOIN-UP"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractInt("ISSTAT", is_stationary);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-UP");

        functions_.push_back(Teuchos::rcp(new FLD::KimMoinUP(mat_id, (bool)is_stationary)));
      }
      else if (function->HaveNamed("KIMMOIN-GRADU"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractInt("ISSTAT", is_stationary);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-GRADU");

        functions_.push_back(Teuchos::rcp(new FLD::KimMoinGradU(mat_id, (bool)is_stationary)));
      }
      else if (function->HaveNamed("KIMMOIN-RHS"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;
        int is_stokes = 0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractInt("ISSTAT", is_stationary);
        function->ExtractInt("ISSTOKES", is_stokes);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-GRADU");

        functions_.push_back(
            Teuchos::rcp(new FLD::KimMoinRHS(mat_id, (bool)is_stationary, (bool)is_stokes)));
      }
      else if (function->HaveNamed("KIMMOIN-STRESS"))
      {
        // read material
        int mat_id = -1;
        int is_stationary = 0;
        double amplitude = 1.0;

        function->ExtractInt("MAT", mat_id);
        function->ExtractInt("ISSTAT", is_stationary);
        function->ExtractDouble("AMPLITUDE", amplitude);

        if (mat_id <= 0) dserror("Please give a (reasonable) 'MAT'/material in KIMMOIN-STRESS");

        functions_.push_back(
            Teuchos::rcp(new FLD::KimMoinStress(mat_id, (bool)is_stationary, amplitude)));
      }
      else if (function->HaveNamed("TURBBOULAYER"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::TurbBouLayerFunction()));
      }
      else if (function->HaveNamed("TURBBOULAYER-BFS"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::TurbBouLayerFunctionBFS()));
      }
      else if (function->HaveNamed("TURBBOULAYERORACLES"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::TurbBouLayerFunctionORACLES()));
      }
      else if (function->HaveNamed("JEFFERY-HAMEL"))
      {
        functions_.push_back(Teuchos::rcp(new FLD::JefferyHamelFlowFunction()));
      }
      else if (function->HaveNamed("ZALESAKSDISK"))
      {
        functions_.push_back(Teuchos::rcp(new ZalesaksDiskFunction()));
      }
      else if (function->HaveNamed("CIRCULARFLAME2"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::CircularFlame2Function()));
      }
      else if (function->HaveNamed("CIRCULARFLAME3"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::CircularFlame3Function()));
      }
      else if (function->HaveNamed("CIRCULARFLAME4"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::CircularFlame4Function()));
      }
      else if (function->HaveNamed("DAMBREAKOBSTACLE"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::DamBreakObstacle()));
      }
      else if (function->HaveNamed("COLLAPSINGWATERCOLUMN"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::CollapsingWaterColumnFunction()));
      }
      else if (function->HaveNamed("COLLAPSINGWATERCOLUMNCOARSE"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::CollapsingWaterColumnFunctionCoarse()));
      }
      else if (function->HaveNamed("IMPACTDROP"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::ImpactFunction()));
      }
      else if (function->HaveNamed("BUBBLES"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::BubbleFunction()));
      }
      else if (function->HaveNamed("ORACLESGFUNC"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::ORACLESGFunction()));
      }
      else if (function->HaveNamed("ROTATINGCONE"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::RotatingConeFunction()));
      }
      else if (function->HaveNamed("LEVELSETCUTTEST"))
      {
        functions_.push_back(Teuchos::rcp(new DRT::UTILS::LevelSetCutTestFunction()));
      }
      else if (function->HaveNamed("FORWARDFACINGSTEP"))
      {
        functions_.push_back(Teuchos::rcp(new GerstenbergerForwardfacingStep()));
      }
      else if (function->HaveNamed("SLIPLENGTHFUNCTION"))
      {
        functions_.push_back(Teuchos::rcp(new SlipLengthLevelSetManipulator()));
      }
      else if (function->HaveNamed("MOVINGLEVELSETCYLINDER"))
      {
        std::vector<double> origin;
        function->ExtractDoubleVector("ORIGIN", origin);

        double radius;
        function->ExtractDouble("RADIUS", radius);

        std::vector<double> direction;
        function->ExtractDoubleVector("DIRECTION", direction);

        double distance;
        function->ExtractDouble("DISTANCE", distance);

        double maxspeed;
        function->ExtractDouble("MAXSPEED", maxspeed);

        functions_.push_back(Teuchos::rcp(
            new MovingLevelSetCylinder(&origin, radius, &direction, distance, maxspeed)));
      }
      else if (function->HaveNamed("TAYLORCOUETTEFLOW"))
      {
        double radius_i;
        function->ExtractDouble("RADIUS_I", radius_i);
        double radius_o;
        function->ExtractDouble("RADIUS_O", radius_o);

        double vel_theta_i;
        function->ExtractDouble("VEL_THETA_I", vel_theta_i);
        double vel_theta_o;
        function->ExtractDouble("VEL_THETA_O", vel_theta_o);

        double sliplength_i;
        function->ExtractDouble("SLIPLENGTH_I", sliplength_i);
        double sliplength_o;
        function->ExtractDouble("SLIPLENGTH_O", sliplength_o);

        double traction_theta_i;
        function->ExtractDouble("TRACTION_THETA_I", traction_theta_i);
        double traction_theta_o;
        function->ExtractDouble("TRACTION_THETA_O", traction_theta_o);

        double viscosity;
        function->ExtractDouble("VISCOSITY", viscosity);

        functions_.push_back(
            Teuchos::rcp(new TaylorCouetteFlow(radius_i, radius_o, vel_theta_i, vel_theta_o,
                sliplength_i, sliplength_o, traction_theta_i, traction_theta_o, viscosity)));
      }
      else if (function->HaveNamed("CONTROLLEDROTATION"))
      {
        std::string fileName;
        function->ExtractString("FILE", fileName);

        std::string type;
        function->ExtractString("TYPE", type);

        std::vector<double> origin;
        function->ExtractDoubleVector("ORIGIN", origin);

        functions_.push_back(Teuchos::rcp(
            new ControlledRotationFunction(fileName, type, origin[0], origin[1], origin[2])));
      }
      else if (function->HaveNamed("ACCELERATIONPROFILE"))
      {
        std::string fileName;
        function->ExtractString("FILE", fileName);

        functions_.push_back(Teuchos::rcp(new AccelerationProfileFunction(fileName)));
      }
      else if (function->HaveNamed("FASTPOLYNOMIAL"))
      {
        std::vector<double> coefficients;
        function->ExtractDoubleVector("COEFF", coefficients);

        functions_.push_back(Teuchos::rcp(new FastPolynomialFunction(&coefficients)));
      }
      else if (function->HaveNamed("TRANSLATEDFUNCTION"))
      {
        int origin, local;
        function->ExtractInt("ORIGIN", origin);
        function->ExtractInt("LOCAL", local);

        if (origin <= 0 or origin >= i)
          dserror(
              "ORIGIN function ID (currently %d) must be positive and smaller than "
              "TRANSLATEDFUNCTION (currently %d).",
              origin, i);
        if (local <= 0 or local >= i)
          dserror(
              "LOCAL function ID (currently %d) must be positive and smaller than "
              "TRANSLATEDFUNCTION (currently %d).",
              local, i);

        Teuchos::RCP<Function> origin_funct = Teuchos::rcpFromRef(Funct(origin - 1));
        Teuchos::RCP<Function> local_funct = Teuchos::rcpFromRef(Funct(local - 1));

        functions_.push_back(Teuchos::rcp(new TranslatedFunction(origin_funct, local_funct)));
      }
      else if (function->HaveNamed("VARFUNCTION"))
      {
        Teuchos::RCP<VariableExprFunction> vecfunc = Teuchos::rcp(new VariableExprFunction());

        std::string component;
        function->ExtractString("VARFUNCTION", component);

        std::vector<std::pair<std::string, double>> constants;
        if (function->HaveNamed("CONSTANTS"))
          function->ExtractPairOfStringAndDoubleVector("CONSTANTS", constants);

        vecfunc->AddExpr(component, constants);
        functions_.push_back(vecfunc);
      }
      else if (function->HaveNamed("POROMULTIPHASESCATRA_FUNCTION"))
      {
        std::string type;
        function->ExtractString("POROMULTIPHASESCATRA_FUNCTION", type);

        std::vector<std::pair<std::string, double>> params;
        if (function->HaveNamed("PARAMS"))
          function->ExtractPairOfStringAndDoubleVector("PARAMS", params);

        Teuchos::RCP<POROMULTIPHASESCATRA::PoroMultiPhaseScaTraFunction> vecfunc = Teuchos::null;
        if (type == "TUMOR_GROWTH_LAW_HEAVISIDE")
          vecfunc = Teuchos::rcp(new POROMULTIPHASESCATRA::TumorGrowthLawHeaviside(params));
        else if (type == "NECROSIS_LAW_HEAVISIDE")
          vecfunc = Teuchos::rcp(new POROMULTIPHASESCATRA::NecrosisLawHeaviside(params));
        else if (type == "OXYGEN_CONSUMPTION_LAW_HEAVISIDE")
          vecfunc = Teuchos::rcp(new POROMULTIPHASESCATRA::OxygenConsumptionLawHeaviside(params));
        else if (type == "TUMOR_GROWTH_LAW_HEAVISIDE_OXY")
          vecfunc = Teuchos::rcp(new POROMULTIPHASESCATRA::TumorGrowthLawHeavisideOxy(params));
        else if (type == "TUMOR_GROWTH_LAW_HEAVISIDE_NECRO")
          vecfunc = Teuchos::rcp(new POROMULTIPHASESCATRA::TumorGrowthLawHeavisideNecro(params));
        else if (type == "OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_CONT")
          vecfunc =
              Teuchos::rcp(new POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawCont(params));
        else if (type == "OXYGEN_TRANSVASCULAR_EXCHANGE_LAW_DISC")
          vecfunc =
              Teuchos::rcp(new POROMULTIPHASESCATRA::OxygenTransvascularExchangeLawDisc(params));
        else
          dserror("Wrong type of POROMULTIPHASESCATRA_FUNCTION");

        functions_.push_back(vecfunc);
      }
      else if (DRT::UTILS::CombustFunctionHaveNamed(function, &functions_))
      {
      }
      else if (DRT::UTILS::XfluidFunctionHaveNamed(function, &functions_))
      {
      }
      else
      {
        // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        // define a new vector of functions
        Teuchos::RCP<ExprFunction> vecfunc = Teuchos::rcp(new ExprFunction());

        // evaluate the maximum component and the number of variables
        int maxcomp = 0;
        int maxvar = -1;
        for (unsigned int n = 0; n < functions.size(); ++n)
        {
          Teuchos::RCP<DRT::INPUT::LineDefinition> functnumcompvar = functions[n];
          functnumcompvar->ExtractInt("COMPONENT", maxcomp);
          functnumcompvar->ExtractInt("VARIABLE", maxvar);
        }

        // evaluate the number of rows used for the definition of the variables
        int numrowsvar;
        numrowsvar = functions.size() - maxcomp - 1;

        // define a vector of strings
        std::vector<std::string> functstring(maxcomp + 1);

        // read each row where the components of the i-th function are defined
        for (int n = 0; n <= maxcomp; ++n)
        {
          // update the current row
          Teuchos::RCP<DRT::INPUT::LineDefinition> functcomp = functions[n];

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
        for (int j = 1; j <= numrowsvar; ++j)
        {
          // update the current row
          Teuchos::RCP<DRT::INPUT::LineDefinition> timevar = functions[maxcomp + j];

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
          struct periodicstruct periodicdata;
          periodicdata.periodic = timevar->FindString("PERIODIC");
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
            // read the number of points
            int numpoints;
            timevar->ExtractInt("NUMPOINTS", numpoints);

            // read times
            std::vector<double> times;
            bool bynum = timevar->FindString("BYNUM");

            if (bynum)  // times defined by number of points
            {
              // read the time range
              std::vector<double> timerange;
              timevar->ExtractDoubleVector("TIMERANGE", timerange);

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
            }
            else  // times defined by vector
            {
              timevar->ExtractDoubleVector("TIMES", times);
            }

            // check if the times are in ascending order
            for (unsigned int k = 1; k < times.size(); ++k)
            {
              if (times[k] <= times[k - 1]) dserror("the TIMES must be in ascending order");
            }

            // read values
            std::vector<double> values;
            timevar->ExtractDoubleVector("VALUES", values);

            vardef =
                Teuchos::rcp(new LinearInterpolationVariable(varname, times, values, periodicdata));
          }
          else if (vartype == "multifunction")
          {
            // read the number of points
            int numpoints;
            timevar->ExtractInt("NUMPOINTS", numpoints);

            // read times
            std::vector<double> times;
            bool bynum = timevar->FindString("BYNUM");

            if (bynum)  // times defined by number of points
            {
              // read the time range
              std::vector<double> timerange;
              timevar->ExtractDoubleVector("TIMERANGE", timerange);

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
            }
            else  // times defined by vector
            {
              timevar->ExtractDoubleVector("TIMES", times);
            }

            // check if the times are in ascending order
            for (unsigned int k = 1; k < times.size(); ++k)
            {
              if (times[k] <= times[k - 1]) dserror("the TIMES must be in ascending order");
            }

            // read descriptions (strings separated with spaces)
            std::vector<std::string> description_vec;
            timevar->ExtractStringVector("DESCRIPTION", description_vec);

            // check if the number of times = number of descriptions + 1
            int numtimes = times.size();
            int numdescriptions = description_vec.size();
            if (numtimes != numdescriptions + 1)
              dserror("the number of TIMES and the number of DESCRIPTIONs must be consistent");

            vardef = Teuchos::rcp(
                new MultiFunctionVariable(varname, times, description_vec, periodicdata));
          }
          else if (vartype == "fourierinterpolation")
          {
            // read the number of points
            int numpoints;
            timevar->ExtractInt("NUMPOINTS", numpoints);

            // read times
            std::vector<double> times;
            bool bynum = timevar->FindString("BYNUM");

            if (bynum)  // times defined by number of points
            {
              // read the time range
              std::vector<double> timerange;
              timevar->ExtractDoubleVector("TIMERANGE", timerange);

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
            }
            else  // times defined by vector
            {
              timevar->ExtractDoubleVector("TIMES", times);
            }

            // check if the times are in ascending order
            for (unsigned int k = 1; k < times.size(); ++k)
            {
              if (times[k] <= times[k - 1]) dserror("the TIMES must be in ascending order");
            }

            // read values
            std::vector<double> values;
            timevar->ExtractDoubleVector("VALUES", values);

            vardef = Teuchos::rcp(
                new FourierInterpolationVariable(varname, times, values, periodicdata));
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

        functions_.push_back(vecfunc);
      }
    }
  }
}


DRT::UTILS::Function& DRT::UTILS::FunctionManager::Funct(int num)
{
  // ensure that desired function is available (prevents segmentation fault)
  if (functions_.size() < (unsigned int)(num + 1) || num < 0)
    dserror("function %d not available", num + 1);

  return *(functions_[num]);
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::ExprFunction::ExprFunction()
{
  dim_ = DRT::Problem::Instance()->NDim();
  expr_.clear();
  exprdd_.clear();
  variables_.clear();
  isparsed_ = false;
}


DRT::UTILS::ExprFunction::~ExprFunction() {}


void DRT::UTILS::ExprFunction::AddExpr(
    std::string buf, std::vector<std::vector<Teuchos::RCP<FunctionVariable>>> variables)
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

  return;
}


double DRT::UTILS::ExprFunction::Evaluate(const int index, const double* x, double t)
{
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  if (not isparsed_) ParseExpressions();

  double index_mod = index;

  if (expr_.size() == 1)
  {
    index_mod = 0;
  }

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
    unsigned int n = 0;
    bool containtime = false;
    while (containtime == false)
    {
      if (n == variables_[i].size())
      {
        dserror("the VARIABLE %d is not defined at time %f", i, t);
      }
      else
      {
        containtime = variables_[i][n]->ContainTime(t);
      }
      if (containtime == false)
      {
        ++n;
      }
    }

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

  double index_mod = index;

  if (expr_.size() == 1)
  {
    index_mod = 0;
  }

  // define Fad object for evaluation
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double>> FAD;

  // define FAD variables
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  FAD xfad(4 + variables_.size(), 0, x[0]);
  FAD yfad(4 + variables_.size(), 1, x[1]);
  FAD zfad(4 + variables_.size(), 2, x[2]);
  FAD tfad(4 + variables_.size(), 3, t);

  xfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 0, x[0]);
  yfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 1, x[1]);
  zfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 2, x[2]);
  tfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 3, t);

  std::vector<FAD> fadvectvars(variables_.size());
  for (unsigned int i = 0; i < variables_.size(); ++i)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int n = 0;
    bool containtime = false;
    while (containtime == false)
    {
      if (n == variables_[i].size())
      {
        dserror("the variable %d is not defined in the time considered", i);
      }
      else
      {
        containtime = variables_[i][n]->ContainTime(t);
      }
      if (containtime == false)
      {
        ++n;
      }
    }
    fadvectvars[i] = FAD(4 + variables_.size(), 4 + i, variables_[i][n]->Value(t));
    fadvectvars[i].val() =
        Sacado::Fad::DFad<double>(4 + variables_.size(), 4 + i, variables_[i][n]->Value(t));
  }
  FAD fdfad;

  // set spatial variables
  switch (dim_)
  {
    case 3:
    {
      exprdd_[index_mod]->SetValue("x", xfad);
      exprdd_[index_mod]->SetValue("y", yfad);
      exprdd_[index_mod]->SetValue("z", zfad);
      break;
    }
    case 2:
    {
      exprdd_[index_mod]->SetValue("x", xfad);
      exprdd_[index_mod]->SetValue("y", yfad);
      exprdd_[index_mod]->SetValue("z", 0);
      break;
    }
    case 1:
    {
      exprdd_[index_mod]->SetValue("x", xfad);
      exprdd_[index_mod]->SetValue("y", 0);
      exprdd_[index_mod]->SetValue("z", 0);
      break;
    }
    default:
      dserror("Problem dimension has to be 1, 2, or 3.");
      break;
  }

  // set temporal variable
  exprdd_[index_mod]->SetValue("t", tfad);

  // set the values of the variables at time t
  for (unsigned int i = 0; i < variables_.size(); ++i)
  {
    exprdd_[index_mod]->SetValue(variables_[i][0]->Name(), fadvectvars[i]);
  }

  // evaluation of derivatives
  fdfad = exprdd_[index_mod]->Evaluate();

  // result vector
  std::vector<double> res(3, 0.0);
  for (unsigned int d = 0; d < 3; ++d)
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

  double index_mod = index;
  if (expr_.size() == 1)
  {
    index_mod = 0;
  }

  // define Fad object for evaluation
  typedef Sacado::Fad::DFad<Sacado::Fad::DFad<double>> FAD;

  // define FAD variables
  // we consider a function of the type F = F ( x, y, z, t, v1(t), ..., vn(t) )
  FAD xfad(4 + variables_.size(), 0, x[0]);
  FAD yfad(4 + variables_.size(), 1, x[1]);
  FAD zfad(4 + variables_.size(), 2, x[2]);
  FAD tfad(4 + variables_.size(), 3, t);

  xfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 0, x[0]);
  yfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 1, x[1]);
  zfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 2, x[2]);
  tfad.val() = Sacado::Fad::DFad<double>(4 + variables_.size(), 3, t);

  std::vector<FAD> fadvectvars(variables_.size());
  for (unsigned int i = 0; i < variables_.size(); ++i)
  {
    // find the right definition of the variable according to the hierarchy
    unsigned int n = 0;
    bool containtime = false;
    while (containtime == false)
    {
      if (n == variables_[i].size())
      {
        dserror("the variable %d is not defined in the time considered", i);
      }
      else
      {
        containtime = variables_[i][n]->ContainTime(t);
      }
      if (containtime == false)
      {
        ++n;
      }
    }
    fadvectvars[i] = FAD(4 + variables_.size(), 4 + i, variables_[i][n]->Value(t));
    fadvectvars[i].val() =
        Sacado::Fad::DFad<double>(4 + variables_.size(), 4 + i, variables_[i][n]->Value(t));
  }
  FAD fdfad;

  // set spatial variables
  switch (dim_)
  {
    case 3:
    {
      exprdd_[index_mod]->SetValue("x", xfad);
      exprdd_[index_mod]->SetValue("y", yfad);
      exprdd_[index_mod]->SetValue("z", zfad);
      break;
    }
    case 2:
    {
      exprdd_[index_mod]->SetValue("x", xfad);
      exprdd_[index_mod]->SetValue("y", yfad);
      exprdd_[index_mod]->SetValue("z", 0);
      break;
    }
    case 1:
    {
      exprdd_[index_mod]->SetValue("x", xfad);
      exprdd_[index_mod]->SetValue("y", 0);
      exprdd_[index_mod]->SetValue("z", 0);
      break;
    }
    default:
      dserror("Problem dimension has to be 1, 2, or 3.");
      break;
  }

  // set temporal variable
  exprdd_[index_mod]->SetValue("t", tfad);

  // set the values of the variables at time t
  for (unsigned int i = 0; i < variables_.size(); ++i)
  {
    exprdd_[index_mod]->SetValue(variables_[i][0]->Name(), fadvectvars[i]);
  }

  // add the 1st time derivative at time t
  if (deg >= 1)
  {
    // evaluation of derivatives
    fdfad = exprdd_[index_mod]->Evaluate();

    // evaluation of dF/dt applying the chain rule:
    // dF/dt = dF*/dt + sum_i(dF/dvi*dvi/dt)
    double fdfad_dt = fdfad.dx(3).val();
    for (unsigned int i = 0; i < variables_.size(); ++i)
    {
      // find the right definition of the variable according to the hierarchy
      unsigned int n = 0;
      bool containtime = false;
      while (containtime == false)
      {
        containtime = variables_[i][n]->ContainTime(t);
        if (n == variables_[i].size())
        {
          dserror("the variable %d is not defined in the time considered", i);
        }
        if (containtime == false)
        {
          ++n;
        }
      }
      fdfad_dt += fdfad.dx(4 + i).val() * variables_[i][n]->TimeDerivativeValue(t);
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

    for (unsigned int i = 0; i < variables_.size(); ++i)
    {
      fdfad_dt2_term[i] = 0;
      // find the right definition of the variable according to the hierarchy
      unsigned int n = 0;
      bool containtime = false;
      while (containtime == false)
      {
        containtime = variables_[i][n]->ContainTime(t);
        if (n == variables_[i].size())
        {
          dserror("the variable %d is not defined in the time considered", i);
        }
        if (containtime == false)
        {
          ++n;
        }
      }
      fdfad_dt2_term[i] += fdfad.dx(3).dx(4 + i);
      fdfad_dt2_term[i] += fdfad.dx(4 + i).dx(3);
      for (unsigned int j = 0; j < variables_.size(); ++j)
      {
        // find the right definition of the variable according to the hierarchy
        unsigned int m = 0;
        bool containtime = false;
        while (containtime == false)
        {
          containtime = variables_[j][m]->ContainTime(t);
          if (m == variables_[j].size())
          {
            dserror("the variable %d is not defined in the time considered", j);
          }
          if (containtime == false)
          {
            ++m;
          }
        }
        fdfad_dt2_term[i] += fdfad.dx(4 + i).dx(4 + j) * variables_[j][m]->TimeDerivativeValue(t);
      }
      fdfad_dt2_term[i] *= variables_[i][n]->TimeDerivativeValue(t);
      fdfad_dt2_term[i] += fdfad.dx(4 + i).val() * variables_[i][n]->TimeDerivativeValue(t, 2);
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
  if (index > (int)expr_.size() - 1 || index < 0)
    dserror("Tried to add a variable to a function in a not available dimension.");

  return expr_[index]->IsVariable(varname);
}


void DRT::UTILS::ExprFunction::AddVariable(
    const int index, const std::string& varname, double varvalue)
{
  if (index > (int)expr_.size() - 1 || index < 0)
    dserror("Tried to add a variable to a function in a not available dimension.");
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
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
DRT::UTILS::VariableExprFunction::VariableExprFunction() : isparsed_(false) {}


DRT::UTILS::VariableExprFunction::~VariableExprFunction() {}


void DRT::UTILS::VariableExprFunction::AddExpr(
    std::string buf, std::vector<std::pair<std::string, double>> constants)
{
  // do the almost same as the expression function (base class) but do not yet parse!

  // build the parser for the function evaluation
  Teuchos::RCP<DRT::PARSER::Parser<double>> parser =
      Teuchos::rcp(new DRT::PARSER::Parser<double>(buf));

  // build the parser for the function derivative evaluation
  Teuchos::RCP<DRT::PARSER::Parser<Sacado::Fad::DFad<double>>> parserd =
      Teuchos::rcp(new DRT::PARSER::Parser<Sacado::Fad::DFad<double>>(buf));

  // add constants
  for (std::vector<std::pair<std::string, double>>::iterator it = constants.begin();
       it != constants.end(); it++)
    parser->AddVariable(it->first, it->second);
  for (std::vector<std::pair<std::string, double>>::iterator it = constants.begin();
       it != constants.end(); it++)
    parserd->AddVariable(it->first, it->second);

  // save the parsers
  expr_.push_back(parser);
  exprd_.push_back(parserd);

  isparsed_ = false;

  return;
}


bool DRT::UTILS::VariableExprFunction::IsVariable(int index, const std::string& varname) const
{
  if (index > (int)expr_.size() - 1 || index < 0)
    dserror("Tried to add a variable to a function in a not available dimension.");

  return expr_[index]->IsVariable(varname);
}


void DRT::UTILS::VariableExprFunction::AddVariable(
    int index, const std::string& varname, double varvalue)
{
  if (index > (int)expr_.size() - 1 || index < 0)
    dserror("Tried to add a variable to a function in a not available dimension.");
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
  // set the values of the variables
  std::vector<std::pair<std::string, double>>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
    expr_[index]->SetValue(it->first, it->second);

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}


double DRT::UTILS::VariableExprFunction::Evaluate(const int index,
    const std::vector<std::pair<std::string, double>>& variables,
    const std::vector<std::pair<std::string, double>>& constants)
{
  if (not isparsed_) ParseExpressions();
  // set the values of the variables
  std::vector<std::pair<std::string, double>>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
    expr_[index]->SetValue(it->first, it->second);
  // set the values of the constants
  for (it = constants.begin(); it != constants.end(); it++)
  {
    if (expr_[index]->IsVariable(it->first)) expr_[index]->SetValue(it->first, it->second);
  }

  // evaluate the function and return the result
  return expr_[index]->Evaluate();
}


std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateDerivative(
    const int index, const std::vector<std::pair<std::string, double>>& variables)
{
  if (not isparsed_) ParseExpressions();

  // Fad object for evaluation
  // sacado data type replaces "double"
  typedef Sacado::Fad::DFad<double> FAD;

  // number of variables
  int numvariables = variables.size();

  // counter for variable numbering
  int counter = 0;

  // set the values of the variables
  std::vector<std::pair<std::string, double>>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
  {
    // for 1st order derivatives
    FAD varfad(numvariables, counter, it->second);
    // set the value in expression
    exprd_[index]->SetValue(it->first, varfad);
    // update counter
    counter++;
  }

  // evaluate the expression
  FAD fdfad = exprd_[index]->Evaluate();

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

  // Fad object for evaluation
  // sacado data type replaces "double"
  typedef Sacado::Fad::DFad<double> FAD;

  // number of variables
  int numvariables = variables.size();

  // counter for variable numbering
  int counter = 0;

  // set the values of the variables
  std::vector<std::pair<std::string, double>>::const_iterator it;
  for (it = variables.begin(); it != variables.end(); it++)
  {
    // for 1st order derivatives
    FAD varfad(numvariables, counter, it->second);
    // set the value in expression
    exprd_[index]->SetValue(it->first, varfad);
    // update counter
    counter++;
  }

  // set the values of the constants
  for (it = constants.begin(); it != constants.end(); it++)
  {
    if (exprd_[index]->IsVariable(it->first))
      // set the value in expression
      exprd_[index]->SetValue(it->first, it->second);
  }

  // evaluate the expression
  FAD fdfad = exprd_[index]->Evaluate();

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
      variables.push_back(std::pair<std::string, double>("x", x[0]));  //-x_[index]
      variables.push_back(std::pair<std::string, double>("y", x[1]));  //-y_[index]
      variables.push_back(std::pair<std::string, double>("z", x[2]));  //-z_[index]
    }
    case 2:
    {
      variables.push_back(std::pair<std::string, double>("x", x[0]));  //-x_[index]
      variables.push_back(std::pair<std::string, double>("y", x[1]));  //-y_[index]
    }
    case 1:
    {
      variables.push_back(std::pair<std::string, double>("x", x[0]));  //-x_[index]
    }
  }

  variables.push_back(std::pair<std::string, double>("t", t));

  return Evaluate(index, variables);
}


std::vector<double> DRT::UTILS::VariableExprFunction::EvaluateSpatialDerivative(
    int index, const double* x, const double t)
{
  std::vector<std::pair<std::string, double>> variables(4);

  variables[0] = std::pair<std::string, double>("x", x[0]);
  variables[1] = std::pair<std::string, double>("y", x[1]);
  variables[2] = std::pair<std::string, double>("z", x[2]);
  variables[3] = std::pair<std::string, double>("t", t);

  return EvaluateDerivative(index, variables);
}
