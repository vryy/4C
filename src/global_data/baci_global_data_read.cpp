/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities to read and fill global data

\level 1

    */
/*----------------------------------------------------------------------*/

#include "baci_global_data_read.H"

#include "baci_comm_utils.H"
#include "baci_contact_constitutivelaw_bundle.H"
#include "baci_contact_constitutivelaw_constitutivelaw_definition.H"
#include "baci_global_data.H"
#include "baci_global_legacy_module.H"
#include "baci_inpar_validconditions.H"
#include "baci_inpar_validcontactconstitutivelaw.H"
#include "baci_inpar_validmaterials.H"
#include "baci_io.H"
#include "baci_io_inputreader.H"
#include "baci_io_linedefinition.H"
#include "baci_io_materialdefinition.H"
#include "baci_lib_conditiondefinition.H"
#include "baci_lib_discret_hdg.H"
#include "baci_lib_discret_xfem.H"
#include "baci_lib_discret_xwall.H"
#include "baci_lib_dofset_independent.H"
#include "baci_lib_elementreader.H"
#include "baci_lib_meshreader.H"
#include "baci_lib_nodereader.H"
#include "baci_lib_utils_createdis.H"
#include "baci_mat_elchmat.H"
#include "baci_mat_elchphase.H"
#include "baci_mat_micromaterial.H"
#include "baci_mat_newman_multiscale.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_scatra_mat_multiscale.H"
#include "baci_nurbs_discret.H"
#include "baci_particle_engine_particlereader.H"
#include "baci_rebalance.H"
#include "baci_utils_function.H"

BACI_NAMESPACE_OPEN

void GLOBAL::ReadFields(GLOBAL::Problem& problem, INPUT::DatFileReader& reader, const bool readmesh)
{
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> xfluiddis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> structaledis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> thermdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> lubricationdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> scatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> scatra_micro_dis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> cellscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> fluidscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> structscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> artscatradis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> arterydis = Teuchos::null;  //_1D_ARTERY_
  Teuchos::RCP<DRT::Discretization> airwaydis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> optidis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> porofluiddis = Teuchos::null;  // fpsi, poroelast
  Teuchos::RCP<DRT::Discretization> elemagdis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> celldis = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> pboxdis = Teuchos::null;

  // decide which kind of spatial representation is required
  const CORE::FE::ShapeFunctionType distype = problem.SpatialApproximationType();
  auto output_control = problem.OutputControlFile();

  // the basic mesh reader. now add desired node and element readers to it!
  INPUT::MeshReader meshreader(reader, "--NODE COORDS");

  switch (problem.GetProblemType())
  {
    case GLOBAL::ProblemType::fsi:
    case GLOBAL::ProblemType::fsi_redmodels:
    case GLOBAL::ProblemType::fsi_lung:
    {
      if (distype == CORE::FE::ShapeFunctionType::nurbs)
      {
        structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
      }
      else if (INPUT::IntegralValue<int>(
                   problem.FluidDynamicParams().sublist("WALL MODEL"), "X_WALL"))
      {
        structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXWall("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }
      else
      {
        structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
        if (INPUT::IntegralValue<bool>(
                problem.XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
          xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      if (xfluiddis != Teuchos::null)
        xfluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("fluid", fluiddis);
      if (xfluiddis != Teuchos::null) problem.AddDis("xfluid", xfluiddis);
      problem.AddDis("ale", aledis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      if (xfluiddis != Teuchos::null)
      {
        meshreader.AddElementReader(INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS"));
      }
      else
        meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::gas_fsi:
    case GLOBAL::ProblemType::ac_fsi:
    case GLOBAL::ProblemType::thermo_fsi:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          dserror("Nurbs discretization not possible for fs3i!");
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1", reader.Comm()));
          structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));
      fluidscatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis, output_control, distype)));
      structscatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structscatradis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("fluid", fluiddis);
      problem.AddDis("ale", aledis);
      problem.AddDis("scatra1", fluidscatradis);
      problem.AddDis("scatra2", structscatradis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(
          INPUT::ElementReader(fluidscatradis, reader, "--TRANSPORT ELEMENTS"));
      meshreader.AddElementReader(
          INPUT::ElementReader(structscatradis, reader, "--TRANSPORT2 ELEMENTS"));

#ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false, false, false);
#endif

      break;
    }
    case GLOBAL::ProblemType::biofilm_fsi:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          dserror("Nurbs discretization not possible for biofilm problems!");
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          structaledis = Teuchos::rcp(new DRT::Discretization("structale", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));
      structaledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structaledis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("fluid", fluiddis);
      problem.AddDis("ale", aledis);
      problem.AddDis("structale", structaledis);


      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

#ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false, false, false);
#endif

      // fluid scatra field
      fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis, output_control, distype)));
      problem.AddDis("scatra1", fluidscatradis);

      // structure scatra field
      structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structscatradis, output_control, distype)));
      problem.AddDis("scatra2", structscatradis);

      break;
    }
    case GLOBAL::ProblemType::fsi_xfem:
    case GLOBAL::ProblemType::fluid_xfem:
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      problem.AddDis("structure", structdis);
      meshreader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.StructuralDynamicParams(), "GEOMETRY"),
          nullptr);

      if (INPUT::IntegralValue<int>(
              problem.XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
        fluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
        problem.AddDis("fluid", fluiddis);

        xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
        xfluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis, output_control, distype)));
        problem.AddDis("xfluid", xfluiddis);

        meshreader.AddElementReader(
            INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID"));
      }
      else
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid", reader.Comm()));
        fluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
        problem.AddDis("fluid", fluiddis);

        meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
            INPUT::IntegralValue<INPAR::GeometryType>(problem.FluidDynamicParams(), "GEOMETRY"),
            nullptr);
      }

      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));
      problem.AddDis("ale", aledis);
      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));
      break;
    }
    case GLOBAL::ProblemType::fpsi_xfem:
    {
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("porofluid", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("porofluid", porofluiddis);
      problem.AddDis("fluid", fluiddis);
      problem.AddDis("ale", aledis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.FluidDynamicParams(), "GEOMETRY"),
          nullptr);

      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::ale:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
          break;
        }
        default:
        {
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.AddDis("ale", aledis);

      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::fluid:
    case GLOBAL::ProblemType::fluid_redmodels:
    {
      if (distype == CORE::FE::ShapeFunctionType::hdg)
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationHDG("fluid", reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }
      else if (distype == CORE::FE::ShapeFunctionType::nurbs)
      {
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));

        // create discretization writer - in constructor set ingto and owned by corresponding
        // discret
        fluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }
      else if (INPUT::IntegralValue<int>(
                   problem.FluidDynamicParams().sublist("WALL MODEL"), "X_WALL"))
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXWall("fluid", reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }
      else
      {
        // fluiddis  = Teuchos::rcp(new DRT::Discretization("fluid",reader.Comm()));
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }

      problem.AddDis("fluid", fluiddis);

      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.FluidDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }
    case GLOBAL::ProblemType::lubrication:
    {
      // create empty discretizations
      lubricationdis = Teuchos::rcp(new DRT::Discretization("lubrication", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      lubricationdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(lubricationdis, output_control, distype)));

      problem.AddDis("lubrication", lubricationdis);

      meshreader.AddElementReader(
          INPUT::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::cardiac_monodomain:
    case GLOBAL::ProblemType::scatra:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra", reader.Comm()));
          break;
        }
        case CORE::FE::ShapeFunctionType::hdg:
        {
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::DiscretizationHDG("scatra", reader.Comm()));
          break;
        }
        default:
        {
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));


      problem.AddDis("fluid", fluiddis);
      problem.AddDis("scatra", scatradis);

      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::sti:
    {
      // safety checks
      if (distype == CORE::FE::ShapeFunctionType::nurbs)
        dserror("Scatra-thermo interaction does not work for nurbs discretizations yet!");

      // create empty discretizations for scalar and thermo fields
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
      thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));

      // create discretization writers
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));
      thermdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(thermdis, output_control, distype)));

      // add empty discretizations to global problem
      problem.AddDis("scatra", scatradis);
      problem.AddDis("thermo", thermdis);

      // add element reader to node reader
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::fluid_ale:
    case GLOBAL::ProblemType::freesurf:
    {
      if (distype == CORE::FE::ShapeFunctionType::hdg)
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationHDG("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }
      else if (distype == CORE::FE::ShapeFunctionType::nurbs)
      {
        fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
      }
      else if (INPUT::IntegralValue<int>(
                   problem.FluidDynamicParams().sublist("WALL MODEL"), "X_WALL"))
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXWall("fluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }
      else
      {
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
        if (INPUT::IntegralValue<bool>(
                problem.XFluidDynamicParams().sublist("GENERAL"), "XFLUIDFLUID"))
          xfluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("xfluid", reader.Comm()));
        aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
      }


      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      if (xfluiddis != Teuchos::null)
      {
        xfluiddis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(xfluiddis, output_control, distype)));
      }
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));


      problem.AddDis("fluid", fluiddis);
      if (xfluiddis != Teuchos::null)
      {
        problem.AddDis("xfluid", xfluiddis);  // xfem discretization on slot 1
      }
      problem.AddDis("ale", aledis);

      if (xfluiddis != Teuchos::null)
      {
        meshreader.AddElementReader(INPUT::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS"));
      }
      else
        meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::tsi:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          thermdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("thermo", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      thermdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(thermdis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("thermo", thermdis);

      meshreader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.StructuralDynamicParams(), "GEOMETRY"),
          nullptr);
      meshreader.AddAdvancedReader(thermdis, reader, "THERMO",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.ThermalDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }
    case GLOBAL::ProblemType::thermo:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          thermdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("thermo", reader.Comm()));
          break;
        }
        default:
        {
          thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      thermdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(thermdis, output_control, distype)));

      problem.AddDis("thermo", thermdis);

      meshreader.AddElementReader(INPUT::ElementReader(thermdis, reader, "--THERMO ELEMENTS"));

      break;
    }

    case GLOBAL::ProblemType::structure:
    {
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));

      problem.AddDis("structure", structdis);

      meshreader.AddAdvancedReader(structdis, reader, "STRUCTURE",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.StructuralDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }

    case GLOBAL::ProblemType::polymernetwork:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      pboxdis = Teuchos::rcp(new DRT::Discretization("boundingbox", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      pboxdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(pboxdis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("boundingbox", pboxdis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(
          INPUT::ElementReader(pboxdis, reader, "--PERIODIC BOUNDINGBOX ELEMENTS"));

      break;
    }

    case GLOBAL::ProblemType::loma:
    {
      // create empty discretizations
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.AddDis("fluid", fluiddis);
      problem.AddDis("scatra", scatradis);

      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }

    case GLOBAL::ProblemType::fluid_xfem_ls:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      if (problem.GetProblemType() == GLOBAL::ProblemType::fluid_xfem_ls)
        fluiddis = Teuchos::rcp(new DRT::DiscretizationXFEM("fluid", reader.Comm()));
      else
        fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("fluid", fluiddis);
      problem.AddDis("scatra", scatradis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.FluidDynamicParams(), "GEOMETRY"),
          nullptr);
      // meshreader.AddElementReader(Teuchos::rcp(new INPUT::ElementReader(fluiddis, reader,
      // "--FLUID ELEMENTS")));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      break;
    }

    case GLOBAL::ProblemType::elch:
    {
      // create empty discretizations
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          fluiddis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("ale", reader.Comm()));
          scatra_micro_dis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra_micro", reader.Comm()));
          break;
        }
        default:
        {
          fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
          aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));
          scatra_micro_dis = Teuchos::rcp(new DRT::Discretization("scatra_micro", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));
      scatra_micro_dis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatra_micro_dis, output_control, distype)));

      problem.AddDis("fluid", fluiddis);
      problem.AddDis("scatra", scatradis);
      problem.AddDis("ale", aledis);
      problem.AddDis("scatra_micro", scatra_micro_dis);

      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));
      meshreader.AddElementReader(
          INPUT::ElementReader(scatra_micro_dis, reader, "--TRANSPORT2 ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::art_net:  // _1D_ARTERY_
    {
      // create empty discretizations
      arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));

      // create empty discretizations
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          dserror("Nurbs Discretization not possible for artery");
          break;
        }
        default:
        {
          scatradis = Teuchos::rcp(new DRT::Discretization("artery_scatra", reader.Comm()));
          break;
        }
      }

      problem.AddDis("artery", arterydis);
      problem.AddDis("artery_scatra", scatradis);

      // create discretization writer - in constructor set into and owned by corresponding discret
      arterydis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(arterydis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      meshreader.AddElementReader(INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::red_airways:  // _reduced D airways
    {
      // create empty discretizations
      airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      airwaydis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(airwaydis, output_control, distype)));

      problem.AddDis("red_airway", airwaydis);

      meshreader.AddElementReader(
          INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::struct_ale:  // structure with ale
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("ale", aledis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::poroelast:
    case GLOBAL::ProblemType::poromultiphase:
    {
      // create empty discretizations
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          porofluiddis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("porofluid", porofluiddis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));

      if (INPUT::IntegralValue<bool>(problem.PoroMultiPhaseDynamicParams(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        arterydis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.AddDis("artery", arterydis);
        meshreader.AddElementReader(INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));
      }

      break;
    }
    case GLOBAL::ProblemType::poromultiphasescatra:
    {
      // create empty discretizations
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("structure", reader.Comm()));
          porofluiddis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("scatra", reader.Comm()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
          porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
          scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("porofluid", porofluiddis);
      problem.AddDis("scatra", scatradis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      if (INPUT::IntegralValue<bool>(
              problem.PoroMultiPhaseScatraDynamicParams(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        arterydis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.AddDis("artery", arterydis);
        meshreader.AddElementReader(INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));

        artscatradis = Teuchos::rcp(new DRT::Discretization("artery_scatra", reader.Comm()));
        artscatradis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(artscatradis, output_control, distype)));
        problem.AddDis("artery_scatra", artscatradis);
        meshreader.AddElementReader(
            INPUT::ElementReader(artscatradis, reader, "--TRANSPORT ELEMENTS"));
      }

      break;
    }
    case GLOBAL::ProblemType::porofluidmultiphase:
    {
      // create empty discretizations
      switch (distype)
      {
        case CORE::FE::ShapeFunctionType::nurbs:
        {
          porofluiddis =
              Teuchos::rcp(new DRT::NURBS::NurbsDiscretization("porofluid", reader.Comm()));
          break;
        }
        default:
        {
          porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));

      problem.AddDis("porofluid", porofluiddis);

      meshreader.AddElementReader(INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));

      if (INPUT::IntegralValue<bool>(problem.PoroFluidMultiPhaseDynamicParams(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        arterydis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.AddDis("artery", arterydis);
        meshreader.AddElementReader(INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));
      }
      break;
    }
    case GLOBAL::ProblemType::fpsi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("porofluid", porofluiddis);
      problem.AddDis("fluid", fluiddis);
      problem.AddDis("ale", aledis);

      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::fbi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("fluid", fluiddis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddAdvancedReader(fluiddis, reader, "FLUID",
          INPUT::IntegralValue<INPAR::GeometryType>(problem.FluidDynamicParams(), "GEOMETRY"),
          nullptr);

      break;
    }
    case GLOBAL::ProblemType::immersed_fsi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("fluid", fluiddis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::fps3i:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
      fluiddis = Teuchos::rcp(new DRT::DiscretizationFaces("fluid", reader.Comm()));
      aledis = Teuchos::rcp(new DRT::Discretization("ale", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      fluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("porofluid", porofluiddis);
      problem.AddDis("fluid", fluiddis);
      problem.AddDis("ale", aledis);


      meshreader.AddElementReader(INPUT::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      // fluid scatra field
      fluidscatradis = Teuchos::rcp(new DRT::Discretization("scatra1", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(fluidscatradis, output_control, distype)));
      problem.AddDis("scatra1", fluidscatradis);

      // poro structure scatra field
      structscatradis = Teuchos::rcp(new DRT::Discretization("scatra2", reader.Comm()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structscatradis, output_control, distype)));
      problem.AddDis("scatra2", structscatradis);

      break;
    }
    case GLOBAL::ProblemType::poroscatra:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      porofluiddis = Teuchos::rcp(new DRT::Discretization("porofluid", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("porofluid", porofluiddis);
      problem.AddDis("scatra", scatradis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      break;
    }
    case GLOBAL::ProblemType::ehl:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      lubricationdis = Teuchos::rcp(new DRT::Discretization("lubrication", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      lubricationdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(lubricationdis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("lubrication", lubricationdis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(
          INPUT::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::ssi:
    case GLOBAL::ProblemType::ssti:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("scatra", scatradis);

      // consider case of additional scatra manifold
      if (INPUT::IntegralValue<bool>(
              problem.SSIControlParams().sublist("MANIFOLD"), "ADD_MANIFOLD"))
      {
        auto scatra_manifold_dis =
            Teuchos::rcp(new DRT::Discretization("scatra_manifold", reader.Comm()));
        scatra_manifold_dis->SetWriter(Teuchos::rcp(
            new IO::DiscretizationWriter(scatra_manifold_dis, output_control, distype)));
        problem.AddDis("scatra_manifold", scatra_manifold_dis);
      }

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      if (problem.GetProblemType() == GLOBAL::ProblemType::ssti)
      {
        thermdis = Teuchos::rcp(new DRT::Discretization("thermo", reader.Comm()));
        thermdis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(thermdis, output_control, distype)));
        problem.AddDis("thermo", thermdis);
        meshreader.AddElementReader(INPUT::ElementReader(thermdis, reader, "--TRANSPORT ELEMENTS"));
      }

      break;
    }
    case GLOBAL::ProblemType::particle:
    case GLOBAL::ProblemType::pasi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));

      problem.AddDis("structure", structdis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      break;
    }
    case GLOBAL::ProblemType::level_set:
    {
      // create empty discretizations
      scatradis = Teuchos::rcp(new DRT::Discretization("scatra", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      scatradis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.AddDis("scatra", scatradis);

      meshreader.AddElementReader(INPUT::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      break;
    }
    case GLOBAL::ProblemType::np_support:
    {
      // no discretizations and nodes needed for supporting procs
      break;
    }
    case GLOBAL::ProblemType::elemag:
    {
      // create empty discretizations
      elemagdis = Teuchos::rcp(new DRT::DiscretizationHDG("elemag", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      elemagdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(elemagdis, output_control, distype)));

      problem.AddDis("elemag", elemagdis);

      std::set<std::string> elemagelementtypes;
      elemagelementtypes.insert("ELECTROMAGNETIC");
      elemagelementtypes.insert("ELECTROMAGNETICDIFF");

      meshreader.AddElementReader(INPUT::ElementReader(
          elemagdis, reader, "--ELECTROMAGNETIC ELEMENTS", elemagelementtypes));

      break;
    }
    case GLOBAL::ProblemType::redairways_tissue:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(new DRT::Discretization("structure", reader.Comm()));
      airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway", reader.Comm()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(structdis, output_control, distype)));
      airwaydis->SetWriter(
          Teuchos::rcp(new IO::DiscretizationWriter(airwaydis, output_control, distype)));

      problem.AddDis("structure", structdis);
      problem.AddDis("red_airway", airwaydis);

      meshreader.AddElementReader(INPUT::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.AddElementReader(
          INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS"));
    }
    break;
    default:
      dserror("Unknown problem type: %d", problem.GetProblemType());
      break;
  }

  // add artery or airways discretizations only for the following problem types
  switch (problem.GetProblemType())
  {
    case GLOBAL::ProblemType::fsi_redmodels:
    case GLOBAL::ProblemType::fsi_lung:
    case GLOBAL::ProblemType::fluid_ale:
    case GLOBAL::ProblemType::fluid_redmodels:
    {
      if (distype == CORE::FE::ShapeFunctionType::polynomial)
      {
        // create empty discretizations
        arterydis = Teuchos::rcp(new DRT::Discretization("artery", reader.Comm()));
        // create discretization writer - in constructor set into and owned by corresponding discret
        arterydis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.AddDis("artery", arterydis);
        meshreader.AddElementReader(INPUT::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));

        airwaydis = Teuchos::rcp(new DRT::Discretization("red_airway", reader.Comm()));
        // create discretization writer - in constructor set into and owned by corresponding discret
        airwaydis->SetWriter(
            Teuchos::rcp(new IO::DiscretizationWriter(airwaydis, output_control, distype)));
        problem.AddDis("red_airway", airwaydis);
        meshreader.AddElementReader(
            INPUT::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS"));
      }
    }
    break;
    default:
      break;
  }

  if (readmesh)  // now read and allocate!
  {
    // we read nodes and elements for the desired fields as specified above
    meshreader.ReadAndPartition();

    CORE::COMM::NestedParallelismType npType =
        GLOBAL::Problem::Instance()->GetCommunicators()->NpType();
    // care for special applications
    switch (problem.GetProblemType())
    {
      case GLOBAL::ProblemType::elch:
      case GLOBAL::ProblemType::fsi:
      case GLOBAL::ProblemType::fsi_redmodels:
      case GLOBAL::ProblemType::fsi_lung:
      case GLOBAL::ProblemType::scatra:
      case GLOBAL::ProblemType::structure:
      {
        // read microscale fields from second, third, ... input file if necessary
        // (in case of multi-scale material models)
        if (npType != CORE::COMM::NestedParallelismType::copy_dat_file)
          ReadMicroFields(problem, reader);
        break;
      }
      case GLOBAL::ProblemType::np_support:
      {
        // read microscale fields from second, third, ... inputfile for supporting processors
        ReadMicrofieldsNPsupport(problem);
        break;
      }
      default:
        break;
    }
  }  // if(readmesh)
}

void GLOBAL::ReadMicroFields(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  // check whether micro material is specified
  const int id_struct =
      GLOBAL::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_struct_multiscale);
  const int id_scatra =
      GLOBAL::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_scatra_multiscale);
  const int id_elch =
      GLOBAL::Problem::Instance()->Materials()->FirstIdByType(INPAR::MAT::m_newman_multiscale);

  // return if no multiscale material is used
  if (id_struct == -1 and id_scatra == -1 and id_elch == -1) return;

  // safety check
  if ((id_struct != -1 and id_scatra != -1) or (id_struct != -1 and id_elch != -1) or
      (id_scatra != -1 and id_elch != -1))
    dserror("Cannot have more than one multi-scale material!");

  // store name of macro-scale discretization in string
  std::string macro_dis_name("");
  if (id_struct != -1)
    macro_dis_name = "structure";
  else
    macro_dis_name = "scatra";

  // fetch communicators
  Teuchos::RCP<Epetra_Comm> lcomm = problem.GetCommunicators()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem.GetCommunicators()->GlobalComm();

  GLOBAL::Problem* macro_problem = GLOBAL::Problem::Instance();
  Teuchos::RCP<DRT::Discretization> macro_dis = macro_problem->GetDis(macro_dis_name);

  // repartition macro problem for a good distribution of elements with micro material
  if (macro_dis_name == "structure")
  {
    // do weighted repartitioning to obtain new row/column maps
    const Teuchos::ParameterList rebalanceParams;
    Teuchos::RCP<const Epetra_CrsGraph> nodeGraph = macro_dis->BuildNodeGraph();
    const auto& [nodeWeights, edgeWeights] = CORE::REBALANCE::BuildWeights(*macro_dis);
    const auto& [rownodes, colnodes] =
        CORE::REBALANCE::RebalanceNodeMaps(nodeGraph, rebalanceParams, nodeWeights, edgeWeights);

    // rebuild the discretization with new maps
    macro_dis->Redistribute(*rownodes, *colnodes, true, true, true);
  }

  // make sure that we read the micro discretizations only on the processors on
  // which elements with the corresponding micro material are evaluated

  std::set<int> my_multimat_IDs;

  // take care also of ghosted elements! -> ElementColMap!
  for (int i = 0; i < macro_dis->ElementColMap()->NumMyElements(); ++i)
  {
    DRT::Element* actele = macro_dis->lColElement(i);
    Teuchos::RCP<MAT::Material> actmat = actele->Material();

    if (id_elch != -1 and actmat->MaterialType() == INPAR::MAT::m_elchmat)
    {
      // extract wrapped material
      auto elchmat = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(actmat);
      auto elchphase =
          Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(elchmat->PhaseById(elchmat->PhaseID(0)));
      actmat = elchphase->MatById(elchphase->MatID(0));
    }

    if ((actmat->MaterialType() == INPAR::MAT::m_struct_multiscale and
            macro_dis_name == "structure") or
        (actmat->MaterialType() == INPAR::MAT::m_scatra_multiscale and
            macro_dis_name == "scatra") or
        (actmat->MaterialType() == INPAR::MAT::m_newman_multiscale and macro_dis_name == "scatra"))
    {
      MAT::PAR::Parameter* actparams = actmat->Parameter();
      my_multimat_IDs.insert(actparams->Id());
    }
  }

  // check which macro procs have an element with micro material
  int foundmicromat = 0;
  int foundmicromatmyrank = -1;
  if (my_multimat_IDs.size() != 0)
  {
    foundmicromat = 1;
    foundmicromatmyrank = lcomm->MyPID();
  }

  // find out how many procs have micro material
  int nummicromat = 0;
  lcomm->SumAll(&foundmicromat, &nummicromat, 1);
  // broadcast number of procs that have micro material
  gcomm->Broadcast(&nummicromat, 1, 0);

  // every proc needs to know which procs have micro material in order to distribute colors
  // array is filled with either its local proc id or -1 when no micro mat was found
  std::vector<int> foundmyranks;
  foundmyranks.resize(lcomm->NumProc(), -1);
  lcomm->GatherAll(&foundmicromatmyrank, foundmyranks.data(), 1);

  // determine color of macro procs with any contribution to micro material, only important for
  // procs with micro material color starts with 0 and is incremented for each group
  int color = -1;
  if (foundmicromat == 1)
  {
    for (int foundmyrank : foundmyranks)
    {
      if (foundmyrank != -1) ++color;
      if (foundmyrank == foundmicromatmyrank) break;
    }
  }
  else
  {
    color = MPI_UNDEFINED;
  }

  // do the splitting of the communicator (macro proc must always be proc in subcomm with lowest
  // key
  // --> 0 is inserted here)
  MPI_Comm mpi_local_comm;
  MPI_Comm_split((Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm, true)->GetMpiComm()), color,
      0 /*important here*/, &mpi_local_comm);

  // sort out macro procs that do not have micro material
  if (foundmicromat == 1)
  {
    // create the sub communicator that includes one macro proc and some supporting procs
    Teuchos::RCP<Epetra_Comm> subgroupcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
    problem.GetCommunicators()->SetSubComm(subgroupcomm);

    // find out how many micro problems have to be solved on this macro proc
    int microcount = 0;
    for (const auto& material_map : *problem.Materials()->Map())
    {
      int matid = material_map.first;
      if (my_multimat_IDs.find(matid) != my_multimat_IDs.end()) microcount++;
    }
    // and broadcast it to the corresponding group of procs
    subgroupcomm->Broadcast(&microcount, 1, 0);

    for (const auto& material_map : *problem.Materials()->Map())
    {
      int matid = material_map.first;

      if (my_multimat_IDs.find(matid) != my_multimat_IDs.end())
      {
        Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);

        // initialize variables storing micro-scale information
        int microdisnum(-1);
        std::string micro_dis_name = "";
        std::string micro_inputfile_name("");
        GLOBAL::Problem* micro_problem(nullptr);

        // structure case
        if (macro_dis_name == "structure")
        {
          // access multi-scale structure material
          auto* micromat = static_cast<MAT::MicroMaterial*>(mat.get());

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->MicroDisNum();
          subgroupcomm->Broadcast(&microdisnum, 1, 0);

          // set name of micro-scale discretization
          micro_dis_name = "structure";

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->MicroInputFileName();

          // instantiate micro-scale problem
          micro_problem = GLOBAL::Problem::Instance(microdisnum);
        }

        // scalar transport case
        else
        {
          // access multi-scale scalar transport material
          MAT::ScatraMultiScale* micromat = nullptr;
          if (id_scatra != -1)
            micromat = dynamic_cast<MAT::ScatraMatMultiScale*>(mat.get());
          else if (id_elch != -1)
            micromat = dynamic_cast<MAT::NewmanMultiScale*>(mat.get());
          else
            dserror("How the heck did you get here?!");

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->MicroDisNum();
          subgroupcomm->Broadcast(&microdisnum, 1, 0);

          // set unique name of micro-scale discretization
          std::stringstream name;
          name << "scatra_multiscale_" << microdisnum;
          micro_dis_name = name.str();

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->MicroInputFileName();

          // instantiate micro-scale problem
          micro_problem = GLOBAL::Problem::Instance(microdisnum);
        }

        if (micro_inputfile_name[0] != '/')
        {
          std::string filename = reader.MyInputfileName();
          std::string::size_type pos = filename.rfind('/');
          if (pos != std::string::npos)
          {
            std::string path = filename.substr(0, pos + 1);
            micro_inputfile_name.insert(micro_inputfile_name.begin(), path.begin(), path.end());
          }
        }

        // broadcast micro input file name
        int length = static_cast<int>(micro_inputfile_name.length());
        subgroupcomm->Broadcast(&length, 1, 0);
        subgroupcomm->Broadcast((const_cast<char*>(micro_inputfile_name.c_str())), length, 0);

        // start with actual reading
        INPUT::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

        Teuchos::RCP<DRT::Discretization> dis_micro =
            Teuchos::rcp(new DRT::Discretization(micro_dis_name, micro_reader.Comm()));

        // replace standard dofset inside micro discretization by independent dofset
        // to avoid inconsistent dof numbering in non-nested parallel settings with more than one
        // micro discretization
        if (problem.GetCommunicators()->NpType() ==
            CORE::COMM::NestedParallelismType::no_nested_parallelism)
          dis_micro->ReplaceDofSet(Teuchos::rcp(new DRT::IndependentDofSet()));

        // create discretization writer - in constructor set into and owned by corresponding
        // discret
        dis_micro->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(dis_micro,
            micro_problem->OutputControlFile(), micro_problem->SpatialApproximationType())));

        micro_problem->AddDis(micro_dis_name, dis_micro);

        ReadParameter(*micro_problem, micro_reader);

        // read materials of microscale
        // CAUTION: materials for microscale cannot be read until
        // micro_reader is activated, since else materials will again be
        // read from macroscale inputfile. Besides, materials MUST be read
        // before elements are read since elements establish a connection
        // to the corresponding material! Thus do not change position of
        // function calls!
        problem.Materials()->SetReadFromProblem(microdisnum);

        ReadMaterials(*micro_problem, micro_reader);

        INPUT::MeshReader micromeshreader(micro_reader, "--NODE COORDS");

        if (micro_dis_name == "structure")
        {
          micromeshreader.AddElementReader(
              INPUT::ElementReader(dis_micro, micro_reader, "--STRUCTURE ELEMENTS"));
        }
        else
          micromeshreader.AddElementReader(
              INPUT::ElementReader(dis_micro, micro_reader, "--TRANSPORT ELEMENTS"));

        micromeshreader.ReadAndPartition();


        {
          CORE::UTILS::FunctionManager function_manager;
          GlobalLegacyModuleCallbacks().AttachFunctionDefinitions(function_manager);
          function_manager.ReadInput(micro_reader);
          micro_problem->SetFunctionManager(std::move(function_manager));
        }

        ReadResult(*micro_problem, micro_reader);
        ReadConditions(*micro_problem, micro_reader);

        // At this point, everything for the microscale is read,
        // subsequent reading is only for macroscale
        dis_micro->FillComplete();

        // broadcast restart information
        int restart_step = problem.Restart();
        subgroupcomm->Broadcast(&restart_step, 1, 0);
        problem.SetRestartStep(restart_step);

        // set the problem number from which to call materials again to zero
        // (i.e. macro problem), cf. MAT::Material::Factory!
        problem.Materials()->ResetReadFromProblem();
      }
    }
    problem.Materials()->ResetReadFromProblem();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadMicrofieldsNPsupport(GLOBAL::Problem& problem)
{
  Teuchos::RCP<Epetra_Comm> lcomm = problem.GetCommunicators()->LocalComm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem.GetCommunicators()->GlobalComm();

  // receive number of procs that have micro material
  int nummicromat = 0;
  gcomm->Broadcast(&nummicromat, 1, 0);

  // prepare the supporting procs for a splitting of gcomm

  // groups should be equally sized
  // in a first step every macro proc that needs support gets procpergroup supporting procs
  int procpergroup = int(floor((lcomm->NumProc()) / nummicromat));
  std::vector<int> supgrouplayout(nummicromat, procpergroup);
  // remaining procs are added to the groups in the beginning
  int remainingProcs = lcomm->NumProc() - procpergroup * nummicromat;
  for (int k = 0; k < remainingProcs; ++k)
  {
    supgrouplayout[k]++;
  }

  // secondly: colors are distributed
  // color starts with 0 and is incremented for each group
  int color = -1;
  int gsum = 0;
  do
  {
    color++;
    gsum += supgrouplayout[color];
  } while (gsum <= lcomm->MyPID());

  // do the splitting of the communicator
  MPI_Comm mpi_local_comm;
  MPI_Comm_split((Teuchos::rcp_dynamic_cast<Epetra_MpiComm>(gcomm, true)->GetMpiComm()), color,
      gcomm->MyPID(), &mpi_local_comm);

  // create the sub communicator that includes one macro proc and some supporting procs
  Teuchos::RCP<Epetra_Comm> subgroupcomm = Teuchos::rcp(new Epetra_MpiComm(mpi_local_comm));
  problem.GetCommunicators()->SetSubComm(subgroupcomm);

  // number of micro problems for this sub group
  int microcount = 0;
  subgroupcomm->Broadcast(&microcount, 1, 0);

  for (int n = 0; n < microcount; n++)
  {
    // broadcast microdis number
    int microdisnum = -1;
    subgroupcomm->Broadcast(&microdisnum, 1, 0);

    GLOBAL::Problem* micro_problem = GLOBAL::Problem::Instance(microdisnum);

    // broadcast micro input file name
    int length = -1;
    std::string micro_inputfile_name;
    subgroupcomm->Broadcast(&length, 1, 0);
    micro_inputfile_name.resize(length);
    subgroupcomm->Broadcast((const_cast<char*>(micro_inputfile_name.c_str())), length, 0);

    // start with actual reading
    INPUT::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

    Teuchos::RCP<DRT::Discretization> structdis_micro =
        Teuchos::rcp(new DRT::Discretization("structure", micro_reader.Comm()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis_micro->SetWriter(Teuchos::rcp(new IO::DiscretizationWriter(structdis_micro,
        micro_problem->OutputControlFile(), micro_problem->SpatialApproximationType())));

    micro_problem->AddDis("structure", structdis_micro);

    ReadParameter(*micro_problem, micro_reader);

    // read materials of microscale
    // CAUTION: materials for microscale cannot be read until
    // micro_reader is activated, since else materials will again be
    // read from macroscale inputfile. Besides, materials MUST be read
    // before elements are read since elements establish a connection
    // to the corresponding material! Thus do not change position of
    // function calls!
    problem.Materials()->SetReadFromProblem(microdisnum);

    ReadMaterials(*micro_problem, micro_reader);

    INPUT::MeshReader micromeshreader(micro_reader, "--NODE COORDS");
    micromeshreader.AddElementReader(
        INPUT::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS"));
    micromeshreader.ReadAndPartition();

    // read conditions of microscale
    // -> note that no time curves and spatial functions can be read!

    ReadConditions(*micro_problem, micro_reader);

    // At this point, everything for the microscale is read,
    // subsequent reading is only for macroscale
    structdis_micro->FillComplete();

    // broadcast restart information
    int restart_step = problem.Restart();
    subgroupcomm->Broadcast(&restart_step, 1, 0);
    problem.SetRestartStep(restart_step);

    // set the problem number from which to call materials again to zero
    // (i.e. macro problem), cf. MAT::Material::Factory!
    problem.Materials()->ResetReadFromProblem();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadParameter(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList("DAT FILE"));

  reader.ReadSection("--DISCRETISATION", *list);
  reader.ReadSection("--PROBLEM SIZE", *list);
  reader.ReadSection("--PROBLEM TYP", *list);
  reader.ReadSection("--BINNING STRATEGY", *list);
  reader.ReadSection("--BOUNDINGVOLUME STRATEGY", *list);
  reader.ReadSection("--IO", *list);
  reader.ReadSection("--IO/EVERY ITERATION", *list);
  reader.ReadSection("--IO/MONITOR STRUCTURE DBC", *list);
  reader.ReadSection("--IO/RUNTIME VTK OUTPUT", *list);
  reader.ReadSection("--IO/RUNTIME VTK OUTPUT/FLUID", *list);
  reader.ReadSection("--IO/RUNTIME VTK OUTPUT/STRUCTURE", *list);
  reader.ReadSection("--IO/RUNTIME VTK OUTPUT/BEAMS", *list);
  reader.ReadSection("--IO/RUNTIME VTP OUTPUT STRUCTURE", *list);
  reader.ReadSection("--STRUCTURAL DYNAMIC", *list);
  reader.ReadSection("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadSection("--STRUCTURAL DYNAMIC/GENALPHA", *list);
  reader.ReadSection("--STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadSection("--STRUCTURAL DYNAMIC/GEMM", *list);
  reader.ReadSection("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY/JOINT EXPLICIT", *list);
  reader.ReadSection("--MORTAR COUPLING", *list);
  reader.ReadSection("--MORTAR COUPLING/PARALLEL REDISTRIBUTION", *list);
  reader.ReadSection("--CONTACT DYNAMIC", *list);
  reader.ReadSection("--CONTACT DYNAMIC/AUGMENTED", *list);
  reader.ReadSection("--CONTACT DYNAMIC/AUGMENTED/COMBO", *list);
  reader.ReadSection("--CONTACT DYNAMIC/AUGMENTED/STEEPESTASCENT", *list);
  reader.ReadSection("--CONTACT DYNAMIC/AUGMENTED/LAGRANGE_MULTIPLIER_FUNCTION", *list);
  reader.ReadSection("--CONTACT DYNAMIC/AUGMENTED/PLOT", *list);
  reader.ReadSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING", *list);
  reader.ReadSection(
      "--CARDIOVASCULAR 0D-STRUCTURE COUPLING/SYS-PUL CIRCULATION PARAMETERS", *list);
  reader.ReadSection("--CARDIOVASCULAR 0D-STRUCTURE COUPLING/RESPIRATORY PARAMETERS", *list);
  reader.ReadSection("--BROWNIAN DYNAMICS", *list);
  reader.ReadSection("--BEAM INTERACTION", *list);
  reader.ReadSection("--BEAM INTERACTION/SPHERE BEAM LINK", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO BEAM CONTACT", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SPHERE CONTACT", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE CONTACT", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SOLID SURFACE/RUNTIME VTK OUTPUT", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING", *list);
  reader.ReadSection("--BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING/RUNTIME VTK OUTPUT", *list);
  reader.ReadSection("--BEAM INTERACTION/CROSSLINKING", *list);
  reader.ReadSection("--THERMAL DYNAMIC", *list);
  reader.ReadSection("--THERMAL DYNAMIC/GENALPHA", *list);
  reader.ReadSection("--THERMAL DYNAMIC/ONESTEPTHETA", *list);
  reader.ReadSection("--TSI DYNAMIC", *list);
  reader.ReadSection("--TSI DYNAMIC/MONOLITHIC", *list);
  reader.ReadSection("--TSI DYNAMIC/PARTITIONED", *list);
  reader.ReadSection("--TSI CONTACT", *list);
  reader.ReadSection("--POROELASTICITY DYNAMIC", *list);
  reader.ReadSection("--POROSCATRA CONTROL", *list);
  reader.ReadSection("--POROFLUIDMULTIPHASE DYNAMIC", *list);
  reader.ReadSection("--POROFLUIDMULTIPHASE DYNAMIC/ARTERY COUPLING", *list);
  reader.ReadSection("--POROMULTIPHASE DYNAMIC", *list);
  reader.ReadSection("--POROMULTIPHASE DYNAMIC/PARTITIONED", *list);
  reader.ReadSection("--POROMULTIPHASE DYNAMIC/MONOLITHIC", *list);
  reader.ReadSection("--POROMULTIPHASESCATRA DYNAMIC", *list);
  reader.ReadSection("--POROMULTIPHASESCATRA DYNAMIC/PARTITIONED", *list);
  reader.ReadSection("--POROMULTIPHASESCATRA DYNAMIC/MONOLITHIC", *list);
  reader.ReadSection("--ELASTO HYDRO DYNAMIC", *list);
  reader.ReadSection("--ELASTO HYDRO DYNAMIC/PARTITIONED", *list);
  reader.ReadSection("--ELASTO HYDRO DYNAMIC/MONOLITHIC", *list);
  reader.ReadSection("--SSI CONTROL", *list);
  reader.ReadSection("--SSI CONTROL/ELCH", *list);
  reader.ReadSection("--SSI CONTROL/MANIFOLD", *list);
  reader.ReadSection("--SSI CONTROL/MONOLITHIC", *list);
  reader.ReadSection("--SSI CONTROL/PARTITIONED", *list);
  reader.ReadSection("--SSTI CONTROL", *list);
  reader.ReadSection("--SSTI CONTROL/MONOLITHIC", *list);
  reader.ReadSection("--SSTI CONTROL/THERMO", *list);
  reader.ReadSection("--FLUID DYNAMIC", *list);
  reader.ReadSection("--FLUID DYNAMIC/RESIDUAL-BASED STABILIZATION", *list);
  reader.ReadSection("--FLUID DYNAMIC/EDGE-BASED STABILIZATION", *list);
  reader.ReadSection("--FLUID DYNAMIC/POROUS-FLOW STABILIZATION", *list);
  reader.ReadSection("--FLUID DYNAMIC/TURBULENCE MODEL", *list);
  reader.ReadSection("--FLUID DYNAMIC/SUBGRID VISCOSITY", *list);
  reader.ReadSection("--FLUID DYNAMIC/WALL MODEL", *list);
  reader.ReadSection("--FLUID DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadSection("--FLUID DYNAMIC/MULTIFRACTAL SUBGRID SCALES", *list);
  reader.ReadSection("--FLUID DYNAMIC/TURBULENT INFLOW", *list);
  reader.ReadSection("--FLUID DYNAMIC/NONLINEAR SOLVER TOLERANCES", *list);
  reader.ReadSection("--LUBRICATION DYNAMIC", *list);
  reader.ReadSection("--SCALAR TRANSPORT DYNAMIC", *list);
  reader.ReadSection("--SCALAR TRANSPORT DYNAMIC/NONLINEAR", *list);
  reader.ReadSection("--SCALAR TRANSPORT DYNAMIC/STABILIZATION", *list);
  reader.ReadSection("--SCALAR TRANSPORT DYNAMIC/S2I COUPLING", *list);
  reader.ReadSection("--SCALAR TRANSPORT DYNAMIC/ARTERY COUPLING", *list);
  reader.ReadSection("--STI DYNAMIC", *list);
  reader.ReadSection("--STI DYNAMIC/MONOLITHIC", *list);
  reader.ReadSection("--STI DYNAMIC/PARTITIONED", *list);
  reader.ReadSection("--FS3I DYNAMIC", *list);
  reader.ReadSection("--FS3I DYNAMIC/PARTITIONED", *list);
  reader.ReadSection("--FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION", *list);
  reader.ReadSection("--FS3I DYNAMIC/AC", *list);
  reader.ReadSection("--ALE DYNAMIC", *list);
  reader.ReadSection("--FSI DYNAMIC", *list);
  reader.ReadSection("--FSI DYNAMIC/CONSTRAINT", *list);
  reader.ReadSection("--FSI DYNAMIC/MONOLITHIC SOLVER", *list);
  reader.ReadSection("--FSI DYNAMIC/PARTITIONED SOLVER", *list);
  reader.ReadSection("--FSI DYNAMIC/TIMEADAPTIVITY", *list);
  reader.ReadSection("--FLUID BEAM INTERACTION", *list);
  reader.ReadSection("--FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING", *list);
  reader.ReadSection("--FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING/RUNTIME VTK OUTPUT", *list);
  reader.ReadSection("--IMMERSED METHOD", *list);
  reader.ReadSection("--IMMERSED METHOD/PARTITIONED SOLVER", *list);
  reader.ReadSection("--FPSI DYNAMIC", *list);
  reader.ReadSection("--ARTERIAL DYNAMIC", *list);
  reader.ReadSection("--REDUCED DIMENSIONAL AIRWAYS DYNAMIC", *list);
  reader.ReadSection("--COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC", *list);
  reader.ReadSection("--SEARCH TREE", *list);
  reader.ReadSection("--XFEM GENERAL", *list);
  reader.ReadSection("--CUT GENERAL", *list);
  reader.ReadSection("--XFLUID DYNAMIC", *list);
  reader.ReadSection("--XFLUID DYNAMIC/GENERAL", *list);
  reader.ReadSection("--XFLUID DYNAMIC/STABILIZATION", *list);
  reader.ReadSection("--XFLUID DYNAMIC/XFPSI MONOLITHIC", *list);
  reader.ReadSection("--LOMA CONTROL", *list);
  reader.ReadSection("--ELCH CONTROL", *list);
  reader.ReadSection("--ELCH CONTROL/DIFFCOND", *list);
  reader.ReadSection("--ELCH CONTROL/SCL", *list);
  reader.ReadSection("--BIOFILM CONTROL", *list);
  reader.ReadSection("--PARTICLE DYNAMIC", *list);
  reader.ReadSection("--PARTICLE DYNAMIC/INITIAL AND BOUNDARY CONDITIONS", *list);
  reader.ReadSection("--PARTICLE DYNAMIC/SPH", *list);
  reader.ReadSection("--PARTICLE DYNAMIC/DEM", *list);
  reader.ReadSection("--PASI DYNAMIC", *list);
  reader.ReadSection("--LEVEL-SET CONTROL", *list);
  reader.ReadSection("--LEVEL-SET CONTROL/REINITIALIZATION", *list);
  reader.ReadSection("--WEAR", *list);
  reader.ReadSection("--BEAM CONTACT", *list);
  reader.ReadSection("--BEAM CONTACT/RUNTIME VTK OUTPUT", *list);
  reader.ReadSection("--BEAM POTENTIAL", *list);
  reader.ReadSection("--BEAM POTENTIAL/RUNTIME VTK OUTPUT", *list);
  reader.ReadSection("--SEMI-SMOOTH PLASTICITY", *list);
  reader.ReadSection("--ELECTROMAGNETIC DYNAMIC", *list);
  reader.ReadSection("--VOLMORTAR COUPLING", *list);
  reader.ReadSection("--CARDIAC MONODOMAIN CONTROL", *list);
  reader.ReadSection("--MOR", *list);
  reader.ReadSection("--MESH PARTITIONING", *list);

  reader.ReadSection("--STRUCT NOX", *list);
  reader.ReadSection("--STRUCT NOX/Direction", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton/Modified", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Newton/Linear Solver", *list);
  reader.ReadSection("--STRUCT NOX/Direction/Steepest Descent", *list);
  reader.ReadSection("--STRUCT NOX/Line Search", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Full Step", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Backtrack", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/Polynomial", *list);
  reader.ReadSection("--STRUCT NOX/Line Search/More'-Thuente", *list);
  reader.ReadSection("--STRUCT NOX/Pseudo Transient", *list);
  reader.ReadSection("--STRUCT NOX/Trust Region", *list);
  reader.ReadSection("--STRUCT NOX/Printing", *list);
  reader.ReadSection("--STRUCT NOX/Status Test", *list);
  reader.ReadSection("--STRUCT NOX/Solver Options", *list);

  // read in solver sections
  // Note: the maximum number of solver blocks in dat files is hardwired here.
  // If you change this do not forget to edit the corresponding parts in
  // validparameters.cpp, too!
  for (int i = 1; i < 10; i++)
  {
    std::stringstream ss;
    ss << "--SOLVER " << i;
    reader.ReadSection(ss.str(), *list);

    // adapt path of XML file if necessary
    Teuchos::ParameterList& sublist = list->sublist(ss.str().substr(2));
    std::vector<std::string> listOfFileNameParameters = {"AMGNXN_XML_FILE", "MUELU_XML_FILE"};

    for (auto& filenameParameter : listOfFileNameParameters)
    {
      auto* xml_filename = sublist.getPtr<std::string>(filenameParameter);
      if (xml_filename != nullptr and *xml_filename != "none")
      {
        // make path relative to input file path if it is not an absolute path
        if ((*xml_filename)[0] != '/')
        {
          std::string filename = reader.MyInputfileName();
          std::string::size_type pos = filename.rfind('/');
          if (pos != std::string::npos)
          {
            std::string tmp = filename.substr(0, pos + 1);
            xml_filename->insert(xml_filename->begin(), tmp.begin(), tmp.end());
          }
        }
      }
    }
  }

  reader.ReadSection("--UMFPACK SOLVER", *list);


  // read in STRUCT NOX/Status Test and modify the xml file name, if there
  // is one.
  if (list->sublist("STRUCT NOX").sublist("Status Test").isParameter("XML File"))
  {
    // adapt path of XML file if necessary
    Teuchos::ParameterList& sublist = list->sublist("STRUCT NOX").sublist("Status Test");
    auto* statustest_xmlfile = sublist.getPtr<std::string>("XML File");
    // make path relative to input file path if it is not an absolute path
    if (((*statustest_xmlfile)[0] != '/') and ((*statustest_xmlfile) != "none"))
    {
      std::string filename = reader.MyInputfileName();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string tmp = filename.substr(0, pos + 1);
        statustest_xmlfile->insert(statustest_xmlfile->begin(), tmp.begin(), tmp.end());
      }
    }
  }  // STRUCT NOX/Status Test

  // check for invalid parameters
  problem.setParameterList(list);

  //---------------------------------------------------------------------
  // Now we have successfully read the whole input file. It's time to access some data

  // 1) get the problem type
  const Teuchos::ParameterList& type = problem.ProblemTypeParams();
  problem.SetProblemType(INPUT::IntegralValue<GLOBAL::ProblemType>(type, "PROBLEMTYP"));

  // 2) get the spatial approximation type
  problem.SetSpatialApproximationType(
      INPUT::IntegralValue<CORE::FE::ShapeFunctionType>(type, "SHAPEFCT"));

  int restart_step = problem.Restart();
  // 3) do the restart business with the four options we support (partially)
  if (restart_step == 0)
  {
    // no restart flag on the command line, so check the restart flag from the input file
    restart_step = type.get<int>("RESTART");
    problem.SetRestartStep(restart_step);
  }
  else  // SetRestartStep() has been called before!
  {
    // There is a non-zero restart flag on the command line, so we ignore the input file.
    // The RESTART flag in the input file should be zero or have the same value!
    const int restartflaginfile = type.get<int>("RESTART");
    if ((restartflaginfile > 0) and (restartflaginfile != restart_step))
      dserror("Restart flags in input file and command line are non-zero and different!");
  }

  // Set restart time based on walltime
  const double restartinterval = problem.IOParams().get<double>("RESTARTWALLTIMEINTERVAL");
  const int restartevry = problem.IOParams().get<int>("RESTARTEVRY");
  problem.RestartManager()->SetupRestartManager(restartinterval, restartevry);

  // 4) set random seed
  // time is in seconds, therefore we add the global processor id to obtain a unique seed on each
  // proc
  {
    int rs = type.get<int>("RANDSEED");
    if (rs < 0)
      rs = static_cast<int>(time(nullptr)) +
           42 * GLOBAL::Problem::Instance(0)->GetCommunicators()->GlobalComm()->MyPID();

    srand((unsigned int)rs);  // Set random seed for stdlibrary. This is deprecated, as it does not
    // produce random numbers on some platforms!
    problem.Random()->SetRandSeed((unsigned int)rs);  // Use this instead.
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadMaterials(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  // create list of known materials
  Teuchos::RCP<std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>> vm = INPUT::ValidMaterials();
  std::vector<Teuchos::RCP<INPUT::MaterialDefinition>>& matlist = *vm;

  // test for each material definition (input file --MATERIALS section)
  // and store in #matmap_
  for (auto& mat : matlist)
  {
    // read material from DAT file of type #matlist[m]
    mat->Read(problem, reader, problem.Materials());
  }

  // check if every material was identified
  const std::string material_section = "--MATERIALS";
  std::vector<const char*> section = reader.Section(material_section);
  if (section.size() > 0)
  {
    for (auto& section_i : section)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(section_i));

      std::string mat;
      std::string number;
      std::string name;
      (*condline) >> mat >> number >> name;
      if ((not(*condline)) or (mat != "MAT"))
        dserror("invalid material line in '%s'", name.c_str());

      // extract material ID
      int matid = -1;
      {
        char* ptr;
        matid = static_cast<int>(strtol(number.c_str(), &ptr, 10));
        if (ptr == number.c_str())
          dserror("failed to read material object number '%s'", number.c_str());
      }

      // processed?
      if (problem.Materials()->Find(matid) == -1)
        dserror("Material 'MAT %d' with name '%s' could not be identified", matid, name.c_str());
    }
  }

  // make fast access parameters
  problem.Materials()->MakeParameters();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadContactConstitutiveLaws(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  // create list of known contact constitutive laws
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>> vm =
      INPUT::ValidContactConstitutiveLaws();
  std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist = *vm;

  // test for each contact constitutive law definition (input file --CONTACT CONSTITUTIVE LAWS
  // section) and store it
  for (auto& m : coconstlawlist)
  {
    // read contact constitutive law from DAT file of type
    m->Read(problem, reader, problem.ContactConstitutiveLaws());
  }

  // check if every contact constitutive law was identified
  const std::string contact_const_laws = "--CONTACT CONSTITUTIVE LAWS";
  std::vector<const char*> section = reader.Section(contact_const_laws);
  if (section.size() > 0)
  {
    for (auto& section_i : section)
    {
      Teuchos::RCP<std::stringstream> condline = Teuchos::rcp(new std::stringstream(section_i));

      std::string coconstlaw;
      std::string number;
      std::string name;
      (*condline) >> coconstlaw >> number >> name;
      if ((not(*condline)) or (coconstlaw != "LAW"))
        dserror("invalid contact constitutive law line in '%s'", name.c_str());

      // extract contact constitutive law ID
      int coconstlawid = -1;
      {
        char* ptr;
        coconstlawid = static_cast<int>(strtol(number.c_str(), &ptr, 10));
        if (ptr == number.c_str())
          dserror("failed to read contact constitutive law object number '%s'", number.c_str());
      }

      // processed?
      if (problem.ContactConstitutiveLaws()->Find(coconstlawid) == -1)
        dserror("Contact constitutive law 'LAW %d' with name '%s' could not be identified",
            coconstlawid, name.c_str());
    }
  }

  // make fast access parameters
  problem.ContactConstitutiveLaws()->MakeParameters();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadCloningMaterialMap(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  INPUT::Lines lines = DRT::UTILS::ValidCloningMaterialMapLines();

  // perform the actual reading and extract the input parameters
  std::vector<INPUT::LineDefinition> input_line_vec = lines.Read(reader);
  for (const auto& input_line : input_line_vec)
  {
    // extract what was read from the input file
    std::string src_field;
    input_line.ExtractString("SRC_FIELD", src_field);
    int src_matid(-1);
    input_line.ExtractInt("SRC_MAT", src_matid);
    std::string tar_field;
    input_line.ExtractString("TAR_FIELD", tar_field);
    int tar_matid(-1);
    input_line.ExtractInt("TAR_MAT", tar_matid);

    // create the key pair
    std::pair<std::string, std::string> fields(src_field, tar_field);

    // enter the material pairing into the map
    std::pair<int, int> matmap(src_matid, tar_matid);
    problem.CloningMaterialMap()[fields].insert(matmap);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadResult(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  problem.GetResultTestManager().ReadInput(reader);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadConditions(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  Teuchos::Time time("", true);
  if (reader.Comm()->MyPID() == 0)
  {
    IO::cout << "Read/generate conditions                          in....";
    IO::cout.flush();
  }

  //--------------------------------------------- read generic node sets
  // read design nodes <-> nodes
  std::vector<std::vector<int>> dnode_fenode;
  reader.ReadDesign("DNODE", dnode_fenode);

  // read design lines <-> nodes
  std::vector<std::vector<int>> dline_fenode;
  reader.ReadDesign("DLINE", dline_fenode);

  // read design surfaces <-> nodes
  std::vector<std::vector<int>> dsurf_fenode;
  reader.ReadDesign("DSURF", dsurf_fenode);

  // read design volumes <-> nodes
  std::vector<std::vector<int>> dvol_fenode;
  reader.ReadDesign("DVOL", dvol_fenode);

  // check for meshfree discretisation to add node set topologies
  std::vector<std::vector<std::vector<int>>*> nodeset(4);
  nodeset[0] = &dnode_fenode;
  nodeset[1] = &dline_fenode;
  nodeset[2] = &dsurf_fenode;
  nodeset[3] = &dvol_fenode;

  // create list of known conditions
  Teuchos::RCP<std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>> vc = INPUT::ValidConditions();
  std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist = *vc;

  // test for each condition definition (input file condition section)
  // - read all conditions that match the definition
  // - add the nodal clouds to the conditions
  // - add the conditions to the appropriate discretizations
  //
  // Note that this will reset (un-FillComplete) the discretizations.
  for (auto& condition : condlist)
  {
    std::multimap<int, Teuchos::RCP<DRT::Condition>> cond;

    // read conditions from dat file
    condition->Read(problem, reader, cond);

    // add nodes to conditions
    std::multimap<int, Teuchos::RCP<DRT::Condition>>::const_iterator curr;
    for (curr = cond.begin(); curr != cond.end(); ++curr)
    {
      switch (curr->second->GType())
      {
        case DRT::Condition::Point:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dnode_fenode.size())
          {
            dserror(
                "DPoint %d not in range [0:%d[\n"
                "DPoint condition on non existent DPoint?",
                curr->first, dnode_fenode.size());
          }
          curr->second->Add("Node Ids", dnode_fenode[curr->first]);
          break;
        case DRT::Condition::Line:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dline_fenode.size())
          {
            dserror(
                "DLine %d not in range [0:%d[\n"
                "DLine condition on non existent DLine?",
                curr->first, dline_fenode.size());
          }
          curr->second->Add("Node Ids", dline_fenode[curr->first]);
          break;
        case DRT::Condition::Surface:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dsurf_fenode.size())
          {
            dserror(
                "DSurface %d not in range [0:%d[\n"
                "DSurface condition on non existent DSurface?",
                curr->first, dsurf_fenode.size());
          }
          curr->second->Add("Node Ids", dsurf_fenode[curr->first]);
          break;
        case DRT::Condition::Volume:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dvol_fenode.size())
          {
            dserror(
                "DVolume %d not in range [0:%d[\n"
                "DVolume condition on non existent DVolume?",
                curr->first, dvol_fenode.size());
          }
          curr->second->Add("Node Ids", dvol_fenode[curr->first]);
          break;
        default:
          dserror("geometry type unspecified");
          break;
      }

      // Iterate through all discretizations and sort the appropriate condition
      // into the correct discretization it applies to

      for (const auto& [name, dis] : problem.DiscretizationRange())
      {
        const std::vector<int>* nodes = curr->second->Nodes();
        if (nodes->size() == 0)
          dserror("%s condition %d has no nodal cloud", condition->Description().c_str(),
              curr->second->Id());

        int foundit = 0;
        for (int node : *nodes)
        {
          foundit = dis->HaveGlobalNode(node);
          if (foundit) break;
        }
        int found = 0;
        dis->Comm().SumAll(&foundit, &found, 1);
        if (found)
        {
          // Insert a copy since we might insert the same condition in many discretizations.
          dis->SetCondition(condition->Name(), Teuchos::rcp(new DRT::Condition(*curr->second)));
        }
      }
    }
  }

  if (reader.Comm()->MyPID() == 0)
  {
    std::cout << time.totalElapsedTime(true) << " secs\n";
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadKnots(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  // get information on the spatial approximation --- we only read knots
  // in the nurbs case
  CORE::FE::ShapeFunctionType distype = problem.SpatialApproximationType();

  // get problem dimension
  int dim = problem.NDim();

  // Iterate through all discretizations and sort the appropriate condition
  // into the correct discretization it applies to

  for (const auto& [name, dis] : problem.DiscretizationRange())
  {
    if (distype == CORE::FE::ShapeFunctionType::nurbs)
    {
      // cast discretisation to nurbs variant to be able
      // to add the knotvector
      auto* nurbsdis = dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

      if (nurbsdis == nullptr)
        dserror("Discretization %s is not a NurbsDiscretization! Panic.", dis->Name().c_str());

      // define an empty knot vector object
      Teuchos::RCP<DRT::NURBS::Knotvector> disknots = Teuchos::null;

      // read the knotvector data from the input
      reader.ReadKnots(dim, dis->Name(), disknots);

      if (disknots == Teuchos::null)
      {
        dserror("Knotvector read failed in Nurbs discretisation\n");
      }

      // make sure atdis is fillcompleted, to be able to call
      // ElementRowMap() on it
      // do not initialize elements, since this would require knot
      // vector values
      if (!dis->Filled())
      {
        dis->FillComplete(false, false, false);
      }

      // the smallest gid in the discretisation determines the access
      // pattern via the element offset
      int smallest_gid_in_dis = dis->ElementRowMap()->MinAllGID();

      // consistency checks
      disknots->FinishKnots(smallest_gid_in_dis);

      // add knots to discretisation
      nurbsdis->SetKnotVector(disknots);
    }
  }  // loop fields
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void GLOBAL::ReadParticles(GLOBAL::Problem& problem, INPUT::DatFileReader& reader)
{
  // no need to read in particles in case of restart
  if (problem.Restart()) return;

  // the basic particle reader
  INPUT::ParticleReader particlereader(reader, "--PARTICLES");

  // do the actual reading of particles
  particlereader.Read(problem.Particles());
}


BACI_NAMESPACE_CLOSE
