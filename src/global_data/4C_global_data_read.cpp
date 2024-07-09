/*----------------------------------------------------------------------*/
/*! \file

\brief Utilities to read and fill global data

\level 1

    */
/*----------------------------------------------------------------------*/

#include "4C_global_data_read.hpp"

#include "4C_comm_utils.hpp"
#include "4C_contact_constitutivelaw_bundle.hpp"
#include "4C_contact_constitutivelaw_constitutivelaw_definition.hpp"
#include "4C_fem_condition_definition.hpp"
#include "4C_fem_discretization_hdg.hpp"
#include "4C_fem_dofset_independent.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fem_nurbs_discretization.hpp"
#include "4C_global_legacy_module.hpp"
#include "4C_inpar_validconditions.hpp"
#include "4C_inpar_validcontactconstitutivelaw.hpp"
#include "4C_inpar_validmaterials.hpp"
#include "4C_io.hpp"
#include "4C_io_dat_file_utils.hpp"
#include "4C_io_elementreader.hpp"
#include "4C_io_geometry_type.hpp"
#include "4C_io_inputreader.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_io_meshreader.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_elchphase.hpp"
#include "4C_mat_materialdefinition.hpp"
#include "4C_mat_micromaterial.hpp"
#include "4C_mat_newman_multiscale.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_scatra_multiscale.hpp"
#include "4C_particle_engine_particlereader.hpp"
#include "4C_rebalance_graph_based.hpp"
#include "4C_xfem_discretization.hpp"
#include "4C_xfem_discretization_utils.hpp"

FOUR_C_NAMESPACE_OPEN

void Global::ReadFields(
    Global::Problem& problem, Core::IO::DatFileReader& reader, const bool read_mesh)
{
  Teuchos::RCP<Core::FE::Discretization> structdis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> fluiddis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> xfluiddis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> aledis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> structaledis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> thermdis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> lubricationdis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> scatradis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> scatra_micro_dis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> cellscatradis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> fluidscatradis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> structscatradis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> artscatradis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> arterydis = Teuchos::null;  //_1D_ARTERY_
  Teuchos::RCP<Core::FE::Discretization> airwaydis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> optidis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> porofluiddis = Teuchos::null;  // fpsi, poroelast
  Teuchos::RCP<Core::FE::Discretization> elemagdis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> celldis = Teuchos::null;
  Teuchos::RCP<Core::FE::Discretization> pboxdis = Teuchos::null;

  // decide which kind of spatial representation is required
  const Core::FE::ShapeFunctionType distype = problem.spatial_approximation_type();
  auto output_control = problem.output_control_file();

  // the basic mesh reader. now add desired node and element readers to it!
  Core::IO::MeshReader meshreader(reader, "--NODE COORDS",
      {.mesh_paritioning_parameters = Problem::instance()->mesh_partitioning_params(),
          .geometric_search_parameters = Problem::instance()->geometric_search_params(),
          .io_parameters = Problem::instance()->io_params()});

  switch (problem.get_problem_type())
  {
    case Core::ProblemType::fsi:
    case Core::ProblemType::fsi_redmodels:
    case Core::ProblemType::fsi_lung:
    {
      if (distype == Core::FE::ShapeFunctionType::nurbs)
      {
        structdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
            "structure", reader.get_comm(), problem.n_dim()));
        fluiddis = Teuchos::rcp(
            new Core::FE::Nurbs::NurbsDiscretization("fluid", reader.get_comm(), problem.n_dim()));
        aledis = Teuchos::rcp(
            new Core::FE::Nurbs::NurbsDiscretization("ale", reader.get_comm(), problem.n_dim()));
      }
      else if (Core::UTILS::IntegralValue<int>(
                   problem.fluid_dynamic_params().sublist("WALL MODEL"), "X_WALL"))
      {
        structdis = Teuchos::rcp(
            new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
        fluiddis = Teuchos::rcp(
            new XFEM::DiscretizationXWall("fluid", reader.get_comm(), problem.n_dim()));
        aledis =
            Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
      }
      else
      {
        structdis = Teuchos::rcp(
            new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
        if (Core::UTILS::IntegralValue<bool>(
                problem.x_fluid_dynamic_params().sublist("GENERAL"), "XFLUIDFLUID"))
          xfluiddis = Teuchos::rcp(
              new XFEM::DiscretizationXFEM("xfluid", reader.get_comm(), problem.n_dim()));
        aledis =
            Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      if (xfluiddis != Teuchos::null)
        xfluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(xfluiddis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("fluid", fluiddis);
      if (xfluiddis != Teuchos::null) problem.add_dis("xfluid", xfluiddis);
      problem.add_dis("ale", aledis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      if (xfluiddis != Teuchos::null)
      {
        meshreader.add_element_reader(
            Core::IO::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS"));
      }
      else
        meshreader.add_element_reader(
            Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case Core::ProblemType::gas_fsi:
    case Core::ProblemType::ac_fsi:
    case Core::ProblemType::thermo_fsi:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          FOUR_C_THROW("Nurbs discretization not possible for fs3i!");
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(
              new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
          fluiddis = Teuchos::rcp(
              new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
          aledis =
              Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
          fluidscatradis = Teuchos::rcp(
              new Core::FE::Discretization("scatra1", reader.get_comm(), problem.n_dim()));
          structscatradis = Teuchos::rcp(
              new Core::FE::Discretization("scatra2", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));
      fluidscatradis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(fluidscatradis, output_control, distype)));
      structscatradis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(structscatradis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("fluid", fluiddis);
      problem.add_dis("ale", aledis);
      problem.add_dis("scatra1", fluidscatradis);
      problem.add_dis("scatra2", structscatradis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(fluidscatradis, reader, "--TRANSPORT ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(structscatradis, reader, "--TRANSPORT2 ELEMENTS"));

#ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false, false, false);
#endif

      break;
    }
    case Core::ProblemType::biofilm_fsi:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          FOUR_C_THROW("Nurbs discretization not possible for biofilm problems!");
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(
              new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
          fluiddis = Teuchos::rcp(
              new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
          aledis =
              Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
          structaledis = Teuchos::rcp(
              new Core::FE::Discretization("structale", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));
      structaledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structaledis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("fluid", fluiddis);
      problem.add_dis("ale", aledis);
      problem.add_dis("structale", structaledis);


      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

#ifdef EXTENDEDPARALLELOVERLAP
      structdis->CreateExtendedOverlap(false, false, false);
#endif

      // fluid scatra field
      fluidscatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra1", reader.get_comm(), problem.n_dim()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(fluidscatradis, output_control, distype)));
      problem.add_dis("scatra1", fluidscatradis);

      // structure scatra field
      structscatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra2", reader.get_comm(), problem.n_dim()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(structscatradis, output_control, distype)));
      problem.add_dis("scatra2", structscatradis);

      break;
    }
    case Core::ProblemType::fsi_xfem:
    case Core::ProblemType::fluid_xfem:
    {
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      problem.add_dis("structure", structdis);
      meshreader.add_advanced_reader(structdis, reader, "STRUCTURE",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.structural_dynamic_params(), "GEOMETRY"),
          nullptr);

      if (Core::UTILS::IntegralValue<int>(
              problem.x_fluid_dynamic_params().sublist("GENERAL"), "XFLUIDFLUID"))
      {
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
        fluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
        problem.add_dis("fluid", fluiddis);

        xfluiddis = Teuchos::rcp(
            new XFEM::DiscretizationXFEM("xfluid", reader.get_comm(), problem.n_dim()));
        xfluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(xfluiddis, output_control, distype)));
        problem.add_dis("xfluid", xfluiddis);

        meshreader.add_element_reader(
            Core::IO::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS", "FLUID"));
      }
      else
      {
        fluiddis =
            Teuchos::rcp(new XFEM::DiscretizationXFEM("fluid", reader.get_comm(), problem.n_dim()));
        fluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
        problem.add_dis("fluid", fluiddis);

        meshreader.add_advanced_reader(fluiddis, reader, "FLUID",
            Core::UTILS::IntegralValue<Core::IO::GeometryType>(
                problem.fluid_dynamic_params(), "GEOMETRY"),
            nullptr);
      }

      aledis =
          Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));
      problem.add_dis("ale", aledis);
      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));
      break;
    }
    case Core::ProblemType::fpsi_xfem:
    {
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      fluiddis =
          Teuchos::rcp(new XFEM::DiscretizationXFEM("fluid", reader.get_comm(), problem.n_dim()));
      porofluiddis = Teuchos::rcp(
          new Core::FE::DiscretizationFaces("porofluid", reader.get_comm(), problem.n_dim()));
      aledis =
          Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("porofluid", porofluiddis);
      problem.add_dis("fluid", fluiddis);
      problem.add_dis("ale", aledis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      meshreader.add_advanced_reader(fluiddis, reader, "FLUID",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.fluid_dynamic_params(), "GEOMETRY"),
          nullptr);

      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case Core::ProblemType::ale:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          aledis = Teuchos::rcp(
              new Core::FE::Nurbs::NurbsDiscretization("ale", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          aledis =
              Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.add_dis("ale", aledis);

      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case Core::ProblemType::fluid:
    case Core::ProblemType::fluid_redmodels:
    {
      if (distype == Core::FE::ShapeFunctionType::hdg)
      {
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationHDG("fluid", reader.get_comm(), problem.n_dim()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }
      else if (distype == Core::FE::ShapeFunctionType::nurbs)
      {
        fluiddis = Teuchos::rcp(
            new Core::FE::Nurbs::NurbsDiscretization("fluid", reader.get_comm(), problem.n_dim()));

        // create discretization writer - in constructor set ingto and owned by corresponding
        // discret
        fluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }
      else if (Core::UTILS::IntegralValue<int>(
                   problem.fluid_dynamic_params().sublist("WALL MODEL"), "X_WALL"))
      {
        fluiddis = Teuchos::rcp(
            new XFEM::DiscretizationXWall("fluid", reader.get_comm(), problem.n_dim()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }
      else
      {
        // fluiddis  = Teuchos::rcp(new Core::FE::Discretization("fluid",reader.Comm()));
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));

        // create discretization writer - in constructor set into and owned by corresponding discret
        fluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      }

      problem.add_dis("fluid", fluiddis);

      meshreader.add_advanced_reader(fluiddis, reader, "FLUID",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.fluid_dynamic_params(), "GEOMETRY"),
          nullptr);

      break;
    }
    case Core::ProblemType::lubrication:
    {
      // create empty discretizations
      lubricationdis = Teuchos::rcp(
          new Core::FE::Discretization("lubrication", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      lubricationdis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(lubricationdis, output_control, distype)));

      problem.add_dis("lubrication", lubricationdis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS"));

      break;
    }
    case Core::ProblemType::cardiac_monodomain:
    case Core::ProblemType::scatra:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          fluiddis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "fluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "scatra", reader.get_comm(), problem.n_dim()));
          break;
        }
        case Core::FE::ShapeFunctionType::hdg:
        {
          fluiddis = Teuchos::rcp(
              new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(
              new Core::FE::DiscretizationHDG("scatra", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          fluiddis = Teuchos::rcp(
              new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(
              new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));


      problem.add_dis("fluid", fluiddis);
      problem.add_dis("scatra", scatradis);

      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }
    case Core::ProblemType::sti:
    {
      // safety checks
      if (distype == Core::FE::ShapeFunctionType::nurbs)
        FOUR_C_THROW("Scatra-thermo interaction does not work for nurbs discretizations yet!");

      // create empty discretizations for scalar and thermo fields
      scatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));
      thermdis =
          Teuchos::rcp(new Core::FE::Discretization("thermo", reader.get_comm(), problem.n_dim()));

      // create discretization writers
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));
      thermdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(thermdis, output_control, distype)));

      // add empty discretizations to global problem
      problem.add_dis("scatra", scatradis);
      problem.add_dis("thermo", thermdis);

      // add element reader to node reader
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }
    case Core::ProblemType::fluid_ale:
    case Core::ProblemType::freesurf:
    {
      if (distype == Core::FE::ShapeFunctionType::hdg)
      {
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationHDG("fluid", reader.get_comm(), problem.n_dim()));
        aledis =
            Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
      }
      else if (distype == Core::FE::ShapeFunctionType::nurbs)
      {
        fluiddis = Teuchos::rcp(
            new Core::FE::Nurbs::NurbsDiscretization("fluid", reader.get_comm(), problem.n_dim()));
        aledis = Teuchos::rcp(
            new Core::FE::Nurbs::NurbsDiscretization("ale", reader.get_comm(), problem.n_dim()));
      }
      else if (Core::UTILS::IntegralValue<int>(
                   problem.fluid_dynamic_params().sublist("WALL MODEL"), "X_WALL"))
      {
        fluiddis = Teuchos::rcp(
            new XFEM::DiscretizationXWall("fluid", reader.get_comm(), problem.n_dim()));
        aledis =
            Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
      }
      else
      {
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
        if (Core::UTILS::IntegralValue<bool>(
                problem.x_fluid_dynamic_params().sublist("GENERAL"), "XFLUIDFLUID"))
          xfluiddis = Teuchos::rcp(
              new XFEM::DiscretizationXFEM("xfluid", reader.get_comm(), problem.n_dim()));
        aledis =
            Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
      }


      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      if (xfluiddis != Teuchos::null)
      {
        xfluiddis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(xfluiddis, output_control, distype)));
      }
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));


      problem.add_dis("fluid", fluiddis);
      if (xfluiddis != Teuchos::null)
      {
        problem.add_dis("xfluid", xfluiddis);  // xfem discretization on slot 1
      }
      problem.add_dis("ale", aledis);

      if (xfluiddis != Teuchos::null)
      {
        meshreader.add_element_reader(
            Core::IO::ElementReader(xfluiddis, reader, "--FLUID ELEMENTS"));
      }
      else
        meshreader.add_element_reader(
            Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case Core::ProblemType::tsi:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "structure", reader.get_comm(), problem.n_dim()));
          thermdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "thermo", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(
              new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
          thermdis = Teuchos::rcp(
              new Core::FE::Discretization("thermo", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      thermdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(thermdis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("thermo", thermdis);

      meshreader.add_advanced_reader(structdis, reader, "STRUCTURE",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.structural_dynamic_params(), "GEOMETRY"),
          nullptr);
      meshreader.add_advanced_reader(thermdis, reader, "THERMO",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.thermal_dynamic_params(), "GEOMETRY"),
          nullptr);

      break;
    }
    case Core::ProblemType::thermo:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          thermdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "thermo", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          thermdis = Teuchos::rcp(
              new Core::FE::Discretization("thermo", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      thermdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(thermdis, output_control, distype)));

      problem.add_dis("thermo", thermdis);

      meshreader.add_element_reader(Core::IO::ElementReader(thermdis, reader, "--THERMO ELEMENTS"));

      break;
    }

    case Core::ProblemType::structure:
    {
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "structure", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(
              new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));

      problem.add_dis("structure", structdis);

      meshreader.add_advanced_reader(structdis, reader, "STRUCTURE",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.structural_dynamic_params(), "GEOMETRY"),
          nullptr);

      break;
    }

    case Core::ProblemType::polymernetwork:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      pboxdis = Teuchos::rcp(
          new Core::FE::Discretization("boundingbox", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      pboxdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(pboxdis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("boundingbox", pboxdis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(pboxdis, reader, "--PERIODIC BOUNDINGBOX ELEMENTS"));

      break;
    }

    case Core::ProblemType::loma:
    {
      // create empty discretizations
      fluiddis = Teuchos::rcp(
          new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
      scatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.add_dis("fluid", fluiddis);
      problem.add_dis("scatra", scatradis);

      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }

    case Core::ProblemType::fluid_xfem_ls:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      if (problem.get_problem_type() == Core::ProblemType::fluid_xfem_ls)
        fluiddis =
            Teuchos::rcp(new XFEM::DiscretizationXFEM("fluid", reader.get_comm(), problem.n_dim()));
      else
        fluiddis = Teuchos::rcp(
            new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
      scatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("fluid", fluiddis);
      problem.add_dis("scatra", scatradis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_advanced_reader(fluiddis, reader, "FLUID",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.fluid_dynamic_params(), "GEOMETRY"),
          nullptr);
      // meshreader.AddElementReader(Teuchos::rcp(new Core::IO::ElementReader(fluiddis, reader,
      // "--FLUID ELEMENTS")));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      break;
    }

    case Core::ProblemType::elch:
    {
      // create empty discretizations
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          fluiddis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "fluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "scatra", reader.get_comm(), problem.n_dim()));
          aledis = Teuchos::rcp(
              new Core::FE::Nurbs::NurbsDiscretization("ale", reader.get_comm(), problem.n_dim()));
          scatra_micro_dis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "scatra_micro", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          fluiddis = Teuchos::rcp(
              new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(
              new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));
          aledis =
              Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));
          scatra_micro_dis = Teuchos::rcp(
              new Core::FE::Discretization("scatra_micro", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));
      scatra_micro_dis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(scatra_micro_dis, output_control, distype)));

      problem.add_dis("fluid", fluiddis);
      problem.add_dis("scatra", scatradis);
      problem.add_dis("ale", aledis);
      problem.add_dis("scatra_micro", scatra_micro_dis);

      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatra_micro_dis, reader, "--TRANSPORT2 ELEMENTS"));

      break;
    }
    case Core::ProblemType::art_net:  // _1D_ARTERY_
    {
      // create empty discretizations
      arterydis =
          Teuchos::rcp(new Core::FE::Discretization("artery", reader.get_comm(), problem.n_dim()));

      // create empty discretizations
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          FOUR_C_THROW("Nurbs discretization not possible for artery");
          break;
        }
        default:
        {
          scatradis = Teuchos::rcp(
              new Core::FE::Discretization("artery_scatra", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      problem.add_dis("artery", arterydis);
      problem.add_dis("artery_scatra", scatradis);

      // create discretization writer - in constructor set into and owned by corresponding discret
      arterydis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(arterydis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      meshreader.add_element_reader(
          Core::IO::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      break;
    }
    case Core::ProblemType::red_airways:  // _reduced D airways
    {
      // create empty discretizations
      airwaydis = Teuchos::rcp(
          new Core::FE::Discretization("red_airway", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      airwaydis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(airwaydis, output_control, distype)));

      problem.add_dis("red_airway", airwaydis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS"));

      break;
    }
    case Core::ProblemType::struct_ale:  // structure with ale
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      aledis =
          Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("ale", aledis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(Core::IO::ElementReader(aledis, reader, "--ALE ELEMENTS"));

      break;
    }
    case Core::ProblemType::poroelast:
    case Core::ProblemType::poromultiphase:
    {
      // create empty discretizations
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "structure", reader.get_comm(), problem.n_dim()));
          porofluiddis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "porofluid", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(
              new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
          porofluiddis = Teuchos::rcp(
              new Core::FE::Discretization("porofluid", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("porofluid", porofluiddis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));

      if (Core::UTILS::IntegralValue<bool>(
              problem.poro_multi_phase_dynamic_params(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(
            new Core::FE::Discretization("artery", reader.get_comm(), problem.n_dim()));
        arterydis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.add_dis("artery", arterydis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));
      }

      break;
    }
    case Core::ProblemType::poromultiphasescatra:
    {
      // create empty discretizations
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          structdis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "structure", reader.get_comm(), problem.n_dim()));
          porofluiddis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "porofluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "scatra", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          structdis = Teuchos::rcp(
              new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
          porofluiddis = Teuchos::rcp(
              new Core::FE::Discretization("porofluid", reader.get_comm(), problem.n_dim()));
          scatradis = Teuchos::rcp(
              new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("porofluid", porofluiddis);
      problem.add_dis("scatra", scatradis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      if (Core::UTILS::IntegralValue<bool>(
              problem.poro_multi_phase_scatra_dynamic_params(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(
            new Core::FE::Discretization("artery", reader.get_comm(), problem.n_dim()));
        arterydis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.add_dis("artery", arterydis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));

        artscatradis = Teuchos::rcp(
            new Core::FE::Discretization("artery_scatra", reader.get_comm(), problem.n_dim()));
        artscatradis->set_writer(Teuchos::rcp(
            new Core::IO::DiscretizationWriter(artscatradis, output_control, distype)));
        problem.add_dis("artery_scatra", artscatradis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(artscatradis, reader, "--TRANSPORT ELEMENTS"));
      }

      break;
    }
    case Core::ProblemType::porofluidmultiphase:
    {
      // create empty discretizations
      switch (distype)
      {
        case Core::FE::ShapeFunctionType::nurbs:
        {
          porofluiddis = Teuchos::rcp(new Core::FE::Nurbs::NurbsDiscretization(
              "porofluid", reader.get_comm(), problem.n_dim()));
          break;
        }
        default:
        {
          porofluiddis = Teuchos::rcp(
              new Core::FE::Discretization("porofluid", reader.get_comm(), problem.n_dim()));
          break;
        }
      }

      // create discretization writer - in constructor set into and owned by corresponding discret
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));

      problem.add_dis("porofluid", porofluiddis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));

      if (Core::UTILS::IntegralValue<bool>(
              problem.poro_fluid_multi_phase_dynamic_params(), "ARTERY_COUPLING"))
      {
        arterydis = Teuchos::rcp(
            new Core::FE::Discretization("artery", reader.get_comm(), problem.n_dim()));
        arterydis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.add_dis("artery", arterydis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));
      }
      break;
    }
    case Core::ProblemType::fpsi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      porofluiddis = Teuchos::rcp(
          new Core::FE::Discretization("porofluid", reader.get_comm(), problem.n_dim()));
      fluiddis = Teuchos::rcp(
          new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
      aledis =
          Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("porofluid", porofluiddis);
      problem.add_dis("fluid", fluiddis);
      problem.add_dis("ale", aledis);

      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      break;
    }
    case Core::ProblemType::fbi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      fluiddis = Teuchos::rcp(
          new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("fluid", fluiddis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_advanced_reader(fluiddis, reader, "FLUID",
          Core::UTILS::IntegralValue<Core::IO::GeometryType>(
              problem.fluid_dynamic_params(), "GEOMETRY"),
          nullptr);

      break;
    }
    case Core::ProblemType::immersed_fsi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      fluiddis = Teuchos::rcp(
          new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("fluid", fluiddis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));

      break;
    }
    case Core::ProblemType::fps3i:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      porofluiddis = Teuchos::rcp(
          new Core::FE::Discretization("porofluid", reader.get_comm(), problem.n_dim()));
      fluiddis = Teuchos::rcp(
          new Core::FE::DiscretizationFaces("fluid", reader.get_comm(), problem.n_dim()));
      aledis =
          Teuchos::rcp(new Core::FE::Discretization("ale", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      fluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(fluiddis, output_control, distype)));
      aledis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(aledis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("porofluid", porofluiddis);
      problem.add_dis("fluid", fluiddis);
      problem.add_dis("ale", aledis);


      meshreader.add_element_reader(Core::IO::ElementReader(fluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      // fluid scatra field
      fluidscatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra1", reader.get_comm(), problem.n_dim()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      fluidscatradis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(fluidscatradis, output_control, distype)));
      problem.add_dis("scatra1", fluidscatradis);

      // poro structure scatra field
      structscatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra2", reader.get_comm(), problem.n_dim()));
      // create discretization writer - in constructor set into and owned by corresponding discret
      structscatradis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(structscatradis, output_control, distype)));
      problem.add_dis("scatra2", structscatradis);

      break;
    }
    case Core::ProblemType::poroscatra:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      porofluiddis = Teuchos::rcp(
          new Core::FE::Discretization("porofluid", reader.get_comm(), problem.n_dim()));
      scatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      porofluiddis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(porofluiddis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("porofluid", porofluiddis);
      problem.add_dis("scatra", scatradis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(porofluiddis, reader, "--FLUID ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      break;
    }
    case Core::ProblemType::ehl:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      lubricationdis = Teuchos::rcp(
          new Core::FE::Discretization("lubrication", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      lubricationdis->set_writer(Teuchos::rcp(
          new Core::IO::DiscretizationWriter(lubricationdis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("lubrication", lubricationdis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(lubricationdis, reader, "--LUBRICATION ELEMENTS"));

      break;
    }
    case Core::ProblemType::ssi:
    case Core::ProblemType::ssti:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      scatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("scatra", scatradis);

      // consider case of additional scatra manifold
      if (Core::UTILS::IntegralValue<bool>(
              problem.ssi_control_params().sublist("MANIFOLD"), "ADD_MANIFOLD"))
      {
        auto scatra_manifold_dis = Teuchos::rcp(
            new Core::FE::Discretization("scatra_manifold", reader.get_comm(), problem.n_dim()));
        scatra_manifold_dis->set_writer(Teuchos::rcp(
            new Core::IO::DiscretizationWriter(scatra_manifold_dis, output_control, distype)));
        problem.add_dis("scatra_manifold", scatra_manifold_dis);
      }

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));

      if (problem.get_problem_type() == Core::ProblemType::ssti)
      {
        thermdis = Teuchos::rcp(
            new Core::FE::Discretization("thermo", reader.get_comm(), problem.n_dim()));
        thermdis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(thermdis, output_control, distype)));
        problem.add_dis("thermo", thermdis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(thermdis, reader, "--TRANSPORT ELEMENTS"));
      }

      break;
    }
    case Core::ProblemType::particle:
    case Core::ProblemType::pasi:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));

      problem.add_dis("structure", structdis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));

      break;
    }
    case Core::ProblemType::level_set:
    {
      // create empty discretizations
      scatradis =
          Teuchos::rcp(new Core::FE::Discretization("scatra", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      scatradis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(scatradis, output_control, distype)));

      problem.add_dis("scatra", scatradis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(scatradis, reader, "--TRANSPORT ELEMENTS"));
      break;
    }
    case Core::ProblemType::np_support:
    {
      // no discretizations and nodes needed for supporting procs
      break;
    }
    case Core::ProblemType::elemag:
    {
      // create empty discretizations
      elemagdis = Teuchos::rcp(
          new Core::FE::DiscretizationHDG("elemag", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      elemagdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(elemagdis, output_control, distype)));

      problem.add_dis("elemag", elemagdis);

      std::set<std::string> elemagelementtypes;
      elemagelementtypes.insert("ELECTROMAGNETIC");
      elemagelementtypes.insert("ELECTROMAGNETICDIFF");

      meshreader.add_element_reader(Core::IO::ElementReader(
          elemagdis, reader, "--ELECTROMAGNETIC ELEMENTS", elemagelementtypes));

      break;
    }
    case Core::ProblemType::redairways_tissue:
    {
      // create empty discretizations
      structdis = Teuchos::rcp(
          new Core::FE::Discretization("structure", reader.get_comm(), problem.n_dim()));
      airwaydis = Teuchos::rcp(
          new Core::FE::Discretization("red_airway", reader.get_comm(), problem.n_dim()));

      // create discretization writer - in constructor set into and owned by corresponding discret
      structdis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis, output_control, distype)));
      airwaydis->set_writer(
          Teuchos::rcp(new Core::IO::DiscretizationWriter(airwaydis, output_control, distype)));

      problem.add_dis("structure", structdis);
      problem.add_dis("red_airway", airwaydis);

      meshreader.add_element_reader(
          Core::IO::ElementReader(structdis, reader, "--STRUCTURE ELEMENTS"));
      meshreader.add_element_reader(
          Core::IO::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS"));
    }
    break;
    default:
      FOUR_C_THROW("Unknown problem type: %d", problem.get_problem_type());
      break;
  }

  // add artery or airways discretizations only for the following problem types
  switch (problem.get_problem_type())
  {
    case Core::ProblemType::fsi_redmodels:
    case Core::ProblemType::fsi_lung:
    case Core::ProblemType::fluid_ale:
    case Core::ProblemType::fluid_redmodels:
    {
      if (distype == Core::FE::ShapeFunctionType::polynomial)
      {
        // create empty discretizations
        arterydis = Teuchos::rcp(
            new Core::FE::Discretization("artery", reader.get_comm(), problem.n_dim()));
        // create discretization writer - in constructor set into and owned by corresponding discret
        arterydis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(arterydis, output_control, distype)));
        problem.add_dis("artery", arterydis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(arterydis, reader, "--ARTERY ELEMENTS"));

        airwaydis = Teuchos::rcp(
            new Core::FE::Discretization("red_airway", reader.get_comm(), problem.n_dim()));
        // create discretization writer - in constructor set into and owned by corresponding discret
        airwaydis->set_writer(
            Teuchos::rcp(new Core::IO::DiscretizationWriter(airwaydis, output_control, distype)));
        problem.add_dis("red_airway", airwaydis);
        meshreader.add_element_reader(
            Core::IO::ElementReader(airwaydis, reader, "--REDUCED D AIRWAYS ELEMENTS"));
      }
    }
    break;
    default:
      break;
  }

  if (read_mesh)  // now read and allocate!
  {
    // we read nodes and elements for the desired fields as specified above
    meshreader.read_and_partition();

    // care for special applications
    switch (problem.get_problem_type())
    {
      case Core::ProblemType::elch:
      case Core::ProblemType::fsi:
      case Core::ProblemType::fsi_redmodels:
      case Core::ProblemType::fsi_lung:
      case Core::ProblemType::scatra:
      case Core::ProblemType::structure:
      {
        // read microscale fields from second, third, ... input file if necessary
        // (in case of multi-scale material models)
        ReadMicroFields(problem, reader);
        break;
      }
      case Core::ProblemType::np_support:
      {
        // read microscale fields from second, third, ... inputfile for supporting processors
        ReadMicrofieldsNPsupport(problem);
        break;
      }
      default:
        break;
    }
  }  // if(read_mesh)
}

void Global::ReadMicroFields(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  // check whether micro material is specified
  const int id_struct = Global::Problem::instance()->materials()->first_id_by_type(
      Core::Materials::m_struct_multiscale);
  const int id_scatra = Global::Problem::instance()->materials()->first_id_by_type(
      Core::Materials::m_scatra_multiscale);
  const int id_elch = Global::Problem::instance()->materials()->first_id_by_type(
      Core::Materials::m_newman_multiscale);

  // return if no multiscale material is used
  if (id_struct == -1 and id_scatra == -1 and id_elch == -1) return;

  // safety check
  if ((id_struct != -1 and id_scatra != -1) or (id_struct != -1 and id_elch != -1) or
      (id_scatra != -1 and id_elch != -1))
    FOUR_C_THROW("Cannot have more than one multi-scale material!");

  // store name of macro-scale discretization in string
  std::string macro_dis_name("");
  if (id_struct != -1)
    macro_dis_name = "structure";
  else
    macro_dis_name = "scatra";

  // fetch communicators
  Teuchos::RCP<Epetra_Comm> lcomm = problem.get_communicators()->local_comm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem.get_communicators()->global_comm();

  Global::Problem* macro_problem = Global::Problem::instance();
  Teuchos::RCP<Core::FE::Discretization> macro_dis = macro_problem->get_dis(macro_dis_name);

  // repartition macro problem for a good distribution of elements with micro material
  if (macro_dis_name == "structure")
  {
    // do weighted repartitioning to obtain new row/column maps
    const Teuchos::ParameterList rebalanceParams;
    Teuchos::RCP<const Epetra_CrsGraph> nodeGraph = macro_dis->build_node_graph();
    const auto& [nodeWeights, edgeWeights] = Core::Rebalance::BuildWeights(*macro_dis);
    const auto& [rownodes, colnodes] =
        Core::Rebalance::RebalanceNodeMaps(nodeGraph, rebalanceParams, nodeWeights, edgeWeights);

    // rebuild the discretization with new maps
    macro_dis->redistribute(*rownodes, *colnodes, true, true, true);
  }

  // make sure that we read the micro discretizations only on the processors on
  // which elements with the corresponding micro material are evaluated

  std::set<int> my_multimat_IDs;

  // take care also of ghosted elements! -> ElementColMap!
  for (int i = 0; i < macro_dis->element_col_map()->NumMyElements(); ++i)
  {
    Core::Elements::Element* actele = macro_dis->l_col_element(i);
    Teuchos::RCP<Core::Mat::Material> actmat = actele->material();

    if (id_elch != -1 and actmat->material_type() == Core::Materials::m_elchmat)
    {
      // extract wrapped material
      auto elchmat = Teuchos::rcp_dynamic_cast<const Mat::ElchMat>(actmat);
      auto elchphase = Teuchos::rcp_dynamic_cast<const Mat::ElchPhase>(
          elchmat->phase_by_id(elchmat->phase_id(0)));
      actmat = elchphase->mat_by_id(elchphase->mat_id(0));
    }

    if ((actmat->material_type() == Core::Materials::m_struct_multiscale and
            macro_dis_name == "structure") or
        (actmat->material_type() == Core::Materials::m_scatra_multiscale and
            macro_dis_name == "scatra") or
        (actmat->material_type() == Core::Materials::m_newman_multiscale and
            macro_dis_name == "scatra"))
    {
      Core::Mat::PAR::Parameter* actparams = actmat->parameter();
      my_multimat_IDs.insert(actparams->id());
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
    problem.get_communicators()->set_sub_comm(subgroupcomm);

    // find out how many micro problems have to be solved on this macro proc
    int microcount = 0;
    for (const auto& material_map : problem.materials()->map())
    {
      int matid = material_map.first;
      if (my_multimat_IDs.find(matid) != my_multimat_IDs.end()) microcount++;
    }
    // and broadcast it to the corresponding group of procs
    subgroupcomm->Broadcast(&microcount, 1, 0);

    for (const auto& material_map : problem.materials()->map())
    {
      int matid = material_map.first;

      if (my_multimat_IDs.find(matid) != my_multimat_IDs.end())
      {
        Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(matid);

        // initialize variables storing micro-scale information
        int microdisnum(-1);
        std::string micro_dis_name = "";
        std::string micro_inputfile_name("");
        Global::Problem* micro_problem(nullptr);

        // structure case
        if (macro_dis_name == "structure")
        {
          // access multi-scale structure material
          auto* micromat = static_cast<Mat::MicroMaterial*>(mat.get());

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->micro_dis_num();
          subgroupcomm->Broadcast(&microdisnum, 1, 0);

          // set name of micro-scale discretization
          micro_dis_name = "structure";

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->micro_input_file_name();

          // instantiate micro-scale problem
          micro_problem = Global::Problem::instance(microdisnum);
        }

        // scalar transport case
        else
        {
          // access multi-scale scalar transport material
          Mat::ScatraMicroMacroCoupling* micromat = nullptr;
          if (id_scatra != -1)
            micromat = dynamic_cast<Mat::ScatraMultiScale*>(mat.get());
          else if (id_elch != -1)
            micromat = dynamic_cast<Mat::NewmanMultiScale*>(mat.get());
          else
            FOUR_C_THROW("How the heck did you get here?!");

          // extract and broadcast number of micro-scale discretization
          microdisnum = micromat->micro_dis_num();
          subgroupcomm->Broadcast(&microdisnum, 1, 0);

          // set unique name of micro-scale discretization
          std::stringstream name;
          name << "scatra_multiscale_" << microdisnum;
          micro_dis_name = name.str();

          // extract name of micro-scale input file
          micro_inputfile_name = micromat->micro_input_file_name();

          // instantiate micro-scale problem
          micro_problem = Global::Problem::instance(microdisnum);
        }

        if (micro_inputfile_name[0] != '/')
        {
          std::string filename = reader.my_inputfile_name();
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
        Core::IO::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

        Teuchos::RCP<Core::FE::Discretization> dis_micro = Teuchos::rcp(
            new Core::FE::Discretization(micro_dis_name, micro_reader.get_comm(), problem.n_dim()));

        // replace standard dofset inside micro discretization by independent dofset
        // to avoid inconsistent dof numbering in non-nested parallel settings with more than one
        // micro discretization
        if (problem.get_communicators()->np_type() ==
            Core::Communication::NestedParallelismType::no_nested_parallelism)
          dis_micro->replace_dof_set(Teuchos::rcp(new Core::DOFSets::IndependentDofSet()));

        // create discretization writer - in constructor set into and owned by corresponding
        // discret
        dis_micro->set_writer(Teuchos::rcp(new Core::IO::DiscretizationWriter(dis_micro,
            micro_problem->output_control_file(), micro_problem->spatial_approximation_type())));

        micro_problem->add_dis(micro_dis_name, dis_micro);

        ReadParameter(*micro_problem, micro_reader);

        // read materials of microscale
        // CAUTION: materials for microscale cannot be read until
        // micro_reader is activated, since else materials will again be
        // read from macroscale inputfile. Besides, materials MUST be read
        // before elements are read since elements establish a connection
        // to the corresponding material! Thus do not change position of
        // function calls!
        problem.materials()->set_read_from_problem(microdisnum);

        ReadMaterials(*micro_problem, micro_reader);

        Core::IO::MeshReader micromeshreader(micro_reader, "--NODE COORDS",
            {.mesh_paritioning_parameters = Problem::instance()->mesh_partitioning_params(),
                .geometric_search_parameters = Problem::instance()->geometric_search_params(),
                .io_parameters = Problem::instance()->io_params()});

        if (micro_dis_name == "structure")
        {
          micromeshreader.add_element_reader(
              Core::IO::ElementReader(dis_micro, micro_reader, "--STRUCTURE ELEMENTS"));
        }
        else
          micromeshreader.add_element_reader(
              Core::IO::ElementReader(dis_micro, micro_reader, "--TRANSPORT ELEMENTS"));

        micromeshreader.read_and_partition();


        {
          Core::UTILS::FunctionManager function_manager;
          GlobalLegacyModuleCallbacks().AttachFunctionDefinitions(function_manager);
          function_manager.read_input(micro_reader);
          micro_problem->set_function_manager(std::move(function_manager));
        }

        ReadResult(*micro_problem, micro_reader);
        ReadConditions(*micro_problem, micro_reader);

        // At this point, everything for the microscale is read,
        // subsequent reading is only for macroscale
        dis_micro->fill_complete();

        // broadcast restart information
        int restart_step = problem.restart();
        subgroupcomm->Broadcast(&restart_step, 1, 0);
        problem.set_restart_step(restart_step);

        // set the problem number from which to call materials again to zero
        // (i.e. macro problem), cf. Mat::Factory!
        problem.materials()->reset_read_from_problem();
      }
    }
    problem.materials()->reset_read_from_problem();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadMicrofieldsNPsupport(Global::Problem& problem)
{
  Teuchos::RCP<Epetra_Comm> lcomm = problem.get_communicators()->local_comm();
  Teuchos::RCP<Epetra_Comm> gcomm = problem.get_communicators()->global_comm();

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
  problem.get_communicators()->set_sub_comm(subgroupcomm);

  // number of micro problems for this sub group
  int microcount = 0;
  subgroupcomm->Broadcast(&microcount, 1, 0);

  for (int n = 0; n < microcount; n++)
  {
    // broadcast microdis number
    int microdisnum = -1;
    subgroupcomm->Broadcast(&microdisnum, 1, 0);

    Global::Problem* micro_problem = Global::Problem::instance(microdisnum);

    // broadcast micro input file name
    int length = -1;
    std::string micro_inputfile_name;
    subgroupcomm->Broadcast(&length, 1, 0);
    micro_inputfile_name.resize(length);
    subgroupcomm->Broadcast((const_cast<char*>(micro_inputfile_name.c_str())), length, 0);

    // start with actual reading
    Core::IO::DatFileReader micro_reader(micro_inputfile_name, subgroupcomm, 1);

    Teuchos::RCP<Core::FE::Discretization> structdis_micro = Teuchos::rcp(
        new Core::FE::Discretization("structure", micro_reader.get_comm(), problem.n_dim()));

    // create discretization writer - in constructor set into and owned by corresponding discret
    structdis_micro->set_writer(Teuchos::rcp(new Core::IO::DiscretizationWriter(structdis_micro,
        micro_problem->output_control_file(), micro_problem->spatial_approximation_type())));

    micro_problem->add_dis("structure", structdis_micro);

    ReadParameter(*micro_problem, micro_reader);

    // read materials of microscale
    // CAUTION: materials for microscale cannot be read until
    // micro_reader is activated, since else materials will again be
    // read from macroscale inputfile. Besides, materials MUST be read
    // before elements are read since elements establish a connection
    // to the corresponding material! Thus do not change position of
    // function calls!
    problem.materials()->set_read_from_problem(microdisnum);

    ReadMaterials(*micro_problem, micro_reader);

    Core::IO::MeshReader micromeshreader(micro_reader, "--NODE COORDS",
        {.mesh_paritioning_parameters = Problem::instance()->mesh_partitioning_params(),
            .geometric_search_parameters = Problem::instance()->geometric_search_params(),
            .io_parameters = Problem::instance()->io_params()});
    micromeshreader.add_element_reader(
        Core::IO::ElementReader(structdis_micro, micro_reader, "--STRUCTURE ELEMENTS"));
    micromeshreader.read_and_partition();

    // read conditions of microscale
    // -> note that no time curves and spatial functions can be read!

    ReadConditions(*micro_problem, micro_reader);

    // At this point, everything for the microscale is read,
    // subsequent reading is only for macroscale
    structdis_micro->fill_complete();

    // broadcast restart information
    int restart_step = problem.restart();
    subgroupcomm->Broadcast(&restart_step, 1, 0);
    problem.set_restart_step(restart_step);

    // set the problem number from which to call materials again to zero
    // (i.e. macro problem), cf. Mat::Factory!
    problem.materials()->reset_read_from_problem();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadParameter(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  Teuchos::RCP<Teuchos::ParameterList> list = Teuchos::rcp(new Teuchos::ParameterList("DAT FILE"));

  reader.read_section("--DISCRETISATION", *list);
  reader.read_section("--PROBLEM SIZE", *list);
  reader.read_section("--PROBLEM TYP", *list);
  reader.read_section("--BINNING STRATEGY", *list);
  reader.read_section("--BOUNDINGVOLUME STRATEGY", *list);
  reader.read_section("--IO", *list);
  reader.read_section("--IO/EVERY ITERATION", *list);
  reader.read_section("--IO/MONITOR STRUCTURE DBC", *list);
  reader.read_section("--IO/RUNTIME VTK OUTPUT", *list);
  reader.read_section("--IO/RUNTIME VTK OUTPUT/FLUID", *list);
  reader.read_section("--IO/RUNTIME VTK OUTPUT/STRUCTURE", *list);
  reader.read_section("--IO/RUNTIME VTK OUTPUT/BEAMS", *list);
  reader.read_section("--IO/RUNTIME VTP OUTPUT STRUCTURE", *list);
  reader.read_section("--STRUCTURAL DYNAMIC", *list);
  reader.read_section("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY", *list);
  reader.read_section("--STRUCTURAL DYNAMIC/GENALPHA", *list);
  reader.read_section("--STRUCTURAL DYNAMIC/ONESTEPTHETA", *list);
  reader.read_section("--STRUCTURAL DYNAMIC/TIMEADAPTIVITY/JOINT EXPLICIT", *list);
  reader.read_section("--MORTAR COUPLING", *list);
  reader.read_section("--MORTAR COUPLING/PARALLEL REDISTRIBUTION", *list);
  reader.read_section("--CONTACT DYNAMIC", *list);
  reader.read_section("--CARDIOVASCULAR 0D-STRUCTURE COUPLING", *list);
  reader.read_section(
      "--CARDIOVASCULAR 0D-STRUCTURE COUPLING/SYS-PUL CIRCULATION PARAMETERS", *list);
  reader.read_section("--CARDIOVASCULAR 0D-STRUCTURE COUPLING/RESPIRATORY PARAMETERS", *list);
  reader.read_section("--BROWNIAN DYNAMICS", *list);
  reader.read_section("--BEAM INTERACTION", *list);
  reader.read_section("--BEAM INTERACTION/SPHERE BEAM LINK", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO BEAM CONTACT", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO SPHERE CONTACT", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO SOLID SURFACE CONTACT", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO SOLID SURFACE MESHTYING", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO SOLID SURFACE", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO SOLID SURFACE/RUNTIME VTK OUTPUT", *list);
  reader.read_section("--BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING", *list);
  reader.read_section(
      "--BEAM INTERACTION/BEAM TO SOLID VOLUME MESHTYING/RUNTIME VTK OUTPUT", *list);
  reader.read_section("--BEAM INTERACTION/CROSSLINKING", *list);
  reader.read_section("--THERMAL DYNAMIC", *list);
  reader.read_section("--THERMAL DYNAMIC/GENALPHA", *list);
  reader.read_section("--THERMAL DYNAMIC/ONESTEPTHETA", *list);
  reader.read_section("--TSI DYNAMIC", *list);
  reader.read_section("--TSI DYNAMIC/MONOLITHIC", *list);
  reader.read_section("--TSI DYNAMIC/PARTITIONED", *list);
  reader.read_section("--TSI CONTACT", *list);
  reader.read_section("--POROELASTICITY DYNAMIC", *list);
  reader.read_section("--POROSCATRA CONTROL", *list);
  reader.read_section("--POROFLUIDMULTIPHASE DYNAMIC", *list);
  reader.read_section("--POROFLUIDMULTIPHASE DYNAMIC/ARTERY COUPLING", *list);
  reader.read_section("--POROMULTIPHASE DYNAMIC", *list);
  reader.read_section("--POROMULTIPHASE DYNAMIC/PARTITIONED", *list);
  reader.read_section("--POROMULTIPHASE DYNAMIC/MONOLITHIC", *list);
  reader.read_section("--POROMULTIPHASESCATRA DYNAMIC", *list);
  reader.read_section("--POROMULTIPHASESCATRA DYNAMIC/PARTITIONED", *list);
  reader.read_section("--POROMULTIPHASESCATRA DYNAMIC/MONOLITHIC", *list);
  reader.read_section("--ELASTO HYDRO DYNAMIC", *list);
  reader.read_section("--ELASTO HYDRO DYNAMIC/PARTITIONED", *list);
  reader.read_section("--ELASTO HYDRO DYNAMIC/MONOLITHIC", *list);
  reader.read_section("--SSI CONTROL", *list);
  reader.read_section("--SSI CONTROL/ELCH", *list);
  reader.read_section("--SSI CONTROL/MANIFOLD", *list);
  reader.read_section("--SSI CONTROL/MONOLITHIC", *list);
  reader.read_section("--SSI CONTROL/PARTITIONED", *list);
  reader.read_section("--SSTI CONTROL", *list);
  reader.read_section("--SSTI CONTROL/MONOLITHIC", *list);
  reader.read_section("--SSTI CONTROL/THERMO", *list);
  reader.read_section("--FLUID DYNAMIC", *list);
  reader.read_section("--FLUID DYNAMIC/RESIDUAL-BASED STABILIZATION", *list);
  reader.read_section("--FLUID DYNAMIC/EDGE-BASED STABILIZATION", *list);
  reader.read_section("--FLUID DYNAMIC/POROUS-FLOW STABILIZATION", *list);
  reader.read_section("--FLUID DYNAMIC/TURBULENCE MODEL", *list);
  reader.read_section("--FLUID DYNAMIC/SUBGRID VISCOSITY", *list);
  reader.read_section("--FLUID DYNAMIC/WALL MODEL", *list);
  reader.read_section("--FLUID DYNAMIC/TIMEADAPTIVITY", *list);
  reader.read_section("--FLUID DYNAMIC/MULTIFRACTAL SUBGRID SCALES", *list);
  reader.read_section("--FLUID DYNAMIC/TURBULENT INFLOW", *list);
  reader.read_section("--FLUID DYNAMIC/NONLINEAR SOLVER TOLERANCES", *list);
  reader.read_section("--LUBRICATION DYNAMIC", *list);
  reader.read_section("--SCALAR TRANSPORT DYNAMIC", *list);
  reader.read_section("--SCALAR TRANSPORT DYNAMIC/NONLINEAR", *list);
  reader.read_section("--SCALAR TRANSPORT DYNAMIC/STABILIZATION", *list);
  reader.read_section("--SCALAR TRANSPORT DYNAMIC/S2I COUPLING", *list);
  reader.read_section("--SCALAR TRANSPORT DYNAMIC/ARTERY COUPLING", *list);
  reader.read_section("--SCALAR TRANSPORT DYNAMIC/EXTERNAL FORCE", *list);
  reader.read_section("--STI DYNAMIC", *list);
  reader.read_section("--STI DYNAMIC/MONOLITHIC", *list);
  reader.read_section("--STI DYNAMIC/PARTITIONED", *list);
  reader.read_section("--FS3I DYNAMIC", *list);
  reader.read_section("--FS3I DYNAMIC/PARTITIONED", *list);
  reader.read_section("--FS3I DYNAMIC/STRUCTURE SCALAR STABILIZATION", *list);
  reader.read_section("--FS3I DYNAMIC/AC", *list);
  reader.read_section("--ALE DYNAMIC", *list);
  reader.read_section("--FSI DYNAMIC", *list);
  reader.read_section("--FSI DYNAMIC/CONSTRAINT", *list);
  reader.read_section("--FSI DYNAMIC/MONOLITHIC SOLVER", *list);
  reader.read_section("--FSI DYNAMIC/PARTITIONED SOLVER", *list);
  reader.read_section("--FSI DYNAMIC/TIMEADAPTIVITY", *list);
  reader.read_section("--FLUID BEAM INTERACTION", *list);
  reader.read_section("--FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING", *list);
  reader.read_section("--FLUID BEAM INTERACTION/BEAM TO FLUID MESHTYING/RUNTIME VTK OUTPUT", *list);
  reader.read_section("--IMMERSED METHOD", *list);
  reader.read_section("--IMMERSED METHOD/PARTITIONED SOLVER", *list);
  reader.read_section("--FPSI DYNAMIC", *list);
  reader.read_section("--ARTERIAL DYNAMIC", *list);
  reader.read_section("--REDUCED DIMENSIONAL AIRWAYS DYNAMIC", *list);
  reader.read_section("--COUPLED REDUCED-D AIRWAYS AND TISSUE DYNAMIC", *list);
  reader.read_section("--SEARCH TREE", *list);
  reader.read_section("--XFEM GENERAL", *list);
  reader.read_section("--CUT GENERAL", *list);
  reader.read_section("--XFLUID DYNAMIC", *list);
  reader.read_section("--XFLUID DYNAMIC/GENERAL", *list);
  reader.read_section("--XFLUID DYNAMIC/STABILIZATION", *list);
  reader.read_section("--XFLUID DYNAMIC/XFPSI MONOLITHIC", *list);
  reader.read_section("--LOMA CONTROL", *list);
  reader.read_section("--ELCH CONTROL", *list);
  reader.read_section("--ELCH CONTROL/DIFFCOND", *list);
  reader.read_section("--ELCH CONTROL/SCL", *list);
  reader.read_section("--BIOFILM CONTROL", *list);
  reader.read_section("--PARTICLE DYNAMIC", *list);
  reader.read_section("--PARTICLE DYNAMIC/INITIAL AND BOUNDARY CONDITIONS", *list);
  reader.read_section("--PARTICLE DYNAMIC/SPH", *list);
  reader.read_section("--PARTICLE DYNAMIC/DEM", *list);
  reader.read_section("--PASI DYNAMIC", *list);
  reader.read_section("--LEVEL-SET CONTROL", *list);
  reader.read_section("--LEVEL-SET CONTROL/REINITIALIZATION", *list);
  reader.read_section("--WEAR", *list);
  reader.read_section("--BEAM CONTACT", *list);
  reader.read_section("--BEAM CONTACT/RUNTIME VTK OUTPUT", *list);
  reader.read_section("--BEAM POTENTIAL", *list);
  reader.read_section("--BEAM POTENTIAL/RUNTIME VTK OUTPUT", *list);
  reader.read_section("--SEMI-SMOOTH PLASTICITY", *list);
  reader.read_section("--ELECTROMAGNETIC DYNAMIC", *list);
  reader.read_section("--VOLMORTAR COUPLING", *list);
  reader.read_section("--CARDIAC MONODOMAIN CONTROL", *list);
  reader.read_section("--MOR", *list);
  reader.read_section("--MESH PARTITIONING", *list);
  reader.read_section("--MULTI POINT CONSTRAINTS", *list);
  reader.read_section("--STRUCT NOX", *list);
  reader.read_section("--STRUCT NOX/Direction", *list);
  reader.read_section("--STRUCT NOX/Direction/Newton", *list);
  reader.read_section("--STRUCT NOX/Direction/Newton/Modified", *list);
  reader.read_section("--STRUCT NOX/Direction/Newton/Linear Solver", *list);
  reader.read_section("--STRUCT NOX/Direction/Steepest Descent", *list);
  reader.read_section("--STRUCT NOX/Line Search", *list);
  reader.read_section("--STRUCT NOX/Line Search/Full Step", *list);
  reader.read_section("--STRUCT NOX/Line Search/Backtrack", *list);
  reader.read_section("--STRUCT NOX/Line Search/Polynomial", *list);
  reader.read_section("--STRUCT NOX/Line Search/More'-Thuente", *list);
  reader.read_section("--STRUCT NOX/Pseudo Transient", *list);
  reader.read_section("--STRUCT NOX/Trust Region", *list);
  reader.read_section("--STRUCT NOX/Printing", *list);
  reader.read_section("--STRUCT NOX/Status Test", *list);
  reader.read_section("--STRUCT NOX/Solver Options", *list);

  // read in solver sections
  // Note: the maximum number of solver blocks in dat files is hardwired here.
  // If you change this do not forget to edit the corresponding parts in
  // validparameters.cpp, too!
  for (int i = 1; i < 10; i++)
  {
    std::stringstream ss;
    ss << "--SOLVER " << i;
    reader.read_section(ss.str(), *list);

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
          std::string filename = reader.my_inputfile_name();
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

  reader.read_section("--UMFPACK SOLVER", *list);


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
      std::string filename = reader.my_inputfile_name();
      std::string::size_type pos = filename.rfind('/');
      if (pos != std::string::npos)
      {
        std::string tmp = filename.substr(0, pos + 1);
        statustest_xmlfile->insert(statustest_xmlfile->begin(), tmp.begin(), tmp.end());
      }
    }
  }  // STRUCT NOX/Status Test

  // check for invalid parameters
  problem.set_parameter_list(list);

  //---------------------------------------------------------------------
  // Now we have successfully read the whole input file. It's time to access some data

  // 1) get the problem type
  const Teuchos::ParameterList& type = problem.problem_type_params();
  problem.set_problem_type(Core::UTILS::IntegralValue<Core::ProblemType>(type, "PROBLEMTYP"));

  // 2) get the spatial approximation type
  problem.set_spatial_approximation_type(
      Core::UTILS::IntegralValue<Core::FE::ShapeFunctionType>(type, "SHAPEFCT"));

  int restart_step = problem.restart();
  // 3) do the restart business with the four options we support (partially)
  if (restart_step == 0)
  {
    // no restart flag on the command line, so check the restart flag from the input file
    restart_step = type.get<int>("RESTART");
    problem.set_restart_step(restart_step);
  }
  else  // SetRestartStep() has been called before!
  {
    // There is a non-zero restart flag on the command line, so we ignore the input file.
    // The RESTART flag in the input file should be zero or have the same value!
    const int restartflaginfile = type.get<int>("RESTART");
    if ((restartflaginfile > 0) and (restartflaginfile != restart_step))
      FOUR_C_THROW("Restart flags in input file and command line are non-zero and different!");
  }

  // Set restart time based on walltime
  const double restartinterval = problem.io_params().get<double>("RESTARTWALLTIMEINTERVAL");
  const int restartevry = problem.io_params().get<int>("RESTARTEVRY");
  problem.restart_manager()->setup_restart_manager(restartinterval, restartevry);

  // 4) set random seed
  // time is in seconds, therefore we add the global processor id to obtain a unique seed on each
  // proc
  {
    int rs = type.get<int>("RANDSEED");
    if (rs < 0)
      rs = static_cast<int>(time(nullptr)) +
           42 * Global::Problem::instance(0)->get_communicators()->global_comm()->MyPID();

    srand((unsigned int)rs);  // Set random seed for stdlibrary. This is deprecated, as it does not
    // produce random numbers on some platforms!
    problem.random()->set_rand_seed((unsigned int)rs);  // Use this instead.
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadMaterials(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  // create list of known materials
  Teuchos::RCP<std::vector<Teuchos::RCP<Mat::MaterialDefinition>>> vm = Input::ValidMaterials();
  std::vector<Teuchos::RCP<Mat::MaterialDefinition>>& matlist = *vm;

  // test for each material definition (input file --MATERIALS section)
  // and store in #matmap_
  auto mmap = problem.materials();
  for (auto& mat : matlist)
  {
    std::vector<std::pair<int, Core::IO::InputParameterContainer>> read_materials =
        mat->read(reader);

    for (const auto& [id, data] : read_materials)
    {
      if (mmap->id_exists(id)) FOUR_C_THROW("More than one material with 'MAT %d'", id);

      // Here we call out to the factory code to create materials from generic input data.
      // We do this via a lambda function which is wrapped inside the LazyPtr. This allows to defer
      // the actual creation of the material until needed. Any other material parameters that are
      // needed during creation, are constructed automatically, when querying them from the list
      // of parameters.
      // Also, this line shows a design flaw where the parameter object needs to know the material
      // id that was chosen in the input file.
      mmap->insert(id, Core::UTILS::LazyPtr<Core::Mat::PAR::Parameter>(
                           [i = id, mat_type = mat->type(), d = data]()
                           { return Mat::make_parameter(i, mat_type, d); }));
    }
  }

  // We have read in all the materials and now we force construction of them all. The LazyPtr
  // ensures that the ordering does not matter. Note that we do not wait any longer for
  // construction, because materials might later be used in code sections that only run on proc 0.
  // Doing anything MPI-parallel inside the material constructors would then fail. Unfortunately,
  // such operations happen in the code base, thus we construct the materials here.
  for (const auto& [id, mat] : mmap->map())
  {
    // This is the point where the material is actually constructed via the side effect that we try
    // to access the material.
    (void)mat.get();
  }


  // check if every material was identified
  const std::string material_section = "--MATERIALS";
  std::vector<const char*> section = reader.section(material_section);
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
        FOUR_C_THROW("invalid material line in '%s'", name.c_str());

      // extract material ID
      int matid = -1;
      {
        char* ptr;
        matid = static_cast<int>(strtol(number.c_str(), &ptr, 10));
        if (ptr == number.c_str())
          FOUR_C_THROW("failed to read material object number '%s'", number.c_str());
      }

      FOUR_C_THROW_UNLESS(problem.materials()->id_exists(matid),
          "Material 'MAT %d' with name '%s' could not be identified", matid, name.c_str());
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadContactConstitutiveLaws(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  // create list of known contact constitutive laws
  Teuchos::RCP<std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>> vm =
      Input::ValidContactConstitutiveLaws();
  std::vector<Teuchos::RCP<CONTACT::CONSTITUTIVELAW::LawDefinition>>& coconstlawlist = *vm;

  // test for each contact constitutive law definition (input file --CONTACT CONSTITUTIVE LAWS
  // section) and store it
  for (auto& m : coconstlawlist)
  {
    // read contact constitutive law from DAT file of type
    m->read(problem, reader, problem.contact_constitutive_laws());
  }

  // check if every contact constitutive law was identified
  const std::string contact_const_laws = "--CONTACT CONSTITUTIVE LAWS";
  std::vector<const char*> section = reader.section(contact_const_laws);
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
        FOUR_C_THROW("invalid contact constitutive law line in '%s'", name.c_str());

      // extract contact constitutive law ID
      int coconstlawid = -1;
      {
        char* ptr;
        coconstlawid = static_cast<int>(strtol(number.c_str(), &ptr, 10));
        if (ptr == number.c_str())
          FOUR_C_THROW(
              "failed to read contact constitutive law object number '%s'", number.c_str());
      }

      // processed?
      if (problem.contact_constitutive_laws()->find(coconstlawid) == -1)
        FOUR_C_THROW("Contact constitutive law 'LAW %d' with name '%s' could not be identified",
            coconstlawid, name.c_str());
    }
  }

  // make fast access parameters
  problem.contact_constitutive_laws()->make_parameters();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadCloningMaterialMap(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  const std::vector<Input::LineDefinition> lines = Core::FE::valid_cloning_material_map_lines();

  // perform the actual reading and extract the input parameters
  std::vector<Input::LineDefinition> input_line_vec =
      Core::IO::DatFileUtils::read_all_lines_in_section(reader, "CLONING MATERIAL MAP", lines);
  for (const auto& input_line : input_line_vec)
  {
    // extract what was read from the input file
    std::string src_field;
    input_line.extract_string("SRC_FIELD", src_field);
    int src_matid(-1);
    input_line.extract_int("SRC_MAT", src_matid);
    std::string tar_field;
    input_line.extract_string("TAR_FIELD", tar_field);
    int tar_matid(-1);
    input_line.extract_int("TAR_MAT", tar_matid);

    // create the key pair
    std::pair<std::string, std::string> fields(src_field, tar_field);

    // enter the material pairing into the map
    std::pair<int, int> matmap(src_matid, tar_matid);
    problem.cloning_material_map()[fields].insert(matmap);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadResult(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  const auto lines = GlobalLegacyModuleCallbacks().valid_result_description_lines();

  // read design nodes <-> nodes, lines <-> nodes, surfaces <-> nodes, volumes <-> nodes
  const auto get_discretization_callback = [](const std::string& name) -> decltype(auto)
  { return *Global::Problem::instance()->get_dis(name); };
  std::vector<std::vector<std::vector<int>>> nodeset(4);
  reader.read_design("DNODE", nodeset[0], get_discretization_callback);
  reader.read_design("DLINE", nodeset[1], get_discretization_callback);
  reader.read_design("DSURF", nodeset[2], get_discretization_callback);
  reader.read_design("DVOL", nodeset[3], get_discretization_callback);
  problem.get_result_test_manager().set_node_set(nodeset);

  problem.get_result_test_manager().set_parsed_lines(
      Core::IO::DatFileUtils::read_all_lines_in_section(reader, "RESULT DESCRIPTION", lines));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadConditions(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  Teuchos::Time time("", true);
  if (reader.get_comm()->MyPID() == 0)
  {
    Core::IO::cout << "Read/generate conditions                          in....";
    Core::IO::cout.flush();
  }

  //--------------------------------------------- read generic node sets
  const auto get_discretization_callback = [](const std::string& name) -> decltype(auto)
  { return *Global::Problem::instance()->get_dis(name); };

  // read design nodes <-> nodes
  std::vector<std::vector<int>> dnode_fenode;
  reader.read_design("DNODE", dnode_fenode, get_discretization_callback);

  // read design lines <-> nodes
  std::vector<std::vector<int>> dline_fenode;
  reader.read_design("DLINE", dline_fenode, get_discretization_callback);

  // read design surfaces <-> nodes
  std::vector<std::vector<int>> dsurf_fenode;
  reader.read_design("DSURF", dsurf_fenode, get_discretization_callback);

  // read design volumes <-> nodes
  std::vector<std::vector<int>> dvol_fenode;
  reader.read_design("DVOL", dvol_fenode, get_discretization_callback);

  // check for meshfree discretisation to add node set topologies
  std::vector<std::vector<std::vector<int>>*> nodeset(4);
  nodeset[0] = &dnode_fenode;
  nodeset[1] = &dline_fenode;
  nodeset[2] = &dsurf_fenode;
  nodeset[3] = &dvol_fenode;

  // create list of known conditions
  Teuchos::RCP<std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>> vc =
      Input::ValidConditions();
  std::vector<Teuchos::RCP<Core::Conditions::ConditionDefinition>>& condlist = *vc;

  // test for each condition definition (input file condition section)
  // - read all conditions that match the definition
  // - add the nodal clouds to the conditions
  // - add the conditions to the appropriate discretizations
  //
  // Note that this will reset (un-fill_complete) the discretizations.
  for (auto& condition : condlist)
  {
    std::multimap<int, Teuchos::RCP<Core::Conditions::Condition>> cond;

    // read conditions from dat file
    condition->read(reader, cond);

    // add nodes to conditions
    std::multimap<int, Teuchos::RCP<Core::Conditions::Condition>>::const_iterator curr;
    for (curr = cond.begin(); curr != cond.end(); ++curr)
    {
      switch (curr->second->g_type())
      {
        case Core::Conditions::geometry_type_point:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dnode_fenode.size())
          {
            FOUR_C_THROW(
                "DPoint %d not in range [0:%d[\n"
                "DPoint condition on non existent DPoint?",
                curr->first, dnode_fenode.size());
          }
          curr->second->set_nodes(dnode_fenode[curr->first]);
          break;
        case Core::Conditions::geometry_type_line:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dline_fenode.size())
          {
            FOUR_C_THROW(
                "DLine %d not in range [0:%d[\n"
                "DLine condition on non existent DLine?",
                curr->first, dline_fenode.size());
          }
          curr->second->set_nodes(dline_fenode[curr->first]);
          break;
        case Core::Conditions::geometry_type_surface:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dsurf_fenode.size())
          {
            FOUR_C_THROW(
                "DSurface %d not in range [0:%d[\n"
                "DSurface condition on non existent DSurface?",
                curr->first, dsurf_fenode.size());
          }
          curr->second->set_nodes(dsurf_fenode[curr->first]);
          break;
        case Core::Conditions::geometry_type_volume:
          if (curr->first < 0 or static_cast<unsigned>(curr->first) >= dvol_fenode.size())
          {
            FOUR_C_THROW(
                "DVolume %d not in range [0:%d[\n"
                "DVolume condition on non existent DVolume?",
                curr->first, dvol_fenode.size());
          }
          curr->second->set_nodes(dvol_fenode[curr->first]);
          break;
        default:
          FOUR_C_THROW("geometry type unspecified");
          break;
      }

      // Iterate through all discretizations and sort the appropriate condition
      // into the correct discretization it applies to

      for (const auto& [name, dis] : problem.discretization_range())
      {
        const std::vector<int>* nodes = curr->second->get_nodes();
        if (nodes->size() == 0)
          FOUR_C_THROW("%s condition %d has no nodal cloud", condition->description().c_str(),
              curr->second->id());

        int foundit = 0;
        for (int node : *nodes)
        {
          foundit = dis->have_global_node(node);
          if (foundit) break;
        }
        int found = 0;
        dis->get_comm().SumAll(&foundit, &found, 1);
        if (found)
        {
          // Insert a copy since we might insert the same condition in many discretizations.
          dis->set_condition(condition->name(), curr->second->copy_without_geometry());
        }
      }
    }
  }

  if (reader.get_comm()->MyPID() == 0)
  {
    std::cout << time.totalElapsedTime(true) << " secs\n";
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadKnots(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  // get information on the spatial approximation --- we only read knots
  // in the nurbs case
  Core::FE::ShapeFunctionType distype = problem.spatial_approximation_type();

  // Iterate through all discretizations and sort the appropriate condition
  // into the correct discretization it applies to

  for (const auto& [name, dis] : problem.discretization_range())
  {
    if (distype == Core::FE::ShapeFunctionType::nurbs)
    {
      // cast discretisation to nurbs variant to be able
      // to add the knotvector
      auto* nurbsdis = dynamic_cast<Core::FE::Nurbs::NurbsDiscretization*>(&(*dis));

      if (nurbsdis == nullptr)
        FOUR_C_THROW("discretization %s is not a NurbsDiscretization! Panic.", dis->name().c_str());

      // define an empty knot vector object
      Teuchos::RCP<Core::FE::Nurbs::Knotvector> disknots = Teuchos::null;

      // read the knotvector data from the input
      reader.read_knots(dis->name(), disknots);

      if (disknots == Teuchos::null)
      {
        FOUR_C_THROW("Knotvector read failed in Nurbs discretisation\n");
      }

      // make sure atdis is fillcompleted, to be able to call
      // ElementRowMap() on it
      // do not initialize elements, since this would require knot
      // vector values
      if (!dis->filled())
      {
        dis->fill_complete(false, false, false);
      }

      // the smallest gid in the discretisation determines the access
      // pattern via the element offset
      int smallest_gid_in_dis = dis->element_row_map()->MinAllGID();

      // consistency checks
      disknots->finish_knots(smallest_gid_in_dis);

      // add knots to discretisation
      nurbsdis->set_knot_vector(disknots);
    }
  }  // loop fields
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Global::ReadParticles(Global::Problem& problem, Core::IO::DatFileReader& reader)
{
  // no need to read in particles in case of restart
  if (problem.restart()) return;

  // the basic particle reader
  Input::ParticleReader particlereader(reader, "--PARTICLES");

  // do the actual reading of particles
  particlereader.read(problem.particles());
}

FOUR_C_NAMESPACE_CLOSE
