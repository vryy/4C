/*----------------------------------------------------------------------*/
/*! \file

 \brief main routine of the main postprocessor filters

\level 1

 */

#include "baci_post_common.H"

#include <Epetra_MpiComm.h>

#include "baci_post_common_writer_base.H"
#include "baci_post_common_single_field_writers.H"

#include "baci_utils_exceptions.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_discret.H"

#include "baci_scatra_ele.H"


void runEnsightVtuFilter(PostProblem& problem)
{
  // each problem type is different and writes different results
  switch (problem.Problemtype())
  {
    case ProblemType::fsi:
    case ProblemType::fsi_redmodels:
    case ProblemType::fsi_lung:
    {
      std::string basename = problem.outname();
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      PostField* alefield = problem.get_discretization(2);
      AleFilter alewriter(alefield, basename);
      alewriter.WriteFiles();
      // 1d artery
      if (problem.num_discr() == 4)
      {
        PostField* field = problem.get_discretization(2);
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();
      }
      if (problem.num_discr() > 2 and problem.get_discretization(2)->name() == "xfluid")
      {
        std::string basename = problem.outname();
        PostField* fluidfield = problem.get_discretization(2);
        FluidFilter xfluidwriter(fluidfield, basename);
        xfluidwriter.WriteFiles();
      }
      break;
    }
    case ProblemType::gas_fsi:
    case ProblemType::ac_fsi:
    case ProblemType::thermo_fsi:
    {
      std::string basename = problem.outname();
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      int numdisc = problem.num_discr();

      for (int i = 0; i < numdisc - 3; ++i)
      {
        PostField* scatrafield = problem.get_discretization(3 + i);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case ProblemType::biofilm_fsi:
    {
      std::string basename = problem.outname();
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      int numdisc = problem.num_discr();

      for (int i = 0; i < numdisc - 4; ++i)
      {
        PostField* scatrafield = problem.get_discretization(3 + i);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case ProblemType::struct_ale:
    {
      PostField* structurefield = problem.get_discretization(0);
      StructureFilter structwriter(
          structurefield, problem.outname(), problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* alefield = problem.get_discretization(1);
      AleFilter alewriter(alefield, problem.outname());
      alewriter.WriteFiles();
      break;
    }
    case ProblemType::structure:
    {
      // Regular solid/structure output
      {
        PostField* structurefield = problem.get_discretization(0);
        StructureFilter structwriter(structurefield, problem.outname(), problem.stresstype(),
            problem.straintype(), problem.optquantitytype());
        structwriter.WriteFiles();
      }

      // Deal with contact / meshtying problems
      if (problem.DoMortarInterfaces())
      {
        /* Loop over all mortar interfaces and process each interface individually
         *
         * Start at i = 1 since discretization '0' is the structure discretization.
         *
         * Note: We assume that there is only one structure discretization and
         * that all other discretizations are mortar interface discretizations.
         */
        for (int i = 1; i < problem.num_discr(); ++i)
        {
          PostField* mortarfield = problem.get_discretization(i);
          MortarFilter mortarwriter(mortarfield, problem.outname());
          mortarwriter.WriteFiles();
        }
      }

      break;
    }
    case ProblemType::polymernetwork:
    {
      int numdiscr = problem.num_discr();
      for (int i = 0; i < numdiscr; ++i)
      {
        std::string disname = problem.get_discretization(i)->name();
        if (disname == "structure")
        {
          PostField* structure = problem.get_discretization(i);
          StructureFilter writer(structure, problem.outname());
          writer.WriteFiles();
        }
        else if (disname == "ia_structure")
        {
          PostField* ia_structure = problem.get_discretization(i);
          StructureFilter iawriter(ia_structure, problem.outname());
          iawriter.WriteFiles();
        }
        else if (disname == "boundingbox")
        {
          PostField* boxdiscret = problem.get_discretization(i);
          StructureFilter boxwriter(boxdiscret, problem.outname());
          boxwriter.WriteFiles();
        }
        else if (disname == "bins")
        {
          PostField* visualizebins = problem.get_discretization(i);
          StructureFilter binwriter(visualizebins, problem.outname());
          binwriter.WriteFiles();
        }
        else
        {
          dserror("unknown discretization for postprocessing of polymer network problem!");
        }
      }
      break;
    }
    case ProblemType::fluid:
    {
      if (problem.num_discr() == 2)
      {
        if (problem.get_discretization(1)->name() == "xfluid")
        {
          std::string basename = problem.outname();
          PostField* fluidfield = problem.get_discretization(1);
          XFluidFilter xfluidwriter(fluidfield, basename);
          xfluidwriter.WriteFiles();
        }
      }
      // don't add break here !!!!
    }
    case ProblemType::fluid_redmodels:
    {
      if (problem.num_discr() == 2)
      {
        // 1d artery
        std::string basename = problem.outname();
        PostField* field = problem.get_discretization(1);
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();
        if (problem.get_discretization(1)->name() == "xfluid")
        {
          std::string basename = problem.outname();
          PostField* fluidfield = problem.get_discretization(1);
          XFluidFilter xfluidwriter(fluidfield, basename);
          xfluidwriter.WriteFiles();
        }
      }
      // don't add break here !!!!
    }
    case ProblemType::fluid_ale:
    case ProblemType::freesurf:
    {
      PostField* field = problem.get_discretization(0);
      FluidFilter writer(field, problem.outname());
      writer.WriteFiles();
      if (problem.num_discr() > 1 and problem.get_discretization(1)->name() == "xfluid")
      {
        std::string basename = problem.outname();
        PostField* fluidfield = problem.get_discretization(1);
        FluidFilter xfluidwriter(fluidfield, basename);
        xfluidwriter.WriteFiles();
      }
      break;
    }
    case ProblemType::particle:
    case ProblemType::pasi:
    {
      int numdiscr = problem.num_discr();
      for (int i = 0; i < numdiscr; ++i)
      {
        std::string disname = problem.get_discretization(i)->name();
        if (disname == "bins")
        {
          PostField* visualizebins = problem.get_discretization(i);
          StructureFilter binwriter(visualizebins, problem.outname());
          binwriter.WriteFiles();
        }
        else if (disname == "structure")
        {
          PostField* structure = problem.get_discretization(i);
          StructureFilter writer(structure, problem.outname(), problem.stresstype(),
              problem.straintype(), problem.optquantitytype());
          writer.WriteFiles();
        }
        else
        {
          dserror("Particle problem has illegal discretization name!");
        }
      }
      break;
    }
    case ProblemType::level_set:
    {
      std::string basename = problem.outname();

      PostField* scatrafield = problem.get_discretization(0);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      break;
    }
    case ProblemType::redairways_tissue:
    {
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, problem.outname(), problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      StructureFilter fluidwriter(
          fluidfield, problem.outname(), problem.stresstype(), problem.straintype());
      fluidwriter.WriteFiles();
      break;
    }
    case ProblemType::ale:
    {
      PostField* field = problem.get_discretization(0);
      AleFilter writer(field, problem.outname());
      writer.WriteFiles();
      break;
    }
    case ProblemType::lubrication:
    {
      PostField* lubricationfield = problem.get_discretization(0);
      LubricationFilter lubricationwriter(lubricationfield, problem.outname());
      lubricationwriter.WriteFiles();
      break;
    }
    case ProblemType::porofluidmultiphase:
    {
      std::string basename = problem.outname();

      PostField* field = problem.get_discretization(0);
      PoroFluidMultiPhaseFilter writer(field, problem.outname());
      writer.WriteFiles();

      // write output for artery
      if (problem.num_discr() == 2)
      {
        PostField* field = problem.get_discretization(1);
        // AnyFilter writer(field, problem.outname());
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();
      }
      break;
    }
    case ProblemType::poromultiphase:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      PoroFluidMultiPhaseFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();
      if (problem.num_discr() == 3)
      {
        // artery
        PostField* field = problem.get_discretization(2);
        // AnyFilter writer(field, problem.outname());
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();
      }
      break;
    }
    case ProblemType::poromultiphasescatra:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      PoroFluidMultiPhaseFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      // no artery discretization
      if (problem.num_discr() == 3)
      {
        PostField* scatrafield = problem.get_discretization(2);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }
      else if (problem.num_discr() == 4)
      {
        // artery
        PostField* field = problem.get_discretization(2);
        // AnyFilter writer(field, problem.outname());
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();

        // scatra
        PostField* scatrafield = problem.get_discretization(3);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }
      else if (problem.num_discr() == 5)
      {
        // artery
        PostField* field = problem.get_discretization(2);
        // AnyFilter writer(field, problem.outname());
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();

        // artery scatra
        PostField* artscatrafield = problem.get_discretization(3);
        ScaTraFilter artscatrawriter(artscatrafield, basename);
        artscatrawriter.WriteFiles();

        // scatra
        PostField* scatrafield = problem.get_discretization(4);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }
      else
        dserror("wrong number of discretizations");
      break;
    }
    case ProblemType::cardiac_monodomain:
    case ProblemType::scatra:
    {
      std::string basename = problem.outname();
      // do we have a fluid discretization?
      int numfields = problem.num_discr();
      if (numfields == 2)
      {
        PostField* fluidfield = problem.get_discretization(0);
        FluidFilter fluidwriter(fluidfield, basename);
        fluidwriter.WriteFiles();

        PostField* scatrafield = problem.get_discretization(1);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }
      else if (numfields == 1)
      {
        PostField* scatrafield = problem.get_discretization(0);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }
      else
        dserror("number of fields does not match: got %d", numfields);

      break;
    }
    case ProblemType::sti:
    {
      // extract label for output files
      std::string basename = problem.outname();

      // safety check
      if (problem.num_discr() != 2)
        dserror("Must have exactly two discretizations for scatra-thermo interaction problems!");

      DRT::ELEMENTS::Transport* transport_element = dynamic_cast<DRT::ELEMENTS::Transport*>(
          problem.get_discretization(0)->discretization()->lRowElement(0));
      if (transport_element == NULL)
        dserror("Elements of unknown type on scalar transport discretization!");

      if (transport_element->ImplType() == INPAR::SCATRA::impltype_elch_electrode_thermo or
          transport_element->ImplType() == INPAR::SCATRA::impltype_elch_diffcond_thermo)
      {
        ElchFilter elchwriter(problem.get_discretization(0), basename);
        elchwriter.WriteFiles();
      }
      else
      {
        dserror("Scatra-thermo interaction not yet implemented for standard scalar transport!");
        ScaTraFilter scatrawriter(problem.get_discretization(0), basename);
        scatrawriter.WriteFiles();
      }

      ScaTraFilter thermowriter(problem.get_discretization(1), basename);
      thermowriter.WriteFiles();

      break;
    }
    case ProblemType::fsi_xfem:
    case ProblemType::fpsi_xfem:
    case ProblemType::fluid_xfem_ls:
    case ProblemType::two_phase_flow:
    {
      std::cout << "|=============================================================================|"
                << std::endl;
      std::cout << "|==  Output for General Problem" << std::endl;

      int numfield = problem.num_discr();

      std::cout << "|==  Number of discretizations: " << numfield << std::endl;

      std::string basename = problem.outname();

      std::cout << "\n|==  Start postprocessing for discretizations:" << std::endl;

      for (int idx_int = 0; idx_int < numfield; idx_int++)
      {
        std::cout
            << "\n|=============================================================================|"
            << std::endl;
        std::string disname = problem.get_discretization(idx_int)->name();
        PostField* field = problem.get_discretization(idx_int);
        if (disname == "structure")
        {
          std::cout << "|==  Structural Field ( " << disname << " )" << std::endl;
          StructureFilter structwriter(field, basename, problem.stresstype(), problem.straintype());
          structwriter.WriteFiles();
        }
        else if (disname == "fluid" or disname == "xfluid" or disname == "porofluid")
        {
          std::cout << "|==    Fluid Field ( " << disname << " )" << std::endl;
          FluidFilter fluidwriter(field, basename);
          fluidwriter.WriteFiles();
        }
        else if (disname == "scatra")
        {
          std::cout << "|==    Scatra Field ( " << disname << " )" << std::endl;
          ScaTraFilter scatrawriter(field, basename);
          scatrawriter.WriteFiles();
        }
        else if (disname == "ale")
        {
          std::cout << "|==    Ale Field ( " << disname << " )" << std::endl;
          //          AleFilter alewriter(field, basename);
          //          alewriter.WriteFiles();
        }
        else if (disname.compare(1, 12, "boundary_of_"))
        {
          std::cout << "|==    Interface Field ( " << disname << " )" << std::endl;
          InterfaceFilter ifacewriter(field, basename);
          ifacewriter.WriteFiles();
        }
        else
          dserror(
              "You try to postprocess a discretization with name %s, maybe you should add it here?",
              disname.c_str());
      }
      std::cout << "|=============================================================================|"
                << std::endl;
      break;
    }
    case ProblemType::fluid_xfem:
    {
      std::cout << "Output FLUID-XFEM Problem" << std::endl;

      int numfield = problem.num_discr();

      std::cout << "Number of discretizations: " << numfield << std::endl;
      for (int i = 0; i < numfield; i++)
      {
        std::cout << "dis-name i=" << i << ": " << problem.get_discretization(i)->name()
                  << std::endl;
      }

      if (numfield == 0) dserror("we expect at least a fluid field, numfield=%i", numfield);
      std::string basename = problem.outname();

      // XFluid in the standard case, embedded fluid for XFF
      std::cout << "  Fluid Field" << std::endl;
      PostField* fluidfield = problem.get_discretization(0);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      // start index for interface discretizations
      int idx_int = 1;
      if (problem.num_discr() > 1 and problem.get_discretization(1)->name() == "xfluid")
      {
        // XFluid for XFF
        std::cout << "  XFluid Field" << std::endl;
        PostField* xfluidfield = problem.get_discretization(1);
        FluidFilter xfluidwriter(xfluidfield, basename);
        xfluidwriter.WriteFiles();
        idx_int += 1;
      }

      // all other fields are interface fields
      for (int i = idx_int; i < numfield; i++)
      {
        std::cout << "  Interface Field ( " << problem.get_discretization(i)->name() << " )"
                  << std::endl;
        PostField* ifacefield = problem.get_discretization(i);
        InterfaceFilter ifacewriter(ifacefield, basename);
        ifacewriter.WriteFiles();
      }
      break;
    }
    case ProblemType::loma:
    {
      std::string basename = problem.outname();

      PostField* fluidfield = problem.get_discretization(0);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      PostField* scatrafield = problem.get_discretization(1);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();
      break;
    }
    case ProblemType::elch:
    {
      std::string basename = problem.outname();
      int numfield = problem.num_discr();

      if (numfield == 3)
      {
        // Fluid, ScaTra and ALE fields are present
        PostField* fluidfield = problem.get_discretization(0);
        FluidFilter fluidwriter(fluidfield, basename);
        fluidwriter.WriteFiles();

        PostField* scatrafield = problem.get_discretization(1);
        ElchFilter elchwriter(scatrafield, basename);
        elchwriter.WriteFiles();

        PostField* alefield = problem.get_discretization(2);
        AleFilter alewriter(alefield, basename);
        alewriter.WriteFiles();
      }
      else if (numfield == 2)
      {
        auto* const dis0 = problem.get_discretization(0);
        auto* const dis1 = problem.get_discretization(1);

        if (dis0->name() == "scatra" and dis1->name() == "scatra_micro")
        {
          PostField* scatrafield = dis0;
          ElchFilter elchwriter(scatrafield, basename);
          elchwriter.WriteFiles();

          PostField* scatrafield_micro = dis1;
          ElchFilter elchwriter_micro(scatrafield_micro, basename);
          elchwriter_micro.WriteFiles();
        }
        else
        {
          // Fluid and ScaTra fields are present
          PostField* fluidfield = dis0;
          FluidFilter fluidwriter(fluidfield, basename);
          fluidwriter.WriteFiles();

          PostField* scatrafield = dis1;
          ElchFilter elchwriter(scatrafield, basename);
          elchwriter.WriteFiles();
        }
        break;
      }
      else if (numfield == 1)
      {
        // only a ScaTra field is present
        PostField* scatrafield = problem.get_discretization(0);
        ElchFilter elchwriter(scatrafield, basename);
        elchwriter.WriteFiles();
      }
      else
        dserror("number of fields does not match: got %d", numfield);
      break;
    }
    case ProblemType::art_net:
    {
      std::string basename = problem.outname();
      PostField* field = problem.get_discretization(0);
      // AnyFilter writer(field, problem.outname());
      StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
      writer.WriteFiles();

      // write output for scatra
      if (problem.num_discr() == 2)
      {
        PostField* scatrafield = problem.get_discretization(1);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case ProblemType::thermo:
    {
      PostField* field = problem.get_discretization(0);
      ThermoFilter writer(field, problem.outname(), problem.heatfluxtype(), problem.tempgradtype());
      writer.WriteFiles();
      break;
    }
    case ProblemType::tsi:
    {
      std::cout << "Output TSI Problem" << std::endl;

      std::string basename = problem.outname();

      PostField* thermfield = problem.get_discretization(0);
      ThermoFilter thermwriter(
          thermfield, basename, problem.heatfluxtype(), problem.tempgradtype());
      thermwriter.WriteFiles();

      PostField* structfield = problem.get_discretization(1);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();
      break;
    }
    case ProblemType::red_airways:
    {
      std::string basename = problem.outname();
      PostField* field = problem.get_discretization(0);
      // AnyFilter writer(field, problem.outname());
      StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
      writer.WriteFiles();
      //      writer.WriteFiles();

      break;
    }
    case ProblemType::poroelast:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();
      break;
    }
    case ProblemType::poroscatra:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      PostField* scatrafield = problem.get_discretization(2);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      break;
    }
    case ProblemType::fpsi:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* porofluidfield = problem.get_discretization(1);
      FluidFilter porofluidwriter(porofluidfield, basename);
      porofluidwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(2);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      break;
    }
    case ProblemType::immersed_fsi:
    case ProblemType::fbi:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      break;
    }
    case ProblemType::fps3i:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* porofluidfield = problem.get_discretization(1);
      FluidFilter porofluidwriter(porofluidfield, basename);
      porofluidwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(2);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      /////////////

      int numdisc = problem.num_discr();

      for (int i = 0; i < numdisc - 4; ++i)
      {
        PostField* scatrafield = problem.get_discretization(4 + i);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case ProblemType::ehl:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* lubricationfield = problem.get_discretization(1);
      LubricationFilter lubricationwriter(lubricationfield, basename);
      lubricationwriter.WriteFiles();

      break;
    }
    case ProblemType::ssi:
    {
      std::string basename = problem.outname();

      const int numfields = problem.num_discr();

      // remark: scalar transport discretization number is one for old structural time integration!
      PostField* scatrafield = problem.get_discretization(0);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      if (numfields == 2)
      {
        // remark: structure discretization number is zero for old structural time integration!
        PostField* structfield = problem.get_discretization(1);
        StructureFilter structwriter(
            structfield, basename, problem.stresstype(), problem.straintype());
        structwriter.WriteFiles();
      }
      else if (numfields == 3)
      {
        // remark: structure discretization number is zero for old structural time integration!
        PostField* structfield = problem.get_discretization(2);
        StructureFilter structwriter(
            structfield, basename, problem.stresstype(), problem.straintype());
        structwriter.WriteFiles();

        PostField* scatra_manifoldfield = problem.get_discretization(1);
        ScaTraFilter scatra_manifoldfieldwriter(scatra_manifoldfield, basename);
        scatra_manifoldfieldwriter.WriteFiles();
      }
      else
        dserror("Unknwon number of solution fields");

      break;
    }

    case ProblemType::ssti:
    {
      std::string basename = problem.outname();

      // remark: scalar transport discretization number is one for old structural time integration!
      PostField* scatrafield = problem.get_discretization(0);

      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      // remark: structure discretization number is zero for old structural time integration!
      PostField* structfield = problem.get_discretization(2);

      StructureFilter structwriter(
          structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* thermofield = problem.get_discretization(1);
      ScaTraFilter thermowriter(thermofield, basename);
      thermowriter.WriteFiles();

      break;
    }
    case ProblemType::elemag:
    {
      PostField* field = problem.get_discretization(0);
      ElemagFilter writer(field, problem.outname());
      writer.WriteFiles();
      break;
    }
    case ProblemType::invana:
    {
      PostField* field = problem.get_discretization(0);
      InvanaFilter writer(field, problem.outname());
      writer.WriteFiles();
      break;
    }
    case ProblemType::none:
    {
      // Special problem type that contains one discretization and any number
      // of vectors. We just want to see whatever there is.
      PostField* field = problem.get_discretization(0);
      AnyFilter writer(field, problem.outname());
      writer.WriteFiles();
      break;
    }
    default:
      dserror("problem type %d not yet supported", problem.Problemtype());
      break;
  }
}



namespace
{
  std::string get_filter(int argc, char** argv)
  {
    Teuchos::CommandLineProcessor clp(false, false, false);
    std::string filter("ensight");
    clp.setOption("filter", &filter, "filter to run [ensight, vtu, vti]");
    // ignore warnings about unrecognized options.
    std::ostringstream warning;
    clp.parse(argc, argv, &warning);
    return filter;
  }
}  // namespace


/*!
 \brief post-processor main routine

 Select the appropriate filter and run!

 \author kronbichler
 \date 03/14
 */
int main(int argc, char** argv)
{
  try
  {
    std::string filter = get_filter(argc, argv);
    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString("Main BACI post-processor\n");

    PostProblem problem(My_CLP, argc, argv);

    if (filter == "ensight" || filter == "vtu" || filter == "vtu_node_based" || filter == "vti")
      runEnsightVtuFilter(problem);
    else
      dserror("Unknown filter %s given, supported filters: [ensight|vtu|vti]", filter.c_str());

  }  // try
  catch (std::runtime_error& err)
  {
    char line[] = "=========================================================================\n";
    std::cout << "\n\n" << line << err.what() << "\n" << line << "\n" << std::endl;

    // proper cleanup
    DRT::Problem::Done();
#ifdef DSERROR_DUMP
    abort();
#endif

    MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
  }  // catch

  // proper cleanup
  DRT::Problem::Done();

  return 0;
}
