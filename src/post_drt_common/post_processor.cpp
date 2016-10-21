/*!
 \file post_processor.cpp

 \brief main routine of the main postprocessor filters

\level 1

 \maintainer Martin Kronbichler
             kronbichler@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15235
 */

#include "post_drt_common.H"

#include "post_writer_base.H"
#include "post_single_field_writers.H"
#include "../post_gid/post_drt_gid.H"

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_scatra_ele/scatra_ele.H"


void runEnsightVtuFilter(PostProblem    &problem)
{
#if 0
    for (int i = 0; i<problem.num_discr(); ++i)
    {
        PostField* field = problem.get_discretization(i);
        StructureFilter writer(field, problem.outname());
        writer.WriteFiles();
    }
#endif

    // each problem type is different and writes different results
    switch (problem.Problemtype())
    {
    case prb_fsi:
    case prb_fsi_redmodels:
    case prb_fsi_lung:
    {
        std::string basename = problem.outname();
        PostField* structfield = problem.get_discretization(0);
        StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
        structwriter.WriteFiles();

        PostField* fluidfield = problem.get_discretization(1);
        FluidFilter fluidwriter(fluidfield, basename);
        fluidwriter.WriteFiles();

        //PostField* alefield = problem.get_discretization(2);
        //AleFilter alewriter(alefield, basename);
        //alewriter.WriteFiles();
#ifdef D_ARTNET
        // 1d artery
        if (problem.num_discr()== 4)
        {
          PostField* field = problem.get_discretization(2);
          StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
          writer.WriteFiles();
        }
#endif //D_ARTNET
        if (problem.num_discr() > 2 and problem.get_discretization(2)->name()=="xfluid")
        {
          std::string basename = problem.outname();
          PostField* fluidfield = problem.get_discretization(2);
          FluidFilter xfluidwriter(fluidfield, basename);
          xfluidwriter.WriteFiles();
        }
        break;
    }
    case prb_gas_fsi:
    case prb_ac_fsi:
    case prb_thermo_fsi:
    {
      std::string basename = problem.outname();
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      int numdisc = problem.num_discr();

      for (int i=0; i<numdisc-3; ++i)
      {
        PostField* scatrafield = problem.get_discretization(3+i);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case prb_biofilm_fsi:
    {
      std::string basename = problem.outname();
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      int numdisc = problem.num_discr();

      for (int i=0; i<numdisc-4; ++i)
      {
        PostField* scatrafield = problem.get_discretization(3+i);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case prb_struct_ale:
    {
      PostField* structurefield = problem.get_discretization(0);
      StructureFilter structwriter(structurefield, problem.outname(), problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* alefield = problem.get_discretization(1);
      AleFilter alewriter(alefield, problem.outname());
      alewriter.WriteFiles();
      break;
    }
    case prb_structure:
    case prb_invana:
    {
        PostField* field = problem.get_discretization(0);
        StructureFilter writer(field, problem.outname(), problem.stresstype(), problem.straintype());
        writer.WriteFiles();
        break;
    }
    case prb_fluid:
    {
      if (problem.num_discr()== 2)
      {
        if (problem.get_discretization(1)->name()=="xfluid")
        {
          std::string basename = problem.outname();
          PostField* fluidfield = problem.get_discretization(1);
          XFluidFilter xfluidwriter(fluidfield, basename);
          xfluidwriter.WriteFiles();
        }
      }
      // don't add break here !!!!
    }
    case prb_fluid_redmodels:
    {
      if (problem.num_discr()== 2)
      {
        // 1d artery
#ifdef D_ARTNET
        std::string basename = problem.outname();
        PostField* field = problem.get_discretization(1);
        StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
        writer.WriteFiles();
#endif //D_ARTNET
        if (problem.get_discretization(1)->name()=="xfluid")
        {
          std::string basename = problem.outname();
          PostField* fluidfield = problem.get_discretization(1);
          XFluidFilter xfluidwriter(fluidfield, basename);
          xfluidwriter.WriteFiles();
        }
      }
      // don't add break here !!!!
    }
    case prb_fluid_ale:
    case prb_freesurf:
    {
      PostField* field = problem.get_discretization(0);
      FluidFilter writer(field, problem.outname());
      writer.WriteFiles();
      if (problem.num_discr()>1 and problem.get_discretization(1)->name()=="xfluid")
      {
        std::string basename = problem.outname();
        PostField* fluidfield = problem.get_discretization(1);
        FluidFilter xfluidwriter(fluidfield, basename);
        xfluidwriter.WriteFiles();
      }
      break;
    }
    case prb_particle:
    {
      switch(problem.num_discr())
      {
        case 1:
        {
          PostField* particlefield = problem.get_discretization(0);
          ParticleFilter  particlewriter(particlefield, problem.outname());
          particlewriter.WriteFiles();

          break;
        }

        case 3:
        {
          PostField* structure = problem.get_discretization(0);
          StructureFilter writer(structure, problem.outname(), problem.stresstype(), problem.straintype());
          writer.WriteFiles();

          PostField* particlewallfield = problem.get_discretization(1);
          StructureFilter wallwriter(particlewallfield, problem.outname());
          wallwriter.WriteFiles();

          PostField* particlefield = problem.get_discretization(2);
          ParticleFilter  particlewriter(particlefield, problem.outname());
          particlewriter.WriteFiles();

          break;
        }

        default:
        {
          dserror("Particle problem has illegal number of discretizations!");
          break;
        }
      }

      break;
    }
    case prb_crack:
    {
      PostField* crackfield = problem.get_discretization(0);
      StructureFilter writer(crackfield, problem.outname(), problem.stresstype(), problem.straintype());
      writer.WriteFilesChangingGeom();
      break;
    }
    case prb_cavitation:
    {
      PostField* fluidfield = problem.get_discretization(0);
      FluidFilter fluidwriter(fluidfield, problem.outname());
      fluidwriter.WriteFiles();

      PostField* particlewallfield = problem.get_discretization(1);
      StructureFilter writer(particlewallfield, problem.outname());
      writer.WriteFiles();

      PostField* particlefield = problem.get_discretization(2);
      ParticleFilter  particlewriter(particlefield, problem.outname());
      particlewriter.WriteFiles();
      break;
    }
    case prb_level_set:
    {
      std::string basename = problem.outname();

      PostField* scatrafield = problem.get_discretization(0);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      // check if we have a particle field
      int numfield = problem.num_discr();
      if(numfield==2)
      {
        PostField* particlefield = problem.get_discretization(1);
        ParticleFilter  particlewriter(particlefield, basename);
        particlewriter.WriteFiles();
      }

      break;
    }
    case prb_redairways_tissue:
    {
        PostField* structfield = problem.get_discretization(0);
        StructureFilter structwriter(structfield, problem.outname(), problem.stresstype(), problem.straintype());
        structwriter.WriteFiles();

        PostField* fluidfield = problem.get_discretization(1);
        StructureFilter fluidwriter(fluidfield, problem.outname(), problem.stresstype(), problem.straintype());
        fluidwriter.WriteFiles();
        break;
    }
    case prb_ale:
    {
        PostField* field = problem.get_discretization(0);
        AleFilter writer(field, problem.outname());
        writer.WriteFiles();
        break;
    }
    case prb_lubrication:
    {
        PostField* lubricationfield = problem.get_discretization(0);
        LubricationFilter lubricationwriter(lubricationfield, problem.outname());
        lubricationwriter.WriteFiles();
        break;
    }
    case prb_porofluidmultiphase:
    {
        PostField* field = problem.get_discretization(0);
        PoroFluidMultiPhaseFilter writer(field, problem.outname());
        writer.WriteFiles();
        break;
    }
    case prb_poromultiphase:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      PoroFluidMultiPhaseFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();
      break;
    }
    case prb_poromultiphasescatra:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      PoroFluidMultiPhaseFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      PostField* scatrafield = problem.get_discretization(2);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();
      break;
    }
    case prb_scatra_endoexocytosis:
    case prb_cardiac_monodomain:
    case prb_scatra:
    {
        std::string basename = problem.outname();
        // do we have a fluid discretization?
        int numfield = problem.num_discr();
        if(numfield==2)
        {
          PostField* fluidfield = problem.get_discretization(0);
          FluidFilter fluidwriter(fluidfield, basename);
          fluidwriter.WriteFiles();

          PostField* scatrafield = problem.get_discretization(1);
          ScaTraFilter scatrawriter(scatrafield, basename);
          scatrawriter.WriteFiles();
        }
        else if (numfield==1)
        {
          PostField* scatrafield = problem.get_discretization(0);
          ScaTraFilter scatrawriter(scatrafield, basename);
          scatrawriter.WriteFiles();
        }
        else
          dserror("number of fields does not match: got %d",numfield);

        break;
    }
    case prb_sti:
    {
      // extract label for output files
      std::string basename = problem.outname();

      // safety check
      if(problem.num_discr() != 2)
        dserror("Must have exactly two discretizations for scatra-thermo interaction problems!");

      DRT::ELEMENTS::Transport* transport_element = dynamic_cast<DRT::ELEMENTS::Transport*>(problem.get_discretization(0)->discretization()->lRowElement(0));
      if(transport_element == NULL)
        dserror("Elements of unknown type on scalar transport discretization!");

      if(transport_element->ImplType() == INPAR::SCATRA::impltype_elch_electrode_thermo or
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
    case prb_fsi_xfem:
    case prb_fsi_crack:
    {
      std::cout << "Output FSI-XFEM Problem" << std::endl;

      int numfield = problem.num_discr();

      std::cout << "Number of discretizations: " << numfield << std::endl;
      for(int i=0; i<numfield;i++)
      {
        std::cout << "dis-name i=" << i << ": " << problem.get_discretization(i)->name() << std::endl;
      }


      std::string basename = problem.outname();

      std::cout << "\nStart postprocessing for discretizations: \n" << std::endl;

      std::cout << "  Structural Field ( "<< problem.get_discretization(0)->name() << " )" << std::endl;
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, problem.outname(), problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      std::cout << "  Fluid Field ( "<< problem.get_discretization(1)->name() << " )" << std::endl;
      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      // start index for interface discretizations
      int idx_int = 2;

      if (numfield > 2 && problem.get_discretization(2)->name()=="xfluid") // currently the name of the embedded fluid!
      {
        std::cout << "  XFluid Field ( "<< problem.get_discretization(2)->name() << " )" << std::endl;
        PostField* xfluidfield = problem.get_discretization(2);
        FluidFilter xfluidwriter(xfluidfield, basename);
        xfluidwriter.WriteFiles();
        idx_int += 1; // increase the index
      }

      // all other fields are interface or ale fields
      for(int i=idx_int; i<numfield;i++)
      {
        if (problem.get_discretization(i)->name() != "ale")
        {
          std::cout << "  Interface Field ( "<< problem.get_discretization(i)->name() << " )" << std::endl;
          PostField* ifacefield = problem.get_discretization(i);
          InterfaceFilter ifacewriter(ifacefield, basename);
          ifacewriter.WriteFiles();
        }
        else
        {
          std::cout << "  Ale Field ( "<< problem.get_discretization(i)->name() << " )" << std::endl;
          PostField* alefield = problem.get_discretization(i);
          AleFilter alewriter(alefield, basename);
          alewriter.WriteFiles();
        }
      }

      break;
    }
    case prb_fpsi_xfem:
    {
      int numfield = problem.num_discr();

      std::string basename = problem.outname();

      std::cout << "  Structural Field ( "<< problem.get_discretization(0)->name() << " )" << std::endl;
      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      std::cout << "  Porofluid Field ( "<< problem.get_discretization(1)->name() << " )" << std::endl;
      PostField* porofluidfield = problem.get_discretization(1);
      FluidFilter porofluidwriter(porofluidfield, basename);
      porofluidwriter.WriteFiles();

      std::cout << "  Ale Field ( "<< problem.get_discretization(2)->name() << " )" << std::endl;
      PostField* alefield = problem.get_discretization(2);
      AleFilter alewriter(alefield, basename);
      alewriter.WriteFiles();

      std::cout << "  Fluid Field ( "<< problem.get_discretization(3)->name() << " )" << std::endl;
      PostField* fluidfield = problem.get_discretization(3);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      // all other fields are interface fields
      for(int i=4; i<numfield;i++)
      {
        std::cout << "  Interface Field ( "<< problem.get_discretization(i)->name() << " )" << std::endl;
        PostField* ifacefield = problem.get_discretization(i);
        InterfaceFilter ifacewriter(ifacefield, basename);
        ifacewriter.WriteFiles();
      }

      break;
    }
    case prb_fluid_xfem:
    {
      std::cout << "Output FLUID-XFEM Problem" << std::endl;

      int numfield = problem.num_discr();

      std::cout << "Number of discretizations: " << numfield << std::endl;
      for(int i=0; i<numfield;i++)
      {
        std::cout << "dis-name i=" << i << ": " << problem.get_discretization(i)->name() << std::endl;
      }

      if (numfield == 0)
        dserror("we expect at least a fluid field, numfield=%i",numfield);
      std::string basename = problem.outname();

      // XFluid in the standard case, embedded fluid for XFF
      std::cout << "  Fluid Field" << std::endl;
      PostField* fluidfield = problem.get_discretization(0);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      // start index for interface discretizations
      int idx_int = 1;
      if (problem.num_discr()>1 and problem.get_discretization(1)->name()=="xfluid")
      {
        // XFluid for XFF
        std::cout << "  XFluid Field" << std::endl;
        PostField* xfluidfield = problem.get_discretization(1);
        FluidFilter xfluidwriter(xfluidfield, basename);
        xfluidwriter.WriteFiles();
        idx_int += 1;
      }

      // all other fields are interface fields
      for(int i=idx_int; i<numfield;i++)
      {
        std::cout << "  Interface Field ( "<< problem.get_discretization(i)->name() << " )" << std::endl;
        PostField* ifacefield = problem.get_discretization(i);
        InterfaceFilter ifacewriter(ifacefield, basename);
        ifacewriter.WriteFiles();
      }
      break;
    }
    case prb_loma:
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
    case prb_two_phase_flow:
    case prb_fluid_xfem_ls:
    {
      std::string basename = problem.outname();

      PostField* fluidfield = problem.get_discretization(0);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      PostField* scatrafield = problem.get_discretization(1);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      // check if we have a particle field
      int numfield = problem.num_discr();
      if(numfield==3)
      {
        PostField* particlefield = problem.get_discretization(2);
        ParticleFilter particlewriter(particlefield, basename);
        particlewriter.WriteFiles();
      }
      break;
    }
    case prb_elch:
    {
      std::string basename = problem.outname();
      int numfield = problem.num_discr();
      if(numfield==3)
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
      else if(numfield==2)
      {
        // Fluid and ScaTra fields are present
        PostField* fluidfield = problem.get_discretization(0);
        FluidFilter fluidwriter(fluidfield, basename);
        fluidwriter.WriteFiles();

        PostField* scatrafield = problem.get_discretization(1);
        ElchFilter elchwriter(scatrafield, basename);
        elchwriter.WriteFiles();
        break;
      }
      else if (numfield==1)
      {
        // only a ScaTra field is present
        PostField* scatrafield = problem.get_discretization(0);
        ElchFilter elchwriter(scatrafield, basename);
        elchwriter.WriteFiles();
      }
      else
        dserror("number of fields does not match: got %d",numfield);
      break;
    }
    case prb_combust:
    {
        std::string basename = problem.outname();

        PostField* fluidfield = problem.get_discretization(0);
        XFluidFilter fluidwriter(fluidfield, basename);
        fluidwriter.WriteFiles();

        PostField* scatrafield = problem.get_discretization(1);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();

        // check if we have a particle field
        int numfield = problem.num_discr();
        if(numfield==3)
        {
          PostField* particlefield = problem.get_discretization(2);
          ParticleFilter particlewriter(particlefield, basename);
          particlewriter.WriteFiles();
        }
        break;
    }
    case prb_art_net:
    {
      std::string basename = problem.outname();
      PostField* field = problem.get_discretization(0);
      //AnyFilter writer(field, problem.outname());
      StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
      writer.WriteFiles();

      break;
    }
    case prb_thermo:
    {
      PostField* field = problem.get_discretization(0);
      ThermoFilter writer(field, problem.outname(), problem.heatfluxtype(), problem.tempgradtype());
      writer.WriteFiles();
      break;
    }
    case prb_tsi:
    case prb_tfsi_aero:
    {
      std::cout << "Output TSI Problem" << std::endl;

      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* thermfield = problem.get_discretization(1);
      ThermoFilter thermwriter(thermfield, basename, problem.heatfluxtype(), problem.tempgradtype());
      thermwriter.WriteFiles();
      break;
    }
    case prb_red_airways:
    {
      std::string basename = problem.outname();
      PostField* field = problem.get_discretization(0);
      //AnyFilter writer(field, problem.outname());
      StructureFilter writer(field, basename, problem.stresstype(), problem.straintype());
      writer.WriteFiles();
      //      writer.WriteFiles();

      break;
    }
    case prb_poroelast:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();
      break;
    }
    case prb_poroscatra:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      PostField* scatrafield = problem.get_discretization(2);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      break;
    }
    case prb_fpsi:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* porofluidfield = problem.get_discretization(1);
      FluidFilter porofluidwriter(porofluidfield, basename);
      porofluidwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(2);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      break;
    }
    case prb_immersed_fsi:
    case prb_immersed_ale_fsi:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(1);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      break;
    }
    case prb_immersed_cell:
    {
      std::string basename = problem.outname();

      for(int field=0; field<problem.num_discr(); field++)
      {
        PostField* postfield = problem.get_discretization(field);
        std::cout<<"Write Field "<<field<<": "<<postfield->name()<<std::endl;
        if(postfield->name() == "cell" or postfield->name() == "structure")
        {
          StructureFilter cellwriter(postfield, basename, problem.stresstype(), problem.straintype());
          cellwriter.WriteFiles();
        }
        else if(postfield->name() == "cellscatra" or postfield->name() == "scatra")
        {
          ScaTraFilter scatrawriter(postfield, basename);
          scatrawriter.WriteFiles();
        }
        else if(postfield->name() == "ale")
        {
          AleFilter alewriter(postfield,basename);
          alewriter.WriteFiles();
        }
        else if(postfield->name() == "porofluid" or postfield->name() == "fluid")
        {
          FluidFilter fluidwriter(postfield, basename);
          fluidwriter.WriteFiles();
        }
        else
        {
          dserror("unknown field name");
        }
      }

      break;
    }
    case prb_fps3i:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* porofluidfield = problem.get_discretization(1);
      FluidFilter porofluidwriter(porofluidfield, basename);
      porofluidwriter.WriteFiles();

      PostField* fluidfield = problem.get_discretization(2);
      FluidFilter fluidwriter(fluidfield, basename);
      fluidwriter.WriteFiles();

      /////////////

      int numdisc = problem.num_discr();

      for (int i=0; i<numdisc-4; ++i)
      {
        PostField* scatrafield = problem.get_discretization(4+i);
        ScaTraFilter scatrawriter(scatrafield, basename);
        scatrawriter.WriteFiles();
      }

      break;
    }
    case prb_ehl:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* lubricationfield = problem.get_discretization(1);
      LubricationFilter lubricationwriter(lubricationfield, basename);
      lubricationwriter.WriteFiles();

      break;
    }
    case prb_ssi:
    {
      std::string basename = problem.outname();

      PostField* structfield = problem.get_discretization(0);
      StructureFilter structwriter(structfield, basename, problem.stresstype(), problem.straintype());
      structwriter.WriteFiles();

      PostField* scatrafield = problem.get_discretization(1);
      ScaTraFilter scatrawriter(scatrafield, basename);
      scatrawriter.WriteFiles();

      break;
    }
    case prb_fluid_topopt:
    {
      std::string basename = problem.outname();

      for (int i=0;i<problem.num_discr();i++)
      {
        std::string disname = problem.get_discretization(i)->discretization()->Name();

        if (disname.compare("fluid") == 0) // 0=true
        {
          PostField* fluidfield = problem.get_discretization(i);
          FluidFilter fluidwriter(fluidfield, basename);
          fluidwriter.WriteFiles();
        }
        else if (disname.compare("opti") == 0) // 0=true
        {
          PostField* optifield = problem.get_discretization(i);
          ScaTraFilter optiwriter(optifield, basename);
          optiwriter.WriteFiles();
        }
        else
          dserror("unknown discretization for postprocessing of topopt problem!");
      }

      break;
    }
    case prb_acou:
    {

      for (int i=0;i<problem.num_discr();i++)
      {
        std::string disname = problem.get_discretization(i)->discretization()->Name();
        if (disname.compare("acou") == 0) // 0=true
        {
          PostField* field = problem.get_discretization(i);
          AcouFilter writer(field, problem.outname());
          writer.WriteFiles();
        }
        else if(disname.compare("scatra") == 0)
        {
          PostField* field1 = problem.get_discretization(i);
          ScaTraFilter writer1(field1, problem.outname());
          writer1.WriteFiles();
        }
        else
          dserror("unknown discretization for postprocessing of acoustical problem!");
      }
      break;
    }
    case prb_uq:
    {
      for (int i=0;i<problem.num_discr();i++)
      {
        std::string disname = problem.get_discretization(i)->discretization()->Name();
        if( (disname.compare("structure") == 0) || (disname.compare("red_airway") == 0 ) ) // 0=true
        {
          PostField* field = problem.get_discretization(i);
          StructureFilter writer(field, problem.outname(), problem.stresstype(), problem.straintype());
          writer.WriteFiles();
        }
        else if(disname.compare("ale") == 0)
        {
          PostField* field = problem.get_discretization(i);
          AleFilter writer(field, problem.outname());
          writer.WriteFiles();
          break;
        }
        else
        {
           dserror("Unknown discretization type for problem type UQ");
        }
      } break;
    }
    case prb_none:
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
  std::string get_filter (int argc, char** argv)
  {
    Teuchos::CommandLineProcessor clp(false, false, false);
    std::string filter ("ensight");
    clp.setOption("filter",&filter,"filter to run [ensight, gid, vtu, vti]");
    // ignore warnings about unrecognized options.
    std::ostringstream warning;
    clp.parse(argc, argv, &warning);
    return filter;
  }
}


/*!
 \brief post-processor main routine

 Select the appropriate filter and run!

 \author kronbichler
 \date 03/14
 */
int main(
        int argc,
        char** argv)
{
  try
  {
    std::string filter = get_filter(argc, argv);
    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString("Main BACI post-processor\n");

    PostProblem problem(My_CLP, argc, argv);

    if (filter == "ensight" || filter == "vtu" || filter == "vti")
      runEnsightVtuFilter(problem);
    else if (filter == "gid")
      PostGid::runGidFilter(problem);
    else
      dserror("Unknown filter %s given, supported filters: [ensight|vtu|vti|gid]", filter.c_str());

  } // try
  catch ( std::runtime_error & err )
  {
    char line[] = "=========================================================================\n";
    std::cout << "\n\n"
        << line
        << err.what()
        << "\n"
        << line
        << "\n" << std::endl;

    // proper cleanup
    DRT::Problem::Done();
#ifdef DSERROR_DUMP
    abort();
#endif

#ifdef PARALLEL
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
    exit(1);
#endif
  } // catch

  // proper cleanup
  DRT::Problem::Done();

  return 0;
}

