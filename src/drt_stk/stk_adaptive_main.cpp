
#include "stk_adaptive_main.H"
#include "stk_fluid.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../stk_lib/stk_discret.H"
#include "../stk_fluid/fluid_implicit.H"

#include "../linalg/linalg_solver.H"

#include "../drt_io/io_control.H"

void adaptive_main()
{
#ifdef STKADAPTIVE
  switch (genprob.probtyp)
  {
    case prb_structure:
      switch (genprob.timetyp)
      {
        case time_static:
          /* nonlinear statics with new discretization */
          //break;
        case time_dynamic:
          //break;
        default:
          dserror("Unspecified time handling");
      }
      break;

    case prb_struct_multi:
    {
      switch (genprob.timetyp)
      {
        case time_dynamic:
          //break;
        case time_static:
          dserror("structural multi-scale algorithm only implemented for dynamic problems");
        default:
          dserror("Unspecified time handling");
      }
    }
    break;

    case prb_fluid_pm:
    case prb_fluid:
    {
      Teuchos::RCP<DRT::Discretization> actdis = null;
      actdis = DRT::Problem::Instance()->Dis( genprob.numff, 0 );
      if ( not actdis->HaveDofs() )
      {
        actdis->FillComplete();
      }

      // -------------------------------------------------------------------
      // create a solver
      // -------------------------------------------------------------------
      Teuchos::RCP<LINALG::Solver> solver =
        rcp(new LINALG::Solver(DRT::Problem::Instance()->FluidSolverParams(),
                               actdis->Comm(),
                               DRT::Problem::Instance()->ErrorFile()->Handle()));
      actdis->ComputeNullSpaceIfNecessary(solver->Params());

#if 1

      STK::Discretization dis( actdis->Comm() );
      STK::FLD::Fluid fluid( dis, solver );

      dis.Setup( *actdis, fluid );

#else
      STK::Fluid fluid( *actdis, solver );
      fluid.SetupSTKMesh();

      // full refinement
      fluid.RefineAll();
      fluid.RefineAll();

      //fluid.RefineHalf();

      //fluid.RefineFirst();

      fluid.SetupDRTMesh();

      fluid.Integrate();

      DRT::Problem::Instance()->AddFieldTest(fluid.CreateFieldTest());
      DRT::Problem::Instance()->TestAll(actdis->Comm());
#endif

      break;
    }
    case prb_scatra:
      //break;
    case prb_fluid_xfem:
      //break;
    case prb_fluid_ale:
      //break;
    case prb_freesurf:
      //break;
    case prb_fsi:
    case prb_pfsi:
    case prb_fsi_lung:
      //break;
    case prb_fsi_xfem:
      //break;
    case prb_ale:
    case prb_thermo:
      //break;
    case prb_tsi:
      //break;
    case prb_loma:
      //break;
    case prb_elch:
      //break;
    case prb_combust:
      //break;
    case prb_art_net:
      //break;
    case prb_red_airways:
      //break;

    default:
      dserror("solution of unknown problemtyp %d requested", genprob.probtyp);
      break;
  }
#endif
}
