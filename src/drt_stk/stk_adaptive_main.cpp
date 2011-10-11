
#include "stk_adaptive_main.H"

#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_dserror.H"

#include "../stk_structure/str_problem.H"
#include "../stk_fluid/fluid_problem.H"
#include "../stk_fsi/fsi_problem.H"

void adaptive_main()
{
#ifdef STKADAPTIVE
  switch (genprob.probtyp)
  {
    case prb_structure:
      switch (genprob.timetyp)
      {
        case time_static:
          dserror( "nonlinear statics with new discretization" );
          break;
        case time_dynamic:
        {
          STK::STR::Problem prb;
          prb.Execute();
          break;
        }
        default:
          dserror("Unspecified time handling");
      }
      break;

    case prb_fluid_pm:
    case prb_fluid:
    {
      STK::FLD::Problem prb;
      prb.Execute();
      break;
    }
    case prb_fsi:
    {
      STK::FSI::Problem prb;
      prb.Execute();
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
    case prb_pfsi:
    case prb_fsi_lung:
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
#else
  dserror("Adaptive mesh is not available (flag STKADAPTIVE not defined)");
#endif
}
