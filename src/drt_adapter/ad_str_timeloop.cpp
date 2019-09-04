/*----------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

\maintainer Anh-Tu Vuong

\level 1

*/
/*----------------------------------------------------------------------*/

#include "ad_str_timeloop.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimeLoop::Integrate()
{
  // error checking variables
  INPAR::STR::ConvergenceStatus convergencestatus = INPAR::STR::conv_success;

  // target time #timen_ and step #stepn_ already set
  // time loop
  while (NotFinished() and (convergencestatus == INPAR::STR::conv_success or
                               convergencestatus == INPAR::STR::conv_fail_repeat))
  {
    // call the predictor
    PrePredict();
    PrepareTimeStep();
    PostPredict();

    // integrate time step, i.e. do corrector steps
    // after this step we hold disn_, etc
    PreSolve();
    convergencestatus = Solve();
    PostSolve();

    // if everything is fine
    if (convergencestatus == INPAR::STR::conv_success)
    {
      // calculate stresses, strains and energies
      // note: this has to be done before the update since otherwise a potential
      // material history is overwritten
      PrepareOutput();

      // update displacements, velocities, accelerations
      // after this call we will have disn_==dis_, etc
      // update time and step
      // update everything on the element level
      PreUpdate();
      Update();
      PostUpdate();

      // write output
      PreOutput();
      Output();
      PostOutput();

      // print info about finished time step
      PrintStep();
    }
    // todo: remove this as soon as old structure time integration is gone
    else if (DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(
                 DRT::Problem::Instance()->StructuralDynamicParams(), "INT_STRATEGY") ==
             INPAR::STR::int_old)
    {
      convergencestatus =
          PerformErrorAction(convergencestatus);  // something went wrong update error code
                                                  // according to chosen divcont action
    }
  }

  PostTimeLoop();

  // that's it say what went wrong
  return convergencestatus;
}
