/*----------------------------------------------------------------------*/
/*!
\file ad_str_timeloop.cpp

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_str_timeloop.H"
#include "../drt_inpar/inpar_structure.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int ADAPTER::StructureTimeLoop::Integrate()
{
  // error checking variables
  INPAR::STR::ConvergenceStatus convergencestatus = INPAR::STR::conv_success;

  // target time #timen_ and step #stepn_ already set
  // time loop
  while ( NotFinished() and (convergencestatus == INPAR::STR::conv_success) )
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
    if(convergencestatus == INPAR::STR::conv_success)
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
    else // something went wrong update error code according to chosen divcont action
    {
      convergencestatus = PerformErrorAction(convergencestatus);
    }
  }

  PostTimeLoop();

  // that's it say what went wrong
  return convergencestatus;
}

