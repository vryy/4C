/*!------------------------------------------------------------------------------------------------*
 \file ssi_utils.cpp

 \brief Utility methods for SSI

\level 1

\maintainer Julia Hoermann
            hoermann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264

 *------------------------------------------------------------------------------------------------*/

#include "ssi_utils.H"

#include <Epetra_Map.h>

#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_particle/binning_strategy.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 09/2014 */
/* Function for checking that the different time steps are a
 multiplicative of each other                                           */

int SSI::Utils::CheckTimeStepping(double dt1, double dt2)
{
  double workdt1 = std::min(dt1, dt2);
  double workdt2 = std::max(dt1, dt2);
  double t1 = 0.0;
  int i=0;

  while(true)
  {  i++;
     t1 = i* workdt1;

    if (std::abs(t1-workdt2) < 10E-10)
      break;

    else
      if (t1 > workdt2)
        dserror("Chosen time steps %f and %f are not a multiplicative of each other", dt1, dt2);
  }
  return i;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
/*                                                        AN,JH 10/2014 */
//Modification of time parameter list for problem with different time step size

void SSI::Utils::ChangeTimeParameter(const Epetra_Comm& comm,
    Teuchos::ParameterList& ssiparams,
    Teuchos::ParameterList& scatradyn,
    Teuchos::ParameterList& sdyn)
{
  bool difftimestep = DRT::INPUT::IntegralValue<int>(ssiparams, "DIFFTIMESTEPSIZE");

  if (difftimestep) //Create subproblems with different time steps
  {
    // Check correct choice of time stepping for single fields
    double scatrastep = scatradyn.get<double>("TIMESTEP");
    double solidstep  = sdyn.get<double>("TIMESTEP");

    SSI::Utils::CheckTimeStepping(scatrastep, solidstep);

    //modify global time step size
    ssiparams.set<double>   ("TIMESTEP"    ,std::min(scatrastep,solidstep));
  }
  else
  {
    // -------------------------------------------------------------------
    // overrule certain parameters for coupled problems
    // -------------------------------------------------------------------
    // the default time step size
    scatradyn.set<double>   ("TIMESTEP"    ,ssiparams.get<double>("TIMESTEP"));
    sdyn.set<double>   ("TIMESTEP"    ,ssiparams.get<double>("TIMESTEP"));
    // maximum simulation time
    scatradyn.set<double>   ("MAXTIME"     ,ssiparams.get<double>("MAXTIME"));
    sdyn.set<double>   ("MAXTIME"     ,ssiparams.get<double>("MAXTIME"));
    // maximum number of timesteps
    scatradyn.set<int>      ("NUMSTEP"     ,ssiparams.get<int>("NUMSTEP"));
    sdyn.set<int>      ("NUMSTEP"     ,ssiparams.get<int>("NUMSTEP"));
  }

  // Check correct input of restart. Code relies that both time value RESTARTEVRYTIME and RESULTSEVRYTIME are
  // given if restart from time is applied
  double restarttime = ssiparams.get<double>("RESTARTEVRYTIME");
  double updatetime  = ssiparams.get<double>("RESULTSEVRYTIME");
  if ((updatetime > 0.0) or (restarttime > 0.0))
    if (!(updatetime > 0.0) and !(restarttime > 0.0))
      dserror("If time controlled output and restart is desired, both parameters RESTARTEVRYTIME and RESULTSEVRYTIME has to be set");

  // set restart params
  int scatrarestart;
  int structurerestart;

  if  (restarttime > 0.0)
  {
    scatrarestart    = SSI::Utils::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), restarttime);
    structurerestart = SSI::Utils::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), restarttime);
  }
  else
  {
    int restart        = ssiparams.get<int>("RESTARTEVRY");
    scatrarestart    = restart;
    structurerestart = restart;

  }

  // set output params
  int scatraupres;
  int structureupres;

  if  (updatetime > 0.0)
  {
    scatraupres      = SSI::Utils::CheckTimeStepping(scatradyn.get<double>("TIMESTEP"), updatetime);
    structureupres   = SSI::Utils::CheckTimeStepping(sdyn.get<double>("TIMESTEP"), updatetime);
  }
  else
  {
    int update       = ssiparams.get<int>("RESULTSEVRY");
    scatraupres      = update;
    structureupres   = update;
  }

  // restart
  scatradyn.set<int> ("RESTARTEVRY" ,scatrarestart);
  sdyn.set<int>      ("RESTARTEVRY" ,structurerestart);
  // solution output
  scatradyn.set<int> ("RESULTSEVRY"       ,scatraupres);
  sdyn.set<int>      ("RESULTSEVRY" ,structureupres);

  if (comm.MyPID() == 0)
  {
    std::cout
        << "====================== Overview of chosen time stepping: ==============================\n"
        << "\t Timestep scatra:           "<< scatradyn.get<double>("TIMESTEP") << "\n"
        << "\t Timestep structure:        "<< sdyn.get<double>("TIMESTEP") << "\n"
        << "\t Result step scatra:        "<< scatradyn.get<int>("RESULTSEVRY") << "\n"
        << "\t Result step structure:     "<< sdyn.get<int>("RESULTSEVRY") << "\n"
        << "\t Restart step scatra:       "<< scatradyn.get<int>("RESTARTEVRY") << "\n"
        << "\t Restart step structure:    "<< sdyn.get<int>("RESTARTEVRY") << "\n"
        << "========================================================================================\n \n";
  }
}


/*----------------------------------------------------------------------*
 |  Redistribute using BinningStrategy                      rauch 08/16 |
 *----------------------------------------------------------------------*/
void SSI::Utils::RedistributeDiscretizationsByBinning(const std::vector<Teuchos::RCP<DRT::Discretization> > vector_of_discretizations)
{
  // safety check
  if (vector_of_discretizations.size()==0)
    dserror("No discretizations provided for binning !");

  // get communicator
  const Epetra_Comm& comm = vector_of_discretizations[0]->Comm();

  // redistribute discr. with help of binning strategy
  if(comm.NumProc()>1)
  {
    if(comm.MyPID() == 0)
    {
      std::cout<<"+---------------------------------------------------------------"<<std::endl;
      std::cout<<"| Redistribute discretizations using Binning Strategy ...                                   "<<std::endl;
      // fist we need to call FillComplete on all discretizations
      for(int dis_num=0; dis_num < (int)(vector_of_discretizations.size()); dis_num++)
      {
        std::cout<<"| Redistribute discretization "<<vector_of_discretizations[dis_num]->Name()<<std::endl;
      }
      std::cout<<"+---------------------------------------------------------------"<<std::endl;
    }

    // fist we need to call FillComplete on all discretizations
    for(int dis_num=0; dis_num < (int)(vector_of_discretizations.size()); dis_num++)
    {
      vector_of_discretizations[dis_num]->FillComplete(false,false,false);
    }

    // create vector of maps
    std::vector<Teuchos::RCP<Epetra_Map> > stdelecolmap;
    std::vector<Teuchos::RCP<Epetra_Map> > stdnodecolmap;

    // binning strategy is created and parallel redistribution is performed
    Teuchos::RCP<BINSTRATEGY::BinningStrategy> binningstrategy =
        Teuchos::rcp(new BINSTRATEGY::BinningStrategy(vector_of_discretizations,stdelecolmap,stdnodecolmap));
  }

  return;
}
