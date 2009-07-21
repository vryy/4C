
#ifdef CCADISCRET

#include <iostream>

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "fsi_algorithmbase.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::AlgorithmBase::AlgorithmBase(Epetra_Comm& comm,
                                  const Teuchos::ParameterList& timeparams)
  : comm_(comm)
{
  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, timeparams);

  step_ = 0;
  time_ = 0.;
  dt_ = timeparams.get<double>("TIMESTEP");
  nstep_ = timeparams.get<int>("NUMSTEP");
  maxtime_ = timeparams.get<double>("MAXTIME");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::AlgorithmBase::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::AlgorithmBase::SetTimeStep(double time, int step)
{
  step_ = step;
  time_ = time;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::AlgorithmBase::PrintHeader()
{
  if (Comm().MyPID()==0)
    std::cout << "\n"
              << method_ << "\n"
              << "TIME:  "    << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_
              << "     STEP = " YELLOW_LIGHT << std::setw(4) << step_ << END_COLOR "/" << std::setw(4) << nstep_
//               << "\n"
//               << NOX::Utils::fill(82)
              << "\n\n";
}

#endif
