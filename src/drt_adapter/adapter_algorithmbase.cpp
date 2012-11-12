


#include <iostream>

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"

#include "adapter_algorithmbase.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::AlgorithmBase::AlgorithmBase(const Epetra_Comm& comm,
                                      const Teuchos::ParameterList& timeparams)
  : comm_(comm),
    printscreen_(DRT::Problem::Instance()->IOParams().get<int>("STDOUTEVRY"))
{
  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, timeparams);

  step_ = 0;
  time_ = 0.;
  dt_ = timeparams.get<double>("TIMESTEP");
  nstep_ = timeparams.get<int>("NUMSTEP");
  maxtime_ = timeparams.get<double>("MAXTIME");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::AlgorithmBase::SetTimeStep(double time, int step)
{
  step_ = step;
  time_ = time;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::AlgorithmBase::PrintHeader()
{
  if (Comm().MyPID()==0 and printscreen_ and (step_%printscreen_==0))
    std::cout << "\n"
              << method_ << "\n"
              << "TIME:  "    << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_
              << "     STEP = " YELLOW_LIGHT << std::setw(4) << step_ << END_COLOR "/" << std::setw(4) << nstep_
//               << "\n"
//               << NOX::Utils::fill(82)
              << "\n\n";
}

