/*----------------------------------------------------------------------*/
/*! \file

\brief Base algorithm for all kinds of coupled problems


\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_adapter_algorithmbase.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_pstream.hpp"
#include "4C_utils_exceptions.hpp"

#include <iostream>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::AlgorithmBase::AlgorithmBase(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : comm_(comm), printscreen_(GLOBAL::Problem::Instance()->IOParams().get<int>("STDOUTEVRY"))
{
  if (comm_.MyPID() == 0) INPUT::PrintDefaultParameters(IO::cout, timeparams);

  step_ = 0;
  time_ = 0.;
  dt_ = timeparams.get<double>("TIMESTEP");
  nstep_ = timeparams.get<int>("NUMSTEP");
  maxtime_ = timeparams.get<double>("MAXTIME");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::AlgorithmBase::SetTimeStep(const double time, const int step)
{
  step_ = step;
  time_ = time;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::AlgorithmBase::print_header()
{
  if (Comm().MyPID() == 0 and printscreen_ and (step_ % printscreen_ == 0))
  {
    IO::cout << "\n"
             << method_ << "\n"
             << "TIME:  " << std::scientific << time_ << "/" << std::scientific << maxtime_
             << "     DT = " << std::scientific << dt_ << "     STEP = " << std::setw(4) << step_
             << "/" << std::setw(4) << nstep_ << "\n\n";
  }
}

FOUR_C_NAMESPACE_CLOSE
