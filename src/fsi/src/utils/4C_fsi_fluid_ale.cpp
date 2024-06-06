/*----------------------------------------------------------------------*/
/*! \file

\brief Solve fluid problems on ALE mesh

\level 2

*/
/*----------------------------------------------------------------------*/

#include "4C_fsi_fluid_ale.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAleAlgorithm::FluidAleAlgorithm(const Epetra_Comm& comm)
    : FluidMovingBoundaryBaseAlgorithm(
          Global::Problem::Instance()->FSIDynamicParams(), "FSICoupling"),
      comm_(comm)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::Instance()->FSIDynamicParams();

  if (comm_.MyPID() == 0) Input::PrintDefaultParameters(Core::IO::cout, fsidyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::Timeloop()
{
  while (not_finished())
  {
    prepare_time_step();
    solve();
    update();
    output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::read_restart(int step)
{
  time_ = MBFluidField()->read_restart(step);
  step_ = step;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;

  if (Comm().MyPID() == 0)
    std::cout << "\n"
              << "TIME:  " << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_ << "     STEP = " << std::setw(4) << step_
              << "/" << std::setw(4) << nstep_ << "\n\n";

  MBFluidField()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::solve() { MBFluidField()->nonlinear_solve(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::update() { MBFluidField()->Update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::output() { MBFluidField()->Output(); }

FOUR_C_NAMESPACE_CLOSE
