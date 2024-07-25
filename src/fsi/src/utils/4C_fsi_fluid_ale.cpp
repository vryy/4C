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
          Global::Problem::instance()->fsi_dynamic_params(), "FSICoupling"),
      comm_(comm)
{
  const Teuchos::ParameterList& fsidyn = Global::Problem::instance()->fsi_dynamic_params();

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::timeloop()
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
  time_ = mb_fluid_field()->read_restart(step);
  step_ = step;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;

  if (get_comm().MyPID() == 0)
    std::cout << "\n"
              << "TIME:  " << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_ << "     STEP = " << std::setw(4) << step_
              << "/" << std::setw(4) << nstep_ << "\n\n";

  mb_fluid_field()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::solve() { mb_fluid_field()->nonlinear_solve(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::update() { mb_fluid_field()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAleAlgorithm::output() { mb_fluid_field()->output(); }

FOUR_C_NAMESPACE_CLOSE
