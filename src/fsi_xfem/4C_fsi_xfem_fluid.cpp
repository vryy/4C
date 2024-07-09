/*----------------------------------------------------------------------*/
/*! \file

\brief ...

\level 2

*/

#include "4C_fsi_xfem_fluid.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_pstream.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidXFEMAlgorithm::FluidXFEMAlgorithm(const Epetra_Comm& comm)
    : FluidMovingBoundaryBaseAlgorithm(
          Global::Problem::instance()->fluid_dynamic_params(), "FSICoupling"),
      comm_(comm)
{
  const Teuchos::ParameterList& fluiddyn = Global::Problem::instance()->fluid_dynamic_params();

  if (comm_.MyPID() == 0) Input::PrintDefaultParameters(Core::IO::cout, fluiddyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fluiddyn.get<double>("TIMESTEP");
  nstep_ = fluiddyn.get<int>("NUMSTEP");
  maxtime_ = fluiddyn.get<double>("MAXTIME");
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::timeloop()
{
  if (Global::Problem::instance()->get_problem_type() == Core::ProblemType::fluid_xfem)
  {
    if (get_comm().MyPID() == 0)
      std::cout << "Integrate routine for MOVING INTERFACES"
                << "\n"
                << std::endl;


    while (not_finished())
    {
      prepare_time_step();
      solve();
      update();
      output();
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::read_restart(int step)
{
  time_ = mb_fluid_field()->read_restart(step);
  step_ = step;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;


  mb_fluid_field()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::solve() { mb_fluid_field()->nonlinear_solve(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::update() { mb_fluid_field()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::output() { mb_fluid_field()->output(); }

FOUR_C_NAMESPACE_CLOSE
