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
          Global::Problem::Instance()->FluidDynamicParams(), "FSICoupling"),
      comm_(comm)
{
  const Teuchos::ParameterList& fluiddyn = Global::Problem::Instance()->FluidDynamicParams();

  if (comm_.MyPID() == 0) Input::PrintDefaultParameters(Core::IO::cout, fluiddyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fluiddyn.get<double>("TIMESTEP");
  nstep_ = fluiddyn.get<int>("NUMSTEP");
  maxtime_ = fluiddyn.get<double>("MAXTIME");
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::Timeloop()
{
  if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::fluid_xfem)
  {
    if (Comm().MyPID() == 0)
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
  time_ = MBFluidField()->read_restart(step);
  step_ = step;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::prepare_time_step()
{
  step_ += 1;
  time_ += dt_;


  MBFluidField()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::solve() { MBFluidField()->nonlinear_solve(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::update() { MBFluidField()->update(); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidXFEMAlgorithm::output() { MBFluidField()->output(); }

FOUR_C_NAMESPACE_CLOSE
