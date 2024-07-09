/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all fluid time integrations


\level 1

*/
/*-----------------------------------------------------------*/

#include "4C_fluid_timint.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fluid.hpp"
#include "4C_io_discretization_visualization_writer_mesh.hpp"
#include "4C_io_visualization_parameters.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN

FLD::TimInt::TimInt(const Teuchos::RCP<Core::FE::Discretization>& discret,
    const Teuchos::RCP<Core::LinAlg::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<Core::IO::DiscretizationWriter>& output)
    : discret_(discret),
      solver_(solver),
      params_(params),
      output_(output),
      runtime_output_writer_(Teuchos::null),
      runtime_output_params_(),
      time_(0.0),
      step_(0),
      dta_(params_->get<double>("time step size")),
      stepmax_(params_->get<int>("max number timesteps")),
      maxtime_(params_->get<double>("total time")),
      itemax_(params_->get<int>("max nonlin iter steps")),
      uprestart_(params_->get("write restart every", -1)),
      upres_(params_->get("write solution every", -1)),
      timealgo_(
          Core::UTILS::GetAsEnum<Inpar::FLUID::TimeIntegrationScheme>(*params_, "time int algo")),
      physicaltype_(Core::UTILS::GetAsEnum<Inpar::FLUID::PhysicalType>(*params_, "Physical Type")),
      myrank_(discret_->get_comm().MyPID()),
      updateprojection_(false),
      projector_(Teuchos::null),
      kspsplitter_(Teuchos::null)
{
  // check for special fluid output which is to be handled by an own writer object
  Teuchos::RCP<const Teuchos::ParameterList> fluid_runtime_output_list =
      Teuchos::rcp(new Teuchos::ParameterList(
          Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT").sublist("FLUID")));

  bool output_fluid =
      (bool)Core::UTILS::IntegralValue<int>(*fluid_runtime_output_list, "OUTPUT_FLUID");

  // create and initialize parameter container object for fluid specific runtime output
  if (output_fluid)
  {
    runtime_output_params_.init(*fluid_runtime_output_list);
    runtime_output_params_.setup();

    // TODO This does not work for restarted simulations as the time_ is not yet correctly set.
    // However, this is called before the restart is read and someone with knowledge on the module
    // has to refactor the code. The only implication is that in restarted simulations the .pvd file
    // does not contain the steps of the simulation that is restarted from
    runtime_output_writer_ = Teuchos::rcp(new Core::IO::DiscretizationVisualizationWriterMesh(
        discret_, Core::IO::VisualizationParametersFactory(
                      Global::Problem::instance()->io_params().sublist("RUNTIME VTK OUTPUT"),
                      *Global::Problem::instance()->output_control_file(), time_)));
  }
}

Teuchos::RCP<const Epetra_Map> FLD::TimInt::dof_row_map(unsigned nds)
{
  return Teuchos::rcp(discretization()->dof_row_map(nds), false);
}


void FLD::TimInt::increment_time_and_step()
{
  step_ += 1;
  time_ += dta_;

  return;
}

FOUR_C_NAMESPACE_CLOSE
