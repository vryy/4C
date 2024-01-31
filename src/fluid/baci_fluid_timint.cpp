/*-----------------------------------------------------------*/
/*! \file

\brief Base class for all fluid time integrations


\level 1

*/
/*-----------------------------------------------------------*/

#include "baci_fluid_timint.H"

#include "baci_fluid_utils_mapextractor.H"
#include "baci_global_data.H"
#include "baci_inpar_fluid.H"
#include "baci_inpar_parameterlist_utils.H"
#include "baci_io_discretization_visualization_writer_mesh.H"
#include "baci_io_visualization_parameters.H"
#include "baci_lib_discret.H"

#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

FLD::TimInt::TimInt(const Teuchos::RCP<DRT::Discretization>& discret,
    const Teuchos::RCP<CORE::LINALG::Solver>& solver,
    const Teuchos::RCP<Teuchos::ParameterList>& params,
    const Teuchos::RCP<IO::DiscretizationWriter>& output)
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
      timealgo_(INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo")),
      physicaltype_(INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type")),
      myrank_(discret_->Comm().MyPID()),
      updateprojection_(false),
      projector_(Teuchos::null),
      kspsplitter_(Teuchos::null)
{
  // check for special fluid output which is to be handled by an own writer object
  Teuchos::RCP<const Teuchos::ParameterList> fluid_runtime_output_list =
      Teuchos::rcp(new Teuchos::ParameterList(
          GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT").sublist("FLUID")));

  bool output_fluid = (bool)INPUT::IntegralValue<int>(*fluid_runtime_output_list, "OUTPUT_FLUID");

  // create and initialize parameter container object for fluid specific runtime output
  if (output_fluid)
  {
    runtime_output_params_.Init(*fluid_runtime_output_list);
    runtime_output_params_.Setup();
    runtime_output_writer_ = Teuchos::rcp(new IO::DiscretizationVisualizationWriterMesh(
        discret_, IO::VisualizationParametersFactory(
                      GLOBAL::Problem::Instance()->IOParams().sublist("RUNTIME VTK OUTPUT"))));
  }
}

Teuchos::RCP<const Epetra_Map> FLD::TimInt::DofRowMap(unsigned nds)
{
  return Teuchos::rcp(Discretization()->DofRowMap(nds), false);
}


void FLD::TimInt::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dta_;

  return;
}

BACI_NAMESPACE_CLOSE
