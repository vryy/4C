/*!----------------------------------------------------------------------
\file fluid_timint.cpp
\brief Baseclass for all fluid time integrations

\level 1

<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#include "fluid_timint.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_inpar/inpar_parameterlist_utils.H"
#include "../drt_inpar/inpar_fluid.H"

#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

FLD::TimInt::TimInt(
  const Teuchos::RCP<DRT::Discretization>&      discret,
  const Teuchos::RCP<LINALG::Solver>&           solver,
  const Teuchos::RCP<Teuchos::ParameterList>&   params,
  const Teuchos::RCP<IO::DiscretizationWriter>& output
):
  discret_     (discret),
  solver_      (solver),
  params_      (params),
  output_      (output),
  time_        (0.0),
  step_        (0),
  dta_         (params_->get<double>("time step size")),
  stepmax_     (params_->get<int>   ("max number timesteps")),
  maxtime_     (params_->get<double>("total time")),
  itemax_      (params_->get<int>("max nonlin iter steps")),
  uprestart_   (params_->get("write restart every", -1)),
  upres_       (params_->get("write solution every", -1)),
  timealgo_    (DRT::INPUT::get<INPAR::FLUID::TimeIntegrationScheme>(*params_, "time int algo")),
  physicaltype_(DRT::INPUT::get<INPAR::FLUID::PhysicalType>(*params_, "Physical Type")),
  myrank_      (discret_->Comm().MyPID()),
  updateprojection_(false),
  projector_(Teuchos::null),
  kspsplitter_(Teuchos::null)
{}

FLD::TimInt::~TimInt()
{
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
