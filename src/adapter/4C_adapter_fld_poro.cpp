/*----------------------------------------------------------------------*/
/*! \file

 \brief Fluid field adapter for poroelasticity




\level 2

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_poro.hpp"

#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_lib_assemblestrategy.hpp"
#include "4C_lib_condition_utils.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"

FOUR_C_NAMESPACE_OPEN

/*======================================================================*/
/* constructor */
ADAPTER::FluidPoro::FluidPoro(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : ADAPTER::FluidFPSI::FluidFPSI(fluid, dis, solver, params, output, isale, dirichletcond)
{
  // make sure

  if (fluid_ == Teuchos::null) FOUR_C_THROW("Failed to create the underlying fluid adapter");

  Discretization()->GetCondition("NoPenetration", nopencond_);
}

/*======================================================================*/
/* evaluate poroelasticity specific constraint*/
void ADAPTER::FluidPoro::EvaluateNoPenetrationCond(Teuchos::RCP<Epetra_Vector> Cond_RHS,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ConstraintMatrix,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> StructVelConstraintMatrix,
    Teuchos::RCP<Epetra_Vector> condVector, Teuchos::RCP<std::set<int>> condIDs,
    POROELAST::Coupltype coupltype)
{
  if (!(Discretization()->Filled())) FOUR_C_THROW("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) FOUR_C_THROW("AssignDegreesOfFreedom() was not called");

  Discretization()->SetState(0, "dispnp", Dispnp());
  Discretization()->SetState(0, "scaaf", Scaaf());

  Teuchos::ParameterList params;

  params.set("timescale", TimeScaling());

  if (coupltype == POROELAST::fluidfluid)
  {
    // first, find out which dofs will be constraint
    params.set<int>("action", FLD::no_penetrationIDs);
    params.set<int>("Physical Type", INPAR::FLUID::poro);
    Discretization()->EvaluateCondition(params, condVector, "NoPenetration");

    // write global IDs of dofs on which the no penetration condition is applied (can vary in time
    // and iteration)
    {
      const int ndim = GLOBAL::Problem::Instance()->NDim();
      const int ndof = ndim + 1;
      const int length = condVector->MyLength();
      const int nnod = length / ndof;
      const Epetra_BlockMap& map = condVector->Map();
      bool isset = false;
      for (int i = 0; i < nnod; i++)
      {
        isset = false;
        for (int j = 0; j < ndof; j++)
        {
          if ((*condVector)[i * ndof + j] != 0.0 and isset == false)
          {
            condIDs->insert(map.GID(i * ndof + j));
            isset = true;
            // break;
          }
          else
            (*condVector)[i * ndof + j] = 0.0;
        }
      }
    }

    // set action for elements
    params.set<int>("action", FLD::no_penetration);
    // params.set<Teuchos::RCP< std::set<int> > >("condIDs",condIDs);
    params.set<POROELAST::Coupltype>("coupling", POROELAST::fluidfluid);
    params.set<int>("Physical Type", INPAR::FLUID::poro);

    DRT::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                  // fluiddofset for column
        ConstraintMatrix,                   // fluid matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

    Discretization()->SetState(0, "condVector", condVector);

    Discretization()->EvaluateCondition(params, fluidstrategy, "NoPenetration");
  }
  else if (coupltype == POROELAST::fluidstructure)
  {
    Discretization()->SetState(0, "velnp", Velnp());
    Discretization()->SetState(0, "gridv", GridVel());

    Discretization()->SetState(0, "condVector", condVector);

    // set action for elements
    params.set<int>("action", FLD::no_penetration);
    params.set<POROELAST::Coupltype>("coupling", POROELAST::fluidstructure);
    params.set<int>("Physical Type", INPAR::FLUID::poro);

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy couplstrategy(0,  // fluiddofset for row
        1,                                  // structdofset for column
        ConstraintMatrix,
        StructVelConstraintMatrix,  // fluid-mechanical matrix
        Cond_RHS, Teuchos::null, Teuchos::null);

    // evaluate the fluid-mechanical system matrix on the fluid element
    Discretization()->EvaluateCondition(params, couplstrategy, "NoPenetration");
  }
  else
    FOUR_C_THROW("unknown coupling type for no penetration BC");

  Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<CORE::LINALG::MapExtractor> ADAPTER::FluidPoro::VelPresSplitter()
{
  return fluidimpl_->VelPresSplitter();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidPoro::Output(const int step, const double time)
{
  // set variables that allow the forced
  // output at this point.

  // we surely want to write now
  int upres_ = 1;
  // poro is always ALE
  bool alefluid_ = true;
  // for immersed we want to write eledata in every step
  bool write_eledata_every_step_ = true;

  // write standard output if no arguments are provided (default -1)
  if (step == -1 and time == -1.0) FluidField()->Output();
  // write extra output for specified step and time
  else
  {
    // print info to screen
    if (FluidField()->Discretization()->Comm().MyPID() == 0)
      std::cout << "\n   Write EXTRA FLUID Output Step=" << step << " Time=" << time << " ...   \n"
                << std::endl;

    // step number and time
    FluidField()->DiscWriter()->NewStep(step, time);

    // time step, especially necessary for adaptive dt
    FluidField()->DiscWriter()->WriteDouble("timestep", FluidField()->Dt());

    // velocity/pressure vector
    FluidField()->DiscWriter()->WriteVector("velnp", FluidField()->Velnp());
    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure =
        FluidField()->GetVelPressSplitter()->ExtractCondVector(FluidField()->Velnp());
    FluidField()->DiscWriter()->WriteVector("pressure", pressure);

    if (alefluid_) FluidField()->DiscWriter()->WriteVector("dispnp", FluidField()->Dispnp());

    // write domain decomposition for visualization (only once!)
    if ((FluidField()->Step() == upres_ or FluidField()->Step() == 0) and
        !write_eledata_every_step_)
      FluidField()->DiscWriter()->WriteElementData(true);
    else
      FluidField()->DiscWriter()->WriteElementData(true);

    return;

  }  // write extra output
}

FOUR_C_NAMESPACE_CLOSE
