/*----------------------------------------------------------------------*/
/*! \file

 \brief Fluid field adapter for poroelasticity




\level 2

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#include "4C_adapter_fld_poro.hpp"

#include "4C_discretization_condition_utils.hpp"
#include "4C_discretization_fem_general_assemblestrategy.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
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

  discretization()->GetCondition("no_penetration", nopencond_);
}

/*======================================================================*/
/* evaluate poroelasticity specific constraint*/
void ADAPTER::FluidPoro::evaluate_no_penetration_cond(Teuchos::RCP<Epetra_Vector> Cond_RHS,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> ConstraintMatrix,
    Teuchos::RCP<CORE::LINALG::SparseMatrix> struct_vel_constraint_matrix,
    Teuchos::RCP<Epetra_Vector> condVector, Teuchos::RCP<std::set<int>> condIDs,
    POROELAST::Coupltype coupltype)
{
  if (!(discretization()->Filled())) FOUR_C_THROW("fill_complete() was not called");
  if (!discretization()->HaveDofs()) FOUR_C_THROW("assign_degrees_of_freedom() was not called");

  discretization()->set_state(0, "dispnp", Dispnp());
  discretization()->set_state(0, "scaaf", Scaaf());

  Teuchos::ParameterList params;

  params.set("timescale", TimeScaling());

  if (coupltype == POROELAST::fluidfluid)
  {
    // first, find out which dofs will be constraint
    params.set<int>("action", FLD::no_penetrationIDs);
    params.set<int>("Physical Type", INPAR::FLUID::poro);
    discretization()->evaluate_condition(params, condVector, "no_penetration");

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

    CORE::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
        0,                                       // fluiddofset for column
        ConstraintMatrix,                        // fluid matrix
        Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

    discretization()->set_state(0, "condVector", condVector);

    discretization()->evaluate_condition(params, fluidstrategy, "no_penetration");
  }
  else if (coupltype == POROELAST::fluidstructure)
  {
    discretization()->set_state(0, "velnp", Velnp());
    discretization()->set_state(0, "gridv", GridVel());

    discretization()->set_state(0, "condVector", condVector);

    // set action for elements
    params.set<int>("action", FLD::no_penetration);
    params.set<POROELAST::Coupltype>("coupling", POROELAST::fluidstructure);
    params.set<int>("Physical Type", INPAR::FLUID::poro);

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of fluid_field:
    // fluiddofset = 0, structdofset = 1
    CORE::FE::AssembleStrategy couplstrategy(0,  // fluiddofset for row
        1,                                       // structdofset for column
        ConstraintMatrix,
        struct_vel_constraint_matrix,  // fluid-mechanical matrix
        Cond_RHS, Teuchos::null, Teuchos::null);

    // evaluate the fluid-mechanical system matrix on the fluid element
    discretization()->evaluate_condition(params, couplstrategy, "no_penetration");
  }
  else
    FOUR_C_THROW("unknown coupling type for no penetration BC");

  discretization()->ClearState();

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
  if (step == -1 and time == -1.0) fluid_field()->Output();
  // write extra output for specified step and time
  else
  {
    // print info to screen
    if (fluid_field()->discretization()->Comm().MyPID() == 0)
      std::cout << "\n   Write EXTRA FLUID Output Step=" << step << " Time=" << time << " ...   \n"
                << std::endl;

    // step number and time
    fluid_field()->DiscWriter()->NewStep(step, time);

    // time step, especially necessary for adaptive dt
    fluid_field()->DiscWriter()->WriteDouble("timestep", fluid_field()->Dt());

    // velocity/pressure vector
    fluid_field()->DiscWriter()->WriteVector("velnp", fluid_field()->Velnp());
    // (hydrodynamic) pressure
    Teuchos::RCP<Epetra_Vector> pressure =
        fluid_field()->GetVelPressSplitter()->ExtractCondVector(fluid_field()->Velnp());
    fluid_field()->DiscWriter()->WriteVector("pressure", pressure);

    if (alefluid_) fluid_field()->DiscWriter()->WriteVector("dispnp", fluid_field()->Dispnp());

    // write domain decomposition for visualization (only once!)
    if ((fluid_field()->Step() == upres_ or fluid_field()->Step() == 0) and
        !write_eledata_every_step_)
      fluid_field()->DiscWriter()->WriteElementData(true);
    else
      fluid_field()->DiscWriter()->WriteElementData(true);

    return;

  }  // write extra output
}

FOUR_C_NAMESPACE_CLOSE
