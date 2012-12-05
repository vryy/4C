/*!----------------------------------------------------------------------*
 \file ad_fld_poro.cpp

 \brief Fluid field adapter for poroelasticity

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*/
#include "ad_fld_poro.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_io/io.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../drt_lib/drt_assemblestrategy.H"

/*======================================================================*/
/* constructor */
ADAPTER::FluidPoro::FluidPoro(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output,
    bool isale,
    bool dirichletcond)
: FluidFSI(fluid,dis,solver,params,output,isale,dirichletcond)
{
  // make sure

  if (fluid_ == Teuchos::null)
    dserror("Failed to create the underlying fluid adapter");

  Discretization()->GetCondition("NoPenetration", nopencond_);
}

/*======================================================================*/
/* evaluate poroelasticity specific constraint*/
void ADAPTER::FluidPoro::EvaluateNoPenetrationCond(Teuchos::RCP<Epetra_Vector> Cond_RHS,
                                              Teuchos::RCP<LINALG::SparseMatrix> ConstraintMatrix,
                                              Teuchos::RCP<LINALG::SparseMatrix> StructVelConstraintMatrix,
                                              Teuchos::RCP< std::set<int> > condIDs,
                                              int coupltype)
{
  if (!(Discretization()->Filled())) dserror("FillComplete() was not called");
  if (!Discretization()->HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  Discretization()->ClearState();
  Discretization()->SetState("dispn", Dispn());
  Discretization()->SetState("dispnp", Dispnp());
  Discretization()->SetState(0,"velnp",Velnp());
  Discretization()->SetState(0,"gridv",GridVel());

  ConstraintMatrix->Zero();

  Teuchos::ParameterList params;

  params.set("timescale",TimeScaling());

  condIDs->clear();

  // set action for elements
  params.set<int>("action",FLD::no_penetration);
  if(coupltype==0)
  {
    params.set("coupling","fluid fluid");

    Teuchos::RCP<Epetra_Vector> condset = Teuchos::rcp(new Epetra_Vector(ConstraintMatrix->RowMap(),true) );

    DRT::AssembleStrategy fluidstrategy(
        0,              // fluiddofset for row
        0,              // fluiddofset for column
        ConstraintMatrix,           // fluid-mechanical matrix
        Teuchos::null,
        condset,
        Teuchos::null,
        Teuchos::null
    );

    Discretization()->EvaluateCondition( params, fluidstrategy,"NoPenetration" );

    //write global IDs of dofs on which the no penetration condition is applied (can vary in time)
    const int length = condset->MyLength();
    const Epetra_BlockMap& map = condset->Map();
    for(int i=0; i<length ; i++)
    {
      if((*condset)[i] != 0.0)
        condIDs->insert(map.GID(i));
    }
  }
  else if(coupltype==1)
  {
    params.set("coupling","fluid structure");
    StructVelConstraintMatrix->Zero();

    // build specific assemble strategy for the fluid-mechanical system matrix
    // from the point of view of FluidField:
    // fluiddofset = 0, structdofset = 1
    DRT::AssembleStrategy couplstrategy(
        0,              // fluiddofset for row
        1,              // structdofset for column
        ConstraintMatrix,
        StructVelConstraintMatrix,           // fluid-mechanical matrix
        Cond_RHS,
        Teuchos::null,
        Teuchos::null
    );

    // evaluate the fluid-mechancial system matrix on the fluid element
    Discretization()->EvaluateCondition( params, couplstrategy,"NoPenetration" );
  }
  else
    dserror("unknown coupling type for no penetration BC");

  Discretization()->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidPoro::EvaluateContiBoundaryCond(
                                              Teuchos::RCP<LINALG::SparseMatrix> FluidFluidMatrix,
                                              Teuchos::RCP<LINALG::SparseMatrix> FluidStructMatrix)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidPoro::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  fluidimpl_->AddDirichCond(maptoadd);
}

/*----------------------------------------------------------------------*/
