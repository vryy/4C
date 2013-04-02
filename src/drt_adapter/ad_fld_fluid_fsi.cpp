/*----------------------------------------------------------------------*/
/*!
\file ad_fld_fluid_fsi.cpp

\brief Fluid field adapter for fsi

Can only be used in conjunction with fluidimplicittimeint!

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fluid_fsi.H"

#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
/*======================================================================*/
/* constructor */
ADAPTER::FluidFSI::FluidFSI(Teuchos::RCP<Fluid> fluid,
    Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output,
    bool isale,
    bool dirichletcond)
: FluidWrapper(fluid),
  dis_(dis),
  params_(params),
  output_(output),
  interface_(Teuchos::rcp(new FLD::UTILS::MapExtractor())),
  meshmap_(Teuchos::rcp(new LINALG::MapExtractor()))
{
  // make sure
  if (fluid_ == Teuchos::null)
    dserror("Failed to create the underlying fluid adapter");

  // cast fluid to fluidimplicit
  fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_);
  if (fluidimpl_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::FluidImplicitTimeInt.");

  interface_->Setup(*dis);
  fluidimpl_->SetSurfaceSplitter(&(*interface_));

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint
  const Teuchos::RCP<const LINALG::MapExtractor> dbcmaps = fluidimpl_->DirichMaps();
  std::vector<Teuchos::RCP<const Epetra_Map> > maps;
  maps.push_back(VelocityRowMap());
  maps.push_back(interface_->OtherMap());
  maps.push_back(dbcmaps->OtherMap());
  innervelmap_ = LINALG::MultiMapExtractor::IntersectMaps(maps);

  if (dirichletcond)
  {
    // mark all interface velocities as dirichlet values
    fluidimpl_->AddDirichCond(interface_->FSICondMap());
  }

  interfaceforcen_ = Teuchos::rcp(new Epetra_Vector(*(interface_->FSICondMap())));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::DofRowMap()
{
  return DofRowMap(0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::DofRowMap(unsigned nds)
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap(nds);
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double ADAPTER::FluidFSI::TimeScaling() const
{
  if (params_->get<bool>("interface second order"))
    return 2./fluidimpl_->Dt();
  else
    return 1./fluidimpl_->Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::Update()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(fluidimpl_->TrueResidual());

  interfaceforcen_ = fluidimpl_->ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);

  fluidimpl_->TimeUpdate();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_->InsertFSICondVector(ivel,relax);
  fluidimpl_->LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(fluidimpl_->TrueResidual());

  return fluidimpl_->ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);
}

/*----------------------------------------------------------------------*
 | Return interface velocity at new time level n+1                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceVelnp()
{
  return interface_->ExtractFSICondVector(fluidimpl_->Velnp());
}


/*----------------------------------------------------------------------*
 | Return interface velocity at old time level n                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceVeln()
{
  return interface_->ExtractFSICondVector(fluidimpl_->Veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractFreeSurfaceVeln()
{
  return Interface()->ExtractFSCondVector(fluidimpl_->Veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  // apply the interface velocities
  interface_->InsertFSICondVector(ivel,fluidimpl_->ViewOfVelnp());

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  if (DRT::INPUT::IntegralValue<int>(fsidyn,"DIVPROJECTION"))
  {
    // project the velocity field into a divergence free subspace
    // (might enhance the linear solver, but we are still not sure.)
    ProjVelToDivZero();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_->InsertCondVector(fluiddisp,fluidimpl_->ViewOfDispnp());

  // new grid velocity
  fluidimpl_->UpdateGridv();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  meshmap_->InsertCondVector(gridvel,fluidimpl_->ViewOfGridVel());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_->Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().PointSameAs(veln->Map()))  { dserror("Maps do not match, but they have to."); }
#endif

  /*
   * Delta u(n+1,i+1) = fac * Delta d(n+1,i+1) - dt * u(n)
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double timescale = TimeScaling();
  fcx->Update(-timescale*fluidimpl_->Dt(),*veln,timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::DisplacementToVelocity(
    Teuchos::RCP<Epetra_Vector> fcx,
    Teuchos::RCP<Epetra_Vector> ddgpred,
    Teuchos::RCP<Epetra_Vector> dugpred
)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().PointSameAs(veln->Map()))    { dserror("Maps do not match, but they have to."); }
  if (! fcx->Map().PointSameAs(ddgpred->Map())) { dserror("Maps do not match, but they have to."); }
  if (! fcx->Map().PointSameAs(dugpred->Map())) { dserror("Maps do not match, but they have to."); }
#endif


  /*
   * Delta u(n+1,i+1) = fac * [ Delta d(n+1,i+1) - dt * u_fluid(n) + Delta d_(predicted) ]
   *
   *                  - Delta u_fluid(predicted)
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double ts = TimeScaling();
  fcx->Update(-Dt()*ts,*veln,ts,*ddgpred,ts);
  fcx->Update(-1.0,*dugpred,1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().PointSameAs(veln->Map()))  { dserror("Maps do not match, but they have to."); }
#endif

  /*
   * Delta d(n+1,i+1) = tau * Delta u(n+1,i+1) + dt * u(n)]
   *
   *             / = dt / 2   if interface time integration is second order
   * with tau = |
   *             \ = dt       if interface time integration is first order
   */
  const double tau = 1./TimeScaling();
  fcx->Update(Dt(), *veln, tau);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::VelocityToDisplacement(
    Teuchos::RCP<Epetra_Vector> fcx,
    Teuchos::RCP<Epetra_Vector> ddgpred,
    Teuchos::RCP<Epetra_Vector> dugpred
)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

#ifdef DEBUG
  // check, whether maps are the same
  if (! fcx->Map().PointSameAs(veln->Map()))    { dserror("Maps do not match, but they have to."); }
  if (! fcx->Map().PointSameAs(ddgpred->Map())) { dserror("Maps do not match, but they have to."); }
  if (! fcx->Map().PointSameAs(dugpred->Map())) { dserror("Maps do not match, but they have to."); }
#endif


  /*
   * Delta d(n+1,i+1) = tau * [ Delta u(n+1,i+1) + Delta u(predicted)]
   *
   *                  + dt * u(n) - Delta d_structure(predicted)
   *
   *             / = dt / 2   if interface time integration is second order
   * with tau = |
   *             \ = dt       if interface time integration is first order
   */
  const double tau = 1./TimeScaling();
  fcx->Update(Dt(), *veln, tau, *dugpred, tau);
  fcx->Update(-1.0, *ddgpred, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::FreeSurfDisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface()->ExtractFSCondVector(Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = TimeScaling();
  fcx->Update(-timescale*Dt(),*veln,timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::FreeSurfVelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<Epetra_Vector> veln = Interface()->ExtractFSCondVector(Veln());

  // We convert Delta u(n+1,i+1) to Delta d(n+1,i+1) here.
  //
  // Delta d(n+1,i+1) = ( theta Delta u(n+1,i+1) + u(n) ) * dt
  //
  double timescale = 1./TimeScaling();
  fcx->Update(Dt(),*veln,timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::IntegrateInterfaceShape()
{
  return interface_->ExtractFSICondVector(fluidimpl_->IntegrateInterfaceShape("FSICoupling"));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::UseBlockMatrix(bool splitmatrix)
{
  Teuchos::RCP<std::set<int> > condelements = Interface()->ConditionedElementMap(*Discretization());
  fluidimpl_->UseBlockMatrix(condelements,*Interface(),*Interface(),splitmatrix);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> ADAPTER::FluidFSI::LinearSolver()
{
  return fluidimpl_->LinearSolver();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ProjVelToDivZero()
{
  // This projection affects also the inner DOFs. Unfortunately, the matrix
  // does not look nice. Hence, the inversion of B^T*B is quite costly and
  // we are not sure yet whether it is worth the effort.

  //   get maps with Dirichlet DOFs and fsi interface DOFs
  std::vector<Teuchos::RCP<const Epetra_Map> > dbcfsimaps;
  dbcfsimaps.push_back(fluid_->GetDBCMapExtractor()->CondMap());
  dbcfsimaps.push_back(interface_->FSICondMap());

  // create a map with all DOFs that have either a Dirichlet boundary condition
  // or are located on the fsi interface
  Teuchos::RCP<Epetra_Map> dbcfsimap = LINALG::MultiMapExtractor::MergeMaps(dbcfsimaps);

  // create an element map with offset
  const int numallele = fluidimpl_->Discretization()->NumGlobalElements();
  const int mapoffset = dbcfsimap->MaxAllGID()+fluidimpl_->Discretization()->ElementRowMap()->MinAllGID() + 1;
  Teuchos::RCP<Epetra_Map> elemap = Teuchos::rcp(new Epetra_Map(numallele,mapoffset,fluidimpl_->Discretization()->Comm()));

  // create the combination of dbcfsimap and elemap
  std::vector<Teuchos::RCP<const Epetra_Map> > domainmaps;
  domainmaps.push_back(dbcfsimap);
  domainmaps.push_back(elemap);
  Teuchos::RCP<Epetra_Map> domainmap = LINALG::MultiMapExtractor::MergeMaps(domainmaps);

  // build the corresponding map extractor
  Teuchos::RCP<LINALG::MapExtractor> domainmapex = Teuchos::rcp(new LINALG::MapExtractor(*domainmap, dbcfsimap));

  const int numofrowentries = 82;
  Teuchos::RCP<LINALG::SparseMatrix> B = Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(),numofrowentries,false));

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  fluid_->Discretization()->ClearState();
  fluid_->Discretization()->SetState("dispnp",fluid_->Dispnp());

  // loop over all fluid elements
  for (int lid = 0; lid < fluid_->Discretization()->NumMyColElements(); lid++)
  {
    // get pointer to current element
    DRT::Element * actele = fluid_->Discretization()->lColElement(lid);

    // get element location vector and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    actele->LocationVector(*fluid_->Discretization(),lm,lmowner,lmstride);

    // get dimension of element matrices and vectors
    const int eledim = (int)lm.size();

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(eledim);

    // set action in order to calculate the integrated divergence operator via an Evaluate()-call
    Teuchos::ParameterList params;
    params.set<int>("action",FLD::calc_divop);

    // call the element specific evaluate method
    actele->Evaluate(params,*fluid_->Discretization(),lm,elematrix1,elematrix2,
            elevector1,elevector2,elevector3);

    // assembly
    std::vector<int> lmcol(1);
    lmcol[0] = actele->Id() + dbcfsimap->MaxAllGID() + 1;
    B->Assemble(actele->Id(),lmstride,elevector1,lm,lmowner,lmcol);
  } // end of loop over all fluid elements

  fluid_->Discretization()->ClearState();

  // insert '1's for all DBC and interface DOFs
  for (int i = 0; i < dbcfsimap->NumMyElements(); i++)
  {
    int rowid = dbcfsimap->GID(i);
    int colid = dbcfsimap->GID(i);
    B->Assemble(1.0,rowid,colid);
  }

  B->Complete(*domainmap,*DofRowMap());

  // Compute the projection operator
  Teuchos::RCP<LINALG::SparseMatrix> BTB = LINALG::Multiply(*B,true,*B,false,true);

  Teuchos::RCP<Epetra_Vector> BTvR = Teuchos::rcp(new Epetra_Vector(*domainmap));
  B->Multiply(true,*fluidimpl_->ViewOfVelnp(),*BTvR);
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*dbcfsimap, true));

  domainmapex->InsertCondVector(zeros,BTvR);

  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*domainmap));

  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  const int simplersolvernumber = fdyn.get<int>("SIMPLER_SOLVER");
  if (simplersolvernumber == (-1))
    dserror("no simpler solver, that is used to solve this system, defined for fluid pressure problem. \nPlease set SIMPLER_SOLVER in FLUID DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(simplersolvernumber),
                           fluidimpl_->Discretization()->Comm(),
                           DRT::Problem::Instance()->ErrorFile()->Handle()));

  if(solver->Params().isSublist("ML Parameters"))
  {
    solver->Params().sublist("ML Parameters").set("PDE equations",1);
    solver->Params().sublist("ML Parameters").set("null space: dimension",1);
    const int plength = BTB->RowMap().NumMyElements();
    Teuchos::RCP<std::vector<double> > pnewns = Teuchos::rcp(new std::vector<double>(plength,1.0));
    solver->Params().sublist("ML Parameters").set("null space: vectors",&((*pnewns)[0]));
    solver->Params().sublist("ML Parameters").remove("nullspace",false); // necessary?
    solver->Params().sublist("Michael's secret vault").set<RCP<std::vector<double> > >("pressure nullspace",pnewns);
  }

  solver->Solve(BTB->EpetraOperator(),x,BTvR,true,true);

  Teuchos::RCP<Epetra_Vector> vmod = Teuchos::rcp(new Epetra_Vector(fluidimpl_->ViewOfVelnp()->Map(),true));
  B->Apply(*x,*vmod);
  fluidimpl_->ViewOfVelnp()->Update(-1.0, *vmod, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  fluidimpl_->RemoveDirichCond(maptoremove);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::Reset(bool completeReset, bool newFiles, int iter)

{
  fluidimpl_->Reset(completeReset, newFiles, iter);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::CalculateError()

{
  fluidimpl_->EvaluateErrorComparedToAnalyticalSol();
  return;
}

