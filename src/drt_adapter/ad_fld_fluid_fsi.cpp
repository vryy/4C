/*----------------------------------------------------------------------*/
/*! \file

\brief Fluid field adapter for fsi. Can only be used in conjunction with FLD::FluidImplicitTimeInt

\level 2

*/
/*----------------------------------------------------------------------*/
#include "ad_fld_fluid_fsi.H"

#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_fluid/fluidimplicitintegration.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid/fluid_utils.H"
#include "../drt_fluid_ele/fluid_ele_action.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_utils_sparse_algebra_assemble.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../linalg/linalg_solver.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_RCP.hpp>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <vector>
#include <set>
/*======================================================================*/
/* constructor */
ADAPTER::FluidFSI::FluidFSI(Teuchos::RCP<Fluid> fluid, Teuchos::RCP<DRT::Discretization> dis,
    Teuchos::RCP<LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<IO::DiscretizationWriter> output, bool isale, bool dirichletcond)
    : FluidWrapper(fluid),
      dis_(dis),
      params_(params),
      output_(output),
      dirichletcond_(dirichletcond),
      interface_(Teuchos::rcp(new FLD::UTILS::MapExtractor())),
      meshmap_(Teuchos::rcp(new LINALG::MapExtractor())),
      locerrvelnp_(Teuchos::null),
      auxintegrator_(INPAR::FSI::timada_fld_none),
      numfsidbcdofs_(0),
      methodadapt_(ada_none)
{
  // make sure
  if (fluid_ == Teuchos::null) dserror("Failed to create the underlying fluid adapter");
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::Init()
{
  // call base class init
  FluidWrapper::Init();

  // cast fluid to fluidimplicit
  if (fluidimpl_ == Teuchos::null)
    fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_);

  if (fluidimpl_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::FluidImplicitTimeInt.");

  // default dofset for coupling
  int nds_master = 0;

  // set nds_master = 2 in case of HDG discretization
  // (nds = 0 used for trace values, nds = 1 used for interior values)
  if (DRT::Problem::Instance()->SpatialApproximationType() == ShapeFunctionType::shapefunction_hdg)
  {
    nds_master = 2;
  }

  // create fluid map extractor
  SetupInterface(nds_master);

  fluidimpl_->SetSurfaceSplitter(&(*interface_));

  // create map of inner velocity dof (no FSI or Dirichlet conditions)
  BuildInnerVelMap();

  if (dirichletcond_)
  {
    // mark all interface velocities as dirichlet values
    fluidimpl_->AddDirichCond(Interface()->FSICondMap());
  }

  interfaceforcen_ = Teuchos::rcp(new Epetra_Vector(*(Interface()->FSICondMap())));

  // time step size adaptivity in monolithic FSI
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const bool timeadapton =
      DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"), "TIMEADAPTON");
  if (timeadapton)
  {
    // extract the type of auxiliary integrator from the input parameter list
    auxintegrator_ = DRT::INPUT::IntegralValue<INPAR::FSI::FluidMethod>(
        fsidyn.sublist("TIMEADAPTIVITY"), "AUXINTEGRATORFLUID");

    if (auxintegrator_ != INPAR::FSI::timada_fld_none)
    {
      // determine type of adaptivity
      if (AuxMethodOrderOfAccuracy() > fluidimpl_->MethodOrderOfAccuracy())
        methodadapt_ = ada_upward;
      else if (AuxMethodOrderOfAccuracy() < fluidimpl_->MethodOrderOfAccuracy())
        methodadapt_ = ada_downward;
      else
        methodadapt_ = ada_orderequal;
    }

    //----------------------------------------------------------------------------
    // Handling of Dirichlet BCs in error estimation
    //----------------------------------------------------------------------------
    // Create intersection of fluid DOFs that hold a Dirichlet boundary condition
    // and are located at the FSI interface.
    std::vector<Teuchos::RCP<const Epetra_Map>> intersectionmaps;
    intersectionmaps.push_back(GetDBCMapExtractor()->CondMap());
    intersectionmaps.push_back(Interface()->FSICondMap());
    Teuchos::RCP<Epetra_Map> intersectionmap =
        LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

    // store number of interface DOFs subject to Dirichlet BCs on structure and fluid side of the
    // interface
    numfsidbcdofs_ = intersectionmap->NumGlobalElements();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::DofRowMap() { return DofRowMap(0); }


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
    return 2. / Dt();
  else
    return 1. / Dt();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::Update()
{
  if (DRT::Problem::Instance()->SpatialApproximationType() !=
      ShapeFunctionType::shapefunction_hdg)  // TODO als fix this!
  {
    Teuchos::RCP<Epetra_Vector> interfaceforcem = Interface()->ExtractFSICondVector(TrueResidual());

    interfaceforcen_ = fluidimpl_->ExtrapolateEndPoint(interfaceforcen_, interfaceforcem);
  }

  FluidWrapper::Update();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap, true);
  Interface()->InsertFSICondVector(ivel, relax);
  fluidimpl_->LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ADAPTER::FluidFSI::InnerVelocityRowMap() { return innervelmap_; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceForces()
{
  Teuchos::RCP<Epetra_Vector> interfaceforcem = Interface()->ExtractFSICondVector(TrueResidual());

  return fluidimpl_->ExtrapolateEndPoint(interfaceforcen_, interfaceforcem);
}

/*----------------------------------------------------------------------*
 | Return interface velocity at new time level n+1                      |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceVelnp()
{
  return Interface()->ExtractFSICondVector(Velnp());
}


/*----------------------------------------------------------------------*
 | Return interface velocity at old time level n                        |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractInterfaceVeln()
{
  return Interface()->ExtractFSICondVector(Veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::ExtractFreeSurfaceVeln()
{
  return Interface()->ExtractFSCondVector(Veln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  // apply the interface velocities
  Interface()->InsertFSICondVector(ivel, fluidimpl_->WriteAccessVelnp());

  const Teuchos::ParameterList& fsipart =
      DRT::Problem::Instance()->FSIDynamicParams().sublist("PARTITIONED SOLVER");
  if (DRT::INPUT::IntegralValue<int>(fsipart, "DIVPROJECTION"))
  {
    // project the velocity field into a divergence free subspace
    // (might enhance the linear solver, but we are still not sure.)
    ProjVelToDivZero();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyInitialMeshDisplacement(
    Teuchos::RCP<const Epetra_Vector> initfluiddisp)
{
  // cast fluid to fluidimplicit
  if (fluidimpl_ == Teuchos::null)
    fluidimpl_ = Teuchos::rcp_dynamic_cast<FLD::FluidImplicitTimeInt>(fluid_);

  if (fluidimpl_ == Teuchos::null)
    dserror("Failed to cast ADAPTER::Fluid to FLD::FluidImplicitTimeInt.");

  meshmap_->InsertCondVector(initfluiddisp, fluidimpl_->CreateDispn());
  meshmap_->InsertCondVector(initfluiddisp, fluidimpl_->CreateDispnp());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
{
  meshmap_->InsertCondVector(fluiddisp, fluidimpl_->WriteAccessDispnp());

  // new grid velocity
  fluidimpl_->UpdateGridv();
}

/*----------------------------------------------------------------------*
 | Update fluid griv velocity via FD approximation           Thon 12/14 |
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::UpdateGridv()
{
  // new grid velocity via FD approximation
  fluidimpl_->UpdateGridv();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
{
  meshmap_->InsertCondVector(gridvel, fluidimpl_->WriteAccessGridVel());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm, const int nds_master)
{
  meshmap_->Setup(
      *dis_->DofRowMap(nds_master), mm, LINALG::SplitMap(*dis_->DofRowMap(nds_master), *mm));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

#ifdef DEBUG
  // check, whether maps are the same
  if (!fcx->Map().PointSameAs(veln->Map()))
  {
    dserror("Maps do not match, but they have to.");
  }
#endif

  /*
   * Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   *
   *             / = 2 / dt   if interface time integration is second order
   * with fac = |
   *             \ = 1 / dt   if interface time integration is first order
   */
  const double timescale = TimeScaling();
  fcx->Update(-timescale * Dt(), *veln, timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
{
  // get interface velocity at t(n)
  const Teuchos::RCP<const Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

#ifdef DEBUG
  // check, whether maps are the same
  if (!fcx->Map().PointSameAs(veln->Map()))
  {
    dserror("Maps do not match, but they have to.");
  }
#endif

  /*
   * Delta d(n+1,i+1) = tau * Delta u(n+1,i+1) + dt * u(n)]
   *
   *             / = dt / 2   if interface time integration is second order
   * with tau = |
   *             \ = dt       if interface time integration is first order
   */
  const double tau = 1. / TimeScaling();
  fcx->Update(Dt(), *veln, tau);
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
  fcx->Update(-timescale * Dt(), *veln, timescale);
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
  double timescale = 1. / TimeScaling();
  fcx->Update(Dt(), *veln, timescale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::IntegrateInterfaceShape()
{
  return Interface()->ExtractFSICondVector(fluidimpl_->IntegrateInterfaceShape("FSICoupling"));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::UseBlockMatrix(bool splitmatrix)
{
  Teuchos::RCP<std::set<int>> condelements = Interface()->ConditionedElementMap(*Discretization());
  fluidimpl_->UseBlockMatrix(condelements, *Interface(), *Interface(), splitmatrix);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::Solver> ADAPTER::FluidFSI::LinearSolver()
{
  return FluidWrapper::LinearSolver();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ProjVelToDivZero()
{
  // This projection affects also the inner DOFs. Unfortunately, the matrix
  // does not look nice. Hence, the inversion of B^T*B is quite costly and
  // we are not sure yet whether it is worth the effort.

  //   get maps with Dirichlet DOFs and fsi interface DOFs
  std::vector<Teuchos::RCP<const Epetra_Map>> dbcfsimaps;
  dbcfsimaps.push_back(GetDBCMapExtractor()->CondMap());
  dbcfsimaps.push_back(Interface()->FSICondMap());

  // create a map with all DOFs that have either a Dirichlet boundary condition
  // or are located on the fsi interface
  Teuchos::RCP<Epetra_Map> dbcfsimap = LINALG::MultiMapExtractor::MergeMaps(dbcfsimaps);

  // create an element map with offset
  const int numallele = Discretization()->NumGlobalElements();
  const int mapoffset = dbcfsimap->MaxAllGID() + Discretization()->ElementRowMap()->MinAllGID() + 1;
  Teuchos::RCP<Epetra_Map> elemap =
      Teuchos::rcp(new Epetra_Map(numallele, mapoffset, Discretization()->Comm()));

  // create the combination of dbcfsimap and elemap
  std::vector<Teuchos::RCP<const Epetra_Map>> domainmaps;
  domainmaps.push_back(dbcfsimap);
  domainmaps.push_back(elemap);
  Teuchos::RCP<Epetra_Map> domainmap = LINALG::MultiMapExtractor::MergeMaps(domainmaps);

  // build the corresponding map extractor
  Teuchos::RCP<LINALG::MapExtractor> domainmapex =
      Teuchos::rcp(new LINALG::MapExtractor(*domainmap, dbcfsimap));

  const int numofrowentries = 82;
  Teuchos::RCP<LINALG::SparseMatrix> B =
      Teuchos::rcp(new LINALG::SparseMatrix(*DofRowMap(), numofrowentries, false));

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  Discretization()->ClearState();
  Discretization()->SetState("dispnp", Dispnp());

  // loop over all fluid elements
  for (int lid = 0; lid < Discretization()->NumMyColElements(); lid++)
  {
    // get pointer to current element
    DRT::Element* actele = Discretization()->lColElement(lid);

    // get element location vector and ownerships
    std::vector<int> lm;
    std::vector<int> lmowner;
    std::vector<int> lmstride;
    actele->LocationVector(*Discretization(), lm, lmowner, lmstride);

    // get dimension of element matrices and vectors
    const int eledim = (int)lm.size();

    // Reshape element matrices and vectors and initialize to zero
    elevector1.Size(eledim);

    // set action in order to calculate the integrated divergence operator via an Evaluate()-call
    Teuchos::ParameterList params;
    params.set<int>("action", FLD::calc_divop);

    // call the element specific evaluate method
    actele->Evaluate(
        params, *Discretization(), lm, elematrix1, elematrix2, elevector1, elevector2, elevector3);

    // assembly
    std::vector<int> lmcol(1);
    lmcol[0] = actele->Id() + dbcfsimap->MaxAllGID() + 1;
    B->Assemble(actele->Id(), lmstride, elevector1, lm, lmowner, lmcol);
  }  // end of loop over all fluid elements

  Discretization()->ClearState();

  // insert '1's for all DBC and interface DOFs
  for (int i = 0; i < dbcfsimap->NumMyElements(); i++)
  {
    int rowid = dbcfsimap->GID(i);
    int colid = dbcfsimap->GID(i);
    B->Assemble(1.0, rowid, colid);
  }

  B->Complete(*domainmap, *DofRowMap());

  // Compute the projection operator
  Teuchos::RCP<LINALG::SparseMatrix> BTB = LINALG::Multiply(*B, true, *B, false, true);

  Teuchos::RCP<Epetra_Vector> BTvR = Teuchos::rcp(new Epetra_Vector(*domainmap));
  B->Multiply(true, *Velnp(), *BTvR);
  Teuchos::RCP<Epetra_Vector> zeros = Teuchos::rcp(new Epetra_Vector(*dbcfsimap, true));

  domainmapex->InsertCondVector(zeros, BTvR);

  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(new Epetra_Vector(*domainmap));

  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  const int simplersolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  if (simplersolvernumber == (-1))
    dserror(
        "no simpler solver, that is used to solve this system, defined for fluid pressure problem. "
        "\nPlease set LINEAR_SOLVER in FLUID DYNAMIC to a valid number!");

  Teuchos::RCP<LINALG::Solver> solver =
      Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(simplersolvernumber),
          Discretization()->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()));

  if (solver->Params().isSublist("ML Parameters"))
  {
    solver->Params().sublist("ML Parameters").set("PDE equations", 1);
    solver->Params().sublist("ML Parameters").set("null space: dimension", 1);
    const int plength = BTB->RowMap().NumMyElements();
    Teuchos::RCP<std::vector<double>> pnewns = Teuchos::rcp(new std::vector<double>(plength, 1.0));
    solver->Params().sublist("ML Parameters").set("null space: vectors", &((*pnewns)[0]));
    solver->Params().sublist("ML Parameters").remove("nullspace", false);  // necessary?
    solver->Params()
        .sublist("Michael's secret vault")
        .set<Teuchos::RCP<std::vector<double>>>("pressure nullspace", pnewns);
  }

  solver->Solve(BTB->EpetraOperator(), x, BTvR, true, true);

  Teuchos::RCP<Epetra_Vector> vmod = Teuchos::rcp(new Epetra_Vector(Velnp()->Map(), true));
  B->Apply(*x, *vmod);
  WriteAccessVelnp()->Update(-1.0, *vmod, 1.0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::Reset(bool completeReset, int numsteps, int iter)

{
  FluidWrapper::Reset(completeReset, numsteps, iter);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::CalculateError()

{
  fluidimpl_->EvaluateErrorComparedToAnalyticalSol();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::TimeStepAuxiliar()
{
  // current state
  Teuchos::RCP<const Epetra_Vector> veln = Teuchos::rcp(new Epetra_Vector(*Veln()));
  Teuchos::RCP<const Epetra_Vector> accn = Teuchos::rcp(new Epetra_Vector(*Accn()));

  // prepare vector for solution of auxiliary time step
  locerrvelnp_ = Teuchos::rcp(new Epetra_Vector(*fluid_->DofRowMap(), true));

  // ---------------------------------------------------------------------------

  // calculate time step with auxiliary time integrator, i.e. the extrapolated solution
  switch (auxintegrator_)
  {
    case INPAR::FSI::timada_fld_none:
    {
      break;
    }
    case INPAR::FSI::timada_fld_expleuler:
    {
      ExplicitEuler(*veln, *accn, *locerrvelnp_);

      break;
    }
    case INPAR::FSI::timada_fld_adamsbashforth2:
    {
      if (Step() >= 1)  // AdamsBashforth2 only if at least second time step
      {
        // Acceleration from previous time step
        Teuchos::RCP<Epetra_Vector> accnm =
            Teuchos::rcp(new Epetra_Vector(*ExtractVelocityPart(Accnm())));

        AdamsBashforth2(*veln, *accn, *accnm, *locerrvelnp_);
      }
      else  // ExplicitEuler as starting algorithm
      {
        ExplicitEuler(*veln, *accn, *locerrvelnp_);
      }

      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integration scheme for fluid field.");
      break;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::ExplicitEuler(
    const Epetra_Vector& veln, const Epetra_Vector& accn, Epetra_Vector& velnp) const
{
  // Do a single explicit Euler step
  velnp.Update(1.0, veln, Dt(), accn, 0.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::AdamsBashforth2(const Epetra_Vector& veln, const Epetra_Vector& accn,
    const Epetra_Vector& accnm, Epetra_Vector& velnp) const
{
  // time step sizes of current and previous time step
  const double dt = Dt();
  const double dto = fluidimpl_->DtPrevious();

  // Do a single Adams-Bashforth 2 step
  velnp.Update(1.0, veln, 0.0);
  velnp.Update((2.0 * dt * dto + dt * dt) / (2 * dto), accn, -dt * dt / (2.0 * dto), accnm, 1.0);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::IndicateErrorNorms(double& err, double& errcond, double& errother,
    double& errinf, double& errinfcond, double& errinfother)
{
  // compute estimation of local discretization error
  if (methodadapt_ == ada_orderequal)
  {
    const double coeffmarch = fluidimpl_->MethodLinErrCoeffVel();
    const double coeffaux = AuxMethodLinErrCoeffVel();
    locerrvelnp_->Update(-1.0, *Velnp(), 1.0);
    locerrvelnp_->Scale(coeffmarch / (coeffaux - coeffmarch));
  }
  else
  {
    // schemes do not have the same order of accuracy
    locerrvelnp_->Update(-1.0, *Velnp(), 1.0);
  }

  // set '0' on all pressure DOFs
  Teuchos::RCP<const Epetra_Vector> zeros =
      Teuchos::rcp(new Epetra_Vector(locerrvelnp_->Map(), true));
  LINALG::ApplyDirichlettoSystem(locerrvelnp_, zeros, *(PressureRowMap()));
  // TODO: Do not misuse ApplyDirichlettoSystem()...works for this purpose here: writes zeros into
  // all pressure DoFs

  // set '0' on Dirichlet DOFs
  zeros = Teuchos::rcp(new Epetra_Vector(locerrvelnp_->Map(), true));
  LINALG::ApplyDirichlettoSystem(locerrvelnp_, zeros, *(GetDBCMapExtractor()->CondMap()));

  // extract the condition part of the full error vector (i.e. only interface velocity DOFs)
  Teuchos::RCP<Epetra_Vector> errorcond =
      Teuchos::rcp(new Epetra_Vector(*Interface()->ExtractFSICondVector(locerrvelnp_)));

  /* in case of structure split: extract the other part of the full error vector
   * (i.e. interior velocity and all pressure DOFs) */
  Teuchos::RCP<Epetra_Vector> errorother =
      Teuchos::rcp(new Epetra_Vector(*Interface()->ExtractOtherVector(locerrvelnp_)));

  // calculate L2-norms of different subsets of temporal discretization error vector
  // (neglect Dirichlet and pressure DOFs for length scaling)
  err = CalculateErrorNorm(*locerrvelnp_,
      GetDBCMapExtractor()->CondMap()->NumGlobalElements() + PressureRowMap()->NumGlobalElements());
  errcond = CalculateErrorNorm(*errorcond, numfsidbcdofs_);
  errother = CalculateErrorNorm(
      *errorother, PressureRowMap()->NumGlobalElements() +
                       (GetDBCMapExtractor()->CondMap()->NumGlobalElements() - numfsidbcdofs_));

  // calculate L-inf-norms of temporal discretization errors
  locerrvelnp_->NormInf(&errinf);
  errorcond->NormInf(&errinfcond);
  errorother->NormInf(&errinfother);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> ADAPTER::FluidFSI::CalculateWallShearStresses()
{
  // get inputs
  Teuchos::RCP<const Epetra_Vector> trueresidual = fluidimpl_->TrueResidual();
  double dt = fluidimpl_->Dt();

  // Get WSSManager
  Teuchos::RCP<FLD::UTILS::StressManager> stressmanager = fluidimpl_->StressManager();

  // Since the WSS Manager cannot be initialized in the FluidImplicitTimeInt::Init()
  // it is not so sure if the WSSManager is jet initialized. So let's be safe here..
  if (stressmanager == Teuchos::null) dserror("Call of StressManager failed!");
  if (not stressmanager->IsInit()) dserror("StressManager has not been initialized jet!");

  // Call StressManager to calculate WSS from residual
  Teuchos::RCP<Epetra_Vector> wss = stressmanager->GetWallShearStresses(trueresidual, dt);

  return wss;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double ADAPTER::FluidFSI::CalculateErrorNorm(const Epetra_Vector& vec, const int numneglect) const
{
  double norm = 1.0e+12;

  vec.Norm2(&norm);

  if (vec.GlobalLength() - numneglect > 0.0)
    norm /= sqrt((double)(vec.GlobalLength() - numneglect));
  else
    norm = 0.0;

  return norm;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int ADAPTER::FluidFSI::AuxMethodOrderOfAccuracy() const
{
  if (auxintegrator_ == INPAR::FSI::timada_fld_none)
    return 0;
  else if (auxintegrator_ == INPAR::FSI::timada_fld_expleuler)
    return 1;
  else if (auxintegrator_ == INPAR::FSI::timada_fld_adamsbashforth2)
    return 2;
  else
  {
    dserror("Unknown auxiliary time integration scheme for fluid field.");
    return 0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double ADAPTER::FluidFSI::AuxMethodLinErrCoeffVel() const
{
  if (auxintegrator_ == INPAR::FSI::timada_fld_none)
    return 0.0;
  else if (auxintegrator_ == INPAR::FSI::timada_fld_expleuler)
    return 0.5;
  else if (auxintegrator_ == INPAR::FSI::timada_fld_adamsbashforth2)
  {
    // time step sizes of current and previous time step
    const double dtc = Dt();
    const double dto = fluidimpl_->DtPrevious();

    // leading error coefficient
    return (2 * dtc + 3 * dto) / (12 * dtc);
  }
  else
  {
    dserror("Unknown auxiliary time integration scheme for fluid field.");
    return 0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double ADAPTER::FluidFSI::GetTimAdaErrOrder() const
{
  if (auxintegrator_ != INPAR::FSI::timada_fld_none)
  {
    if (methodadapt_ == ada_upward)
      return fluidimpl_->MethodOrderOfAccuracyVel();
    else
      return AuxMethodOrderOfAccuracy();
  }
  else
  {
    dserror(
        "Cannot return error order for adaptive time integration, since"
        "no auxiliary scheme has been chosen for the fluid field.");
    return 0.0;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const std::string ADAPTER::FluidFSI::GetTimAdaMethodName() const
{
  switch (auxintegrator_)
  {
    case INPAR::FSI::timada_fld_none:
    {
      return "none";
      break;
    }
    case INPAR::FSI::timada_fld_expleuler:
    {
      return "ExplicitEuler";
      break;
    }
    case INPAR::FSI::timada_fld_adamsbashforth2:
    {
      return "AdamsBashfort2";
      break;
    }
    default:
    {
      dserror("Unknown auxiliary time integration scheme for fluid field.");
      return "";
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::SetupInterface(const int nds_master)
{
  interface_->Setup(*dis_, false, false, nds_master);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::BuildInnerVelMap()
{
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;
  maps.push_back(FluidWrapper::VelocityRowMap());
  maps.push_back(Interface()->OtherMap());
  maps.push_back(GetDBCMapExtractor()->OtherMap());
  innervelmap_ = LINALG::MultiMapExtractor::IntersectMaps(maps);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ADAPTER::FluidFSI::UpdateSlaveDOF(Teuchos::RCP<Epetra_Vector>& f)
{
  fluidimpl_->UpdateSlaveDOF(f);
}
