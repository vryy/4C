#ifdef CCADISCRET

#include "fsi_monolithiclagrange.H"
#include "fsi_statustest.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicLagrange::MonolithicLagrange(Epetra_Comm& comm)
  : BlockMonolithic(comm)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  SetDefaultParameters(fsidyn,NOXParameterList());
  linearsolverstrategy_ = Teuchos::getIntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

  // right now we use matching meshes at the interface

  // set up the coupling objects
  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // structure (master) to fluid (slave)
  coupsf.SetupConditionCoupling(*StructureField().Discretization(),
                                StructureField().Interface(),
                                *FluidField().Discretization(),
                                FluidField().Interface(),
                                "FSICoupling");

  // structure (master) to ale (slave)
  coupsa.SetupConditionCoupling(*StructureField().Discretization(),
                                StructureField().Interface(),
                                *AleField().Discretization(),
                                AleField().Interface(),
                                "FSICoupling");

  // assure that both coupling objects use the same map
  // at their common master side (i.e. the structure)
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("Structure interface dof maps do not match!");
  // and that it is in fact existant
  if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // In the following we can assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.

  // fluid (master) to ale (slave)
  // note: the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();
  coupfa.SetupCoupling(*FluidField().Discretization(),
                       *AleField().Discretization(),
                       *fluidnodemap,
                       *alenodemap);

  FluidField().SetMeshMap(coupfa.MasterDofMap());

  // fluid (master) to ale (slave) only at the interface
  icoupfa_.SetupConditionCoupling(*FluidField().Discretization(),
                                  FluidField().Interface(),
                                  *AleField().Discretization(),
                                  AleField().Interface(),
                                  "FSICoupling");

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField().DofRowMap());
  vecSpaces.push_back(FluidField()    .DofRowMap());
  vecSpaces.push_back(AleField()      .Interface().OtherMap());
  vecSpaces.push_back(UTILS::ShiftMap(StructureField().Interface().CondMap(),vecSpaces));

  SetDofRowMaps(vecSpaces);

  // create block system matrix
  systemmatrix_ = Teuchos::rcp(new LagrangianBlockMatrix(Extractor(),
                                                         StructureField().LinearSolver(),
                                                         FluidField().LinearSolver(),
                                                         AleField().LinearSolver()));

  // Switch fluid to interface split block matrix
  // TODO do we need splitting up natrices? (guess not)
  /*
  FluidField().UseBlockMatrix(FluidField().Interface(),
                              FluidField().Interface(),
                              true);
  */

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  // setup Lagrangian coupling terms
  coupsf.SetupCouplingMatrices(*vecSpaces[3],
                               *StructureField().DofRowMap(),
                               *FluidField().DofRowMap());
                               //*FluidField().Interface().CondMap());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::SetupRHS");

  SetupVector(f,
              StructureField().RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  if (firstcall)
  {
    // additional rhs term for ALE equations
    // -dt Aig u(n)
    //
    //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
    //
    // And we are concerned with the u(n) part here.

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();
    if (a==Teuchos::null)
      dserror("expect ale block matrix");

    LINALG::SparseMatrix& aig = a->Matrix(0,1);

    // extract u,gamma,n (and translate it to the other maps)
    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    Teuchos::RCP<Epetra_Vector> sveln = FluidToStruct(fveln);
    Teuchos::RCP<Epetra_Vector> aveln = StructToAle(sveln);

    Teuchos::RCP<Epetra_Vector> rhs;

    // ale (#2): A,ig * u,gamma,n
    rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
    aig.Apply(*aveln,*rhs);

    rhs->Scale(-1.*Dt());
    Extractor().AddVector(*rhs,2,f);

    // lagrangian multiplier (#3): -C,fs * u,gamma,n
    Teuchos::RCP<Epetra_CrsMatrix> cfs = StructureFluidCoupling().SlaveToMasterMat();
    rhs = Teuchos::rcp(new Epetra_Vector(cfs->RowMap()));
    cfs->Apply(*FluidField().Veln(),*rhs);

    rhs->Scale( 1.*Dt());
    Extractor().AddVector(*rhs,3,f);

    // shape derivatives (#1): (F^G,gg + F^G,ig) * u,gamma,n
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().MeshMoveMatrix();
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
      LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
      fmig.Apply(*fveln,*rhs);
      Teuchos::RCP<Epetra_Vector> veln = FluidField().Interface().InsertOtherVector(rhs);

      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
      fmgg.Apply(*fveln,*rhs);
      FluidField().Interface().InsertCondVector(rhs,veln);

      veln->Scale(-1.*Dt());
      Extractor().AddVector(*veln,1,f);
    }
  }

  // NOX expects a different sign here
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::SetupSystemMatrix");

  // extract Jacobian matrices and put them into composite system

  // get all the prerequisites...

  // coupling objects
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  // (unused) const ADAPTER::Coupling& coupsa = StructureAleCoupling();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // system matrices
  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();
  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);

  // scaling factors (needed for coupling matrices)
  // note: these scaling factors are not guaranteed to be available until the
  // elements have been evaluated.
  // note: delta = delta_t * theta
  double timescale = FluidField().TimeScaling();	// delta_t
  double resscale  = FluidField().ResidualScaling();	// theta

  // mesh move matrix (needed for shape derivatives)
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().MeshMoveMatrix();

  // ... and assemble the matrices into the combined block matrix
  // note: we have to make sure that the maps of the blocks we want to insert
  // here match the maps of the combined block matrix

  // structure
  mat.Assign(0,0,View,*s);

  // fluid
  mat.Assign(1,1,View,*f);

  // ale
  mat.Assign(2,2,View,aii);

  aigtransform_(*a,             // full src
                aig,            // src
                1./timescale,   // scale
                ADAPTER::Coupling::SlaveConverter(icoupfa_), // converter
                mat.Matrix(2,1)); // dst

#if 1

  // lagrangian coupling matrices
  mat.Matrix(3,0).Add(*coupsf.MasterToMasterMat(),false,1.0,0.0);
  mat.Matrix(3,1).Add(*coupsf.SlaveToMasterMat() ,false,-1./timescale,0.0);
  //mat.Matrix(3,1).Add(*coupsf.SlaveToMasterMat() ,false,-1.,0.0);

  mat.Matrix(0,3).Add(*coupsf.MasterToMasterMatTrans(),false,1.0,0.0);
  mat.Matrix(1,3).Add(*coupsf.SlaveToMasterMatTrans(),false,-1./resscale,0.0);
  //mat.Matrix(1,3).Add(*coupsf.SlaveToMasterMatTrans() ,false,-1.,0.0);

#else

  Epetra_Vector v(*Extractor().Map(3));
  v.PutScalar(1.);
  LINALG::SparseMatrix d(v);
  d.Complete();
  mat.Matrix(3,3).Add(d,false,1.0,0.0);

#endif

  // add optional fluid linearization with respect to mesh motion block
  if (mmm != Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    // add fmgg, fmig to (1,1)-block (in addition to fluid)
    mat.Matrix(1,1).Add(fmgg,false,1./timescale,1.0);
    mat.Matrix(1,1).Add(fmig,false,1./timescale,1.0);

    // insert fmgi, fmii into (1,2)-block using columnmap transformations between fluid and ale
    fmgitransform_(*mmm,        // full src
                   fmgi,        // src
                   1.,          // scale
                   ADAPTER::Coupling::MasterConverter(coupfa), // converter
                   mat.Matrix(1,2), // dst
                   false,       // exactmatch
                   false);      // addmatrix
    fmiitransform_(*mmm,        // full src
                   fmii,        // src
                   1.,          // scale
                   ADAPTER::Coupling::MasterConverter(coupfa), // converter
                   mat.Matrix(1,2), // dst
                   false,       // exactmatch
                   true);       // addmatrix
  }

  // Done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::InitialGuess");

  SetupVector(*ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  //should we scale the system?
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const bool scaling_infnorm = (bool)getIntegralValue<int>(fsidyn,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    scolsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    acolsum_ = rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,3).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(3,2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

    Extractor().InsertVector(*sx,0,b);
    Extractor().InsertVector(*ax,2,b);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  const bool scaling_infnorm = (bool)getIntegralValue<int>(fsidyn,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = Extractor().ExtractVector(x,0);
    Teuchos::RCP<Epetra_Vector> ay = Extractor().ExtractVector(x,2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0))
      dserror("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0))
      dserror("ale scaling failed");

    Extractor().InsertVector(*sy,0,x);
    Extractor().InsertVector(*ay,2,x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b,0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b,2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0))
      dserror("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0))
      dserror("ale scaling failed");

    Extractor().InsertVector(*sx,0,b);
    Extractor().InsertVector(*ax,2,b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,3).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(3,2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::SetupVector(Epetra_Vector &f,
                                          Teuchos::RCP<const Epetra_Vector> sv,
                                          Teuchos::RCP<const Epetra_Vector> fv,
                                          Teuchos::RCP<const Epetra_Vector> av,
                                          double fluidscale)
{
  f.PutScalar(0.);

  Extractor().InsertVector(*sv,0,f);

  // TODO do we need fluidscale here?
//   if (fluidscale != 0)
//   {
//     Teuchos::RCP<Epetra_Vector> modfv = Teuchos::rcp(new Epetra_Vector(fv->Map()));
//     modfv->Update(fluidscale,*fv,0.0);

//     Extractor().InsertVector(*modfv,1,f);
//   }
//   else
  {
    Extractor().InsertVector(*fv,1,f);
  }

  Teuchos::RCP<Epetra_Vector> aov = AleField().Interface().ExtractOtherVector(av);
  Extractor().InsertVector(*aov,2,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MonolithicLagrange::CreateLinearSystem(ParameterList& nlParams,
                                            NOX::Epetra::Vector& noxSoln,
                                            Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  NOX::Epetra::Interface::Jacobian* iJac = this;
  NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP< Epetra_Operator > J = systemmatrix_;
  const Teuchos::RCP< Epetra_Operator > M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
    linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                               lsParams,
                                                               Teuchos::rcp(iJac,false),
                                                               J,
                                                               Teuchos::rcp(iPrec,false),
                                                               M,
                                                               noxSoln));
    break;
  default:
    dserror("unsupported linear block solver strategy: %d", linearsolverstrategy_);
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::MonolithicLagrange::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                          Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));


  // tests for structural displacements

  std::vector<Teuchos::RCP<const Epetra_Map> > structdisp;
  structdisp.push_back(StructureField().Interface().OtherMap());
  structdisp.push_back(Teuchos::null);
  LINALG::MultiMapExtractor structdispextract(*DofRowMap(),structdisp);

  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp =
    Teuchos::rcp(new NOX::FSI::PartialNormF("displacement",
                                            structdispextract,0,
                                            nlParams.get("Norm abs disp", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("displacement update",
                                                 Extractor(),0,
                                                 nlParams.get("Norm abs disp", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(structureDisp);

  structcombo->addStatusTest(structureDisp);
  structcombo->addStatusTest(structureDispUpdate);

  converged->addStatusTest(structcombo);

  // test for interface forces

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  Teuchos::RCP<ADAPTER::Coupling::Converter> converter =
    Teuchos::rcp(new ADAPTER::Coupling::SlaveConverter(coupsf));

  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialSumNormF> interfaceForce =
    Teuchos::rcp(new NOX::FSI::PartialSumNormF("interface",
                                               LINALG::MapExtractor(*DofRowMap(),
                                                                    StructureField().Interface().CondMap(),
                                                                    Teuchos::null),
                                               1.0,
                                               LINALG::MapExtractor(*DofRowMap(),
                                                                    FluidField().Interface().CondMap(),
                                                                    Teuchos::null),
                                               FluidField().ResidualScaling(),
                                               converter,
                                               nlParams.get("Norm abs disp", 1.0e-6),
                                               NOX::FSI::PartialNormF::Scaled));

  AddStatusTest(interfaceForce);
  interfacecombo->addStatusTest(interfaceForce);
  converged->addStatusTest(interfacecombo);

  // tests for fluid velocities

  std::vector<Teuchos::RCP<const Epetra_Map> > fluidvel;
  fluidvel.push_back(FluidField().InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);

  LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(),fluidvel);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel =
    Teuchos::rcp(new NOX::FSI::PartialNormF("velocity",
                                            fluidvelextract,0,
                                            nlParams.get("Norm abs vel", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("velocity update",
                                                 fluidvelextract,0,
                                                 nlParams.get("Norm abs vel", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(innerFluidVel);

  fluidvelcombo->addStatusTest(innerFluidVel);
  fluidvelcombo->addStatusTest(innerFluidVelUpdate);

  converged->addStatusTest(fluidvelcombo);

  // setup tests for fluid pressure

  std::vector<Teuchos::RCP<const Epetra_Map> > fluidpress;
  fluidpress.push_back(FluidField().PressureRowMap());
  fluidpress.push_back(Teuchos::null);

  LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(),fluidpress);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress =
    Teuchos::rcp(new NOX::FSI::PartialNormF("pressure",
                                            fluidpressextract,0,
                                            nlParams.get("Norm abs pres", 1.0e-6),
                                            NOX::Abstract::Vector::TwoNorm,
                                            NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate =
    Teuchos::rcp(new NOX::FSI::PartialNormUpdate("pressure update",
                                                 fluidpressextract,0,
                                                 nlParams.get("Norm abs pres", 1.0e-6),
                                                 NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(fluidPress);

  fluidpresscombo->addStatusTest(fluidPress);
  fluidpresscombo->addStatusTest(fluidPressUpdate);

  converged->addStatusTest(fluidpresscombo);

  // tests for interface forces
  // tests for lagrangian multipliers

//   Teuchos::RCP<NOX::StatusTest::Combo> lagrangecombo =
//     Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

//   Teuchos::RCP<NOX::FSI::PartialNormUpdate> lagrangeUpdate =
//     Teuchos::rcp(new NOX::FSI::PartialNormUpdate("lagrange update",
//                                                  Extractor(),3,
//                                                  nlParams.get("Norm abs lagrange", 1.0e-6),
//                                                  NOX::FSI::PartialNormUpdate::Scaled));

//   lagrangecombo->addStatusTest(lagrangeUpdate);
//   converged->addStatusTest(lagrangecombo);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicLagrange::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                  Teuchos::RCP<const Epetra_Vector>& sx,
                                                  Teuchos::RCP<const Epetra_Vector>& fx,
                                                  Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicLagrange::ExtractFieldVectors");

  sx = Extractor().ExtractVector(x,0);
  fx = Extractor().ExtractVector(x,1);

  Teuchos::RCP<Epetra_Vector> fcx = FluidField().Interface().ExtractCondVector(fx);
  FluidField().VelocityToDisplacement(fcx);
  Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);

  //ax = Extractor().ExtractVector(x,2);

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface().InsertOtherVector(aox);
  AleField().Interface().InsertCondVector(acx, a);
  ax = a;
}


#endif // ifdef ccadiscret
