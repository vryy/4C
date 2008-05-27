#ifdef CCADISCRET

#include "fsi_monolithicoverlap.H"
#include "fsi_statustest.H"
#include "fsi_nox_aitken.H"
#include "fsi_nox_newton.H"
#include "fsi_debugwriter.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#define FLUIDBLOCKMATRIX

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicOverlap::MonolithicOverlap(Epetra_Comm& comm)
  : Monolithic(comm),
    havefluidstructcolmap_(false),
    havepermfluidstructcolmap_(false),
    havealestructcolmap_(false)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  SetDefaultParameters(fsidyn,NOXParameterList());

  // right now we use matching meshes at the interface

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField().Discretization(),
                                StructureField().Interface(),
                                *FluidField().Discretization(),
                                FluidField().Interface(),
                                "FSICoupling");

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField().Discretization(),
                                StructureField().Interface(),
                                *AleField().Discretization(),
                                AleField().Interface(),
                                "FSICoupling");

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // init transfer from interface to field
  //StructureField().SetInterfaceMap(coupsf.MasterDofMap());
  //FluidField()    ->SetInterfaceMap(coupsf.SlaveDofMap());
  //AleField()      ->SetInterfaceMap(coupsa.SlaveDofMap());

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField().Discretization(),
                       *AleField().Discretization(),
                       *fluidnodemap,
                       *alenodemap);

  FluidField().SetMeshMap(coupfa.MasterDofMap());

  // create combined map

  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField().DofRowMap());
  vecSpaces.push_back(FluidField()    .Interface().OtherMap());
  vecSpaces.push_back(AleField()      .Interface().OtherMap());

  SetDofRowMaps(vecSpaces);

  // create block system matrix

  systemmatrix_ = Teuchos::rcp(new OverlappingBlockMatrix(Extractor(),
                                                          StructureField().LinearSolver(),
                                                          FluidField().LinearSolver(),
                                                          AleField().LinearSolver()));

  // Complete the empty matrix. The non-empty blocks will be inserted (in
  // Filled() state) later.
  systemmatrix_->Complete();

#ifdef FLUIDBLOCKMATRIX
  /*----------------------------------------------------------------------*/
  // Switch fluid to interface split block matrix
  FluidField().UseBlockMatrix(FluidField().Interface(),
                              FluidField().Interface());
#endif

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetDefaultParameters(const Teuchos::ParameterList& fsidyn,
                                                  Teuchos::ParameterList& list)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Preconditioner", "None");
  //nlParams.set("Norm abs F", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", fsidyn.get<int>("ITEMAX"));

  nlParams.set("Norm abs pres", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel",  fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsidyn.get<double>("CONVTOL"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");

  dirParams.set<std::string>("Method","User Defined");
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this,false);
  dirParams.set("User Defined Direction Factory",newtonfactory);

  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  //Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  // status tests are expensive, but instructive
  //solverOptions.set<std::string>("Status Test Check Type","Minimal");
  solverOptions.set<std::string>("Status Test Check Type","Complete");

  // be explicit about linear solver parameters
  lsParams.set<std::string>("Aztec Solver","GMRES");
  lsParams.set<int>("Size of Krylov Subspace",25);
  lsParams.set<std::string>("Preconditioner","User Defined");
  lsParams.set<int>("Output Frequency",AZ_all);
  lsParams.set<bool>("Output Solver Details",true);

  // adaptive tolerance settings
  lsParams.set<double>("base tolerance",fsidyn.get<double>("BASETOL"));

#if 0
  // add Aitken relaxation to Newton step
  // there is nothing to be gained...
  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
    Teuchos::rcp(new NOX::FSI::AitkenFactory());
  lineSearchParams.set("Method","User Defined");
  lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

  //lineSearchParams.sublist("Aitken").set("max step size",
  //fsidyn.get<double>("MAXOMEGA"));
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic>
FSI::MonolithicOverlap::buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
                                       Teuchos::ParameterList& params) const
{
  Teuchos::RCP<NOX::FSI::Newton> newton = Teuchos::rcp(new NOX::FSI::Newton(gd,params));
  for (unsigned i=0; i<statustests_.size(); ++i)
  {
    statustests_[i]->SetNewton(newton);
  }
  return newton;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::computeF");
  Evaluate(Teuchos::rcp(&x,false));
  SetupRHS(F);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::computeJacobian");
  Evaluate(Teuchos::rcp(&x,false));
  LINALG::BlockSparseMatrixBase& mat = Teuchos::dyn_cast<LINALG::BlockSparseMatrixBase>(Jac);
  SetupSystemMatrix(mat);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computePreconditioner(const Epetra_Vector &x,
                                                   Epetra_Operator &M,
                                                   Teuchos::ParameterList *precParams)
{
  // Create preconditioner operator.
  // The blocks are already there. And this is the perfect place to initialize
  // the block preconditioners.
  systemmatrix_->SetupBlockPrecond();
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupRHS(Epetra_Vector& f)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::SetupRHS");

  SetupVector(f,
              StructureField().RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::SetupSystemMatrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupsa = StructureAleCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();

  //LINALG::PrintSparsityToPostscript(*s);

  // split fluid matrix

#ifdef FLUIDBLOCKMATRIX
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();
#else
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();
#endif

  // map between fluid interface and structure column map

  if (not havefluidstructcolmap_)
  {
#ifdef FLUIDBLOCKMATRIX
    DRT::Exporter ex(blockf->FullRowMap(),blockf->FullColMap(),blockf->Comm());
#else
    DRT::Exporter ex(f->RowMap(),f->ColMap(),f->Comm());
#endif
    coupsf.FillSlaveToMasterMap(fluidstructcolmap_);
    ex.Export(fluidstructcolmap_);
    havefluidstructcolmap_ = true;
  }

#ifndef FLUIDBLOCKMATRIX
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf =
    f->Split<LINALG::DefaultBlockMatrixStrategy>(FluidField().Interface(),
                                                 FluidField().Interface());
  blockf->Complete();
#endif

  LINALG::SparseMatrix& fii = blockf->Matrix(0,0);
  LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
  LINALG::SparseMatrix& fgi = blockf->Matrix(1,0);
  LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);

  /*----------------------------------------------------------------------*/

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

  if (a==Teuchos::null)
    dserror("expect ale block matrix");

  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);

  // map between ale interface and structure column map

  if (not havealestructcolmap_)
  {
    DRT::Exporter ex(a->FullRowMap(),a->FullColMap(),a->Comm());
    coupsa.FillSlaveToMasterMap(alestructcolmap_);
    ex.Export(alestructcolmap_);
    havealestructcolmap_ = true;
  }

  double scale     = FluidField().ResidualScaling();
  double timescale = FluidField().TimeScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  AddFluidInterface(scale*timescale,fgg,*s);

  mat.Assign(0,0,View,*s);
  mat.Assign(0,1,View,*ConvertFgiRowmap(fgi,scale,s->RowMap()));

  mat.Assign(1,0,View,*ConvertFigColmap(fig,fluidstructcolmap_,s->DomainMap(),timescale));
  mat.Assign(1,1,View,fii);

  mat.Assign(2,0,View,*ConvertFigColmap(aig,alestructcolmap_,StructureField().DomainMap(),1.));
  mat.Assign(2,2,View,aii);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::InitialGuess");

  SetupVector(*ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupVector(Epetra_Vector &f,
                                         Teuchos::RCP<const Epetra_Vector> sv,
                                         Teuchos::RCP<const Epetra_Vector> fv,
                                         Teuchos::RCP<const Epetra_Vector> av,
                                         double fluidscale)
{

  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> fov = FluidField()    .Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      .Interface().ExtractOtherVector(av);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = FluidField().Interface().ExtractCondVector(fv);
    Teuchos::RCP<Epetra_Vector> modsv = StructureField().Interface().InsertCondVector(FluidToStruct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    Extractor().InsertVector(*modsv,0,f);
  }
  else
  {
    Extractor().InsertVector(*sv,0,f);
  }

  Extractor().InsertVector(*fov,1,f);
  Extractor().InsertVector(*aov,2,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MonolithicOverlap::CreateLinearSystem(ParameterList& nlParams,
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

  linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams,
                                                             lsParams,
                                                             Teuchos::rcp(iJac,false),
                                                             J,
                                                             Teuchos::rcp(iPrec,false),
                                                             M,
                                                             noxSoln));

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::MonolithicOverlap::CreateStatusTest(Teuchos::ParameterList& nlParams,
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

  // setup tests for structural displacements

  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<PartialNormF> structureDisp =
    Teuchos::rcp(new PartialNormF("displacement",
                                  Extractor(),0,
                                  nlParams.get("Norm abs disp", 1.0e-6),
                                  PartialNormF::Scaled));
  Teuchos::RCP<PartialNormUpdate> structureDispUpdate =
    Teuchos::rcp(new PartialNormUpdate("displacement update",
                                       Extractor(),0,
                                       nlParams.get("Norm abs disp", 1.0e-6),
                                       PartialNormUpdate::Scaled));

  statustests_.push_back(structureDisp);
  structcombo->addStatusTest(structureDisp);
  //structcombo->addStatusTest(structureDispUpdate);

  converged->addStatusTest(structcombo);

  // setup tests for fluid velocities

  std::vector<Teuchos::RCP<const Epetra_Map> > fluidvel;
  fluidvel.push_back(FluidField().InnerVelocityRowMap());
  fluidvel.push_back(LINALG::SplitMap(*DofRowMap(),
                                      *FluidField().InnerVelocityRowMap()));
  LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(),fluidvel);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<PartialNormF> innerFluidVel =
    Teuchos::rcp(new PartialNormF("velocity",
                                  fluidvelextract,0,
                                  nlParams.get("Norm abs vel", 1.0e-6),
                                  PartialNormF::Scaled));
  Teuchos::RCP<PartialNormUpdate> innerFluidVelUpdate =
    Teuchos::rcp(new PartialNormUpdate("velocity update",
                                       fluidvelextract,0,
                                       nlParams.get("Norm abs vel", 1.0e-6),
                                       PartialNormUpdate::Scaled));

  statustests_.push_back(innerFluidVel);
  fluidvelcombo->addStatusTest(innerFluidVel);
  //fluidvelcombo->addStatusTest(innerFluidVelUpdate);

  converged->addStatusTest(fluidvelcombo);

  // setup tests for fluid pressure

  std::vector<Teuchos::RCP<const Epetra_Map> > fluidpress;
  fluidpress.push_back(FluidField().InnerVelocityRowMap());
  fluidpress.push_back(LINALG::SplitMap(*DofRowMap(),
                                        *FluidField().InnerVelocityRowMap()));
  LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(),fluidpress);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo =
    Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<PartialNormF> fluidPress =
    Teuchos::rcp(new PartialNormF("pressure",
                                  fluidpressextract,0,
                                  nlParams.get("Norm abs pres", 1.0e-6),
                                  PartialNormF::Scaled));
  Teuchos::RCP<PartialNormUpdate> fluidPressUpdate =
    Teuchos::rcp(new PartialNormUpdate("pressure update",
                                       fluidpressextract,0,
                                       nlParams.get("Norm abs pres", 1.0e-6),
                                       PartialNormUpdate::Scaled));

  statustests_.push_back(fluidPress);
  fluidpresscombo->addStatusTest(fluidPress);
  //fluidpresscombo->addStatusTest(fluidPressUpdate);

  converged->addStatusTest(fluidpresscombo);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                 Teuchos::RCP<const Epetra_Vector>& sx,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::ExtractFieldVectors");

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x,0);
  Teuchos::RCP<const Epetra_Vector> scx = StructureField().Interface().ExtractCondVector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);

  FluidField().ConvertInterfaceUnknown(fcx);

  Teuchos::RCP<Epetra_Vector> f = FluidField().Interface().InsertOtherVector(fox);
  FluidField().Interface().InsertCondVector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface().InsertOtherVector(aox);
  AleField().Interface().InsertCondVector(acx, a);
  ax = a;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::AddFluidInterface(double scale,
                                               const LINALG::SparseMatrix& fgg,
                                               LINALG::SparseMatrix& s)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::AddFluidInterface");

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> pfgg = coupsf.SlaveToPermSlave(fgg);
  Teuchos::RCP<Epetra_CrsMatrix> perm_fgg = pfgg->EpetraMatrix();

  const Epetra_Map& perm_rowmap = perm_fgg->RowMap();
  const Epetra_Map& perm_colmap = perm_fgg->ColMap();

  if (not havepermfluidstructcolmap_)
  {
    // fgg->RowMap()
    DRT::Exporter ex(*coupsf.SlaveDofMap(),perm_colmap,perm_fgg->Comm());
    coupsf.FillSlaveToMasterMap(permfluidstructcolmap_);
    ex.Export(permfluidstructcolmap_);
    havepermfluidstructcolmap_ = true;
  }

#define INTERFACE_MATRIX_UNCOMPLETE

#ifdef INTERFACE_MATRIX_UNCOMPLETE
  // uncomplete because the fluid interface can have more connections than the
  // structural one. (Tet elements in fluid can cause this.) We should do
  // this just once...
  s.UnComplete();
#endif

  int rows = perm_fgg->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = perm_fgg->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    int gid = perm_rowmap.GID(i);
    std::map<int,int>::iterator iter = permfluidstructcolmap_.find(gid);
    if (iter==permfluidstructcolmap_.end())
      dserror("gid %d not found", gid);
    int structgid = iter->second;

    for (int j=0; j<NumEntries; ++j)
    {
      int perm_gid = perm_colmap.GID(Indices[j]);
      iter = permfluidstructcolmap_.find(perm_gid);
      if (iter==permfluidstructcolmap_.end())
        dserror("gid %d not found", perm_gid);
      int index = iter->second;

      // There might be zeros on Dirichlet lines that are not included in the
      // structure matrix graph. Ignore them.
      if (Values[j]!=0)
      {
        double value = Values[j]*scale;
        s.Assemble(value, structgid, index);
      }
    }
  }

#ifdef INTERFACE_MATRIX_UNCOMPLETE
  s.Complete();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
FSI::MonolithicOverlap::ConvertFigColmap(const LINALG::SparseMatrix& fig,
                                         const std::map<int,int>& convcolmap,
                                         const Epetra_Map& domainmap,
                                         double scale)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::ConvertFigColmap");

  const Epetra_Map& rowmap = fig.RowMap();
  const Epetra_Map& colmap = fig.ColMap();

  Teuchos::RCP<Epetra_CrsMatrix> convfig = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,
                                                                             fig.MaxNumEntries()));

  Teuchos::RCP<Epetra_CrsMatrix> efig = fig.EpetraMatrix();

  int rows = efig->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = efig->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    std::vector<int> structindices(NumEntries);
    for (int j=0; j<NumEntries; ++j)
    {
      int gid = colmap.GID(Indices[j]);
      std::map<int,int>::const_iterator iter = convcolmap.find(gid);
      if (iter==convcolmap.end())
        dserror("gid %d not found in map for lid %d at %d", gid, Indices[j], j);
      structindices[j] = iter->second;
    }

    err = convfig->InsertGlobalValues(rowmap.GID(i), NumEntries, Values, &structindices[0]);
    if (err)
      dserror("InsertGlobalValues error: %d", err);
  }

  convfig->FillComplete(domainmap,rowmap);
  convfig->Scale(scale);

  return Teuchos::rcp(new LINALG::SparseMatrix(convfig,fig.ExplicitDirichlet(),fig.SaveGraph()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
FSI::MonolithicOverlap::ConvertFgiRowmap(const LINALG::SparseMatrix& fgi,
                                         double scale,
                                         const Epetra_Map& structrowmap)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicOverlap::ConvertFgiRowmap");

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  // redistribute fluid matrix to match distribution of structure matrix

  Teuchos::RCP<LINALG::SparseMatrix> pfgi = coupsf.SlaveToPermSlave(fgi);
  Teuchos::RCP<Epetra_CrsMatrix> efgi = pfgi->EpetraMatrix();

  // create new matrix with structure row map and fill it

  Teuchos::RCP<Epetra_CrsMatrix> sfgi = Teuchos::rcp(new Epetra_CrsMatrix(Copy,structrowmap,pfgi->MaxNumEntries()));

  const Epetra_Map& colmap = efgi->ColMap();

  int rows = efgi->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = efgi->ExtractMyRowView(i, NumEntries, Values, Indices);
    if (err!=0)
      dserror("ExtractMyRowView error: %d", err);

    // pull indices back to global
    std::vector<int> structindices(NumEntries);
    for (int j=0; j<NumEntries; ++j)
    {
      structindices[j] = colmap.GID(Indices[j]);
    }

    // put row into matrix with structure row map
    err = sfgi->InsertGlobalValues(coupsf.MasterDofMap()->GID(i), NumEntries, Values, &structindices[0]);
    if (err)
      dserror("InsertGlobalValues error: %d", err);
  }

  sfgi->FillComplete(efgi->DomainMap(),structrowmap);
  sfgi->Scale(scale);

  return Teuchos::rcp(new LINALG::SparseMatrix(sfgi,fgi.ExplicitDirichlet(),fgi.SaveGraph()));
}


#endif
