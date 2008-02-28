#ifdef CCADISCRET

#include "fsi_monolithicoverlap.H"
#include "fsi_statustest.H"
#include "fsi_overlapprec.H"

#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicOverlap::MonolithicOverlap(Epetra_Comm& comm)
  : Monolithic(comm),
    havefluidstructcolmap_(false),
    havepermfluidstructcolmap_(false)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  SetDefaultParameters(fsidyn,NOXParameterList());

  // right now we use matching meshes at the interface

  FSI::Coupling& coupsf = StructureFluidCoupling();
  FSI::Coupling& coupsa = StructureAleCoupling();
  FSI::Coupling& coupfa = FluidAleCoupling();

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

  /*----------------------------------------------------------------------*/
  sysmattimer_         = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::SetupSysMat");
  igtimer_             = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::InitialGuess");
  rhstimer_            = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::SetupRHS");
  exctracttimer_       = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::ExtractFieldVectors");
  fluidinterfacetimer_ = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::AddFluidInterface");
  figcolmaptimer_      = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::ConvertFigColmap");
  fgirowmaptimer_      = Teuchos::TimeMonitor::getNewTimer("FSI::MonolithicOverlap::ConvertFgiRowmap");

  /*----------------------------------------------------------------------*/
  // Assume linear ALE. Prepare ALE system matrix once and for all.

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();
  const LINALG::SparseMatrix* aii = AleField().InteriorMatrixBlock();
  const LINALG::SparseMatrix* aig = AleField().InterfaceMatrixBlock();

  if (a==Teuchos::null or aii==NULL or aig==NULL)
    dserror("expect ale block matrix");

  // map between ale interface and structure column map

  DRT::Exporter ex(a->FullRowMap(),a->FullColMap(),a->Comm());
  std::map<int,int> alestructcolmap;
  coupsa.FillSlaveToMasterMap(alestructcolmap);
  ex.Export(alestructcolmap);

  aig_ = ConvertFigColmap(*aig,alestructcolmap,StructureField().DomainMap(),1.);
  aii_ = Teuchos::rcp(new LINALG::SparseMatrix(*aii,View));
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
  nlParams.set("Norm abs F", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", fsidyn.get<int>("ITEMAX"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");

  dirParams.set<std::string>("Method","Newton");

  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  lsParams.set<std::string>("Convergence Test","r0");
  lsParams.set<double>("Tolerance",1e-6);
  lsParams.set<std::string>("Preconditioner","User Defined");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  Evaluate(Teuchos::rcp(&x,false));
  SetupRHS(F);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
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
  // We don't need to do anything special here. The blocks are already there.
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupRHS(Epetra_Vector& f) const
{
  Teuchos::TimeMonitor monitor(*rhstimer_);

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
  Teuchos::TimeMonitor monitor(*sysmattimer_);

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const FSI::Coupling& coupsf = StructureFluidCoupling();
  //const FSI::Coupling& coupsa = StructureAleCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> s = StructureField().SystemMatrix();

  //LINALG::PrintSparsityToPostscript(*s);

  // split fluid matrix

  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();

  // map between fluid interface and structure column map

  if (not havefluidstructcolmap_)
  {
    DRT::Exporter ex(f->RowMap(),f->ColMap(),f->Comm());
    coupsf.FillSlaveToMasterMap(fluidstructcolmap_);
    ex.Export(fluidstructcolmap_);
    havefluidstructcolmap_ = true;
  }

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf =
    f->Split<LINALG::DefaultBlockMatrixStrategy>(FluidField().Interface(),
                                                 FluidField().Interface());
  blockf->Complete();

  LINALG::SparseMatrix& fii = blockf->Matrix(0,0);
  LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
  LINALG::SparseMatrix& fgi = blockf->Matrix(1,0);
  LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);

  double scale = FluidField().ResidualScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  AddFluidInterface(scale/Dt(),fgg,*s);

  mat.Assign(0,0,View,*s);
  mat.Assign(0,1,View,*ConvertFgiRowmap(fgi,scale,s->RowMap()));

  mat.Assign(1,0,View,*ConvertFigColmap(fig,fluidstructcolmap_,s->DomainMap(),1./Dt()));
  mat.Assign(1,1,View,fii);

  mat.Assign(2,0,View,*aig_);
  mat.Assign(2,2,View,*aii_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  Teuchos::TimeMonitor monitor(*igtimer_);

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
                                         double fluidscale) const
{

  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> fov = FluidField()    .Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      .Interface().ExtractOtherVector(av);

  Teuchos::RCP<Epetra_Vector> fcv = FluidField()    .Interface().ExtractCondVector(fv);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
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
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
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

  Teuchos::RCP<PartialNormF> structureDisp =
    Teuchos::rcp(new PartialNormF("displacement",
                                  *DofRowMap(),
                                  *StructureField().DofRowMap(),
                                  nlParams.get("Norm abs disp", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(structureDisp);

  Teuchos::RCP<PartialNormF> innerFluidVel =
    Teuchos::rcp(new PartialNormF("velocity",
                                  *DofRowMap(),
                                  *FluidField().InnerVelocityRowMap(),
                                  nlParams.get("Norm abs vel", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(innerFluidVel);

  Teuchos::RCP<PartialNormF> fluidPress =
    Teuchos::rcp(new PartialNormF("pressure",
                                  *DofRowMap(),
                                  *FluidField().PressureRowMap(),
                                  nlParams.get("Norm abs pres", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(fluidPress);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                 Teuchos::RCP<const Epetra_Vector>& sx,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax) const
{
  Teuchos::TimeMonitor monitor(*exctracttimer_);

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x,0);
  Teuchos::RCP<const Epetra_Vector> scx = StructureField().Interface().ExtractCondVector(sx);

  // process fluid unknowns
  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);

  // get interface displacement at t(n)
  Teuchos::RCP<Epetra_Vector> dispn = StructureField().Interface().ExtractCondVector(StructureField().Dispn());

  // get interface velocity at t(n)
  Teuchos::RCP<Epetra_Vector> veln = FluidField().Interface().ExtractCondVector(FluidField().Veln());

  // We convert Delta d(n+1,i+1) to Delta u(n+1,i+1) here.
  //fcx->Update(-1./Dt(),*StructToFluid(dispn),1./Dt());
  fcx->Update(-1.,*veln,1./Dt());

  Teuchos::RCP<Epetra_Vector> f = FluidField().Interface().InsertOtherVector(fox);
  FluidField().Interface().InsertCondVector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  // Here we have to add the current structural interface displacement because
  // we solve for the displacement increments on the structural side but for
  // the absolute displacements on the ale side.
  //
  // We have x = d(n+1,i+1) - d(n) here, this makes things quite easy.
  acx->Update(1.0,*StructToAle(dispn),1.0);

  Teuchos::RCP<Epetra_Vector> a = AleField().Interface().InsertOtherVector(aox);
  AleField().Interface().InsertCondVector(acx, a);
  ax = a;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::AddFluidInterface(double scale,
                                               const LINALG::SparseMatrix& fgg,
                                               LINALG::SparseMatrix& s) const
{
  Teuchos::TimeMonitor monitor(*fluidinterfacetimer_);

  const FSI::Coupling& coupsf = StructureFluidCoupling();

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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix>
FSI::MonolithicOverlap::ConvertFigColmap(const LINALG::SparseMatrix& fig,
                                         const std::map<int,int>& convcolmap,
                                         const Epetra_Map& domainmap,
                                         double scale) const
{
  Teuchos::TimeMonitor monitor(*figcolmaptimer_);

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

    for (int j=0; j<NumEntries; ++j)
    {
      int gid = colmap.GID(Indices[j]);
      std::map<int,int>::const_iterator iter = convcolmap.find(gid);
      if (iter==convcolmap.end())
        dserror("gid %d not found in map", gid);
      Indices[j] = iter->second;
    }

    err = convfig->InsertGlobalValues(rowmap.GID(i), NumEntries, Values, Indices);
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
                                         const Epetra_Map& structrowmap) const
{
  Teuchos::TimeMonitor monitor(*fgirowmaptimer_);

  const FSI::Coupling& coupsf = StructureFluidCoupling();

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
    for (int j=0; j<NumEntries; ++j)
    {
      Indices[j] = colmap.GID(Indices[j]);
    }

    // put row into matrix with structure row map
    err = sfgi->InsertGlobalValues(coupsf.MasterDofMap()->GID(i), NumEntries, Values, Indices);
    if (err)
      dserror("InsertGlobalValues error: %d", err);
  }

  sfgi->FillComplete(efgi->DomainMap(),structrowmap);
  sfgi->Scale(scale);

  return Teuchos::rcp(new LINALG::SparseMatrix(sfgi,fgi.ExplicitDirichlet(),fgi.SaveGraph()));
}


#endif
