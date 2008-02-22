#ifdef CCADISCRET

#include "fsi_monolithicoverlap.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicOverlap::MonolithicOverlap(Epetra_Comm& comm)
  : Monolithic(comm),
    havefluidstructcolmap_(false),
    havepermfluidstructcolmap_(false)
{
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

  AleField().Evaluate(Teuchos::null);

  // split ale matrix

  Teuchos::RCP<Epetra_CrsMatrix> aii;
  Teuchos::RCP<Epetra_CrsMatrix> aig;

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();
  const LINALG::SparseMatrix* mat_aii = AleField().InteriorMatrixBlock();
  const LINALG::SparseMatrix* mat_aig = AleField().InterfaceMatrixBlock();

  if (a==Teuchos::null or mat_aii==NULL or mat_aig==NULL)
    dserror("expect ale block matrix");

  aii = mat_aii->EpetraMatrix();
  aig = mat_aig->EpetraMatrix();

  // map between ale interface and structure column map

  DRT::Exporter ex(a->FullRowMap(),a->FullColMap(),a->Comm());
  std::map<int,int> alestructcolmap;
  coupsa.FillSlaveToMasterMap(alestructcolmap);
  ex.Export(alestructcolmap);

  aig_ = ConvertFigColmap(aig,alestructcolmap,StructureField().DomainMap(),1.);
  aii_ = aii;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicOverlap::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupRHS(Teuchos::RCP<Epetra_Vector> f) const
{
  Teuchos::TimeMonitor monitor(*rhstimer_);

  SetupVector(f,
              StructureField().RHS(),
              FluidField().RHS(),
              AleField().RHS(),
              FluidField().ResidualScaling());

  // NOX expects a different sign here.
  f->Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupSystemMatrix(Teuchos::RCP<LINALG::BlockSparseMatrixBase> mat) const
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  Teuchos::TimeMonitor monitor(*igtimer_);

  SetupVector(ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicOverlap::SetupVector(Teuchos::RCP<Epetra_Vector> &f,
                                         Teuchos::RCP<const Epetra_Vector> sv,
                                         Teuchos::RCP<const Epetra_Vector> fv,
                                         Teuchos::RCP<const Epetra_Vector> av,
                                         double fluidscale) const
{

  // extract the inner and boundary dofs of all three fields

  Teuchos::RCP<Epetra_Vector> fov = FluidField()    .Interface().ExtractOtherVector(fv);
  Teuchos::RCP<Epetra_Vector> aov = AleField()      .Interface().ExtractOtherVector(av);

  Teuchos::RCP<Epetra_Vector> fcv = FluidField()    .Interface().ExtractCondVector(fv);

  // add fluid interface values to structure vector
  Teuchos::RCP<Epetra_Vector> modsv = StructureField().Interface().InsertCondVector(FluidToStruct(fcv));
  modsv->Update(1.0, *sv, fluidscale);

  Extractor().InsertVector(sv,0,f);
  Extractor().InsertVector(fov,1,f);
  Extractor().InsertVector(aov,2,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::MonolithicOverlap::CreateLinearSystem(ParameterList& nlParams,
                                           const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
                                           NOX::Epetra::Vector& noxSoln,
                                           Teuchos::RCP<NOX::Utils> utils)
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::MonolithicOverlap::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                         Teuchos::RCP<NOX::Epetra::Group> grp)
{
  return Teuchos::null;
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
                                               Teuchos::RCP<Epetra_CrsMatrix> fgg,
                                               Teuchos::RCP<Epetra_CrsMatrix> s) const
{
  Teuchos::TimeMonitor monitor(*fluidinterfacetimer_);

  const FSI::Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<Epetra_CrsMatrix> perm_fgg = coupsf.SlaveToPermSlave(fgg);

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
        err = s->SumIntoGlobalValues(structgid, 1, &value, &index);
        if (err)
          dserror("SumIntoGlobalValues error: %d", err);
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix>
FSI::MonolithicOverlap::ConvertFigColmap(Teuchos::RCP<Epetra_CrsMatrix> fig,
                                         const std::map<int,int>& convcolmap,
                                         const Epetra_Map& domainmap,
                                         double scale) const
{
  Teuchos::TimeMonitor monitor(*figcolmaptimer_);

  const Epetra_Map& rowmap = fig->RowMap();
  const Epetra_Map& colmap = fig->ColMap();

  Teuchos::RCP<Epetra_CrsMatrix> convfig = Teuchos::rcp(new Epetra_CrsMatrix(Copy,rowmap,
                                                                             fig->MaxNumEntries()));

  int rows = fig->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = fig->ExtractMyRowView(i, NumEntries, Values, Indices);
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

  return convfig;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_CrsMatrix>
FSI::MonolithicOverlap::ConvertFgiRowmap(Teuchos::RCP<Epetra_CrsMatrix> fgi,
                                         double scale,
                                         const Epetra_Map& structrowmap) const
{
  Teuchos::TimeMonitor monitor(*fgirowmaptimer_);

  const FSI::Coupling& coupsf = StructureFluidCoupling();

  // redistribute fluid matrix to match distribution of structure matrix

  fgi = coupsf.SlaveToPermSlave(fgi);

  // create new matrix with structure row map and fill it

  Teuchos::RCP<Epetra_CrsMatrix> sfgi = Teuchos::rcp(new Epetra_CrsMatrix(Copy,structrowmap,fgi->MaxNumEntries()));

  //const Epetra_Map& rowmap = fgi->RowMap();
  const Epetra_Map& colmap = fgi->ColMap();

  int rows = fgi->NumMyRows();
  for (int i=0; i<rows; ++i)
  {
    int NumEntries;
    double *Values;
    int *Indices;
    int err = fgi->ExtractMyRowView(i, NumEntries, Values, Indices);
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

  sfgi->FillComplete(fgi->DomainMap(),structrowmap);
  sfgi->Scale(scale);
  return sfgi;
}


#endif
