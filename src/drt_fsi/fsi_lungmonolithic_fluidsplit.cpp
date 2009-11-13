#ifdef CCADISCRET

#include "fsi_lungmonolithic_fluidsplit.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_adapter/adapter_structure_lung.H"
#include "../drt_adapter/adapter_fluid_lung.H"

#define FLUIDSPLITAMG

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::LungMonolithicFluidSplit::LungMonolithicFluidSplit(Epetra_Comm& comm)
  : LungMonolithic(comm)
{
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupSystem()
{
  GeneralSetup();

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  //-----------------------------------------------------------------------------
  // create combined map
  //-----------------------------------------------------------------------------

  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField().DofRowMap());
#ifdef FLUIDSPLITAMG
  vecSpaces.push_back(FluidField().DofRowMap());
#else
  ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());
  vecSpaces.push_back(fluidfield.FSIInterface().OtherMap());
#endif
  vecSpaces.push_back(AleField().Interface().OtherMap());
  vecSpaces.push_back(ConstrMap_);

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  FluidField().UseBlockMatrix(true);

  // build ale system matrix in splitted system
  AleField().BuildSystemMatrix(false);

  // get the PCITER from inputfile
  vector<int> pciter;
  vector<double> pcomega;
  vector<int> spciter;
  vector<double> spcomega;
  vector<int> fpciter;
  vector<double> fpcomega;
  vector<int> apciter;
  vector<double> apcomega;
  {
    int    word1;
    double word2;
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"PCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"PCOMEGA"));
      while (pciterstream >> word1)
        pciter.push_back(word1);
      while (pcomegastream >> word2)
        pcomega.push_back(word2);
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"STRUCTPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"STRUCTPCOMEGA"));
      while (pciterstream >> word1)
        spciter.push_back(word1);
      while (pcomegastream >> word2)
        spcomega.push_back(word2);
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"FLUIDPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"FLUIDPCOMEGA"));
      while (pciterstream >> word1)
        fpciter.push_back(word1);
      while (pcomegastream >> word2)
        fpcomega.push_back(word2);
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsidyn,"ALEPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsidyn,"ALEPCOMEGA"));
      while (pciterstream >> word1)
        apciter.push_back(word1);
      while (pcomegastream >> word2)
        apcomega.push_back(word2);
    }
  }

  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------

  switch(linearsolverstrategy_)
  {
  case INPAR::FSI::PreconditionedKrylov:
    systemmatrix_ = Teuchos::rcp(new LungOverlappingBlockMatrix(Extractor(),
                                                                StructureField(),
                                                                FluidField(),
                                                                AleField(),
                                                                false,
                                                                Teuchos::getIntegralValue<int>(fsidyn,"SYMMETRICPRECOND"),
                                                                pcomega[0],
                                                                pciter[0],
                                                                spcomega[0],
                                                                spciter[0],
                                                                fpcomega[0],
                                                                fpciter[0],
                                                                apcomega[0],
                                                                apciter[0],
                                                                DRT::Problem::Instance()->ErrorFile()->Handle()));

  break;
  case INPAR::FSI::FSIAMG:
  default:
    dserror("Unsupported type of monolithic solver");
  break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupRHS");

  Teuchos::RCP<Epetra_Vector> structureRHS = rcp(new Epetra_Vector(*StructureField().Discretization()->DofRowMap()));
  structureRHS->Update(1.0, *StructureField().RHS(), 1.0, *AddStructRHS_, 0.0);

  Teuchos::RCP<Epetra_Vector> fluidRHS = rcp(new Epetra_Vector(*FluidField().Discretization()->DofRowMap()));
  fluidRHS->Update(1.0, *FluidField().RHS(), 1.0, *AddFluidRHS_, 0.0);

  double scale = FluidField().ResidualScaling();

  SetupVector(f,
              structureRHS,
              fluidRHS,
              AleField().RHS(),
              ConstrRHS_,
              scale);

  if (firstcall)
  {
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();

    LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
    LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);
    LINALG::SparseMatrix& fGg = blockf->Matrix(3,1);

    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    double timescale = FluidField().TimeScaling();

    ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(*fluidfield.InnerSplit().FullMap()));
    Teuchos::RCP<Epetra_Vector> tmp = Teuchos::rcp(new Epetra_Vector(fig.RowMap()));
    fig.Apply(*fveln,*tmp);
    fluidfield.InnerSplit().InsertOtherVector(tmp,rhs);
    tmp = Teuchos::rcp(new Epetra_Vector(fGg.RowMap()));
    fGg.Apply(*fveln,*tmp);
    fluidfield.InnerSplit().InsertCondVector(tmp,rhs);
    rhs->Scale(timescale*Dt());
#ifdef FLUIDSPLITAMG
    rhs = fluidfield.FSIInterface().InsertOtherVector(rhs);
#endif
    Extractor().AddVector(*rhs,1,f);

    rhs = Teuchos::rcp(new Epetra_Vector(fgg.RowMap()));
    fgg.Apply(*fveln,*rhs);
    rhs->Scale(scale*timescale*Dt());
    rhs = FluidToStruct(rhs);
    rhs = StructureField().Interface().InsertFSICondVector(rhs);
    Extractor().AddVector(*rhs,0,f);
  }

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::SetupSystemMatrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupsa = StructureAleCoupling();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

  /*----------------------------------------------------------------------*/

  double scale     = FluidField().ResidualScaling();
  double timescale = FluidField().TimeScaling();

  /*----------------------------------------------------------------------*/
  // structure part

  LINALG::SparseMatrix s = *StructureField().SystemMatrix();
  s.UnComplete();
  s.Add(AddStructConstrMatrix_->Matrix(0,0), false, 1.0, 1.0);
  //s.Complete();
  mat.Assign(0,0,View,s);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  LINALG::SparseMatrix scon = AddStructConstrMatrix_->Matrix(0,1);
  mat.Assign(0,3,View,scon);

  /*----------------------------------------------------------------------*/
  // fluid part

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField().BlockSystemMatrix();

  LINALG::SparseMatrix& fii = blockf->Matrix(0,0);
  LINALG::SparseMatrix& fig = blockf->Matrix(0,1);
  LINALG::SparseMatrix& fiG = blockf->Matrix(0,3);
  LINALG::SparseMatrix& fgi = blockf->Matrix(1,0);
  LINALG::SparseMatrix& fgg = blockf->Matrix(1,1);
  LINALG::SparseMatrix& fgG = blockf->Matrix(1,3);
  LINALG::SparseMatrix& fGi = blockf->Matrix(3,0);
  LINALG::SparseMatrix& fGg = blockf->Matrix(3,1);
  LINALG::SparseMatrix& fGG = blockf->Matrix(3,3);

  //mat.Matrix(0,0).UnComplete();
  fggtransform_(fgg,
                scale*timescale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                ADAPTER::Coupling::SlaveConverter(coupsf),
                mat.Matrix(0,0),
                true,
                true);
  mat.Matrix(0,0).Complete();

  fgitransform_(fgi,
                scale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                mat.Matrix(0,1));
  fgGtransform_(fgG,
                scale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                mat.Matrix(0,1),
                true);

  figtransform_(blockf->FullRowMap(),
                blockf->FullColMap(),
                fig,
                timescale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                mat.Matrix(1,0));
  fGgtransform_(blockf->FullRowMap(),
                blockf->FullColMap(),
                fGg,
                timescale,
                ADAPTER::Coupling::SlaveConverter(coupsf),
                mat.Matrix(1,0),
                true,
                true);

  mat.Matrix(1,1).Add(fii, false, 1.0, 0.0);
  mat.Matrix(1,1).Add(fiG, false, 1.0, 1.0);
  mat.Matrix(1,1).Add(fGi, false, 1.0, 1.0);
  mat.Matrix(1,1).Add(fGG, false, 1.0, 1.0);

#ifdef FLUIDSPLITAMG
  Teuchos::RCP<LINALG::SparseMatrix> eye = LINALG::Eye(*FluidField().Interface().FSICondMap());
  mat.Matrix(1,1).Add(*eye,false,1.0,1.0);
#endif

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();

  if (mmm!=Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    LINALG::SparseMatrix& fmiG = mmm->Matrix(0,3);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
    LINALG::SparseMatrix& fmgG = mmm->Matrix(1,3);
    LINALG::SparseMatrix& fmGi = mmm->Matrix(3,0);
    LINALG::SparseMatrix& fmGg = mmm->Matrix(3,1);
    LINALG::SparseMatrix& fmGG = mmm->Matrix(3,3);

    // We cannot copy the pressure value. It is not used anyway. So no exact
    // match here.
    fmiitransform_(mmm->FullRowMap(),
                   mmm->FullColMap(),
                   fmii,
                   1.,
                   ADAPTER::Coupling::MasterConverter(coupfa),
                   mat.Matrix(1,2),
                   false,
                   false);
    fmGitransform_(mmm->FullRowMap(),
                   mmm->FullColMap(),
                   fmGi,
                   1.,
                   ADAPTER::Coupling::MasterConverter(coupfa),
                   mat.Matrix(1,2),
                   false,
                   true);

    fmgitransform_(fmgi,
                   scale,
                   ADAPTER::Coupling::SlaveConverter(coupsf),
                   ADAPTER::Coupling::MasterConverter(coupfa),
                   mat.Matrix(0,2),
                   false,
                   false);

    fmigtransform_(mmm->FullRowMap(),
                   mmm->FullColMap(),
                   fmig,
                   1.,
                   ADAPTER::Coupling::SlaveConverter(coupsf),
                   mat.Matrix(1,0),
                   false,
                   true);
    fmGgtransform_(mmm->FullRowMap(),
                   mmm->FullColMap(),
                   fmGg,
                   1.,
                   ADAPTER::Coupling::SlaveConverter(coupsf),
                   mat.Matrix(1,0),
                   false,
                   true);

    fmggtransform_(fmgg,
                   scale,
                   ADAPTER::Coupling::SlaveConverter(coupsf),
                   ADAPTER::Coupling::SlaveConverter(coupsf),
                   mat.Matrix(0,0),
                   false,
                   true);

    fmiGtransform_(mmm->FullRowMap(),
                   mmm->FullColMap(),
                   fmiG,
                   1.,
                   ADAPTER::Coupling::SlaveConverter(coupfsout_),
                   mat.Matrix(1,0),
                   false,
                   true);
    fmGGtransform_(mmm->FullRowMap(),
                   mmm->FullColMap(),
                   fmGG,
                   1.,
                   ADAPTER::Coupling::SlaveConverter(coupfsout_),
                   mat.Matrix(1,0),
                   false,
                   true);

    fmgGtransform_(fmgG,
                   scale,
                   ADAPTER::Coupling::SlaveConverter(coupsf),
                   ADAPTER::Coupling::SlaveConverter(coupfsout_),
                   mat.Matrix(0,0),
                   false,
                   true);

    LINALG::SparseMatrix& addfmgg = AddFluidShapeDerivMatrix_->Matrix(1,1);
    LINALG::SparseMatrix& addfmgG = AddFluidShapeDerivMatrix_->Matrix(1,3);

    LINALG::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->Matrix(3,1);
    LINALG::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->Matrix(3,3);

    addfmgGtransform_(addfmgG,
                      scale,
                      ADAPTER::Coupling::SlaveConverter(coupsf),
                      ADAPTER::Coupling::MasterConverter(coupfsout_),
                      mat.Matrix(0,0),
                      false,
                      true);

    addfmGGtransform_(AddFluidShapeDerivMatrix_->FullRowMap(),
                      AddFluidShapeDerivMatrix_->FullColMap(),
                      addfmGG,
                      1.,
                      ADAPTER::Coupling::MasterConverter(coupfsout_),
                      mat.Matrix(1,0),
                      false,
                      true);

    addfmggtransform_(addfmgg,
                      scale,
                      ADAPTER::Coupling::SlaveConverter(coupsf),
                      ADAPTER::Coupling::MasterConverter(coupsf),
                      mat.Matrix(0,0),
                      false,
                      true);

    addfmGgtransform_(AddFluidShapeDerivMatrix_->FullRowMap(),
                      AddFluidShapeDerivMatrix_->FullColMap(),
                      addfmGg,
                      1.,
                      ADAPTER::Coupling::MasterConverter(coupsf),
                      mat.Matrix(1,0),
                      false,
                      true);
  }

  /*----------------------------------------------------------------------*/
  // fluid constraint part

  // split in two blocks according to inner and fsi dofs

  ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());

  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,FluidField().Discretization()->Comm()));
  LINALG::MapExtractor extractor;
  extractor.Setup(*ConstrMap_,emptymap,ConstrMap_);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> fluidconstrblocks =
    FluidConstrMatrix_->Split<LINALG::DefaultBlockMatrixStrategy>(extractor,
                                                                  fluidfield.FSIInterface());

  fluidconstrblocks->Complete();

  LINALG::SparseMatrix& fcii = fluidconstrblocks->Matrix(0,0);
  LINALG::SparseMatrix& fcgi = fluidconstrblocks->Matrix(1,0);

  mat.Matrix(1,3).Zero();
  mat.Matrix(1,3).Add(fcii,false,1.0,0.0);

  // add interface part to structure block

  mat.Matrix(0,3).UnComplete();
  fcgitransform_(fcgi,
                 scale,
                 ADAPTER::Coupling::SlaveConverter(coupsf),
                 mat.Matrix(0,3),
                 true);
  mat.Matrix(0,3).Complete(*ConstrMap_,*StructureField().Discretization()->DofRowMap());

  /*----------------------------------------------------------------------*/
  // ale part

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

  if (a==Teuchos::null)
    dserror("expect ale block matrix");

  a->Complete();

  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);
  LINALG::SparseMatrix& aiG = a->Matrix(0,3);

  mat.Assign(2,2,View,aii);

  aigtransform_(a->FullRowMap(),
                a->FullColMap(),
                aig,
                1.,
                ADAPTER::Coupling::SlaveConverter(coupsa),
                mat.Matrix(2,0));

  aiGtransform_(a->FullRowMap(),
                a->FullColMap(),
                aiG,
                1.,
                ADAPTER::Coupling::SlaveConverter(coupsaout_),
                mat.Matrix(2,0),
                false,
                true);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure

  mat.Assign(3,0,View,AddStructConstrMatrix_->Matrix(1,0));

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid

  // split in two blocks according to inner and fsi structure dofs

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> constrfluidblocks =
    ConstrFluidMatrix_->Split<LINALG::DefaultBlockMatrixStrategy>(fluidfield.FSIInterface(),
                                                                  extractor);
  constrfluidblocks->Complete();

  LINALG::SparseMatrix& cfii = constrfluidblocks->Matrix(0,0);
  LINALG::SparseMatrix& cfig = constrfluidblocks->Matrix(0,1);

  mat.Matrix(3,1).Zero();
  mat.Matrix(3,1).Add(cfii,false,1.0,0.0);

  cfigtransform_(constrfluidblocks->FullRowMap(),
                 constrfluidblocks->FullColMap(),
                 cfig,
                 timescale,
                 ADAPTER::Coupling::MasterConverter(coupsf),
                 mat.Matrix(3,0),
                 false,
                 true);

  /*----------------------------------------------------------------------*/
  // constraint part -> "ale"

  LINALG::SparseMatrix& caig = ConstrAleMatrix_->Matrix(0,1);

  caigtransform_(*coupsf.MasterDofMap(),
                 caig.ColMap(),
                 caig,
                 1.0,
                 ADAPTER::Coupling::MasterConverter(coupsf),
                 mat.Matrix(3,0),
                 false,
                 true);

  LINALG::SparseMatrix& caiG = ConstrAleMatrix_->Matrix(0,3);
  caiGtransform_(*coupsaout_.SlaveDofMap(),
                 caiG.ColMap(),
                 caiG,
                 1.0,
                 ADAPTER::Coupling::SlaveConverter(coupsaout_),
                 mat.Matrix(3,0),
                 false,
                 true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFluidSplit::InitialGuess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess = rcp(new Epetra_Vector(*ConstrMap_,true));

  SetupVector(*ig,
              StructureField().InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              ConstraintInitialGuess,
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::SetupVector(Epetra_Vector &f,
                                                    Teuchos::RCP<const Epetra_Vector> sv,
                                                    Teuchos::RCP<const Epetra_Vector> fv,
                                                    Teuchos::RCP<const Epetra_Vector> av,
                                                    Teuchos::RCP<const Epetra_Vector> cv,
                                                    double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());
  Teuchos::RCP<Epetra_Vector> fov = fluidfield.FSIInterface().ExtractOtherVector(fv);
#ifdef FLUIDSPLITAMG
  fov = fluidfield.FSIInterface().InsertOtherVector(fov);
#endif
  Teuchos::RCP<Epetra_Vector> aov = AleField().Interface().ExtractOtherVector(av);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> fcv = FluidField().Interface().ExtractFSICondVector(fv);
    Teuchos::RCP<Epetra_Vector> modsv = StructureField().Interface().InsertFSICondVector(FluidToStruct(fcv));
    modsv->Update(1.0, *sv, fluidscale);

    Extractor().InsertVector(*modsv,0,f);
  }
  else
  {
    Extractor().InsertVector(*sv,0,f);
  }

  Extractor().InsertVector(*fov,1,f);
  Extractor().InsertVector(*aov,2,f);
  Extractor().InsertVector(*cv,3,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicFluidSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                        Teuchos::RCP<const Epetra_Vector>& sx,
                                                        Teuchos::RCP<const Epetra_Vector>& fx,
                                                        Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::LungMonolithicFluidSplit::ExtractFieldVectors");

  ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());

  // We have overlap at the interface. Thus we need the interface part of the
  // structure vector and append it to the fluid and ale vector. (With the
  // right translation.)

  sx = Extractor().ExtractVector(x,0);
  Teuchos::RCP<const Epetra_Vector> scx = StructureField().Interface().ExtractFSICondVector(sx);

  // process fluid unknowns

  Teuchos::RCP<const Epetra_Vector> fox = Extractor().ExtractVector(x,1);
#ifdef FLUIDSPLITAMG
  fox = fluidfield.FSIInterface().ExtractOtherVector(fox);
#endif
  Teuchos::RCP<Epetra_Vector> fcx = StructToFluid(scx);

  FluidField().DisplacementToVelocity(fcx);

  Teuchos::RCP<Epetra_Vector> f = fluidfield.FSIInterface().InsertOtherVector(fox);
  FluidField().Interface().InsertFSICondVector(fcx, f);
  fx = f;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);
  Teuchos::RCP<Epetra_Vector> a = AleField().Interface().InsertOtherVector(aox);
  AleField().Interface().InsertVector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox = StructureField().Interface().ExtractLungASICondVector(sx);
  Teuchos::RCP<Epetra_Vector> acox = StructToAleOutflow(scox);
  AleField().Interface().InsertVector(acox, 3, a);

  ax = a;
}


#endif
