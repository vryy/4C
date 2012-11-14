/*----------------------------------------------------------------------*/
/*!
\file fsi_lungmonolithic_structuresplit.cpp
\brief Volume-coupled FSI (structure-split)

<pre>
Maintainer: Lena Yoshihara
            yoshihara@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15303
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "fsi_lungmonolithic_structuresplit.H"
#include "fsi_matrixtransform.H"
#include "fsi_lung_overlapprec.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/ad_str_lung.H"
#include "../drt_adapter/ad_fld_lung.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_io/io_control.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::LungMonolithicStructureSplit::LungMonolithicStructureSplit(const Epetra_Comm& comm,
                                                                const Teuchos::ParameterList& timeparams)
  : LungMonolithic(comm,timeparams)
{
  sggtransform_ = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  sgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  sigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  aigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  aiGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmGitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmGGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  addfmGGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  scgitransform_ = Teuchos::rcp(new UTILS::MatrixRowTransform);
  csigtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  caiGtransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupSystem()
{
  GeneralSetup();

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield = Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  //-----------------------------------------------------------------------------
  // create combined map
  //-----------------------------------------------------------------------------

  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(structfield->FSIInterface()->OtherMap());
  vecSpaces.push_back(FluidField().DofRowMap());
  // remaining (not coupled) dofs of ale field
  vecSpaces.push_back(AleField().Interface()->Map(0));
  // additional volume constraints
  vecSpaces.push_back(ConstrMap_);

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);

  FluidField().UseBlockMatrix(false);

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
                                                                *StructureField(),
                                                                FluidField(),
                                                                AleField(),
                                                                true,
                                                                DRT::INPUT::IntegralValue<int>(fsidyn,"SYMMETRICPRECOND"),
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
void FSI::LungMonolithicStructureSplit::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupRHS");

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield = Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());
  ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());

  Teuchos::RCP<Epetra_Vector> structureRHS = Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap()));
  structureRHS->Update(1.0, *StructureField()->RHS(), 1.0, *AddStructRHS_, 0.0);

  double fluidscale = FluidField().ResidualScaling();

  Teuchos::RCP<Epetra_Vector> fluidRHS = Teuchos::rcp(new Epetra_Vector(*FluidField().Discretization()->DofRowMap()));
  fluidRHS->Update(1.0, *FluidField().RHS(), 1.0, *AddFluidRHS_, 0.0);

  SetupVector(f,
              structureRHS,
              fluidRHS,
              AleField().RHS(),
              ConstrRHS_,
              fluidscale);

  if (firstcall)
  {
    // additional rhs term for ALE equations
    // -dt Aig u(n)
    //
    //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
    //
    // And we are concerned with the u(n) part here.

    //--------------------------------------------------------------------------------
    // ale
    //--------------------------------------------------------------------------------
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();
    if (a==Teuchos::null)
      dserror("expect ale block matrix");

    LINALG::SparseMatrix& aig = a->Matrix(0,1);

    Teuchos::RCP<Epetra_Vector> fveln = FluidField().ExtractInterfaceVeln();
    Teuchos::RCP<Epetra_Vector> sveln = FluidToStruct(fveln);
    Teuchos::RCP<Epetra_Vector> aveln = StructToAle(sveln);
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
    aig.Apply(*aveln,*rhs);

    rhs->Scale(-1.*Dt());

    Extractor().AddVector(*rhs,2,f);

    //--------------------------------------------------------------------------------
    // structure
    //--------------------------------------------------------------------------------
    Teuchos::RCP<Epetra_Vector> veln = StructureField()->Interface()->InsertFSICondVector(sveln);
    rhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));

    Teuchos::RCP<LINALG::SparseMatrix> s = StructureField()->SystemMatrix();
    s->Apply(*veln,*rhs);

    Teuchos::RCP<Epetra_Vector> addrhs = Teuchos::rcp(new Epetra_Vector(veln->Map()));
    AddStructConstrMatrix_->Matrix(0,0).Apply(*veln,*addrhs);

    rhs->Update(1.0, *addrhs, 1.0);
    rhs->Scale(-1.*Dt());

    veln = structfield->FSIInterface()->ExtractOtherVector(rhs);
    Extractor().AddVector(*veln,0,f);

    veln = StructureField()->Interface()->ExtractFSICondVector(rhs);
    veln = FluidField().Interface()->InsertFSICondVector(StructToFluid(veln));

    veln->Scale(1./fluidscale);
    Extractor().AddVector(*veln,1,f);

    //--------------------------------------------------------------------------------
    // constraint structure
    //--------------------------------------------------------------------------------
    // split in two blocks according to inner and fsi structure dofs
    Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,StructureField()->Discretization()->Comm()));
    LINALG::MapExtractor extractor;
    extractor.Setup(*ConstrMap_,emptymap,ConstrMap_);

    Teuchos::RCP<LINALG::BlockSparseMatrixBase> constrstructblocks =
      AddStructConstrMatrix_->Matrix(1,0).Split<LINALG::DefaultBlockMatrixStrategy>(*structfield->FSIInterface(),
                                                                                    extractor);
    constrstructblocks->Complete();

    LINALG::SparseMatrix& csig = constrstructblocks->Matrix(0,1);

    rhs = Teuchos::rcp(new Epetra_Vector(csig.RowMap()));
    csig.Apply(*sveln,*rhs);
    rhs->Scale(-1.*Dt());
    Extractor().AddVector(*rhs,3,f);

    //--------------------------------------------------------------------------------
    // constraint ale
    //--------------------------------------------------------------------------------
    LINALG::SparseMatrix& caig = ConstrAleMatrix_->Matrix(0,1);
    caig.Apply(*fveln,*rhs);
    rhs->Scale(-1.*Dt());
    Extractor().AddVector(*rhs,3,f);

    //--------------------------------------------------------------------------------
    // shape derivatives
    //--------------------------------------------------------------------------------
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = fluidfield.ShapeDerivatives();

    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
      LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
      LINALG::SparseMatrix& fmGg = mmm->Matrix(3,1);

      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
      fmig.Apply(*fveln,*rhs);
      veln = FluidField().Interface()->InsertVector(rhs,0);

      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
      fmgg.Apply(*fveln,*rhs);
      FluidField().Interface()->InsertVector(rhs,1,veln);

      rhs = Teuchos::rcp(new Epetra_Vector(fmGg.RowMap()));
      fmGg.Apply(*fveln,*rhs);
      FluidField().Interface()->InsertVector(rhs,3,veln);

      veln->Scale(-1.*Dt());
      Extractor().AddVector(*veln,1,f);

      LINALG::SparseMatrix& afmgg = AddFluidShapeDerivMatrix_->Matrix(1,1);
      LINALG::SparseMatrix& afmGg = AddFluidShapeDerivMatrix_->Matrix(3,1);

      rhs = Teuchos::rcp(new Epetra_Vector(afmgg.RowMap()));
      afmgg.Apply(*fveln,*rhs);
      veln = FluidField().Interface()->InsertVector(rhs,1);

      rhs = Teuchos::rcp(new Epetra_Vector(afmGg.RowMap()));
      afmGg.Apply(*fveln,*rhs);
      FluidField().Interface()->InsertVector(rhs,3,veln);

      veln->Scale(-1.*Dt());
      Extractor().AddVector(*veln,1,f);
    }
  }

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::SetupSystemMatrix");

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here. Extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();

  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField().SystemMatrix();

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  double scale     = FluidField().ResidualScaling();
  double timescale = FluidField().TimeScaling();

  /*----------------------------------------------------------------------*/
  // fluid part

  mat.Assign(1,1,View,*f);

  // fluid linearization with respect to mesh motion block
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField().ShapeDerivatives();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();

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

    LINALG::SparseMatrix& addfmGg = AddFluidShapeDerivMatrix_->Matrix(3,1);
    LINALG::SparseMatrix& addfmGG = AddFluidShapeDerivMatrix_->Matrix(3,3);

    mat.Matrix(1,1).Add(fmig,false,1./timescale,1.0);
    mat.Matrix(1,1).Add(fmgg,false,1./timescale,1.0);
    mat.Matrix(1,1).Add(fmGg,false,1./timescale,1.0);

    (*fmgitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmgi,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false,
                      false);
    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmii,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false,
                      true);
    (*fmGitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmGi,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      mat.Matrix(1,2),
                      false,
                      true);

    (*fmgGtransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmgG,
                      1.,
                      ADAPTER::CouplingMasterConverter(*coupfsout_),
                      mat.Matrix(1,0),
                      true,
                      false);
    (*fmiGtransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmiG,
                      1.,
                      ADAPTER::CouplingMasterConverter(*coupfsout_),
                      mat.Matrix(1,0),
                      true,
                      true);
    (*fmGGtransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmGG,
                      1.,
                      ADAPTER::CouplingMasterConverter(*coupfsout_),
                      mat.Matrix(1,0),
                      true,
                      true);

    mat.Matrix(1,1).Add(addfmGg,false,1./timescale,1.0);

    (*addfmGGtransform_)(AddFluidShapeDerivMatrix_->FullRowMap(),
                         AddFluidShapeDerivMatrix_->FullColMap(),
                         addfmGG,
                         1.,
                         ADAPTER::CouplingMasterConverter(*coupfsout_),
                         mat.Matrix(1,0),
                         true,
                         true);
  }

  /*----------------------------------------------------------------------*/
  // fluid constraint part

  mat.Assign(1,3,View,*FluidConstrMatrix_);

  /*----------------------------------------------------------------------*/
  // structure and additional structure part

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield = Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  LINALG::SparseMatrix s = *StructureField()->SystemMatrix();
  s.UnComplete();
  s.Add(AddStructConstrMatrix_->Matrix(0,0), false, 1.0, 1.0);
  s.Complete();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks =
    s.Split<LINALG::DefaultBlockMatrixStrategy>(*structfield->FSIInterface(),
                                                *structfield->FSIInterface());

  blocks->Complete();

  LINALG::SparseMatrix& sii = blocks->Matrix(0,0);
  LINALG::SparseMatrix& sig = blocks->Matrix(0,1);
  LINALG::SparseMatrix& sgi = blocks->Matrix(1,0);
  LINALG::SparseMatrix& sgg = blocks->Matrix(1,1);

  mat.Assign(0,0,View,sii);

  (*sigtransform_)(blocks->FullRowMap(),
                   blocks->FullColMap(),
                   sig,
                   1./timescale,
                   ADAPTER::CouplingMasterConverter(coupsf),
                   mat.Matrix(0,1));

  (*sggtransform_)(sgg,
                   1./(scale*timescale),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   mat.Matrix(1,1),
                   true,
                   true);

  (*sgitransform_)(sgi,
                   1./scale,
                   ADAPTER::CouplingMasterConverter(coupsf),
                   mat.Matrix(1,0),
                   true);

  /*----------------------------------------------------------------------*/
  // structure constraint part

  // split in two blocks according to inner and fsi structure dofs

  Teuchos::RCP<Epetra_Map> emptymap = Teuchos::rcp(new Epetra_Map(-1,0,NULL,0,StructureField()->Discretization()->Comm()));
  LINALG::MapExtractor extractor;
  extractor.Setup(*ConstrMap_,emptymap,ConstrMap_);

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> structconstrblocks =
    AddStructConstrMatrix_->Matrix(0,1).Split<LINALG::DefaultBlockMatrixStrategy>(extractor,
                                                                                  *structfield->FSIInterface());


  structconstrblocks->Complete();

  LINALG::SparseMatrix& scii = structconstrblocks->Matrix(0,0);
  LINALG::SparseMatrix& scgi = structconstrblocks->Matrix(1,0);

  mat.Assign(0,3,View,scii);

  // add interface part to fluid block

  mat.Matrix(1,3).UnComplete();
  (*scgitransform_)(scgi,
                    1./scale,
                    ADAPTER::CouplingMasterConverter(coupsf),
                    mat.Matrix(1,3),
                    true);

  /*----------------------------------------------------------------------*/
  // ale part

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField().BlockSystemMatrix();

  if (a==Teuchos::null)
    dserror("expect ale block matrix");

  a->Complete();

  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);
  LINALG::SparseMatrix& aiG = a->Matrix(0,3);

  (*aiGtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aiG,
                   1.,
                   ADAPTER::CouplingSlaveConverter(*coupsaout_),
                   mat.Matrix(2,0));
  (*aigtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aig,
                   1./timescale,
                   ADAPTER::CouplingSlaveConverter(*icoupfa_),
                   mat.Matrix(2,1));
  mat.Assign(2,2,View,aii);

  /*----------------------------------------------------------------------*/
  // constraint part -> fluid

  mat.Assign(3,1,View,*ConstrFluidMatrix_);

  /*----------------------------------------------------------------------*/
  // constraint part -> structure
  // split in two blocks according to inner and fsi structure dofs

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> constrstructblocks =
    AddStructConstrMatrix_->Matrix(1,0).Split<LINALG::DefaultBlockMatrixStrategy>(*structfield->FSIInterface(),
                                                                                  extractor);

  constrstructblocks->Complete();

  LINALG::SparseMatrix& csii = constrstructblocks->Matrix(0,0);
  LINALG::SparseMatrix& csig = constrstructblocks->Matrix(0,1);

  mat.Assign(3,0,View,csii);

  mat.Matrix(3,1).UnComplete();
  (*csigtransform_)(*coupsf.MasterDofMap(),
                    csig.ColMap(),
                    csig,
                    1./timescale,
                    ADAPTER::CouplingMasterConverter(coupsf),
                    mat.Matrix(3,1),
                    true,
                    true);

  /*----------------------------------------------------------------------*/
  // constraint part -> ale

  LINALG::SparseMatrix& caiG = ConstrAleMatrix_->Matrix(0,3);
  (*caiGtransform_)(*coupfsout_->MasterDofMap(),
                    caiG.ColMap(),
                    caiG,
                    1.0,
                    ADAPTER::CouplingMasterConverter(*coupfsout_),
                    mat.Matrix(3,0),
                    true,
                    true);

  /*----------------------------------------------------------------------*/
  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicStructureSplit::InitialGuess");

  Teuchos::RCP<Epetra_Vector> ConstraintInitialGuess = Teuchos::rcp(new Epetra_Vector(*ConstrMap_,true));

  SetupVector(*ig,
              StructureField()->InitialGuess(),
              FluidField().InitialGuess(),
              AleField().InitialGuess(),
              ConstraintInitialGuess,
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::SetupVector(Epetra_Vector &f,
                                                    Teuchos::RCP<const Epetra_Vector> sv,
                                                    Teuchos::RCP<const Epetra_Vector> fv,
                                                    Teuchos::RCP<const Epetra_Vector> av,
                                                    Teuchos::RCP<const Epetra_Vector> cv,
                                                    double fluidscale)
{
  // extract the inner and boundary dofs of all three fields

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield = Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  Teuchos::RCP<Epetra_Vector> sov = structfield->FSIInterface()->ExtractOtherVector(sv);

  Teuchos::RCP<Epetra_Vector> aov = AleField().Interface()->ExtractVector(av, 0);

  if (fluidscale!=0)
  {
    // add fluid interface values to structure vector
    Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);
    Teuchos::RCP<Epetra_Vector> modfv = FluidField().Interface()->InsertFSICondVector(StructToFluid(scv));
    modfv->Update(1.0, *fv, 1./fluidscale);

    Extractor().InsertVector(*modfv,1,f);
  }
  else
  {
    Extractor().InsertVector(*fv,1,f);
  }

  Extractor().InsertVector(*sov,0,f);
  Extractor().InsertVector(*aov,2,f);
  Extractor().InsertVector(*cv,3,f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithicStructureSplit::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                            Teuchos::RCP<const Epetra_Vector>& sx,
                                                            Teuchos::RCP<const Epetra_Vector>& fx,
                                                            Teuchos::RCP<const Epetra_Vector>& ax)
{
  const Teuchos::RCP<ADAPTER::StructureLung>& structfield = Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  fx = Extractor().ExtractVector(x,1);

  // process structure unknowns

  Teuchos::RCP<Epetra_Vector> fcx = FluidField().Interface()->ExtractFSICondVector(fx);
  FluidField().VelocityToDisplacement(fcx);
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);

  Teuchos::RCP<Epetra_Vector> s = structfield->FSIInterface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);
  Teuchos::RCP<Epetra_Vector> a = AleField().Interface()->InsertVector(aox,0);
  AleField().Interface()->InsertVector(acx, 1, a);

  Teuchos::RCP<Epetra_Vector> scox = StructureField()->Interface()->ExtractLungASICondVector(sx);
  Teuchos::RCP<Epetra_Vector> acox = StructToAleOutflow(scox);
  AleField().Interface()->InsertVector(acox, 3, a);

  ax = a;
}
