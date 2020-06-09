/*----------------------------------------------------------------------*/
/*! \file
\brief Volume-coupled FSI (base class)

\level 3

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/
#include "fsi_lungmonolithic.H"
#include "fsi_lung_overlapprec.H"
#include "fsi_overlapprec_amgnxn.H"
#include "../drt_adapter/ad_str_lung.H"
#include "../drt_adapter/ad_fld_lung.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "fsi_statustest.H"
#include "../drt_io/io_control.H"
#include "fsi_monolithic_linearsystem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_constraint/constraintdofset.H"
#include "../drt_structure/stru_aux.H"
#include "../linalg/linalg_utils_densematrix_communication.H"
#include "../linalg/linalg_utils_sparse_algebra_manipulation.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_ale_fsi.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::LungMonolithic::LungMonolithic(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : BlockMonolithic(comm, timeparams)
{
  icoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsaout_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfsout_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfaout_ = Teuchos::rcp(new ADAPTER::Coupling());

  //-----------------------------------------------------------------------------
  // additional fluid-structure volume constraints
  //-----------------------------------------------------------------------------

  // Since the current design of the constraint manager and fsi
  // algorithms complicates the neat combination of both
  // (e.g. concerning the question of who owns what actually) in this
  // special application, the general functionality of the constraint
  // manager is included here on the algorithm level (as far as
  // needed).

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());
  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  // consistency check: all dofs contained in ale(fluid)-structure coupling need to
  // be part of the structure volume constraint, too. this needs to be checked because during
  // SetupSystemMatrix, we rely on this information!
  const Teuchos::RCP<const Epetra_Map> asimap = StructureField()->Interface()->LungASICondMap();
  for (int i = 0; i < asimap->NumMyElements(); ++i)
  {
    if (structfield->LungConstrMap()->LID(asimap->GID(i)) == -1)
      dserror("dof of asi coupling is not contained in enclosing boundary");
  }

  std::set<int> FluidLungVolConIDs;
  std::set<int> StructLungVolConIDs;
  int FluidMinLungVolConID;
  int StructMinLungVolConID;

  structfield->ListLungVolCons(StructLungVolConIDs, StructMinLungVolConID);
  fluidfield->ListLungVolCons(FluidLungVolConIDs, FluidMinLungVolConID);

  // We want to be sure that both fluid and structure fields hold the
  // same constraint IDs. After all, every airway outlet needs to be
  // coupled to one structural volume. Therefore, merely comparing the
  // overall number and the minimum constraint ID is not sufficient here.

  for (std::set<int>::iterator iter = FluidLungVolConIDs.begin(); iter != FluidLungVolConIDs.end();
       ++iter)
  {
    if (StructLungVolConIDs.find(*iter) == StructLungVolConIDs.end())
      dserror("No matching in fluid and structure lung volume constraints");
  }

  NumConstrID_ = FluidLungVolConIDs.size();

  ConstrDofSet_ = Teuchos::rcp(new ::UTILS::ConstraintDofSet());
  ConstrDofSet_->AssignDegreesOfFreedom(FluidField()->Discretization(), NumConstrID_, 0);

  // The "OffsetID" is used during the evaluation of constraints on
  // the element level. For assembly of the constraint parts, the gid
  // of the constraint dof (= Lagrange multiplier) needs to be known.
  //
  // gid = current constraint ID - minimum constraint ID + first gid of all constraints
  //                             \__________________________  ________________________/
  //                                                        \/
  //                                                   - OffsetID_
  //
  // By including the minimum constraint ID, one allows also to define
  // a set of constraints not starting from 1 in the input file.
  // Since the "OffsetID" is subtracted later on, we save its negative
  // value here.

  OffsetID_ = FluidMinLungVolConID - ConstrDofSet_->FirstGID();
  ConstrMap_ = Teuchos::rcp(new Epetra_Map(*(ConstrDofSet_->DofRowMap())));

  // build an all reduced version of the constraintmap, since sometimes all processors
  // have to know all values of the constraints and Lagrange multipliers
  RedConstrMap_ = LINALG::AllreduceEMap(*ConstrMap_);

  // create importer
  ConstrImport_ = Teuchos::rcp(new Epetra_Export(*RedConstrMap_, *ConstrMap_));

  // initialize associated matrices and vectors

  // NOTE: everything that is determined in the fluid adapter needs to
  // be based on the fluid dofmap, i.e. also the ale related stuff!
  // Corresponding matrices then need to be transformed to the ale
  // dofmap using corresponding matrix transformators.

  LagrMultVec_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));
  LagrMultVecOld_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));
  IncLagrMultVec_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  // build merged structure dof map
  Teuchos::RCP<Epetra_Map> FullStructDofMap =
      LINALG::MergeMap(*StructureField()->DofRowMap(), *ConstrMap_, false);
  LINALG::MapExtractor StructConstrExtractor(
      *FullStructDofMap, ConstrMap_, StructureField()->DofRowMap());

  AddStructConstrMatrix_ =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          StructConstrExtractor, StructConstrExtractor, 81, false, true));

  AddFluidShapeDerivMatrix_ =
      Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
          *FluidField()->Interface(), *FluidField()->Interface(), 108, false, true));
  FluidConstrMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *FluidField()->Discretization()->DofRowMap(), NumConstrID_, false, true));
  ConstrFluidMatrix_ = Teuchos::rcp(new LINALG::SparseMatrix(
      *ConstrMap_, FluidField()->Discretization()->DofRowMap()->NumGlobalElements(), false, true));

  // additional "ale" matrices filled in the fluid elements
  Teuchos::RCP<Epetra_Map> emptymap =
      Teuchos::rcp(new Epetra_Map(-1, 0, NULL, 0, FluidField()->Discretization()->Comm()));
  LINALG::MapExtractor constrextractor;
  constrextractor.Setup(*ConstrMap_, emptymap, ConstrMap_);
  AleConstrMatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      constrextractor, *FluidField()->Interface(), 108, false, true));
  ConstrAleMatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      *FluidField()->Interface(), constrextractor, 108, false, true));

  AddStructRHS_ =
      Teuchos::rcp(new Epetra_Vector(*StructureField()->Discretization()->DofRowMap(), true));
  AddFluidRHS_ =
      Teuchos::rcp(new Epetra_Vector(*FluidField()->Discretization()->DofRowMap(), true));
  ConstrRHS_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  OldVols_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));
  CurrVols_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));
  SignVolsRed_ = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_, true));
  dVstruct_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  OldFlowRates_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));
  CurrFlowRates_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));
  dVfluid_ = Teuchos::rcp(new Epetra_Vector(*ConstrMap_, true));

  // time integration factor for flow rates
  theta_ = 0.5;

  // determine initial volumes of parenchyma balloons
  Teuchos::RCP<Epetra_Vector> OldVolsRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));

  structfield->InitializeVolCon(OldVolsRed, SignVolsRed_, OffsetID_);

  OldVols_->PutScalar(0.0);
  OldVols_->Export(*OldVolsRed, *ConstrImport_, Add);

  // determine initial flow rates at outlets
  Teuchos::RCP<Epetra_Vector> OldFlowRatesRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));

  fluidfield->InitializeVolCon(OldFlowRatesRed, OffsetID_);

  OldFlowRates_->PutScalar(0.0);
  OldFlowRates_->Export(*OldFlowRatesRed, *ConstrImport_, Add);

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::GeneralSetup()
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  linearsolverstrategy_ =
      DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  SetDefaultParameters(fsidyn, NOXParameterList());

  // ToDo: Set more detailed convergence tolerances like in standard FSI
  // additionally set tolerance for volume constraint
  NOXParameterList().set("Norm abs vol constr", fsimono.get<double>("CONVTOL"));

  //-----------------------------------------------------------------------------
  // ordinary fsi coupling
  //-----------------------------------------------------------------------------

  // right now we use matching meshes at the interface

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();

  const int ndim = DRT::Problem::Instance()->NDim();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->FSICondMap(), *FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa_->SetupConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->FSICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = AleField()->Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField()->Discretization(), *AleField()->Discretization(),
      *fluidnodemap, *alenodemap, ndim);

  FluidField()->SetMeshMap(coupfa.MasterDofMap());

  aleresidual_ = Teuchos::rcp(new Epetra_Vector(*AleField()->Interface()->Map(0)));

  //-----------------------------------------------------------------------------
  // additional coupling of structure and ale field at the outflow boundary
  //-----------------------------------------------------------------------------

  // coupling of structure and ale dofs at airway outflow
  coupsaout_->SetupConstrainedConditionCoupling(*StructureField()->Discretization(),
      StructureField()->Interface()->LungASICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->LungASICondMap(), "StructAleCoupling", "FSICoupling", ndim);
  if (coupsaout_->MasterDofMap()->NumGlobalElements() == 0)
    dserror("No nodes in matching structure ale interface. Empty coupling condition?");

  // coupling of fluid and structure dofs at airway outflow
  coupfsout_->SetupConstrainedConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->LungASICondMap(), *StructureField()->Discretization(),
      StructureField()->Interface()->LungASICondMap(), "StructAleCoupling", "FSICoupling", ndim);
  if (coupfsout_->MasterDofMap()->NumGlobalElements() == 0)
    dserror("No nodes in matching structure ale/fluid interface. Empty coupling condition?");

  // coupling of fluid and ale dofs at airway outflow
  coupfaout_->SetupConstrainedConditionCoupling(*FluidField()->Discretization(),
      FluidField()->Interface()->LungASICondMap(), *AleField()->Discretization(),
      AleField()->Interface()->LungASICondMap(), "StructAleCoupling", "FSICoupling", ndim);
  if (coupfaout_->MasterDofMap()->NumGlobalElements() == 0)
    dserror("No nodes in matching ale fluid ouflow interface. Empty coupling condition?");

  //-----------------------------------------------------------------------------
  // enable output of changes in volumes in text file
  //-----------------------------------------------------------------------------

  if (Comm().MyPID() == 0)
  {
    std::string outputprefix = DRT::Problem::Instance()->OutputControlFile()->NewOutputFileName();
    std::string dfluidfilename;
    std::string dstructfilename;
    std::string absstructfilename;
    std::string absfluidfilename;
    size_t posn = outputprefix.rfind('-');
    if (posn != std::string::npos)
    {
      std::string number = outputprefix.substr(posn + 1);
      std::string prefix = outputprefix.substr(0, posn);
      std::ostringstream sf;
      sf << prefix << "_dVfluid"
         << "-" << number << ".txt";
      dfluidfilename = sf.str();
      std::ostringstream ss;
      ss << prefix << "_dVstruct"
         << "-" << number << ".txt";
      dstructfilename = ss.str();
      std::ostringstream sas;
      sas << prefix << "_absVstruct"
          << "-" << number << ".txt";
      absstructfilename = sas.str();
    }
    else
    {
      std::ostringstream sf;
      sf << outputprefix << "_dVfluid.txt";
      dfluidfilename = sf.str();
      std::ostringstream ss;
      ss << outputprefix << "_dVstruct.txt";
      dstructfilename = ss.str();
      std::ostringstream sas;
      sas << outputprefix << "_absVstruct.txt";
      absstructfilename = sas.str();
    }

    outfluiddvol_.open(dfluidfilename.c_str());
    outstructdvol_.open(dstructfilename.c_str());
    outstructabsvol_.open(absstructfilename.c_str());
  }

  writerestartevery_ = fsidyn.get<int>("RESTARTEVRY");

  // ToDo: Setup the monolithic DBC map extractor and use only this to handle DBCs in Matrix and RHS
  dbcmaps_ = Teuchos::null;

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::LungMonolithic::StructToAleOutflow(
    Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsaout_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  //-----------------------------------------------------------------------------
  // evaluation of all fields
  //-----------------------------------------------------------------------------

  FSI::Monolithic::Evaluate(x);

  //-----------------------------------------------------------------------------
  // evaluation of lung volume constraints
  //-----------------------------------------------------------------------------

  if (x != Teuchos::null)
  {
    // extract sum of all iterative increments of lagrange multipliers
    // in this time step (this is what we get from NOX)
    IncLagrMultVec_ = Extractor().ExtractVector(x, 3);

    // update current lagrange multipliers
    LagrMultVec_->Update(1.0, *LagrMultVecOld_, 1.0, *IncLagrMultVec_, 0.0);
  }

  //-----------------------------------------------------------------------------
  // structure part

  // create redundant vectors
  Teuchos::RCP<Epetra_Vector> LagrMultVecRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
  LINALG::Export(*LagrMultVec_, *LagrMultVecRed);
  Teuchos::RCP<Epetra_Vector> CurrVolsRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());
  CurrVolsRed->PutScalar(0.0);
  AddStructRHS_->PutScalar(0.0);
  AddStructConstrMatrix_->Zero();

  structfield->EvaluateVolCon(
      AddStructConstrMatrix_, AddStructRHS_, CurrVolsRed, SignVolsRed_, LagrMultVecRed, OffsetID_);

  // Export redundant vector into distributed one
  CurrVols_->PutScalar(0.0);
  CurrVols_->Export(*CurrVolsRed, *ConstrImport_, Add);


  // negative sign (for shift to rhs) is already taken into account!
  dVstruct_->Update(1.0, *CurrVols_, -1.0, *OldVols_, 0.0);
  ConstrRHS_->Update(-1.0, *dVstruct_, 0.0);

  //-----------------------------------------------------------------------------
  // fluid/ale part

  Teuchos::RCP<ADAPTER::FluidLung> fluidfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::FluidLung>(FluidField());

  // create redundant vector
  Teuchos::RCP<Epetra_Vector> CurrFlowRatesRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));

  CurrFlowRatesRed->PutScalar(0.0);
  AddFluidRHS_->PutScalar(0.0);
  AddFluidShapeDerivMatrix_->Zero();
  FluidConstrMatrix_->Zero();
  ConstrFluidMatrix_->Zero();
  AleConstrMatrix_->Zero();
  ConstrAleMatrix_->Zero();

  const double dt = Dt();
  const double dttheta = dt * theta_;

  fluidfield->EvaluateVolCon(AddFluidShapeDerivMatrix_, FluidConstrMatrix_, ConstrFluidMatrix_,
      AleConstrMatrix_, ConstrAleMatrix_, AddFluidRHS_, CurrFlowRatesRed, LagrMultVecRed, OffsetID_,
      dttheta);

  // Export redundant vector into distributed one
  CurrFlowRates_->PutScalar(0.0);
  CurrFlowRates_->Export(*CurrFlowRatesRed, *ConstrImport_, Add);

  // negative sign (for shift to rhs) is already taken into account!
  dVfluid_->Update(dt * theta_, *CurrFlowRates_, dt * (1.0 - theta_), *OldFlowRates_, 0.0);
  ConstrRHS_->Update(1.0, *dVfluid_, 1.0);

  //   std::cout << "CurrFlowRates_:\n" << *CurrFlowRates_ << std::endl;
  //   std::cout << "CurrVols_:\n" << *CurrVols_ << std::endl;
  //   std::cout << "LagrMultVec_:\n" << *LagrMultVec_ << std::endl;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  // should we scale the system?
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3, 0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 3).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(3, 2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b, 2);

    if (sx->Multiply(1.0, *srowsum_, *sx, 0.0)) dserror("structure scaling failed");
    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0)) dserror("ale scaling failed");

    Extractor().InsertVector(*sx, 0, b);
    Extractor().InsertVector(*ax, 2, b);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::UnscaleSolution(
    LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> sy = Extractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ay = Extractor().ExtractVector(x, 2);

    if (sy->Multiply(1.0, *scolsum_, *sy, 0.0)) dserror("structure scaling failed");
    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) dserror("ale scaling failed");

    Extractor().InsertVector(*sy, 0, x);
    Extractor().InsertVector(*ay, 2, x);

    Teuchos::RCP<Epetra_Vector> sx = Extractor().ExtractVector(b, 0);
    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b, 2);

    if (sx->ReciprocalMultiply(1.0, *srowsum_, *sx, 0.0)) dserror("structure scaling failed");
    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) dserror("ale scaling failed");

    Extractor().InsertVector(*sx, 0, b);
    Extractor().InsertVector(*ax, 2, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0, 0).EpetraMatrix();
    srowsum_->Reciprocal(*srowsum_);
    scolsum_->Reciprocal(*scolsum_);
    if (A->LeftScale(*srowsum_) or A->RightScale(*scolsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0, 3).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(3, 0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2, 2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(2, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2, 3).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1, 2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(3, 2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");
  }

  // very simple hack just to see the linear solution

  Epetra_Vector r(b.Map());
  mat.Apply(x, r);
  r.Update(1., b, 1.);

  Teuchos::RCP<Epetra_Vector> sr = Extractor().ExtractVector(r, 0);
  Teuchos::RCP<Epetra_Vector> fr = Extractor().ExtractVector(r, 1);
  Teuchos::RCP<Epetra_Vector> ar = Extractor().ExtractVector(r, 2);
  Teuchos::RCP<Epetra_Vector> cr = Extractor().ExtractVector(r, 3);

  // increment additional ale residual
  aleresidual_->Update(-1., *ar, 0.);

  std::ios_base::fmtflags flags = Utils()->out().flags();

  double n, ns, nf, na, nc;
  r.Norm2(&n);
  sr->Norm2(&ns);
  fr->Norm2(&nf);
  ar->Norm2(&na);
  cr->Norm2(&nc);
  Utils()->out() << std::scientific << "\nlinear solver quality:\n"
                 << "L_2-norms:\n"
                 << END_COLOR "   |r|=" YELLOW << n << END_COLOR "   |rs|=" YELLOW << ns
                 << END_COLOR "   |rf|=" YELLOW << nf << END_COLOR "   |ra|=" YELLOW << na
                 << END_COLOR "   |rc|=" YELLOW << nc << END_COLOR "\n";
  r.NormInf(&n);
  sr->NormInf(&ns);
  fr->NormInf(&nf);
  ar->NormInf(&na);
  cr->NormInf(&nc);
  Utils()->out() << "L_inf-norms:\n"
                 << END_COLOR "   |r|=" YELLOW << n << END_COLOR "   |rs|=" YELLOW << ns
                 << END_COLOR "   |rf|=" YELLOW << nf << END_COLOR "   |ra|=" YELLOW << na
                 << END_COLOR "   |rc|=" YELLOW << nc << END_COLOR "\n";

  Utils()->out().flags(flags);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> FSI::LungMonolithic::CreateLinearSystem(
    Teuchos::ParameterList& nlParams, NOX::Epetra::Vector& noxSoln, Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList* lsParams = NULL;

  // in case of nonlinCG the linear solver list is somewhere else
  if (dirParams.get("Method", "User Defined") == "User Defined")
    lsParams = &(newtonParams.sublist("Linear Solver"));
  else if (dirParams.get("Method", "User Defined") == "NonlinearCG")
    lsParams = &(dirParams.sublist("Nonlinear CG").sublist("Linear Solver"));
  else
    dserror("Unknown nonlinear method");

  NOX::Epetra::Interface::Jacobian* iJac = this;
  NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = systemmatrix_;
  const Teuchos::RCP<Epetra_Operator> M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
    case INPAR::FSI::PreconditionedKrylov:
    case INPAR::FSI::AMGnxn:
      linSys = Teuchos::rcp(new  // NOX::Epetra::LinearSystemAztecOO(
          FSI::MonolithicLinearSystem(printParams, *lsParams, Teuchos::rcp(iJac, false), J,
              Teuchos::rcp(iPrec, false), M, noxSoln));
      break;
    case INPAR::FSI::FSIAMG:
    default:
      dserror("unsupported linear block solver strategy: fsiamg");
      break;
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo> FSI::LungMonolithic::CreateStatusTest(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  //   Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
  //   Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
  //   combo->addStatusTest(update);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  //   combo->addStatusTest(update);
  combo->addStatusTest(maxiters);

  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));

  // setup tests for structural displacements

  Teuchos::RCP<NOX::StatusTest::Combo> structcombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> structureDisp = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "displacement", Extractor(), 0, nlParams.get<double>("Norm abs disp"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> structureDispUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("displacement update", Extractor(), 0,
          nlParams.get<double>("Norm abs disp"), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(structureDisp);
  structcombo->addStatusTest(structureDisp);
  // structcombo->addStatusTest(structureDispUpdate);

  converged->addStatusTest(structcombo);

  // setup tests for interface

  std::vector<Teuchos::RCP<const Epetra_Map>> interface;
  interface.push_back(FluidField()->Interface()->FSICondMap());
  interface.push_back(Teuchos::null);
  LINALG::MultiMapExtractor interfaceextract(*DofRowMap(), interface);

  Teuchos::RCP<NOX::StatusTest::Combo> interfacecombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "interface", interfaceextract, 0, nlParams.get<double>("Norm abs vel"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("interface update", interfaceextract, 0,
          nlParams.get<double>("Norm abs vel"), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(interfaceTest);
  interfacecombo->addStatusTest(interfaceTest);
  // interfacecombo->addStatusTest(interfaceTestUpdate);

  converged->addStatusTest(interfacecombo);

  // setup tests for fluid velocities

  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(FluidField()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidvelextract(*DofRowMap(), fluidvel);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidvelcombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "velocity", fluidvelextract, 0, nlParams.get<double>("Norm abs vel"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("velocity update", fluidvelextract, 0,
          nlParams.get<double>("Norm abs vel"), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(innerFluidVel);
  fluidvelcombo->addStatusTest(innerFluidVel);
  // fluidvelcombo->addStatusTest(innerFluidVelUpdate);

  converged->addStatusTest(fluidvelcombo);

  // setup tests for fluid pressure

  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(FluidField()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidpressextract(*DofRowMap(), fluidpress);

  Teuchos::RCP<NOX::StatusTest::Combo> fluidpresscombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "pressure", fluidpressextract, 0, nlParams.get<double>("Norm abs pres"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("pressure update", fluidpressextract, 0,
          nlParams.get<double>("Norm abs pres"), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(fluidPress);
  fluidpresscombo->addStatusTest(fluidPress);
  // fluidpresscombo->addStatusTest(fluidPressUpdate);

  converged->addStatusTest(fluidpresscombo);


  // setup tests for volume constraint

  std::vector<Teuchos::RCP<const Epetra_Map>> volconstr;
  volconstr.push_back(ConstrMap_);
  volconstr.push_back(Teuchos::null);
  LINALG::MultiMapExtractor volconstrextract(*DofRowMap(), volconstr);

  Teuchos::RCP<NOX::StatusTest::Combo> volconstrcombo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> VolConstr = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "volume constraint", volconstrextract, 0, nlParams.get<double>("Norm abs vol constr"),
      NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> VolConstrUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("volume constraint update", volconstrextract, 0,
          nlParams.get<double>("Norm abs vol constr"), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(VolConstr);
  volconstrcombo->addStatusTest(VolConstr);
  // volconstrcombo->addStatusTest(VolConstrUpdate);


  converged->addStatusTest(volconstrcombo);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::Update()
{
  FSI::BlockMonolithic::Update();

  // update fluid flow rates and structure volumes and lagrange multipliers
  OldVols_->Update(1.0, *CurrVols_, 0.0);
  OldFlowRates_->Update(1.0, *CurrFlowRates_, 0.0);
  LagrMultVecOld_->Update(1.0, *LagrMultVec_, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField()->Output();

  // additional output of volume constraint related forces
  //   ADAPTER::StructureLung& structfield =
  //   dynamic_cast<ADAPTER::StructureLung&>(StructureField());
  //   structfield.OutputForces(AddStructRHS_);

  // Write history vectors in case of restart
  // This is done using the structure DiscretizationWriter, hence it
  // is placed in between the output of the single fields.
  if (writerestartevery_ and (Step() % writerestartevery_ == 0))
  {
    const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
        Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());
    Teuchos::RCP<Epetra_Vector> OldFlowRatesRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
    LINALG::Export(*OldFlowRates_, *OldFlowRatesRed);
    Teuchos::RCP<Epetra_Vector> OldVolsRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
    LINALG::Export(*OldVols_, *OldVolsRed);
    Teuchos::RCP<Epetra_Vector> LagrMultVecOldRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
    LINALG::Export(*LagrMultVecOld_, *LagrMultVecOldRed);
    structfield->WriteVolConRestart(OldFlowRatesRed, OldVolsRed, LagrMultVecOldRed);
  }

  FluidField()->Output();

  // additional output of volume constraint related forces
  //   ADAPTER::FluidLung& fluidfield = dynamic_cast<ADAPTER::FluidLung&>(FluidField());
  //   fluidfield->OutputForces(AddFluidRHS_);

  AleField()->Output();

  // output of volumes for visualization (e.g. gnuplot)

  Teuchos::RCP<Epetra_Vector> dVfluidRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
  LINALG::Export(*dVfluid_, *dVfluidRed);

  if (Comm().MyPID() == 0)
  {
    outfluiddvol_ << Step();
    for (int i = 0; i < dVfluidRed->MyLength(); ++i)
    {
      outfluiddvol_ << "\t" << (*dVfluidRed)[i];
    }
    outfluiddvol_ << "\n" << std::flush;
  }

  Teuchos::RCP<Epetra_Vector> dVstructRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
  LINALG::Export(*dVstruct_, *dVstructRed);

  if (Comm().MyPID() == 0)
  {
    outstructdvol_ << Step();
    for (int i = 0; i < dVstructRed->MyLength(); ++i)
    {
      outstructdvol_ << "\t" << (*dVstructRed)[i];
    }
    outstructdvol_ << "\n" << std::flush;
  }

  Teuchos::RCP<Epetra_Vector> VstructRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
  LINALG::Export(*CurrVols_, *VstructRed);

  if (Comm().MyPID() == 0)
  {
    outstructabsvol_ << Step();
    for (int i = 0; i < VstructRed->MyLength(); ++i)
    {
      outstructabsvol_ << "\t" << (*VstructRed)[i];
    }
    outstructabsvol_ << "\n" << std::flush;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::ReadRestart(int step)
{
  FSI::Monolithic::ReadRestart(step);

  const Teuchos::RCP<ADAPTER::StructureLung>& structfield =
      Teuchos::rcp_dynamic_cast<ADAPTER::StructureLung>(StructureField());

  Teuchos::RCP<Epetra_Vector> OldFlowRatesRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
  Teuchos::RCP<Epetra_Vector> OldVolsRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));
  Teuchos::RCP<Epetra_Vector> OldLagrMultRed = Teuchos::rcp(new Epetra_Vector(*RedConstrMap_));

  structfield->ReadVolConRestart(step, OldFlowRatesRed, OldVolsRed, OldLagrMultRed);

  // Export redundant vector into distributed one
  OldVols_->PutScalar(0.0);
  OldVols_->Export(*OldVolsRed, *ConstrImport_, Insert);
  CurrVols_->Update(1.0, *OldVols_, 0.0);
  OldFlowRates_->PutScalar(0.0);
  OldFlowRates_->Export(*OldFlowRatesRed, *ConstrImport_, Insert);
  CurrFlowRates_->Update(1.0, *OldFlowRates_, 0.0);
  LagrMultVecOld_->PutScalar(0.0);
  LagrMultVecOld_->Export(*OldLagrMultRed, *ConstrImport_, Insert);
  LagrMultVec_->Update(1.0, *LagrMultVecOld_, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::PrepareTimeStep()
{
  FSI::BlockMonolithic::PrepareTimeStep();

  // additional lung volume constraint stuff

  // Update of Lagrange multipliers, current volumes and flow rates is
  // not necessary here, since these values are already equal to the
  // "old" ones (cf. Update()). Note that we assume a constant
  // predictor here!

  IncLagrMultVec_->PutScalar(0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> FSI::LungMonolithic::SystemMatrix() const
{
  return systemmatrix_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::LungMonolithic::CreateSystemMatrix(bool structuresplit)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // get the PCITER from inputfile
  std::vector<int> pciter;
  std::vector<double> pcomega;
  std::vector<int> spciter;
  std::vector<double> spcomega;
  std::vector<int> fpciter;
  std::vector<double> fpcomega;
  std::vector<int> apciter;
  std::vector<double> apcomega;
  {
    std::string word1;
    std::string word2;
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "PCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "PCOMEGA"));
      while (pciterstream >> word1) pciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) pcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "STRUCTPCITER"));
      std::istringstream pcomegastream(
          Teuchos::getNumericStringParameter(fsimono, "STRUCTPCOMEGA"));
      while (pciterstream >> word1) spciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) spcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCOMEGA"));
      while (pciterstream >> word1) fpciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) fpcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "ALEPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "ALEPCOMEGA"));
      while (pciterstream >> word1) apciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) apcomega.push_back(std::atof(word2.c_str()));
    }
  }

  //-----------------------------------------------------------------------------
  // create block system matrix
  //-----------------------------------------------------------------------------

  switch (linearsolverstrategy_)
  {
    case INPAR::FSI::PreconditionedKrylov:
      systemmatrix_ = Teuchos::rcp(
          new LungOverlappingBlockMatrix(Extractor(), *StructureField(), *FluidField(), *AleField(),
              structuresplit, DRT::INPUT::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"),
              pcomega[0], pciter[0], spcomega[0], spciter[0], fpcomega[0], fpciter[0], apcomega[0],
              apciter[0], DRT::Problem::Instance()->ErrorFile()->Handle()));
      break;
    case INPAR::FSI::AMGnxn:
    {
      // Parse BLOCKSMOOTHER list
      std::vector<std::string> blocksmoother;
      std::string word;
      std::istringstream blocksmootherstream(
          Teuchos::getNumericStringParameter(fsimono, "BLOCKSMOOTHER"));
      while (blocksmootherstream >> word) blocksmoother.push_back(word);
      // We assume that the xml file is given in the first position of the BLOCKSMOOTHER list
      std::string amgnxn_xml = "none";
      if ((int)blocksmoother.size() > 0)
        amgnxn_xml = blocksmoother[0];
      else
        dserror("Not found xml file in the first position of the BLOCKSMOOTHER list");
      systemmatrix_ = Teuchos::rcp(new OverlappingBlockMatrixAMGnxn(Extractor(), *StructureField(),
          *FluidField(), *AleField(), structuresplit, amgnxn_xml,
          DRT::Problem::Instance()->ErrorFile()->Handle(), "LungFSI"));
    }
    break;
    case INPAR::FSI::FSIAMG:
    default:
      dserror("Unsupported type of monolithic solver");
      break;
  }
}
