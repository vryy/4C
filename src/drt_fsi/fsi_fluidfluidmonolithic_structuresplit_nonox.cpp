/*!----------------------------------------------------------------------
\file fsi_fluidfluidmonolithic_structuresplit_nonox.cpp
\brief Control routine for monolithic fluid-fluid-fsi (structuresplit)
using XFEM

\level 3

<pre>
\maintainer  Benedikt Schott
             schott@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15241
</pre>

*----------------------------------------------------------------------*/
#include "fsi_fluidfluidmonolithic_structuresplit_nonox.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "fsi_matrixtransform.H"
#include "fsi_debugwriter.H"
#include "fsi_statustest.H"
#include "fsi_monolithic_linearsystem.H"

#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_ale/ale_utils_mapextractor.H"
#include "../drt_adapter/ad_fld_fluid_fluid_fsi.H"
#include "../drt_adapter/ad_ale_xffsi.H"
#include "../drt_inpar/inpar_xfem.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_io/io_pstream.H"

#include "../drt_inpar/inpar_ale.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidFluidMonolithicStructureSplitNoNOX::FluidFluidMonolithicStructureSplitNoNOX(const Epetra_Comm& comm,
                                                                                      const Teuchos::ParameterList& timeparams)
  : MonolithicNoNOX(comm,timeparams)
{
  // Throw an error if there are DBCs on structural interface DOFs.
  std::vector<Teuchos::RCP<const Epetra_Map> > intersectionmaps;
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  Teuchos::RCP<Epetra_Map> intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);

  if (intersectionmap->NumGlobalElements() != 0)
  {
    std::stringstream errormsg;
    errormsg  << "  +---------------------------------------------------------------------------------------------+" << std::endl
              << "  |                DIRICHLET BOUNDARY CONDITIONS ON SLAVE SIDE OF FSI INTERFACE                 |" << std::endl
              << "  +---------------------------------------------------------------------------------------------+" << std::endl
              << "  | NOTE: The slave side of the interface is not allowed to carry Dirichlet boundary conditions.|" << std::endl
              << "  |                                                                                             |" << std::endl
              << "  | This is a structure split scheme. Hence, master and slave field are chosen as follows:      |" << std::endl
              << "  |     MASTER  = FLUID                                                                         |" << std::endl
              << "  |     SLAVE   = STRUCTURE                                                                     |" << std::endl
              << "  |                                                                                             |" << std::endl
              << "  | Dirichlet boundary conditions were detected on slave interface degrees of freedom. Please   |" << std::endl
              << "  | remove Dirichlet boundary conditions from the slave side of the FSI interface.              |" << std::endl
              << "  | Only the master side of the FSI interface is allowed to carry Dirichlet boundary conditions.|" << std::endl
              << "  +---------------------------------------------------------------------------------------------+" << std::endl;

    dserror(errormsg.str());
  }

#ifdef DEBUG
  // check if removing Dirichlet conditions was successful
  intersectionmaps.resize(0);
  intersectionmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  intersectionmaps.push_back(StructureField()->Interface()->FSICondMap());
  intersectionmap = LINALG::MultiMapExtractor::IntersectMaps(intersectionmaps);
  if (intersectionmap->NumGlobalElements() != 0)
    dserror("Could not remove structural interface Dirichlet conditions from structure DBC map.");
#endif

  sggtransform_  = Teuchos::rcp(new UTILS::MatrixRowColTransform);
  sgitransform_  = Teuchos::rcp(new UTILS::MatrixRowTransform);
  sigtransform_  = Teuchos::rcp(new UTILS::MatrixColTransform);
  aigtransform_  = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new UTILS::MatrixColTransform);
  fsaigtransform_= Teuchos::rcp(new UTILS::MatrixColTransform);
  fsmgitransform_= Teuchos::rcp(new UTILS::MatrixColTransform);

  // Recovering of Lagrange multiplier happens on structure field
  lambda_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap()));
  ddiinc_ = Teuchos::null;
  solipre_= Teuchos::null;
  ddginc_ = Teuchos::null;
  solgpre_= Teuchos::null;
  ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->Interface()->FSICondMap(),true));

  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupSystem()
{
  // setup coupling
  FSI::MonolithicNoNOX::SetupSystem();

  // create combined map
  CreateCombinedDofRowMap();

  // Use normal matrix for fluid equations but build (splitted) mesh movement
  // linearization (if requested in the input file)
  FluidField()->UseBlockMatrix(false);

  // Use splitted structure matrix
  StructureField()->UseBlockMatrix();

  // build ale system matrix in splitted system
  AleField()->CreateSystemMatrix(AleField()->Interface());

  /*----------------------------------------------------------------------*/
  // initialize systemmatrix_
  systemmatrix_ = Teuchos::rcp(
    new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      Extractor(),
      Extractor(),
      81,
      false,
      true
      )
    );
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupRHS");

  SetupVector(f,
              StructureField()->RHS(),
              FluidField()->RHS(),
              AleField()->RHS(),
              FluidField()->ResidualScaling());

  if (firstcall)
  {
     // get structure matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocks = StructureField()->BlockSystemMatrix();
    if (blocks==Teuchos::null)
      dserror("expect structure block matrix");

    // get fluid shape derivatives matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();

    //get ale matrix
    Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocka = AleField()->BlockSystemMatrix();
    if (blocka==Teuchos::null)
      dserror("expect ale block matrix");

    // extract structure and ale submatrices
    LINALG::SparseMatrix& sig = blocks->Matrix(0,1);  // S_{I\Gamma}
    LINALG::SparseMatrix& sgg = blocks->Matrix(1,1);  // S_{\Gamma\Gamma}
    LINALG::SparseMatrix& aig = blocka->Matrix(0,1);  // A_{I\Gamma}

    // get time integration parameters of structure an fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = StructureField()->TimIntParam();
    const double ftiparam = FluidField()->TimIntParam();

    // some scaling factors for fluid
    const double scale     = FluidField()->ResidualScaling();
    const double dt        = FluidField()->Dt();

    // some often re-used vectors
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;

    // ---------- inner structural DOFs
    /* The following terms are added to the inner structural DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  S_{I\Gamma} * \Delta d_{\Gamma,p}
     *
     * (2)  - dt * S_{I\Gamma} * u^{n}_{\Gamma}
     *
     *  Remarks on all terms:
     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
     *
     */

    // ----------addressing term 1
    rhs = Teuchos::rcp(new Epetra_Vector(sig.RowMap(),true));
    sig.Apply(*ddgpred_,*rhs);

    Extractor().AddVector(*rhs,0,f);

    // ----------addressing term 2
    rhs = Teuchos::rcp(new Epetra_Vector(sig.RowMap(),true));
    Teuchos::RCP<Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();
    sig.Apply(*FluidToStruct(fveln),*rhs);
    rhs->Scale(-dt);

    Extractor().AddVector(*rhs,0,f);

    // ---------- end of inner structural DOFs

    // we need a temporary vector with the whole fluid dofs where we
    // can insert the embedded dofrowmap into it
    // kruse 30.04. --> we don't, see adapter
    //Teuchos::RCP<Epetra_Vector> fluidfluidtmp = LINALG::CreateVector(*FluidField()->DofRowMap(),true);

    // ---------- inner fluid DOFs
    /* The following terms are added to the inner fluid DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  - dt * F^{G}_{I\Gamma} * u^{n}_{\Gamma}
     *
     */
    // ----------addressing term 1
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap(),true));
      fmig.Apply(*fveln,*rhs);
      rhs->Scale(-dt);

      rhs = FluidField()->Interface()->InsertOtherVector(rhs);
      Extractor().AddVector(*rhs,1,f);
    }

    // ---------- end of inner fluid DOFs

    // ---------- interface fluid DOFs
    /* The following terms are added to the interface fluid DOFs of right hand side:
     *
     * rhs_firstnewtonstep =
     *
     * (1)  - dt * F^{G}_{\Gamma\Gamma} * u^{n}_{\Gamma}
     *
     * (2)  - dt * (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * u^{n}_{\Gamma}
     *
     * (3)  + (1-ftiparam)/(1-stiparam) * S_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
     *
     * Remarks on all terms:
     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
     *
     */

    // ----------addressing term 1
    if (mmm!=Teuchos::null)
    {
      LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap(),true));
      fmgg.Apply(*fveln,*rhs);
      rhs->Scale(-dt);
      rhs = FluidField()->Interface()->InsertFSICondVector(rhs);
      Extractor().AddVector(*rhs,1,f);
    }

    // ----------addressing term 2
    rhs = Teuchos::rcp(new Epetra_Vector(sgg.RowMap(),true));
    sgg.Apply(*FluidToStruct(fveln),*rhs);
    rhs->Scale(-dt * (1.-ftiparam) / ((1-stiparam) * scale));

    rhs = StructToFluid(rhs);
    rhs = FluidField()->Interface()->InsertFSICondVector(rhs);
    Extractor().AddVector(*rhs,1,f);

    // ----------addressing term 3
    rhs = Teuchos::rcp(new Epetra_Vector(sgg.RowMap(),true));
    sgg.Apply(*ddgpred_,*rhs);
    rhs->Scale((1.-ftiparam) / ((1-stiparam) * scale));

    rhs = StructToFluid(rhs);
    rhs = FluidField()->Interface()->InsertFSICondVector(rhs);
    Extractor().AddVector(*rhs,1,f);

    // ----------addressing term 1
    rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap(),true));

    aig.Apply(*FluidToAleInterface(fveln),*rhs);
    rhs->Scale(-dt);
    Extractor().AddVector(*rhs,2,f);
    // ---------- end of inner ale DOFs
  }

  // -----------------------------------------------------
  // Reset quantities for previous iteration step since they still store values from the last time step
  ddiinc_ = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true);
  solipre_ = Teuchos::null;
  ddginc_ = LINALG::CreateVector(*StructureField()->Interface()->FSICondMap(),true);
  solgpre_ = Teuchos::null;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupSystemMatrix()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupSystemMatrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  const ADAPTER::Coupling& coupfa = FluidAleCoupling();
  const ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> s = StructureField()->BlockSystemMatrix();
  if (s==Teuchos::null)
    dserror("expect structure block matrix");
  Teuchos::RCP<LINALG::SparseMatrix> f = FluidField()->SystemMatrix();
  if (f==Teuchos::null)
    dserror("expect fluid matrix");
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> a = AleField()->BlockSystemMatrix();
  if (a==Teuchos::null)
    dserror("expect ale block matrix");

  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);

  // scaling factors for fluid
  const double scale     = FluidField()->ResidualScaling();
  const double timescale = FluidField()->TimeScaling();

  // get time integration parameters of structure an fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // store parts of structural matrix to know them in the next iteration as previous iteration matricesu
  sgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,0)));
  sggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(s->Matrix(1,1)));

  /*----------------------------------------------------------------------*/
  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  systemmatrix_->Assign(0,0,LINALG::View,s->Matrix(0,0));

  (*sigtransform_)(s->FullRowMap(),
                   s->FullColMap(),
                   s->Matrix(0,1),
                   1./timescale,
                   ADAPTER::CouplingMasterConverter(coupsf),
                   systemmatrix_->Matrix(0,1));
  (*sggtransform_)(s->Matrix(1,1),
                   ((1.0-ftiparam)/(1.0-stiparam))*(1./(scale*timescale)),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   *f,
                   true,
                   true);

  Teuchos::RCP<LINALG::SparseMatrix> lsgi = Teuchos::rcp(new LINALG::SparseMatrix(f->RowMap(),81,false));
  (*sgitransform_)(s->Matrix(1,0),
                   ((1.0-ftiparam)/(1.0-stiparam))*(1./scale),
                   ADAPTER::CouplingMasterConverter(coupsf),
                   *lsgi);

  lsgi->Complete(s->Matrix(1,0).DomainMap(),f->RangeMap());

  //systemmatrix_->Assign(1,1,View,*f);
  systemmatrix_->Assign(1,0,LINALG::View,*lsgi);

  (*aigtransform_)(a->FullRowMap(),
                   a->FullColMap(),
                   aig,
                   1./timescale,
                   ADAPTER::CouplingSlaveConverter(icoupfa),
                   systemmatrix_->Matrix(2,1));

                   systemmatrix_->Assign(2,2,LINALG::View,aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();
  if (mmm!=Teuchos::null)
  {
    LINALG::SparseMatrix& fmii = mmm->Matrix(0,0);
    LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
    LINALG::SparseMatrix& fmgi = mmm->Matrix(1,0);
    LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);

    systemmatrix_->Matrix(1,1).Add(fmgg,false,1./timescale,1.0);
    systemmatrix_->Matrix(1,1).Add(fmig,false,1./timescale,1.0);

    Teuchos::RCP<LINALG::SparseMatrix> lfmgi = Teuchos::rcp(new LINALG::SparseMatrix(f->RowMap(),81,false));
    (*fmgitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmgi,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      //systemmatrix_->Matrix(1,2),
                      *lfmgi,
                      false,
                      false);

    (*fmiitransform_)(mmm->FullRowMap(),
                      mmm->FullColMap(),
                      fmii,
                      1.,
                      ADAPTER::CouplingMasterConverter(coupfa),
                      *lfmgi,
                      false,
                      true);

    lfmgi->Complete(aii.DomainMap(),f->RangeMap());

    systemmatrix_->Assign(1,2,LINALG::View,*lfmgi);
  }

  f->Complete();

  // finally assign fluid block
  systemmatrix_->Assign(1,1,LINALG::View,*f);

  // done. make sure all blocks are filled.
  systemmatrix_->Complete();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicStructureSplitNoNOX::InitialGuess");

  SetupVector(*ig,
              StructureField()->InitialGuess(),
              FluidField()->InitialGuess(),
              AleField()->InitialGuess(),
              0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::ScaleSystem(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  //should we scale the system?
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsimono,"INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(0,0).EpetraMatrix();
    srowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    scolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*srowsum_);
    A->InvColSums(*scolsum_);
    if (A->LeftScale(*srowsum_) or
        A->RightScale(*scolsum_) or
        mat.Matrix(0,1).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->LeftScale(*srowsum_) or
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(),false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
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

/*----------------------------------------------------------------------*
 |  map containing the dofs with Dirichlet BC
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FSI::FluidFluidMonolithicStructureSplitNoNOX::CombinedDBCMap()
{
  //Create a combined map vector with the 3 field DBC maps
  std::vector<Teuchos::RCP<const Epetra_Map> > alldbcmaps;

  // structure DBC
  alldbcmaps.push_back(StructureField()->GetDBCMapExtractor()->CondMap());
  // fluid DBC
  alldbcmaps.push_back(FluidField()->GetDBCMapExtractor()->CondMap());
  // ALE-DBC
  std::vector<Teuchos::RCP<const Epetra_Map> > aleintersectionmaps;
  aleintersectionmaps.push_back(AleField()->GetDBCMapExtractor()->CondMap());
  aleintersectionmaps.push_back(AleField()->Interface()->OtherMap());
  Teuchos::RCP<Epetra_Map> aleintersectionmap = LINALG::MultiMapExtractor::IntersectMaps(aleintersectionmaps);
  alldbcmaps.push_back(aleintersectionmap);

  // Merge the maps
  Teuchos::RCP<Epetra_Map> alldbcmap = LINALG::MultiMapExtractor::MergeMaps(alldbcmaps);

  return alldbcmap;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::UnscaleSolution(LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)DRT::INPUT::IntegralValue<int>(fsimono,"INFNORMSCALING");

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
        mat.Matrix(1,0).EpetraMatrix()->RightScale(*scolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->RightScale(*scolsum_))
      dserror("structure scaling failed");

    A = mat.Matrix(2,2).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or
        A->RightScale(*acolsum_) or
        mat.Matrix(2,0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(2,1).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0,2).EpetraMatrix()->RightScale(*acolsum_) or
        mat.Matrix(1,2).EpetraMatrix()->RightScale(*acolsum_))
      dserror("ale scaling failed");
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupVector(Epetra_Vector &f,
                                                          Teuchos::RCP<const Epetra_Vector> sv,
                                                          Teuchos::RCP<const Epetra_Vector> fv,
                                                          Teuchos::RCP<const Epetra_Vector> av,
                                                          double fluidscale)
{
  // get time integration parameters of structure an fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  // structure inner
  Teuchos::RCP<Epetra_Vector> sov = StructureField()->Interface()->ExtractOtherVector(sv);

  // ale inner
  Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractOtherVector(av);

  // add fluid interface values to structure vector
  // scv: structure fsi dofs
  Teuchos::RCP<Epetra_Vector> scv = StructureField()->Interface()->ExtractFSICondVector(sv);

  if (fabs(fluidscale) > EPS15)
  {
    // modfv: whole embedded fluid map but entries at fsi dofs
    Teuchos::RCP<Epetra_Vector> modfv = FluidField()->Interface()->InsertFSICondVector(StructToFluid(scv));

    // modfv = modfv * 1/fluidscale * d/b
    modfv->Scale( (1./fluidscale)*(1.0-ftiparam)/(1.0-stiparam));

    // add contribution of Lagrange multiplier from previous time step
    if (lambda_ != Teuchos::null)
    {
      Teuchos::RCP<Epetra_Vector> lambdaglobal = FluidField()->Interface()->InsertFSICondVector(StructToFluid(lambda_));
      modfv->Update((-ftiparam+stiparam*(1.0-ftiparam)/(1.0-stiparam))/fluidscale, *lambdaglobal, 1.0);
    }

    modfv->Update(1.0,*fv,1.0);
    Extractor().InsertVector(*modfv,1,f);
  }
  else
  {
    Extractor().InsertVector(*fv,1,f);
  }

  Extractor().InsertVector(*sov,0,f);
  Extractor().InsertVector(*aov,2,f);
}

/*----------------------------------------------------------------------
* - Called from Evaluate() method in Newton-loop with x=x_sum_
*   (increment sum)
* - Field contributions sx,fx,ax are recovered from x
----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                                       Teuchos::RCP<const Epetra_Vector>& sx,
                                                                       Teuchos::RCP<const Epetra_Vector>& fx,
                                                                       Teuchos::RCP<const Epetra_Vector>& ax)
{
#ifdef DEBUG
  if(ddgpred_ == Teuchos::null) { dserror("Vector 'ddgpred_' has not been initialized properly."); }
#endif

  // ----------------------
  // process fluid unknowns
  fx = Extractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> fcx = FluidField()->Interface()->ExtractFSICondVector(fx);

  // ----------------------
  // process ale unknowns
  Teuchos::RCP<Epetra_Vector> fcxforale = Teuchos::rcp(new Epetra_Vector(*fcx));
  FluidField()->VelocityToDisplacement(fcxforale);
  Teuchos::RCP<Epetra_Vector> acx = FluidToStruct(fcxforale);
  acx = StructToAle(acx);

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,2);
  // Teuchos::RCP<Epetra_Vector> acx = StructToAle(scx);

  Teuchos::RCP<Epetra_Vector> a = AleField()->Interface()->InsertOtherVector(aox);
  AleField()->Interface()->InsertFSICondVector(acx, a);

  ax = a;

  // ----------------------
  // process structure unknowns
  Teuchos::RCP<const Epetra_Vector> sox = Extractor().ExtractVector(x,0);

  //Teuchos::RCP<Epetra_Vector> scx = FluidToStruct(fcx);
  // convert ALE interface displacements to structure interface displacements
  Teuchos::RCP<Epetra_Vector> scx = AleToStruct(acx);
  scx->Update(-1.0, *ddgpred_, 1.0);

  Teuchos::RCP<Epetra_Vector> s = StructureField()->Interface()->InsertOtherVector(sox);
  StructureField()->Interface()->InsertFSICondVector(scx, s);
  sx = s;

  // IMPORTANT:
  // you can get these increments in a similar way like in fluidsplit, without saving the
  // previous variables. This is important if you are considering the whole term of lagrange-multiplier.
  // ----------------------
  // Store field vectors to know them later on as previous quantities
  if (solipre_ != Teuchos::null)
    ddiinc_->Update(1.0, *sox, -1.0, *solipre_, 0.0);  // compute current iteration increment
  else
    ddiinc_ = Teuchos::rcp(new Epetra_Vector(*sox));   // first iteration increment

  solipre_ = sox;                                      // store current step increment

  if (solgpre_ != Teuchos::null)
    ddginc_->Update(1.0, *scx, -1.0, *solgpre_, 0.0);  // compute current iteration increment
  else
    ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx));   // first iteration increment

  solgpre_ = scx;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::Output()
{
  StructureField()->Output();

  // output Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = StructureField()->Interface()->InsertFSICondVector(lambda_);

    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    const int uprestart = fsidyn.get<int>("RESTARTEVRY");
    const int upres = fsidyn.get<int>("RESULTSEVRY");
    if((uprestart != 0 && FluidField()->Step() % uprestart == 0) || FluidField()->Step() % upres == 0)
      StructureField()->DiscWriter()->WriteVector("fsilambda", lambdafull);
  }

  FluidField()->    Output();
  AleField()->      Output();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if(Comm().MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::ReadRestart(int step)
{
  // read Lagrange multiplier
  {
    Teuchos::RCP<Epetra_Vector> lambdafull = Teuchos::rcp(new Epetra_Vector(*StructureField()->DofRowMap(),true));
    IO::DiscretizationReader reader = IO::DiscretizationReader(StructureField()->Discretization(),step);
    reader.ReadVector(lambdafull, "fsilambda");
    lambda_ = StructureField()->Interface()->ExtractFSICondVector(lambdafull);
  }

  StructureField()->ReadRestart(step);
  FluidField()->ReadRestart(step);
  AleField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(),FluidField()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::CreateCombinedDofRowMap()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicStructureSplitNoNOX::SetupNewSystem()");

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;
  vecSpaces.push_back(StructureField()->Interface()->OtherMap());
  vecSpaces.push_back(FluidField()->DofRowMap());
  vecSpaces.push_back(AleField()->Interface()->OtherMap());

  if (vecSpaces[0]->NumGlobalElements()==0)
    dserror("No inner structural equations. Splitting not possible. Panic.");

  SetDofRowMaps(vecSpaces);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidFluidMonolithicStructureSplitNoNOX::BuildConvergenceNorms()
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicStructureSplitNoNOX::BuildCovergenceNorms()");
  //----------------------------
  // build residual norms
  rhs_->Norm2(&normrhs_);

  // structural Dofs
  StructureField()->Interface()->ExtractOtherVector(StructureField()->RHS())->Norm2(&normstrrhsL2_);
  StructureField()->Interface()->ExtractOtherVector(StructureField()->RHS())->NormInf(&normstrrhsInf_);

  // extract fluid Dofs
  Teuchos::RCP<const Epetra_Vector> rhs = Extractor().ExtractVector(rhs_,1);

  FluidField()->Interface()->ExtractFSICondVector(rhs)->Norm2(&norminterfacerhsL2_);
  FluidField()->Interface()->ExtractFSICondVector(rhs)->NormInf(&norminterfacerhsInf_);

  // inner fluid velocity Dofs: inner velocity dofs without db-condition
  std::vector<Teuchos::RCP<const Epetra_Map> > innerfluidvel;
  innerfluidvel.push_back(FluidField()->InnerVelocityRowMap());
  innerfluidvel.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidvelextract(*(FluidField()->DofRowMap()),innerfluidvel);
  fluidvelextract.ExtractVector(FluidField()->RHS(),0)->Norm2(&normflvelrhsL2_);
  fluidvelextract.ExtractVector(FluidField()->RHS(),0)->NormInf(&normflvelrhsInf_);

  // fluid pressure Dofs: pressure dofs (at interface and inner) without db-condition
  std::vector<Teuchos::RCP<const Epetra_Map> > fluidpres;
  fluidpres.push_back(FluidField()->PressureRowMap());
  fluidpres.push_back(Teuchos::null);
  LINALG::MultiMapExtractor fluidpresextract(*(FluidField()->DofRowMap()),fluidpres);
  fluidpresextract.ExtractVector(FluidField()->RHS(),0)->Norm2(&normflpresrhsL2_);
  fluidpresextract.ExtractVector(FluidField()->RHS(),0)->NormInf(&normflpresrhsInf_);

  // ale
  AleField()->RHS()->Norm2(&normalerhsL2_);
  //-------------------------------
  // build solution increment norms

  // build increment norm
  iterinc_->Norm2(&norminc_);

  // structural Dofs
  Extractor().ExtractVector(iterinc_,0)->Norm2(&normstrincL2_);
  Extractor().ExtractVector(iterinc_,0)->NormInf(&normstrincInf_);

  // interface
  Teuchos::RCP<const Epetra_Vector> inc = Extractor().ExtractVector(iterinc_,1);
  FluidField()->Interface()->ExtractFSICondVector(inc)->Norm2(&norminterfaceincL2_);
  FluidField()->Interface()->ExtractFSICondVector(inc)->NormInf(&norminterfaceincInf_);

  // inner fluid velocity Dofs
  fluidvelextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),0)->Norm2(&normflvelincL2_);
  fluidvelextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),0)->NormInf(&normflvelincInf_);

  // fluid pressure Dofs
  fluidpresextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),0)->Norm2(&normflpresincL2_);
  fluidpresextract.ExtractVector(Extractor().ExtractVector(iterinc_,1),0)->NormInf(&normflpresincInf_);

  // ale
  Extractor().ExtractVector(iterinc_,2)->Norm2(&normaleincL2_);

  //get length of the structural, fluid and ale vector
  ns_ = (*(StructureField()->RHS())).GlobalLength(); //structure
  ni_ = (*(FluidField()->Interface()->ExtractFSICondVector(rhs))).GlobalLength(); //fluid fsi
  nf_ = (*(FluidField()->RHS())).GlobalLength(); //fluid inner
  nfv_ = (*(fluidvelextract.ExtractVector(FluidField()->RHS(),0))).GlobalLength(); //fluid velocity
  nfp_ = (*(fluidpresextract.ExtractVector(FluidField()->RHS(),0))).GlobalLength();//fluid pressure
  na_ = (*(AleField()->RHS())).GlobalLength(); //ale
  nall_ = (*rhs_).GlobalLength(); //all

}

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface                     */
/*----------------------------------------------------------------------*/
 void FSI::FluidFluidMonolithicStructureSplitNoNOX::RecoverLagrangeMultiplier()
 {
   TEUCHOS_FUNC_TIME_MONITOR("FSI::FluidFluidMonolithicStructureSplitNoNOX::RecoverLagrangeMultiplier");

   // get time integration parameters of structural time integrator
   // to enable consistent time integration among the fields
   const double stiparam = StructureField()->TimIntParam();

  // some often re-used vectors
  Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;     // stores intermediate result of terms (3)-(8)
  Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
  Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

  /* Recovery of Lagrange multiplier lambda^{n+1} is done by the following
   * condensation expression:
   *
   * lambda^{n+1} =
   *
   * (1)  - stiparam / (1.-stiparam) * lambda^{n}
   *
   * (2)  + 1. / (1.-stiparam) * tmpvec
   *
   * with tmpvec =
   *
   * (3)    r_{\Gamma}^{S,n+1}
   *
   * (4)  + S_{\Gamma I} * \Delta d_{I}^{S,n+1}
   *
   * (5)  + tau * S_{\Gamma\Gamma} * \Delta u_{\Gamma}^{F,n+1}
   *
   * (6)  + dt * S_{\Gamma\Gamma} * u_{\Gamma}^n]
   *
   * Remark on term (6):
   * Term (6) has to be considered only in the first Newton iteration.
   * Hence, it will usually not be computed since in general we need more
   * than one nonlinear iteration until convergence.
   *
   * Remarks on all terms:
   * +  Division by (1.0 - stiparam) will be done in the end
   *    since this is common to all terms
   * +  tau: time scaling factor for interface time integration (tau = 1/FluidField()->TimeScaling())
   * +  neglecting terms (4)-(5) should not alter the results significantly
   *    since at the end of the time step the solution increments tend to zero.
   *
   *                                                 Matthias Mayr (10/2012)
   */

  // ---------Addressing term (1)
  lambda_->Update(-stiparam,*lambda_,0.0);
  // ---------End of term (1)

  // ---------Addressing term (3)
  Teuchos::RCP<Epetra_Vector> structureresidual = StructureField()->Interface()->ExtractFSICondVector(StructureField()->RHS());
  structureresidual->Scale(-1.0); // invert sign to obtain residual, not rhs
  tmpvec = Teuchos::rcp(new Epetra_Vector(*structureresidual));
  // ---------End of term (3)

  /* You might want to comment out terms (4) to (6) since they tend to
   * introduce oscillations in the Lagrange multiplier field for certain
   * material properties of the structure.
   *                                                    Matthias Mayr 11/2012
  // ---------Addressing term (4)
  auxvec = Teuchos::rcp(new Epetra_Vector(sgicur_->RangeMap(),true));
  sgicur_->Apply(*ddiinc_,*auxvec);
  tmpvec->Update(1.0,*auxvec,1.0);
  // ---------End of term (4)

  // ---------Addressing term (5)
  auxvec = Teuchos::rcp(new Epetra_Vector(sggcur_->RangeMap(),true));
  sggcur_->Apply(*ddginc_,*auxvec);
  tmpvec->Update(1.0/timescale,*auxvec,1.0);
  // ---------End of term (5)

  //---------Addressing term (6)
  if (firstcall_)
  {
    auxvec = Teuchos::rcp(new Epetra_Vector(sggprev_->RangeMap(),true));
    sggprev_->Apply(*FluidToStruct(FluidField()->ExtractInterfaceVeln()),*auxvec);
    tmpvec->Update(Dt(),*auxvec,1.0);
  }
  // ---------End of term (6)
  *
  */

  // ---------Addressing term (2)
  lambda_->Update(1.0,*tmpvec,1.0);
  // ---------End of term (2)

  // finally, divide by -(1.-stiparam) which is common to all terms
  lambda_->Scale(1./(1.0-stiparam));

  // Finally, the Lagrange multiplier 'lambda_' is recovered here.
  // It represents nodal forces acting onto the structure.


//------ old version changed on 5/12/12  --------------
//   // compute the product S_{\Gamma I} \Delta d_I
//   Teuchos::RCP<Epetra_Vector> sgiddi = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true); // store the prodcut 'S_{\GammaI} \Delta d_I^{n+1}' in here
//   (sgicur_->EpetraMatrix())->Multiply(false, *ddiinc_, *sgiddi);

//    // compute the product S_{\Gamma\Gamma} \Delta d_\Gamma
//    Teuchos::RCP<Epetra_Vector> sggddg = LINALG::CreateVector(*StructureField()->Interface()->OtherMap(),true); // store the prodcut '\Delta t / 2 * S_{\Gamma\Gamma} \Delta u_\Gamma^{n+1}' in here
//    (sggcur_->EpetraMatrix())->Multiply(false, *ddginc_, *sggddg);

//    // Update the Lagrange multiplier:
//    /* \lambda^{n+1} =  - a/b*\lambda^n - f_\Gamma^S
//     *                  - S_{\Gamma I} \Delta d_I - S_{\Gamma\Gamma} \Delta d_\Gamma
//     */
//    lambda_->Update(1.0, *fgcur_, -stiparam);
//    //lambda_->Update(-1.0, *sgiddi, -1.0, *sggddg, 1.0);
//    lambda_->Scale(1/(1.0-stiparam)); //entire Lagrange multiplier it divided by (1.-stiparam)
//

   return;
}

void FSI::FluidFluidMonolithicStructureSplitNoNOX::HandleFluidDofMapChangeInNewton()
{
  if (Comm().MyPID() == 0)
    IO::cout << " New Map!! " << IO::endl;

  // save the old x_sum
  Teuchos::RCP<Epetra_Vector> x_sum_n =  LINALG::CreateVector(*DofRowMap(), true);
  *x_sum_n = *x_sum_;
  Teuchos::RCP<const Epetra_Vector> sx_n;
  Teuchos::RCP<const Epetra_Vector> ax_n;
  sx_n = Extractor().ExtractVector(x_sum_n,0);
  ax_n = Extractor().ExtractVector(x_sum_n,2);

  // rebuild combined dof-map
  CreateCombinedDofRowMap();
  // re-initialize systemmatrix_
  systemmatrix_ = Teuchos::rcp(
    new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
      Extractor(),
      Extractor(),
      81,
      false,
      true
      )
    );

  rhs_ = LINALG::CreateVector(*DofRowMap(), true);
  iterinc_ = LINALG::CreateVector(*DofRowMap(), true);
  zeros_ = LINALG::CreateVector(*DofRowMap(), true);
  x_sum_ = LINALG::CreateVector(*DofRowMap(),true);

  // build the new iter_sum
  Extractor().InsertVector(sx_n,0,x_sum_);
  Extractor().InsertVector(FluidField()->Stepinc(),1,x_sum_);
  Extractor().InsertVector(ax_n,2,x_sum_);
  nf_ = (*(FluidField()->RHS())).GlobalLength();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::FluidFluidMonolithicStructureSplitNoNOX::HasFluidDofMapChanged(
        const Epetra_BlockMap& fluidincrementmap)
{
  return ! FluidField()->DofRowMap()->SameAs(fluidincrementmap);
}
