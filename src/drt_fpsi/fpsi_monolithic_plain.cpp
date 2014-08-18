/*----------------------------------------------------------------------*/
/*!
\file fpsi_monolithic_plain.cpp

\brief Solve FPSI problem with matching grids using a monolithic scheme
       in its plain form. Only interface ale displacements are condensed.

<pre>
Maintainer: Andreas Rauch
            rauch@lnm.mw.tum.de

</pre>
*/

/*----------------------------------------------------------------------*/
//GENERAL
#include <Teuchos_TimeMonitor.hpp>

//FPSI includes
#include "fpsi_monolithic_plain.H"
#include "fpsi_monolithic.H"

//POROELAST includes
#include "../drt_poroelast/poro_base.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_poroelast/poroelast_monolithic.H"

//FSI includes
#include "../drt_fsi/fsi_matrixtransform.H"
#include "../drt_fsi/fsi_debugwriter.H"
#include "../drt_fsi/fsi_statustest.H"
#include "../drt_fsi/fsi_overlapprec_fsiamg.H"
#include "../drt_fsi/fsi_monolithic_linearsystem.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"

//FLUID includes
#include "../drt_fluid/fluid_utils_mapextractor.H"

//STRUCTURE includes
#include "../drt_structure/stru_aux.H"

//ALE includes
#include "../drt_ale/ale_utils_mapextractor.H"

//LINALG includes
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"

//IO includes
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"

//OTHERS
#include "../drt_constraint/constraint_manager.H"
#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::Monolithic_Plain::Monolithic_Plain(const Epetra_Comm& comm,
                                         const Teuchos::ParameterList& fpsidynparams,
                                         const Teuchos::ParameterList& poroelastdynparams)
  : Monolithic(comm,fpsidynparams,poroelastdynparams)
{
  // create transformation objects for coupling terms
  couplingrowtransform_     = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform2_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform3_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform4_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform5_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingcoltransform_     = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  couplingcoltransform2_    = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  couplingrowcoltransform_  = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  couplingrowcoltransform2_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);


  // create transformation object for the ale condensation
  aigtransform_           = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  kpatransform_           = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);


  fggtransform_           = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  fggtransform2_           = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  fgitransform_           = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  fgitransform2_           = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  figtransform_           = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  figtransform2_           = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  aigtransform2_          = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  fmiitransform_          = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  fmgitransform_          = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);

  // Recovery of Lagrange multiplier happens on fluid field
  lambda_   = Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FSICondMap(),true));
  fmgiprev_ = Teuchos::null;
  fmgicur_  = Teuchos::null;
  fmggprev_ = Teuchos::null;
  fmggcur_  = Teuchos::null;
  fgiprev_  = Teuchos::null;
  fgicur_   = Teuchos::null;
  fggprev_  = Teuchos::null;
  fggcur_   = Teuchos::null;

#ifdef DEBUG
  // check whether allocation was successful
  if (fggtransform_   == Teuchos::null) { dserror("Allocation of 'fggtransform_' failed."); }
  if (fggtransform2_   == Teuchos::null) { dserror("Allocation of 'fggtransform2_' failed."); }
  if (fgitransform_   == Teuchos::null) { dserror("Allocation of 'fgitransform_' failed."); }
  if (fgitransform2_   == Teuchos::null) { dserror("Allocation of 'fgitransform2_' failed."); }
  if (figtransform_   == Teuchos::null) { dserror("Allocation of 'figtransform_' failed."); }
  if (figtransform2_   == Teuchos::null) { dserror("Allocation of 'figtransform2_' failed."); }
  if (aigtransform_   == Teuchos::null) { dserror("Allocation of 'aigtransform_' failed."); }
  if (fmiitransform_  == Teuchos::null) { dserror("Allocation of 'fmiitransform_' failed."); }
  if (fmgitransform_  == Teuchos::null) { dserror("Allocation of 'fmgitransform_' failed."); }
  if (lambda_         == Teuchos::null) { dserror("Allocation of 'lambda_' failed."); }
#endif

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupSystem()
{

  const Teuchos::ParameterList& fpsidynparams   = DRT::Problem::Instance()->FPSIDynamicParams();

  SetDefaultParameters(fpsidynparams);

  // call SetupSystem in base classes
  PoroField()->SetupSystem();
  FPSI::Monolithic::SetupSystem();

  // create combined map
  std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces;

  vecSpaces.push_back(PoroField()->DofRowMap());

  vecSpaces.push_back(FluidField()->DofRowMap());

  vecSpaces.push_back(AleField()      ->Interface()->OtherMap());

  //Modify block_numbers manually according to vecSpaces.push_back order!
  poro_block_               = 0;
  fluid_block_              = 1;
  ale_i_block_              = 2;

  if(vecSpaces[poro_block_]->NumGlobalElements()==0)
    dserror("Poro: Hey I just met you and this is crazy, but here's my number "
            "so call me maybe ... !\n (0) 89 289 15300");
  if(vecSpaces[fluid_block_]->NumGlobalElements()==0)
    dserror("No inner fluid equations. Splitting not possible. Roundhouse-Kick!");
  if(vecSpaces[ale_i_block_]->NumGlobalElements()==0)
    dserror("ALE ?! Roundhouse-Kick!");
  // merge maps and create full monolithic FPSI-DofRowMap
  SetDofRowMaps(vecSpaces);

  // switch fluid to interface split block matrix
  FluidField()->UseBlockMatrix(true);

  // build ale system matrix in splitted system
  AleField()->BuildSystemMatrix(false);

  // initialize FPSI-systemmatrix_
  systemmatrix_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<
      LINALG::DefaultBlockMatrixStrategy>(Extractor(), Extractor(), 81, false,
      true));

  // create off-diagonal coupling matrices

  // Block Porous medium - Fluid
  k_pf_ = Teuchos::rcp(new LINALG::SparseMatrix(
                  *(PoroField()->DofRowMap()), 81, true, true));

  // Block Fluid - Porous medium

  k_fp_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(Extractor().Map(fluid_block_)), 81, true, true)); //Changed it -- is done!! ChrAg

  k_pf_struct_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(PoroField()->StructureField()->DofRowMap()), 81, true, true));

  k_pf_porofluid_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(FluidField()->DofRowMap()),  81, true, true));//Should also be changed? ChrAg

  k_fp_porofluid_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(PoroField()->FluidField()->Discretization()->DofRowMap()),  81, true, true));

  // Sub Block PoroFluid - Structure
  k_pfs_ = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(PoroField()->FluidField()->Discretization()->DofRowMap()),
                          81, true, true));

  // create matrices needed for ale condensation

  k_pa_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<
              LINALG::DefaultBlockMatrixStrategy>(*AleField()->Interface(),*FluidField()->FPSIInterface(), 81, false, true));

  k_ap_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(AleField()->Interface()->OtherMap()), 81, true, true));
  k_aa_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(AleField()->Interface()->OtherMap()), 81, true, true));

  k_fa_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*AleField()->Interface(),*FluidField()->FPSIInterface(), 81, false, true));

  // create additional map extractors!

  //   create combined map
   std::vector<Teuchos::RCP<const Epetra_Map> > vecSpaces2;

  Teuchos::RCP<Epetra_Map> Fluidcombined = LINALG::MergeMap(FluidField()->Interface()->OtherMap(),FluidField()->Interface()->FPSICondMap());
  vecSpaces2.push_back(Fluidcombined);
  vecSpaces2.push_back(FluidField()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi));

  SetFluidFPSI_FSIDofRowMaps(vecSpaces2);


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetDofRowMaps(const std::vector<Teuchos::RCP<
    const Epetra_Map> >& maps)
{
  fullmap_ = LINALG::MultiMapExtractor::MergeMaps(maps);

  // full FPSI-blockmap
  blockrowdofmap_.Setup(*fullmap_, maps);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetFluidFPSI_FSIDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  Teuchos::RCP<Epetra_Map> fullfluidfpsifsimap = LINALG::MultiMapExtractor::MergeMaps(maps);

  // fullfluidfpsifsimap
  fluidfpsi_fsi_rowdofmap_.Setup(*fullfluidfpsifsimap, maps);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::SetDefaultParameters(const Teuchos::ParameterList& fpsidynparams)
{
// to do
return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupRHS(bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic_Plain::SetupRHS");
  // create full monolithic rhs vector

  rhs_ = Teuchos::rcp(new Epetra_Vector(*DofRowMap(), true));

  firstcall_ = firstcall;

  PoroField()->SetupRHS(firstcall_);

  SetupVector(*rhs_,
              PoroField() ->RHS(),
              FluidField()->RHS(),
              AleField()  ->RHS(),
              FluidField()->ResidualScaling());

  if (FSI_Interface_exists_)
  {
    SetupRHSLambda(*rhs_);
    if (firstcall)
    {
      SetupRHSFirstIter(*rhs_);
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic_Plain::SetupSystemMatrix");
    mat.UnComplete();

    // get single field block matrices
    Teuchos::RCP<LINALG::BlockSparseMatrixBase>       pbm    = PoroField() -> SystemBlockMatrix();

    // Todo: use only blockmatrix pbm instead of sparse matrix
    // vuong 08/14
    Teuchos::RCP<LINALG::SparseMatrix>                p      = PoroField() -> SystemSparseMatrix();

    const Teuchos::RCP<LINALG::SparseMatrix>          f      = FluidField()-> SystemSparseMatrix();
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase> fbm    = FluidField()-> BlockSystemMatrix();
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a      = AleField()  -> BlockSystemMatrix();

    //Get Idx of fluid and ale field map extractors
    const int &fidx_other = FLD::UTILS::MapExtractor::cond_other;
    const int &fidx_fsi = FLD::UTILS::MapExtractor::cond_fsi;
    const int &fidx_fpsi = FLD::UTILS::MapExtractor::cond_fpsi;

    const int &aidx_other = ALE::UTILS::MapExtractor::cond_other;
    const int &aidx_fsi = ALE::UTILS::MapExtractor::cond_fsi;
    const int &aidx_fpsi = ALE::UTILS::MapExtractor::cond_fpsi;

    // FPSI Couplings
    const ADAPTER::Coupling& coupsa_fpsi        =   StructureAleCoupling_FPSI();
    const ADAPTER::Coupling& coupsf_fpsi        =   StructureFluidCoupling_FPSI();
    const ADAPTER::Coupling& smallcoupsf_fpsi   =   SmallStructureFluidCoupling_FPSI();

    //General Couplings
    const ADAPTER::Coupling& coupfa             =   FluidAleCoupling();

    //FSI Couplings
    const ADAPTER::Coupling& coupsf_fsi         =   StructureFluidCoupling_FSI();
    const ADAPTER::Coupling& coupsa_fsi         =   StructureAleCoupling_FSI();

    // extract ale submatrices
    Teuchos::RCP<LINALG::SparseMatrix> k_pa = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pa_);
    Teuchos::RCP<LINALG::SparseMatrix> k_ap = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_ap_);
    Teuchos::RCP<LINALG::SparseMatrix> k_aa = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_aa_);
    k_ap -> UnComplete();
    k_aa -> UnComplete();

    // scaling factors for fluid
    const double scale     = FluidField()->ResidualScaling();
    const double timescale = FluidField()->TimeScaling();

    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = PoroField()->StructureField()->TimIntParam();
    const double ftiparam = FluidField()->TimIntParam();

  /*----------------------------------------------------------------------*/
  // build block matrix

     /*----------------------------------------------------------------------*/
     //insert poro
     p->UnComplete();
     mat.Assign(poro_block_,poro_block_,View, *p);

     //Assign fii + Coupling Parts
     f->UnComplete();
     mat.Assign(fluid_block_,fluid_block_,View,*f);

  //Insert condensed Fluid Blocks: Fgg and Fgi (+ Fg_gFPSI) --> g is on the FSI-Interface

     if (FSI_Interface_exists_)
     {
         // extract fluid submatrices -- use block matrices just for fsi boundary matrices as they have to be condensed!
         // --> others will be assigned directly by the fluid sparse matrix
         LINALG::SparseMatrix fgi      =   fbm->Matrix(fidx_fsi,fidx_other);
         LINALG::SparseMatrix fg_gfpsi =   fbm->Matrix(fidx_fsi,fidx_fpsi);
         LINALG::SparseMatrix& fgg     =   fbm->Matrix(fidx_fsi,fidx_fsi);

        (*fggtransform_)( fgg,
                         (1.0-stiparam)/(1.0-ftiparam)*scale*timescale*1,
                         ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                         ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                         *p,
                         true,
                         true);

        (*fgitransform_)(fgi,
                        (1.0-stiparam)/(1.0-ftiparam)*scale,
                        ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                        mat.Matrix(poro_block_, fluid_block_),
                        false); //Addmatrix

        (*fgitransform2_)(fg_gfpsi,
                        (1.0-stiparam)/(1.0-ftiparam)*scale,
                        ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                        mat.Matrix(poro_block_, fluid_block_),
                        true); //Addmatrix

    // Insert ale: Aii, Aig and Ai_gFPSI--> g is on the FSI-Interface

       LINALG::SparseMatrix& aig = a->Matrix(aidx_other,aidx_fsi);

      (*aigtransform2_)( a->FullRowMap(),
                       a->FullColMap(),
                       aig,
                       1.,
                       ADAPTER::CouplingSlaveConverter(coupsa_fsi),
                       mat.Matrix(ale_i_block_,poro_block_));
    }

     LINALG::SparseMatrix& aii = a->Matrix(aidx_other,aidx_other);
     LINALG::SparseMatrix& ai_gfpsi = a->Matrix(aidx_other,aidx_fpsi);

    Teuchos::RCP<LINALG::SparseMatrix>  laig = Teuchos::rcp(new LINALG::SparseMatrix(ai_gfpsi.RangeMap(),81,false));
    (*aigtransform_)( a->FullRowMap(),
                    a->FullColMap(),
                    ai_gfpsi,
                    1.,
                    ADAPTER::CouplingSlaveConverter(coupsa_fpsi),
                    *laig);

    laig->Complete(p->DomainMap(),laig->RangeMap());

    k_aa -> Add(  aii,false,1.0,0.0);
    k_ap -> Add(*laig,false,1.0,0.0);

    k_aa -> Complete(aii.DomainMap(),aii.RangeMap());
    k_ap -> Complete(p->DomainMap(),laig->RangeMap());
    mat.Assign(ale_i_block_,ale_i_block_,View,*k_aa);
    mat.Assign(ale_i_block_,poro_block_,View,*k_ap);


  // insert coupling terms

  Teuchos::RCP<LINALG::SparseMatrix> k_pf           = PoroFluidCouplingMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> k_fp           = FluidPoroCouplingMatrix();

  ApplyCouplingTerms(k_pf,k_fp,p,f,k_pa_, k_fa_);

  Teuchos::RCP<LINALG::SparseMatrix> aa  = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()),81,false));
  aa->Assign(View,Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_pa_)->Matrix(fidx_fpsi,aidx_other));
  aa->Complete(*AleField()->Interface()->OtherMap(),*FluidField()->FPSIInterface()->FPSICondMap());

  Teuchos::RCP<LINALG::SparseMatrix> aaa = Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->DofRowMap()),81,false));
  (*couplingrowtransform5_)( *aa,
                              1.0,
                              ADAPTER::CouplingSlaveConverter(coupsf_fpsi), // important to use slave converter
                             *aaa,
                              false);

  aaa -> Complete(*AleField()->Interface()->OtherMap(),*(PoroField()->DofRowMap()));

  mat.Matrix(poro_block_, ale_i_block_).Add(*aaa,false,1.0,1.0);

  Teuchos::RCP<LINALG::SparseMatrix> aa2 = Teuchos::rcp(new LINALG::SparseMatrix(*(Extractor().Map(fluid_block_)),81,false));
  aa2 -> Zero();
  Teuchos::RCP<LINALG::SparseMatrix> fluid_ale_gamma = Teuchos::rcp(new LINALG::SparseMatrix(*(Extractor().Map(fluid_block_)),81,false));
  Teuchos::RCP<LINALG::SparseMatrix> fluid_ale_omega = Teuchos::rcp(new LINALG::SparseMatrix(*(Extractor().Map(fluid_block_)),81,false));
  fluid_ale_omega -> Assign(View,Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_fa_)->Matrix(fidx_other,aidx_other));
  fluid_ale_omega -> Complete(*AleField()->Interface()->OtherMap(),*(Extractor().Map(fluid_block_)));
  fluid_ale_gamma -> Assign(View,Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_fa_)->Matrix(fidx_fpsi,aidx_other));
  fluid_ale_gamma -> Complete(*AleField()->Interface()->OtherMap(),*(Extractor().Map(fluid_block_)));  //map of fluid_block_shouldnt be right...?

  aa2 -> Add(*fluid_ale_omega,false,1.0,1.0);
  aa2 -> Add(*fluid_ale_gamma,false,1.0,1.0);

  aa2 -> Complete(AleField()->BlockSystemMatrix()->DomainMap(aidx_other),*(Extractor().Map(fluid_block_)));
  //mat.Assign(fluid_block_,ale_i_block_,View,*aa2); // assign fluid_ale coupling matrix (due to Neumann Integration)
  mat.Matrix(fluid_block_,ale_i_block_).Add(*aa2,false,1.0,1.0);


//  //////////////////////////////////////////////
//  //////                                  //////
//  //////    Linearization of FluidField   //////
//  //////    with respect to ale mesh      //////
//  //////             motion               //////
//  //////                                  //////
//  //////////////////////////////////////////////
  const Teuchos::ParameterList& fpsidynparams   = DRT::Problem::Instance()->FPSIDynamicParams();
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> fluidalematrix = FluidField()->ShapeDerivatives();
  if(Teuchos::getIntegralValue<int>(fpsidynparams,"USESHAPEDERIVATIVES"))
  {
    k_fp_->UnComplete();
//
    // extract submatrices
//
//
    if (fluidalematrix != Teuchos::null)
    {
      //There is no fpsi-fsi overlap in the block matrixes, all dofs which are on both interfaces belong
      //to the fsi-block matrix!

            LINALG::SparseMatrix& fluidalematrix_ii = fluidalematrix->Matrix(fidx_other,fidx_other);

            LINALG::SparseMatrix& fluidalematrix_gg = fluidalematrix->Matrix(fidx_fpsi,fidx_fpsi);
            LINALG::SparseMatrix& fluidalematrix_gi = fluidalematrix->Matrix(fidx_fpsi,fidx_other);
            LINALG::SparseMatrix& fluidalematrix_ig = fluidalematrix->Matrix(fidx_other,fidx_fpsi);
//
//      // add fluid_ale block ii and gi
//      // those two blocks are not condensed since they belong to the columns of the inner ale dofs
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_ii,
          1.0,
          ADAPTER::CouplingMasterConverter(coupfa), // row converter: important to use slave converter
          mat.Matrix(fluid_block_,ale_i_block_),
          false, // bool exactmatch = true (default)
          true);

      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_gi,
          1.0,
          ADAPTER::CouplingMasterConverter(coupfa), // row converter: important to use slave converter
          mat.Matrix(fluid_block_,ale_i_block_),
          false, // bool exactmatch = true (default)
          true);
//
      // reuse k_fp to assign linearization with respect to condensed interface dofs to the structural block
      Teuchos::RCP<LINALG::SparseMatrix> tempale1 = Teuchos::rcp(new LINALG::SparseMatrix(fluidalematrix_ig.RowMap(),81,false));
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_ig,
          1.0,
          ADAPTER::CouplingMasterConverter(smallcoupsf_fpsi), // row converter: important to use slave converter
          *tempale1,
          false,  // bool exactmatch = true (default)
          false); // bool add

      Teuchos::RCP<LINALG::SparseMatrix> tempale2 = Teuchos::rcp(new LINALG::SparseMatrix(fluidalematrix_gg.RowMap(),81,false));
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_gg,
          1.0,
          ADAPTER::CouplingMasterConverter(smallcoupsf_fpsi), // row converter: important to use slave converter
          *tempale2,
          false,  // bool exactmatch = true (default)
          false); // bool add

      tempale1  -> Complete(PoroField()->StructureRangeMap(),f->RangeMap());
      tempale2  -> Complete(PoroField()->StructureRangeMap(),f->RangeMap());

      k_fp -> Add(*tempale1,false,1.0,1.0);
      k_fp -> Add(*tempale2,false,1.0,1.0);
      k_fp -> Complete(p->RangeMap(),f->RangeMap());

      if (FSI_Interface_exists_)
      {
        LINALG::SparseMatrix& fluidalematrix_gg_fsi = fluidalematrix->Matrix(fidx_fsi,fidx_fsi);
        LINALG::SparseMatrix& fluidalematrix_gi_fsi = fluidalematrix->Matrix(fidx_fsi,fidx_other);
        LINALG::SparseMatrix& fluidalematrix_ig_fsi = fluidalematrix->Matrix(fidx_other,fidx_fsi);

        LINALG::SparseMatrix& fluidalematrix_gfsigfpsi = fluidalematrix->Matrix(fidx_fsi,fidx_fpsi);
        LINALG::SparseMatrix& fluidalematrix_gfpsigfsi = fluidalematrix->Matrix(fidx_fpsi,fidx_fsi);


        (*figtransform_)(FluidField()->BlockSystemMatrix()->FullRowMap(),
                  FluidField()->BlockSystemMatrix()->FullColMap(),
                 fluidalematrix_ig_fsi,
                 1.0,
                 ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                 mat.Matrix(fluid_block_,poro_block_),
                 false, // bool exactmatch = true (default)
                 true);

        (*figtransform2_)(FluidField()->BlockSystemMatrix()->FullRowMap(),
                  FluidField()->BlockSystemMatrix()->FullColMap(),
                 fluidalematrix_gfpsigfsi,
                 1.0,
                 ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                 mat.Matrix(fluid_block_,poro_block_),
                 false,
                 true);

        (*fggtransform_)(fluidalematrix_gg_fsi,
                 (1.0-stiparam)/(1.0-ftiparam)*scale,
                 ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                 ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                 *p,
                 false,
                 true);

        (*fggtransform2_)(fluidalematrix_gfsigfpsi,
                 (1.0-stiparam)/(1.0-ftiparam)*scale,
                 ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                 ADAPTER::CouplingMasterConverter(smallcoupsf_fpsi),
                 *p,
                 false,
                 true);

        (*fmgitransform_)(fluidalematrix_gi_fsi,
                (1.0-stiparam)/(1.0-ftiparam)*scale,
                ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                ADAPTER::CouplingMasterConverter(coupfa),
                mat.Matrix(poro_block_,ale_i_block_),
                false,
                true);
      }
    }
    else // if shapederivatives = no in FluidDynamics section in dat-file
    {
      std::cout<<"WARNING: Linearization with respect to mesh motion of fluid subproblem is switched off!"<<std::endl;
    }
  } // if useshapederivatives


      //////////////////////////////////////////////////////
     ///////                                        ///////
    ///////   FILL FSI/FPSI - INTERFACE OVERLAP    ///////
   ///////                                        ///////
  //////////////////////////////////////////////////////

  if (FSI_Interface_exists_)
  {
    (*figtransform_)(p->RowMap(),
                     f->ColMap(),
                     *k_pf,
                     timescale,
                     ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                     mat.Matrix(poro_block_,poro_block_),
                     false,//no exactmatch! (just FSI Part should be extracted)
                     true); //Add


    //Add (Tau * Fig) to Structural Column (condensation of fluid velocities) ... done here to catch also FPSI-Coupling terms (for overlapping FSI/FPSI Interfaces)
    (*figtransform_)(f->RowMap(),
                     f->ColMap(),
                     *f,
                     timescale,
                     ADAPTER::CouplingSlaveConverter(coupsf_fsi),
                     mat.Matrix(fluid_block_,poro_block_),
                     false, //no exactmatch! (just FSI Part should be extracted)
                     false); //Assign

    mat.Matrix(poro_block_,fluid_block_).Add(*k_pf,false,1.0,1.0);// add poro_fluid coupling matrix
    mat.Matrix(fluid_block_,poro_block_).Add(*k_fp,false,1.0,1.0);// add fluid_poro coupling matrix
  }
  else //Add blocks as already matrices have been added to these blocks!
  {
    mat.Assign(poro_block_,fluid_block_, View,*k_pf);
    mat.Assign(fluid_block_,poro_block_, View,*k_fp);
  }

//    //+++ part in the f...Matrix from overlap will be removed by ApplyDBC of condensed velocity DOFS in LinearSolve()!

  // done. make sure all blocks are filled.

  mat.Complete();

  //if (FSI_Interface_exists_)
  {
    // store parts of fluid matrix to know them in the next iteration as previous iteration matrices
    fgiprev_ = fgicur_;
    fggprev_ = fggcur_;
    fgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(fbm->Matrix(FLD::UTILS::MapExtractor::cond_fsi,FLD::UTILS::MapExtractor::cond_other)));
    fggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(fbm->Matrix(FLD::UTILS::MapExtractor::cond_fsi,FLD::UTILS::MapExtractor::cond_fsi)));

    // store parts of fluid shape derivative matrix to know them in the next iteration as previous iteration matrices
    fmgiprev_ = fmgicur_;
    fmggprev_ = fmggcur_;
    if (fluidalematrix!=Teuchos::null)
    {
      fmgicur_ = Teuchos::rcp(new LINALG::SparseMatrix(fluidalematrix->Matrix(FLD::UTILS::MapExtractor::cond_fsi,FLD::UTILS::MapExtractor::cond_other)));
      fmggcur_ = Teuchos::rcp(new LINALG::SparseMatrix(fluidalematrix->Matrix(FLD::UTILS::MapExtractor::cond_fsi,FLD::UTILS::MapExtractor::cond_fsi)));
    }
  }

//  // Write Matrix to Matlab format
//  static int globindex = 0;
//  ++globindex;
//
//  // print to file in matlab format
//  std::ostringstream filename;
//  std::string filebase = "blockmatrix";
//  filename << "/scratch/rauch/workspace/baci/baci_debug/matlab_output/" << filebase << "_" << globindex << ".mtl";
//  //std::cout<<"filename: "<< filename.str() <<endl;
//  LINALG::PrintBlockMatrixInMatlabFormat(filename.str().c_str(),mat);

  //std::cout<<" !!! Matrices Built !!! "<<endl;
  //dserror(" !!!  Matrices Built !!! ");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupVector(Epetra_Vector &f,
                                         Teuchos::RCP<const Epetra_Vector> sv,
                                         Teuchos::RCP<const Epetra_Vector> fv,
                                         Teuchos::RCP<const Epetra_Vector> av,
                                         double fluidscale)
{

  //Get FluidField Block Matrix
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> fbm    = FluidField()-> BlockSystemMatrix();
  const LINALG::SparseMatrix fgg = fbm->Matrix(FLD::UTILS::MapExtractor::cond_fsi,FLD::UTILS::MapExtractor::cond_fsi);

  // get time integration parameters of structure and fluid time integrators
  // to enable consistent time integration among the fields
  const double stiparam = PoroField()->StructureField()->TimIntParam();
  const double ftiparam = FluidField()->TimIntParam();

  Extractor().InsertVector(*fv,fluid_block_,f);  // add fluid contributions to 'f'

  if (FSI_Interface_exists_) //in case FSI interface exists, add term from condensation to RHS
  {
  // add fluid interface values to structure vector
  Teuchos::RCP<Epetra_Vector> fcvgfsi = FluidField()->Interface()->ExtractFSICondVector(fv);

  Teuchos::RCP<Epetra_Vector> modsv = PoroField()->StructureField()->Interface()->InsertFSICondVector(FluidToStruct_FSI(fcvgfsi)); //(fvg)fg -> (fvg)sg -> (fvg)s
  Teuchos::RCP<Epetra_Vector> mmodsv = PoroField()->Extractor()->InsertVector(modsv,0); //(fvg)s -> (fvg)p

  mmodsv->Update(1.0, *sv, (1.0-stiparam)/(1.0-ftiparam)*fluidscale);

  Extractor().InsertVector(*mmodsv,poro_block_,f);  // add poroelast contributions to 'f'
  }
  else
    Extractor().InsertVector(*sv,poro_block_,f);  // add poroelast contributions to 'f'

  Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractOtherVector(av);
  Extractor().InsertVector(*aov,ale_i_block_,f);  // add ALE contributions to 'f'


}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupRHSLambda(Epetra_Vector& f)
{
  if (lambda_ != Teuchos::null) //FSI - Interface?
  {
    // get time integration parameters of structure and fluid time integrators
    // to enable consistent time integration among the fields
    const double stiparam = PoroField()->StructureField()->TimIntParam();
    const double ftiparam = FluidField()->TimIntParam();

    // project Lagrange multiplier field onto the master interface DOFs and consider temporal scaling
    Teuchos::RCP<Epetra_Vector> lambdafull = PoroField()->StructureField()->Interface()->InsertFSICondVector(FluidToStruct_FSI(lambda_)); //(lambda)fg -> (lambda)sg -> (lambda)s
    Teuchos::RCP<Epetra_Vector> lambdafull_p = PoroField()->Extractor()->InsertVector(lambdafull,0); //(lambda)s -> (lambda)p
    lambdafull_p->Scale(stiparam-(ftiparam*(1.0-stiparam))/(1.0-ftiparam));

    // add Lagrange multiplier
    Extractor().AddVector(*lambdafull_p,poro_block_,f);
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupRHSFirstIter(Epetra_Vector& f)
{
  // This method is directly take from FSImonolithic_fluidsplit. As at the moment no predictors are considered (does not improve the convergence a lot),
  // all terms coming from predictors are commented out (for later implementation).

  // get time integration parameters of structure and fluid time integrators
   // to enable consistent time integration among the fields
   const double stiparam = PoroField()->StructureField()->TimIntParam();
   const double ftiparam = FluidField()->TimIntParam();

   // some scaling factors for fluid
   const double timescale = FluidField()->TimeScaling();
   const double scale     = FluidField()->ResidualScaling();

   // old interface velocity of fluid field (FSI Cond Vector)
   const Teuchos::RCP<const Epetra_Vector> fveln = FluidField()->ExtractInterfaceVeln();

   // get fluid matrix
   const Teuchos::RCP<LINALG::BlockSparseMatrixBase> blockf = FluidField()->BlockSystemMatrix();

   // get fluid shape derivatives matrix
   const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();

   // get ale matrix
   Teuchos::RCP<LINALG::BlockSparseMatrixBase> blocka = AleField()->BlockSystemMatrix();

 #ifdef DEBUG
   if (blockf==Teuchos::null)  { dserror("Expected Teuchos::rcp to fluid block matrix."); }
   if (blocka==Teuchos::null)  { dserror("Expected Teuchos::rcp to ale block matrix."); }
 #endif

   // extract fluid and ale submatrices
   const LINALG::SparseMatrix& fig = blockf->Matrix(FLD::UTILS::MapExtractor::cond_other,FLD::UTILS::MapExtractor::cond_fsi); // F_{I\Gamma}
   const LINALG::SparseMatrix& fgg = blockf->Matrix(FLD::UTILS::MapExtractor::cond_fsi,FLD::UTILS::MapExtractor::cond_fsi); // F_{\Gamma\Gamma}
   //const LINALG::SparseMatrix& aig = blocka->Matrix(ALE::UTILS::MapExtractor::cond_other,ALE::UTILS::MapExtractor::cond_fsi); // A_{I\Gamma}

   // some often re-used vectors
   Teuchos::RCP<Epetra_Vector> rhs = Teuchos::null;

   // Different contributions/terms to the rhs are separated by the following comment line
   // ---------- structural interface DOFs
   /* The following terms are added to the structural interface DOFs of right hand side:
    *
    * rhs_firstnewtonstep =
    *
    * (1)  + (1-stiparam)/(1-ftiparam) * dt / tau * F_{\Gamma\Gamma} * u^{n}_{\Gamma}
    *
    * (2)  - (1-stiparam)/(1-ftiparam) / tau * F_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
    *
    * (3)  - (1-stiparam)/(1-ftiparam) * F^{G}_{\Gamma\Gamma} * \Delta d_{\Gamma,p}
    *
    * Remarks on all terms:
    * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
    *
    */
   // ----------addressing term 1
   rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(),true));

   fgg.Apply(*fveln,*rhs);

   rhs->Scale(scale * (1.-stiparam)/(1.-ftiparam) * Dt() * timescale);
   rhs = FluidToStruct_FSI(rhs);
   rhs = PoroField()->StructureField()->Interface()->InsertFSICondVector(rhs);
   rhs = PoroField()->Extractor()->InsertVector(rhs,0); //s->p

   if (PoroField()->StructureField()->GetSTCAlgo() == INPAR::STR::stc_currsym) //??ChrAg
   {
     Teuchos::RCP<LINALG::SparseMatrix> stcmat = PoroField()->StructureField()->GetSTCMat();
     stcmat->Multiply(true,*rhs,*rhs);
   }

   Extractor().AddVector(*rhs,0,f);
   // ----------end of term 1

//   // ----------addressing term 2:
//   rhs = Teuchos::rcp(new Epetra_Vector(fgg.RangeMap(),true));
//
//   fgg.Apply(*StructToFluid(ddgpred_), *rhs);
//
//   rhs->Scale(-scale * (1.-stiparam) / (1.-ftiparam) * timescale);
//   rhs = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(rhs));
//
//   Extractor().AddVector(*rhs,0,f);
//   // ----------end of term 2
//
//   // ----------addressing term 3:
//   if (mmm != Teuchos::null)
//   {
//     // extract F^{G}_{\Gamma\Gamma}
//     const LINALG::SparseMatrix& fmgg = mmm->Matrix(1,1);
//
//     rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RangeMap(),true));
//
//     fmgg.Apply(*StructToFluid(ddgpred_), *rhs);
//
//     rhs->Scale(-(1.-stiparam) / (1.-ftiparam));
//     rhs = StructureField()->Interface()->InsertFSICondVector(FluidToStruct(rhs));
//
//     Extractor().AddVector(*rhs,0,f);
//   }
//   // ----------end of term 3
   // ----------end of structural interface DOFs

   // ---------- inner fluid DOFs
   /* The following terms are added to the inner fluid DOFs of right hand side:
    *
    * rhs_firstnewtonstep =
    *
    * (1)  + dt / tau * F_{I \Gamma} * u^{n}_{\Gamma}
    *
    * (2)  - 1 / tau F_{I \Gamma} * \Delta d_{\Gamma,p}
    *
    * (3)  - F^{G}_{I \Gamma} * \Delta d_{\Gamma,p}
    *
    * Remarks on all terms:
    * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
    *
    */
   // ----------addressing term 1
   rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));

   fig.Apply(*fveln,*rhs);

   rhs->Scale(Dt() * timescale);

 #ifdef FLUIDSPLITAMG
   rhs = FluidField().Interface()->InsertOtherVector(rhs);
 #endif

   rhs = FluidField()->Interface()->InsertOtherVector(rhs);
   Extractor().AddVector(*rhs,1,f);
   // ----------end of term 1

//   // ----------addressing term 2
//   rhs = Teuchos::rcp(new Epetra_Vector(fig.RangeMap(),true));
//
//   fig.Apply(*StructToFluid(ddgpred_),*rhs);
//
//   rhs->Scale(-timescale);
//
// #ifdef FLUIDSPLITAMG
//   rhs = FluidField().Interface()->InsertOtherVector(rhs);
// #endif
//
//   Extractor().AddVector(*rhs,1,f);
//   // ----------end of term 2
//
//   // ----------addressing term 3
//   if(mmm != Teuchos::null)
//   {
//     // extract F^{G}_{I \Gamma}
//     const LINALG::SparseMatrix& fmig = mmm->Matrix(0,1);
//
//     rhs = Teuchos::rcp(new Epetra_Vector(fmig.RangeMap(),true));
//
//     fmig.Apply(*StructToFluid(ddgpred_),*rhs);
//
//     rhs->Scale(-1.);
//
// #ifdef FLUIDSPLITAMG
//   rhs = FluidField().Interface()->InsertOtherVector(rhs);
// #endif
//
//     Extractor().AddVector(*rhs,1,f);
//   }
//   // ----------end of term 3
//   // ----------end of inner fluid DOFs

   // ---------- inner ale DOFs
   /* The following terms are added to the inner ale DOFs of right hand side:
    *
    * rhs_firstnewtonstep =
    *
    * (1)  - A_{I \Gamma} * \Delta d_{\Gamma,p}
    *
    */
//   // ----------addressing term 1
//   rhs = Teuchos::rcp(new Epetra_Vector(aig.RangeMap(),true));
//
//   aig.Apply(*StructToAle(ddgpred_),*rhs);
//   rhs->Scale(-1.0);
//
//   Extractor().AddVector(*rhs,2,f);
//   // ----------end of term 1
   // ---------- end of inner ale DOFs

   // Reset quantities of previous iteration step since they still store values from the last time step
   ddginc_ = LINALG::CreateVector(*PoroField()->StructureField()->Interface()->FSICondMap(),true);
   duiinc_ = LINALG::CreateVector(*FluidField()->Interface()->OtherMap(),true);
   ddialeinc_ = LINALG::CreateVector(*AleField()->Interface()->OtherMap(),true);
   soliprev_ = Teuchos::null;
   solgprev_ = Teuchos::null;
   fgicur_ = Teuchos::null;
   fggcur_ = Teuchos::null;
   return;

  }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                 Teuchos::RCP<const Epetra_Vector>& px,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax,  ///< ale displacements
                                                 bool firstiter_)                       ///< firstiteration?
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic_Plain::ExtractFieldVectors");

  // porous medium
  px = Extractor().ExtractVector(x,poro_block_);
  Teuchos::RCP<Epetra_Vector> sx = PoroField()->Extractor()->ExtractVector(px,0); //Extract structural dofs

  // extract inner ALE solution increment
   Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x,ale_i_block_);

   // convert structure solution increment to ALE solution increment at the interface

   Teuchos::RCP<Epetra_Vector> scx_fsi = PoroField()->StructureField()->Interface()->ExtractFSICondVector(sx);

   Teuchos::RCP<const Epetra_Vector> acx_fsi = StructToAle_FSI(scx_fsi);

   // put inner --- ALE solution together
   Teuchos::RCP<Epetra_Vector> a = AleField()->Interface()->InsertOtherVector(aox);
   //AleField()->Interface()->InsertFPSICondVector(acx_fpsi, a); //Already done by Ale().ApplyInterfaceDisplacements()
   //AleField()->Interface()->InsertFSICondVector(acx_fsi, a); //Already done by Ale().ApplyInterfaceDisplacements()
   ax = a;

   // ---------------------------------------------------------------------------
   // process fluid unknowns
   // ---------------------------------------------------------------------------
   // extract inner fluid solution increment from NOX increment
   Teuchos::RCP<Epetra_Vector> f = Extractor().ExtractVector(x,fluid_block_);
   #ifdef FLUIDSPLITAMG
     fox = FluidField()->Interface()->ExtractOtherVector(fox);
   #endif

   // convert ALE solution increment to fluid solution increment at the interface
   Teuchos::RCP<Epetra_Vector> fcx_fsi = AleToFluidInterface_FSI(acx_fsi);

   if (firstiter_)
     FluidField()->DisplacementToVelocity(fcx_fsi); // Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1) - dt * u(n))
   else
     fcx_fsi->Scale(FluidField()->TimeScaling()); // Delta u(n+1,i+1) = fac * (Delta d(n+1,i+1)

   FluidField()->Interface()->InsertFSICondVector(fcx_fsi, f);
   Teuchos::RCP<Epetra_Vector> fox = Extractor_FluidFPSI_FSI().ExtractVector(f,0);

   fx = f;


   // ---------------------------------------------------------------------------

   // Store field vectors to know them later on as previous quantities
   // inner ale displacement increment
   // interface structure displacement increment
   if (disgprev_ != Teuchos::null)
     ddginc_->Update(1.0, *scx_fsi, -1.0, *disgprev_, 0.0);        // compute current iteration increment
   else
     ddginc_ = Teuchos::rcp(new Epetra_Vector(*scx_fsi));                   // first iteration increment

   disgprev_ = scx_fsi;                                            // store current step increment
   // ------------------------------------

   // inner ale displacement increment
   if (solialeprev_ != Teuchos::null)
     ddialeinc_->Update(1.0, *aox, -1.0, *solialeprev_, 0.0);  // compute current iteration increment
   else
     ddialeinc_ = Teuchos::rcp(new Epetra_Vector(*aox));                // first iteration increment

   solialeprev_ = aox;                                         // store current step increment
   // ------------------------------------

   // fluid solution increment
   if (soliprev_ != Teuchos::null)                             // compute current iteration increment
     duiinc_->Update(1.0, *fox, -1.0, *soliprev_, 0.0);
   else                                                        // first iteration increment
     duiinc_ =  Teuchos::rcp(new Epetra_Vector(*fox));
                                                               // store current step increment
   soliprev_ = fox;
   // ------------------------------------

}

/*----------------------------------------------------------------------*/
/* Recover the Lagrange multiplier at the interface   ... FPSI adapted version of mayr.mt (03/2012) */
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::RecoverLagrangeMultiplier()
{
  // This method is directly take from FSImonolithic_fluidsplit. As at the moment no predictors are considered (does not improve the convergence a lot),
  // all terms coming from predictors are commented out (for later implementation).

  if (lambda_ != Teuchos::null) //FSI - Interface?
  {
    // get time integration parameter of fluid time integrator
    // to enable consistent time integration among the fields
    const double ftiparam = FluidField()->TimIntParam();

    // some scaling factors for fluid
    const double timescale  = FluidField()->TimeScaling();
    const double scale      = FluidField()->ResidualScaling();

    // get fluid shape derivative matrix
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase> mmm = FluidField()->ShapeDerivatives();

    // some often re-used vectors
    Teuchos::RCP<Epetra_Vector> tmpvec = Teuchos::null;     // stores intermediate result of terms (3)-(8)
    Teuchos::RCP<Epetra_Vector> auxvec = Teuchos::null;     // just for convenience
    Teuchos::RCP<Epetra_Vector> auxauxvec = Teuchos::null;  // just for convenience

    /* Recovery of Lagrange multiplier \lambda_^{n+1} is done by the following
     * condensation expression:
     *
     * lambda_^{n+1} =
     *
     * (1)  - ftiparam / (1.-ftiparam) * lambda^{n}
     *
     * (2)  - 1. / (1.-ftiparam) * tmpvec
     *
     * with tmpvec =
     *
     * (3)    r_{\Gamma}^{F,n+1}
     *
     * (4)  + 1 / tau * F_{\Gamma\Gamma} * \Delta d_{\Gamma}^{S,n+1}
     *
     * (5)  + F_{\Gamma\Gamma}^{G} * \Delta d_{\Gamma}^{S,n+1}
     *
     * (6)  + F_{\Gamma I} * \Delta u_{I}^{F,n+1}
     *
     * (7)  + F_{\Gamma I}^{G} * \Delta d_{I}^{G,n+1}
     *
     * (8)  + dt / tau * F_{\Gamma\Gamma} * u_{\Gamma}^n]
     *
     * Remark on term (8):
     * Term (8) has to be considered only in the first Newton iteration.
     * Hence, it will usually not be computed since in general we need more
     * than one nonlinear iteration until convergence.
     *
     * Remarks on all terms:
     * +  Division by -(1.0 - ftiparam) will be done in the end
     *    since this is common to all terms
     * +  tau: time scaling factor for interface time integration (tau = 1/FluidField().TimeScaling())
     * +  neglecting terms (4)-(8) should not alter the results significantly
     *    since at the end of the time step the solution increments tend to zero.
     *
     *                                                 Matthias Mayr (10/2012)
     */

    // ---------Addressing term (1)
    lambda_->Update(ftiparam,*lambda_,0.0);
    // ---------End of term (1)

    // ---------Addressing term (3)
    Teuchos::RCP<Epetra_Vector> fluidresidual = FluidField()->Interface()->ExtractFSICondVector(FluidField()->RHS());
    fluidresidual->Scale(-1.0); // invert sign to obtain residual, not rhs
    tmpvec = Teuchos::rcp(new Epetra_Vector(*fluidresidual));
    // ---------End of term (3)

    // ---------Addressing term (4)
    auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(),true));

    fggprev_->Apply(*StructToFluid_FSI(ddginc_), *auxvec);
    tmpvec->Update(timescale,*auxvec,1.0);
    // ---------End of term (4)

    // ---------Addressing term (5)
    if (fmggprev_!=Teuchos::null)
    {
      auxvec = Teuchos::rcp(new Epetra_Vector(fmggprev_->RangeMap(),true));
      fmggprev_->Apply(*StructToFluid_FSI(ddginc_),*auxvec);
      tmpvec->Update(1.0,*auxvec,1.0);
    }
    // ---------End of term (5)

    // ---------Addressing term (6)
    auxvec = Teuchos::rcp(new Epetra_Vector(fgiprev_->RangeMap(),true));
    fgiprev_->Apply(*duiinc_,*auxvec);
    tmpvec->Update(1.0,*auxvec,1.0);
    // ---------End of term (6)

    // ---------Addressing term (7)
    if (fmgiprev_!=Teuchos::null)
    {
      /* For matrix-vector-product, the DomainMap() of the matrix and the Map() of the vector
       * have to match. DomainMap() contains inner velocity DOFs and all pressure DOFs.
       * The inner ale displacement increment is converted to the fluid map using AleToFluid().
       * This results in a map that contains all velocity but no pressure DOFs.
       *
       * We have to circumvent some trouble with Epetra_BlockMaps since we cannot split
       * an Epetra_BlockMap into inner and interface DOFs.
       *
       * We create a map extractor 'velothermap' in order to extract the inner velocity
       * DOFs after calling AleToFluid(). Afterwards, a second map extractor 'velotherpressuremapext'
       * is used to append pressure DOFs filled with zeros.
       *
       * Finally, maps match and matrix-vector-multiplication can be done.
       */

      // extract inner velocity DOFs after calling AleToFluid()
      Teuchos::RCP<Epetra_Map> velothermap = LINALG::SplitMap(*FluidField()->VelocityRowMap(),*InterfaceFluidAleCoupling_FSI().MasterDofMap());
      LINALG::MapExtractor velothermapext = LINALG::MapExtractor(*FluidField()->VelocityRowMap(),velothermap,false);
      auxvec = Teuchos::rcp(new Epetra_Vector(*velothermap, true));
      velothermapext.ExtractOtherVector(AleToFluid(AleField()->Interface()->InsertOtherVector(ddialeinc_)),auxvec);

      // add pressure DOFs
      LINALG::MapExtractor velotherpressuremapext = LINALG::MapExtractor(fmgiprev_->DomainMap(),velothermap);
      auxauxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->DomainMap(), true));
      velotherpressuremapext.InsertCondVector(auxvec,auxauxvec);

      // prepare vector to store result of matrix-vector-product
      auxvec = Teuchos::rcp(new Epetra_Vector(fmgiprev_->RangeMap(),true));

      // Now, do the actual matrix-vector-product
      fmgiprev_->Apply(*auxauxvec,*auxvec);
      tmpvec->Update(1.0,*auxvec,1.0);
    }
    // ---------End of term (7)

    // ---------Addressing term (8)
    if (firstcall_)
    {
      auxvec = Teuchos::rcp(new Epetra_Vector(fggprev_->RangeMap(),true));
      fggprev_->Apply(*FluidField()->ExtractInterfaceVeln(),*auxvec);
      tmpvec->Update(Dt()*timescale,*auxvec,1.0);
    }
    // ---------End of term (8)

    // ---------Addressing term (2)
    lambda_->Update(scale,*tmpvec,1.0); // scale with ResidualScaling() to get [N/m^2]
    // ---------End of term (2)

    // Finally, divide by (1.0-ftiparam) which is common to all terms
    lambda_->Scale(-1.0/(1.0-ftiparam));

    // Finally, the Lagrange multiplier 'lambda_' is recovered here.
    // It represents nodal forces acting onto the structure.
  }
  return;
}


