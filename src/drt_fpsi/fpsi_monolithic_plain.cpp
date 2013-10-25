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
#include "../drt_adapter/ad_str_fsiwrapper.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::Monolithic_Plain::Monolithic_Plain(const Epetra_Comm& comm,
                                         const Teuchos::ParameterList& fpsidynparams,
                                         const Teuchos::ParameterList& poroelastdynparams)
  : Monolithic(comm,fpsidynparams,poroelastdynparams)
{
  // create transformation object for the ale condensation
  aigtransform_           = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  kpatransform_           = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);

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

  vecSpaces.push_back(PoroField()     ->DofRowMap());
  //std::cout<<"PoroField()"<<*(PoroField()->DofRowMap())<<endl;
  vecSpaces.push_back(FluidField()    ->DofRowMap());
  //std::cout<<"FluidField()"<<*(FluidField()->Discretization()->DofRowMap())<<endl;
  vecSpaces.push_back(AleField()      ->Interface()->OtherMap());
//  std::cout<<"AleField()"<<*(AleField()->Interface()->OtherMap())<<endl;

  if(vecSpaces[0]->NumGlobalElements()==0)
    dserror("Hey I just met you and this is crazy, but here's my number "
            "so call me maybe ... !\n (0) 89 289 15300");
  if(vecSpaces[1]->NumGlobalElements()==0)
    dserror("No inner fluid equations. Splitting not possible. Roundhouse-Kick!");
  if(vecSpaces[2]->NumGlobalElements()==0)
    dserror("ALE ?! Roundhouse-Kick!");

  // merge maps and create full monolithic FPSI-DofRowMap
  SetDofRowMaps(vecSpaces);

  // switch fluid to interface split block matrix
  FluidField()->UseBlockMatrix(true);

  // build ale system matrix in splitted system
  AleField()->BuildSystemMatrix(false);

  // initialize FPSI-systemmatrix_
  systemmatrix_ = rcp(new LINALG::BlockSparseMatrix<
      LINALG::DefaultBlockMatrixStrategy>(Extractor(), Extractor(), 81, false,
      true));

  // create off-diagonal coupling matrices

  // Block Porous medium - Fluid
  k_pf_ = Teuchos::rcp(new LINALG::SparseMatrix(
                  *(PoroField()->DofRowMap()), 81, true, true));

  // Block Fluid - Porous medium
  k_fp_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(FluidField()->Discretization()->DofRowMap()), 81, true, true));

  k_pf_struct_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(PoroField()->StructureField()->DofRowMap()), 81, true, true));

  k_pf_porofluid_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(FluidField()->DofRowMap()),  81, true, true));

  k_fp_porofluid_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(PoroField()->FluidField()->Discretization()->DofRowMap()),  81, true, true));


  // Sub Block PoroFluid - Structure
  k_pfs_ = Teuchos::rcp(new LINALG::SparseMatrix(
                        *(PoroField()->FluidField()->Discretization()->DofRowMap()),
                          81, true, true));

  // create matrices needed for ale condensation

  k_pa_ = rcp(new LINALG::BlockSparseMatrix<
              LINALG::DefaultBlockMatrixStrategy>(*AleField()->Interface(),*FluidField()->FPSIInterface(), 81, false, true));

  k_ap_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(AleField()->Interface()->OtherMap()), 81, true, true));
  k_aa_ = Teuchos::rcp(new LINALG::SparseMatrix(
                      *(AleField()->Interface()->OtherMap()), 81, true, true));

  k_fa_ = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*AleField()->Interface(),*FluidField()->FPSIInterface(), 81, false, true));

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

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::SetupSystemMatrix(LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic_Plain::SetupSystemMatrix");

  //PoroField()->SetupSystemMatrix();

  const ADAPTER::Coupling& coupsa  = StructureAleCoupling();
  const ADAPTER::Coupling& coupsf  = StructureFluidCoupling();
  const ADAPTER::Coupling& smallcoupsf  = SmallStructureFluidCoupling();
  const ADAPTER::Coupling& coupfa  = FluidAleCoupling();


  // get single field block matrices
        Teuchos::RCP<LINALG::SparseMatrix>          p      = PoroField() -> SystemMatrix();
  const Teuchos::RCP<LINALG::SparseMatrix>          f      = FluidField()-> SystemSparseMatrix();
  const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a      = AleField()  -> BlockSystemMatrix();
  /*----------------------------------------------------------------------*/
  // build block matrix

  // insert poroelast structure
  p->UnComplete();
  mat.Assign(0,0,View,*p);

  //std::cout<<"p ColMap():   \n" <<p->ColMap()<<endl;
  //std::cout<<"mat ColMap(): \n" << mat.Matrix(0,0).ColMap()<<endl;

  // insert fluid
  f->UnComplete();
  mat.Assign(1,1,View,*f);

  // insert ale
  Teuchos::RCP<LINALG::SparseMatrix> k_pa = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pa_);
  Teuchos::RCP<LINALG::SparseMatrix> k_ap = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_ap_);
  Teuchos::RCP<LINALG::SparseMatrix> k_aa = Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_aa_);
  k_ap -> UnComplete();
  k_aa -> UnComplete();

    // condense ale
  LINALG::SparseMatrix& aii = a->Matrix(0,0);
  LINALG::SparseMatrix& aig = a->Matrix(0,1);

    Teuchos::RCP<LINALG::SparseMatrix>  laig = Teuchos::rcp(new LINALG::SparseMatrix(aii.RowMap(),81,false));
    (*aigtransform_)( a->FullRowMap(),
                      a->FullColMap(),
                      aig,
                      1.,
                      ADAPTER::CouplingSlaveConverter(coupsa),
                      *laig);

    laig->Complete(p->DomainMap(),laig->RangeMap());


    k_aa -> Add(  aii,false,1.0,1.0);
    k_ap -> Add(*laig,false,1.0,1.0);

    k_aa -> Complete(aii.DomainMap(),aii.RangeMap());
    k_ap -> Complete(p->DomainMap(),laig->RangeMap());
    mat.Assign(2,2,View,*k_aa);
    mat.Assign(2,0,View,*k_ap);

  // insert coupling terms

  Teuchos::RCP<LINALG::SparseMatrix> k_pf           = PoroFluidCouplingMatrix();
  Teuchos::RCP<LINALG::SparseMatrix> k_fp           = FluidPoroCouplingMatrix();

  ApplyCouplingTerms(k_pf,k_fp,p,f,k_pa_, k_fa_);

  Teuchos::RCP<LINALG::SparseMatrix> aa  = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()),81,false));
  aa->Assign(View,Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_pa_)->Matrix(1,0));
  aa->Complete(*AleField()->Interface()->OtherMap(),*FluidField()->FPSIInterface()->FSICondMap());

  Teuchos::RCP<LINALG::SparseMatrix> aaa = Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->DofRowMap()),81,false));
  (*couplingrowtransform5_)( *aa,
                              1.0,
                              ADAPTER::CouplingSlaveConverter(coupsf), // important to use slave converter
                             *aaa,
                              false);

  aaa -> Complete(*AleField()->Interface()->OtherMap(),p->RangeMap());
  mat.Assign(0, 2, View, *(aaa)); // assign poro_ale coupling matrix


  Teuchos::RCP<LINALG::SparseMatrix> aa2 = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()),81,false));
  aa2 -> Zero();
  Teuchos::RCP<LINALG::SparseMatrix> fluid_ale_gamma = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()),81,false));
  Teuchos::RCP<LINALG::SparseMatrix> fluid_ale_omega = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()),81,false));
  fluid_ale_omega -> Assign(View,Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_fa_)->Matrix(0,0));
  fluid_ale_omega -> Complete(*AleField()->Interface()->OtherMap(),f->RangeMap());
  fluid_ale_gamma -> Assign(View,Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(k_fa_)->Matrix(1,0));
  fluid_ale_gamma -> Complete(*AleField()->Interface()->OtherMap(),f->RangeMap());
  //std::cout<<"fluid_ale_omega: "<<*fluid_ale_omega<<endl;
  //std::cout<<"fluid_ale_gamma: "<<*fluid_ale_gamma<<endl;

  aa2 -> Add(*fluid_ale_omega,false,1.0,1.0);
  aa2 -> Add(*fluid_ale_gamma,false,1.0,1.0);

  aa2 -> Complete(AleField()->BlockSystemMatrix()->DomainMap(0),f->RangeMap());
  mat.Assign(1,2,View,*aa2); // assign fluid_ale coupling matrix (due to Neumann Integration)


  //////////////////////////////////////////////
  ///////                                 //////
  //////    Linearization of FluidField   //////
  //////    with respect to ale mesh      //////
  //////             motion               //////
  //////                                  //////
  //////////////////////////////////////////////
  const Teuchos::ParameterList& fpsidynparams   = DRT::Problem::Instance()->FPSIDynamicParams();
  if(Teuchos::getIntegralValue<int>(fpsidynparams,"USESHAPEDERIVATIVES"))
  {
    k_fp_->UnComplete();
    const Teuchos::RCP<LINALG::BlockSparseMatrixBase> fluidalematrix = FluidField()->ShapeDerivatives();

    if (fluidalematrix != Teuchos::null)
    {
      LINALG::SparseMatrix& fluidalematrix_ii = fluidalematrix->Matrix(0,0);
      LINALG::SparseMatrix& fluidalematrix_gg = fluidalematrix->Matrix(1,1);
      LINALG::SparseMatrix& fluidalematrix_gi = fluidalematrix->Matrix(1,0);
      LINALG::SparseMatrix& fluidalematrix_ig = fluidalematrix->Matrix(0,1);
      //std::cout<<fluidalematrix_ig<<endl;

      // add fluid_ale block ii and gi
      // those two blocks are not condensed since they belong to the columns of the inne ale dofs
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_ii,
          1.0,
          ADAPTER::CouplingMasterConverter(coupfa), // row converter: important to use slave converter
          mat.Matrix(1,2),
          false,
          true); // bool exactmatch = true (default)
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_gi,
          1.0,
          ADAPTER::CouplingMasterConverter(coupfa), // row converter: important to use slave converter
          mat.Matrix(1,2),
          false,
          true); // bool exactmatch = true (default)


      // reuse k_fp to assign linearization with respect to condensed interface dofs to the structural block
      Teuchos::RCP<LINALG::SparseMatrix> tempale1 = Teuchos::rcp(new LINALG::SparseMatrix(fluidalematrix_ig.RowMap(),81,false));
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_ig,
          1.0,
          ADAPTER::CouplingMasterConverter(smallcoupsf), // row converter: important to use slave converter
          *tempale1,
          false,  // bool exactmatch = true (default)
          false); // bool add
      Teuchos::RCP<LINALG::SparseMatrix> tempale2 = Teuchos::rcp(new LINALG::SparseMatrix(fluidalematrix_gg.RowMap(),81,false));
      (*couplingcoltransform_)( FluidField()->BlockSystemMatrix()->FullRowMap(),
          FluidField()->BlockSystemMatrix()->FullColMap(),
          fluidalematrix_gg,
          1.0,
          ADAPTER::CouplingMasterConverter(smallcoupsf), // row converter: important to use slave converter
          *tempale2,
          false,  // bool exactmatch = true (default)
          false); // bool add

      tempale1  -> Complete(*PoroField()->StructureRangeMap(),f->RangeMap());
      tempale2  -> Complete(*PoroField()->StructureRangeMap(),f->RangeMap());

      //Teuchos::RCP<LINALG::SparseMatrix> tempale12 = Teuchos::rcp(new LINALG::SparseMatrix((FluidField()->SystemSparseMatrix()->RowMap()),81,false));
      //tempale12->Add(*tempale2,false,1.0,1.0);
      //tempale12-> Complete(PoroField()->StructureField()->SystemMatrix()->RangeMap(),FluidField()->SystemSparseMatrix()->RangeMap());
      //std::cout<<*tempale12<<endl;


      k_fp -> Add(*tempale1,false,1.0,1.0);
      k_fp -> Add(*tempale2,false,1.0,1.0);
      k_fp -> Complete(p->RangeMap(),f->RangeMap());
    }
    else // if shapederivatives = no in FluidDynamics section in dat-file
    {
      std::cout<<"WARNING: Linearization with respect to mesh motion of fluid subproblem is switched off!"<<std::endl;
    }
  } // if useshapederivatives


  mat.Assign(0, 1, View, *(k_pf)); // assign poro_fluid coupling matrix
  mat.Assign(1, 0, View, *(k_fp)); // assign fluid_poro coupling matrix
  // done. make sure all blocks are filled.
  mat.Complete();

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
//
//  std::ostringstream filename2;
//  filename2 << "/scratch/rauch/workspace/baci/baci_debug/matlab_output/" << filebase << "_" << globindex << "_rhs";
//  LINALG::PrintVectorInMatlabFormat(filename2.str().c_str(),*rhs_);

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
  Extractor().InsertVector(*sv,0,f);  // add poroelast contributions to 'f'
  Extractor().InsertVector(*fv,1,f);  // add fluid contributions to 'f'

  Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractOtherVector(av);
  Extractor().InsertVector(*aov,2,f);  // add ALE contributions to 'f'
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic_Plain::ExtractFieldVectors(Teuchos::RCP<const Epetra_Vector> x,
                                                 Teuchos::RCP<const Epetra_Vector>& sx,
                                                 Teuchos::RCP<const Epetra_Vector>& fx,
                                                 Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic_Plain::ExtractFieldVectors");

  // porous medium
  sx = Extractor().ExtractVector(x,0);
  // fluid
  fx = Extractor().ExtractVector(x, 1);
  // ale
  ax = Extractor().ExtractVector(x, 2);

}
