/*----------------------------------------------------------------------*/
/*!
\file fpsi_coupling.cpp

\brief FPSI Coupling Object: Holds all objects on the Fluid-Poro-Interface and evaluates the Fluid-Poro-Coupling Matrixes!

<pre>
Maintainer: Ager Christoph & Andreas Rauch
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*/
///*----------------------------------------------------------------------*/
//GENERAL includes
//#include <sstream>
#include <Teuchos_TimeMonitor.hpp>
//#include <Teuchos_Time.hpp>
#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>
//
// POROELAST includes
#include "../drt_poroelast/poroelast_monolithic.H"
//
//// FSI includes
#include "../drt_fsi/fsi_matrixtransform.H"
//
// LINALG includes
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_blocksparsematrix.H"
//
// INPAR includes
#include "../drt_inpar/inpar_fpsi.H"
//
// LIB includes
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_assemblestrategy.H"
//
// ADAPTER includes
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_ale_fpsi.H"
//
// STRUCTURE includes
#include "../drt_structure/stru_aux.H"
//
// FLUID includes
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

// Header Include
#include "fpsi_coupling.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FPSI::FPSICoupling::FPSICoupling(Teuchos::RCP<POROELAST::Monolithic> poro,
                                 Teuchos::RCP< ::ADAPTER::Fluid> fluid,
                                 Teuchos::RCP< ::ADAPTER::AleFpsiWrapper> ale,
                                 Teuchos::RCP<std::map<int,int> > Fluid_PoroFluid_InterfaceMap,
                                 Teuchos::RCP<std::map<int,int> > PoroFluid_Fluid_InterfaceMap):
  poro_(poro),
  fluid_(fluid),
  ale_(ale),
  isfirstcall_(true),
  Fluid_PoroFluid_InterfaceMap_(Fluid_PoroFluid_InterfaceMap), //to be removed later
  PoroFluid_Fluid_InterfaceMap_(PoroFluid_Fluid_InterfaceMap),
  conductivity_(0.0)
{
  SetupInterfaceCoupling();
  InitCouplingMatrixesRHS();
  ReInitCouplingMatrixTransform();
  return;
}
/*----------------------------------------------------------------------/
| Initialize Coupling Matrixes and Coupling RHS              ager 12/14 |
/----------------------------------------------------------------------*/
void FPSI::FPSICoupling::InitCouplingMatrixesRHS()
{
  c_pp_ = Teuchos::RCP<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*poro_->DofRowMap(), 81, true, true));
  c_ff_ = Teuchos::RCP<LINALG::BlockSparseMatrixBase>(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*fluid_->Interface(),*fluid_->Interface(),81,true,true));
  c_pf_ = Teuchos::RCP<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*poro_->DofRowMap(), 81, true, true));
  c_fp_ = Teuchos::RCP<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*fluid_->DofRowMap(), 81, true, true));
  c_pa_ = Teuchos::RCP<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*poro_->DofRowMap(), 81, true, true));
  c_fa_ = Teuchos::RCP<LINALG::SparseMatrix>(new LINALG::SparseMatrix(*fluid_->DofRowMap(), 81, true, true));

  c_rhs_p_ = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*poro_->DofRowMap(),true));
  c_rhs_f_ = Teuchos::RCP<Epetra_Vector>(new Epetra_Vector(*fluid_->DofRowMap(),true));

  return;
}

/*----------------------------------------------------------------------/
| Setup the Coupling Object                                  ager 12/14 |
/----------------------------------------------------------------------*/
void FPSI::FPSICoupling::SetupInterfaceCoupling()
{
  const int ndim = DRT::Problem::Instance()->NDim();


  // porous fluid to fluid
  icoup_pf_f_  = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_pf_f_->SetupConditionCoupling(*PoroField()->FluidField()->Discretization(),
                               PoroField() ->FluidField()->FPSIInterface()->FPSICondMap(),
                              *FluidField()->Discretization(),
                               FluidField()->FPSIInterface()->FPSICondMap(),
                              "FPSICoupling",
                               ndim+1,
                               false);

// porous structure to fluid
  icoup_ps_f_  = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_ps_f_->SetupConditionCoupling(*PoroField()->StructureField()->Discretization(),
                               PoroField() ->StructureField()->Interface()->FPSICondMap(),
                              *FluidField()->Discretization(),
                               FluidField()->Interface()->FPSICondMap(),
                              "FPSICoupling",
                               ndim,
                               false);

// porous structure to ale
  icoup_ps_a_  = Teuchos::rcp(new ADAPTER::Coupling());
  icoup_ps_a_->SetupConditionCoupling(*PoroField()->StructureField()->Discretization(),
                               PoroField() ->StructureField()->Interface()->FPSICondMap(),
                              *AleField()  ->Discretization(),
                               AleField()  ->Interface()->FPSICondMap(),
                               "FPSICoupling",
                               ndim,
                               false);

  return;
}

/*-----------------------------------------------------------------------/
| Method reinitializes the matrix transformation objects      ager 12/14 |
/-----------------------------------------------------------------------*/
void FPSI::FPSICoupling::ReInitCouplingMatrixTransform()
{
  // create transformation objects for coupling terms
  couplingrowtransform_     = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform2_     = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform3_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform4_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingrowtransform5_    = Teuchos::rcp(new FSI::UTILS::MatrixRowTransform);
  couplingcoltransform_     = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  couplingcoltransform2_     = Teuchos::rcp(new FSI::UTILS::MatrixColTransform);
  couplingrowcoltransform_  = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
  couplingrowcoltransform2_ = Teuchos::rcp(new FSI::UTILS::MatrixRowColTransform);
}

/*-------------------------------------------------------------------------------/
| Evaluate Coupling Matrixes and Coupling RHS    orig. rauch / modif. ager 12/14 |
/-------------------------------------------------------------------------------*/
void FPSI::FPSICoupling::EvaluateCouplingMatrixesRHS()
{
  //Evaluates all Coupling Matrixes ...
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::FPSICoupling::EvaluateCouplingMatrixesRHS");

  Teuchos::RCP<LINALG::SparseMatrix> k_fp_porofluid = Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->FluidField()->DofRowMap()),81,true,true));
  Teuchos::RCP<LINALG::SparseMatrix> k_pf_porofluid = Teuchos::rcp(new LINALG::SparseMatrix(*(FluidField()->DofRowMap()),81,true,true));

  //set all coupling matrixes to zero!!
  c_pp_->Zero();
  c_ff_->Zero();
  c_pf_->Zero();
  c_fp_->Zero();
  c_fa_->Zero();
  c_pa_->Zero();

  k_pf_porofluid ->Zero();

  const ADAPTER::Coupling& couppff_fpsi     = *icoup_pf_f_;
  const ADAPTER::Coupling& coupsf_fpsi      = *icoup_ps_f_;
  const ADAPTER::Coupling& coup_ps_a_fpsi     = *icoup_ps_a_;

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  INPAR::FPSI::PartitionedCouplingMethod method = DRT::INPUT::IntegralValue<INPAR::FPSI::PartitionedCouplingMethod>(fpsidynparams,"PARTITIONED");

  if(method != INPAR::FPSI::nocoupling)
  {

  // set general vector values needed by elements

   PoroField()->FluidField()->Discretization()->ClearState();

   PoroField()->FluidField()->Discretization()->SetState(0,"dispnp" ,PoroField()->FluidField()->Dispnp());

   PoroField()->FluidField()->Discretization()->SetState(0,"gridv"  ,PoroField()->FluidField()->GridVel());
   PoroField()->FluidField()->Discretization()->SetState(0,"dispn"  ,PoroField()->FluidField()->Dispn());
   PoroField()->FluidField()->Discretization()->SetState(0,"veln"   ,PoroField()->FluidField()->Veln());
   PoroField()->FluidField()->Discretization()->SetState(0,"velaf"  ,PoroField()->FluidField()->Velnp());
   PoroField()->FluidField()->Discretization()->SetState(0,"velnp"  ,PoroField()->FluidField()->Velnp());

   FluidField()->Discretization()->ClearState();

   FluidField()->Discretization()->SetState(0,"dispnp",FluidField()->Dispnp());
   FluidField()->Discretization()->SetState(0,"gridv" ,FluidField()->GridVel());
   FluidField()->Discretization()->SetState(0,"dispn" ,FluidField()->Dispn());
   FluidField()->Discretization()->SetState(0,"veln"  ,FluidField()->Veln());
   FluidField()->Discretization()->SetState(0,"velaf" ,FluidField()->Velnp());
   FluidField()->Discretization()->SetState(0,"velnp" ,FluidField()->Velnp());


  // create the parameters for the discretization
  Teuchos::ParameterList fparams;

  // action for elements
  fparams.set<int>("action", FLD::fpsi_coupling);
  fparams.set("timescale",PoroField()->FluidField()->ResidualScaling());

  fparams.set("dt",fpsidynparams.get<double>("TIMESTEP"));
  fparams.set<int>("physical type", PoroField()->FluidField()->PhysicalType());


  if (method == INPAR::FPSI::monolithic)
  {
    fparams.set<std::string>("fillblock","Porofluid_Freefluid");
    fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
    DRT::AssembleStrategy fluidstrategy(
        0,              // fluiddofset for row
        0,              // fluiddofset for column
        k_pf_porofluid, // coupling matrix with fluid rowmap
        Teuchos::null,  // no other matrix or vectors
        Teuchos::null ,
        Teuchos::null,
        Teuchos::null
    );

    // what's the current problem type? Is it a fps3i problem?
    PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();

    if(probtype==prb_fps3i)
    {
      if(conductivity_==0.0)
      {
        dserror("In the case of FPS3I, a positive conductivity must be set in DESIGN SCATRA COUPLING SURF CONDITIONS");
      }
      fparams.set("membrane conductivity",conductivity_);
    }


    FluidField()->Discretization()->EvaluateCondition( fparams, fluidstrategy,"FPSICoupling" );
    k_pf_porofluid -> Complete(*FluidField()->DofRowMap(),*FluidField()->DofRowMap());
    Teuchos::RCP<LINALG::SparseMatrix> temp2 = Teuchos::rcp(new LINALG::SparseMatrix(*(PoroField()->FluidField()->DofRowMap()),81,false)); // temp matrix with porofluid rowmap

    {
    TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
    (*couplingrowtransform_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid),
                              1.0,
                              ADAPTER::CouplingSlaveConverter(couppff_fpsi),
                             *temp2,
                              false);
    }

    temp2          -> Complete(*FluidField()->DofRowMap(),PoroField()->FluidRangeMap());
    c_pf_->Add(*temp2,         false,1.0,0.0);

    fparams.set<std::string>("fillblock","Porofluid_Structure");
    fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
    k_pf_porofluid -> Zero();
    DRT::AssembleStrategy fluidstrategy21(
        0,                   // porofluiddofset for row
        0,                   // structuredofset for column
        k_pf_porofluid,     // coupling matrix with fluid rowmap
        Teuchos::null,       // no other matrix or vectors
        Teuchos::null,
        Teuchos::null,
        Teuchos::null
    );

    FluidField()    -> Discretization()->EvaluateCondition( fparams, fluidstrategy21, "FPSICoupling" );
    k_pf_porofluid -> Complete(*FluidField()->DofRowMap(),*FluidField()->DofRowMap());
    Teuchos::RCP<LINALG::SparseOperator> temp51 = Teuchos::rcp(new LINALG::SparseMatrix((*PoroField()->DofRowMap()),81,false));
    {
    TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
    (*couplingrowcoltransform2_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid),
        1.0,
        ADAPTER::CouplingSlaveConverter(couppff_fpsi), // row converter: important to use slave converter
        ADAPTER::CouplingSlaveConverter(coupsf_fpsi), //  col converter: important to use slave converter
        *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp51),
        false); // bool exactmatch = true (default)
    }
    temp51  -> Complete(*PoroField()->DofRowMap(),*PoroField()->DofRowMap());
    c_pp_->Add(*temp51,false,1.0,1.0);

    fparams.set<std::string>("fillblock","Fluid_Porofluid");
    fparams.set("InterfaceFacingElementMap", PoroFluid_Fluid_InterfaceMap_);
    k_fp_porofluid->Zero();
    DRT::AssembleStrategy porofluidstrategy(
        0,                   // porofluiddofset for row
        0,                   // porofluiddofset for column
        k_fp_porofluid,      // porofluid-structure coupling matrix
        Teuchos::null ,      // no other matrix or vectors
        Teuchos::null ,
        Teuchos::null ,
        Teuchos::null
    );


  PoroField()     -> FluidField()->Discretization()->EvaluateCondition( fparams, porofluidstrategy, "FPSICoupling" );
  k_fp_porofluid  -> Complete(PoroField()->FluidDomainMap(),PoroField()->FluidRangeMap());

  Teuchos::RCP<LINALG::SparseMatrix> temp = Teuchos::rcp(new LINALG::SparseMatrix(*FluidField()->DofRowMap(),81,false));
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingrowtransform2_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_fp_porofluid),
                             1.0,
                             ADAPTER::CouplingMasterConverter(couppff_fpsi),
                            *temp,
                             false); //true
  }
  temp  -> Complete(PoroField()->FluidDomainMap(),*FluidField()->DofRowMap());

  c_fp_ -> Add(*temp,false,1.0,1.0);


  fparams.set<std::string>("fillblock","Fluid_Structure");
  fparams.set("InterfaceFacingElementMap", PoroFluid_Fluid_InterfaceMap_);


  Teuchos::RCP<LINALG::SparseMatrix> k_pfs_ = Teuchos::rcp(new LINALG::SparseMatrix( //move me somewhere else
                        *(PoroField()->FluidField()->Discretization()->DofRowMap()),
                          81, true, true));

  k_pfs_ -> UnComplete();

  DRT::AssembleStrategy structurestrategy(
      0,                   // porofluiddofset for row
      1,                   // structuredofset for column
      k_pfs_,              // coupling matrix with porofluid rowmap
      Teuchos::null,       // no other matrix or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
  );

  PoroField() -> FluidField()->Discretization()->EvaluateCondition( fparams, structurestrategy, "FPSICoupling" );
  k_pfs_      -> Complete(PoroField()->StructureDomainMap(),PoroField()->FluidRangeMap());

  Teuchos::RCP<LINALG::SparseMatrix> temp3 = Teuchos::rcp(new LINALG::SparseMatrix(*FluidField()->DofRowMap(),81,false));
  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingrowtransform3_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pfs_),
                             1.0,
                             ADAPTER::CouplingMasterConverter(couppff_fpsi),
                            *temp3,
                             true);
  }

  temp3  -> Complete(PoroField()->StructureDomainMap(),*FluidField()->DofRowMap());

  c_fp_ -> Add(*temp3,false,1.0,1.0);

  ///// Fluid_Structure (fluid part / linearization of tangentials with respect to displacements)
  k_pf_porofluid->Reset();
  k_pf_porofluid -> UnComplete();

  DRT::AssembleStrategy structurestrategy2(
      0,                   // porofluiddofset for row
      0,                   // fluiddofset for column
      k_pf_porofluid,     // coupling matrix with porofluid rowmap
      Teuchos::null,       // no other matrix or vectors
      Teuchos::null,
      Teuchos::null,
      Teuchos::null
   );

  fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
  FluidField()->Discretization()->EvaluateCondition( fparams, structurestrategy2, "FPSICoupling" );

  k_pf_porofluid -> Complete(*FluidField()->DofRowMap(),*FluidField()->DofRowMap());
  Teuchos::RCP<LINALG::SparseMatrix> temp32 = Teuchos::rcp(new LINALG::SparseMatrix(*FluidField()->DofRowMap(),81,false));

  {
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
  (*couplingcoltransform_)( *FluidField()->DofRowMap(),
                            *FluidField()->DofRowMap(),
                            *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid),
                             1.0,
                             ADAPTER::CouplingSlaveConverter(coupsf_fpsi), // row converter: important to use slave converter
                            *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp32),
                             false); // bool exactmatch = true (default)
  }

  temp32  -> Complete(PoroField()->StructureDomainMap(),*FluidField()->DofRowMap());
  c_fp_ -> Add(*temp32,false,1.0,1.0);

  fparams.set<std::string>("fillblock","Fluid_Fluid");
  fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
    k_pf_porofluid -> Zero();
    k_pf_porofluid -> UnComplete();

    DRT::AssembleStrategy fluidfluidstrategy(
        0,                   // fluiddofset for row
        0,                   // fluiddofset for column
        c_ff_,     // porofluid-structure coupling matrix
        Teuchos::null,       // no other matrix or vectors
        Teuchos::null,
        Teuchos::null,
        Teuchos::null
    );

      FluidField()   ->Discretization()->EvaluateCondition( fparams, fluidfluidstrategy, "FPSICoupling" );

      fparams.set<std::string>("fillblock","Structure_Fluid");
      fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
      k_pf_porofluid -> Zero();
      k_pf_porofluid -> UnComplete();

      DRT::AssembleStrategy structurefluidstrategy(
          0,                   // fluid dofset for row
          0,                   // fluid dofset for column
          k_pf_porofluid ,    // coupling matrix with fluid rowmap
          Teuchos::null,       // no other matrix or vectors
          Teuchos::null,
          Teuchos::null,
          Teuchos::null
      );

      FluidField()   ->Discretization()->EvaluateCondition( fparams, structurefluidstrategy, "FPSICoupling" );

      k_pf_porofluid->Complete(*FluidField()->DofRowMap(), *FluidField()->DofRowMap());
      Teuchos::RCP<LINALG::SparseMatrix> temp4 = Teuchos::rcp(new LINALG::SparseMatrix((*PoroField()->StructureField()->DofRowMap()),81,false));
      {
      TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
      (*couplingrowtransform4_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid),
                                 1.0,
                                 ADAPTER::CouplingSlaveConverter(coupsf_fpsi), // important to use slave converter
                                *temp4,
                                 false);
      }

      temp4  -> Complete(*FluidField()->DofRowMap(),PoroField()->StructureRangeMap());


      c_pf_ -> Add(*temp4,false,1.0,1.0);

      fparams.set<std::string>("fillblock","Structure_Structure");
      fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
          k_pf_porofluid -> Zero();
          k_pf_porofluid -> UnComplete();

          DRT::AssembleStrategy structurestructurestrategy(
              0,                   // fluid dofset for row
              0,                   // fluid dofset for column
              k_pf_porofluid,     // coupling matrix with fluid rowmap
              Teuchos::null,       // no other matrix or vectors
              Teuchos::null,
              Teuchos::null,
              Teuchos::null
          );

          FluidField()   ->Discretization()->EvaluateCondition( fparams, structurestructurestrategy, "FPSICoupling" );
          // condense linearization with respect to the ale mesh motion (interface structural displacements = interface ale displacements)
          k_pf_porofluid->Complete(*FluidField()->DofRowMap(), *FluidField()->DofRowMap());

          Teuchos::RCP<LINALG::SparseOperator> temp5 = Teuchos::rcp(new LINALG::SparseMatrix((*PoroField()->DofRowMap()),81,false));
          {
          TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
          (*couplingrowcoltransform_)(*Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(k_pf_porofluid),
                                     1.0,
                                     ADAPTER::CouplingSlaveConverter(coupsf_fpsi), // row converter: important to use slave converter
                                     ADAPTER::CouplingSlaveConverter(coupsf_fpsi), // col converter: important to use slave converter
                                    *Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(temp5),
                                     false); // bool exactmatch = true (default)
          }

          temp5  -> Complete(*PoroField()->DofRowMap(),*PoroField()->DofRowMap());
          c_pp_      -> Add(*temp5,false,1.0,1.0);

      // Process inner ale dofs
          fparams.set<std::string>("fillblock","Structure_Ale");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);

          //temporal matrix
          //todo (initialization should be avoided in every iteration...)
          Teuchos::RCP<LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> > temp6 =
              Teuchos::rcp(new LINALG::BlockSparseMatrix<
                  LINALG::DefaultBlockMatrixStrategy>(*AleField()->Interface(),*FluidField()->Interface(), 81, false, false)); //Use FluidField()->Interface =

          //assemble into fluid row and column dofs -> need to transform rows to structure dofs and cols to ale dofs
          DRT::AssembleStrategy structurealestrategy(
              0,                   // fluid dofset for row
              1,                   // ale dofset for column
              temp6,               // coupling matrix with fluid rowmap
              Teuchos::null,       // no other matrix or vectors
              Teuchos::null,
              Teuchos::null,
              Teuchos::null
          );

          // evaluate coupling terms
          FluidField()   ->Discretization()->EvaluateCondition( fparams, structurealestrategy, "FPSICoupling" );
          temp6->Complete(); //for row transform!

          {
          TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::transform");
          (*couplingrowtransform5_)(temp6->Matrix(FLD::UTILS::MapExtractor::cond_fpsi,ALE::UTILS::MapExtractor::cond_other),
                                     1.0,
                                     ADAPTER::CouplingSlaveConverter(coupsf_fpsi), // important to use slave converter
                                    *c_pa_,
                                     false);
          }

  } // if monolithic
  else
  {
    if (isfirstcall_)
    std::cout<<"LINEARIZATION OF FPSI INTERFACE IS TURNED OFF \n"
               "IF YOU THINK THAT SUCKS. SET 'PARITIONED' to 'monolithic' \n"
               "IN THE FPSI SECTION OF YOUR INPUT FILE !!! "<<std::endl;
    isfirstcall_ = false;
  }


        //////////////////////////////////////
        ///////      __       ___       //////
        ///////      |_|  |_| |__       //////
        ///////      | \  | |  __|      //////
        //////                          //////
        //////////////////////////////////////

        Teuchos::RCP<Epetra_Vector> temprhs  = Teuchos::null;
        Teuchos::RCP<Epetra_Vector> temprhs2 = Teuchos::null;;

        fparams.set<std::string>("fillblock","conti");
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
        //temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(),true));
        //temprhs  = Teuchos::rcp(new Epetra_Vector(*Extractor().Map(fluid_block_),true)); // to get the full original map!
        temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->DofRowMap(),true)); // to get the full original map!
        temprhs2 = Teuchos::rcp(new Epetra_Vector(*PoroField() ->DofRowMap(),true));
        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        DRT::AssembleStrategy rhscontistrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        FluidField()   ->Discretization()->EvaluateCondition( fparams, rhscontistrategy, "FPSICoupling" );

        // extract FPSI part of the fluid field
        temprhs = FluidField()->FPSIInterface()->ExtractFPSICondVector(temprhs); //temprhs always has to be of full fluid dofrowmap

        // replace global fluid interface dofs through porofluid interface dofs
        temprhs = iFluidToPorofluid(temprhs);
        // insert porofluid interface entries into vector with full porofield length (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
        PoroField()->Interface().InsertVector(temprhs,3,temprhs2);
        // add vector with full porofield length to global rhs

        c_rhs_p_->Update(1.0,*temprhs2,0.0);

        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        fparams.set<std::string>("fillblock","structure");
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->Interface()->FullMap(),true)); //
        temprhs2 = Teuchos::rcp(new Epetra_Vector(*PoroField() ->Interface().FullMap(),true));

        DRT::AssembleStrategy rhsstructurestrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        FluidField()   ->Discretization()->EvaluateCondition( fparams, rhsstructurestrategy, "FPSICoupling" );

        // extract FPSI part of the fluid field
        temprhs = FluidField()->Interface()->ExtractFPSICondVector(temprhs); //
        // replace global fluid interface dofs through porofluid interface dofs
        temprhs = iFluidToPorostruct(temprhs); //
        // insert porofluid interface entries into vector with full porofield length (0: inner dofs of structure, 1: interface dofs of structure, 2: inner dofs of porofluid, 3: interface dofs of porofluid )
        PoroField()->Interface().InsertVector(temprhs,1,temprhs2);
        // add vector with full porofield length to global rhs

        c_rhs_p_->Update(1.0,*temprhs2,1.0);

        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        fparams.set<std::string>("fillblock","fluid");
        fparams.set("InterfaceFacingElementMap", PoroFluid_Fluid_InterfaceMap_);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*PoroField()->FluidField()->FPSIInterface()->FullMap(),true)); //
        temprhs2 = Teuchos::rcp(new Epetra_Vector(*FluidField()->FPSIInterface()->FullMap(),true));

        DRT::AssembleStrategy rhsfluidstrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        PoroField()->FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidstrategy, "FPSICoupling" );
        // extract FPSI part of the fluid field
        temprhs = PoroField()->FluidField()->FPSIInterface()->ExtractFPSICondVector(temprhs); //

        // replace global fluid interface dofs through porofluid interface dofs
        temprhs = iPorofluidToFluid(temprhs);
        // insert porofluid interface entries into vector with full fluidfield length
        FluidField()->FPSIInterface()->InsertVector(temprhs,FLD::UTILS::MapExtractor::cond_fpsi,temprhs2);
        // add vector with full porofield length to global rhs
        c_rhs_f_->Update(1.0,*temprhs2,0.0);

        temprhs ->PutScalar(0.0);
        temprhs2->PutScalar(0.0);

        fparams.set<std::string>("fillblock","fluidfluid"); // (wot,tangentialfac*uot) part
        fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
        temprhs  = Teuchos::rcp(new Epetra_Vector(*FluidField()->FPSIInterface()->FullMap(),true)); //


        DRT::AssembleStrategy rhsfluidfluidstrategy(
            0,                   // fluid dofset for row
            0,                   // fluid dofset for column
            Teuchos::null,
            Teuchos::null,
            temprhs,             // rhs vector
            Teuchos::null,
            Teuchos::null
        );

        FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy, "FPSICoupling" );

        c_rhs_f_->Update(1.0,*temprhs,1.0);

        temprhs ->PutScalar(0.0);
        temprhs2 ->PutScalar(0.0);

        //////////////////////////////////////
        //////                          //////
        //////   NEUMANN INTEGRATION    //////
        //////                          //////
        //////////////////////////////////////
        if (method == INPAR::FPSI::monolithic or method == INPAR::FPSI::RobinNeumann)
        {
          fparams.set<std::string>("fillblock","NeumannIntegration");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);

          Teuchos::RCP<LINALG::BlockSparseMatrixBase> tmp_c_ff = Teuchos::RCP<LINALG::BlockSparseMatrixBase>(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(*fluid_->Interface(),*fluid_->Interface(),81,true,true));

          DRT::AssembleStrategy rhsfluidfluidstrategy2(
              0,                   // fluid dofset for row
              0,                   // fluid dofset for column
              tmp_c_ff ,    // coupling matrix with fluid rowmap
              Teuchos::null,
              temprhs,             // rhs vector
              Teuchos::null,
              Teuchos::null
          );
          temprhs->PutScalar(0.0);
          FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy2, "NeumannIntegration" );

          c_rhs_f_->Update(1.0,*temprhs,1.0);

          tmp_c_ff->Complete();

          c_ff_ -> Add(*tmp_c_ff ,false,1.0,1.0);
          {
            fparams.set<std::string>("fillblock","NeumannIntegration_Ale");
            fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);

            Teuchos::RCP< LINALG::BlockSparseMatrixBase> tmp_c_fa = Teuchos::rcp(new LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>(
                *AleField()->Interface(),*FluidField()->FPSIInterface(), 81, false, true));

            DRT::AssembleStrategy rhsfluidfluidstrategy3(
                0,                   // fluid dofset for row
                1,                   // ale dofset for column
                tmp_c_fa,            // coupling matrix with fluid rowmap
                Teuchos::null,
                Teuchos::null,
                Teuchos::null,
                Teuchos::null
                );

            FluidField()->Discretization()->EvaluateCondition(fparams, rhsfluidfluidstrategy3, "NeumannIntegration");
            tmp_c_fa->Complete();

            //Add all inner ale parts to c_fa_ directly
            c_fa_->Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_fpsi, ALE::UTILS::MapExtractor::cond_other),false,1.0,0.0);
            c_fa_->Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_other),false,1.0,1.0);
            c_fa_->Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_other),false,1.0,1.0);

            //-->now transform ale fpsi block to structure (is condensed)!!!
            LINALG::SparseMatrix tmp_c_fp = LINALG::SparseMatrix(*FluidField()->DofRowMap(),81); //still in ale domain map!!!

            //Add all condensed parts to tmp_c_fa...
            tmp_c_fp.Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_fpsi, ALE::UTILS::MapExtractor::cond_fpsi),false,1.0,0.0);
            tmp_c_fp.Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_fpsi),false,1.0,1.0);
            tmp_c_fp.Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_fpsi),false,1.0,1.0);
            tmp_c_fp.Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_fpsi, ALE::UTILS::MapExtractor::cond_fsi),false,1.0,1.0);
            tmp_c_fp.Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_other, ALE::UTILS::MapExtractor::cond_fsi),false,1.0,1.0);
            tmp_c_fp.Add(tmp_c_fa->Matrix(FLD::UTILS::MapExtractor::cond_fsi, ALE::UTILS::MapExtractor::cond_fsi),false,1.0,1.0);
            tmp_c_fp.Complete(*AleField()->Interface()->FPSICondMap(),*FluidField()->DofRowMap());

            //For Ale Condensation ==> AleColumns to StructuralColumns
           (*couplingcoltransform2_)(AleField()->BlockSystemMatrix()->FullRowMap(),
                                     AleField()->BlockSystemMatrix()->FullColMap(),
                                     tmp_c_fp,
                                      1.0,
                                      ADAPTER::CouplingSlaveConverter(coup_ps_a_fpsi), // row converter: important to use slave converter
                                     *c_fp_,
                                      false,  // bool exactmatch = true (default)
                                      true); //bool add = false (default)

          }
        }
        else // only fill rhs
        {
          fparams.set<std::string>("fillblock","NeumannIntegration");
          fparams.set("InterfaceFacingElementMap", Fluid_PoroFluid_InterfaceMap_);
          temprhs->PutScalar(0.0);
          temprhs2->PutScalar(0.0);
          FluidField()->Discretization()->EvaluateCondition( fparams, rhsfluidfluidstrategy, "NeumannIntegration" );

          c_rhs_f_->Update(1.0,*temprhs,1.0);
        }

        ////////////////////////////
        // DONE -> CLEAR STATES
        PoroField() ->FluidField()->Discretization()->ClearState();
        FluidField()->Discretization()->ClearState();

  } // if not nocoupling

  ////////////////////////////
  // DONE -> Complete Coupling Matrixes
  c_ff_->Complete();
  c_pp_->Complete();

}

/*----------------------------------------------------------------------/
| set hydraulic conductivity
/----------------------------------------------------------------------*/
void FPSI::FPSICoupling::SetConductivity(double conduct)
{
  conductivity_ = conduct;
}
