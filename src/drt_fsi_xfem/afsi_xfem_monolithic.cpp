/*!----------------------------------------------------------------------
\file afsi_xfem_monolithic.cpp
\brief  Control routine for monolithic eXtendedAleFluidStructureInteraction (XAFSI) solved via a classical Newton scheme
        taking into account changing fluid dofsets

\level 3

<pre>
\maintainer Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*----------------------------------------------------------------------*/

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_xfsi.H"
#include "../drt_adapter/ad_ale_fpsi.H"

#include "../drt_fsi/fsi_debugwriter.H"

#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"

#include "../drt_inpar/inpar_solver.H"
#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_structure/stru_aux.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
#include "../drt_constraint/constraint_manager.H"
#include "../drt_io/io_pstream.H"

#include "fsi_xfem_monolithic.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>

//poro ... & contact
#include "../drt_adapter/ad_str_fpsiwrapper.H"

#include "../drt_contact/contact_poro_lagrange_strategy.H"
#include "../drt_mortar/mortar_manager_base.H"
#include "../drt_contact/meshtying_contact_bridge.H"

#include "../drt_adapter/ad_str_poro_wrapper.H"


#include "../drt_fsi_xfem/afsi_xfem_monolithic.H"

#include "coupling_comm_manager.H"
/*----------------------------------------------------------------------*/
// constructor
/*----------------------------------------------------------------------*/
FSI::MonolithicAFSI_XFEM::MonolithicAFSI_XFEM(const Epetra_Comm& comm,
                                    const Teuchos::ParameterList& timeparams) // FSIDynamicParams
  : MonolithicXFEM(comm, timeparams)
{
  if (!HaveAle()) dserror("MonolithicAFSI_XFEM: For AFSI always Ale Fluid is required!");

  const int ndim = DRT::Problem::Instance()->NDim();

  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_->SetupCoupling(*FluidField()->Discretization(),
                       *AleField()->Discretization(),
                       *FluidField()->Discretization()->NodeRowMap(),
                       *AleField()->Discretization()->NodeRowMap(),
                       ndim,
                       false);

   Ale_Struct_coupling_ = Teuchos::rcp(new XFEM::Coupling_Comm_Manager(*AleField()->Discretization(),*StructureField()->Discretization(),"StructAleCoupling"));

  //build ale system matrix in splitted system
   AleField()->CreateSystemMatrix(AleField()->Interface());

    //create interfaces (mapextractors which are required for the coupling objects in the FPSI Coupl()!!!)
    //initialize xFluid now! - Coupling required filleddofset! -- Here first cut is done!!!
    FluidField()->Init();
    {
      // set initial field by given function
      // we do this here, since we have direct access to all necessary parameters
      const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
      INPAR::FLUID::InitialField initfield = DRT::INPUT::IntegralValue<INPAR::FLUID::InitialField>(fdyn,"INITIALFIELD");
      if(initfield != INPAR::FLUID::initfield_zero_field)
      {
        int startfuncno = fdyn.get<int>("STARTFUNCNO");
        if (initfield != INPAR::FLUID::initfield_field_by_function and
            initfield != INPAR::FLUID::initfield_disturbed_field_from_function)
        {
          startfuncno=-1;
        }
        FluidField()->SetInitialFlowField(initfield,startfuncno);
      }
    }

    meshmap_   = Teuchos::rcp(new LINALG::MapExtractor());
    meshmap_->Setup(*FluidField()->DiscretisationXFEM()->InitialDofRowMap(),FluidAleCoupling().MasterDofMap(),
        LINALG::SplitMap(*FluidField()->DiscretisationXFEM()->InitialDofRowMap(),*FluidAleCoupling().MasterDofMap()));

  //This has to be done after Init()!
  SetupCouplingObjects();
  return;
}

//-------------------------------------------------------------------------
// Assign Ale Sysmat
//-------------------------------------------------------------------------
 void FSI::MonolithicAFSI_XFEM::AddCouplingSysmat(
     Teuchos::RCP<LINALG::BlockSparseMatrixBase>& sysmat,
     Teuchos::RCP<LINALG::SparseMatrix>& s ,
     Teuchos::RCP<LINALG::SparseMatrix>& f,
     const double scaling_S,
     const double scaling_F)
 {
   //Get Idx of fluid and ale field map extractors
   const int &aidx_other = ALE::UTILS::MapExtractor::cond_other;

   const Teuchos::RCP<LINALG::BlockSparseMatrixBase> a      = AleField()  -> BlockSystemMatrix();

   //ALE Condensation
   LINALG::SparseMatrix& aii = a->Matrix(aidx_other,aidx_other);

   sysmat->Assign(ale_i_block_,ale_i_block_,LINALG::View,aii);

   //  //////////////////////////////////////////////
   //  //////                                  //////
   //  //////    Linearization of FluidField   //////
   //  //////    with respect to ale mesh      //////
   //  //////             motion               //////
   //  //////                                  //////
   //  //////////////////////////////////////////////

   // TODO: THIS IS STILL MISSING, BUT USUALLY DOES NOT HAVE A BIG INFLUENCE ON THE CONVERGENCE!!!
   return;

 }

 //-------------------------------------------------------------------------
 // Add Ale RHS
 //-------------------------------------------------------------------------
  void FSI::MonolithicAFSI_XFEM::AddCouplingRHS(Teuchos::RCP<Epetra_Vector>& rhs)
  {
    Teuchos::RCP<const Epetra_Vector> av = AleField()->RHS();
    Teuchos::RCP<Epetra_Vector> aov = AleField()->Interface()->ExtractOtherVector(av);
    Extractor().InsertVector(*aov,ale_i_block_,*rhs);  // add ALE contributions to 'rhs'
  }

 //-------------------------------------------------------------------------
 // Set Ale Interface Displacements
 //-------------------------------------------------------------------------
 void FSI::MonolithicAFSI_XFEM::SetAleInterfaceDisplacements()
 {
     //Sets structural conditioned Dispnp onto Ale
     Ale_Struct_coupling_->InsertVector(1,StructureField()->Dispnp(),0,AleField()->WriteAccessDispnp(),XFEM::Coupling_Comm_Manager::full_to_full);
     Teuchos::RCP<const Epetra_Vector> aledisplacements = AleToFluid(AleField()->Dispnp());
     //Set Fluid Dispnp
     meshmap_->InsertCondVector(aledisplacements,FluidField()->WriteAccessDispnp());

     // new grid velocity
     FluidField()->UpdateGridv();

     return;
 }

 //-------------------------------------------------------------------------
 // Evaluate Ale
 //-------------------------------------------------------------------------
 void FSI::MonolithicAFSI_XFEM::EvaluateAle(Teuchos::RCP<const Epetra_Vector> ax)
 {
   // ale field
   Epetra_Time ta(Comm());

   if (ax != Teuchos::null) //set just ale displacements, which are not set by a dirichlet boundary condition!!!
   {
     Teuchos::RCP<Epetra_Vector> DispnpAle = Teuchos::rcp(new Epetra_Vector(*AleField()->DofRowMap()),true);
     DispnpAle->Update(1.0,*AleField()->Interface()->InsertOtherVector(ax),1.0,*AleField()->Dispn(),0.0); //update ale disp here...
     AleField()->GetDBCMapExtractor()->InsertOtherVector(
         AleField()->GetDBCMapExtractor()->ExtractOtherVector(DispnpAle),AleField()->WriteAccessDispnp()); //just update displacements which are not on dbc condition
   }

   AleField()->Evaluate(); //Evaluate with Teuchos::null as dispnp_ in AleField() was just set before!
   SetAleInterfaceDisplacements();

   IO::cout << "ale time: " << ta.ElapsedTime() << IO::endl;
   return;
 }


 /*----------------------------------------------------------------------*
  | setup of the monolithic XFSI system,                                 |
  | setup a new combined block row map and a new block matrix            |
  *----------------------------------------------------------------------*/
 void FSI::MonolithicAFSI_XFEM::SetupSystem()
 {
   //Call Base Class
   MonolithicXFEM::SetupSystem();
 }

 /*----------------------------------------------------------------------*/
 // read restart data for monolithic XFSI system
 /*----------------------------------------------------------------------*/
 void FSI::MonolithicAFSI_XFEM::ReadRestart(int step)
 {
   //Call Base Class
   MonolithicXFEM::ReadRestart(step);

   //--------------------------------
   // read ale field
   AleField()->ReadRestart(step);
   FluidField()->CreateInitialState();
   FluidField()->ReadRestart(step);

 }

 /*----------------------------------------------------------------------*
  | prepare the time step for fluid and structure                        |
  *----------------------------------------------------------------------*/
 void FSI::MonolithicAFSI_XFEM::PrepareTimeStep()
 {
   //Call Base Class
   MonolithicXFEM::PrepareTimeStep();

   AleField()->PrepareTimeStep();
 }

 /*----------------------------------------------------------------------*
  | recover Lagrange multiplier (structural forces) needed for rhs in    |
  | next time step and update single fields                              |
  *----------------------------------------------------------------------*/
 void FSI::MonolithicAFSI_XFEM::Update()
 {
   //Call Base Class
   MonolithicXFEM::Update();

   AleField()->Update();
 }

 /*----------------------------------------------------------------------*
  | write output
  *----------------------------------------------------------------------*/
 void FSI::MonolithicAFSI_XFEM::Output()
 {
   //Call Base Class
   MonolithicXFEM::Output();

   AleField()->Output();
 }

 Teuchos::RCP<Epetra_Vector> FSI::MonolithicAFSI_XFEM::FluidToAle(Teuchos::RCP<const Epetra_Vector> iv) const
 {
   return coupfa_->MasterToSlave(iv);
 }
 Teuchos::RCP<Epetra_Vector> FSI::MonolithicAFSI_XFEM::AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
 {
   return coupfa_->SlaveToMaster(iv);
 }



