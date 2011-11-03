//#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_fsi/fsi_dyn.H"

#include "../drt_fsi/fsi_monolithicfluidsplit.H"
#include "../drt_fsi/fsi_monolithiclagrange.H"
#include "../drt_fsi/fsi_monolithicstructuresplit.H"
#include "../drt_fsi/fsi_utils.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"

#include "../drt_fsi/fs_monolithic.H"

#include "../drt_fsi/fsi_nox_aitken.H"
#include "../drt_fsi/fsi_nox_group.H"
#include "../drt_fsi/fsi_nox_newton.H"
#include "../drt_fsi/fsi_statustest.H"

#include <NOX_Epetra_Interface_Jacobian.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_scatra/scatra_utils.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"


#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "fs3i_biofilm_growth.H"
#include "../drt_adapter/adapter_structure_bio.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


//#define SCATRABLOCKMATRIXMERGE


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::BiofilmGrowth::BiofilmGrowth(
	Epetra_Comm& comm)
:FS3I_1WC(comm)/*,
 ADAPTER::StructureBio(comm,	///< communicator
					   prbdyn, 	///< problem-specific parameters
					   1,  		///< we need an ALE formulation of the structure
					   0, 		///< scatra discretization number
					   "FSICoupling")*/
{	  // make sure that initial time derivative of concentration is not calculated
	  // automatically (i.e. field-wise)
	  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
	  if (DRT::INPUT::IntegralValue<int>(scatradyn,"SKIPINITDER")==false)
		  dserror("Initial time derivative of phi must not be calculated automatically -> set SKIPINITDER to false");
	  //fsi parameters
	  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
	  dt_fsi = fsidyn.get<double>("TIMESTEP");
	  nstep_fsi = fsidyn.get<int>("NUMSTEP");
	  maxtime_fsi = fsidyn.get<double>("MAXTIME");
	  step_fsi = 0;
	  time_fsi = 0.;
	  //surface growth parameters
	  const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();
	  dt_bio= biofilmcontrol.get<double>("BIOTIMESTEP");
	  nstep_bio= biofilmcontrol.get<int>("BIONUMSTEP");
	  grownvolume_ = biofilmcontrol.get<double>("GROWNVOLUME");
	  step_bio=0;
	  time_bio = 0.;

	  //total time
	  time_ = 0.;

	  idispn_= fsi_->FluidField().ExtractInterfaceVeln();
	  idispnp_= fsi_->FluidField().ExtractInterfaceVeln();
	  iveln_= fsi_->FluidField().ExtractInterfaceVeln();

	  idispn_->PutScalar(0.0);
	  idispnp_->PutScalar(0.0);
	  iveln_->PutScalar(0.0);


	  const string condname = "FSICoupling";

	  // set up ale-fluid couplings
	  icoupfa_.SetupConditionCoupling(*(fsi_->FluidField().Discretization()),
									  (fsi_->FluidField().Interface().FSICondMap()),
									  *(fsi_->AleField().Discretization()),
									  (fsi_->AleField().Interface().FSICondMap()),
									  condname,
	                                  genprob.ndim);

	  // the fluid-ale coupling always matches
	  const Epetra_Map* fluidnodemap = fsi_->FluidField().Discretization()->NodeRowMap();
	  const Epetra_Map* alenodemap   = fsi_->AleField().Discretization()->NodeRowMap();

	  coupfa_.SetupCoupling(*(fsi_->FluidField().Discretization()),
 							*(fsi_->AleField().Discretization()),
							*fluidnodemap,
							*alenodemap,
	                        genprob.ndim);

	  /// do we need this. What's for???
	  fsi_->FluidField().SetMeshMap(coupfa_.MasterDofMap());

	  Teuchos::RCP<ADAPTER::StructureBio> structbio = Teuchos::rcp(new ADAPTER::StructureBio(comm, fsidyn, 1, 0, "FSICoupling"));

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmGrowth::Timeloop()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  //outer loop for surface growth
  while (step_bio < nstep_bio)
  {
	  printf(" surface growth step = %3d   \n",step_bio);
	  printf(" Total time = %3f   \n",time_);

	  fsi_->SetupSystem();

	  // inner loop for fsi and scatra
	  InnerTimeloop();

	  // compute interface displacement and velocity
	  ComputeInterfaceVectors(idispnp_,iveln_);

	  if (idispnp_!=Teuchos::null)
	    {
	      // if we have values at the interface we need to apply them
	      fsi_->AleField().ApplyInterfaceDisplacements(FluidToAle(idispnp_));
	    }

	   // Note: We do not look for moving ale boundaries (outside the coupling
	   // interface) on the fluid side. Thus if you prescribe time variable ale
	   // Dirichlet conditions the according fluid Dirichlet conditions will not
	   // notice.
	  fsi_->AleField().BuildSystemMatrix();
	  fsi_->AleField().Solve();


	  //it is not necessary to apply the discplacement to the fluid and scatra field
	  //because during the fsi and scatra step (in the inner timeloop) they will first
	  //ask the ale field for the displacement

	  //Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluidField(fsi_->AleField().ExtractDisplacement());

	  // apply the displacement to the fluid field
	  //fsi_->FluidField().ApplyMeshDisplacement(fluiddisp);
	  //fsi_->FluidField().NonlinearSolve();

	  // apply the displacement to the scatra field
	  /*for (unsigned i=0; i<scatravec_.size(); ++i)  //// perhaps it has to be applied only to scatravec_[0] now
		  {
		  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
	   	  scatra->ScaTraField().ApplyMeshMovement(
	   			  fluiddisp,
	   			  fsi_->FluidField().Discretization()
				  );
		  }
	  LinearSolveScatra();*/

	  // do all the settings and solve the structure on a deforming mesh: still not properly implemented!!!
	  //StructAleSolve(idispnp_);

	  structbio->StructAleSolve(idispnp_);


	  step_bio++;
	  time_bio+=dt_bio;
	  time_ = time_bio + time_fsi;
  }

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmGrowth::InnerTimeloop()
{
  fsi_->PrepareTimeloop();

  double t=0.;
  step_fsi=0;

  while (step_fsi < nstep_fsi and t+1e-10*dt_fsi < maxtime_fsi)
  {
	DoFsiStep();
	DoScatraStep();

    step_fsi++;
    t+=dt_fsi;
  }
  time_fsi+=t;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmGrowth::DoScatraStep()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  if (comm.MyPID()==0)
  {
    cout<<"\n***********************\n GAS TRANSPORT SOLVER \n***********************\n";
  }

  // first scatra field is associated with fluid, second scatra field is
  // associated with structure

  bool stopnonliniter=false;
  int itnum = 0;

  PrepareTimeStep();

  while (stopnonliniter==false)
  {
    SetMeshDisp();
    SetVelocityFields();

    EvaluateScatraFields();

    SetupCoupledScatraSystem();

    stopnonliniter = AbortScatraNonlinIter(itnum);
    if (stopnonliniter)
      break;

    // transfer moving mesh data
    for (unsigned i=0; i<scatravec_.size(); ++i)
      {
        Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
        scatra->ScaTraField().ApplyMeshMovement(
        		fsi_->FluidField().Dispnp(),
        		fsi_->FluidField().Discretization()
        );
      }


    LinearSolveScatra();
    FieldUpdateIter();

    itnum++;
  }

  UpdateAndOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmGrowth::ComputeInterfaceVectors(
    RCP<Epetra_Vector> idispnp,
    RCP<Epetra_Vector> iveln)
{
	  // calculate normal flux vector field only at FSICoupling boundaries (no output to file)
	  vector<std::string> condnames(1);
	  condnames[0] = "FSICoupling";

	  Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[0];
	  Teuchos::RCP<Epetra_MultiVector> flux = scatra->ScaTraField().CalcFluxAtBoundary(condnames,false);

	  // access discretizations
	  RCP<DRT::Discretization> fluiddis = fsi_->FluidField().Discretization();
	  RCP<DRT::Discretization> scatradis = scatra->ScaTraField().Discretization();

	  // no support for growth caused by the mass transfer of multiple species
	  // id of the species causing the growth
	  int reactingspeciesid = 0;

	  const Epetra_BlockMap& ivelmap = iveln->Map();

	  // loop over all local nodes of fluid discretization
	  for (int lnodeid=0; lnodeid < fluiddis->NumMyRowNodes(); lnodeid++)
	  {
	    // Here we rely on the fact that the scatra discretization
	    // is a clone of the fluid mesh. => a scatra node has the same
	    // local (and global) ID as its corresponding fluid node!

		// get the processor's local fluid node with the same lnodeid
	    DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
	    // get the degrees of freedom associated with this fluid node
	    vector<int> fluidnodedofs = fluiddis->Dof(fluidlnode);

	    if(ivelmap.MyGID(fluidnodedofs[0])) // is this GID (implies: node) relevant for iveln_?
	    {
	      // determine number of space dimensions (numdof - pressure dof)
	      const int numdim = ((int) fluidnodedofs.size()) -1;
	      // number of dof per node in ScaTra
	      int numscatradof = scatradis->NumDof(scatradis->lRowNode(lnodeid));

	      vector<double> Values(numdim);
	      for(int index=0;index<numdim;++index)
	      {
	        const int pos = lnodeid*numscatradof+reactingspeciesid;
	        // interface growth has opposite direction respect to mass flow -> minus sign
	        Values[index] = (-grownvolume_)*((*flux)[index])[pos];
	      }

	      // now insert only the first numdim entries (pressure dof is not inserted!)
	      int error = iveln_->ReplaceGlobalValues(numdim,&Values[0],&fluidnodedofs[0]);
	      if (error > 0) dserror("Could not insert values into vector iveln_: error %d",error);
	    }
	  }

	  // have to compute an approximate displacement from given interface velocity
	  // id^{n+1} = id^{n} + \delta t vel_i
	  idispnp->Update(1.0,*idispn_,0.0);
	  idispnp->Update(dt_bio,*iveln_,1.0);
	  cout<<"grownvolume="<<grownvolume_<<endl;
	  cout<<"idispnp="<<*idispnp<<endl;
	return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmGrowth::FluidToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FS3I::BiofilmGrowth::AleToFluidField(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::BiofilmGrowth::UpdateAndOutput()
{
  for (unsigned i=0; i<scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField().Update();

    // perform time shift of interface displacement
    idispn_->Update(1.0, *idispnp_ , 0.0);

    // in order to do not have artificial velocity
    idispnp_=idispn_;

    scatra->ScaTraField().Output();
  }
}
