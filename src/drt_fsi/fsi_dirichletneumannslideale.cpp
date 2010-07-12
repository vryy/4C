
#ifdef CCADISCRET

#include "fsi_dirichletneumannslideale.H"
#include "fsi_debugwriter.H"
#include "fsi_utils.H"
#include "../drt_geometry/searchtree.H"
#include "../drt_mortar/mortar_interface.H"
#include "../drt_inpar/inpar_fsi.H"

extern struct _GENPROB genprob;


#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannSlideale::DirichletNeumannSlideale(Epetra_Comm& comm)
  : DirichletNeumann(comm)
{
  
  displacementcoupling_ = 
      DRT::Problem::Instance()->FSIDynamicParams().get<std::string>("COUPVARIABLE") == "Displacement";
  
  // get discretizations
	RCP<DRT::Discretization> masterdis = (StructureField().Discretization());
	RCP<DRT::Discretization> slavedis = MBFluidField().Discretization();
  
	slideale_ = rcp(new FSI::UTILS::SlideAleUtils(masterdis,slavedis,StructureFluidCouplingMortar()));

	islave_ = Teuchos::rcp(new Epetra_Vector(*StructureFluidCouplingMortar().SlaveDofRowMap(),true));
	
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannSlideale::FSIOp(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  if (displacementcoupling_)
  {
    
    const Teuchos::RCP<Epetra_Vector> idispn = rcp(new Epetra_Vector(x));
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("idispn",*idispn);

    const Teuchos::RCP<Epetra_Vector> iforce = FluidOp(idispn, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("iforce",*iforce);

    const Teuchos::RCP<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("idispnp",*idispnp);
   
    F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);
  }
  else
  {
    const Teuchos::RCP<Epetra_Vector> iforcen = rcp(new Epetra_Vector(x));
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("iforcen",*iforcen);

    const Teuchos::RCP<Epetra_Vector> idisp = StructOp(iforcen, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("idisp",*idisp);

    const Teuchos::RCP<Epetra_Vector> iforcenp = FluidOp(idisp, fillFlag);
    if (MyDebugWriter()!=Teuchos::null)
      MyDebugWriter()->WriteVector("iforcenp",*iforcenp);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannSlideale::Remeshing()
{
	
	//dispn and dispnp of structure, used for surface integral and velocity of the fluid in the interface 
	Teuchos::RCP<Epetra_Vector> idispn = StructureField().ExtractInterfaceDispn();
	Teuchos::RCP<Epetra_Vector> idisptotal = StructureField().ExtractInterfaceDispnp();
	Teuchos::RCP<Epetra_Vector> idispstep = StructureField().ExtractInterfaceDispnp();
	idispstep->Update(-1.0, *idispn, 1.0);
	
  INPAR::FSI::SlideALEProj aletype = 
      Teuchos::getIntegralValue<INPAR::FSI::SlideALEProj>(DRT::Problem::Instance()->FSIDynamicParams(),"SLIDEALEPROJ");
	
	slideale_->Remeshing(idisptotal,
                        idispstep,
                        (StructureField().Discretization()),
                        MBFluidField().Discretization(),
                        islave_,
                        StructureFluidCouplingMortar(),
                        Comm(),
                        aletype);
	
	RCP<Epetra_Map> masterdofrowmap = StructureFluidCouplingMortar().MasterDofRowMap();
	
	//interface velocity for fluid
	RCP<Epetra_Vector> ivel = LINALG::CreateVector(*masterdofrowmap,true);
	ivel->Update(1./Dt(), *idispstep, 0.0);
	
	Teuchos::RCP<Epetra_Vector> ivelfluid = StructToFluid(ivel);
	//solve ale and fluid again with known correct solution
	MBFluidField().NonlinearSolve(islave_,ivelfluid); 
	//solve of structure not needed because no changes were applied

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::DirichletNeumannSlideale::FluidOp(Teuchos::RCP<Epetra_Vector> idispcurr,
		const FillType fillFlag)
{
	FSI::Partitioned::FluidOp(idispcurr,fillFlag);

	if (fillFlag==User)
	{
		dserror("not implemented");
		// SD relaxation calculation
		return FluidToStruct(MBFluidField().RelaxationSolve(StructToFluid(idispcurr),Dt()));
	}
	else
	{
		// normal fluid solve

		// the displacement -> velocity conversion at the interface
		const Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idispcurr);

		// A rather simple hack. We need something better!
		const int itemax = MBFluidField().Itemax();
		if (fillFlag==MF_Res and mfresitemax_ > 0)
			MBFluidField().SetItemax(mfresitemax_ + 1);

		//new Epetra_Vector for aledisp in interface
		Teuchos::RCP<Epetra_Vector> iale = Teuchos::rcp(new Epetra_Vector(*(StructureFluidCouplingMortar().MasterDofRowMap()),true)); 

		Teuchos::RCP<Epetra_Vector> idispn = StructureField().ExtractInterfaceDispn();

		iale->Update(1.0, *idispcurr, 0.0);

		//iale reduced by old displacement dispn and instead added the real last displacements
		iale->Update(1.0, *FTStemp_, -1.0, *idispn, 1.0);
		
		MBFluidField().NonlinearSolve(StructToFluid(iale),StructToFluid(ivel));

		MBFluidField().SetItemax(itemax);

		return FluidToStruct(MBFluidField().ExtractInterfaceForces());
	}
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::DirichletNeumannSlideale::StructOp(Teuchos::RCP<Epetra_Vector> iforce,
		const FillType fillFlag)
{
	FSI::Partitioned::StructOp(iforce,fillFlag);

	if (fillFlag==User)
	{
		// SD relaxation calculation
		return StructureField().RelaxationSolve(iforce);
	}
	else
	{
		// normal structure solve
		StructureField().ApplyInterfaceForces(iforce);
		StructureField().Solve();
		return StructureField().ExtractInterfaceDispnp();
	}
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannSlideale::InitialGuess()
{
	if (displacementcoupling_)
	{
		//real displacement of slave side at time step begin on master side --> for calcualtion of FluidOp 
		FTStemp_ = FluidToStruct(islave_);
		// predict displacement
		return StructureField().PredictInterfaceDispnp();
	}
	else
	{
		const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
		if (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR")!=1)
		{
			dserror("unknown interface force predictor '%s'",
					fsidyn.get<string>("PREDICTOR").c_str());
		}
		return InterfaceForce();
	}
}

#endif
