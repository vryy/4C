
#ifdef CCADISCRET

#include "fsi_dirichletneumannslideale.H"
#include "fsi_debugwriter.H"
#include "../drt_geometry/searchtree.H"

extern struct _GENPROB genprob;


#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::DirichletNeumannSlideale::DirichletNeumannSlideale(Epetra_Comm& comm)
  : DirichletNeumann(comm)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  displacementcoupling_ = fsidyn.get<std::string>("COUPVARIABLE") == "Displacement";
  
  	//master objects in interface
	RCP<DRT::Discretization> masterdis = (StructureField().Discretization());
	map<int, RefCountPtr<DRT::Element> > masterelements;
	map<int, DRT::Node*> masternodes;
  map<int, DRT::Node*> mastergnodes;
  
	DRT::UTILS::FindConditionObjects(*masterdis, masternodes, mastergnodes, masterelements,"FSICoupling");
	imasternodes_ = mastergnodes;
	imastereles_ = masterelements;
	
  	//slave objects in interface
	RCP<DRT::Discretization> slavedis = MBFluidField().Discretization();
	map<int, RefCountPtr<DRT::Element> > slaveelements;
	map<int, DRT::Node*> slavenodes;
  map<int, DRT::Node*> slavegnodes;
	DRT::UTILS::FindConditionObjects(*slavedis, slavenodes, slavegnodes, slaveelements,"FSICoupling");
	islavenodes_ = slavenodes;
	
	//useful displacement vectors
	RCP<Epetra_Map> masterdofmap = StructureFluidCouplingMortar().MasterDofMap();
	RCP<Epetra_Map> slavecolmap = StructureFluidCouplingMortar().SlaveColMap();
	RCP<Epetra_Map> slavedofmap = StructureFluidCouplingMortar().SlaveDofMap();
	RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*masterdofmap,*slavedofmap, true);
	idispms_ = LINALG::CreateVector(*dofrowmap, true);
	islave_ = Teuchos::rcp(new Epetra_Vector(*slavedofmap,true));
	
	centerdisptotal_.resize(genprob.ndim);
	
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
	
  const Teuchos::ParameterList& input = DRT::Problem::Instance()->FSIDynamicParams();
  INPAR::FSI::SlideALEProj aletype = Teuchos::getIntegralValue<INPAR::FSI::SlideALEProj>(input,"SLIPALEPROJ");
	const int dim = genprob.ndim;
	
	//dispn and dispnp of structure, used for surface integral and velocity of the fluid in the interface 
	Teuchos::RCP<Epetra_Vector> idispn = StructureField().ExtractInterfaceDispn();
	Teuchos::RCP<Epetra_Vector> idisptotal = StructureField().ExtractInterfaceDispnp();
	Teuchos::RCP<Epetra_Vector> idispstep = StructureField().ExtractInterfaceDispnp();
	idispstep->Update(-1.0, *idispn, 1.0);
	
	map<int, RCP<DRT::Element> >& masterelements = imastereles_;
	RCP<DRT::Discretization> masterdis = (StructureField().Discretization());
	RCP<DRT::Discretization> slavedis = MBFluidField().Discretization();

	//for evaluation (surface integral)
	masterdis->SetState("displacementtotal",idisptotal); 		 
	masterdis->SetState("displacementincr",idispstep); 			
	Teuchos::ParameterList params;
	
	// define element matrices and vectors
	Epetra_SerialDenseMatrix elematrix1;
	Epetra_SerialDenseMatrix elematrix2;
	Epetra_SerialDenseVector elevector1;
	Epetra_SerialDenseVector elevector2;
	Epetra_SerialDenseVector elevector3;


	vector<double> mycenterdisp(dim);
	vector<double> centerdisp(dim);
	//variabel for length (2D) or area (3D) of the interface
	double mylengthcirc = 0.0;
	double lengthcirc = 0.0;

	
	//calculating the center displacement by loop over all struct elements
	map<int, RefCountPtr<DRT::Element> >::const_iterator elemiter;
	for (elemiter = masterelements.begin(); elemiter != masterelements.end(); ++elemiter)
	{
		RefCountPtr<DRT::Element> iele = elemiter->second;

		vector<int> lm;
		vector<int> lmowner;
		iele->LocationVector(*masterdis,lm,lmowner);
		elevector2.Size(1);		//length of circ with gaussinteg
		elevector3.Size(dim);		//centerdisp part of ele	

		params.set<string>("action","calc_struct_centerdisp");	    
		int err = iele->Evaluate(params,*masterdis,lm,elematrix1,elematrix2,elevector1,elevector2,elevector3);
		if (err) dserror("error while evaluating elements");
		mylengthcirc += elevector2[0];

		//disp of the interface
		for (int i=0; i<dim ;i++)
		{
			mycenterdisp[i] += elevector3[i]; 	
		}
	} //end of ele loop
	//masterdis->ClearState();

	//Assemble
	Comm().SumAll(&mylengthcirc, &lengthcirc, 1);
	Comm().SumAll(&mycenterdisp[0], &centerdisp[0], dim);

	if (lengthcirc == 0.0)
		dserror("no masterelements in the interface or wrong calculation");
		
	//calculating the final disp of the interface and summation over all time steps
	for (int i=0; i<dim ;i++)
	{
		centerdisp[i] = centerdisp[i] / lengthcirc;
		centerdisptotal_[i] +=centerdisp[i];
	}

	RCP<Epetra_Map> slavecolmap = StructureFluidCouplingMortar().SlaveColMap();
	RCP<Epetra_Map> slavedofmap = StructureFluidCouplingMortar().SlaveDofMap();
	RCP<Epetra_Map> masterdofmap = StructureFluidCouplingMortar().MasterDofMap();
	
	RCP<Epetra_Vector> islavestep = LINALG::CreateVector(*slavedofmap,true);
	
	//filling of islavestep with disp of the interface in every direction ("Schwerpunktsverschiebung")
	if (aletype==INPAR::FSI::ALEprojection_curr)
	{		
		//filling of islavestep for "Updated Lagrange Projektion"
		for (int iter=0; iter < (islavestep->MyLength()); iter++)
		{
			int a = iter % dim;
			int err = islavestep->ReplaceMyValues(1, &centerdisp[a], &iter);
			if (err == 1) dserror("error while replace values");
		}
		islave_->Update(1.0, *islavestep, 1.0);
	}
	else
	{
		//filling of islavestep for "Total Lagrange Projektion"
		for (int iter=0; iter < (islave_->MyLength()); iter++)
		{
			int a = iter % dim;
			int err = islave_->ReplaceMyValues(1, &centerdisptotal_[a], &iter);
			if (err == 1) dserror("error while replace values");
		}
	}
	
	//projection method 2D and 3D
	
	map<int, DRT::Node*> masternodes = imasternodes_ ;

	//currentpositions of masternodes for the search tree (always 3 coordinates)
	std::map<int,LINALG::Matrix<3,1> > currentpositions;
	currentpositions.clear();

	RefCountPtr<const Epetra_Vector> disp = masterdis->GetState("displacementtotal");
	
	DRT::Discretization& interfacedis = StructureFluidCouplingMortar().Interface()->Discret();
	RCP<Epetra_Map> msfullnodemap =  StructureFluidCouplingMortar().Interface()->MasterFullDofs();
	RCP<Epetra_Map> msfullelemap =  StructureFluidCouplingMortar().Interface()->MasterFullElements();

	RCP<Epetra_Import> interimpo = rcp (new Epetra_Import(*msfullnodemap,*masterdofmap));

	RefCountPtr<Epetra_Vector> reddisp = LINALG::CreateVector(*msfullnodemap,true);
	
	reddisp -> Import(*idisptotal,*interimpo,Add);
	
	map<int, RCP<DRT::Element> > masterreduelements;
	
	for (int eleind = 0; eleind<msfullelemap->NumMyElements(); eleind++)
	{
	  DRT::Element* tmpele = interfacedis.lColElement(eleind);
	  masterreduelements[tmpele->Id()]= rcp(tmpele,false);
	  
	  const int* n = tmpele->NodeIds();
	  
    for (int j=0; j < tmpele->NumNode(); j++)
    {
      const int gid = n[j];
      const DRT::Node* node = interfacedis.gNode(gid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      interfacedis.Dof(node, lm);
      vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;
      
      DRT::UTILS::ExtractMyValues(*reddisp,mydisp,lm);
      
      for (int a=0; a<3; a++)
      {
        currpos(a,0) = node->X()[a] + mydisp[a];
      }
      currentpositions[node->Id()] = currpos;
    }
	}
	//init of search tree
	Teuchos::RCP<GEO::SearchTree> searchTree = rcp(new GEO::SearchTree(8));
	const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofEles(masterreduelements, currentpositions);

	if(dim==2)
		searchTree->initializeTreeSlideALE(rootBox, masterreduelements, GEO::TreeType(GEO::QUADTREE));
	else if(dim==3)
		searchTree->initializeTreeSlideALE(rootBox, masterreduelements, GEO::TreeType(GEO::OCTTREE));
	else dserror("wrong dimension");
	
	
	map<int, DRT::Node*> slavenodes = islavenodes_;
	
	map<int, DRT::Node*>::const_iterator nodeiter;
	for (nodeiter = slavenodes.begin(); nodeiter != slavenodes.end(); ++nodeiter)
	{
		DRT::Node* node = nodeiter->second;
		vector<int> dofs = slavedis->Dof(node);		//gids
		vector<int> place(dim);
		for(int p=0; p<dim; p++)
			place[p] = slavedofmap->LID(dofs[p]);	//lids of gids

		Epetra_SerialDenseVector alenodecurr(dim);  // current coord of observed "alenode"

		//either "Updated Lagrange Projektion" or "Total Lagrange Projektion"
		if (aletype==INPAR::FSI::ALEprojection_curr)
		{
			for(int p=0; p<dim; p++)
				alenodecurr[p] =   node->X()[p]  + (*islave_)[(place[p])];	//access into islave_ only with lids
		}
		else
		{
			for(int p=0; p<dim; p++)
				alenodecurr[p] =   node->X()[p]  + centerdisptotal_[p];				
		}
		

		vector <double> finaldxyz(dim);

		//searchtree
		LINALG::Matrix<3,1>  querypoint;
		querypoint(0,0) = alenodecurr[0];
		querypoint(1,0) = alenodecurr[1];
		querypoint(2,0) = alenodecurr[2];

		//search for near elements next to the query point
		std::map<int,std::set<int> > 	closeeles = searchTree->searchElementsSlideALE(masterreduelements, currentpositions, querypoint);

		
		//search for the nearest point to project on
		if(dim == 2)
		{
			LINALG::Matrix<3,1> minDistCoords;
			GEO::nearestObjectInNode(masternodes,	masterreduelements, currentpositions, closeeles, querypoint, minDistCoords);
			finaldxyz[0] = minDistCoords(0,0) - alenodecurr[0]; 
			finaldxyz[1] = minDistCoords(1,0) - alenodecurr[1];
		}
		else
		{
			LINALG::Matrix<3,1> minDistCoords;
			GEO::nearestObjectInNode(rcp(&interfacedis,false), masterreduelements, currentpositions, closeeles,
					querypoint, minDistCoords);
			finaldxyz[0] = minDistCoords(0,0) - alenodecurr[0]; 
			finaldxyz[1] = minDistCoords(1,0) - alenodecurr[1];
			finaldxyz[2] = minDistCoords(2,0) - alenodecurr[2];
			
		}
		
		//refill of islavestep with values of the projection
		for(int p=0; p<dim; p++)
		{
			int err = islavestep->ReplaceMyValues(1, &finaldxyz[p], &place[p]);
			if (err == 1) dserror("error while replacing values");
		}

	}

	//summation in islave_
	islave_->Update(1.0, *islavestep, 1.0);

	//filling of idispms_ for mortar
	idispms_->Scale(0.0);
	
	RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*masterdofmap,*slavedofmap, true);
	RCP<Epetra_Import> msimpo = rcp (new Epetra_Import(*dofrowmap,*masterdofmap));
	RCP<Epetra_Import> slimpo = rcp (new Epetra_Import(*dofrowmap,*slavedofmap));

	idispms_ -> Import(*idisptotal,*msimpo,Add);
	idispms_ -> Import(*islave_,*slimpo,Add);
	
	//in here is the remeshing --> new D,M,Dinv out of disp of master and slave side
	StructureFluidCouplingMortar().Evaluate(idispms_);
	
	//interface velocity for fluid
	RCP<Epetra_Vector> ivel = LINALG::CreateVector(*masterdofmap,true);
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
		Teuchos::RCP<Epetra_Vector> iale = Teuchos::rcp(new Epetra_Vector(*(StructureFluidCouplingMortar().MasterDofMap()),true)); 

		Teuchos::RCP<Epetra_Vector> idispn = StructureField().ExtractInterfaceDispn();

		iale->Update(1.0, *idispcurr, 0.0);

		//iale reduced by old displacement dispn and instead added the real last displacements
		//iale->Update(1.0, *FTStemp_, -1.0, *idispn, 1.0);
		
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

		//FluidToStruct is dependent on the current configuration,
		//so only one time in time step is enough to calculate
	  
	  RCP<Epetra_Map> slavecolmap = StructureFluidCouplingMortar().SlaveColMap();
	  RCP<Epetra_Map> slavedofmap = StructureFluidCouplingMortar().SlaveDofMap();
  
	  
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
