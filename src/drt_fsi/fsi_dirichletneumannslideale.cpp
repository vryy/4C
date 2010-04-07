
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
  
  displacementcoupling_ = 
      DRT::Problem::Instance()->FSIDynamicParams().get<std::string>("COUPVARIABLE") == "Displacement";
  
  // declare master objects in interface
	RCP<DRT::Discretization> masterdis = (StructureField().Discretization());
	map<int, RefCountPtr<DRT::Element> > masterelements;
	map<int, DRT::Node*> masternodes;
  map<int, DRT::Node*> mastergnodes;
  
  // initialize master objects in interface
	DRT::UTILS::FindConditionObjects(*masterdis, masternodes, mastergnodes, masterelements,"FSICoupling");
	imastergnodes_ = mastergnodes;
	imastereles_ = masterelements;
	
  // declare slave objects in interface
	RCP<DRT::Discretization> slavedis = MBFluidField().Discretization();
	map<int, RefCountPtr<DRT::Element> > slaveelements;
	map<int, DRT::Node*> slavenodes;  // complete map of slave nodes 
	map<int, DRT::Node*> slavemnodes; // partial map of sliding slave nodes 
  map<int, DRT::Node*> slavegnodes; // dummy map
  
  
  //initialize slave objects in interface
	DRT::UTILS::FindConditionedNodes(*slavedis, "FSICoupling", slavenodes);
	DRT::UTILS::FindConditionedNodes(*slavedis, "FSICouplingNoSlide", slavemnodes);
	islaveconfnodes_ = slavemnodes;
	islaveslidnodes_ = slavenodes;
	map<int, DRT::Node*>::iterator it;
	for ( it=islaveconfnodes_.begin() ; it != islaveconfnodes_.end(); it++ )
	{
	  int err = islaveslidnodes_.erase((*it).first);
	  if (!err)
	    dserror("Non sliding interface has to be a subset of FSI-interface or empty");
	}
	
	// useful displacement vectors
	RCP<Epetra_Map> masterdofrowmap = StructureFluidCouplingMortar().MasterDofRowMap();
	RCP<Epetra_Map> slavedofrowmap = StructureFluidCouplingMortar().SlaveDofRowMap();
	RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*masterdofrowmap,*slavedofrowmap, true);
	idispms_ = LINALG::CreateVector(*dofrowmap, true);
	islave_ = Teuchos::rcp(new Epetra_Vector(*slavedofrowmap,true));
	
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
	
	const int dim = genprob.ndim;
	
	//dispn and dispnp of structure, used for surface integral and velocity of the fluid in the interface 
	Teuchos::RCP<Epetra_Vector> idispn = StructureField().ExtractInterfaceDispn();
	Teuchos::RCP<Epetra_Vector> idisptotal = StructureField().ExtractInterfaceDispnp();
	Teuchos::RCP<Epetra_Vector> idispstep = StructureField().ExtractInterfaceDispnp();
	idispstep->Update(-1.0, *idispn, 1.0);
	
	//calculate structural interface center of gravity
	vector<double> centerdisp = Centerdisp(idisptotal, idispstep);	
	
	// Replace ALE disp by average interface translation (rot free)
	RCP<Epetra_Map> slavedofrowmap = StructureFluidCouplingMortar().SlaveDofRowMap();
	RCP<Epetra_Map> masterdofrowmap = StructureFluidCouplingMortar().MasterDofRowMap();
	
	RCP<Epetra_Vector> islavestep = LINALG::CreateVector(*slavedofrowmap,true);
	
  // We need the master elements on every processor for the projection of the slave nodes.
  // Furthermore we need the current position of the masternodes on every processor.
  // Elements provided by interface discretization, necessary maps provided by interface.
  DRT::Discretization& interfacedis = StructureFluidCouplingMortar().Interface()->Discret();
  RCP<Epetra_Map> msfullnodemap =  StructureFluidCouplingMortar().Interface()->MasterFullDofs();
  RCP<Epetra_Map> msfullelemap =  StructureFluidCouplingMortar().Interface()->MasterFullElements();

  // Redistribute displacement of masternodes on the interface to all processors.
  RCP<Epetra_Import> interimpo = rcp (new Epetra_Import(*msfullnodemap,*masterdofrowmap));
  RCP<Epetra_Vector> reddisp = LINALG::CreateVector(*msfullnodemap,true);
  reddisp -> Import(*idisptotal,*interimpo,Add);
  
  // map with fully reduced master element distribution 
  map<int, RCP<DRT::Element> > masterreduelements;
  for (int eleind = 0; eleind<msfullelemap->NumMyElements(); eleind++)
  {
    DRT::Element* tmpele = interfacedis.lColElement(eleind);
    masterreduelements[tmpele->Id()]= rcp(tmpele,false);
  }
	
	//currentpositions of master nodes for the search tree (always 3 coordinates)
	std::map<int,LINALG::Matrix<3,1> > currentpositions = 
	    CurrentMasterPos(idisptotal,reddisp,interfacedis,msfullelemap);
	
	//project sliding slave nodes onto master interface surface
	islavestep->PutScalar(0.0);
	SlideProjection(islavestep,interfacedis,centerdisp,masterreduelements,currentpositions);

	//For the NON sliding ALE Nodes, use structure displacements
	RCP<Epetra_Vector> istrudis = StructToFluid(idisptotal);
	
	map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = islaveconfnodes_.begin(); nodeiter != islaveconfnodes_.end(); ++nodeiter)
  {
    DRT::Node* node = nodeiter->second;
    vector<int> lids(dim);
    for(int p=0; p<dim; p++)
      //lids of gids of node
      lids[p] = slavedofrowmap->LID((MBFluidField().Discretization()->Dof(node))[p]);

    // current coord of ale node = ref coord + islave_
    vector <double> finaldxyz(dim); 

    for(int p=0; p<dim; p++)
      finaldxyz[p] =   (*istrudis)[(lids[p])]; 
      
    int err = islavestep->ReplaceMyValues(dim, &finaldxyz[0], &lids[0]);
    if (err == 1) dserror("error while replacing values");

  }
  
  //put everything together
  islave_->Update(1.0, *islavestep, 1.0);
  
  
	//merge displacement values of interface nodes (master+slave) into idispms_ for mortar
	idispms_->Scale(0.0);
	
	RCP<Epetra_Map> dofrowmap = LINALG::MergeMap(*masterdofrowmap,*slavedofrowmap, true);
	RCP<Epetra_Import> msimpo = rcp (new Epetra_Import(*dofrowmap,*masterdofrowmap));
	RCP<Epetra_Import> slimpo = rcp (new Epetra_Import(*dofrowmap,*slavedofrowmap));

	idispms_ -> Import(*idisptotal,*msimpo,Add);
	idispms_ -> Import(*islave_,*slimpo,Add);
	
	//new D,M,Dinv out of disp of master and slave side
	StructureFluidCouplingMortar().Evaluate(idispms_);
	
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
vector<double> FSI::DirichletNeumannSlideale::Centerdisp
(
    Teuchos::RCP<Epetra_Vector> idisptotal,
    Teuchos::RCP<Epetra_Vector>idispstep
)
{
  const int dim = genprob.ndim;
  // get structure and fluid discretizations  and set stated for element evaluation 
  RCP<DRT::Discretization> masterdis = (StructureField().Discretization());
  
  masterdis->SetState("displacementtotal",idisptotal);     
  masterdis->SetState("displacementincr",idispstep);
  
  //define stuff needed by the elements
  Teuchos::ParameterList params;
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;

  //prepare variables for length (2D) or area (3D) of the interface
  vector<double> mycenterdisp(dim);
  vector<double> centerdisp(dim); 
  double mylengthcirc = 0.0;
  double lengthcirc = 0.0;
  
  //calculating the center displacement by evaluating structure interface elements
  map<int, RefCountPtr<DRT::Element> >::const_iterator elemiter;
  for (elemiter = imastereles_.begin(); elemiter != imastereles_.end(); ++elemiter)
  {
    
    RefCountPtr<DRT::Element> iele = elemiter->second;
    vector<int> lm;
    vector<int> lmowner;
    iele->LocationVector(*masterdis,lm,lmowner);
    elevector2.Size(1);   //length of circ with gaussinteg
    elevector3.Size(dim);   //centerdisp part of ele  

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
  masterdis->ClearState();

  //Communicate to 'assemble' length and center displacements
  Comm().SumAll(&mylengthcirc, &lengthcirc, 1);
  Comm().SumAll(&mycenterdisp[0], &centerdisp[0], dim);

  if (lengthcirc <= 1.0E-6)
    dserror("Zero interface length!");
    
  //calculating the final disp of the interface and summation over all time steps
  for (int i=0; i<dim ;i++)
  {
    centerdisp[i] = centerdisp[i] / lengthcirc;
    centerdisptotal_[i] +=centerdisp[i];
  }
  
  return centerdisp;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<int,LINALG::Matrix<3,1> > FSI::DirichletNeumannSlideale::CurrentMasterPos
(
    Teuchos::RCP<Epetra_Vector> idisptotal,
    Teuchos::RCP<Epetra_Vector> reddisp,
    DRT::Discretization& interfacedis,
    RCP<Epetra_Map> msfullelemap
)
{
  std::map<int,LINALG::Matrix<3,1> > currentpositions;

  // map with fully reduced master element distribution 
  map<int, RCP<DRT::Element> > masterreduelements;
  for (int eleind = 0; eleind<msfullelemap->NumMyElements(); eleind++)
  {
    DRT::Element* tmpele = interfacedis.lColElement(eleind);
    masterreduelements[tmpele->Id()]= rcp(tmpele,false);
    
    const int* n = tmpele->NodeIds();
    
    // fill currentpositions
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
  return currentpositions;
}

void FSI::DirichletNeumannSlideale::SlideProjection
(
    Teuchos::RCP<Epetra_Vector> islavestep,
    DRT::Discretization& interfacedis,
    vector<double> centerdisp,
    std::map<int, RCP<DRT::Element> > masterreduelements,
    std::map<int,LINALG::Matrix<3,1> > currentpositions
)
{
  INPAR::FSI::SlideALEProj aletype = 
      Teuchos::getIntegralValue<INPAR::FSI::SlideALEProj>(DRT::Problem::Instance()->FSIDynamicParams(),"SLIDEALEPROJ");
  const int dim = genprob.ndim;
  
  // Project slave nodes onto the master interface
  //init of search tree
  Teuchos::RCP<GEO::SearchTree> searchTree = rcp(new GEO::SearchTree(8));
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofEles(masterreduelements, currentpositions);

  if(dim==2)
    searchTree->initializeTreeSlideALE(rootBox, masterreduelements, GEO::TreeType(GEO::QUADTREE));
  else if(dim==3)
    searchTree->initializeTreeSlideALE(rootBox, masterreduelements, GEO::TreeType(GEO::OCTTREE));
  else dserror("wrong dimension");
  
  // translation + projection
  map<int, DRT::Node*>::const_iterator nodeiter;
  for (nodeiter = islaveslidnodes_.begin(); nodeiter != islaveslidnodes_.end(); ++nodeiter)
  {
    
    DRT::Node* node = nodeiter->second;
    vector<int> lids(dim);
    for(int p=0; p<dim; p++)
      //lids of gids of node
      lids[p] = (StructureFluidCouplingMortar().SlaveDofRowMap())->LID((MBFluidField().Discretization()->Dof(node))[p]);  

    // translate ale by centerdisp
    if (aletype==INPAR::FSI::ALEprojection_curr)
    {   
      //filling of islavestep for "Updated Lagrange Projektion"
      int err = islave_->SumIntoMyValues(dim, &centerdisp[0], &lids[0]);
      if (err == 1) dserror("error while adding values");
    }
    else
    {
      //filling of islavestep for "Total Lagrange Projektion"
      int err = islave_->ReplaceMyValues(dim, &centerdisptotal_[0], &lids[0]);
      if (err == 1) dserror("error while replacing values");
    }
    
    // current coord of ale node = ref coord + islave_
    LINALG::Matrix<3,1> alenodecurr;  

    for(int p=0; p<dim; p++)
      alenodecurr(p,0) =   node->X()[p]  + (*islave_)[(lids[p])]; 
      
    // coordinates to project to
    vector <double> finaldxyz(dim);

    //search for near elements next to the query point
    std::map<int,std::set<int> >  closeeles = 
        searchTree->searchElementsSlideALE(masterreduelements, currentpositions, alenodecurr);

    
    //search for the nearest point to project on
    if(dim == 2)
    {
      LINALG::Matrix<3,1> minDistCoords;
      GEO::nearestObjectInNode(imastergnodes_,  masterreduelements, currentpositions, 
          closeeles, alenodecurr, minDistCoords);
      finaldxyz[0] = minDistCoords(0,0) - alenodecurr(0,0); 
      finaldxyz[1] = minDistCoords(1,0) - alenodecurr(1,0);
    }
    else
    {
      LINALG::Matrix<3,1> minDistCoords;
      GEO::nearestObjectInNode(rcp(&interfacedis,false), masterreduelements, currentpositions, 
          closeeles, alenodecurr, minDistCoords);
      finaldxyz[0] = minDistCoords(0,0) - alenodecurr(0,0); 
      finaldxyz[1] = minDistCoords(1,0) - alenodecurr(1,0);
      finaldxyz[2] = minDistCoords(2,0) - alenodecurr(2,0);
      
    }
    
    //refill of islavestep with values of the projection
    int err = islavestep->ReplaceMyValues(dim, &finaldxyz[0], &lids[0]);
    if (err == 1) dserror("error while replacing values");

  }
}


#endif
