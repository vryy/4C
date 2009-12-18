/*!-------------------------------------------------------------------
\file drt_potential_volume.cpp

\brief  Class controlling volume stresses due to potential forces
        between interfaces of mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_potential_volume.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include <cstdlib>



/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 05/09|
 *-------------------------------------------------------------------*/
POTENTIAL::VolumePotential::VolumePotential(
    Teuchos::RCP<DRT::Discretization>     discretRCP,
    DRT::Discretization&                  discret,
    const GEO::TreeType&                  treetype):
    Potential(discretRCP, discret),
    treetype_(treetype)
{
  //treetype_ = treetype;
#ifdef PARALLEL
  double cutoff =  (discretRCP_->GetCondition("Potential"))->GetDouble("cutOff");
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDisPar(*discretRCP_,cutoff);
#else
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*discretRCP_);
#endif
  DRT::UTILS::CollectElementsByConditionLabel(*discretRCP_, elementsByLabel_,"Potential" );
  searchTree_->initializeTree(rootBox, elementsByLabel_, treetype_);
  
  const Epetra_Map* dofcolmap = discretRCP_->DofColMap();
  disp_col_ = LINALG::CreateVector(*dofcolmap, true);

  std::cout << "Volume potential constructor done" << endl;
}

    

/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::EvaluatePotential( ParameterList& p,
                                                    RefCountPtr<Epetra_Vector> disp,
                                                    RefCountPtr<Epetra_Vector> fint,
                                                    RefCountPtr<LINALG::SparseMatrix> stiff)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  disp_col_ = discret_.GetState("displacement");
  
  // update displacement for volume discretization discretRCP
  UpdateDisplacementsOfPotentialDiscretization(disp_col_);

  EvaluatePotentialCondition(p,stiff,null,fint,null,null,"Potential");
  return;
}



/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| for volume elements                                                |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule3D&  gaussrule,
    ParameterList&                  params,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize Lennard Jones potential constant variables
  RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null);
  const int     label     = cond->GetInt("label");
  const double  cutOff    = cond->GetDouble("cutOff");
  // initialize time variables
  double curvefac = GetTimeCurveFactor(params);

  std::map<int,std::set<int> > localEleIds;
  std::map<int,std::vector<PotentialElementContainer> > nonlocalPecs;

  // tree search
  TreeSearch(discretRCP_, elemXAABBList_, element, localEleIds, nonlocalPecs, cutOff, label); 

  // compute internal force and stiffness matrix
  // TODO if poteles empty don t do assembly
  ComputeFandK(element, gaussrule, localEleIds, nonlocalPecs, lm, K_surf, F_int, cond, label, curvefac);
  // cout << "stiffness stop" << endl;
  return;
}



/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    ParameterList&                  params,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize potential condition variables
  RefCountPtr<DRT::Condition> cond      = params.get<RefCountPtr<DRT::Condition> >("condition",null);
  const int                   label     = cond->GetInt("label");
  const double                cutOff    = cond->GetDouble("cutOff");
  // initialize time variables
  double curvefac = GetTimeCurveFactor(params);
  
  // tree search
  std::map<int,std::set<int> > localEleIds;
  std::map<int,std::vector<PotentialElementContainer> > nonlocalPecs;
  TreeSearch(discretRCP_, elemXAABBList_, element, localEleIds, nonlocalPecs, cutOff, label);  

  // compute internal force and stiffness matrix
  // TODO if poteles empty don t do assembly
  ComputeFandK(element, gaussrule, localEleIds, nonlocalPecs, lm, K_surf, F_int, cond, label, curvefac);
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| update displacements in redundant boundary discretization          |
| from solid discretization                                          |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::UpdateDisplacementsOfPotentialDiscretization(
    Teuchos::RCP<const Epetra_Vector>     idisp_solid
    )
{
  currentpositions_.clear();
  // run over volume discretization
  {
    for (int lid = 0; lid < discretRCP_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = discretRCP_->lColNode(lid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      discretRCP_->Dof(node, lm);
   
      vector<double> mydisp(3);
      LINALG::Matrix<3,1> currpos;
      DRT::UTILS::ExtractMyValues(*idisp_solid,mydisp,lm);
      currpos(0) = node->X()[0] + mydisp[0];
      currpos(1) = node->X()[1] + mydisp[1];
      currpos(2) = node->X()[2] + mydisp[2];
      currentpositions_[node->Id()] = currpos;
    }
  }

  // reinitialize search tree
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*discretRCP_, currentpositions_);
  searchTree_->initializeTree(rootBox, elementsByLabel_, treetype_);
  
  //build boxes around every element
  for (int lid = 0; lid < discretRCP_->NumMyColElements(); ++lid)
  {
  	const DRT::Element*  element = discretRCP_->lColElement(lid);
  	const LINALG::SerialDenseMatrix xyze(GEO::getCurrentNodalPositions(element,currentpositions_));
	  GEO::EleGeoType eleGeoType(GEO::HIGHERORDER);
	  GEO::checkRoughGeoType(element, xyze, eleGeoType);
	  elemXAABBList_[element->Id()]=GEO::computeFastXAABB(element->Shape(), xyze, eleGeoType);	
  }
  return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::TreeSearch(
    const Teuchos::RCP<DRT::Discretization>     		        potentialdis,
    std::map<int,LINALG::Matrix<3,2> >&         		        elemXAABBList,
    const DRT::Element*                         		        element,
    std::map<int,std::set<int> >&               		        localEleIds,
    std::map<int,std::vector<PotentialElementContainer> >&  nonlocalPecs,
    const double                                		        radius,
    const int                                   		        label)
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  
  // get AABB of element
  LINALG::Matrix<3,2> eleXAABB = elemXAABBList[element->Id()];
  for(int dim = 0; dim < 3; dim++)
  {
    eleXAABB(dim,0) = eleXAABB(dim,0) - radius;
    eleXAABB(dim,1) = eleXAABB(dim,1) + radius;
  }
  
  // local search
  searchTree_->queryPotentialElements(elemXAABBList, eleXAABB, localEleIds, label);
  
  const int numprocs = comm.NumProc();
  // if only one processor is running no parallel search has to be performed
  if(numprocs==1)
    return;
  
#ifdef PARALLEL
  
  // parallel search
  const int myrank = comm.MyPID();
  int forward = 0;
  int backward = 0;
  int tag_send = myrank+1;
  int tag_recv = myrank;
  
  //where is my destination and source
  if (myrank==0)              //first processor
  {
    forward=myrank+1;
    backward=numprocs-1;
  }
  else if (myrank==(numprocs-1))   //last processor
  {
    forward=0;
    backward=myrank-1;
    tag_send = 0;
  }
  else if (myrank!=0 && myrank!=(numprocs-1))  //neither first nor last processor
  {
    forward=myrank+1;
    backward=myrank-1;
  }
  else
    dserror("something is wrong with myrank and numprocs");
  
  
  // enlarge box by cut off radius
  std::vector<double> eleXAABB_send(6,0.0);
  std::vector<double> eleXAABB_recv(6,0.0);
  for(int dim = 0; dim < 3; dim++)
  {
    eleXAABB_send[dim] = eleXAABB(dim,0);
    eleXAABB_send[dim+3] = eleXAABB(dim,1);
  }
   
  // initialize with dummy character in case that no elements are found on one proc
  int numEle_send = 0;
  int numEle_recv = 0;
  std::vector<char> data_send(1,'d');
  std::vector<char> data_recv;
  MPI_Status status;

  for(int i_proc = 0; i_proc < numprocs; i_proc++)
  {  
    // send XAABB 
    int err = MPI_Sendrecv( &eleXAABB_send[0], 6, MPI_DOUBLE, forward,   tag_send, 
                            &eleXAABB_recv[0], 6, MPI_DOUBLE, backward,  tag_recv, 
                            MPI_COMM_WORLD, &status);

    if (err != 0)
      dserror("mpi sendrecv error %d", err);

    eleXAABB_send = eleXAABB_recv;

    // search for gids of potential elements on the current proc
    std::map<int,std::set<int> > potEleGids;
    LINALG::Matrix<3,2> eleXAABB_LG;
    for(int dim = 0; dim < 3; dim++)
    {
      eleXAABB_LG(dim, 0) = eleXAABB_recv[dim];
      eleXAABB_LG(dim, 1) = eleXAABB_recv[dim+3];
      cout << "eleXAABB_LG(dim, 0) = " << eleXAABB_LG(dim, 0) << "myrank = " << myrank <<  endl;
      cout << "eleXAABB_LG(dim, 1) = " << eleXAABB_LG(dim, 1) << "myrank = " << myrank <<  endl;
    }
    searchTree_->queryPotentialElements(elemXAABBList, eleXAABB_LG, potEleGids, label);
    
    // loop over labelIter->first
    for(std::map<int, std::set<int> >::const_iterator labelIter = potEleGids.begin(); labelIter != potEleGids.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
        cout << "eleId" <<  *eleIter << endl;
        DRT::Element* element = discretRCP_->gElement(*eleIter);
        vector<int> lmowner;
        vector<int> lm;
        element->LocationVector(*discretRCP_,lm,lmowner);

        RCP<PotentialElementContainer> pec = rcp( new PotentialElementContainer(
            element->Id(),
            element->Owner(),
            element->Shape(),
            labelIter->first,
            SpatialConfiguration(currentpositions_,element,prob_dim_),
            ReferenceConfiguration(element, prob_dim_),
            lm));

        pec->Pack(data_send);
        numEle_send++;
      }

    potEleGids.clear();
    
    // send number of elements
    MPI_Sendrecv( &numEle_send, 1, MPI_INT, forward,  tag_send, 
                  &numEle_recv, 1, MPI_INT, backward, tag_recv, 
                  MPI_COMM_WORLD, &status);
    numEle_send = numEle_recv;

    // let the next proc now about the data size he will receive
    int data_size_send =  (int) data_send.size();
    int data_size_recv = 0;
    MPI_Sendrecv(&data_size_send, 1, MPI_INT, forward,  tag_send, 
                 &data_size_recv, 1, MPI_INT, backward, tag_recv, 
                 MPI_COMM_WORLD, &status);
    data_recv.resize(data_size_recv);
    
    // 2. Send data to the next processor and receive incoming data from previous processor
    MPI_Sendrecv(&data_send[0], data_size_send, MPI_CHAR, forward, tag_send, 
                 &data_recv[0], data_size_recv, MPI_CHAR, backward,tag_recv, 
                 MPI_COMM_WORLD, &status);
    data_send = data_recv;
    data_recv.clear();
  }

  // unpack elements after returning to the starting processor
  // store in potentialElements 
  std::set<int> pecIds;
  int position = 1;  // jump over dummy
  for(int i_ele = 0; i_ele < numEle_send; i_ele++)
  {
    PotentialElementContainer pec;
    pec.Unpack(data_send, position);
    // std::map<int, std::set<PotentialElementContainer>  set because some of the
    // elements are sends a few times since loop over col elements
    if(pecIds.find(pec.Id()) == pecIds.end()) // not yet store store pec
    {
      nonlocalPecs[pec.Body_label()].push_back(pec);
      pecIds.insert(pec.Id());
    }
  }

  if(position != (int) data_send.size())
    dserror("something is wrong with the data vector");
#endif
  
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix (volume)        |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::ComputeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule3D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   vector<int>&                     lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   RefCountPtr<DRT::Condition>      cond,
   const int                        label,
   const double                     curvefac)
{

  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, potElements, lm); //MINE potentialdis_
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta");	
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(3,numnode);
    LINALG::Matrix<3,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);

    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule3D rule_pot = DRT::UTILS::intrule3D_undefined;
         GetGaussRule3D(element_pot->Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints3D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(3,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
         
           const int numdof = 3;
	
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 3; dim++)
			     F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
	
	         // computation of stiffness matrix (possibly with non-local values)
	         for (int inode = 0;inode < numnode; ++inode)
	           for(int dim = 0; dim < 3; dim++)
	           {
	             // k,ii
	             for (int jnode = 0; jnode < numnode; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
	                   funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
	
	             // k,ij
	             for (int jnode = 0;jnode < numnode_pot; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
	                    funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);
	
	            }
	         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix (volume)        |
| for parallel search
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::ComputeFandK(
   const DRT::Element*                                    actEle,
   const DRT::UTILS::GaussRule3D&                         gaussrule,
   std::map<int,std::set<int> >&                          localEleIds,
   std::map<int,std::vector<PotentialElementContainer> >& nonlocalPecs,
   vector<int>&                                           lm,
   Epetra_SerialDenseMatrix&                              K_surf,
   Epetra_SerialDenseVector&                              F_int,
   RefCountPtr<DRT::Condition>                            cond,
   const int                                              label,
   const double                                           curvefac)
{

  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, localEleIds, nonlocalPecs, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta");	
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(3,numnode);
    LINALG::Matrix<3,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);

    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = localEleIds.begin(); labelIter != localEleIds.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule3D rule_pot = DRT::UTILS::intrule3D_undefined;
         GetGaussRule3D(element_pot->Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints3D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(3,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
         
           const int numdof = 3;
	
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 3; dim++)
			     F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
	
	         // computation of stiffness matrix (possibly with non-local values)
	         for (int inode = 0;inode < numnode; ++inode)
	           for(int dim = 0; dim < 3; dim++)
	           {
	             // k,ii
	             for (int jnode = 0; jnode < numnode; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
	                   funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
	
	             // k,ij
	             for (int jnode = 0;jnode < numnode_pot; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
	                    funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);
	
	            }
	         } // loop over all gauss points of the potential element
      } // loop over all potential elements

      //----------------------------------------------------------------------
    // run over all influencing elements from different procs
    //----------------------------------------------------------------------
    for(std::map<int, std::vector<POTENTIAL::PotentialElementContainer> >::iterator labelIter = nonlocalPecs.begin(); 
          labelIter != nonlocalPecs.end(); labelIter++)
      for(std::vector<POTENTIAL::PotentialElementContainer>::iterator pecIter = (labelIter->second).begin(); 
          pecIter != (labelIter->second).end(); pecIter++)
      {
         // obtain current potential dofs
         vector<int> lmpot = (*pecIter).GetLm();
       
         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule3D rule_pot = DRT::UTILS::intrule3D_undefined;
         GetGaussRule3D((*pecIter).Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints3D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = (*pecIter).NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(3,numnode_pot);
           LINALG::Matrix<3,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(*pecIter, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<3,1>  potderiv1;
           LINALG::Matrix<3,3>  potderiv2;

           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
         
           const int numdof = 3;
	
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 3; dim++)
			     F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
	
	         // computation of stiffness matrix (possibly with non-local values)
	         for (int inode = 0;inode < numnode; ++inode)
	           for(int dim = 0; dim < 3; dim++)
	           {
	             // k,ii
	             for (int jnode = 0; jnode < numnode; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
	                   funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
	
	             // k,ij
	             for (int jnode = 0;jnode < numnode_pot; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
	                    funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);
	
	            }
	         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}


/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::ComputeFandK(
   const DRT::Element*                      actEle,
   const DRT::UTILS::GaussRule2D&           gaussrule,
   std::map<int,std::set<int> >&            potElements,
   vector<int>&                             lm,
   Epetra_SerialDenseMatrix&                K_surf,
   Epetra_SerialDenseVector&                F_int,
   RefCountPtr<DRT::Condition>              cond,
   const int                                label,
   const double                             curvefac)
{
  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(2,numnode);
    LINALG::Matrix<2,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
         GetGaussRule2D(element_pot->Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(2,numnode_pot);
           LINALG::Matrix<2,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<2,1>  potderiv1;
           LINALG::Matrix<2,2>  potderiv2;
           
           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
            
           const int numdof = 2;
  
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 2; dim++)
               F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
  
           // computation of stiffness matrix (possibly with non-local values)
           for (int inode = 0;inode < numnode; ++inode)
             for(int dim = 0; dim < 2; dim++)
             {
               // k,ii
               for (int jnode = 0; jnode < numnode; ++jnode)
                 for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                   K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                     funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
  
               // k,ij
               for (int jnode = 0;jnode < numnode_pot; ++jnode)
                 for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                   K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                      funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

            }
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}




/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix (surface)       |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::ComputeFandK(
    const DRT::Element*                                         actEle,
    const DRT::UTILS::GaussRule2D&                              gaussrule,
    std::map<int,std::set<int> >&                               localEleIds,
    std::map<int,std::vector<PotentialElementContainer> >&      nonlocalPecs,
    vector<int>&                                                lm,
    Epetra_SerialDenseMatrix&                                   K_surf,
    Epetra_SerialDenseVector&                                   F_int,
    RefCountPtr<DRT::Condition>                                 cond,
    const int                                                   label,
    const double                                                curvefac)
{
  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, localEleIds, nonlocalPecs, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(2,numnode);
    LINALG::Matrix<2,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = localEleIds.begin(); labelIter != localEleIds.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner);

         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
         GetGaussRule2D(element_pot->Shape(), rule_pot);
         const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(2,numnode_pot);
           LINALG::Matrix<2,1>        x_pot_gp(true);

           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot,
                                                gp_pot, x_pot_gp, curvefac);

           // evaluate Lennard Jones potential and its derivatives
           LINALG::Matrix<2,1>  potderiv1;
           LINALG::Matrix<2,2>  potderiv2;
           
           EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
            
           const int numdof = 2;
  
           // computation of internal forces (possibly with non-local values)
           for (int inode = 0; inode < numnode; inode++)
             for(int dim = 0; dim < 2; dim++)
               F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);
  
           // computation of stiffness matrix (possibly with non-local values)
           for (int inode = 0;inode < numnode; ++inode)
             for(int dim = 0; dim < 2; dim++)
             {
               // k,ii
               for (int jnode = 0; jnode < numnode; ++jnode)
                 for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                   K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                     funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);
  
               // k,ij
               for (int jnode = 0;jnode < numnode_pot; ++jnode)
                 for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                   K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                      funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

            }
         } // loop over all gauss points of the potential element
      } // loop over all potential elements
    
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::vector<PotentialElementContainer> >::iterator labelIter = nonlocalPecs.begin(); labelIter != nonlocalPecs.end(); labelIter++)
      for(std::vector<PotentialElementContainer>::iterator pecIter = (labelIter->second).begin(); pecIter != (labelIter->second).end(); pecIter++)
      {
        // obtain current potential dofs
        vector<int> lmpot = (*pecIter).GetLm();
      
        // obtain Gaussrule and integration points
        DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
        GetGaussRule2D((*pecIter).Shape(), rule_pot);
        const DRT::UTILS::IntegrationPoints2D intpoints_pot(rule_pot);
        //----------------------------------------------------------------------
        // run over all gauss points of a influencing element
        //----------------------------------------------------------------------
        for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
        {
          // compute func, deriv, x_gp and factor
          const int numnode_pot = (*pecIter).NumNode();
          LINALG::SerialDenseVector  funct_pot(numnode_pot);
          LINALG::SerialDenseMatrix  deriv_pot(2,numnode_pot);
          LINALG::Matrix<2,1>        x_pot_gp(true);

          const double fac_pot = ComputeFactor(*pecIter, funct_pot, deriv_pot, intpoints_pot,
              gp_pot, x_pot_gp, curvefac);

          // evaluate Lennard Jones potential and its derivatives
          LINALG::Matrix<2,1>  potderiv1;
          LINALG::Matrix<2,2>  potderiv2;

          EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);

          const int numdof = 2;

          // computation of internal forces (possibly with non-local values)
          for (int inode = 0; inode < numnode; inode++)
            for(int dim = 0; dim < 2; dim++)
              F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta*potderiv1(dim)*fac_pot);

          // computation of stiffness matrix (possibly with non-local values)
          for (int inode = 0;inode < numnode; ++inode)
            for(int dim = 0; dim < 2; dim++)
            {
              // k,ii
              for (int jnode = 0; jnode < numnode; ++jnode)
                for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                  K_surf(inode*numdof+dim, jnode*numdof+dim_pot) +=
                    funct(inode)*beta*fac*(beta*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

              // k,ij
              for (int jnode = 0;jnode < numnode_pot; ++jnode)
                for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                  K_surf(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                    funct(inode)*beta*fac*(beta*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

            }
        } // loop over all gauss points of the potential element
      } // loop over all potential elements
  } // loop over all gauss points of the actual element

  return;
}


/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for volume elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    const DRT::Element*                     element,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints3D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];	
  const double e1 = intpoints.qxg[gp][1];
  const double e2 = intpoints.qxg[gp][2];
  
  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_3D(funct,e0,e1,e2,element->Shape());
  DRT::UTILS::shape_function_3D_deriv1(deriv,e0,e1,e2,element->Shape());

  
  LINALG::SerialDenseMatrix Jacobi(3,3);
  LINALG::SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);		
  LINALG::SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(currentpositions_,element,x,3);	
  Jacobi.Multiply('N','N',1.0,deriv,X,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = Jacobi(0,0)*(   Jacobi(1,1)*Jacobi(2,2) 
                                    - Jacobi(2,1)*Jacobi(1,2)  ) -
                      Jacobi(0,1)*(   Jacobi(1,0)*Jacobi(2,2) 
                                    - Jacobi(2,0)*Jacobi(1,2)  ) +
                      Jacobi(0,2)*(   Jacobi(1,0)*Jacobi(2,1) 
                                    - Jacobi(2,0)*Jacobi(1,1)  );                 		
                            		
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for volume elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    PotentialElementContainer&	            pec,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints3D&  intpoints,
    const int                               gp,
    LINALG::Matrix<3,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = pec.NumNode();
  const double e0 = intpoints.qxg[gp][0];	
  const double e1 = intpoints.qxg[gp][1];
  const double e2 = intpoints.qxg[gp][2];
  
  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_3D(funct,e0,e1,e2,pec.Shape());
  DRT::UTILS::shape_function_3D_deriv1(deriv,e0,e1,e2,pec.Shape());

  LINALG::SerialDenseMatrix Jacobi(3,3);
  LINALG::SerialDenseMatrix X = pec.GetReferenceConfiguration();
  LINALG::SerialDenseMatrix x = pec.GetSpatialConfiguration();
  Jacobi.Multiply('N','N',1.0,deriv,X,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = Jacobi(0,0)*(   Jacobi(1,1)*Jacobi(2,2) 
                                    - Jacobi(2,1)*Jacobi(1,2)  ) -
                      Jacobi(0,1)*(   Jacobi(1,0)*Jacobi(2,2) 
                                    - Jacobi(2,0)*Jacobi(1,2)  ) +
                      Jacobi(0,2)*(   Jacobi(1,0)*Jacobi(2,1) 
                                    - Jacobi(2,0)*Jacobi(1,1)  );                 		
  
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}


/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for surface elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    const DRT::Element*                     element,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    LINALG::Matrix<2,1>&                    x_gp,
    const double                            curve_fac)
{

  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];

  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());

  LINALG::SerialDenseMatrix Jacobi(2,2);
  LINALG::SerialDenseMatrix X(numnode,2);
  ReferenceConfiguration(element,X,2);
  LINALG::SerialDenseMatrix x(numnode,2);
  SpatialConfiguration(currentpositions_,element,x,2);
  Jacobi.Multiply('N','N',1.0,deriv,X,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = Jacobi(0,0)*Jacobi(1,1) - Jacobi(1,0)*Jacobi(0,1);
  		  				               		
                              		
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 2; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}




/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp, for surface elements  u.may 08/08|
 *----------------------------------------------------------------------*/
double POTENTIAL::VolumePotential::ComputeFactor(
    PotentialElementContainer&              pec,
    LINALG::SerialDenseVector&              funct,
    LINALG::SerialDenseMatrix&              deriv,
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    LINALG::Matrix<2,1>&                    x_gp,
    const double                            curve_fac)
{
  const int numnode = pec.NumNode();
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];

  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,pec.Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,pec.Shape());

  LINALG::SerialDenseMatrix Jacobi(2,2);
  LINALG::SerialDenseMatrix X = pec.GetReferenceConfiguration();
  LINALG::SerialDenseMatrix x = pec.GetSpatialConfiguration();
  Jacobi.Multiply('N','N',1.0,deriv,X,0.0);

  // detA maps the reference configuration to the parameter space domain
  const double detA = Jacobi(0,0)*Jacobi(1,1) - Jacobi(1,0)*Jacobi(0,1);                                                                 
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 2; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);

  return factor;
}



/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 09/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::GetGaussRule3D(
    const DRT::Element::DiscretizationType& distype,
    DRT::UTILS::GaussRule3D&                rule_pot)
{
  // obtain Gaussrule and integration points
  switch (distype)
  {
    case DRT::Element::hex8:
      rule_pot = DRT::UTILS::intrule_hex_8point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return;
}



/*----------------------------------------------------------------------*
 |  get Gauss rule and integraion points                     u.may 04/09|
 *----------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::GetGaussRule2D(
    const DRT::Element::DiscretizationType& distype,
    DRT::UTILS::GaussRule2D&                rule_pot)
{
  // obtain Gaussrule and integration points
  switch (distype)
  {
    case DRT::Element::quad4:
      rule_pot = DRT::UTILS::intrule_quad_4point;
      break;
    case DRT::Element::quad8: case DRT::Element::quad9:
      rule_pot = DRT::UTILS::intrule_quad_9point;
      break;
    case DRT::Element::tri3:
      rule_pot = DRT::UTILS::intrule_tri_3point;
      break;
    case DRT::Element::tri6:
      rule_pot = DRT::UTILS::intrule_tri_6point;
      break;
    default:
      dserror("unknown number of nodes for gaussrule initialization");
  }
  return;
}




#endif




