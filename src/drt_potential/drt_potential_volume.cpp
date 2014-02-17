/*!-------------------------------------------------------------------
\file drt_potential_volume.cpp

\brief  Class controlling volume stresses due to intermolecular interaction forces
        between interfaces of mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/

#include "drt_potential_volume.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_geometry/searchtree_geometry_service.H"
#include "../drt_geometry/intersection_service_templates.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_exporter.H"
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
#ifdef PARALLEL
      // TODO check for differetn potential cutoffs
  double cutoff =  (discretRCP_->GetCondition("Potential"))->GetDouble("cutOff");
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDisPar(*discretRCP_,cutoff);
#else
  const LINALG::Matrix<3,2> rootBox = GEO::getXAABBofDis(*discretRCP_);
#endif
  DRT::UTILS::CollectElementsByConditionLabel(*discretRCP_, elementsByLabel_,"Potential" );
  searchTree_->initializeTree(rootBox, elementsByLabel_, treetype_);
  InvertElementsPerLabel(elementsByLabel_,labelByElement_);

  const Epetra_Map* dofcolmap = discretRCP_->DofColMap();
  disp_col_ = LINALG::CreateVector(*dofcolmap, true);

  std::cout << "Volume potential constructor done" << std::endl;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::EvaluatePotential( Teuchos::ParameterList& p,
                                                    Teuchos::RCP<Epetra_Vector> disp,
                                                    Teuchos::RCP<Epetra_Vector> fint,
                                                    Teuchos::RCP<LINALG::SparseMatrix> stiff)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  disp_col_ = discret_.GetState("displacement");

  // update displacement for volume discretization discretRCP
  UpdateDisplacementsOfPotentialDiscretization(disp_col_);

  EvaluateVolumePotentialCondition(p,stiff,Teuchos::null,fint,Teuchos::null,Teuchos::null,"Potential");
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| evaluate potential conditions based on a Epetra_FecrsMatrix        |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::EvaluateVolumePotentialCondition(
    Teuchos::ParameterList&                 params,
    Teuchos::RCP<LINALG::SparseMatrix>      systemmatrix1,
    Teuchos::RCP<LINALG::SparseMatrix>      systemmatrix2,
    Teuchos::RCP<Epetra_Vector>             systemvector1,
    Teuchos::RCP<Epetra_Vector>             systemvector2,
    Teuchos::RCP<Epetra_Vector>             systemvector3,
    const std::string&                      condstring)
{
  if (!discret_.Filled()) dserror("FillComplete() was not called");
  if (!discret_.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time < 0.0) usetime = false;

  if(time < 0.0)
    std::cout <<  "no time curve set " << std::endl;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

   // get conditions
   std::vector<DRT::Condition*> potentialcond;
   discret_.GetCondition(condstring, potentialcond);

   int num_local_ele = 0;
   for(std::vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++ condIter)
   {
     std::map<int,Teuchos::RCP<DRT::Element> >& geom = (*condIter)->Geometry();
     num_local_ele += (int) geom.size();
   }

   int num_max_ele = 0;
   discret_.Comm().MaxAll(&num_local_ele,&num_max_ele,1 );
   int num_dummy_ele = num_max_ele - num_local_ele;

  //----------------------------------------------------------------------
  // loop through potential conditions and evaluate them
  //----------------------------------------------------------------------
   for(std::vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++ condIter)
   {
     std::map<int,Teuchos::RCP<DRT::Element> >& geom = (*condIter)->Geometry();
     // if (geom.empty()) dserror("evaluation of condition with empty geometry");
     // no check for empty geometry here since in parallel computations
     // can exist processors which do not own a portion of the elements belonging
     // to the condition geometry
     std::map<int,Teuchos::RCP<DRT::Element> >::iterator curr;

     // Evaluate Loadcurve if defined. Put current load factor in parameterlist
     const std::vector<int>*    curve  = (*condIter)->Get<std::vector<int> >("curve");
     int                   curvenum = -1;
     if (curve) curvenum = (*curve)[0];
     double                curvefac = 1.0;
     if (curvenum>=0 && usetime)
       curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);

     params.set("LoadCurveFactor",curvefac);

     params.set<Teuchos::RCP<DRT::Condition> >("condition", Teuchos::rcp(*condIter,false));

     // define element matrices and vectors
     Epetra_SerialDenseMatrix elematrix1;
     Epetra_SerialDenseMatrix elematrix2;
     Epetra_SerialDenseVector elevector1;
     Epetra_SerialDenseVector elevector2;
     Epetra_SerialDenseVector elevector3;

     for (curr=geom.begin(); curr!=geom.end(); ++curr)
     {
       // get element location vector and ownerships
       std::vector<int> lm;
       std::vector<int> lmrowowner;
       std::vector<int> lmstride;
       curr->second->LocationVector(discret_,lm,lmrowowner,lmstride);
       const int rowsize = lm.size();

       // tree search has to be called before evaluate for each element
       TreeSearch((curr->second).get(), (*condIter)->GetDouble("cutOff"), (*condIter)->GetInt("label"));

       // Reshape element matrices and vectors and init to zero in element->evaluate
       // call the element specific evaluate method
       int err = curr->second->Evaluate( params,discret_,lm, elematrix1,elematrix2,
           elevector1,elevector2,elevector3);

       if (err) dserror("error while evaluating elements");

       // specify lm row and lm col
       std::vector<int> lmrow;
       std::vector<int> lmcol;
       // only local values appeared
       if((int) lm.size() == rowsize)
       {
         lmrow.resize(lm.size());
         lmcol.resize(lm.size());
         lmrow = lm;
         lmcol = lm;
       }
       // non-local values appeared
       else if((int) lm.size() > rowsize)
       {
         for(int i = 0; i < rowsize; i++)
           lmrow.push_back(lm[i]);

         lmcol.resize(lm.size());
         lmcol = lm;
       }
       else
         dserror("lm is not properly filled");

       // assembly
       int eid = curr->second->Id();
       if (assemblemat1) systemmatrix1->FEAssemble(eid,elematrix1,lmrow,lmrowowner,lmcol);
       if (assemblemat2) systemmatrix2->FEAssemble(eid,elematrix2,lmrow,lmrowowner,lmcol);

       if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lmrow,lmrowowner);
       if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lmrow,lmrowowner);
       if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lmrow,lmrowowner);
     }
   } // for(std::vector<DRT::Condition*>::iterator condIter = potentialcond->begin() ; condIter != potentialcond->end(); ++ condIter)

   // tree search over dummy elements to allow for send_recv
   for(int i_dummy = 0; i_dummy < num_dummy_ele; i_dummy++)
     TreeSearchDummy();

  return;
} // end of EvaluateVolumePotentialCondition




/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| for volume elements                                                |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::StiffnessAndInternalForcesPotential(
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule3D&  gaussrule,
    Teuchos::ParameterList&         params,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       K_stiff,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize Lennard Jones potential constant variables
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition",Teuchos::null);

  // initialize time variables
  double curvefac = GetTimeCurveFactor(params);

  // tree search is performed in evaluate potential and
  // local and nonlocal pecs are filled

  // compute internal force and stiffness matrix
  ComputeFandK(element, gaussrule, lm, K_stiff, F_int, cond, curvefac);
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
    Teuchos::ParameterList&         params,
    std::vector<int>&               lm,
    Epetra_SerialDenseMatrix&       K_stiff,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize potential condition variables
  Teuchos::RCP<DRT::Condition> cond      = params.get<Teuchos::RCP<DRT::Condition> >("condition",Teuchos::null);

  // initialize time variables
  double curvefac = GetTimeCurveFactor(params);

  // compute internal force and stiffness matrix
  ComputeFandK(element, gaussrule, lm, K_stiff, F_int, cond, curvefac);
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
  for (int lid = 0; lid < discretRCP_->NumMyColNodes(); ++lid)
  {
    const DRT::Node* node = discretRCP_->lColNode(lid);
    std::vector<int> lm;
    lm.reserve(3);
    // extract global dof ids
    discretRCP_->Dof(node, lm);

    std::vector<double> mydisp(3);
    LINALG::Matrix<3,1> currpos;
    DRT::UTILS::ExtractMyValues(*idisp_solid,mydisp,lm);
    currpos(0) = node->X()[0] + mydisp[0];
    currpos(1) = node->X()[1] + mydisp[1];
    currpos(2) = node->X()[2] + mydisp[2];
    currentpositions_[node->Id()] = currpos;
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
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::TreeSearch(
    const DRT::Element*                                     element,
    const double                                            radius,
    const int                                               label)
{
   // get AABB of element and enlarge by cut off radius
   LINALG::Matrix<3,2> eleXAABB = elemXAABBList_[element->Id()];
   for(int dim = 0; dim < 3; dim++)
   {
     eleXAABB(dim,0) = eleXAABB(dim,0) - radius;
     eleXAABB(dim,1) = eleXAABB(dim,1) + radius;
   }

   // local ids and nonlocal pecs
   TreeSearchElement(eleXAABB, label, false);

   return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
| serial version of search method                                    |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::TreeSearchDummy()
{
   // create dummy ele XAABB
   LINALG::Matrix<3,2> eleXAABB(true);

   // local ids and nonlocal pecs
   TreeSearchElement(eleXAABB, -1, true);
   return;
}



/*-------------------------------------------------------------------*
| (protected)                                             umay  09/09|
| tree  search                                                       |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::TreeSearchElement(
    const LINALG::Matrix<3,2>&                              eleXAABB,
    const int                                   		        label,
    const bool                                              dummy)
{
  // local search
  localEleIds_.clear();
  if(label != -1) // dummy
    searchTree_->queryPotentialElements(elemXAABBList_, eleXAABB, localEleIds_, label);

  const int numprocs = discret_.Comm().NumProc();
  // if only one processor is running no parallel search has to be performed
  if(numprocs==1)
    return;

#ifdef PARALLEL

  const Epetra_MpiComm* comm = dynamic_cast<const Epetra_MpiComm*>(&(discret_.Comm()));

  std::set<int> pecIds;
  for(std::map<int, std::set<int> >::const_iterator labelIter = localEleIds_.begin(); labelIter != localEleIds_.end(); labelIter++)
    for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      pecIds.insert(*eleIter);

  // parallel search
  const int myrank = discret_.Comm().MyPID();
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
  else if(myrank==(numprocs-1))   //last processor
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
  std::vector<char> data_recv(1,'d');
  MPI_Status status;
  int label_send = label;
  int label_recv = 0;

  // send XAABB to other procs and collect pecs
  for(int i_proc = 0; i_proc < (numprocs-1); i_proc++)
  {
    // send XAABB
    int err = MPI_Sendrecv( &eleXAABB_send[0], 6, MPI_DOUBLE, forward,   tag_send,
                            &eleXAABB_recv[0], 6, MPI_DOUBLE, backward,  tag_recv,
                            comm->Comm(), &status);

    if (err != 0)
      dserror("mpi sendrecv error %d", err);
    eleXAABB_send = eleXAABB_recv;


    MPI_Sendrecv( &label_send, 1, MPI_INT, forward,   tag_send,
                  &label_recv, 1, MPI_INT, backward,  tag_recv,
                  comm->Comm(), &status);
    label_send = label_recv;


    // search for gids of potential elements on the current proc
    LINALG::Matrix<3,2> eleXAABB_LG;
    for(int dim = 0; dim < 3; dim++)
    {
      eleXAABB_LG(dim, 0) = eleXAABB_recv[dim];
      eleXAABB_LG(dim, 1) = eleXAABB_recv[dim+3];
    }

    std::map<int,std::set<int> > potEleGids;
    if(label_recv != -1)
      searchTree_->queryPotentialElements(elemXAABBList_, eleXAABB_LG, potEleGids, label_recv);

    std::vector<Teuchos::RCP<PotentialElementContainer> > vpec;

    // loop over labelIter->first
    for(std::map<int, std::set<int> >::const_iterator labelIter = potEleGids.begin(); labelIter != potEleGids.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
        DRT::Element* element = discretRCP_->gElement(*eleIter);
        std::vector<int> lmowner;
        std::vector<int> lm;
        std::vector<int> lmstride;
        element->LocationVector(*discretRCP_,lm,lmowner,lmstride);
        const double beta = GetAtomicDensity(element->Id(), "Potential", labelByElement_);

        Teuchos::RCP<PotentialElementContainer> pec = Teuchos::rcp( new PotentialElementContainer(
            element->Id(),
            element->Shape(),
            labelIter->first,
            beta,
            SpatialConfiguration(currentpositions_,element,prob_dim_),
            ReferenceConfiguration(element, prob_dim_),
            lm));

        //pec->Pack(data_send);
        vpec.push_back( pec );
        numEle_send++;
      }

    DRT::PackBuffer data;
    for ( std::vector<Teuchos::RCP<PotentialElementContainer> >::iterator i=vpec.begin(); i!=vpec.end(); ++i )
    {
      Teuchos::RCP<PotentialElementContainer> pec = *i;
      pec->Pack( data );
    }
    data.StartPacking();
    for ( std::vector<Teuchos::RCP<PotentialElementContainer> >::iterator i=vpec.begin(); i!=vpec.end(); ++i )
    {
      Teuchos::RCP<PotentialElementContainer> pec = *i;
      pec->Pack( data );
    }
    swap( data_send, data() );

    // send number of elements
    MPI_Sendrecv( &numEle_send, 1, MPI_INT, forward,  tag_send,
                  &numEle_recv, 1, MPI_INT, backward, tag_recv,
                  comm->Comm(), &status);
    numEle_send = numEle_recv;

    // let the next proc now about the data size he will receive
    int data_size_send =  (int) data_send.size();
    int data_size_recv = 0;
    MPI_Sendrecv(&data_size_send, 1, MPI_INT, forward,  tag_send,
                 &data_size_recv, 1, MPI_INT, backward, tag_recv,
                 comm->Comm(), &status);

    data_recv.clear();
    data_recv.resize(data_size_recv);

    // 2. Send data to the next processor and receive incoming data from previous processor
    MPI_Sendrecv(&data_send[0], data_size_send, MPI_CHAR, forward, tag_send,
                 &data_recv[0], data_size_recv, MPI_CHAR, backward,tag_recv,
                 comm->Comm(), &status);

    data_send = data_recv;
  } // loop over numprocs


  // unpack elements after returning to the starting processor
  if(label == -1 && numEle_send != 0)
    dserror("dummy element obtained potential elements");

  nonlocalPecs_.clear();
  std::vector<char>::size_type position = 0;
  for(int i_ele = 0; i_ele < numEle_send; ++i_ele)
  {
    Teuchos::RCP<PotentialElementContainer> pec = Teuchos::rcp( new PotentialElementContainer());
    pec->Unpack(data_send, position);
    // std::map<int, std::set<PotentialElementContainer>  set because some of the
    // elements are sends a few times since loop over col elements
    if(pecIds.find(pec->Id()) == pecIds.end()) // not yet store store pec
    {
      //pec->Print();
      nonlocalPecs_[pec->Body_label()].push_back(*pec);
      pecIds.insert(pec->Id());
    }
  }

  if(position != data_send.size())
    dserror("something is wrong with the data vector on proc %d", myrank);

  if(((label == -1) && (int) nonlocalPecs_.size() != 0) )//|| (label != label_recv))
    dserror("something is wrong with the dummy search %d", myrank);


#endif

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
   std::vector<int>&                                      lm,
   Epetra_SerialDenseMatrix&                              K_stiff,
   Epetra_SerialDenseVector&                              F_int,
   Teuchos::RCP<DRT::Condition>                            cond,
   const double                                           curvefac)
{

  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, localEleIds_, nonlocalPecs_, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_stiff.Shape(ndofrow, ndofcol);

  if(localEleIds_.empty() && nonlocalPecs_.empty())
   return;

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta   = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int                   numnode = actEle->NumNode();
    Epetra_SerialDenseVector    funct(numnode);
    Epetra_SerialDenseMatrix    deriv(3,numnode);
    LINALG::Matrix<3,1>         x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = localEleIds_.begin(); labelIter != localEleIds_.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

         // obtain current potential dofs
         std::vector<int> lmpot;
         std::vector<int> lmowner;
         std::vector<int> lmstride;
         element_pot->LocationVector(*discretRCP_,lmpot,lmowner,lmstride);

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
           const int                  numnode_pot = element_pot->NumNode();
           const double beta_pot =    GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);
           Epetra_SerialDenseVector   funct_pot(numnode_pot);
           Epetra_SerialDenseMatrix   deriv_pot(3,numnode_pot);
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
			     F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta_pot*potderiv1(dim)*fac_pot);

	         // computation of stiffness matrix (possibly with non-local values)
	         for (int inode = 0;inode < numnode; ++inode)
	           for(int dim = 0; dim < 3; dim++)
	           {
	             // k,ii
	             for (int jnode = 0; jnode < numnode; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_stiff(inode*numdof+dim, jnode*numdof+dim_pot) +=
	                   funct(inode)*beta*fac*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

	             // k,ij
	             for (int jnode = 0;jnode < numnode_pot; ++jnode)
	               for(int dim_pot = 0; dim_pot < 3; dim_pot++)
	                 K_stiff(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
	                    funct(inode)*beta*fac*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

	            }
	         } // loop over all gauss points of the potential element
      } // loop over all potential elements


    //----------------------------------------------------------------------
    // run over all influencing elements from different procs
    //----------------------------------------------------------------------
    for(std::map<int, std::vector<POTENTIAL::PotentialElementContainer> >::iterator labelIter = nonlocalPecs_.begin();
    labelIter != nonlocalPecs_.end(); labelIter++)
      for(std::vector<POTENTIAL::PotentialElementContainer>::iterator pecIter = (labelIter->second).begin();
      pecIter != (labelIter->second).end(); pecIter++)
      {

        // get atom density
        const double beta_pot = (*pecIter).Beta();
        // obtain current potential dofs
        std::vector<int> lmpot = (*pecIter).GetLm();

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
          const int                   numnode_pot = (*pecIter).NumNode();
          Epetra_SerialDenseVector    funct_pot(numnode_pot);
          Epetra_SerialDenseMatrix    deriv_pot(3,numnode_pot);
          LINALG::Matrix<3,1>         x_pot_gp(true);

          const double fac_pot = ComputeFactor(*pecIter, funct_pot, deriv_pot, intpoints_pot,
                                              gp_pot, x_pot_gp, curvefac);

          // evaluate Lennard Jones potential and its derivatives
          LINALG::Matrix<3,1>  potderiv1;
          LINALG::Matrix<3,3>  potderiv2;

          EvaluatePotentialfromCondition(cond, x_gp, x_pot_gp, potderiv1, potderiv2);
          //potderiv1.Print(std::cout);
          //potderiv2.Print(std::cout);

          const int numdof = 3;
          // computation of internal forces (possibly with non-local values)
          for (int inode = 0; inode < numnode; inode++)
            for(int dim = 0; dim < 3; dim++)
              F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta_pot*potderiv1(dim)*fac_pot);

          //F_int.Print(std::cout);
          // computation of stiffness matrix (possibly with non-local values)
          for (int inode = 0;inode < numnode; ++inode)
            for(int dim = 0; dim < 3; dim++)
            {
              // k,ii
              for (int jnode = 0; jnode < numnode; ++jnode)
                for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                  K_stiff(inode*numdof+dim, jnode*numdof+dim_pot) +=
                    funct(inode)*beta*fac*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

              // k,ij
              for (int jnode = 0;jnode < numnode_pot; ++jnode)
                for(int dim_pot = 0; dim_pot < 3; dim_pot++)
                  K_stiff(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                    funct(inode)*beta*fac*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

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
    std::vector<int>&                                           lm,
    Epetra_SerialDenseMatrix&                                   K_stiff,
    Epetra_SerialDenseVector&                                   F_int,
    Teuchos::RCP<DRT::Condition>                                 cond,
    const double                                                curvefac)
{
  // determine global row indices (lmrow) and global colum indices (lm)
  std::vector<int> lmrow = lm;
  CollectLmcol(discretRCP_, localEleIds_, nonlocalPecs_, lm);

  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();
  const int ndofcol    = lm.size();
  F_int.Size(ndofrow);
  K_stiff.Shape(ndofrow, ndofcol);

  if(localEleIds_.empty() && nonlocalPecs_.empty())
   return;

  // number of atoms (~0.2 nm) per volume area in reference configuration
  // here equal for all bodies in n/µm^2
  const double  beta    = cond->GetDouble("beta");
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor
    const int                 numnode = actEle->NumNode();
    Epetra_SerialDenseVector  funct(numnode);
    Epetra_SerialDenseMatrix  deriv(2,numnode);
    LINALG::Matrix<2,1>       x_gp(true);

    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::set<int> >::const_iterator labelIter = localEleIds_.begin(); labelIter != localEleIds_.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
        const DRT::Element* element_pot = discretRCP_->gElement(*eleIter);

        // obtain current potential dofs
        std::vector<int> lmpot;
        std::vector<int> lmowner;
        std::vector<int> lmstride;
        element_pot->LocationVector(*discretRCP_,lmpot,lmowner,lmstride);

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
          const int                  numnode_pot = element_pot->NumNode();
          const double beta_pot =    GetAtomicDensity(element_pot->Id(), "Potential", labelByElement_);;
          Epetra_SerialDenseVector   funct_pot(numnode_pot);
          Epetra_SerialDenseMatrix   deriv_pot(2,numnode_pot);
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
              F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta_pot*potderiv1(dim)*fac_pot);

          // computation of stiffness matrix (possibly with non-local values)
          for (int inode = 0;inode < numnode; ++inode)
            for(int dim = 0; dim < 2; dim++)
            {
              // k,ii
              for (int jnode = 0; jnode < numnode; ++jnode)
                for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                  K_stiff(inode*numdof+dim, jnode*numdof+dim_pot) +=
                    funct(inode)*beta*fac*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

              // k,ij
              for (int jnode = 0;jnode < numnode_pot; ++jnode)
                for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                  K_stiff(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                    funct(inode)*beta*fac*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);
            }
        } // loop over all gauss points of the potential element
      } // loop over all potential elements

    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    for(std::map<int, std::vector<PotentialElementContainer> >::iterator labelIter = nonlocalPecs_.begin(); labelIter != nonlocalPecs_.end(); labelIter++)
      for(std::vector<PotentialElementContainer>::iterator pecIter = (labelIter->second).begin(); pecIter != (labelIter->second).end(); pecIter++)
      {

        const double beta_pot = (*pecIter).Beta();
        // obtain current potential dofs
        std::vector<int> lmpot = (*pecIter).GetLm();

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
          const int                 numnode_pot = (*pecIter).NumNode();
          Epetra_SerialDenseVector  funct_pot(numnode_pot);
          Epetra_SerialDenseMatrix  deriv_pot(2,numnode_pot);
          LINALG::Matrix<2,1>       x_pot_gp(true);

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
              F_int[inode*numdof+dim] += funct(inode)*beta*fac*(beta_pot*potderiv1(dim)*fac_pot);

          // computation of stiffness matrix (possibly with non-local values)
          for (int inode = 0;inode < numnode; ++inode)
            for(int dim = 0; dim < 2; dim++)
            {
              // k,ii
              for (int jnode = 0; jnode < numnode; ++jnode)
                for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                  K_stiff(inode*numdof+dim, jnode*numdof+dim_pot) +=
                    funct(inode)*beta*fac*(beta_pot*potderiv2(dim,dim_pot)*funct(jnode)*fac_pot);

              // k,ij
              for (int jnode = 0;jnode < numnode_pot; ++jnode)
                for(int dim_pot = 0; dim_pot < 2; dim_pot++)
                  K_stiff(inode*numdof+dim, GetLocalIndex(lm,lmpot[jnode*numdof+dim_pot]) ) +=
                    funct(inode)*beta*fac*(beta_pot*(-1)*potderiv2(dim,dim_pot)*funct_pot(jnode)*fac_pot);

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
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
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

  Epetra_SerialDenseMatrix Jacobi(3,3);
  Epetra_SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);
  Epetra_SerialDenseMatrix x(numnode,3);
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
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
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

  Epetra_SerialDenseMatrix Jacobi(3,3);
  Epetra_SerialDenseMatrix X = pec.GetReferenceConfiguration();
  Epetra_SerialDenseMatrix x = pec.GetSpatialConfiguration();

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
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
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

  Epetra_SerialDenseMatrix Jacobi(2,2);
  Epetra_SerialDenseMatrix X(numnode,2);
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
    Epetra_SerialDenseVector&               funct,
    Epetra_SerialDenseMatrix&               deriv,
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

  Epetra_SerialDenseMatrix Jacobi(2,2);
  Epetra_SerialDenseMatrix X = pec.GetReferenceConfiguration();
  Epetra_SerialDenseMatrix x = pec.GetSpatialConfiguration();

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
    case DRT::Element::hex20:
    case DRT::Element::hex27:
      rule_pot = DRT::UTILS::intrule_hex_27point;
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



///////////////////////// test potential //////////////////////////////////////

/*-------------------------------------------------------------------*
| (public)                                                 umay 01/10|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/
void POTENTIAL::VolumePotential::TestEvaluatePotential(
  Teuchos::ParameterList&                      p,
  Teuchos::RCP<Epetra_Vector>          disp,
  Teuchos::RCP<Epetra_Vector>          fint,
  Teuchos::RCP<LINALG::SparseMatrix>   stiff,
  const double                        time,
  const int                           step)
{
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);

  EvaluateVolumePotentialCondition(p,stiff,Teuchos::null,fint,Teuchos::null,Teuchos::null,"Potential");

  Teuchos::RCP<const Epetra_Vector>        disp_col = discret_.GetState("displacement");
  // compute test results
  std::map<int, std::set<int> > empty_set;

  if( p.get<std::string>("solution type") == "Sphere" )
    computeTestVanDerWaalsSpheres(Teuchos::null, elementsByLabel_, empty_set, disp_col, fint,
                                time, step, p.get("vdw_radius", 0.0), p.get("n_offset", 0.0));
  else if( p.get<std::string>("solution type") == "Membrane" )
    computeTestVanDerWaalsMembranes(Teuchos::null, elementsByLabel_, empty_set, disp_col, fint,
                                    time, step, p.get("vdw_radius", 0.0), p.get("n_offset", 0.0),
                                    p.get("thickness", 0.0));
  else
    dserror("specify proper solution type");
  return;
}





