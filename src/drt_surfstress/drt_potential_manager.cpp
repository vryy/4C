 /*!-------------------------------------------------------------------
\file drt_potential_manager.cpp

\brief Class controlling surface stresses due to potential forces
between interfaces


<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*--------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <blitz/array.h>
#include "drt_potential_manager.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_timecurve.H"
#include <cstdlib>


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
UTILS::PotentialManager::PotentialManager(
    Teuchos::RCP<DRT::Discretization> discretRCP,
    DRT::Discretization& discret):
    discretRCP_(discretRCP),
    discret_(discret)
{
  // create discretization from potential boundary elements and distribute them
  // on every processor
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("Potential");
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(discretRCP_, "Potential", "PotBoundary", "BELE3", conditions_to_copy);
  dsassert(boundarydis_->NumGlobalNodes() > 0, "empty discretization detected. Potential conditions applied?");

  // set new dof set
  RCP<UTILS::PotentialDofSet> pdofset = rcp(new UTILS::PotentialDofSet(discretRCP_));
  (*boundarydis_).ReplaceDofSet(pdofset);
  (*boundarydis_).FillComplete(false, false, false);
 
  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *(boundarydis_->NodeRowMap());
  const Epetra_Map elemrowmap = *(boundarydis_->ElementRowMap());
  
  // put all boundary nodes and elements onto all processors
  const Teuchos::RCP<const Epetra_Map> nodecolmap = LINALG::AllreduceEMap(*(boundarydis_->NodeRowMap()));
  const Teuchos::RCP<const Epetra_Map> elemcolmap = LINALG::AllreduceEMap(*(boundarydis_->ElementRowMap()));
  
  // redistribute nodes and elements to column (ghost) map
  boundarydis_->ExportColumnNodes(*nodecolmap);
  boundarydis_->ExportColumnElements(*elemcolmap);
 
  // Now we are done. :)
  // fill complete calls the assgin degrees of fredom method
  // provided by PotentialDofSet class
  const int err = boundarydis_->FillComplete();
  if (err) dserror("FillComplete() returned err=%d",err);
  
  // split dof vector of soliddiscretization into a dof vector of potential boundary condition and 
  // remaining dofs
  DRT::UTILS::SetupNDimExtractor(*discretRCP_ ,"Potential", potboundary_);
  
  // create potential surface dof row map using the solid parallel distribution
  // for all potential surface elements belonging to the solid discretization on a single 
  // processor
  const Teuchos::RCP< const Epetra_Map > potsurface_condmap_onOneProc = potboundary_.CondMap();
  if(! potsurface_condmap_onOneProc->UniqueGIDs() ) dserror("cond_map not unique");
  if(! (potsurface_condmap_onOneProc->PointSameAs(*(*boundarydis_).DofRowMap())))
    dserror("maps are not point equal");
  if(! (potsurface_condmap_onOneProc->SameAs(*(*boundarydis_).DofRowMap())))
    dserror("maps are not equal");
  // create disp vector for boundary discretization corresponding to the 
  // potential condition elements stored on one proc
  idisp_onproc_    = LINALG::CreateVector(*potsurface_condmap_onOneProc,true);
  
  // create potential surface dof col map based on the boundayr dis
  // that is distributed redundantly on all processors
  const Epetra_Map* potsurface_dofcolmap_total = boundarydis_->DofColMap();
  // create disp vector for boundary discretization redundant on all procs
  idisp_total_    = LINALG::CreateVector(*potsurface_dofcolmap_total,true);
  
  // create importer
  importer_ = rcp(new Epetra_Import(idisp_total_->Map(),idisp_onproc_->Map()));

  // set up tree
  const BlitzMat3x2 rootBox = GEO::getXAABBofDis(*boundarydis_);
  octTree_      = rcp(new GEO::SearchTree(8));
  CollectElementsByPotentialCouplingLabel(elementsByLabel_);
  octTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  
  // std::cout << "Potential manager constructor done" << endl;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/

void UTILS::PotentialManager::EvaluatePotential(  ParameterList& p,
                                                  RefCountPtr<Epetra_Vector> disp,
                                                  RefCountPtr<Epetra_Vector> fint,
                                                  RefCountPtr<LINALG::SparseMatrix> stiff)
{   
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  UpdateDisplacementsOfBoundaryDiscretization(disp);
  
  EvaluatePotentialCondition(p,stiff,null,fint,null,null,"Potential");
  
  return;
}



/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::StiffnessAndInternalForcesLJPotential( 
    const DRT::Element*             element,
    const DRT::UTILS::GaussRule2D&  gaussrule,
    ParameterList&                  params,
    vector<int>&                    lm,
    Epetra_SerialDenseMatrix&       K_surf,
    Epetra_SerialDenseVector&       F_int)
{
  // initialize Lennard Jones potential constant variables
  RefCountPtr<DRT::Condition> cond = params.get<RefCountPtr<DRT::Condition> >("condition",null); 

  // find nodal ids influencing a given element
  const int     label     = cond->Getint("label");
  const double  cutOff    = cond->GetDouble("cutOff");
  std::map<int,std::set<int> > potentialElementIds;
  for(int i = 0; i < DRT::UTILS::getNumberOfElementCornerNodes(element->Shape()); i++)
  {
    const DRT::Node* node = element->Nodes()[i];
    const BlitzVec3 x_node = currentboundarypositions_.find(node->Id())->second;
    // octtree search
    octTreeSearchElementsInCutOffRadius(x_node, potentialElementIds, cutOff, label);
    // serial search
    // searchElementsInCutOffRadius(eleId, potentialElementIds, cutOff);
  } 
  
  // compute internal force and stiffness matrix
  // initialize potential variables
  const double  depth     = cond->GetDouble("depth");
  const double  rootDist  = cond->GetDouble("rootDist");
  
  // initialize time variables
  const int    curvenum = cond->Getint("curve");
  const double time     = params.get<double>("total time",-1.0);
  //const double dt     = params.get<double>("delta time",0.0);
  const double t_end    = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).end();
  double curvefac       = 1.0;
  // apply potential forces gradually
  if (time <= t_end)
    curvefac      = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
  
  // compute internal force and stiffness matrix
  // TODO if poteles empty don t do assembly
  computeFandK(element, gaussrule, potentialElementIds, lm, K_surf, F_int, label, depth, rootDist, curvefac);
  // cout << "stiffness stop" << endl;
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| update displacements in redundant boundary discretization          |
| from solid discretization                                          |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::UpdateDisplacementsOfBoundaryDiscretization(
    Teuchos::RCP<Epetra_Vector>     idisp_solid
    )
{
  // extract displacements associated with potential elements from solid discretization
  idisp_onproc_ = potboundary_.ExtractCondVector(idisp_solid);
  idisp_total_->Scale(0.0);
 
 // LINALG::Export(*idisp_onproc_,*idisp_total_);
  
  // import 
  int err = idisp_total_->Import((*idisp_onproc_), (*importer_),Insert);
  if(err) dserror("Import using importer returned err=%d",err);

  currentboundarypositions_.clear();
  {
    for (int lid = 0; lid < boundarydis_->NumMyColNodes(); ++lid)
    {
      const DRT::Node* node = boundarydis_->lColNode(lid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      boundarydis_->Dof(node, lm);
      vector<double> mydisp(3);
      BlitzVec3 currpos(3);
      DRT::UTILS::ExtractMyValues(*idisp_total_,mydisp,lm);
      currpos(0) = node->X()[0] + mydisp[0];
      currpos(1) = node->X()[1] + mydisp[1];
      currpos(2) = node->X()[2] + mydisp[2];
      currentboundarypositions_[node->Id()] = currpos;
    }
  }
  
  // reinitialize search tree
  const BlitzMat3x2 rootBox = GEO::getXAABBofDis(*boundarydis_, currentboundarypositions_);
  octTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
}



/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::searchElementsInCutOffRadius(
    const BlitzVec3&                  point,
    std::map<int,std::set<int> >&     potentialElementIds,
    const double                      radius)
{

const double tol =  1e-7;

// only for testing no labels are distinguished
  for (int lid = 0; lid < boundarydis_->NumMyColNodes(); ++lid)
  {
    // compute distance between point and node
    const BlitzVec3 currpos = currentboundarypositions_[lid];
    double distance = 0;
    for(int k = 0; k < 3; k++)
    {
      const double comp_diff = currpos(k) - point(k);
      const double diff_square = comp_diff*comp_diff;
      distance += diff_square; 
    }
    distance = sqrt(distance);
    
    // if node lies within cutoff radius -> collect element id in set
    if(fabs(distance) - fabs(radius) < tol)
    {
      const DRT::Node* node = boundarydis_->lColNode(lid);
      const DRT::Element *const* elements = node->Elements();
      for(int k = 0; k <  node->NumElement(); k++)
        potentialElementIds[0].insert(elements[k]->Id());
    }
  }
}



/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| serial version of search method                                    |
| method runs over all nodes of the boundary discretization and      |
| checks, if a node lies within a given cut off radius               |
| in this case the element ids of adjacent elements are stored       |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::octTreeSearchElementsInCutOffRadius(
    const BlitzVec3&                  point,
    std::map<int,std::set<int> >&     potentialElementIds,
    const double                      radius,
    const int                         label)
{
  potentialElementIds = 
    octTree_->searchElementsInRadius(*boundarydis_, currentboundarypositions_, point, radius, label);
  
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| collects all emement ids that belong to certain label              |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::CollectElementsByPotentialCouplingLabel(
    std::map<int,std::set<int> >&        elementsByLabel)
{
  // Reset
  elementsByLabel.clear();
  
  // get condition
  vector< DRT::Condition * >  potConditions;
  boundarydis_->GetCondition("Potential", potConditions);
  
  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = potConditions.begin(); conditer!=potConditions.end(); ++conditer)
  {
    DRT::Condition* potCondition = *conditer;
    const int label = potCondition->Getint("label");
    // here the node have to be visitied because the element geometry map doesn t 
    // have ghosted elements
    const vector<int> geometryMap = *potCondition->Nodes();
    vector<int>::const_iterator iterNode;
    for(iterNode = geometryMap.begin(); iterNode != geometryMap.end(); ++iterNode )
    {
      const int nodegid = *iterNode;
      const DRT::Node* node = boundarydis_->gNode(nodegid);
      const DRT::Element*const* elements = node->Elements();
      for (int iele=0; iele < node->NumElement(); ++iele)
      {
        const DRT::Element* boundaryele = elements[iele];
        elementsByLabel[label].insert(boundaryele->Id());
      }
    }
  }
  
  int numOfCollectedIds = 0;
  for(unsigned int i = 0; i < elementsByLabel.size(); i++)
    numOfCollectedIds += elementsByLabel[i].size();
  
  if (boundarydis_->NumMyColElements() != numOfCollectedIds)
    dserror("not all elements collected.");
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| evaluate potential conditions based on a Epetra_FecrsMatrix        |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::EvaluatePotentialCondition(
    ParameterList&                          params,
    RefCountPtr<LINALG::SparseMatrix>       systemmatrix1,
    RefCountPtr<LINALG::SparseMatrix>       systemmatrix2,
    RefCountPtr<Epetra_Vector>              systemvector1,
    Teuchos::RCP<Epetra_Vector>             systemvector2,
    Teuchos::RCP<Epetra_Vector>             systemvector3,
    const string&                           condstring)
{
  if (!discret_.Filled()) dserror("FillComplete() was not called");
  if (!discret_.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time < 0.0) usetime = false;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through potential conditions and evaluate them 
  //----------------------------------------------------------------------
  vector<DRT::Condition*> potentialcond;
  discret_.GetCondition(condstring, potentialcond);
  for(vector<DRT::Condition*>::iterator condIter = potentialcond.begin() ; condIter != potentialcond.end(); ++ condIter)
  {
    map<int,RefCountPtr<DRT::Element> >& geom = (*condIter)->Geometry();
    // if (geom.empty()) dserror("evaluation of condition with empty geometry");
    // no check for empty geometry here since in parallel computations
    // can exist processors which do not own a portion of the elements belonging
    // to the condition geometry
    map<int,RefCountPtr<DRT::Element> >::iterator curr;

    // Evaluate Loadcurve if defined. Put current load factor in parameterlist
    const vector<int>*    curve  = (*condIter)->Get<vector<int> >("curve");
    int                   curvenum = -1;
    if (curve) curvenum = (*curve)[0];
    double                curvefac = 1.0;
    if (curvenum>=0 && usetime)
      curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);
      
    params.set("LoadCurveFactor",curvefac);

    params.set<RefCountPtr<DRT::Condition> >("condition", Teuchos::rcp(*condIter,false));

    // define element matrices and vectors
    Epetra_SerialDenseMatrix elematrix1;
    Epetra_SerialDenseMatrix elematrix2;
    Epetra_SerialDenseVector elevector1;
    Epetra_SerialDenseVector elevector2;
    Epetra_SerialDenseVector elevector3;

    for (curr=geom.begin(); curr!=geom.end(); ++curr)
    {
      // get element location vector and ownerships
      vector<int> lm;
      vector<int> lmowner;
      curr->second->LocationVector(discret_,lm,lmowner);
      const int rowsize = lm.size();

      // Reshape element matrices and vectors and init to zero in element->evaluate
      // call the element specific evaluate method
      int err = curr->second->Evaluate( params,discret_,lm, elematrix1,elematrix2,
                                        elevector1,elevector2,elevector3);
      if (err) dserror("error while evaluating elements");
   
      // specify lm row and lm col
      vector<int> lmrow;
      vector<int> lmcol;
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
      if (assemblemat1) systemmatrix1->FEAssemble(eid,elematrix1,lmrow,lmcol);
      if (assemblemat2) systemmatrix2->FEAssemble(eid,elematrix2,lmrow,lmcol);
      
      if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lmrow,lmowner);
      if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lmrow,lmowner);
      if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lmrow,lmowner);
    }
  } // for(vector<DRT::Condition*>::iterator condIter = potentialcond->begin() ; condIter != potentialcond->end(); ++ condIter)
  return;
} // end of EvaluatePotentialCondition



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| compute internal force vector and stiffness matrix                 |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::computeFandK(
   const DRT::Element*              actEle,
   const DRT::UTILS::GaussRule2D&   gaussrule,
   std::map<int,std::set<int> >&    potElements,
   vector<int>&                     lm,
   Epetra_SerialDenseMatrix&        K_surf,
   Epetra_SerialDenseVector&        F_int,
   const int                        label,
   const double                     depth,
   const double                     rootDist,
   const double                     curvefac)
{
  
  // determine global row indices (lmrow) and global colum indices (lm)
  vector<int> lmrow = lm;
  CollectLmcol(potElements, lm);
  // resize matrix and vector and zero out
  const int ndofrow    = lmrow.size();  
  const int ndofcol    = lm.size(); 
  F_int.Size(ndofrow);
  K_surf.Shape(ndofrow, ndofcol);

  // number of atoms (~0.2 nm) per surface area in reference configuration
  // here equal for all bodies in n/µm^2
  const double beta = 25000000;
  const DRT::UTILS::IntegrationPoints2D intpoints = getIntegrationPoints2D(gaussrule);
  
  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute func, deriv, x_gp and factor 
    const int numnode = actEle->NumNode();
    LINALG::SerialDenseVector   funct(numnode);
    LINALG::SerialDenseMatrix   deriv(2,numnode);
    BlitzVec3                   x_gp = 0.0;
   
    const double fac = ComputeFactor(actEle, funct, deriv, intpoints, gp, x_gp, curvefac);
    //----------------------------------------------------------------------
    // run over all influencing elements called element_pot
    //----------------------------------------------------------------------
    
    for(std::map<int, std::set<int> >::const_iterator labelIter = potElements.begin(); labelIter != potElements.end(); labelIter++)
      for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
      {
         const DRT::Element* element_pot = boundarydis_->gElement(*eleIter);
         
         // obtain current potential dofs
         vector<int> lmpot;
         vector<int> lmowner;
         element_pot->LocationVector(*boundarydis_,lmpot,lmowner);
         
         // obtain Gaussrule and integration points
         DRT::UTILS::GaussRule2D rule_pot = DRT::UTILS::intrule2D_undefined;
         switch (element_pot->Shape())
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
         const DRT::UTILS::IntegrationPoints2D intpoints_pot = getIntegrationPoints2D(rule_pot);
         
         //----------------------------------------------------------------------
         // run over all gauss points of a influencing element
         //----------------------------------------------------------------------
         for (int gp_pot = 0; gp_pot < intpoints_pot.nquad; gp_pot++)
         {
           // compute func, deriv, x_gp and factor 
           const int numnode_pot = element_pot->NumNode();
           LINALG::SerialDenseVector  funct_pot(numnode_pot);
           LINALG::SerialDenseMatrix  deriv_pot(2,numnode_pot);
           BlitzVec3                  x_pot_gp = 0.0;
           
           const double fac_pot = ComputeFactor(element_pot, funct_pot, deriv_pot, intpoints_pot, 
                                                gp_pot, x_pot_gp, curvefac);
           
            // evaluate Lennard Jones potential and its derivatives    
            BlitzVec3     potderiv1;
            BlitzMat3x3   potderiv2;
            EvaluateLennardJonesPotential(depth, rootDist, x_gp, x_pot_gp, potderiv1, potderiv2);
            //cout << "potderiv1 = " << potderiv1 << endl;
            //cout << "potderiv2 = " << potderiv2 << endl;
          
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
| (private)                                               umay  07/08|
|                                                                    |
| update displacements in redundant boundary discretization          |
| from solid discretization                                          |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::EvaluateLennardJonesPotential(
    const double        depth,
    const double        rootDist,
    const BlitzVec3&    x,
    const BlitzVec3&    y,
    BlitzVec3&          potderiv1,
    BlitzMat3x3&        potderiv2)
{

  // evaluate distance related stuff
  double          distance      = 0.0;
  BlitzVec3       distance_vec  = 0.0;
  BlitzVec3       distance_unit = 0.0;
  BlitzMat3x3     du_tensor_du;
  computeDistance(x,y, du_tensor_du, distance_vec, distance_unit, distance);

  //----------------------------------------------------------------------
  // evaluate 1.derivative dphi/du_i
  //----------------------------------------------------------------------
  const double dpotdr = ((12*depth)/rootDist)*((-1)*pow((double)(rootDist/distance),13) + pow((double)(rootDist/distance), 7));
  for(int i = 0; i < 3; i++)
    potderiv1(i) = dpotdr*distance_unit(i);

  //----------------------------------------------------------------------
  // evaluate 2.derivative dphi/du_i d_uiI  (this is not a mistake !!!!)
  //----------------------------------------------------------------------
  const double dpotdrdr = ((12*depth)/(rootDist*rootDist))*(13*pow((double)(rootDist/distance),14) - 7*pow((double)(rootDist/distance), 8));
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      potderiv2(i,j) = 0.0;
  
  for(int i = 0; i < 3; i++)
  {
      potderiv2(i,i) += dpotdr/distance;
      for(int j = 0; j < 3; j++)
        potderiv2(i,j) += (dpotdrdr - (dpotdr/distance))*du_tensor_du(i,j);
  }
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| computes distance vector, distance, the distance unit vector       |
| and r_unit tensorproduct r_unit                                    |
*--------------------------------------------------------------------*/
void UTILS::PotentialManager::computeDistance(
    const BlitzVec3&  x,
    const BlitzVec3&  y,
    BlitzMat3x3&      du_tensor_du,
    BlitzVec3&        dist_vec,
    BlitzVec3&        dist_unit,
    double&           distance)
{
  // compute distance vector
  dist_vec(0) = x(0) - y(0);
  dist_vec(1) = x(1) - y(1);
  dist_vec(2) = x(2) - y(2);

  // compute distance
  distance = GEO::Norm2(dist_vec);

  // compute distance unit vector
  dist_unit(0) = dist_vec(0)/distance;
  dist_unit(1) = dist_vec(1)/distance;
  dist_unit(2) = dist_vec(2)/distance;

  // compute r_unit tensorproduct tensor_product
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      du_tensor_du(i,j) = dist_unit(i)*dist_unit(j); 

}



/*----------------------------------------------------------------------*
 |  determine global column indices for the stiffness matrix u.may 08/08|
 |  with nonlocal entries                                               |  
 *----------------------------------------------------------------------*/
void UTILS::PotentialManager::CollectLmcol(
    std::map<int,std::set<int> >&     potentialElementIds,
    vector<int>&                      lmcol)
{
  for(std::map<int, std::set<int> >::const_iterator labelIter = potentialElementIds.begin(); labelIter != potentialElementIds.end(); labelIter++)
    for(std::set<int>::const_iterator eleIter = (labelIter->second).begin(); eleIter != (labelIter->second).end(); eleIter++)
    {
      const DRT::Element* element = boundarydis_->gElement(*eleIter);
      vector<int> lmowner;
      vector<int> lm;
      element->LocationVector(*boundarydis_,lm,lmowner);
      
      for(int i = 0; i < (int) lm.size(); i++)
      {
        bool alreadyInserted = false;
        for(int j = 0; j < (int) lmcol.size(); j++)
        {
          if(lm[i]==lmcol[j])
          {
            alreadyInserted = true;
            break;
          } 
        }
        if(!alreadyInserted)
          lmcol.push_back(lm[i]);
      }
    }
  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| return local index in stiffness matrix of                          |
*--------------------------------------------------------------------*/
int UTILS::PotentialManager::GetLocalIndex(
    vector<int>&    lmcol,
    int             value)
{
  int localindex = -1;
  for(unsigned int i = 0; i < lmcol.size(); i++)
    if(lmcol[i] == value)
    {
      localindex = i;
      break;
    }
  
  if(localindex == -1)
    dserror("no local index found");
        
  return localindex;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in reference configuration         u.may 08/08|
 *----------------------------------------------------------------------*/
void UTILS::PotentialManager::ReferenceConfiguration(
    const DRT::Element*         element,
    LINALG::SerialDenseMatrix&  X, 
    const int                   numdim) const
{
  const int numnode = element->NumNode();
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j) 
      X(i,j) = element->Nodes()[i]->X()[j];

  return;
}



/*----------------------------------------------------------------------*
 |  get nodal coordinates in spatial configuration           u.may 08/08|
 *----------------------------------------------------------------------*/
void UTILS::PotentialManager::SpatialConfiguration(
    const DRT::Element*         element,
    LINALG::SerialDenseMatrix&  x, 
    const int                   numdim) const
{
  const int numnode = element->NumNode();
  for (int i=0; i<numnode; ++i)
    for (int j = 0; j < numdim; ++j) 
      x(i,j) = currentboundarypositions_.find(element->NodeIds()[i])->second(j);
  return;
}



/*----------------------------------------------------------------------*
 |  compute factor funct, deriv, x_gp,                       u.may 08/08|
 *----------------------------------------------------------------------*/
double UTILS::PotentialManager::ComputeFactor(
    const DRT::Element*                     element, 
    LINALG::SerialDenseVector&              funct, 
    LINALG::SerialDenseMatrix&              deriv, 
    const DRT::UTILS::IntegrationPoints2D&  intpoints,
    const int                               gp,
    BlitzVec3&                              x_gp,
    const double                            curve_fac)
{
  
  const int numnode = element->NumNode();
  const double e0 = intpoints.qxg[gp][0];
  const double e1 = intpoints.qxg[gp][1];
  
  // get shape functions and derivatives of the element
  DRT::UTILS::shape_function_2D(funct,e0,e1,element->Shape());
  DRT::UTILS::shape_function_2D_deriv1(deriv,e0,e1,element->Shape());
  
  LINALG::SerialDenseMatrix dXYZdrs(2,3);
  LINALG::SerialDenseMatrix X(numnode,3);
  ReferenceConfiguration(element,X,3);
  LINALG::SerialDenseMatrix x(numnode,3);
  SpatialConfiguration(element,x,3);
  dXYZdrs.Multiply('N','N',1.0,deriv,X,0.0);
  LINALG::SerialDenseMatrix  metrictensor(2,2);
  metrictensor.Multiply('N','T',1.0,dXYZdrs,dXYZdrs,0.0);
  
  // detA maps the reference configuration to the parameter space domain
  const double detA = sqrt(  metrictensor(0,0)*metrictensor(1,1)
                            -metrictensor(0,1)*metrictensor(1,0));
  double factor = intpoints.qwgt[gp] * detA * curve_fac;
  
  x_gp = 0.0;
  // compute gauss point in physical coordinates
  for (int inode = 0; inode < numnode; inode++)
    for(int dim = 0; dim < 3; dim++)
      x_gp(dim) += funct(inode)*x(inode,dim);
  
  return factor;
}


#endif

