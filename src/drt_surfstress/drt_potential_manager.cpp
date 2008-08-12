 /*!-------------------------------------------------------------------
\file drt_potential_manager.cpp

\brief Class controlling surface stresses due to potential forces
between interfaces
and containing all necessary history data

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
#include "../drt_lib/drt_utils.H"
#include <cstdlib>


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
DRT::PotentialManager::PotentialManager(
    Teuchos::RCP<DRT::Discretization> discretRCP,
    DRT::Discretization& discret):
discretRCP_(discretRCP),
discret_(discret)
{
  surfrowmap_ = DRT::UTILS::GeometryElementMap(discret_, "Potential", false);
  RCP<Epetra_Map> surfcolmap = DRT::UTILS::GeometryElementMap(discret_, "Potential", true);
 //vector< DRT::Condition * > potentialConds;

  // get elements with potential condition
  /*discret.GetCondition("LJ_Potential", potentialConds);
  if(potentialConds.size()==0)
    dserror("number of potential conditions = 0");
  
  for(unsigned int i=0; i<potentialConds.size(); i++)
  {
    surfaceEleMap_.push_back( potentialConds[i]->Geometry() );
  }
  */
  
  // create discretization from potential boundary elements and distribute them
  // on every processor
  vector<string> conditions_to_copy;
  conditions_to_copy.push_back("Potential");
  
  //.Discretization()
  boundarydis_ = DRT::UTILS::CreateDiscretizationFromCondition(discretRCP_, "Potential", "PotBoundary", "BELE3", conditions_to_copy);
  dsassert(boundarydis_->NumGlobalNodes() > 0, "empty discretization detected. Potential conditions applied?");

  // split dof vector of soliddiscretization into a dof vector of potential boundary condition and 
  // remaining dofs
  DRT::UTILS::SetupNDimExtractor(*discretRCP_ ,"Potential", potboundary_);
  
  // create potential surface dof row map using the solid parallel distribution
  // for all potential surface elements belonging to the solid discretization on a single 
  // processor
  const Teuchos::RCP< const Epetra_Map > potsurface_dofrowmap_onOneProc = potboundary_.CondMap();
    
  cout << "numglobal  " << boundarydis_->NumGlobalElements() << endl;
  cout << "nummyrow   " << boundarydis_->NumMyRowElements() << endl;
  cout << "nummycol   " << boundarydis_->NumMyColElements() << endl;

  for(int i  = 0; i < boundarydis_->NumGlobalElements(); i++)
	boundarydis_->gElement(i)->Print(cout);

   cout << endl;
  // create node and element distribution with elements and nodes ghosted on all processors
  const Epetra_Map noderowmap = *boundarydis_->NodeRowMap();
  std::cout << "noderowmap->UniqueGIDs(): " << noderowmap.UniqueGIDs() << endl;
  std::cout << noderowmap << endl;

  for(int i  = 0; i < boundarydis_->NumGlobalElements(); i++)
	boundarydis_->gElement(i)->Print(cout);

  Teuchos::RCP<Epetra_Map> newnodecolmap = LINALG::AllreduceEMap(noderowmap);
  std::cout << *newnodecolmap << endl;

  // create redundant boundary discretization on all procs
  DRT::UTILS::RedistributeWithNewNodalDistribution(*boundarydis_, noderowmap, *newnodecolmap);

  // create potential surface dof row map based on the boundayr dis
  // that is distributed redundantly on all processors
  const Epetra_Map* potsurface_dofrowmap_total = boundarydis_->DofRowMap();

  
  // create disp vector for boundary discretization corresponding to the 
  // potential condition elements stored on one proc
  idisp_onproc_    = LINALG::CreateVector(*potsurface_dofrowmap_onOneProc,true);
  // create disp vector for boundary discretization redundant on all procs
  idisp_total_    = LINALG::CreateVector(*potsurface_dofrowmap_total,true);

  // we start with interfacial area (A_old_) and concentration
  // (con_quot_) = 0. this is wrong but does not make a difference
  // since we apply the equilibrium concentration gradually, thus we
  // do not need these history variables needed for the dynamic model
  // in the beginning.
  A_old_temp_   = rcp(new Epetra_Vector(*surfcolmap,true));
  A_old_        = rcp(new Epetra_Vector(*surfcolmap,true));
  
  // determine depth + aabb
  const BlitzMat3x2 rootBox = GEO::getXAABBofDis(*boundarydis_);
  std::map<int,std::set<int> > elementsByLabel;
  CollectElementsByPotentialCouplingLabel(elementsByLabel);
  octTree_      = rcp(new GEO::SearchTree(8));
  octTree_->initializeTree(rootBox, elementsByLabel_, GEO::TreeType(GEO::OCTTREE));
  
  std::cout << "Potential manager constructor done" << endl;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| Call discretization to evaluate additional contributions due to    |
| potential forces                                                   |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::EvaluatePotential(  ParameterList& p,
                                                RefCountPtr<Epetra_Vector> disp,
                                                RefCountPtr<Epetra_Vector> fint,
                                                RefCountPtr<LINALG::SparseMatrix> stiff)
{   
  // has to be called before Evaluate condition !!
  UpdateDisplacementsOfBoundaryDiscretization(disp);
  // action for elements
  p.set("action","calc_potential_stiff");

  discret_.ClearState();
  discret_.SetState("displacement",disp);
  // write conditionon your own copy from discret condition
  //discret_.EvaluateCondition(p,stiff,null,fint,null,null,"Potential");
  
  // set condid -1 so there is no effect
  //EvaluatePotentialCondition(p,stiff,null,fint,null,null,"Potential", -1);

  return;
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| update surface area                                                |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::Update()
{
  A_old_->Update(1.0, *A_old_temp_, 0.0);
}



/*-------------------------------------------------------------------*
| (public)                                                 umay 06/08|
|                                                                    |
| write restart                                                      |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::GetHistory(RCP<Epetra_Vector> A_old_temp_row)
{
  // Note that the temporal vectors need to be written since in the
  // final ones we still have the data of the former step. The column
  // map based vector used for calculations is exported to a row map
  // based one needed for writing.

  LINALG::Export(*A_old_temp_, *A_old_temp_row);
}



/*-------------------------------------------------------------------*
| (public)                                                umay  06/08|
|                                                                    |
| Calculate additional internal forces and corresponding stiffness   |
| on element level for Lennard-Jones potential interaction forces    |
*--------------------------------------------------------------------*/

void DRT::PotentialManager::StiffnessAndInternalForces( const DRT::Element* element,
                                                        const int curvenum,
                                                        const double& A,
                                                        Epetra_SerialDenseVector& fint,
                                                        Epetra_SerialDenseMatrix& K_surf,
                                                        const int ID,
                                                        const double time,
                                                        const double dt,
                                                        const double label,
                                                        const double depth,
                                                        const double rootDist,
                                                        const double radius)
{
  double traction;
  int LID = A_old_->Map().LID(ID);
  double t_end = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).end();
  
  
  // find nodal ids influencing a given element
  std::set<int> potentialElementIds;
  // use corner nodes
  // TODO check for other elements
  for(int i = 0; i < 4; i++)
  {
    const DRT::Node* node = element->Nodes()[i];
    const BlitzVec3 x_node = currentboundarypositions_.find(node->Id())->second;
    octTreeSearchElementsInCutOffRadius(x_node, potentialElementIds, radius);
  }
  
  //const TinyVector<double,3> eleCenter = getLocalCenterPosition(element->Shape());     
  
  
  //searchElementsInCutOffRadius(eleId, potentialElementIds, cutOff);
  printf("not yet fully implemented");
  /*------------------------------------------------- initialization */
  (*A_old_temp_)[LID] = A;

  // find influencing element surfaces
  // later with the help of a tree 
  // so far elemnts are collected if the center lies in a sphere around
  // the center of the current element with a radius = "cutOff"
  // get displacement
  // 
  /*
     potentialElementIDs = collectPotentialElements(surfrowmap_, currentdisp, cutOff);
     
  */
  // compute stresses acting on the current element due to potential forces
  
  /*-----------calculation of current surface stress and its partial
   *-----------------derivative with respect to the interfacial area */

  if (time <= t_end)         /* gradual application of surface stress */
  {
    traction = 1;
  }
  else
  {
    traction = 2;
  }

  double curvefac = 1.;

  /*------------gradual application of surface stresses via time curve*/
  if (time <= t_end)
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

  double ndof = 5; //Adiff.Length();

  for (int i=0;i<ndof;++i)
    for (int j=0;j<ndof;++j)
      K_surf(i,j) = curvefac;
  
  /*------calculation of current internal force due to surface energy*/
  for (int i=0;i<ndof;++i)
    fint[i] = curvefac;


  return;
}



/*-------------------------------------------------------------------*
| (private)                                               umay  07/08|
|                                                                    |
| update displacements in redundant boundary discretization          |
| from solid discretization                                          |
*--------------------------------------------------------------------*/
void DRT::PotentialManager::UpdateDisplacementsOfBoundaryDiscretization(
    Teuchos::RCP<Epetra_Vector>     idisp_solid
    )
{
  // extract displacements associated with potential elements from solid discretization
  idisp_onproc_ = potboundary_.ExtractCondVector(idisp_solid);
  // map to fluid parallel distribution
  LINALG::Export(*idisp_onproc_,*idisp_total_);

  currentboundarypositions_.clear();
  {
    // row and column nodes are identical here, because boundary dis is redundant on all
    // procs
    for (int lid = 0; lid < boundarydis_->NumMyRowNodes(); ++lid)
    {
      const DRT::Node* node = boundarydis_->lRowNode(lid);
      vector<int> lm;
      lm.reserve(3);
      // extract global dof ids
      boundarydis_->Dof(node, lm);
      vector<double> mydisp(3);
      BlitzVec3 currpos(3);
      DRT::UTILS::ExtractMyValues(*idisp_total_,mydisp,lm);
      currpos[0] = node->X()[0] + mydisp[0];
      currpos[1] = node->X()[1] + mydisp[1];
      currpos[2] = node->X()[2] + mydisp[2];
      currentboundarypositions_[node->Id()] = currpos;
    }
  }
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
void DRT::PotentialManager::searchElementsInCutOffRadius(
    const BlitzVec3&        point,
    std::set<int>&          potentialElementIds,
    const double            radius)
{

const double tol =  1e-7;

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
        potentialElementIds.insert(elements[k]->Id());
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
void DRT::PotentialManager::octTreeSearchElementsInCutOffRadius(
    const BlitzVec3&        point,
    std::set<int>&          potentialElementIds,
    const double            radius)
{

const double tol =  1e-7;

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
        potentialElementIds.insert(elements[k]->Id());
    }
  } 
}



/*-------------------------------------------------------------------*
| (private)                                               umay  08/08|
|                                                                    |
| collects all emement ids that belong to certain label              |
*--------------------------------------------------------------------*/
void DRT::PotentialManager::CollectElementsByPotentialCouplingLabel(
    std::map<int,std::set<int> >&        elementsByLabel)
{
  // Reset
  elementsByLabel.clear();
  
  // get condition
  vector< DRT::Condition * >  potConditions;
  boundarydis_->GetCondition ("Potential", potConditions);
  
  // collect elements by xfem coupling label
  for(vector<DRT::Condition*>::const_iterator conditer = potConditions.begin(); conditer!=potConditions.end(); ++conditer)
  {
    DRT::Condition* potCondition = *conditer;
    const int label = potCondition->Getint("label");
    const map<int, RCP<DRT::Element > > geometryMap = potCondition->Geometry();
    std::set< DRT::Element* > boundaryElements;
    // find all elements of this condition
    map<int, RCP<DRT::Element > >::const_iterator iterGeo;
    for(iterGeo = geometryMap.begin(); iterGeo != geometryMap.end(); ++iterGeo )
    {
      DRT::Element*  cutterElement = iterGeo->second.get();
      elementsByLabel[label].insert(cutterElement->Id());
    }
  }
}


/*
void DRT::PotentialManager::EvaluatePotentialCondition(
    ParameterList&                          params,
    RefCountPtr<LINALG::SparseOperator>     systemmatrix1,
    RefCountPtr<LINALG::SparseOperator>     systemmatrix2,
    RefCountPtr<Epetra_Vector>              systemvector1,
    Teuchos::RCP<Epetra_Vector>             systemvector2,
    Teuchos::RCP<Epetra_Vector>             systemvector3,
    const string&                           condstring,
    const int                               condid)
{
  if (!discret_.Filled()) dserror("FillComplete() was not called");
  if (!discret_.HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");

  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time < 0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;

  const bool assemblemat1 = systemmatrix1!=Teuchos::null;
  const bool assemblemat2 = systemmatrix2!=Teuchos::null;
  const bool assemblevec1 = systemvector1!=Teuchos::null;
  const bool assemblevec2 = systemvector2!=Teuchos::null;
  const bool assemblevec3 = systemvector3!=Teuchos::null;

  //----------------------------------------------------------------------
  // loop through conditions and evaluate them if they match the criterion
  //----------------------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first == condstring)
    {
      DRT::Condition& cond = *(fool->second);
      if (condid == -1 || condid ==cond.Getint("ConditionID"))
      {
        map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
        // if (geom.empty()) dserror("evaluation of condition with empty geometry");
        // no check for empty geometry here since in parallel computations
        // can exist processors which do not own a portion of the elements belonging
        // to the condition geometry
        map<int,RefCountPtr<DRT::Element> >::iterator curr;

        // Evaluate Loadcurve if defined. Put current load factor in parameterlist
        const vector<int>*    curve  = cond.Get<vector<int> >("curve");
        int curvenum = -1;
        if (curve) curvenum = (*curve)[0];
        double curvefac = 1.0;
        if (curvenum>=0 && usetime)
          curvefac = UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

        // Get ConditionID of current condition if defined and write value in parameterlist
        const vector<int>*    CondIDVec  = cond.Get<vector<int> >("ConditionID");
        if (CondIDVec)
        {
          params.set("ConditionID",(*CondIDVec)[0]);
          char factorname[30];
          sprintf(factorname,"LoadCurveFactor %d",(*CondIDVec)[0]);
          params.set(factorname,curvefac);
        }
        else
        {
          params.set("LoadCurveFactor",curvefac);
        }
        params.set<RefCountPtr<DRT::Condition> >("condition", fool->second);

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
          curr->second->LocationVector(*this,lm,lmowner);

        // get dimension of element matrices and vectors
        // Reshape element matrices and vectors and init to zero
        const int eledim = (int)lm.size();
        if (assemblemat1) elematrix1.Shape(eledim,eledim);
        if (assemblemat2) elematrix2.Shape(eledim,eledim);
        if (assemblevec1) elevector1.Size(eledim);
        if (assemblevec2) elevector2.Size(eledim);
        if (assemblevec3) elevector3.Size(eledim);

          // call the element specific evaluate method
          int err = curr->second->Evaluate(params,*this,lm,elematrix1,elematrix2,
                                           elevector1,elevector2,elevector3);
          if (err) dserror("error while evaluating elements");

          // assembly
          int eid = curr->second->Id();
          if (assemblemat1) systemmatrix1->Assemble(eid,elematrix1,lm,lmowner);
          if (assemblemat2) systemmatrix2->Assemble(eid,elematrix2,lm,lmowner);
          if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
          if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
          if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);
        }
      }
    }
  } //for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  return;
} // end of DRT::Discretization::EvaluateCondition
*/


#endif

