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
#include "../drt_xfem/xfsi_searchtree.H"
#include "../drt_adapter/adapter_utils.H"
#include <cstdlib>


/*-------------------------------------------------------------------*
 |  ctor (public)                                          umay 06/08|
 *-------------------------------------------------------------------*/
DRT::PotentialManager::PotentialManager(Teuchos::RCP<DRT::Discretization> discretRCP,
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
  ADAPTER::UTILS::SetupNDimExtractor(*discretRCP_ ,"Potential", potboundary_);
  
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

  
  std::cout << "Potential manager constructor done" << endl;
  
  
  // we start with interfacial area (A_old_) and concentration
  // (con_quot_) = 0. this is wrong but does not make a difference
  // since we apply the equilibrium concentration gradually, thus we
  // do not need these history variables needed for the dynamic model
  // in the beginning.
  A_old_temp_   = rcp(new Epetra_Vector(*surfcolmap,true));
  A_old_        = rcp(new Epetra_Vector(*surfcolmap,true));
  xTree_        = rcp(new XFEM::XSearchTree());
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
  discret_.EvaluateCondition(p,stiff,null,fint,null,null,"Potential");

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

void DRT::PotentialManager::StiffnessAndInternalForces( const int curvenum,
                                                        const double& A,
                                                        Epetra_SerialDenseVector& fint,
                                                        Epetra_SerialDenseMatrix& K_surf,
                                                        const int ID,
                                                        const double time,
                                                        const double dt,
                                                        const double label,
                                                        const double depth,
                                                        const double rootDist,
                                                        const double cutOff)
{

  double traction;
  int LID = A_old_->Map().LID(ID);
  double t_end = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).end();
  set<int> potentialElementIds;
  
  //const TinyVector<double,3> eleCenter = getLocalCenterPosition(       
  //    const DRT::Element::DiscretizationType  distype)      
  
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
    Teuchos::RCP<Epetra_Vector> idisp_solid
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
      vector<double> currpos(3);
      DRT::UTILS::ExtractMyValues(*idisp_total_,mydisp,lm);
      currpos[0] = node->X()[0] + mydisp[0];
      currpos[1] = node->X()[1] + mydisp[1];
      currpos[2] = node->X()[2] + mydisp[2];
      currentboundarypositions_[node->Id()] = currpos;
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
void DRT::PotentialManager::searchElementsInCutOffRadius(
    const vector<double>&   point,
    set< int >&             potentialElementIds, 
    const double            cutOff)
{

const double tol =  1e-7;

  for (int lid = 0; lid < boundarydis_->NumMyColNodes(); ++lid)
  {
    // compute distance between point and node
    const vector<double> currpos = currentboundarypositions_[lid];
    double distance = 0;
    for(int k = 0; k < 3; k++)
    {
      const double comp_diff = currpos[k] - point[k];
      const double diff_square = comp_diff*comp_diff;
      distance += diff_square; 
    }
    distance = sqrt(distance);
    
    // if node lies within cutoff radius -> collect element id in set
    if(fabs(distance) - fabs(cutOff) < tol)
    {
      const DRT::Node* node = boundarydis_->lColNode(lid);
      const DRT::Element *const* elements = node->Elements();
      for(int k = 0; k <  node->NumElement(); k++)
        potentialElementIds.insert(elements[k]->Id());
    }
  }
}



#endif

