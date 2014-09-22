/*!-------------------------------------------------------------------------*\
 * \file drt_meshfree_discret_fillcomplete.cpp
 *
 * \brief some functions to meshfree disrcetisations
 *
 * <pre>
 * Maintainer: Keijo Nissen (nis)
 *             nissen@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15253
 * </pre>
 *
\*--------------------------------------------------------------------------*/

#include "drt_meshfree_discret.H"
#include "drt_meshfree_node.H"
#include "drt_meshfree_cell.H"
#include "../drt_fem_general/drt_utils_maxent_basisfunctions.H"

/*--------------------------------------------------------------------------*
 | unassigns all node information in all cells           (public) nis Apr12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Unassign()
{
  for (size_t i=0; i<elecolptr_.size(); i++)
    dynamic_cast<DRT::MESHFREE::Cell*>(elecolptr_[i])->DeleteNodes();
  assigned_ = false;
  return;
}

/*--------------------------------------------------------------------------*
 | Reset all maps and conditions                         (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::Reset(bool killdofs, bool killcond)
{
  assigned_ = false;

  // delete all maps and pointers concerning points
  pointrowmap_ = Teuchos::null;
  pointcolmap_ = Teuchos::null;
  pointrowptr_.clear();
  pointcolptr_.clear();

  // call Reset() of base class Discretisation
  DRT::Discretization::Reset(killdofs, killcond);

  return;
}

/*--------------------------------------------------------------------------*
 | Finalize construction                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::FillComplete(
  bool assigndegreesoffreedom,
  bool initelements,
  bool doboundaryconditions)
{
  // this rather puzzling return style is copied from base class Discretization
  int ret = 0;

  if (!Filled())
  {
    // order of fillcompleting is crucial:
    // -> vertices -> line(s) -> surface(s) -> volume
    // (re)assign the nodes' face dimension if boundary conditions are rebuild
    for(int i=0; i<4; ++i)
      for(unsigned j=0; j<domaintopology_[i].size(); ++j)
        domaintopology_[i][j].FillComplete(solutionapprox_->GetRange(), comm_->MyPID(), doboundaryconditions);

    // fill parent discretization
    ret = DRT::Discretization::FillComplete(false,false,false);

    // (re)build map of points pointrowmap_, pointcolmap_, pointrowptr, and pointcolptr
    BuildPointMaps();

    // (re)construct element -> node pointers
    BuildElementToPointPointers();

    // (re)construct node -> element pointers
    BuildPointToElementPointers();

    // assign nodes to geometry
    AssignNodesToPointsAndCells();

    // Assign degrees of freedom to elements and nodes
    if (assigndegreesoffreedom) AssignDegreesOfFreedom(0);

    // call element routines to initialize
    if (initelements) InitializeElements();

    // (Re)build the geometry of the boundary conditions
    if (doboundaryconditions) BoundaryConditionsGeometry();
  }

  return ret;
}

/*--------------------------------------------------------------------------*
 | Build pointrowmap_ and pointcolmap_                  (private) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildPointMaps()
{
  const int myrank = Comm().MyPID();
  int numrowpoints     = 0;
  std::map<int,Teuchos::RCP<DRT::Node> >::iterator curr;
  for (curr=point_.begin(); curr != point_.end(); ++curr)
    if (curr->second->Owner() == myrank)
      ++numrowpoints;
  std::vector<int> pointrowids(numrowpoints);
  pointrowptr_.resize(numrowpoints);

  int numcolpoints = (int)point_.size();
  std::vector<int> pointcolids(numcolpoints);
  pointcolptr_.resize(numcolpoints);

  int countrow=0;
  int countcol=0;
  for (curr=point_.begin(); curr != point_.end(); ++curr){
    if (curr->second->Owner() == myrank)
    {
      pointrowids[countrow] = curr->second->Id();
      pointrowptr_[countrow] = curr->second.get();
      pointcolids[countcol] = pointrowids[countrow];
      pointcolptr_[countcol] = pointrowptr_[countrow];
      ++countrow;
    }
    else
    {
      pointcolids[countcol] = curr->second->Id();
      pointcolptr_[countcol] = curr->second.get();
    }
    curr->second->SetLID(countcol);
    ++countcol;
  }

  if (countrow != numrowpoints) dserror("Mismatch in no. of rowpoints");
  if (countcol != numcolpoints) dserror("Mismatch in no. of colpoints");
  pointcolmap_ = Teuchos::rcp(new Epetra_Map(-1,numcolpoints,&pointcolids[0],0,Comm()));
  pointrowmap_ = Teuchos::rcp(new Epetra_Map(-1,numrowpoints,&pointrowids[0],0,Comm()));
  return;
}

/*--------------------------------------------------------------------------*
 | Build ptrs element -> point                           (private) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildElementToPointPointers()
{
  std::map<int,Teuchos::RCP<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = Teuchos::rcp_dynamic_cast<DRT::MESHFREE::Cell>(elecurr->second,true)->BuildPointPointers(point_);
    if (!success)
      dserror("Building element <-> point topology failed");
  }
  return;
}

/*--------------------------------------------------------------------------*
 | Build ptrs point -> element                           (private) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildPointToElementPointers()
{
  std::map<int,Teuchos::RCP<DRT::Node> >::iterator pointcurr;
  for (pointcurr=point_.begin(); pointcurr != point_.end(); ++pointcurr)
    pointcurr->second->ClearMyElementTopology();

  std::map<int,Teuchos::RCP<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    const int  npoint = elecurr->second->NumPoint();
    const int* points = elecurr->second->PointIds();
    for (int j=0; j<npoint; ++j)
    {
      DRT::Node* point = gPoint(points[j]);
      if (!point) dserror("Point %d is not on this proc %d",j,Comm().MyPID());
      else point->AddElementPtr(elecurr->second.get());
    }
  }

  return;
}
