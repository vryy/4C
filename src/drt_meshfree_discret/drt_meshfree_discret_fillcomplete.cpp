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
void DRT::MESHFREE::MeshfreeDiscretization::Reset(bool killdofs)
{
  assigned_ = false;

  // delete all maps and pointers concerning knots
  knotrowmap_ = Teuchos::null;
  knotcolmap_ = Teuchos::null;
  knotrowptr_.clear();
  knotcolptr_.clear();
  knotghoptr_.clear();

  // call Reset() of base class Discretisation
  DRT::Discretization::Reset(killdofs);

  return;
}


/*--------------------------------------------------------------------------*
 | Finalize construction                                 (public) nis Jan12 |
 *--------------------------------------------------------------------------*/
int DRT::MESHFREE::MeshfreeDiscretization::FillComplete(bool assigndegreesoffreedom,
                                                        bool initelements,
                                                        bool doboundaryconditions)
{
  // this rather puzzling return style is copied from base class discret
  int ret = 0;

  if (!filled_)
  {
    // this rather puzzling return style is copied from base class discret
    ret = DRT::Discretization::FillComplete(false,false,false);

    // (re)build map of knots knotrowmap_, knotcolmap_, knotrowptr, knotcolptr, and knotghoptr
    BuildKnotMaps();

    // (re)construct element -> node pointers
    BuildElementToKnotPointers();

    // (re)construct node -> element pointers
    BuildKnotToElementPointers();

    // second assign nodes to geometry if not assigned
    AssignNodesToKnotsAndCells();

    // Assign degrees of freedom to elements and nodes
    if (assigndegreesoffreedom) AssignDegreesOfFreedom(0);

    // call element routines to initialize
    if (initelements) InitializeElements();

    // (Re)build the geometry of the boundary conditions
    if (doboundaryconditions) BoundaryConditionsGeometry();
  }
//  else
//  {
//    if (!assigned_)
//    {
//      // second assign nodes to geometry if not assigned
//      AssignNodesToKnotsAndCells();
//
//      // Assign degrees of freedom to elements and nodes
//      AssignDegreesOfFreedom(0);
//    }
//  }

  return ret;
}

/*--------------------------------------------------------------------------*
 | Build knotrowmap_, knotcolmap_, knotghoptr_          (private) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildKnotMaps()
{
  const int myrank = Comm().MyPID();
  int numrowknots     = 0;
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::iterator curr;
  for (curr=knot_.begin(); curr != knot_.end(); ++curr)
    if (curr->second->Owner() == myrank)
      ++numrowknots;
  std::vector<int> knotrowids(numrowknots);
  knotrowptr_.resize(numrowknots);

  int numcolknots = (int)knot_.size();
  std::vector<int> knotcolids(numcolknots);
  knotcolptr_.resize(numcolknots);

  std::vector<int> knotghoids(numcolknots-numrowknots);
  knotghoptr_.resize(numcolknots-numrowknots);

  int countrow=0;
  int countcol=0;
  int countgho=0;
  for (curr=knot_.begin(); curr != knot_.end(); ++curr){
    if (curr->second->Owner() == myrank)
    {
      knotrowids[countrow] = curr->second->Id();
      knotrowptr_[countrow] = curr->second.get();
      knotcolids[countcol] = knotrowids[countrow];
      knotcolptr_[countcol] = knotrowptr_[countrow];
      ++countrow;
    }
    else
    {
      knotghoids[countgho] = curr->second->Id();
      knotghoptr_[countgho] = curr->second.get();
      knotcolids[countcol] = knotghoids[countgho];
      knotcolptr_[countcol] = knotghoptr_[countgho];
      ++countgho;
    }
    curr->second->SetLID(countcol);
    ++countcol;
  }

  if (countrow != numrowknots) dserror("Mismatch in no. of rowknots");
  if (countcol != numcolknots) dserror("Mismatch in no. of colknots");
  if (countgho != (numcolknots-numrowknots)) dserror("Mismatch in no. of ghoknots");
  knotcolmap_ = Teuchos::rcp(new Epetra_Map(-1,numcolknots,&knotcolids[0],0,Comm()));
  knotrowmap_ = Teuchos::rcp(new Epetra_Map(-1,numrowknots,&knotrowids[0],0,Comm()));
  return;
}

/*--------------------------------------------------------------------------*
 | Build ptrs element -> knot                           (private) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildElementToKnotPointers()
{
  std::map<int,Teuchos::RCP<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    bool success = Teuchos::rcp_dynamic_cast<DRT::MESHFREE::Cell>(elecurr->second,true)->BuildKnotPointers(knot_);
    if (!success)
      dserror("Building element <-> knot topology failed");
  }
  return;
}

/*--------------------------------------------------------------------------*
 | Build ptrs knot -> element                           (private) nis Jan12 |
 *--------------------------------------------------------------------------*/
void DRT::MESHFREE::MeshfreeDiscretization::BuildKnotToElementPointers()
{
  std::map<int,Teuchos::RCP<DRT::MESHFREE::MeshfreeNode> >::iterator knotcurr;
  for (knotcurr=knot_.begin(); knotcurr != knot_.end(); ++knotcurr)
    knotcurr->second->ClearMyElementTopology();

  std::map<int,Teuchos::RCP<DRT::Element> >::iterator elecurr;
  for (elecurr=element_.begin(); elecurr != element_.end(); ++elecurr)
  {
    Teuchos::RCP<DRT::MESHFREE::Cell> cell = Teuchos::rcp_dynamic_cast<DRT::MESHFREE::Cell>(elecurr->second,true);
    const int  nknot = cell->NumKnot();
    const int* knots = cell->KnotIds();
    for (int j=0; j<nknot; ++j)
    {
      DRT::Node* knot = gKnot(knots[j]);
      if (!knot) dserror("Knot %d is not on this proc %d",j,Comm().MyPID());
      else knot->AddElementPtr(elecurr->second.get());
    }
  }

  return;
}
