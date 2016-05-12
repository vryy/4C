/*!----------------------------------------------------------------------
\file drt_discret_xfem.cpp

\brief Implementation a class to manage one discretization with changing dofs

<pre>
\brief Implementation
\level 1
\maintainer Christoph Ager
            http://www.lnm.mw.tum.de
            089 - 289 - 15249
</pre>
*/

/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "drt_discret_xfem.H"
#include "../drt_xfem/xfem_dofset.H"
#include "../linalg/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                             ager 11/14|
 |  comm             (in)  a communicator object                        |
 *----------------------------------------------------------------------*/
DRT::DiscretizationXFEM::DiscretizationXFEM(const std::string name, Teuchos::RCP<Epetra_Comm> comm) :
DiscretizationFaces(name, comm),
initialized_(false),
initialfulldofrowmap_(Teuchos::null),
initialpermdofrowmap_(Teuchos::null)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Finalize construction (public)                             ager 11/14|
 *----------------------------------------------------------------------*/
int DRT::DiscretizationXFEM::InitialFillComplete(bool assigndegreesoffreedom,
                                      bool initelements,
                                      bool doboundaryconditions)
{
  //Call from BaseClass
  int val = DRT::Discretization::FillComplete(assigndegreesoffreedom, initelements, doboundaryconditions);

  if (!assigndegreesoffreedom)
    dserror("DiscretizationXFEM: Call InitialFillComplete() with assigndegreesoffreedom = true!");

  //Store initial dofs of the discretisation
  StoreInitialDofs();
  return val;
}

/*----------------------------------------------------------------------*
 |  checks if Discretization is initialized (protected)  ager 11/14|
 *----------------------------------------------------------------------*/
 bool DRT::DiscretizationXFEM::Initialized() const
{
  if (!initialized_) dserror("DiscretizationXFEM is not initialized! - Call InitialFillComplete() once!");
  return initialized_;
}

/*----------------------------------------------------------------------*
 |  Store Initial Dofs (private)                               ager 11/14|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationXFEM::StoreInitialDofs()
{
  // store copy of initial dofset
  initialdofsets_.clear();
  for (unsigned int dofset = 0; dofset < dofsets_.size(); ++dofset)
  {
    initialdofsets_.push_back(dofsets_[dofset]->Clone());
  }

  // store map required for export to active dofs
  if (initialdofsets_.size() > 1)
    dserror("DiscretizationXFEM: At the moment just one initial dofset is supported by DiscretisationXFEM!");

  Teuchos::RCP<DRT::FixedSizeDofSet> fsds = Teuchos::rcp_dynamic_cast<DRT::FixedSizeDofSet>(initialdofsets_[0]);
  if (fsds == Teuchos::null)
    dserror("DiscretizationXFEM: Cast to DRT::FixedSizeDofSet failed!");

  Teuchos::RCP<XFEM::XFEMDofSet> xfds = Teuchos::rcp_dynamic_cast<XFEM::XFEMDofSet>(initialdofsets_[0]);
  if (xfds != Teuchos::null)
    dserror("DiscretizationXFEM: Initial Dofset shouldn't be a XFEM::XFEMDofSet!");

  int numdofspernode = 0;
  fsds->GetReservedMaxNumDofperNode(numdofspernode);

  if(NumMyColNodes() == 0 ) dserror("no column node on this proc available!");
  int numdofspernodedofset = fsds->NumDof(lColNode(0));
  int numdofsetspernode = 0;

  if (numdofspernode%numdofspernodedofset)
    dserror("DiscretizationXFEM: Dividing numdofspernode / numdofspernodedofset failed!");
  else
    numdofsetspernode = numdofspernode/numdofspernodedofset;

  initialfulldofrowmap_ = ExtendMap(fsds->DofRowMap(),numdofspernodedofset,numdofsetspernode,true);
  initialpermdofrowmap_ = ExtendMap(fsds->DofRowMap(),numdofspernodedofset,numdofsetspernode,false);

  initialized_ = true;

  return;
}

/*------------------------------------------------------------------------------*
 * Export Vector with initialdofrowmap (all nodes have one dofset) - to Vector  |
* with all active dofs (public)                                       ager 11/14|
 *  *---------------------------------------------------------------------------*/
void DRT::DiscretizationXFEM::ExportInitialtoActiveVector(
    Teuchos::RCP<const Epetra_Vector>& initialvec,
    Teuchos::RCP< Epetra_Vector>& activevec)
{
  // Is the discretization initialized?
  Initialized();

  Teuchos::RCP<Epetra_Vector> fullvec = Teuchos::rcp(new Epetra_Vector(*initialpermdofrowmap_,true));

  { //Export manually as target.Map().UniqueGIDs() gives = true, although this shouldn't be the case
    //(UniqueGIDs() just checks if gid occurs on more procs!)
    if (initialvec->Comm().NumProc() == 1 && activevec->Comm().NumProc() == 1) //for one proc , Export works fine!
    {
      LINALG::Export(*initialvec,*fullvec);
    }
    else
    {
      Epetra_Import importer(fullvec->Map(), initialvec->Map());
      int err = fullvec->Import(*initialvec, importer, Insert);
      if (err)
        dserror("Export using exporter returned err=%d", err);
    }
  }
  fullvec->ReplaceMap(*initialfulldofrowmap_); ///replace |1 2 3 4|1 2 3 4| -> |1 2 3 4|5 6 7 8|
  LINALG::Export(*fullvec,*activevec);
}

/*----------------------------------------------------------------------*
 |  get dof row map (public)                                 ager 11/14 |
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DiscretizationXFEM::InitialDofRowMap(unsigned nds) const
{
  Initialized();
  dsassert(nds<initialdofsets_.size(),"undefined initial dof set");

  return initialdofsets_[nds]->DofRowMap();
}


/*----------------------------------------------------------------------*
 |  get dof column map (public)                              ager 11/14 |
 *----------------------------------------------------------------------*/
const Epetra_Map* DRT::DiscretizationXFEM::InitialDofColMap(unsigned nds) const
{
  Initialized();
  dsassert(nds<initialdofsets_.size(),"undefined initial dof set");

  return initialdofsets_[nds]->DofColMap();
}

/*---------------------------------------------------------------------------*
 * Takes DofRowMap with just one xfem-Dofset and duplicates                  |
 * the dof gids for export to active dofs                          ager 11/14|
 *---------------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> DRT::DiscretizationXFEM::ExtendMap(const Epetra_Map* srcmap,
                                                            int numdofspernodedofset,
                                                            int numdofsets,
                                                            bool uniquenumbering)
{
  int numsrcelements = srcmap->NumMyElements();
  const int* srcgids = srcmap->MyGlobalElements();
  std::vector<int> dstgids;
  for (int i = 0; i < numsrcelements; i += numdofspernodedofset)
  {
    if (numsrcelements < i + numdofspernodedofset) dserror("ExtendMap(): Check your srcmap!");
    for (int dofset = 0; dofset < numdofsets; ++dofset)
    {
      for (int dof = 0; dof < numdofspernodedofset; ++dof)
      {
        dstgids.push_back(srcgids[i+dof] + uniquenumbering*dofset*numdofspernodedofset);
      }
    }
  }

  return  Teuchos::rcp(new Epetra_Map(-1,dstgids.size(), &dstgids[0],0,srcmap->Comm()));
}

/*----------------------------------------------------------------------*
 |  set a reference to a data vector (public)                mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::DiscretizationXFEM::SetInitialState(unsigned nds,const std::string& name,Teuchos::RCP<const Epetra_Vector> state)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::DiscretizationXFEM::SetInitialState");

  if (!HaveDofs()) dserror("FillComplete() was not called");
  const Epetra_Map* colmap = InitialDofColMap(nds);
  const Epetra_BlockMap& vecmap = state->Map();

  if (state_.size()<=nds)
    state_.resize(nds+1);

  // if it's already in column map just set a reference
  // This is a rough test, but it might be ok at this place. It is an
  // error anyway to hand in a vector that is not related to our dof
  // maps.
  if (vecmap.PointSameAs(*colmap))
  {
    state_[nds][name] = state;
  }
  else // if it's not in column map export and allocate
  {
#ifdef DEBUG
    if (not InitialDofRowMap(nds)->SameAs(state->Map()))
    {
      dserror("row map of discretization and state vector %s are different. This is a fatal bug!",name.c_str());
    }
#endif
    Teuchos::RCP<Epetra_Vector> tmp = LINALG::CreateVector(*colmap,false);
    LINALG::Export(*state,*tmp);
    state_[nds][name] = tmp;
  }
  return;
}
