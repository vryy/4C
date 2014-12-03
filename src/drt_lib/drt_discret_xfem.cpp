/*!----------------------------------------------------------------------
\file drt_discret_xfem.cpp

\brief a class to manage one discretization with changing dofs

<pre>
Maintainer: Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 - 15249
</pre>
*/

/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>

#include "drt_discret_xfem.H"
#include "../drt_xfem/xfem_fluiddofset.H"
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

  Teuchos::RCP<XFEM::FluidDofSet> xfds = Teuchos::rcp_dynamic_cast<XFEM::FluidDofSet>(initialdofsets_[0]);
  if (xfds != Teuchos::null)
    dserror("DiscretizationXFEM: Initial Dofset shouldn't be a XFEM::FluidDofSet!");

  int numdofspernode = 0;
  fsds->GetReservedMaxNumDofperNode(numdofspernode);
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
void DRT::DiscretizationXFEM::ExportInitialtoActiveVector(Teuchos::RCP<Epetra_Vector>& initialvec,
                                                    Teuchos::RCP<Epetra_Vector>& activevec)
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
  Initialized()
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
