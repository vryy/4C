/*----------------------------------------------------------------------*/
/*!
\file xfluid_timeInt.cpp

\brief provides the xfluid timeIntegration,
       maps vectors from old interface position to vectors at new interface position,
       determines the reconstruction method for missing and unreasonable ghost and standard values

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>

*/
/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>


#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_gmsh.H"

#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_element.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_cutwizard.H"


#include "../drt_inpar/inpar_xfem.H"

#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "xfem_dofset.H"
#include "xfem_condition_manager.H"
#include "xfluid_timeInt.H"

#include <iostream>


//#define DEBUG_TIMINT


// -------------------------------------------------------------------
// constructor
// -------------------------------------------------------------------
XFEM::XFluidTimeInt::XFluidTimeInt(
    const bool                                                                  is_newton_increment_transfer, /// monolithic newton increment transfer or time step transfer?
    const Teuchos::RCP<DRT::Discretization> &                                   dis,                    /// discretization
    const Teuchos::RCP<XFEM::ConditionManager> &                                condition_manager,      /// condition manager
    const Teuchos::RCP<GEO::CutWizard> &                                        wizard_old,             /// cut wizard at t^n
    const Teuchos::RCP<GEO::CutWizard> &                                        wizard_new,             /// cut wizard at t^(n+1)
    const Teuchos::RCP<XFEM::XFEMDofSet> &                                      dofset_old,             /// XFEM dofset at t^n
    const Teuchos::RCP<XFEM::XFEMDofSet> &                                      dofset_new,             /// XFEM dofset at t^(n+1)
    const INPAR::XFEM::XFluidTimeIntScheme                                      xfluid_timintapproach,  /// xfluid_timintapproch
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> > &                   node_to_reconstr_method,/// reconstruction map for nodes and its dofsets
    std::map<INPAR::XFEM::XFluidTimeInt, std::map<int,std::set<int> > > &       reconstr_method_to_node,/// inverse reconstruction map for nodes and its dofsets
    const int                                                                   step,                   /// global time step
    const bool                                                                  xfluid_timint_check_interfacetips ///check interfacetips?
) :
is_newton_increment_transfer_(is_newton_increment_transfer),
dis_(dis),
condition_manager_(condition_manager),
wizard_old_(wizard_old),
wizard_new_(wizard_new),
dofset_old_(dofset_old),
dofset_new_(dofset_new),
timeint_scheme_ (xfluid_timintapproach),
node_to_reconstr_method_(node_to_reconstr_method),
reconstr_method_to_node_(reconstr_method_to_node),
step_(step),
xfluid_timint_check_interfacetips_(xfluid_timint_check_interfacetips)
{

  myrank_  = dis->Comm().MyPID();
  numproc_ = dis->Comm().NumProc();

  permutation_map_ = Teuchos::rcp(new std::map<int,int>);

  return;
} // end constructor


// -------------------------------------------------------------------
// set and print reconstruction status for nodes
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::SetAndPrintStatus(const bool screenout)
{

  reconstr_counts_.clear();
  for(std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >::const_iterator node_it=node_to_reconstr_method_.begin();
      node_it!= node_to_reconstr_method_.end();
      node_it++)
  {
    const std::vector<INPAR::XFEM::XFluidTimeInt> & nodesets = node_it->second;

    for(std::vector<INPAR::XFEM::XFluidTimeInt>::const_iterator sets= nodesets.begin();
        sets != nodesets.end();
        sets++)
    {
      std::map<INPAR::XFEM::XFluidTimeInt,int>::iterator it = reconstr_counts_.find(*sets);

      if(it!=reconstr_counts_.end())
        (it->second)++; // increase counter
      else
        reconstr_counts_.insert(std::pair<INPAR::XFEM::XFluidTimeInt,int>(*sets,1)); // initialize counter with 1
    }
  }

  int nummethods=INPAR::XFEM::Xf_TimeInt_undefined+1; // has to be larger than the maximum of the enum INPAR::XFEM::XFluidTimeInt

  Epetra_SerialDenseVector cpu_methods(nummethods);
  Epetra_SerialDenseVector glob_methods(nummethods);

  for (int i=0; i<nummethods; ++i)
  {
    cpu_methods(i) = 0.0;
    glob_methods(i) = 0.0;
  }

  for(std::map<INPAR::XFEM::XFluidTimeInt,int>::iterator reconstrMethod = reconstr_counts_.begin();
      reconstrMethod != reconstr_counts_.end();
      reconstrMethod++)
  {
    int index = (int)(reconstrMethod->first);
    cpu_methods(index) = reconstrMethod->second;
  }

  // reduce and sum over all procs
  dis_->Comm().SumAll(cpu_methods.Values(), glob_methods.Values(), nummethods);

  if(screenout)
  {
    IO::cout << "\n+-------------------------------------------------------+"
             << "\nProc " << myrank_ << ": XFEM::XFluidTimeInt::Status:" << IO::endl;

    for(int method_idx=0; method_idx< nummethods; method_idx++)
    {
      printf("\n%s:\t #dofsets(P%i/allprocs):\t(%i/%i)",
             MapMethodEnumToString(INPAR::XFEM::XFluidTimeInt(method_idx)).c_str(),
             myrank_,
             (int)cpu_methods(method_idx),
             (int)glob_methods(method_idx) );
    }
    IO::cout << "\n+-------------------------------------------------------+\n" << IO::endl;
  }

  return;
}


/*----------------------------------------------------------------------*
| returns matching std::string for each reconstruction method   schott 07/12 |
*----------------------------------------------------------------------*/
std::string XFEM::XFluidTimeInt::MapMethodEnumToString
(
   const enum INPAR::XFEM::XFluidTimeInt term
)
{
  // length of return std::string is 14 due to usage in formated screen output
  switch (term)
  {
  case INPAR::XFEM::Xf_TimeInt_STD_by_SL:
    return "Semi-Lagrange: STD(n+1)                    ";
    break;
  case INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST:
    return "Copy Dofset:   GHOST(n)   -> STD(n+1)      ";
    break;
  case INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD:
    return "Copy Dofset:   STD(n)     -> STD(n+1)      ";
    break;
  case INPAR::XFEM::Xf_TimeInt_GHOST_by_GP:
    return "Ghost-Penalty: GHOST(n+1)                  ";
    break;
  case INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST:
    return "Copy Dofset:   GHOST(n)   -> GHOST(n+1)    ";
    break;
  case INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_STD:
    return "Copy Dofset:   STD(n)     -> GHOST(n+1)    ";
    break;
  case INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS:
    return "Copy Dofset:   STD_EMB(n) -> STD/GHOST(n+1)";
    break;
  case INPAR::XFEM::Xf_TimeInt_undefined:
    return "undefined:                                 ";
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }

  return "";
} // ScaTraTimIntImpl::MapTimIntEnumToString


// -------------------------------------------------------------------
// transfer standard and ghost dofs to new map as far as possible and
// mark all dofs for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::TransferDofsToNewMap(
  const std::vector<Teuchos::RCP<const Epetra_Vector> >&   oldRowStateVectors,  /// row map based vectors w.r.t old interface position
  const std::vector<Teuchos::RCP<Epetra_Vector> >&         newRowStateVectors,  /// row map based vectors w.r.t new interface position
  const Teuchos::RCP<std::set<int> >                       dbcgids              /// set of dof gids that must not be changed by ghost penalty reconstruction
)
{
  //------------------------------------------------------------
  // loop background fluid row nodes
  const Epetra_Map* noderowmap = dis_->NodeRowMap();
  const int numrownode = noderowmap->NumMyPoints();

  for (int lid=0; lid<numrownode; ++lid)
  {
    // get global id of a node
    int gid = noderowmap->GID(lid);
    TransferNodalDofsToNewMap(
        oldRowStateVectors,
        newRowStateVectors,
        dbcgids,
        gid);
  }

  // export the reconstruction method information to other procs
  // this has to be done, when SL-nodes mark non-row surrounding nodes (and its some of its dofsets) as Ghost-penalty dofsets
  ExportMethods(newRowStateVectors, dbcgids);
}

// -------------------------------------------------------------------
// transfer standard and ghost dofs to new map as far as possible and
// mark all dofs for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::TransferDofsToNewMap(
  const std::vector<Teuchos::RCP<const Epetra_Vector> >&   oldRowStateVectors,  /// row map based vectors w.r.t old interface position
  const std::vector<Teuchos::RCP<Epetra_Vector> >&         newRowStateVectors,  /// row map based vectors w.r.t new interface position
  const Teuchos::RCP<std::set<int> >                       dbcgids,             /// set of dof gids that must not be changed by ghost penalty reconstruction
  const std::vector<int>&                                  node_gids            /// vector of node gids
)
{
  for (size_t i=0; i<node_gids.size(); ++i)
  {
    // get global id of a node
    TransferNodalDofsToNewMap(
        oldRowStateVectors,
        newRowStateVectors,
        dbcgids,
        node_gids[i]);
  }
  // export the reconstruction method information to other procs
  // this has to be done, when SL-nodes mark non-row surrounding nodes (and its some of its dofsets) as Ghost-penalty dofsets
  ExportMethods(newRowStateVectors, dbcgids);
}

// -------------------------------------------------------------------
// transfer standard and ghost dofs to new map as far as possible
// and mark nodal dofs for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::TransferNodalDofsToNewMap(
  const std::vector<Teuchos::RCP<const Epetra_Vector> >&   oldRowStateVectors,  /// row map based vectors w.r.t old interface position
  const std::vector<Teuchos::RCP<Epetra_Vector> >&         newRowStateVectors,  /// row map based vectors w.r.t new interface position
  const Teuchos::RCP<std::set<int> >                       dbcgids,             /// set of dof gids that must not be changed by ghost penalty reconstruction
  int                                                      gid                  /// nodal gid
)
{
#ifdef DEBUG_TIMINT
  const int numdofpernode = 4;
#endif

  if(oldRowStateVectors.size() != newRowStateVectors.size()) dserror("TransferDofsToNewMap: not equal number of old and new vectors");

  // get the node
  DRT::Node* node = dis_->gNode(gid);

  // get cut nodes with respect to cut wizards at t^n and t^(n+1)
  GEO::CUT::Node * n_new = wizard_new_->GetNode(gid);
  GEO::CUT::Node * n_old = wizard_old_->GetNode(gid);


  //---------------------------------------------------------------------------------
  // switch over different cases dependent of surrounding elements are cut or not at t^n and t^(n+1)
  //     case A: surrounding elements not cut at t^n AND t^(n+1) => copy dofs
  //     case B: at least one surrounding element cut at t^n, uncut elements at t^(n+1)
  //     case C: uncut elements at t^n, but at least one surrounding element cut at t^(n+1)
  //     case D: surrounding elements cut at old and new timestep t^n and t^(n+1)
  //---------------------------------------------------------------------------------

  // initialize to true, as in case that no cut information we have standard FE and the dofsets are unique standard dofsets
  bool unique_std_uncut_np = true;
  bool unique_std_uncut_n = true;

  // check for unique std dofset and surrounding uncut elements at t^(n+1)
  if(n_new != NULL)
  {
    const int numDofSets_new = n_new->NumDofSets(); //= dof_cellsets_new.size()

    // just one dofset at new timestep t^(n+1)
    if(numDofSets_new == 1) unique_std_uncut_np = NonIntersectedElements(node, wizard_new_);
    else                    unique_std_uncut_np = false;
  }

  // check for unique std dofset and surrounding uncut elements at t^n
  if(n_old != NULL)
  {
    const int numDofSets_old = n_old->NumDofSets(); //= dof_cellsets_old.size()

    // just one dofset at new timestep t^n
    if(numDofSets_old == 1) unique_std_uncut_n = NonIntersectedElements(node, wizard_old_);
    else                    unique_std_uncut_n = false;
  }

  //===========================================================
  // switch cases A-D
  //===========================================================
  if(   (n_new == NULL       and n_old == NULL     )
      or(n_new == NULL       and unique_std_uncut_n)
      or(unique_std_uncut_np and n_old == NULL     )
      or(unique_std_uncut_np and unique_std_uncut_n)
      )
  {
    //---------------------------------------------------------------------------------
    // case A: surrounding elements not cut at t^n AND t^(n+1) => copy dofs
    //---------------------------------------------------------------------------------

    // holes are included within the bounding box
    // REMARK: we assume that the structure does not move more than one element and
    //         nodes do not change the side w.r.t structure (not to large movement!!!)

    // just one dofset available at t^n and t^(n+1)
    const int nds_new = 0;
    const int nds_old = 0;

    // copy dofs from this node from t^n -> t^(n+1)
    CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);

  }//end case A
  else if(n_new == NULL or unique_std_uncut_np)
  {
    if(n_old == NULL) dserror("you should call case A here");

    //---------------------------------------------------------------------------------
    // case B: at least one surrounding element cut at t^(n), uncut elements at t^(n+1)
    //---------------------------------------------------------------------------------

    // case a): a similar uncut position at old and new timestep (but node in bounding box at old timestep)
    // case b): at least one std-dof at t^n (  -> std-dof can be copied, std-dof should not have changed the side within one timestep!)
    //                                         -> all ghost-dofs have to be computed newly)
    // case c): only ghost dofs at t^n (-> to much structural movement within one timestep)

    //-------------------------------
    // t^(n+1)
    // assume one unique dofset
    const int nds_new = 0;
    GEO::CUT::Point::PointPosition pos_new = GEO::CUT::Point::undecided;

    if(n_new == NULL)
      pos_new = GEO::CUT::Point::outside; // by default for nodes outside the cut-boundingbox
    else
    {
      if (n_new->NodalDofSets().size() == 1)
        pos_new = n_new->NodalDofSets()[0]->Position();
      else
        dserror("XFEM::XFluidTimeInt::TransferDofsToNewMap() - Case B: n_new->NodalDofSets().size() != 1");
    }

    if(pos_new != GEO::CUT::Point::outside and pos_new != GEO::CUT::Point::inside)
      dserror("position of unique std-dofset is not inside and not outside, can this happen?");

    std::vector<int> dofs_new;
    dofset_new_->Dof(dofs_new, node, nds_new );

#ifdef DEBUG_TIMINT
    const int numdofs_new = (int)(dofs_new.size());
    if(numdofs_new != numdofpernode)
    {
      dserror("XFLUID-TIMINIT CASE B: node %d,\t %d dofs for node without nodehandle at timestep t^(n+1) ?!", gid, numdofs_new);
    }
#endif

    //-------------------------------
    // t^(n)

    // how many sets of dofs?
    const int numDofSets_old = n_old->NumDofSets(); //= dof_cellsets_old.size()

    //-------------------------------------
    // just one dofset at old timestep t^n
    //-------------------------------------
    if(numDofSets_old == 1 and unique_std_uncut_n)
    {
      dserror("here you should call case A!");
    }
    else if(numDofSets_old == 1 and !unique_std_uncut_n)
    {
      // check if there is a standard dofset
      int nds_old = n_old->GetStandardNodalDofSet(pos_new);

      if(nds_old != -1) // case b)
      {
        // copy values
        CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
      }
      else // case c)
      {
        //case: NEWLY CREATED NODEHANDLE: just one ghost dofset at old time!
        //      structural movement more than one element near node, (here SEMILAGRANGE possible!)

        if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
          dserror("structural movement more than one element near node %d", gid);
        else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                 timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
        }
        else if ( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS,dbcgids);
        }
        else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

      }
    } // some elements cut
    else if(numDofSets_old == 0 )
    {
      // case: "XFLUID-TIMINIT CASE B: node %d,\t no dofset at t^n available,
      //        structural movement more than one element for node ? (here SEMILAGRANGE possible!", gid);

      if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
        dserror("no dofset at t^n available, structural movement more than one element for node %d", gid);
      else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
               timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
      {
        MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
      }
      else if ( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
      {
        MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS,dbcgids);
      }
      else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

    }
    else
    {
      //-------------------------------------
      // more than one dofset at old timestep
      //-------------------------------------

      // check if there is a standard dofset
      int nds_old = n_old->GetStandardNodalDofSet(pos_new);

      if( nds_old == -1 )
      {
        // all dofsets at t^n are not std-dofsets, structural movement more than one element

        if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
          dserror("all dofset at t^n are ghost-dofsets, structural movement more than one element for node %d", gid);
        else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                 timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
        }
        else if ( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS,dbcgids);
        }
        else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

      }
      else
      {
        // copy values
        CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
      }
    } //end more than one dofset
  }//end case B
  else if(n_old == NULL or unique_std_uncut_n)
  {
    if(n_new == NULL) dserror("you should call case A here");

    // how many sets of dofs?
    const int numDofSets_new = n_new->NumDofSets(); //= dof_cellsets_new.size()

    if(numDofSets_new == 0 ) // no values to set
    {
      return; // do nothing for this node
    }

    //---------------------------------------------------------------------------------
    // case C: C: uncut elements at t^n, but at least one surrounding element cut at t^(n+1)
    //---------------------------------------------------------------------------------

    // case a): a similar uncut position at old and new timestep (but node now in bounding box)
    // case b): at least one std-dof at t^(n+1)(  -> std-dof can be copied, std-dof should not have changed the side within one timestep!)
    //                                            -> all ghost-dofs have to be computed newly)
    // case c): only ghost dofs at t^(n+1) (-> to much structural movement within one timestep)

    //-------------------------------
    // t^n
    // assume one unique dofset
    const int nds_old = 0;

    std::vector<int> dofs_old;
    dofset_old_->Dof(dofs_old, node, nds_old );

#ifdef DEBUG_TIMINT
    const int numdofs_old = (int)(dofs_old.size());
    if(numdofs_old != numdofpernode)
    {
      dserror("XFLUID-TIMINIT CASE C: node %d,\t %d dofs for node without nodehandle at timestep t^n ?!", gid, numdofs_old);
    }
#endif

    //-------------------------------
    // t^(n+1)
    const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > & dof_cellsets_new = n_new->NodalDofSets();



    //------------------------------------------
    // just one dofset at new timestep t^(n+1)
    //------------------------------------------
    if(numDofSets_new == 1 and unique_std_uncut_np)
    {
      dserror("here you should call case A!");
    }
    else if(numDofSets_new == 1 and !unique_std_uncut_np)
    {
      const int nds_new = 0;

      // get the unique cellset
      const GEO::CUT::NodalDofSet* nodaldofset = &*(dof_cellsets_new[nds_new]);

      if(nodaldofset->Is_Standard_DofSet()) // case b)
      {
        // copy values or SL
        if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP or
            timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
            timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
        {
          CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
        }
        else if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
        }
      }
      else // case c)
      {
        // REMARK: Ghost penalty reconstruction is more safe than the copy method GHOST_by_COPY_from_STD
        // the ghost dof can be also newly created w.r.t a second other interface,
        // e.g. at t^n a fluid node between two contacting thin structures can become a ghost node w.r.t to a neighbored independent fluid domain
        // (F|S|F -node- F|S|F -> F|S-node-S|F -- F|S|F)
        MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
      }
    } // just one dofset at new timestep
    else
    {
      //------------------------------------------
      // more than one dofset at new timestep t^(n+1)
      //------------------------------------------

      // loop new dofsets
      for(std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> >::const_iterator sets=dof_cellsets_new.begin();
          sets!=dof_cellsets_new.end();
          sets ++)
      {
        Teuchos::RCP<GEO::CUT::NodalDofSet> nodaldofset_new = *sets;
        int nds_new = sets-dof_cellsets_new.begin();

        if(nodaldofset_new->Is_Standard_DofSet()) // first dofset (has been checked to be a std-dofset)
        {
          // copy values or SL
          if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP or
              timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
              timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
          {
            CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
          }
          else if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
          }
        }
        else
        {
          // for newly created ghost dofsets we have to reconstruct ghost dofs
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
        }
      }
    } // more than one dofset
  } //end case C
  else
  {
    if(n_new == NULL or n_old == NULL) dserror("this case should be done before");

    // how many sets of dofs?
    const int numDofSets_new = n_new->NumDofSets(); //= dof_cellsets_new.size()

    if(numDofSets_new == 0 ) // no values to set
    {
      return; // do nothing for this node
    }

    //---------------------------------------------------------------------------------
    // D: surrounding elements cut at old and new timestep t^n and t^(n+1)
    //---------------------------------------------------------------------------------

    // case a): std-set at t^(n+1) -> try to identify vc-connection
    //                                    std found: - Check special case of tips with changed side!
    //                                                  -> Yes: Copy value
    //                                                  -> No:  Semilagrange
    //                                    ghost found and only one ghost dof:
    //                                                   -> copy (or Semilagrange)
    //                                    ghost found but not unique
    //                                                   -> Semilagrange
    // case b): ghost dofset at t^(n+1) -> try to identify vc-connection
    //                                    std or ghost found -> Copy
    //                                    not found -> mark for Ghost Penalty!

    //-------------------------------
    // t^n
    const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > & dof_cellsets_old = n_old->NodalDofSets();


    //-------------------------------
    // t^(n+1)
    const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > & dof_cellsets_new = n_new->NodalDofSets();

    // loop new dofsets
    for(std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> >::const_iterator sets=dof_cellsets_new.begin();
              sets!=dof_cellsets_new.end();
              sets ++)
    {

      const GEO::CUT::NodalDofSet* cell_set = &**sets;

      int nds_new = sets-dof_cellsets_new.begin(); // nodal dofset counter

      //is current set at t^(n+1) std or ghost or dofset
      bool is_std_set_np = cell_set->Is_Standard_DofSet();

      //-----------------------------------------------------------------
      // identify cellsets at t^n with current dofset at t^(n+1)
      int nds_old = IdentifyOldSets(n_old, n_new, dof_cellsets_old, cell_set); // get nds number of corresponding old dofset
      //-----------------------------------------------------------------

      if (timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
      {
        // check, if we already have an entry for that node!
        if (node_to_reconstr_method_.find(gid) != node_to_reconstr_method_.end())
        {
          // disable safety check in case of double cut if we project from other discretization
          nds_old = 0;

          if (nds_old < 0)
            dserror("Projection failed and there is no old dofset to be copied for node %d", gid);
        }
      }

      if(nds_old < 0) // no set or not a unique set found
      {
#ifdef DEBUG_TIMINT
        IO::cout << "XFLUID-TIMINIT CASE D: node " << gid << ",\t no corresponding dofset found at time t^n for dofset " << nds_new << " at time t^(n+1)" << IO::endl;
#endif
        if(is_std_set_np) // std at t^(n+1)
        {
          if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
            dserror("no corresponding dofset at t^n, choose a SL-based-approach for node %d", gid);
          else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                   timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
          }
          else if (timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS,dbcgids);
          }
          else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");
        }
        else // ghost at t^(n+1)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
        }
      }
      else // unique set found
      {
        bool is_std_set_n = dof_cellsets_old[nds_old]->Is_Standard_DofSet();

#ifdef DEBUG_TIMINT
        IO::cout << "XFLUID-TIMINIT CASE D: node " << gid << ",\t corresponding std(yes/no: " << is_std_set_n << ") dofset "
                                                   << nds_old << " found at time t^n for dofset " << nds_new << " at time t^(n+1)" << IO::endl;
#endif

        if(is_std_set_np and is_std_set_n) // std at t^n and t^(n+1)
        {
          // copy values or use SL
          if (timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP or
              timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
              timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
          {
            CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
          }
          else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
          }
          else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");
        }
        else if(is_std_set_np and !is_std_set_n) // ghost at t^n and std at t^(n+1)
        {
          // check number of ghost-dofsets
          // get the number of ghost dofsets -> this check is done in IdentifyOldSets
          bool ghost_set_unique = true;

          if(ghost_set_unique) // ghost dofset is unique
          {
            // copy values

            if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP )
            {
              CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST, newRowStateVectors, oldRowStateVectors,dbcgids);
            }
            else if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP)
            {
#if(1) // is a copy from ghost->std reasonable or not?
            CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST, newRowStateVectors, oldRowStateVectors,dbcgids);
#else
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
#endif
            }
            else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP )
            {
              MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
            }
            else if ( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP)
            {
              MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS,dbcgids);
            }
            else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

          }
          else
          {
            dserror("case should not be called");
            // not unique ghost dofset
            //MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_SemiLagrange,dbcgids);
          }
        }
        else if(!is_std_set_np and is_std_set_n) // ghost at t^(n+1)
        {
          // copy values, should be reasonable since dofsets identified
          CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
//            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
        }
        else if(!is_std_set_np and !is_std_set_n)
        {
          // copy values, should be reasonable since dofsets identified
          //dserror("Copy for node %d", gid);
          CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST, newRowStateVectors, oldRowStateVectors,dbcgids);
//            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
        }
      } // end found set at t^n

    } // loop new dofsets
  }//end case D (node handles at t^n and t^(n+1))

  return;
}

// -------------------------------------------------------------------
// get nodes and dofsets for given reconstruction method
// -------------------------------------------------------------------
std::map<int,std::set<int> > & XFEM::XFluidTimeInt::Get_NodeToDofMap_For_Reconstr(INPAR::XFEM::XFluidTimeInt reconstr)
{
  // returns empty map if reconstruction method is not present
  return reconstr_method_to_node_[reconstr];
}

// -------------------------------------------------------------------
// all surrounding elements non-intersected ?
// -------------------------------------------------------------------
bool XFEM::XFluidTimeInt::NonIntersectedElements(DRT::Node* n, const Teuchos::RCP<GEO::CutWizard> wizard)
{
  const int numele = n->NumElement();

  DRT::Element** elements= n->Elements();

  // loop surrounding elements
  for(int i=0; i< numele; i++)
  {
    DRT::Element * e = elements[i];

    // we have to check elements and its sub-elements in case of quadratic elements
    GEO::CUT::ElementHandle* ehandle = wizard->GetElement(e);

    // elements which do not have an element-handle are non-intersected anyway
    if(ehandle == NULL) continue;

    // check if the element is intersected or not.
    // If the element is not intersected or is just fully or partially touched at a facet
    // a unique volume-cell for the element (each sub-element for quadratic elementhandles holding the same position info)
    // has been produced
    if(ehandle->IsIntersected())
    {
      return false; // at least one element is intersected by a cut-side
    }

  }

  return true; //all surrounding elements uncut
}


// -------------------------------------------------------------------
// find all ghost dofsets around this node and its std-dofset
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::FindSurroundingGhostDofsets(
    std::map<int,std::set<int> >&          ghostDofsets,   /// surrounding ghost dofsets to be filled, map of ghost nodes and correponding ghost nds index w.r.t given std nodal dofset
    const DRT::Node*                       node,           /// node
    const int                              nds_new         /// dofset of node used for finding the surrounding ghost dofs
    )
{
  if(nds_new != 0 ) dserror("do you really want to find ghost dofsets surrounding a non-std node with id %i?", node->Id());

  GEO::CUT::Node * n = wizard_new_->GetNode(node->Id());

  if(n == NULL)
  {
    // it can happen that the node was std at t^n, however cut elements around it (nodehandle available at t^n and SL called)
    // and at new time t^(n+1) it is a std-node without a nodehandle, then there are no surrounding ghost dofs at t^(n+1)
    return;
  }

  const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > & dof_cellsets = n->NodalDofSets();


  // get the corresponding cellset
  const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> & cellset = dof_cellsets[nds_new]->VolumeCellComposite();

  // get for each plain_volumecell_set of the surrounding elements the ghost dofs

  // loop the surrounding (sub-)elements
  for(std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>::const_iterator e_vcset=cellset.begin();
      e_vcset != cellset.end();
      e_vcset++)
  {
    const GEO::CUT::plain_volumecell_set & vcs= *e_vcset;

    // get the element, ask the first vc
    int peid = (*vcs.begin())->ParentElement()->GetParentId();
    DRT::Element * ele = dis_->gElement(peid);
    DRT::Node* * nodes = ele->Nodes();

    const std::vector<int> & nds = (*vcs.begin())->NodalDofSet();

    // which dofset is the corresponding to the current vc-connection and is a ghost dofset?
    for(int n_it=0; n_it < ele->NumNode(); n_it++)
    {
      const DRT::Node* ghost_node = nodes[n_it];
      const int ghost_nid = ghost_node->Id();

      GEO::CUT::Node* ghost_node_cut = wizard_new_->GetNode(ghost_nid);

      if(ghost_node_cut == NULL) continue; // this node is then a standard node or not on this proc

      // check if the neighbored node is a ghost node w.r.t to the cellset of the SL-node
      // if ghost node w.r.t. the cellset then add it to ghostDofsets with the corresponding nds-number
      if(! (dof_cellsets[nds_new]->Contains(ghost_node_cut->point())) )
      {
        std::map<int, std::set<int> >::iterator map_it = ghostDofsets.find(ghost_nid);
        if(map_it==ghostDofsets.end())
        {
          std::set<int> tmp_map;

          tmp_map.insert(nds[n_it]); // insert in map, use map as for hex20 elements, the subelements
          ghostDofsets.insert(std::pair<int, std::set<int> >(ghost_nid,tmp_map));
        }
        else
        {
          map_it->second.insert(nds[n_it]);
        }
      }
    }
  }
}



// -------------------------------------------------------------------
// copy dofs from old vectors to new vector for all row vectors
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::CopyDofs(
    const DRT::Node*                                          node,               /// drt node
    const int                                                 nds_new,            /// nodal dofset at t^(n+1)
    const int                                                 nds_old,            /// nodal dofset at t^n
    const INPAR::XFEM::XFluidTimeInt                          method,             /// reconstruction method
    const std::vector<Teuchos::RCP<Epetra_Vector> >&          newRowStateVectors, /// row map based state vectors at t^(n+1)
    const std::vector<Teuchos::RCP<const Epetra_Vector> >&    oldRowStateVectors, /// row map based state vectors at t^n
    const Teuchos::RCP<std::set<int> >                        dbcgids             /// set of DBC global ids)
    )
{

  if(method != INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST and
     method != INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_STD   and
     method != INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST   and
     method != INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD) dserror("don't call CopyDofs for non-copy reconstruction type");


  std::vector<int> dofs_old;
  std::vector<int> dofs_new;

  dofset_new_->Dof(dofs_new, node, nds_new );
  dofset_old_->Dof(dofs_old, node, nds_old );

#ifdef DEBUG_TIMINT
  if(dofs_new.size() != dofs_old.size())
  {
    dserror("XFLUID-TIMINIT: unequal number of dofs for node %d at old timestep (set %d, %d dofs) and new timestep (set %d, %d dofs)", node->Id(), nds_old, dofs_old.size(), nds_new, dofs_new.size());
  }
#endif

  // try to set the reconstruction method, in case that it has been already set to Ghost-penalty, do not overwrite it
  if(!SetReconstrMethod(node, nds_new, method)) return;

  int vec_count = 0;

  // copy values for all vectors
  for(std::vector<Teuchos::RCP<const Epetra_Vector> >::const_iterator it=oldRowStateVectors.begin(); it!=oldRowStateVectors.end(); it++)
  {

    Teuchos::RCP<Epetra_Vector>       vec_new = newRowStateVectors[vec_count];
    Teuchos::RCP<const Epetra_Vector> vec_old = oldRowStateVectors[vec_count];

    // copy values from old vector to new vector
    for(size_t i=0; i<dofs_new.size(); i++)
    {
      int dof_lid_new = vec_new->Map().LID(dofs_new[i]);
      if(dof_lid_new == -1) dserror("new dof %d not local on this proc!", dofs_new[i]);

      int dof_lid_old = vec_old->Map().LID(dofs_old[i]);
      if(dof_lid_old == -1) dserror("old dof %d not local on this proc!", dofs_old[i]);

      (*vec_new)[dof_lid_new] = (*vec_old)[dof_lid_old];

      // set Dirichlet BC for ghost penalty reconstruction
      int gid = dofs_new[i];

      if (dbcgids != Teuchos::null)
        (*dbcgids).insert(gid);

    }

    vec_count++;
  }

  // create permutation cycles
  if( ( method == INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST
     or method == INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_STD
     or method == INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST
     or method == INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD // no permutation for std-to-std as std is always the first, nevertheless std-to-std reasonable for some cases?
     )
     and nds_new != nds_old)
  {
    //std::cout << "copying from ghost to ghost for node " << node->Id() << " set old: " << nds_old << " set new: " << nds_new << std::endl;

    // copy values from old vector to new vector
    for(size_t i=0; i<dofs_new.size(); i++)
    {
      int dof_gid_new = dofs_new[i];
      int dof_gid_old = dofs_old[i];

      permutation_map_->insert(std::pair<int,int>(dof_gid_old, dof_gid_new));
    }
  }
  return;
}


// -------------------------------------------------------------------
// mark nodal dofs of vector w.r.t new interface position for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::MarkDofs(
    const DRT::Node*                                   node,                 /// drt node
    const int                                          nds_new,              /// nodal dofset at t^(n+1)
    const std::vector<Teuchos::RCP<Epetra_Vector> >&   newRowStateVectors,   /// row map based state vectors at t^(n+1)
    const INPAR::XFEM::XFluidTimeInt                   method,               /// reconstruction method
    const Teuchos::RCP<std::set<int> >                 dbcgids               /// set of dof gids that must not be changed by ghost penalty reconstruction
)
{

  // try to set the reconstruction method, in case that it has been already set to Ghost-penalty, do not overwrite it
  if(!SetReconstrMethod(node, nds_new, method)) return;


  //-------------------------------------
  // for the Semi-Lagrangean algorithm applied to single nodes we reconstruct the surrounding ghost nodes with the Ghost-penalty approach
  // this is done to avoid kinks in the reconstructed field when only std-values are computed using SL but ghost-values would be copied
  // in case of the full SL-algorithm in the boundary zone, all ghost-values will be reconstructed using the GP-approach (this is done in the TransferDofsToNewMap directly)
  if( (timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
       timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP ) and
       method == INPAR::XFEM::Xf_TimeInt_STD_by_SL)
  {
    // map of global nid of a neighbored ghost node and the corresponding ghost-dofset-number which corresponds to the given std dofset at the node
    std::map<int, std::set<int> > ghostDofsets;
    FindSurroundingGhostDofsets(ghostDofsets, node, nds_new);

    for(std::map<int,std::set<int > >::iterator it=ghostDofsets.begin(); it!=ghostDofsets.end(); it++)
    {
      // mark the ghost dofset in case that the node is a row node on this proc
      // otherwise mark it for export at the end of TransferDofsToNewMap
      int nid=it->first;
      std::set<int>& nds = it->second;
      for(std::set<int>::iterator nds_it=nds.begin(); nds_it!=nds.end(); nds_it++)
      {
        const int dofset = *nds_it; // dofset number for ghost node

        // comment this block if we do not want to reconstruct ghost values of ghost nodes around SL-standard nodes
#if(1)
        if(dis_->NodeRowMap()->LID(nid) != -1) // is row node on this proc
        {
          DRT::Node* n = dis_->gNode(nid);
          MarkDofs(n, dofset, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
        }
        else
        {
          // in parallel we have to export the info to other procs
          MarkDofsForExport(nid, dofset, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP);
        }
#endif
      }
    }
  }
  //-------------------------------------

  // get nodal dofs for current dofset w.r.t new interface position
  std::vector<int> dofs_new;
  dofset_new_->Dof(dofs_new, node, nds_new );

  // loop vectors
  for(std::vector<Teuchos::RCP<Epetra_Vector> >::const_iterator it=newRowStateVectors.begin(); it!=newRowStateVectors.end(); it++)
  {

    Teuchos::RCP<Epetra_Vector> vec_new = *it;

    // set a dummy value for dofs in the vector
    for(size_t i=0; i<dofs_new.size(); i++)
    {
      int dof_lid_new = vec_new->Map().LID(dofs_new[i]);
      if(dof_lid_new == -1) dserror("new dof %d not local on this proc!", dofs_new[i]);

      (*vec_new)[dof_lid_new] = 0.0; // zero value as start value for newton for gradient reconstruction

      int gid = dofs_new[i];

      // set Dirichlet BC for ghost penalty reconstruction
      if(method != INPAR::XFEM::Xf_TimeInt_GHOST_by_GP)
      {
        if (dbcgids != Teuchos::null)
          (*dbcgids).insert(gid);
      }
      else
      {
        // find if the dbc has been already set by a previous set non-ghost-penalty reconstruction
        // and ensure that no dbc is set
        std::set<int>::iterator it = (*dbcgids).find(gid);
        if(it!=(*dbcgids).end()) (*dbcgids).erase(gid); // remove the already set dbc
      }
    } // dofs
  } // state vectors

  return;
}

// -------------------------------------------------------------------
// mark one specific nodal dofset with used for export
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::MarkDofsForExport(
    const int                        nid,         /// node id
    const int                        dofset,      /// ghost dofset number
    const INPAR::XFEM::XFluidTimeInt method       /// reconstruction method used for marking the nodal dofset
    )
{
  std::map<int, std::map<int, int> >::iterator it = dofset_marker_export_.find(nid);

  if(it!=dofset_marker_export_.end())
  {
    std::map<int, int > & dofsets_map = it->second;
    std::map<int, int >::iterator nds_it = dofsets_map.find(dofset);
    if(nds_it != dofsets_map.end())
    {
      if(nds_it->second != (int)method) dserror("you want to mark one dofset twice with different methods! nid: %d, dofset: %d, set method: %d, new method: %d", nid, dofset, (int)nds_it->second, (int)method);
      else return; // already set
    }
    else
    {
      dofsets_map.insert(std::pair<int,int>(dofset,(int)method));
    }
  }
  else
  {
    std::map<int, int> dofsets_map;
    dofsets_map.insert(std::pair<int,int>(dofset,(int)method));

    dofset_marker_export_.insert(std::pair<int,std::map<int, int> >(nid, dofsets_map));
  }
}


// -------------------------------------------------------------------
// set the reconstruction method for current nodal dofset, return if set
// -------------------------------------------------------------------
bool XFEM::XFluidTimeInt::SetReconstrMethod(
    const DRT::Node*                     node,                 /// drt node
    const int                            nds_new,              /// nodal dofset w.r.t new interface position
    const INPAR::XFEM::XFluidTimeInt     method                /// which type of reconstruction method
    )
{
  GEO::CUT::Node* n = wizard_new_->GetNode(node->Id());

  int numdofsets = -1;
  if(n!=NULL){ numdofsets = n->NumDofSets();}
  else numdofsets = 1; // one std dofset

  // find the node in map
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >::iterator it;
  it= node_to_reconstr_method_.find(node->Id());

  if(it != node_to_reconstr_method_.end())
  {
    // assume that we call this function for all nodal dofsets in the right order
    if(nds_new >= numdofsets) dserror("no valid nds %d for node %d", nds_new, n->Id());

    //---------------------------------------------------------
    // check if the current set method can be overwritten

    // do not overwrite Ghost-penalty reconstruction method which have been already marked by surrounding SL-dofsets
    if(it->second[nds_new] == INPAR::XFEM::Xf_TimeInt_GHOST_by_GP) return false;

    // do not overwrite already values, expect that method is Xf_TimeInt_GHOST_by_GP
    if(it->second[nds_new] != INPAR::XFEM::Xf_TimeInt_undefined and    // overwriting
       it->second[nds_new] != method and                               // methods not equal
       method != INPAR::XFEM::Xf_TimeInt_GHOST_by_GP and               // overwriting with non-GP approach
       it->second[nds_new] != INPAR::XFEM::Xf_TimeInt_by_PROJ_from_DIS // overwriting projection after previous failure (label correction)
       )
    {
      dserror("inconsistency in reconstruction method, why do want to replace reconstruction method %d with %d for node=%d and nds=%d",  it->second[nds_new], method, node->Id(), nds_new );
    }
    //---------------------------------------------------------

    // remove the node from the old method in the reconstr. method-to-node map
    std::map<INPAR::XFEM::XFluidTimeInt, std::map<int,std::set<int> > >::iterator itm;
    itm = reconstr_method_to_node_.find(it->second[nds_new]);

    if (itm != reconstr_method_to_node_.end())
    {
      std::map<int,std::set<int> >::iterator itn = itm->second.find(node->Id());
      if (itn != itm->second.end())
        itn->second.erase(nds_new);

      if (itn->second.empty())
        itm->second.erase(node->Id());
    }

    (it->second)[nds_new] = method;
    reconstr_method_to_node_[method][node->Id()].insert(nds_new);

    return true;
  }
  else
  {
    std::vector<INPAR::XFEM::XFluidTimeInt> vec(numdofsets,INPAR::XFEM::Xf_TimeInt_undefined); // initialize with undefined status
    vec[nds_new] = method;
    node_to_reconstr_method_.insert(std::pair<int,std::vector<INPAR::XFEM::XFluidTimeInt> >(node->Id(), vec ));
    reconstr_method_to_node_[method][node->Id()].insert(nds_new);
    return true;
  }

  return false;
}


// -------------------------------------------------------------------
// identify cellsets at time t^n with cellsets at time t^(n+1)
// -------------------------------------------------------------------
int XFEM::XFluidTimeInt::IdentifyOldSets(
    const GEO::CUT::Node *                                         n_old,               /// node w.r.t to old wizard
    const GEO::CUT::Node *                                         n_new,               /// node w.r.t to new wizard
    const std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> > &      dof_cellsets_old,    /// all dofcellsets at t^n
    const GEO::CUT::NodalDofSet*                                   cell_set_new        /// dofcellset at t^(n+1) which has to be identified
)
{

  //is current set at t^(n+1) std or ghost or dofset
  bool is_std_set_np = cell_set_new->Is_Standard_DofSet();

  const GEO::CUT::Point::PointPosition pos_new = cell_set_new->Position();

  //--------------------------------------------------------
  // t^(n+1)
  // set of side-ids involved in cutting the current connection of volumecells at t^(n+1)
  // get all side-ids w.r.t to all volumecells contained in current new set around the current node

  GEO::CUT::plain_int_set cutsides_new;
  cell_set_new->CollectCutSides( cutsides_new );

  //--------------------------------------------------------
  // t^n

  std::map<int, std::vector<int> > identified_old_sets; // map of possible nds_sets and related identified sides

  //--------------------------------------------------------
  // PRESELECTION via common cutting sides (for level-sets all sets have the same level-set side and are possible sets)
  //--------------------------------------------------------
  //check each old dofset for identification with new dofset
  for(std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> >::const_iterator old_sets=dof_cellsets_old.begin();
      old_sets!=dof_cellsets_old.end();
      old_sets++)
  {
    const int setnumber = old_sets - dof_cellsets_old.begin();

    GEO::CUT::plain_int_set cutsides_old;
    (*old_sets)->CollectCutSides( cutsides_old );

    //--------------------------------------------------------------
    // check if identification of sets is possible
    // -> check if any side is involved in both cuts (find side-Id involved at time t^(n+1) in set of involved side-Ids at t^n)
    //--------------------------------------------------------------

    // get the common sides of both sets!
    std::vector<int> common_sides;
    common_sides.reserve(std::max(cutsides_new.size(), cutsides_old.size()));

    // sorted vectors are already sorted
    if(cutsides_old.size()> 0 and cutsides_new.size()>0)
    std::set_intersection (cutsides_new.begin(), cutsides_new.end(), cutsides_old.begin(), cutsides_old.end(), std::back_inserter(common_sides)); // back_inserter uses pushback


    if(common_sides.size() > 0)
      identified_old_sets[setnumber] = common_sides; // [] creates a new vector of common side-ids if not created yet
  }

  // "rotation" of interface around a point outside the volumecell composite, such that no common cut side is available for the two sets of cut-sides
  // Nevertheless the standard dofset (if there exists one) could be still a good choice
  // check if there is a standard dofset which is reasonable nevertheless
  if(identified_old_sets.size() == 0 // no identification via common cutsides possible
     and is_std_set_np // it is a standard dofset at the new time step, for which the alternative would be only to use SemiLagrangean
     and timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP // just when we need to copy data, otherwise simulation will stop
  )
  {
    // look again in the old sets if there was a standard dofset
    for(std::vector<Teuchos::RCP<GEO::CUT::NodalDofSet> >::const_iterator old_sets=dof_cellsets_old.begin();
        old_sets!=dof_cellsets_old.end();
        old_sets++)
    {
      if((*old_sets)->Is_Standard_DofSet()) // this set is a unique standard set at tn, might it be a good choice nevertheless?
      {
        const int setnumber = old_sets - dof_cellsets_old.begin();

        bool use_old_std_set = false; // shall the unique standard set be used for further checks?

#if(1) // use the standard dofset, however this might be unsafe
        use_old_std_set = true; // use the standard set just in case that there is at least a common node which is shared by the two sets of cutsides

        // NOTE: an alternative would be to check if there is at least one common node or a common edge between the two sets of cutsides, this might be more safe!
#endif

        if(use_old_std_set)
        {
          IO::cout << "WARNING (xfluid_timeInt.cpp): for std-dofset (t^(n+1)) at node " << n_new->Id()
              << " no std dofset at t^n could be identified via common sides."
              << " However, there is a std-dofset (t^n) which could be used nevertheless! This might be unsafe! Be aware of that!" << IO::endl;

          std::vector<int> dummy; // dummy for common sides ( actually no common sides available! )
          identified_old_sets[setnumber] = dummy; // [] creates a new vector of common side-ids if not created yet
        }
      }
    }
  }

  //---------------------------------------
  // try to find a unique candidate and do safety check if there is only one possible set
  //---------------------------------------

  std::set<int> invalid_old_sets; // find invalid cell-sets and erase them from identified_old_sets afterwards

  if(identified_old_sets.size() > 1)
  {

    bool unique_set_found = false;

    for(std::map<int, std::vector<int> >::iterator it= identified_old_sets.begin();
        it!=identified_old_sets.end();
        it++)
    {
      const int nds_old = it->first;

      if(unique_set_found)
      {
        // make the other sets invalid
        invalid_old_sets.insert(nds_old);
        continue;
      }

      const GEO::CUT::Point::PointPosition pos_old = dof_cellsets_old[nds_old]->Position();

      if(pos_old == GEO::CUT::Point::undecided or
          pos_old == GEO::CUT::Point::oncutsurface) dserror("why is the cellcet position undecided or oncutsurface %i, something wrong", pos_old);

      // dof-cellsets at new and old time have to correspond to each other w.r.t position of their cellsets (same fluid phase)
      if(pos_old == pos_new)
      {
        bool is_std_set_n = dof_cellsets_old[nds_old]->Is_Standard_DofSet();

        if(is_std_set_n) // standard sets are unique, therefore accept the unique set and neglect the others
          unique_set_found = true;
        else
          unique_set_found = false;
      }
      else // changing dofset positions, not allowed here
      {
        // make the other sets invalid
        invalid_old_sets.insert(nds_old);
        continue;
      }
    } // loop old sets
  } // identified_dofsets > 1

  //---------------------------------------
  // remove invalid sets
  //---------------------------------------
  for(std::set<int>::iterator set_it=invalid_old_sets.begin();
      set_it!=invalid_old_sets.end();
      set_it++)
  {
    identified_old_sets.erase(*set_it);
  }

  //---------------------------------------
  // no unique reasonable candidate found
  //---------------------------------------
  if(identified_old_sets.size() > 1)
  {
#ifdef DEBUG_TIMINT
    IO::cout << "Warning: found dofset at t^n not unique, found " << identified_old_sets.size() << " dofsets, set status to not-found!" << IO::endl;
#endif
    return -1;
  }

  //---------------------------------------
  // no reasonable candidate available
  //---------------------------------------
  if(identified_old_sets.size() <= 0)
  {
#ifdef DEBUG_TIMINT
    IO::cout << "Warning: no found dofset at t^n, set status to not-found!" << IO::endl;
#endif
    return -1;
  }

  //---------------------------------------
  // special check for the unique candidate / can we really accept the value?
  //---------------------------------------
  if (identified_old_sets.size() == 1)
  {
    const int nds_old = identified_old_sets.begin()->first;

#ifdef DEBUG_TIMINT
    IO::cout << "Exactly one dofset found: nds= " << nds_old << " perform further special checks:" <<IO::endl;
#endif

    bool is_std_set_n  = dof_cellsets_old[nds_old]->Is_Standard_DofSet();

  //---------------------------------------
  // do special checks
  //---------------------------------------


    bool did_node_change_side = false;
    bool successful_check = false;

    //--------------------------------------
    // special check whether node slides on cut surface
    successful_check = SpecialCheck_SlidingOnSurface(did_node_change_side, n_old, n_new);

    if(!successful_check or (successful_check and did_node_change_side))
      return -1; // do not accept the old value

    if(xfluid_timint_check_interfacetips_)
    {
      //--------------------------------------
      // special check for interface tips if the node has changed the side w.r.t identified sides at t^n and t^(n+1)
      if( (is_std_set_np and is_std_set_n) or (!is_std_set_np and !is_std_set_n))
      {
        std::vector<int> & identified_sides = identified_old_sets[nds_old];
        successful_check = SpecialCheck_InterfaceTips(did_node_change_side, identified_sides, n_old,n_new);

        if(!successful_check or (successful_check and did_node_change_side))
          return -1; // do not accept the old value
      }
    }

    //---------------------------------------
    // if the unique candidate passed all checks we accept the value
    //---------------------------------------
#ifdef DEBUG_TIMINT
    IO::cout << "\t ACCEPT the unique dofset!" <<IO::endl;
#endif

    return nds_old; // all special tests passed -> accept the old value
  }

  return -1;
}


bool XFEM::XFluidTimeInt::SpecialCheck_SlidingOnSurface(
    bool &                                                         changed_side,
    const GEO::CUT::Node *                                         n_old,               /// node w.r.t to old wizard
    const GEO::CUT::Node *                                         n_new                /// node w.r.t to new wizard
)
{
  if(n_old->Id() != n_new->Id()) dserror("XFLUID-TIMINT: not the same node for CheckChangingSide");

  //--------------------------------------------------
  // do some simple checks for changing side based on point positions to reduce checks based on space-time sides

  // check when point moves on the surface of the same side
  GEO::CUT::Point* p_old = n_old->point();
  GEO::CUT::Point* p_new = n_new->point();

  GEO::CUT::Point::PointPosition pos_old = p_old->Position();
  GEO::CUT::Point::PointPosition pos_new = p_new->Position();

  if( (pos_old == GEO::CUT::Point::oncutsurface and pos_new == GEO::CUT::Point::oncutsurface) )
  {
    // first case: point slides on the cut surface (within one side or from 'within the side' to point or edge and vice versa)
    // second case: point moves from one side to another side (-> decide if the sides are neighbors and so on)

    const GEO::CUT::plain_side_set & cut_sides_old  = p_old->CutSides();
    const GEO::CUT::plain_side_set & cut_sides_new  = p_new->CutSides();


#ifdef DEBUG_TIMINT
    IO::cout << "!!WARNING point is on cut surface at time t^n and t^(n+1)" << IO::endl;
#endif


    // check if there is at least one common cutside at both times
    if(cut_sides_old.size() == 0 or cut_sides_new.size() == 0) dserror("there are no cutsides but point is 'oncutsurface' at both times");

    // REMARK: we have to check this based on the real sides not based on the triangulated cutsides
    // get the side-Ids

    std::set<int> on_cut_sides_old;
    std::set<int> on_cut_sides_new;

    // loop cutsides and extract side ids
    for(GEO::CUT::plain_side_set::const_iterator side_it= cut_sides_old.begin();
        side_it!=cut_sides_old.end();
        side_it++)
    {
      int sid=(*side_it)->Id();
      if(sid != -1) on_cut_sides_old.insert(sid); // do not insert sides of the background element
    }

    // loop cutsides and extract side ids
    for(GEO::CUT::plain_side_set::const_iterator side_it= cut_sides_new.begin();
        side_it!=cut_sides_new.end();
        side_it++)
    {
      int sid=(*side_it)->Id();
      if(sid != -1) on_cut_sides_new.insert(sid); // do not insert sides of the background element
    }


    // check if there is at least one common cut side the point lies on
    // this checks the case when point moves from inside the side to an edge or cornerpoint of the same side
    //------------------------------------------------------
    // find common side

    while(true)
    {
      // if at least one vector of nodes is checked
      if((int)on_cut_sides_new.size() == 0 or (int)on_cut_sides_old.size() == 0) break;

      if(*(on_cut_sides_new.begin()) < *(on_cut_sides_old.begin()))      on_cut_sides_new.erase(on_cut_sides_new.begin());
      else if(*(on_cut_sides_new.begin()) > *(on_cut_sides_old.begin())) on_cut_sides_old.erase(on_cut_sides_old.begin());
      else // on_cut_sides_new.begin() == on_cut_sides_old.begin()
      {
#ifdef DEBUG_TIMINT
        IO::cout << "point moves within side " << *(on_cut_sides_new.begin()) << " oncutsurface!" << IO::endl;
#endif

        // no changing side
        changed_side = false;
        return true;
      }
    }

    //------------------------------------------------------
    // if this check was not successful, check if two sides are at least neighbors
    dserror("point moves oncutsurface from t^n to t^(n+1), but remains not within the same side, check if a 'changed side' maybe based on side neighbors would be possible");
    //------------------------------------------------------

    changed_side = false;
    return true;
  }

  changed_side = false;
  return true;
}


// -------------------------------------------------------------------
// check if the node has changed the side w.r.t identified sides at t^n and t^(n+1)
// return if check was successful
// -------------------------------------------------------------------
bool XFEM::XFluidTimeInt::SpecialCheck_InterfaceTips(
    bool&                                                          changed_side,        /// did the node change the side ?
    std::vector<int> &                                             identified_sides,    /// side Id of identified side
    const GEO::CUT::Node *                                         n_old,               /// node w.r.t to old wizard
    const GEO::CUT::Node *                                         n_new                /// node w.r.t to new wizard
    )
{
  bool successful_check = true;

  changed_side = false;

  if(n_old->Id() != n_new->Id()) dserror("XFLUID-TIMINT: not the same node for CheckChangingSide");

  //--------------------------------------------------
  // do some simple checks for changing side based on point positions to reduce checks based on space-time sides

  // check when point moves on the surface of the same side
  GEO::CUT::Point* p_old = n_old->point();
  GEO::CUT::Point* p_new = n_new->point();


  GEO::CUT::Point::PointPosition pos_old = p_old->Position();
  GEO::CUT::Point::PointPosition pos_new = p_new->Position();

  if(pos_old == GEO::CUT::Point::undecided or pos_new == GEO::CUT::Point::undecided) dserror("at least one node position undecided!");

  if(    (pos_old == GEO::CUT::Point::outside and pos_new == GEO::CUT::Point::outside)
      or (pos_old == GEO::CUT::Point::inside  and pos_new == GEO::CUT::Point::inside)
  )
  {
    // continue with check based on space time sides
    // REMARK: the point could have moved through an inside/outside volumecell

    // do nothing and continue with space-time side check
  }
  else // otherwise check not required, continue with other special checks
  {
    changed_side = false;
    return true;
  }


  //------------------------------------

  LINALG::Matrix<3,1> n_coord_old(true);
  n_old->Coordinates(&n_coord_old(0,0));

  LINALG::Matrix<3,1> n_coord_new(true);
  n_new->Coordinates(&n_coord_new(0,0));

  // check if moving node (ALE case)
  LINALG::Matrix<3,1> n_diff(true);
  n_diff.Update(1.0, n_coord_new, -1.0, n_coord_old);

  // TODO: for ALE we have to check whether the path of the point crosses at least one space-time side element
  if(n_diff.Norm2()>1e-14)
  {
    // TODO: currently we expect that node did not change the side around a tip
    // TODO: USE the dserror, at the moment we just throw a warning
    IO::cout << "WARNING: node " << n_old->Id() << " seems to move (background ALE?). Check the SpecialCheck_InterfaceTips-Check for ALE!" << IO::endl;

    changed_side = false;
    return true;

    //dserror("node %i seems to move. Check the SpaceTimeChangingSide-Check for ALE!", n_old->Id());
  }



  // check the change of the node w.r.t to sides using space-time side elements
  // create space-time side elements using the same side at t^n and t^(n+1)
  // check if the node is contained in such space-time surface elements at any time between t^n and t^(n+1)

  int side_count = 0;

  // loop sides
  for(std::vector<int>::iterator sides=identified_sides.begin(); sides!= identified_sides.end(); sides++)
  {
    side_count++;

#ifdef DEBUG_TIMINT
    IO::cout <<  "\t CheckChangingSide for node " << n_old->Id() << " w.r.t side " << *sides << IO::endl;
#endif

    const int coup_sid = *sides; // side id used within the cut

    if(condition_manager_->IsMeshCoupling(coup_sid))
    {
      DRT::Element * side = condition_manager_->GetSide(coup_sid);

      successful_check = SpecialCheck_InterfaceTips_SpaceTime(changed_side, side, coup_sid, n_coord_old);
    }
    else
    {
      if(n_diff.Norm2()>1e-14)
        dserror("background fluid ALE with level-sets interface??? Think about that, does Scatra support this?");

      successful_check = SpecialCheck_InterfaceTips_Levelset(changed_side);
    }

    if(!successful_check) return false; // return non-successful check

    if(changed_side)
    {
#ifdef DEBUG_TIMINT
      IO::cout << "\t\t node " << n_new->Id()<< " changed interface side, node within space-time side element for at least one side" << IO::endl;
#endif
      break; // stop loop over sides
    }
  }

  return successful_check;
}

bool XFEM::XFluidTimeInt::SpecialCheck_InterfaceTips_Levelset(
    bool& changed_side        /// did the node change the side ?
)
{
  // interface tips moving/rotating around a node not supported yet
  // TODO: if necessary, implement a check based on the gradients of the level-set field or the normal vectors of the boundary cells
  // possibility: check the different directions of the normal vectors up to a certain tolerance

  // assume that the interface did not change the side and accept the value
  changed_side = false;
  return true;
}

bool XFEM::XFluidTimeInt::SpecialCheck_InterfaceTips_SpaceTime(
    bool&                         changed_side,              /// did the node change the side ?
    DRT::Element *                side,
    const int                     coup_sid,
    const LINALG::Matrix<3,1>&    n_coord                    /// node coodinates
)
{
  bool node_within_Space_Time_Side = false;

  DRT::Element::DiscretizationType side_distype = side->Shape();

  bool successful_check = false;

  switch(side_distype)
  {
    case DRT::Element::tri3:
    {
      successful_check = WithinSpaceTimeSide<DRT::Element::tri3,DRT::Element::wedge6>(node_within_Space_Time_Side, side, coup_sid, n_coord);
      break;
    }
    case DRT::Element::quad4:
    {
      successful_check = WithinSpaceTimeSide<DRT::Element::quad4,DRT::Element::hex8>(node_within_Space_Time_Side, side, coup_sid, n_coord);
      break;
    }
    case DRT::Element::quad8:
    {
      successful_check = WithinSpaceTimeSide<DRT::Element::quad8,DRT::Element::hex16>(node_within_Space_Time_Side, side, coup_sid, n_coord);
      break;
    }
    case DRT::Element::quad9:
    {
      successful_check = WithinSpaceTimeSide<DRT::Element::quad9,DRT::Element::hex18>(node_within_Space_Time_Side, side, coup_sid, n_coord);
      break;
    }
    default: dserror("side-distype %s not handled", DRT::DistypeToString(side_distype).c_str()); break;
  }

  changed_side = node_within_Space_Time_Side;

  return successful_check;
}

// -------------------------------------------------------------------
// check if the node is within the space time side
// -------------------------------------------------------------------
template<DRT::Element::DiscretizationType side_distype, DRT::Element::DiscretizationType space_time_distype>
bool XFEM::XFluidTimeInt::WithinSpaceTimeSide(
    bool&                         within_space_time_side,    /// within the space time side
    DRT::Element *                side,
    const int                     coup_sid,
    const LINALG::Matrix<3,1>&    n_coord                    /// node coodinates
)
{
  // get the right cutter discretization for the given side
  Teuchos::RCP<DRT::Discretization> cutter_dis = condition_manager_->GetCutterDis( coup_sid );

  std::string state_new = "idispnp";
  std::string state_old = "";

  if(is_newton_increment_transfer_)
    state_old = "idispnpi"; // get displacements from last newton increment
  else
    state_old = "idispn";   // get displacements from last time step

  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> idisp_new = cutter_dis->GetState(state_new);
  Teuchos::RCP<const Epetra_Vector> idisp_old = cutter_dis->GetState(state_old);


  const int numnode_space_time = DRT::UTILS::DisTypeToNumNodePerEle<space_time_distype>::numNodePerElement;

  const int numnode_side = DRT::UTILS::DisTypeToNumNodePerEle<space_time_distype>::numNodePerElement/2;

  // space time side coordinates
  LINALG::Matrix<3,numnode_space_time> xyze_st;

  const int numnode = side->NumNode();
  DRT::Node ** nodes = side->Nodes();


  Epetra_SerialDenseMatrix xyze_old(3, numnode);
  Epetra_SerialDenseMatrix xyze_new(3, numnode);

  for (int i = 0; i < numnode; ++i)
  {
    DRT::Node & node = *nodes[i];

    LINALG::Matrix<3, 1> x_old(node.X());
    LINALG::Matrix<3, 1> x_new(node.X());

    std::vector<int> lm;
    std::vector<double> mydisp_old;
    std::vector<double> mydisp_new;

    cutter_dis->Dof(&node, lm);

    DRT::UTILS::ExtractMyValues(*idisp_old, mydisp_old, lm);
    DRT::UTILS::ExtractMyValues(*idisp_new, mydisp_new, lm);

    // add displacements
    x_old(0) += mydisp_old.at(0);
    x_old(1) += mydisp_old.at(1);
    x_old(2) += mydisp_old.at(2);

    x_new(0) += mydisp_new.at(0);
    x_new(1) += mydisp_new.at(1);
    x_new(2) += mydisp_new.at(2);

    std::copy(x_old.A(), x_old.A() + 3, &xyze_old(0, i));
    std::copy(x_new.A(), x_new.A() + 3, &xyze_new(0, i));
  }

  //------------------------------------------------------------------

  for( int i=0; i<numnode_side; ++i )
  {
    for(int j=0; j< 3; j++)
    {
      xyze_st(j,i) = xyze_old(j,i);
      xyze_st(j,i+numnode_side) = xyze_new(j,i);
    }
  }

  // check the space time side volume
  bool successful_check = CheckSTSideVolume<space_time_distype,numnode_space_time>(xyze_st);

  if(!successful_check) // apply small artificial movement in normal direction(translation)
  {
    double perturbation = 1e-008;

#ifdef DEBUG_TIMINT
    IO::cout << "\t\t\t Not Successful! -> Apply small artificial translation in normal direction at t^(n+1) of " << perturbation << IO::endl;
#endif

    //------------------------------------------------------------------
    // get normal vector on side at new interface position at center of side

    LINALG::Matrix<3,1> normal(true);

    LINALG::Matrix<2,1> xi_side(true);

    if(numnode_side == 3)
    {
      xi_side(0) = 0.3333333333333333;
      xi_side(1) = 0.3333333333333333;
    }
    else if(numnode_side == 4 or numnode_side == 8 or numnode_side == 9)
    {
      xi_side(0) = 0.0;
      xi_side(1) = 0.0;
    }
    else dserror("unknown side type with %d nodes", numnode_side);

    // Initialization
    LINALG::Matrix<2,numnode_side> deriv(true);      // derivatives dr, ds

    LINALG::Matrix<3,2> derxy(true);
    LINALG::Matrix<3,1> dx_dr(true);
    LINALG::Matrix<3,1> dx_ds(true);

    // get current values
    DRT::UTILS::shape_function_2D_deriv1( deriv, xi_side( 0 ), xi_side( 1 ), side_distype );

    LINALG::Matrix<3,numnode_side> xyz_side_new (xyze_new);

    derxy.MultiplyNT(xyz_side_new, deriv);


    // set dx_dr and dx_ds
    for (int i=0; i< 3; i++)
    {
      dx_dr(i) = derxy(i,0);
      dx_ds(i) = derxy(i,1);
    }


    normal(0) = dx_dr(1)*dx_ds(2)-dx_ds(1)*dx_dr(2);
    normal(1) = dx_dr(2)*dx_ds(0)-dx_ds(2)*dx_dr(0);
    normal(2) = dx_dr(0)*dx_ds(1)-dx_ds(0)*dx_dr(1);

    normal.Scale(1.0/normal.Norm2());

    //------------------------------------------------------------------

    for( int i=0; i<numnode_side; ++i )
    {
      for(int j=0; j< 3; j++)
      {
        xyze_st(j,i) = xyze_old(j,i);
        xyze_st(j,i+numnode_side) = xyze_new(j,i)+perturbation*normal(j,0);
      }
    }

    successful_check = CheckSTSideVolume<space_time_distype,numnode_space_time>(xyze_st);

#ifdef DEBUG_TIMINT
    if(successful_check) IO::cout << "\t\t\t Successful check!" << IO::endl;
    else                 IO::cout << "\t\t\t Again not successful check!" << IO::endl;
#endif
  }

  if(!successful_check) return successful_check;

  Teuchos::RCP<GEO::CUT::Position> pos =
        GEO::CUT::PositionFactory::BuildPosition<3,space_time_distype>(xyze_st,n_coord);
  within_space_time_side = pos->Compute();

#ifdef DEBUG_TIMINT
  LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t space time element (r,s,t !!!)
  pos->LocalCoordinates(rst);

  if(within_space_time_side)
  {
    IO::cout << "rst " << rst << IO::endl;
    IO::cout << "\t\t\t changing side found!" << "new: " << xyze_new <<  " old: " << xyze_old << " \n global coords " << n_coord << " \n local coords " << rst << IO::endl;
  }
#endif

  return successful_check;
}


// -------------------------------------------------------------------
// check the volume of the space time side, distorted space-time side ?
// -------------------------------------------------------------------
template<DRT::Element::DiscretizationType space_time_distype, const int numnode_space_time>
bool XFEM::XFluidTimeInt::CheckSTSideVolume( const LINALG::Matrix<3,numnode_space_time>&  xyze_st )
{
  bool successful = true;

  const int nsd = DRT::UTILS::DisTypeToDim<space_time_distype>::dim;

  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<nsd> intpoints_stab(DRT::ELEMENTS::DisTypeToStabGaussRule<space_time_distype>::rule);

  LINALG::Matrix<nsd,1> xsi(true);

  // coordinates of the current integration point
  const double* gpcoord = (intpoints_stab.IP().qxg)[0];
  for (int idim=0;idim<nsd;idim++)
  {
    xsi(idim) = gpcoord[idim];
  }


  LINALG::Matrix<nsd,numnode_space_time> deriv(true);
  LINALG::Matrix<nsd,nsd> xjm(true);

  DRT::UTILS::shape_function_deriv1<space_time_distype>(xsi,deriv);

  xjm.MultiplyNT(deriv,xyze_st);

  double det = 0.0;
  det = xjm.Determinant();
  if ( abs(det) < 1E-14 )
  {
    successful = false;
#ifdef DEBUG_TIMINT
    IO::cout << "\t\t\t det for space-time side element (det == 0.0) => distorted flat st-element" << IO::endl;
#endif
    return successful;
  }

#ifdef DEBUG_TIMINT
  const double wquad = intpoints_stab.IP().qwgt[0];
  std::cout << "\t\t\t space-time side element volume " << wquad*det << std::endl;
#endif

  if (det < 1E-14 and det > -1E-14)
  {
#ifdef DEBUG_TIMINT
    IO::cout << "\t\t\t det for space-time side element (det == " << det << ") => distorted flat st-element" << IO::endl;
#endif
    successful = false;
    return successful;
  }


  return successful;
}


/*------------------------------------------------------------------------------------------------*
 * export data about reconstruction method to neighbor proc and receive data from previous proc schott 04/13 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFluidTimeInt::ExportMethods(
    const std::vector<Teuchos::RCP<Epetra_Vector> >&       newRowStateVectors,               /// row map based vectors w.r.t new interface position
    const Teuchos::RCP<std::set<int> >                     dbcgids                           /// set of dof gids that must not be changed by ghost penalty reconstruction
)
{

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest higher neighbour...)
  for (int dest=(myrank_+1)%numproc_;dest!=myrank_;dest=(dest+1)%numproc_) // dest is the target processor
  {
    // Initialization of sending
    DRT::PackBuffer dataSend; // vector including all data that has to be send to dest proc

    // Initialization
    int source = myrank_-(dest-myrank_); // source proc (sends (dest-myrank_) far and gets from (dest-myrank_) earlier)
    if (source<0)
      source+=numproc_;
    else if (source>=numproc_)
      source -=numproc_;

    //---------------------------------------------------------------------------------------------------------------
    //--------------------- send data to next proc and receive from previous proc -----------------------------------
    //---------------------------------------------------------------------------------------------------------------
    // send current DofSetData to next proc and receive a new map from previous proc
    {
      DRT::PackBuffer dataSend; // data to be sent

      // packing the data
      for (std::map<int, std::map<int, int > >::iterator node_it=dofset_marker_export_.begin();
           node_it!=dofset_marker_export_.end();
           node_it++)
      {
        DRT::ParObject::AddtoPack(dataSend, node_it->first);
        DRT::ParObject::AddtoPack(dataSend, node_it->second);
      }

      dataSend.StartPacking();

      // packing the data
      for (std::map<int, std::map<int, int > >::iterator node_it=dofset_marker_export_.begin();
           node_it!=dofset_marker_export_.end();
           node_it++)
      {
        DRT::ParObject::AddtoPack(dataSend, node_it->first);
        DRT::ParObject::AddtoPack(dataSend, node_it->second);
      }

      std::vector<char> dataRecv;
      sendData(dataSend,dest,source,dataRecv);

      // pointer to current position of group of cells in global std::string (counts bytes)
      std::vector<char>::size_type posinData = 0;

      // unpack received data
      while (posinData < dataRecv.size())
      {
        // unpack volumecell
        int nid = -1;                         // node id
        std::map<int,int> dofset_map;         // dofset map <nds, Method>

        // unpack reconstruction method data
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, nid);
        DRT::ParObject::ExtractfromPack(posinData,dataRecv, dofset_map);

        // distribute the received information on this proc if the info is required on this node

        // set reconstruction method on this proc only for received row node information
        int lid = dis_->NodeRowMap()->LID(nid);


        if(lid == -1) continue; // this is not a row node on this proc

        // get the node
        DRT::Node* node = dis_->gNode(nid);

        for(std::map<int,int>::iterator it=dofset_map.begin(); it!=dofset_map.end(); it++)
        {
          MarkDofs(node, it->first, newRowStateVectors, (INPAR::XFEM::XFluidTimeInt)it->second, dbcgids);
        }
      }

      dis_->Comm().Barrier(); // processors wait for each other
    }

    //---------------------------------------------------------------------------------------------------------------

  } // end loop over procs

  return;
} // end ExportMethods


/*------------------------------------------------------------------------------------------------*
 * basic function sending data to dest and receiving data from source                schott 03/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFluidTimeInt::sendData(
    DRT::PackBuffer&      dataSend,
    int&                  dest,
    int&                  source,
    std::vector<char>&    dataRecv
) const
{

  std::vector<int> lengthSend(1,0);
  lengthSend[0] = dataSend().size();
  int size_one = 1;

#ifdef DEBUG
  std::cout << "--- sending "<< lengthSend[0] << " bytes: from proc " << myrank_ << " to proc " << dest << std::endl;
#endif

  // exporter for sending
  DRT::Exporter exporter(dis_->Comm());

  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter.ISend(myrank_, dest, &(lengthSend[0]) , size_one, length_tag, req_length_data);
  // ... and receive length
  std::vector<int> lengthRecv(1,0);
  exporter.Receive(source, length_tag, lengthRecv, size_one);
  exporter.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);

  // ... and receive data
  dataRecv.clear(); dataRecv.resize(lengthRecv[0]);
  exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
  exporter.Wait(req_data);

#ifdef DEBUG
  std::cout << "--- receiving "<< lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc " << source << std::endl;
#endif
} // end sendData



// -------------------------------------------------------------------
// timint output for reconstruction methods
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::Output()
{

  int step_diff = 500;

  // output for all dofsets of nodes
  const std::string filename = IO::GMSH::GetNewFileNameAndDeleteOldFiles("TIMINT_Method", step_, step_diff, true, dis_->Comm().MyPID());
  std::ofstream gmshfilecontent(filename.c_str());
  gmshfilecontent.setf(std::ios::scientific,std::ios::floatfield);
  gmshfilecontent.precision(16);
  {
    gmshfilecontent << "View \" " << "Reconstr-Method \" {\n";

    for (int i=0; i<dis_->NumMyRowNodes(); ++i)
    {
      const DRT::Node* actnode = dis_->lRowNode(i);
      const LINALG::Matrix<3,1> pos(actnode->X());

      std::map<int,std::vector<INPAR::XFEM::XFluidTimeInt> >::const_iterator it = node_to_reconstr_method_.find(actnode->Id());

      if(it == node_to_reconstr_method_.end())
      {
        continue;
        // this is a node in the void domain
      }

      // time integration reconstruction methods for the node's different dofsets
      const std::vector<INPAR::XFEM::XFluidTimeInt> & nds_methods = it->second;

      for(size_t j=0; j<nds_methods.size(); j++ )
      {
        IO::GMSH::cellWithScalarToStream(DRT::Element::point1, (int)nds_methods[j], pos, gmshfilecontent);
      }
    }
    gmshfilecontent << "};\n";
  }

  gmshfilecontent.close();

  if(myrank_==0) IO::cout << IO::endl;

}






