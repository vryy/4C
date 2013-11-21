/*!
\file xfluid_timeInt.cpp

\brief provides the xfluid timeIntegration,
       maps vectors from old interface position to vectors at new interface position,
       determines the reconstruction method for missing and unreasonable ghost and standard values

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>


#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_io/io_gmsh.H"

#include "../drt_cut/cut_boundingbox.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_point.H"
#include "../drt_cut/cut_element.H"
#include "../drt_cut/cut_volumecell.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_fluid_ele/fluid_ele_utils.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_interface.H"
#include "../drt_fluid_ele/fluid_ele_factory.H"

#include "xfem_fluidwizard.H"
#include "xfem_fluiddofset.H"

#include "xfluid_timeInt.H"

#include <iostream>


//#define DEBUG_TIMINT


// -------------------------------------------------------------------
// constructor
// -------------------------------------------------------------------
XFEM::XFluidTimeInt::XFluidTimeInt(
    const Teuchos::RCP<DRT::Discretization> dis,                               /// discretization
    const Teuchos::RCP<XFEM::FluidWizard>   wizard_old,                        /// fluid wizard at t^n
    const Teuchos::RCP<XFEM::FluidWizard>   wizard_new,                        /// fluid wizard at t^(n+1)
    const Teuchos::RCP<XFEM::FluidDofSet>   dofset_old,                        /// dofset at t^n
    const Teuchos::RCP<XFEM::FluidDofSet>   dofset_new,                        /// dofset at t^(n+1)
    const Teuchos::ParameterList&           params,                            /// parameter list
    INPAR::XFEM::XFluidTimeIntScheme        xfluid_timintapproach,             /// xfluid_timintapproch
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method,  /// reconstruction map for nodes and its dofsets
    int                                     step                               /// global time step
  ) :
  dis_(dis),
  wizard_old_(wizard_old),
  wizard_new_(wizard_new),
  dofset_old_(dofset_old),
  dofset_new_(dofset_new),
  params_ (params),
  timeint_scheme_ (xfluid_timintapproach),
  reconstr_method_(reconstr_method),
  step_(step)
  {

    myrank_  = dis->Comm().MyPID();
    numproc_ = dis->Comm().NumProc();

    return;
  } // end constructor


// -------------------------------------------------------------------
// set and print reconstruction status for nodes
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::SetAndPrintStatus(bool screenout)
{


  for(std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >::iterator node_it=reconstr_method_.begin();
      node_it!= reconstr_method_.end();
      node_it++)
  {
    std::vector<INPAR::XFEM::XFluidTimeInt> nodesets = node_it->second;

    for(std::vector<INPAR::XFEM::XFluidTimeInt>::iterator sets= nodesets.begin();
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
    return "Semi-Lagrange: STD(n+1)              ";
    break;
  case INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST:
    return "Copy Dofset:   GHOST(n) -> STD(n+1)  ";
    break;
  case INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD:
    return "Copy Dofset:   STD(n)   -> STD(n+1)  ";
    break;
  case INPAR::XFEM::Xf_TimeInt_GHOST_by_GP:
    return "Ghost-Penalty: GHOST(n+1)            ";
    break;
  case INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST:
    return "Copy Dofset:   GHOST(n) -> GHOST(n+1)";
    break;
  case INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_STD:
    return "Copy Dofset:   STD(n)   -> GHOST(n+1)";
    break;
  case INPAR::XFEM::Xf_TimeInt_undefined:
    return "undefined:                           ";
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }

  return "";
} // ScaTraTimIntImpl::MapTimIntEnumToString


// -------------------------------------------------------------------
// transfer standard and ghost dofs to new map as far as possible and mark dofs for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::TransferDofsToNewMap(
    const Epetra_Map&                       olddofrowmap,                     /// dof row map w.r.t old interface position
    const Epetra_Map&                       olddofcolmap,                     /// dof col map w.r.t old interface position
    std::vector<RCP<const Epetra_Vector> >& oldRowStateVectors,               /// row map based vectors w.r.t old interface position
    std::vector<RCP<Epetra_Vector> >&       newRowStateVectors,               /// row map based vectors w.r.t new interface position
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method, /// reconstruction map for nodes and its dofsets
    Teuchos::RCP<std::set<int> >            dbcgids                           /// set of dof gids that must not be changed by ghost penalty reconstruction
    )
{

#ifdef DEBUG_TIMINT
  const int numdofpernode = 4;
#endif


  if(oldRowStateVectors.size() != newRowStateVectors.size()) dserror("TransferDofsToNewMap: not equal number of old and new vectors");

  //------------------------------------------------------------
  // loop background fluid row nodes
  const Epetra_Map* noderowmap = dis_->NodeRowMap();
  const int numrownode = noderowmap->NumMyPoints();

  for (int lid=0; lid<numrownode; ++lid)
  {
    // get global id of a node
    int gid = noderowmap->GID(lid);

    // get the node
    DRT::Node* node = dis_->gNode(gid);

    // get cut nodes with respect to fluid wizards at t^n and t^(n+1)
    GEO::CUT::Node * n_new = wizard_new_->GetNode(gid);
    GEO::CUT::Node * n_old = wizard_old_->GetNode(gid);


    //---------------------------------------------------------------------------------
    // switch over different cases dependent of surrounding elements are cut or not at t^n and t^(n+1)
    //     case A: surrounding elements not cut at t^n AND t^(n+1) => copy dofs
    //     case B: at least one surrounding element cut at t^n, uncut elements at t^(n+1)
    //     case C: uncut elements at t^n, but at least one surrounding element cut at t^(n+1)
    //     case D: surrounding elements cut at old and new timestep t^n and t^(n+1)
    //---------------------------------------------------------------------------------


    bool unique_std_uncut_np = false;
    bool unique_std_uncut_n = false;

    // check for unique std dofset and surrounding uncut elements at t^(n+1)
    if(n_new != NULL)
    {
      const int numDofSets_new = n_new->NumDofSets(); //= dof_cellsets_new.size()

      // just one dofset at new timestep t^(n+1)
      if(numDofSets_new == 1) unique_std_uncut_np = UncutEles(n_new);
      else                    unique_std_uncut_np = false;
    }

    // check for unique std dofset and surrounding uncut elements at t^n
    if(n_old != NULL)
    {
      const int numDofSets_old = n_old->NumDofSets(); //= dof_cellsets_old.size()

      // just one dofset at new timestep t^n
      if(numDofSets_old == 1) unique_std_uncut_n = UncutEles(n_old);
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

      std::vector<int> dofs_new;
      dofset_new_->Dof(*node, nds_new, dofs_new );

#ifdef DEBUG_TIMINT
      const int numdofs_new = (int)(dofs_new.size());
      if(numdofs_new != numdofpernode)
      {
        dserror("XFLUID-TIMINIT CASE B: node %d,\t %d dofs for node without nodehandle at timestep t^(n+1) ?!", gid, numdofs_new);
      }
#endif

      //-------------------------------
      // t^(n)
      const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >& dof_cellsets_old = n_old->DofCellSets();

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
        int nds_old = FindStdDofSet(n_old, dof_cellsets_old);

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
            dserror("structural movement more than one element near node %d", n_new->Id());
          else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                   timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
          }
          else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

        }
      } // some elements cut
      else if(numDofSets_old == 0 )
      {
        // case: "XFLUID-TIMINIT CASE B: node %d,\t no dofset at t^n available,
        //        structural movement more than one element for node ? (here SEMILAGRANGE possible!", gid);

        if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
          dserror("no dofset at t^n available, structural movement more than one element for node %d", n_new->Id());
        else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                 timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
        {
          MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
        }
        else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

      }
      else
      {
        //-------------------------------------
        // more than one dofset at old timestep
        //-------------------------------------

        // check if there is a standard dofset
        int nds_old = FindStdDofSet(n_old, dof_cellsets_old);

        if( nds_old == -1 )
        {
          // all dofsets at t^n are not std-dofsets, structural movement more than one element

          if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
            dserror("all dofset at t^n are ghost-dofsets, structural movement more than one element for node %d", n_new->Id());
          else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                   timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
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
        continue; // do nothing for this node
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
      dofset_old_->Dof(*node, nds_old, dofs_old );

#ifdef DEBUG_TIMINT
      const int numdofs_old = (int)(dofs_old.size());
      if(numdofs_old != numdofpernode)
      {
        dserror("XFLUID-TIMINIT CASE C: node %d,\t %d dofs for node without nodehandle at timestep t^n ?!", gid, numdofs_old);
      }
#endif

      //-------------------------------
      // t^(n+1)
      const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >& dof_cellsets_new = n_new->DofCellSets();


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
        const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set = dof_cellsets_new[nds_new];

        if(Is_Std_CellSet(n_new, cell_set)) // case b)
        {
          // copy values or SL
          if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP or
              timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP)
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

        // REMARK: check if the first set is a std set
        int nds_new = FindStdDofSet(n_new, dof_cellsets_new);

        if( nds_new == -1 )
        {
          dserror("XFLUID-TIMINIT CASE C: node %d,\t first dofset is not a std dofset, structural movement more than one element for node?", gid);
        }

        // loop new dofsets
        for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >::const_iterator sets=dof_cellsets_new.begin();
            sets!=dof_cellsets_new.end();
            sets ++)
        {

          if(nds_new == 0) // first dofset (has been checked to be a std-dofset)
          {
            // copy values or SL
            if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP or
                timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP)
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

          nds_new ++;
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
        continue; // do nothing for this node
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
      const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >& dof_cellsets_old = n_old->DofCellSets();

      //-------------------------------
      // t^(n+1)
      const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >& dof_cellsets_new = n_new->DofCellSets();

      int nds_new = 0; // nodal dofset counter

      // loop new dofsets
      for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp > >::const_iterator sets=dof_cellsets_new.begin();
                sets!=dof_cellsets_new.end();
                sets ++)
      {

        const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set = *sets;

        //is current set at t^(n+1) std or ghost or dofset
        bool is_std_set_np = Is_Std_CellSet(n_new, cell_set);

        //-----------------------------------------------------------------
        // identify cellsets at t^n with current dofset at t^(n+1)
        int nds_old = -1; // nds number of corresponding old dofset
        std::vector<int> identified_sides;
        IdentifyOldSets(nds_old, identified_sides, dof_cellsets_old, cell_set);
        //-----------------------------------------------------------------

        if(nds_old == -1) // not found
        {
#ifdef DEBUG_TIMINT
          IO::cout << "XFLUID-TIMINIT CASE D: node " << gid << ",\t no corresponding dofset found at time t^n for dofset " << nds_new << " at time t^(n+1)" << IO::endl;
#endif
          if(is_std_set_np) // std at t^(n+1)
          {
            if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP)
              dserror("no correspoinding dofset at t^n, choose a SL-based-approach for node %d", n_new->Id());
            else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                     timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
            {
              MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
            }
            else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");
          }
          else // ghost at t^(n+1)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
          }
        }
        else // found set
        {
          bool is_std_set_n = Is_Std_CellSet(n_old, dof_cellsets_old[nds_old]);

#ifdef DEBUG_TIMINT
          IO::cout << "XFLUID-TIMINIT CASE D: node " << gid << ",\t corresponding std(yes/no: " << is_std_set_n << ") dofset "
                                                     << nds_old << " found at time t^n for dofset " << nds_new << " at time t^(n+1)" << IO::endl;
#endif

          if(is_std_set_np and is_std_set_n) // std at t^n and t^(n+1)
          {
            // check special case of tips ( change of side w.r.t identified side of volumecell at t^n and t^(n+1) )
            bool std_node_changed_side = false;
            bool successful_check = false;
            successful_check = CheckChangingSide(std_node_changed_side, identified_sides, n_old,n_new, dof_cellsets_old[nds_old], dof_cellsets_old[nds_new] );

            if(successful_check)
            {

              if(!std_node_changed_side)
              {
                // copy values or use SL

                if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP or
                   timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP)
                {
                  CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD, newRowStateVectors, oldRowStateVectors,dbcgids);
                }
                else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
                {
                  MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
                }
                else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

              }
              else
              {
                // only semilagrange (SL) reasonable

                if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP )
                  dserror("std-node changed the side, choose a SL-based-approach for node %d", n_new->Id());
                else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                         timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP )
                {
                  MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
                }
                else dserror("unknwon INPAR::XFEM::Xf_TimIntScheme");

              }
            }
            else
            {
              if(timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP )
              {
                dserror("Changing side check not successful, you should use a combined semi-lagrangean algorithm here for node %d", n_new->Id());
              }
              else if( timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP or
                       timeint_scheme_ == INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP)
              {
                MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_STD_by_SL,dbcgids);
              }
            }
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
            CopyDofs(node, nds_new, nds_old, INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST, newRowStateVectors, oldRowStateVectors,dbcgids);
//            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
          }
        } // end found set at t^n

        nds_new++;
      } // loop new dofsets
    }//end case D (node handles at t^n and t^(n+1))

  } // loop row nodes

  // export the reconstruction method information to other procs
  // this has to be done, when SL-nodes mark non-row surrounding nodes (and its some of its dofsets) as Ghost-penalty dofsets
  ExportMethods(newRowStateVectors, dbcgids);

  return;
}


// -------------------------------------------------------------------
// all surrounding elements uncut ?
// -------------------------------------------------------------------
bool XFEM::XFluidTimeInt::UncutEles(GEO::CUT::Node* n)
{
  // surrounding elements
  const GEO::CUT::plain_element_set& adj_eles = n->Elements();

  // assume all elements are uncut
  bool all_eles_uncut = true;

  // loop surrounding elements
  for(GEO::CUT::plain_element_set::const_iterator eles = adj_eles.begin(); eles!=adj_eles.end(); eles++)
  {
    GEO::CUT::Element* e = *eles;
    if(e->IsCut())
    {
      all_eles_uncut = false;
      return all_eles_uncut; // at least one element is cut
    }
  }

  return all_eles_uncut;
}


// -------------------------------------------------------------------
// find all ghost dofsets around this node and its std-dofset
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::FindSurroundingGhostDofsets(
    std::map<int,std::vector<int> >& ghostDofsets,   /// surrounding ghost dofsets to be filled
    DRT::Node*                       node,           /// node
    int                              nds_new         /// dofset of node used for finding the surrounding ghost dofs
    )
{
  if(nds_new != 0 ) dserror("do you really want to find ghost dofsets surrounding a non-std node?");

  GEO::CUT::Node * n = wizard_new_->GetNode(node->Id());

  if(n == NULL)
  {
    // it can happen that the node was std at t^n, however cut elements around it (nodehandle available at t^n and SL called)
    // and at new time t^(n+1) it is a std-node without a nodehandle, then there are no surrounding ghost dofs at t^(n+1)
    return;
  }

  const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >& dof_cellsets = n->DofCellSets();

  // get the corresponding cellset
  const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> & cellset = dof_cellsets[nds_new];

  // get for each plain_volumecell_set of the surrounding elements the ghost dofs

  // loop the surrounding elements
  for(std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>::const_iterator e_vcset=cellset.begin();
      e_vcset != cellset.end();
      e_vcset++)
  {
    const GEO::CUT::plain_volumecell_set & vcs= *e_vcset;

    // get the element, ask the first vc
    GEO::CUT::Element* e = vcs[0]->ParentElement();
    const std::vector<GEO::CUT::Node*> nodes = e->Nodes();

    const std::vector<int> & nds = vcs[0]->NodalDofSet();

    // which dofset is the corresponding to the current vc-connection and is a ghost dofset?
    for(std::vector<GEO::CUT::Node*>::const_iterator n_it=nodes.begin();
        n_it!=nodes.end(); n_it++)
    {
      // check if the neighbored node is a ghost node w.r.t to the cellset of the SL-node
      // if ghost node w.r.t. the cellset then add it to ghostDofsets with the corresponding nds-number
      if(!Is_Std_CellSet(*n_it, cellset))
      {
        std::map<int, std::vector<int> >::iterator map_it = ghostDofsets.find((*n_it)->Id());
        if(map_it==ghostDofsets.end())
        {
          std::vector<int> tmp_vec;

          tmp_vec.push_back(nds[n_it-nodes.begin()]);
          ghostDofsets.insert(std::pair<int, std::vector<int> >((*n_it)->Id(),tmp_vec));
        }
        else
        {
          map_it->second.push_back(nds[n_it-nodes.begin()]);
        }
      }
    }
  }
}



// -------------------------------------------------------------------
// copy dofs from old vectors to new vector for all row vectors
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::CopyDofs(DRT::Node*                              node,               /// drt node
                                   const int                               nds_new,            /// nodal dofset at w.r.t new interface
                                   const int                               nds_old,            /// the corresponding nodal dofset at w.r.t old interface
                                   INPAR::XFEM::XFluidTimeInt              method,               /// reconstruction method
                                   std::vector<RCP<Epetra_Vector> >&       newRowStateVectors, /// row map based state vectors at new interface
                                   std::vector<RCP<const Epetra_Vector> >& oldRowStateVectors, /// row map based state vectors at old interface
                                   Teuchos::RCP<std::set<int> >            dbcgids             /// set of DBC global ids
)
{

  if(method != INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_GHOST and
     method != INPAR::XFEM::Xf_TimeInt_GHOST_by_COPY_from_STD   and
     method != INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_GHOST   and
     method != INPAR::XFEM::Xf_TimeInt_STD_by_COPY_from_STD) dserror("don't call CopyDofs for non-copy reconstruction type");


  std::vector<int> dofs_old;
  std::vector<int> dofs_new;

  dofset_new_->Dof(*node, nds_new, dofs_new );
  dofset_old_->Dof(*node, nds_old, dofs_old );

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
  for(std::vector<RCP<const Epetra_Vector> >::iterator it=oldRowStateVectors.begin(); it!=oldRowStateVectors.end(); it++)
  {

    RCP<Epetra_Vector>       vec_new = newRowStateVectors[vec_count];
    RCP<const Epetra_Vector> vec_old = oldRowStateVectors[vec_count];

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

  return;
}


// -------------------------------------------------------------------
// mark nodal dofs of vector w.r.t new interface position for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::MarkDofs(
    DRT::Node*                          node,                 /// drt node
    const int                           nds_new,              /// nodal dofset at t^(n+1)
    std::vector<RCP<Epetra_Vector> >&   newRowStateVectors,   /// row map based state vectors at t^(n+1)
    INPAR::XFEM::XFluidTimeInt          method,               /// reconstruction method
    Teuchos::RCP<std::set<int> >        dbcgids               /// set of dof gids that must not be changed by ghost penalty reconstruction
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
    // map of global nid and the corresponding ghost-dofset-number
    std::map<int,std::vector<int > > ghostDofsets;
    FindSurroundingGhostDofsets(ghostDofsets, node, nds_new);

    for(std::map<int,std::vector<int > >::iterator it=ghostDofsets.begin(); it!=ghostDofsets.end(); it++)
    {
      // mark the ghost dofset in case that the node is a row node on this proc
      // otherwise mark it for export at the end of TransferDofsToNewMap
      int nid=it->first;
      std::vector<int>& nds = it->second;
      for(std::vector<int>::iterator nds_it=nds.begin(); nds_it!=nds.end(); nds_it++)
      {
        int dofset = *nds_it;
        if(dis_->NodeRowMap()->LID(nid) != -1)
        {
          DRT::Node* n = dis_->gNode(nid);
          // TODO:
          // uncomment when no ghost-values shall be copied instead of reconstructed
          MarkDofs(n, dofset, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP,dbcgids);
        }
        else
        {
          //TODO:
          // uncomment when no ghost-values shall be copied instead of reconstructed
          MarkDofsForExport(nid, dofset, INPAR::XFEM::Xf_TimeInt_GHOST_by_GP);
        }
      }
    }
  }
  //-------------------------------------

  // get nodal dofs for current dofset w.r.t new interface position
  std::vector<int> dofs_new;
  dofset_new_->Dof(*node, nds_new, dofs_new );

  // loop vectors
  for(std::vector<RCP<Epetra_Vector> >::iterator it=newRowStateVectors.begin(); it!=newRowStateVectors.end(); it++)
  {

    RCP<Epetra_Vector> vec_new = *it;

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
    int                        nid,         /// node id
    int                        dofset,      /// ghost dofset number
    INPAR::XFEM::XFluidTimeInt method       /// reconstruction method used for marking the nodal dofset
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
    DRT::Node*                     node,                 /// drt node
    const int                      nds_new,              /// nodal dofset w.r.t new interface position
    INPAR::XFEM::XFluidTimeInt     method                /// which type of reconstruction method
    )
{
  GEO::CUT::Node* n = wizard_new_->GetNode(node->Id());

  int numdofsets = -1;
  if(n!=NULL){ numdofsets = n->NumDofSets();}
  else numdofsets = 1; // one std dofset

  // find the node in map
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >::iterator it;
  it= reconstr_method_.find(node->Id());

  if(it != reconstr_method_.end())
  {
    // assume that we call this function for all nodal dofsets in the right order
    if(nds_new >= numdofsets) dserror("no valid nds %d for node %d", nds_new, n->Id());

    //---------------------------------------------------------
    // check if the current set method can be overwritten

    // do not overwrite Ghost-penalty reconstruction method which have been already marked by surrounding SL-dofsets
    if(it->second[nds_new] == INPAR::XFEM::Xf_TimeInt_GHOST_by_GP) return false;

    // do not overwrite already values, expect that method is Xf_TimeInt_GHOST_by_GP
    if(it->second[nds_new] != INPAR::XFEM::Xf_TimeInt_undefined and   // overwriting
       it->second[nds_new] != method and                              // methods not equal
       method != INPAR::XFEM::Xf_TimeInt_GHOST_by_GP                  // overwriting with non-GP approach
       )
    {
      dserror("inconsistency in reconstruction method, why do want to replace reconstruction method %d with %d for node=%d and nds=%d",  it->second[nds_new], method, node->Id(), nds_new );
    }
    //---------------------------------------------------------

    (it->second)[nds_new] = method;
    return true;
  }
  else
  {
    std::vector<INPAR::XFEM::XFluidTimeInt> vec(numdofsets,INPAR::XFEM::Xf_TimeInt_undefined); // initialize with undefined status
    vec[nds_new] = method;
    reconstr_method_.insert(std::pair<int,std::vector<INPAR::XFEM::XFluidTimeInt> >(node->Id(), vec ));
    return true;
  }

  return false;
}

// -------------------------------------------------------------------
// find the standard dofset, return the dofset number of std dofset
// -------------------------------------------------------------------
int XFEM::XFluidTimeInt::FindStdDofSet(
    GEO::CUT::Node*                                                               node,          /// cut node
    const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp > >& dof_cell_sets  /// dofcellsets of node
    )
{
  int count = 0;

  for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp > >::const_iterator it = dof_cell_sets.begin();
      it!= dof_cell_sets.end();
      it++)
  {
    if(Is_Std_CellSet(node,*it)) return count;

    count ++;
  }

  return -1;
}

// -------------------------------------------------------------------
// is this node a standard or ghost node w.r.t current set
// -------------------------------------------------------------------
bool XFEM::XFluidTimeInt::Is_Std_CellSet(
    GEO::CUT::Node*                                                  node,     /// cut node
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp >&  cell_set  /// set of volumecells
    )
{
  GEO::CUT::Point* p = node->point();

  // at least one vc has to contain the node
  for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator sets=cell_set.begin(); sets!=cell_set.end(); sets++)
  {
    const GEO::CUT::plain_volumecell_set& set = *sets;

    for(GEO::CUT::plain_volumecell_set::const_iterator vcs=set.begin(); vcs!=set.end(); vcs++)
    {
      if((*vcs)->Contains(p))
      {
        // return if at least one vc contains this point
        return true;
      }
    }
  }

  return false;
}

// -------------------------------------------------------------------
// identify cellsets at time t^n with cellsets at time t^(n+1)
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::IdentifyOldSets(
    int &                                                                         nds_old,            /// set identified nodal dofset at t^n
    std::vector<int> &                                                            identified_sides,   /// set identified using sides (side-Ids)
    const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >&  dof_cellsets_old,   /// all dofcellsets at t^n
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>&                cell_set_new        /// dofcellset at t^(n+1) which has to be identified
    )
{
  // set of side-ids involved in cutting the current connection of volumecells at t^(n+1)
  std::map<int, std::vector<GEO::CUT::BoundaryCell*> >  bcells_new;

  //--------------------------------------------------------
  // t^(n+1)
  // get all side-ids w.r.t to all volumecells contained in current new set around the current node
  for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator adj_eles = cell_set_new.begin(); adj_eles!=cell_set_new.end(); adj_eles++)
  {
    const GEO::CUT::plain_volumecell_set ele_vc = *adj_eles;

    for(GEO::CUT::plain_volumecell_set::const_iterator vcs=ele_vc.begin(); vcs!=ele_vc.end(); vcs++)
    {
      GEO::CUT::VolumeCell* vc = *vcs;

      // get sides involved in creation boundary cells (std::map<sideId,bcells>)
      vc->GetBoundaryCells(bcells_new);
    }
  }


  //--------------------------------------------------------
  // t^n
  int num_found_old_cellsets = 0; // how many dofsets found

  int set_number = 0;

  //check each old dofset for identification with new dofset
  for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >::const_iterator old_sets=dof_cellsets_old.begin();
      old_sets!=dof_cellsets_old.end();
      old_sets++)
  {

    bool identified_set = false;

    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> & old_set = *old_sets;

    // for each set a new map of bcs
    std::map<int, std::vector<GEO::CUT::BoundaryCell*> >  bcells_old;

    // get all side-ids w.r.t to all volumecells contained in this set around the current node
    for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator adj_eles = old_set.begin(); adj_eles!=old_set.end(); adj_eles++)
    {
      const GEO::CUT::plain_volumecell_set ele_vc = *adj_eles;

      for(GEO::CUT::plain_volumecell_set::const_iterator vcs=ele_vc.begin(); vcs!=ele_vc.end(); vcs++)
      {
        GEO::CUT::VolumeCell* vc = *vcs;

        // get sides involved in creation boundary cells (std::map<sideId,bcells>)
        vc->GetBoundaryCells(bcells_old);
      }
    }

    //--------------------------------------------------------------
    // check if identification of sets is possible
    // -> check if any side is involved in both cuts (find side-Id involved at time t^(n+1) in set of involved side-Ids at t^n)
    //--------------------------------------------------------------
    for(std::map<int, std::vector<GEO::CUT::BoundaryCell*> >::iterator old_sides=bcells_old.begin();
        old_sides!=bcells_old.end();
        old_sides++)
    {
      if(bcells_new.find((old_sides->first)) != bcells_new.end())
      {
        identified_set = true;
        identified_sides.push_back(old_sides->first);
      }
    }

    if(identified_set == true)
    {
      nds_old = set_number;
      num_found_old_cellsets++;
    }

    // next possible set
    set_number++;
  }

  if(num_found_old_cellsets > 1)
  {
#ifdef DEBUG_TIMINT
    IO::cout << "Warning: found dofset at t^n not unique, set status to unfound!" << IO::endl;
#endif
    nds_old = -1;
  }



  return;
}

// -------------------------------------------------------------------
// check if the node has changed the side w.r.t identified sides at t^n and t^(n+1)
// return if check was successful
// -------------------------------------------------------------------
bool XFEM::XFluidTimeInt::CheckChangingSide(
    bool&                                                          changed_side,        /// did the node change the side ?
    std::vector<int> &                                             identified_sides,    /// side Id of identified side
    GEO::CUT::Node *                                               n_old,               /// node w.r.t to old wizard
    GEO::CUT::Node *                                               n_new,               /// node w.r.t to new wizard
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_old,        /// dofcellset at t^n
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_new         /// dofcellset at t^(n+1)
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

  if(   (pos_old == GEO::CUT::Point::inside  and pos_new == GEO::CUT::Point::outside)
     or (pos_old == GEO::CUT::Point::outside and pos_new == GEO::CUT::Point::inside )
     )
  {
    // changed from inside->outside or outside->inside
    changed_side = true;
    return true;
  }
  else if(     (pos_old == GEO::CUT::Point::outside and pos_new == GEO::CUT::Point::outside)
            or (pos_old == GEO::CUT::Point::inside  and pos_new == GEO::CUT::Point::inside)
         )
  {
    // continue with check based on space time sides
    // REMARK: the point could have moved through an inside/outside volumecell

    // do nothing and continue with space-time side check
  }
  else if(pos_old == GEO::CUT::Point::oncutsurface and pos_new == GEO::CUT::Point::inside)
  {
    //dserror("point moved from %d position at t^n to %d at t^(n+1), is this a side change?", pos_old, pos_new);
    changed_side = true;
    return true;
  }
  else if(pos_old == GEO::CUT::Point::inside and pos_new == GEO::CUT::Point::oncutsurface)
  {
    //dserror("point moved from %d position at t^n to %d at t^(n+1), is this a side change?", pos_old, pos_new);
    changed_side = true;
    return true;
  }
  else if(pos_old == GEO::CUT::Point::oncutsurface and pos_new == GEO::CUT::Point::outside)
  {
    //dserror("point moved from %d position at t^n to %d at t^(n+1), is this a side change?", pos_old, pos_new);
    changed_side = false;
    return true;
  }
  else if(pos_old == GEO::CUT::Point::outside and pos_new == GEO::CUT::Point::oncutsurface)
  {
    //dserror("point moved from %d position at t^n to %d at t^(n+1), is this a side change?", pos_old, pos_new);
    changed_side = false;
    return true;
  }
  else if( (pos_old == GEO::CUT::Point::oncutsurface and pos_new == GEO::CUT::Point::oncutsurface) )
  {
    // first case: point slides on the cut surface (within one side or from 'within the side' to point or edge and vice versa)
    // second case: point moves from one side to another side (-> decide if the sides are neighbors and so on)

    const GEO::CUT::plain_side_set & cut_sides_old  = p_old->CutSides();
    const GEO::CUT::plain_side_set & cut_sides_new  = p_new->CutSides();

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

    changed_side = false;
    return true;
  }


  if(pos_old == GEO::CUT::Point::oncutsurface and pos_new == GEO::CUT::Point::oncutsurface)
  {
#ifdef DEBUG_TIMINT
    IO::cout << "!!WARNING point is on cut surfac at time t^n and t^(n+1)" << IO::endl;
#endif

    const GEO::CUT::plain_side_set & cut_sides_old  = p_old->CutSides();
    const GEO::CUT::plain_side_set & cut_sides_new  = p_new->CutSides();

    if(cut_sides_old.size() != 1) dserror("more than one cut side");
    if(cut_sides_new.size() != 1) dserror("more than one cut side");

  }
  else if(pos_old == GEO::CUT::Point::oncutsurface or pos_new == GEO::CUT::Point::oncutsurface)
  {
    dserror("point is on cut surface at least at one time t(n) or t(n+1), CheckChangingSide in the right way?");
  }



  LINALG::Matrix<3,1> n_coord(true);
  n_old->Coordinates(&n_coord(0,0));

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

    int sid = *sides;
    GEO::CUT::SideHandle* side_old = wizard_old_->GetCutSide(sid, 0);
    GEO::CUT::SideHandle* side_new = wizard_new_->GetCutSide(sid, 0);

    if(side_old == NULL or side_new == NULL) dserror("no sidehandles available for side %d", sid);


    bool node_within_Space_Time_Side = false;

    if(side_old->Shape() != side_new->Shape()) dserror("the same side has different shapes at tn and t(n+1) -> not possible!");

    DRT::Element::DiscretizationType side_distype = side_old->Shape();

    switch(side_distype)
    {
      case DRT::Element::tri3:
      {
        successful_check = WithinSpaceTimeSide<DRT::Element::tri3,DRT::Element::wedge6>(node_within_Space_Time_Side,side_old, side_new, n_coord);
        break;
      }
      case DRT::Element::quad4:
      {
        successful_check = WithinSpaceTimeSide<DRT::Element::quad4,DRT::Element::hex8>(node_within_Space_Time_Side,side_old, side_new, n_coord);
        break;
      }
      default: dserror("side-distype %s not handled", DRT::DistypeToString(side_distype).c_str()); break;
    }

    if(!successful_check) return false; // return non-successful check

    if(node_within_Space_Time_Side)
    {

#ifdef DEBUG_TIMINT
      IO::cout << "\t\t node " << n_new->Id()<< " changed interface side, node within space-time side element for at least one side" << IO::endl;
#endif

      changed_side = true;
      break; // stop loop over sides
    }

  }

  return successful_check;
}


// -------------------------------------------------------------------
// check if the node is within the space time side
// -------------------------------------------------------------------
template<DRT::Element::DiscretizationType side_distype, DRT::Element::DiscretizationType space_time_distype>
bool XFEM::XFluidTimeInt::WithinSpaceTimeSide(
    bool&                   within_space_time_side,    /// within the space time side
    GEO::CUT::SideHandle*   side_old,                  /// side w.r.t old interface
    GEO::CUT::SideHandle*   side_new,                  /// side w.r.t new interface
    LINALG::Matrix<3,1>&    n_coord                    /// node coodinates
)
{

  const int numnode_space_time = DRT::UTILS::DisTypeToNumNodePerEle<space_time_distype>::numNodePerElement;

  const int numnode_side = DRT::UTILS::DisTypeToNumNodePerEle<space_time_distype>::numNodePerElement/2;

  Epetra_SerialDenseMatrix xyze_old;
  Epetra_SerialDenseMatrix xyze_new;

  side_old->Coordinates(xyze_old);
  side_new->Coordinates(xyze_new);


  // space time side coordinates
  LINALG::Matrix<3,numnode_space_time> xyze_st;

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
    else if(numnode_side == 4)
    {
      xi_side(0) = 0.0;
      xi_side(1) = 0.0;
    }
    else dserror("unknown side type with %d nodes", numnode_side);

    // Initialization
    LINALG::Matrix<2,numnode_side> deriv(true);      // derivatives dr, ds

    LINALG::Matrix<3,2> derxy (true);
    LINALG::Matrix<3,1> dx_dr (true);
    LINALG::Matrix<3,1> dx_ds (true);

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


  GEO::CUT::Position<space_time_distype> pos( xyze_st, n_coord );
  within_space_time_side = pos.Compute();

#ifdef DEBUG_TIMINT
  LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t space time element (r,s,t !!!)
  rst = pos.LocalCoordinates();

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
bool XFEM::XFluidTimeInt::CheckSTSideVolume( LINALG::Matrix<3,numnode_space_time>&  xyze_st )
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
    std::vector<RCP<Epetra_Vector> >&       newRowStateVectors,               /// row map based vectors w.r.t new interface position
    Teuchos::RCP<std::set<int> >            dbcgids                           /// set of dof gids that must not be changed by ghost penalty reconstruction
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

      std::map<int,std::vector<INPAR::XFEM::XFluidTimeInt> >::iterator it = reconstr_method_.find(actnode->Id());

      if(it == reconstr_method_.end())
      {
        continue;
        // this is a node in the void domain
      }

      // time integration reconstruction methods for the node's different dofsets
      std::vector<INPAR::XFEM::XFluidTimeInt>& nds_methods = it->second;

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



