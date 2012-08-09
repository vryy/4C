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


#define DEBUG_TIMINT


// -------------------------------------------------------------------
// constructor
// -------------------------------------------------------------------
XFEM::XFluidTimeInt::XFluidTimeInt(
    const Teuchos::RCP<DRT::Discretization> dis,                               /// discretization
    const Teuchos::RCP<DRT::Discretization> boundarydis,                       /// boundary discretization
    const Teuchos::RCP<XFEM::FluidWizard>   wizard_old,                        /// fluid wizard at t^n
    const Teuchos::RCP<XFEM::FluidWizard>   wizard_new,                        /// fluid wizard at t^(n+1)
    const Teuchos::RCP<XFEM::FluidDofSet>   dofset_old,                        /// dofset at t^n
    const Teuchos::RCP<XFEM::FluidDofSet>   dofset_new,                        /// dofset at t^(n+1)
    const int                               step,                              /// timestep
    const Teuchos::ParameterList&           params,                            /// parameter list
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method   /// reconstruction map for nodes and its dofsets
  ) :
  dis_(dis),
  boundarydis_(boundarydis),
  wizard_old_(wizard_old),
  wizard_new_(wizard_new),
  dofset_old_(dofset_old),
  dofset_new_(dofset_new),
  step_(step),
  params_ (params),
  reconstr_method_(reconstr_method)
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
        reconstr_counts_.insert(pair<INPAR::XFEM::XFluidTimeInt,int>(*sets,1)); // initialize counter with 1
    }
  }

  if(screenout)
  {
    cout << "\nXFEM::XFluidTimeInt::Status:\n" << std::flush;

    for(std::map<INPAR::XFEM::XFluidTimeInt,int>::iterator reconstrMethod = reconstr_counts_.begin();
        reconstrMethod != reconstr_counts_.end();
        reconstrMethod++)
    {
      cout << MapMethodEnumToString(reconstrMethod->first) << ":\t #dofsets:\t" << reconstrMethod->second << endl;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
| returns matching string for each reconstruction method   schott 07/12 |
*----------------------------------------------------------------------*/
std::string XFEM::XFluidTimeInt::MapMethodEnumToString
(
   const enum INPAR::XFEM::XFluidTimeInt term
)
{
  // length of return string is 14 due to usage in formated screen output
  switch (term)
  {
  case INPAR::XFEM::Xf_TimeInt_Copy:
    return "Copy Dofset";
    break;
  case INPAR::XFEM::Xf_TimeInt_GhostPenalty:
    return "Ghost-Penalty";
    break;
  case INPAR::XFEM::Xf_TimeInt_SemiLagrange:
    return "Semi-Lagrange";
    break;
  default :
    dserror("Cannot cope with name enum %d", term);
    return "";
    break;
  }
} // ScaTraTimIntImpl::MapTimIntEnumToString


// -------------------------------------------------------------------
// transfer standard and ghost dofs to new map as far as possible and mark dofs for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::TransferDofsToNewMap(
    const Epetra_Map&                       olddofrowmap,                     /// dof row map w.r.t old interface position
    const Epetra_Map&                       olddofcolmap,                     /// dof col map w.r.t old interface position
    vector<RCP<const Epetra_Vector> >&      oldRowStateVectors,               /// row map based vectors w.r.t old interface position
    vector<RCP<Epetra_Vector> >&            newRowStateVectors,               /// row map based vectors w.r.t new interface position
    std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >& reconstr_method, /// reconstruction map for nodes and its dofsets
    Teuchos::RCP<std::set<int> >            dbcgids                           /// set of dof gids that must not be changed by ghost penalty reconstruction
    )
{
#ifdef DEBUG_TIMINT
  const int numdofpernode = 4;
#endif

  // output map <node, vec<reconstr_method for each dofset>> (assume all vel and pressure dofs are reconstructed in the same way)
  output_reconstr_.clear();

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

#ifdef DEBUG_TIMINT
    // cout if nodehandles available at t^n and t^(n+1)
    //cout << "node: " << gid << "\t handle at tn: " << !(n_old == NULL) << "\t handle at t(n+1): " << !(n_new == NULL) << endl;
#endif

    //---------------------------------------------------------------------------------
    // switch over different cases dependent of surrounding elments are cut or not at t^n and t^(n+1)
    //     case A: surrounding elements not cut at t^n AND t^(n+1) => copy dofs
    //     case B: at least one surrounding element cut at t^(n+1), uncut elements at t^n
    //     case C: uncut elements at t^(n+1), but at least one surrounding element cut at t^n
    //     case D: surrounding elements cut at old and new timestep t^n and t^(n+1)
    //---------------------------------------------------------------------------------


    //------------------------------------------------
    {
      int num_nds = -1;

      if(n_new != NULL) num_nds = n_new->NumDofSets();
      else num_nds = 1;

      output_reconstr_[node->Id()].reserve(num_nds);
    }
    //------------------------------------------------

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


    //switch basic cases
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
      CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);

    }//end case A
    else if(n_new == NULL or unique_std_uncut_np)
    {
      if(n_old == NULL) dserror("you should call case A here");

      //---------------------------------------------------------------------------------
      // case B: at least one surrounding element cut at t^(n+1), uncut elements at t^n
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
          CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
        }
        else // case c)
        {
          dserror("case: NEWLY CREATED NODEHANDLE: just ghost dofsets at old time! structural movement more than one element near node, (here SEMILAGRANGE possible!) ", gid);
        }
      } // some elements cut
      else if(numDofSets_old == 0 )
      {
        dserror("XFLUID-TIMINIT CASE B: node %d,\t no dofset at t^n available, structural movement more than one element for node ? (here SEMILAGRANGE possible!", gid);
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
          dserror("XFLUID-TIMINIT CASE B: node %d,\t all dofset at t^n are not std-dofsets, structural movement more than one element for node? (here SEMILAGRANGE possible!", gid);
        }
        else
        {

          // copy values
          CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
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
      // case C: C: uncut elements at t^(n+1), but at least one surrounding element cut at t^n
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
          // copy values
          CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
        }
        else // case c)
        {
          dserror("XFLUID-TIMINIT CASE B: node %d,\t unique new dofset is a ghost dofset! structural movement more than one element near this node", gid);
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
          dserror("XFLUID-TIMINIT CASE B: node %d,\t first dofset is not a std dofset, structural movement more than one element for node?", gid);
        }

        // loop new dofsets
        for(std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >::const_iterator sets=dof_cellsets_new.begin();
            sets!=dof_cellsets_new.end();
            sets ++)
        {

          if(nds_new == 0) // first dofset (has been checked to be a std-dofset)
          {
            CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
          }
          else
          {
            // for newly created ghost dofsets we have to reconstruct ghost dofs
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GhostPenalty,dbcgids);
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
          cout << "XFLUID-TIMINIT CASE D: node " << gid << ",\t no corresponding dofset found at time t^n for dofset " << nds_new << " at time t^(n+1)" << endl;
#endif
          if(is_std_set_np) // std at t^(n+1)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_SemiLagrange,dbcgids);
          }
          else // ghost at t^(n+1)
          {
            MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_GhostPenalty,dbcgids);
          }
        }
        else // found set
        {
          bool is_std_set_n = Is_Std_CellSet(n_old, dof_cellsets_old[nds_old]);

#ifdef DEBUG_TIMINT
          cout << "XFLUID-TIMINIT CASE D: node " << gid << ",\t corresponding std(yes/no: " << is_std_set_n << ") dofset "
                                                 << nds_old << " found at time t^n for dofset " << nds_new << " at time t^(n+1)" << endl;
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
                // copy values
                CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
              }
              else
              {
                // semilagrange
                MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_SemiLagrange,dbcgids);
              }
            }
            else
            {
              MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_SemiLagrange,dbcgids);
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
#if(1)
              CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
#else
              MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_SemiLagrange,dbcgids);
#endif
            }
            else
            {
              dserror("case should not be called");
              // not unique ghost dofset
              MarkDofs(node, nds_new, newRowStateVectors, INPAR::XFEM::Xf_TimeInt_SemiLagrange,dbcgids);
            }
          }
          else if(!is_std_set_np) // ghost at t^(n+1)
          {
            // copy values
            CopyDofs(node, nds_new, nds_old, newRowStateVectors, oldRowStateVectors,dbcgids);
          }
        } // end found set at t^n

        nds_new++;
      } // loop new dofsets
    }//end case D (node handles at t^n and t^(n+1))

  } // loop row nodes



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
// copy dofs from old vectors to new vector for all row vectors
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::CopyDofs(DRT::Node*                         node,               /// drt node
                                   const int                          nds_new,            /// nodal dofset at w.r.t new interface
                                   const int                          nds_old,            /// the corresponding nodal dofset at w.r.t old interface
                                   vector<RCP<Epetra_Vector> >&       newRowStateVectors, /// row map based state vectors at new interface
                                   vector<RCP<const Epetra_Vector> >& oldRowStateVectors, /// row map based state vectors at old interface
                                   Teuchos::RCP<std::set<int> >       dbcgids
)
{

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

  int vec_count = 0;

  // copy values for all vectors
  for(vector<RCP<const Epetra_Vector> >::iterator it=oldRowStateVectors.begin(); it!=oldRowStateVectors.end(); it++)
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

      //cout << "copy value" << (*vec_old)[dof_lid_old] << endl;

      // set Dirichlet BC for ghost penalty reconstruction
      int gid = dofs_new[i];

      if (dbcgids != Teuchos::null)
        (*dbcgids).insert(gid);

    }

    vec_count++;
  }

  SetReconstrMethod(node, nds_new, INPAR::XFEM::Xf_TimeInt_Copy);

  // set reconstr method for output
  output_reconstr_[node->Id()].push_back(-(int)INPAR::XFEM::Xf_TimeInt_Copy-6);



  return;
}


// -------------------------------------------------------------------
// mark nodal dofs of vector w.r.t new interface position for reconstruction
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::MarkDofs(
    DRT::Node*                     node,                 /// drt node
    const int                      nds_new,              /// nodal dofset at t^(n+1)
    vector<RCP<Epetra_Vector> >&   newRowStateVectors,   /// row map based state vectors at t^(n+1)
    INPAR::XFEM::XFluidTimeInt     method,               /// reconstruction method
    Teuchos::RCP<std::set<int> >   dbcgids               /// set of dof gids that must not be changed by ghost penalty reconstruction
)
{
  // get nodal dofs for current dofset w.r.t new interface position
  std::vector<int> dofs_new;
  dofset_new_->Dof(*node, nds_new, dofs_new );

  // loop vectors
  for(vector<RCP<Epetra_Vector> >::iterator it=newRowStateVectors.begin(); it!=newRowStateVectors.end(); it++)
  {

    RCP<Epetra_Vector> vec_new = *it;

    // TODO: look at this loop -> get better startvalues, e.g. 0.0 for ghost penalty
    // set a dummy value for dofs in the vector
    for(size_t i=0; i<dofs_new.size(); i++)
    {
      int dof_lid_new = vec_new->Map().LID(dofs_new[i]);
      if(dof_lid_new == -1) dserror("new dof %d not local on this proc!", dofs_new[i]);

      (*vec_new)[dof_lid_new] = 0.0; // zero value as start value for newton for gradient reconstruction

      int gid = dofs_new[i];

      // set Dirichlet BC for ghost penalty reconstruction
      if(method != INPAR::XFEM::Xf_TimeInt_GhostPenalty)
      {
        if (dbcgids != Teuchos::null)
          (*dbcgids).insert(gid);
      }
    } // dofs
  } // state vectors

  // set the reconstruction method for current nodal dofset
  SetReconstrMethod(node, nds_new, method);

  // set reconstr method for output
  output_reconstr_[node->Id()].push_back(-((int)method)-6);
  return;
}

// -------------------------------------------------------------------
// set the reconstruction method for current nodal dofset
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::SetReconstrMethod(
    DRT::Node*                     node,                 /// drt node
    const int                      nds_new,              /// nodal dofset w.r.t new interface position
    INPAR::XFEM::XFluidTimeInt     method                /// reconstruction method
    )
{
  std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt> >::iterator it;

  // find the node in map
  it= reconstr_method_.find(node->Id());

  if(it != reconstr_method_.end())
  {
    // assume that we call this function for all nodal dofsets in the right order
    if(nds_new != (int)(it->second).size()) dserror("wrong order to set reconstruction method for all dofsets");

    (it->second).push_back(method);
//    (it->second)[nds_new] = method;
  }
  else
  {
    std::vector<INPAR::XFEM::XFluidTimeInt> vec;
    vec.push_back(method);
    reconstr_method_.insert(pair<int,std::vector<INPAR::XFEM::XFluidTimeInt> >(node->Id(), vec ));
  }
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
    GEO::CUT::Node*                                  node,     /// cut node
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp >&  cell_set  /// set of volumecells
    )
{
  // assume non-standard set
  bool is_std_set = false;

  GEO::CUT::Point* p = node->point();

  GEO::CUT::Point::PointPosition pos = p->Position();

  if(pos == GEO::CUT::Point::oncutsurface)
  {
    cout << "!!WARNING point is on cut surface, Is_Std_CellSet decision?!" << endl;
//    dserror("point is on cut surface, Is_Std_CellSet decision?!");
  }

  // at least one vc has to contain the node
  for(std::set<GEO::CUT::plain_volumecell_set>::const_iterator sets=cell_set.begin(); sets!=cell_set.end(); sets++)
  {
    // break the outer loop if at least one vc contains this point
    if(is_std_set == true) break;

    const GEO::CUT::plain_volumecell_set& set = *sets;

    for(GEO::CUT::plain_volumecell_set::const_iterator vcs=set.begin(); vcs!=set.end(); vcs++)
    {
      if((*vcs)->Contains(p))
      {
        // return if at least one vc contains this point
        is_std_set=true;
        return true;
      }
    }
  }

  return is_std_set;
}

// -------------------------------------------------------------------
// identify cellsets at time t^n with cellsets at time t^(n+1)
// -------------------------------------------------------------------
void XFEM::XFluidTimeInt::IdentifyOldSets(
    int &                                                          nds_old,            /// set identified nodal dofset at t^n
    std::vector<int> &                                             identified_sides,   /// set identified using sides (side-Ids)
    const std::vector<std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp> >&  dof_cellsets_old,   /// all dofcellsets at t^n
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_new        /// dofcellset at t^(n+1) which has to be identified
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

      // get sides involved in creation boundary cells (map<sideId,bcells>)
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

        // get sides involved in creation boundary cells (map<sideId,bcells>)
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
    cout << "Warning: found dofset at t^n not unique, set status to unfound!" << endl;
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
        cout << "point moves within side " << *(on_cut_sides_new.begin()) << " oncutsurface!" << endl;
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
    cout << "!!WARNING point is on cut surfac at time t^n and t^(n+1)" << endl;
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
    cout <<  "\t CheckChangingSide for node " << n_old->Id() << " w.r.t side " << *sides << endl;
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
      default: dserror("side-distype %s not handled", DRT::DistypeToString(side_distype).c_str());
    }

    if(!successful_check) return false; // return non-successful check

    if(node_within_Space_Time_Side)
    {

#ifdef DEBUG_TIMINT
      cout << "\t\t node " << n_new->Id()<< " changed interface side, node within space-time side element for at least one side" << endl;
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
    cout << "\t\t\t Not Successful! -> Apply small artificial translation in normal direction at t^(n+1) of " << perturbation << endl;
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
    if(successful_check) cout << "\t\t\t Successful check!" << endl;
    else                 cout << "\t\t\t Again not successful check!" << endl;
#endif
  }

  if(!successful_check) return successful_check;

  LINALG::Matrix<3,1> rst(true); // local coordinates w.r.t space time element (r,s,t !!!)

  GEO::CUT::Position<space_time_distype> pos( xyze_st, n_coord );
  pos.Compute();

  rst = pos.LocalCoordinates();

  within_space_time_side = false;

  if(pos.WithinLimits())
  {
    within_space_time_side=true;

#ifdef DEBUG_TIMINT
    cout << "\t\t\t changing side found!" << xyze_new << " \n local coords " << rst << endl;
#endif
  }

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
  LINALG::Matrix<nsd,nsd> xji(true);

  DRT::UTILS::shape_function_deriv1<space_time_distype>(xsi,deriv);

  xjm.MultiplyNT(deriv,xyze_st);

  double det = 0.0;

  try{
    det = xji.Invert(xjm);
  }
  catch (std::runtime_error) {
     // code that executes when exception-declaration is thrown
     // in the try block
#ifdef DEBUG_TIMINT
     cout << "\t\t\t det for space-time side element (det == 0.0) => distorted flat st-element" << endl;
#endif
     successful = false;
     return successful;
  };

#ifdef DEBUG_TIMINT
  const double wquad = intpoints_stab.IP().qwgt[0];
  cout << "\t\t\t space-time side element volume " << wquad*det << endl;
#endif

  if (det < 1E-14 and det > -1E-14)
  {
#ifdef DEBUG_TIMINT
    cout << "\t\t\t det for space-time side element (det == " << det << ") => distorted flat st-element" << endl;
#endif
    successful = false;
    return successful;
  }


  return successful;
}
