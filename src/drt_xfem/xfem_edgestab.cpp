/*----------------------------------------------------------------------*/
/*! \file

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented,
continuous interior penalty) scheme

\level 1

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
*/
/*----------------------------------------------------------------------*/



#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret_faces.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_cutwizard.H"


#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_intfaces_calc.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_inpar/inpar_fluid.H"

// Needed for material check.
#include "../drt_xfem/xfem_utils.H"
// Needed for safety check.
#include "../drt_mat/material.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"

#include "xfem_edgestab.H"

/*----------------------------------------------------------------------*
 |  prepares edge based stabilization and ghost penaly in case of XFEM  |
 |  and calls evaluate routine                             schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::EvaluateEdgeStabGhostPenalty(
    Teuchos::ParameterList& eleparams,                ///< element parameter list
    Teuchos::RCP<DRT::Discretization> discret,        ///< discretization
    DRT::ELEMENTS::FluidIntFace* faceele,             ///< face element
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector,         ///< systemvector
    Teuchos::RCP<GEO::CutWizard> wizard,              ///< cut wizard
    bool include_inner,        ///< stabilize also facets with inside position
    bool include_inner_faces,  ///< stabilize also faces with inside position if possible
    bool gmsh_eos_out          ///< stabilization gmsh output
)
{
  //====================================================================================================
  // implementation of edge-based fluid stabilization and ghost penalty
  //====================================================================================================

  // EDGE-based fluid stabilization and EDGE-based Ghost-penalty stabilization
  // REMARK: the current implementation of edge-based stabilization is based on the
  // DiscretizationXFEM extension
  //         using additional information about faces between two volume elements
  // * fluid stabilization has to be integrated for all internal faces
  // * ghost penalty has to be integrated if there is at least one cut element
  //   (because all faces between two elements for that at least one element is cut by the interface
  //   has to be stabilized) NOTE: the limit case that a cut side just touches a node or if the cut
  //   side touches an element side completely, the check
  //         e->IsIntersected() returns false and we do not stabilize the face.
  //         This avoids over-stabilization as it e.g. does not switch on the ghost-penalties for
  //         standard FEM situations. by a more appropriate check which tells you if the neighboring
  //         volumecell is equal to the element itself.
  //   NOTE: it might be helpful and might lead to better results when weighting ghost-penalties by
  //   e.g. volume-fractions.
  //         In that case it still has to be guaranteed not to loose coercivity. To guarantee weak
  //         consistency the scalings have to be bounded by h. Such scalings are not available yet.
  //
  // we distinguish different stabilization cases
  //  1. the master element and slave element (connected via current side)
  //     do not have an elementhandle (standard fluid case)
  //     -> standard fluid stabilization
  //                               => EOS(fluid): YES         GHOST-PENALTY: NO
  //  2. element handles for both parent elements
  //     -> stabilization for each facet and corresponding volumecells of parent elements
  //                               => EOS(fluid): YES         GHOST-PENALTY: Yes (if at least one
  //                               parent element is cut)
  //                                                                         NO  (if both parent
  //                                                                         elements are uncut)
  //  3. just one elementhandle available (at limit of bounding box)
  //     -> stabilization for each facet and corresponding volumecells of parent elements
  //                               => EOS(fluid): YES         GHOST-PENALTY: Yes (if at least one
  //                               parent element is cut)
  //                                                                         NO  (if both parent
  //                                                                         elements are uncut)


  Teuchos::RCP<DRT::DiscretizationFaces> xdiscret =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


  // get the parent fluid elements
  DRT::ELEMENTS::Fluid* p_master = faceele->ParentMasterElement();
  DRT::ELEMENTS::Fluid* p_slave = faceele->ParentSlaveElement();

  // get corresponding element handles if available
  GEO::CUT::ElementHandle* p_master_handle = wizard->GetElement(p_master);
  GEO::CUT::ElementHandle* p_slave_handle = wizard->GetElement(p_slave);

  size_t p_master_numnode = p_master->NumNode();
  size_t p_slave_numnode = p_slave->NumNode();

  // get the parent element
  int p_master_id = p_master->Id();

  std::vector<int> nds_master;
  nds_master.reserve(p_master_numnode);

  std::vector<int> nds_slave;
  nds_slave.reserve(p_slave_numnode);

  INPAR::XFEM::FaceType face_type;

  int num_edgestab = 0;      // how often to stabilize this face for edgebased stabilizations
  int num_ghostpenalty = 0;  // how often to stabilize this face for ghost penalty stabilizations


  // Provide material at both sides:
  Teuchos::RCP<MAT::Material> matptr_m;
  Teuchos::RCP<MAT::Material> matptr_s;
  matptr_m = p_master->Material();
  matptr_s = p_slave->Material();

  //------------------------------------------------------------------------------
  // simplest case: no element handles for both parent elements
  // two uncut elements / standard fluid case
  // problems cut with levelset will not enter here!
  //------------------------------------------------------------------------------
  if (p_master_handle == NULL and p_slave_handle == NULL)
  {
    num_edgestab++;

    if (matptr_m->MaterialType() == INPAR::MAT::m_matlist)
      dserror("The edgebased algo can not handle matlist at the moment, for this entry!");

    face_type = INPAR::XFEM::face_type_std;

    {
      TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");

      for (size_t i = 0; i < p_master_numnode; i++) nds_master.push_back(0);

      for (size_t i = 0; i < p_slave_numnode; i++) nds_slave.push_back(0);
    }

    //--------------------------------------------------------------------------------------------

    // call evaluate and assemble routine
    AssembleEdgeStabGhostPenalty(eleparams, face_type, faceele, matptr_m, matptr_s, nds_master,
        nds_slave, *xdiscret, systemmatrix, systemvector);

    //--------------------------------------------------------------------------------------------
  }
  //------------------------------------------------------------------------------
  // second case: element handles for both parent elements
  // two elements that are maybe cut
  //------------------------------------------------------------------------------
  else if (p_master_handle != NULL and p_slave_handle != NULL)
  {
    // linear elements
    if (p_master->Shape() == DRT::Element::hex8 or p_master->Shape() == DRT::Element::tet4 or
        p_master->Shape() == DRT::Element::wedge6 or p_master->Shape() == DRT::Element::pyramid5)
    {
      GEO::CUT::SideHandle* side = GetFace(faceele, wizard);

      //-------------------------------- loop facets of this side -----------------------------
      // facet of current side
      std::vector<GEO::CUT::Facet*> facets;
      side->Facets(facets);

      if (facets.size() == 0) dserror("there is no facet between two elements with elementhandle!");

      // each facet should have 2 volumecells
      for (std::vector<GEO::CUT::Facet*>::const_iterator f = facets.begin(); f != facets.end(); f++)
      {
        if ((*f)->Position() == GEO::CUT::Point::outside or
            ((*f)->Position() == GEO::CUT::Point::inside and
                (include_inner || include_inner_faces)))
        {
          GEO::CUT::plain_volumecell_set vcs = (*f)->Cells();

          // how many volumecells found?
          if (vcs.size() ==
              2)  // standard XFEM case (facet between two vcs of two neighbouring cut elements
          {
            GEO::CUT::plain_volumecell_set::iterator vc_it = vcs.begin();

            GEO::CUT::VolumeCell* vc1 = *(vc_it);
            vc_it++;
            GEO::CUT::VolumeCell* vc2 = *(vc_it);


            // get the parent element
            int vc_ele1_id = vc1->ParentElement()->Id();
            int vc_ele2_id = vc2->ParentElement()->Id();

            bool all_dofs = (facets.size() == 1 && include_inner_faces);
            if ((*f)->Position() == GEO::CUT::Point::outside || include_inner)
            {
              //------------------------ create nodal dof sets
              TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");
              // which element is the parent element
              if (vc_ele1_id == p_master_id)
              {
                nds_master = vc1->NodalDofSet();
                nds_slave = vc2->NodalDofSet();
              }
              else if (vc_ele2_id == p_master_id)
              {  // switch ele 1 <-> ele 2
                nds_master = vc2->NodalDofSet();
                nds_slave = vc1->NodalDofSet();
              }
              else
                dserror("no element (ele1 and ele2) is the parent element!!! WHY?");
            }
            else if ((*f)->Position() == GEO::CUT::Point::inside && all_dofs)
            {
              for (uint n = 0; n < vc1->ParentElement()->Nodes().size(); ++n)
              {
                if (!vc1->ParentElement()->Nodes()[n]->NodalDofSets().size())
                {
                  all_dofs = false;
                  break;
                }
              }
              if (all_dofs)
                for (uint n = 0; n < vc2->ParentElement()->Nodes().size(); ++n)
                {
                  if (!vc2->ParentElement()->Nodes()[n]->NodalDofSets().size())
                  {
                    all_dofs = false;
                    break;
                  }
                }
              if (all_dofs)
              {
                nds_master.clear();
                nds_slave.clear();
                if (vc1->ParentElement()->NumNodes() == vc2->ParentElement()->NumNodes())
                {
                  //------------------------ create nodal dof sets
                  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");
                  for (uint n = 0; n < vc2->ParentElement()->Nodes().size(); ++n)
                  {
                    nds_master.push_back(0);
                    nds_slave.push_back(0);
                  }
                }
                else
                  dserror("Number of Nodes different between Master and Slave Element!");
              }
            }

            if ((*f)->Position() == GEO::CUT::Point::inside && !include_inner && !all_dofs)
              continue;
            //------------------------

            num_edgestab++;

            // at least one element has to be cut
            if (p_master_handle->IsIntersected() or p_slave_handle->IsIntersected())
            {
              num_ghostpenalty++;

              face_type = INPAR::XFEM::face_type_ghost_penalty;
            }
            else
              face_type = INPAR::XFEM::face_type_std;

            XFEM::UTILS::GetVolumeCellMaterial(p_master, matptr_m, (*f)->Position());
            XFEM::UTILS::GetVolumeCellMaterial(p_slave, matptr_s, (*f)->Position());

            //--------------------------------------------------------------------------------------------

            // call evaluate and assemble routine
            AssembleEdgeStabGhostPenalty(eleparams, face_type, faceele, matptr_m, matptr_s,
                nds_master, nds_slave, *xdiscret, systemmatrix, systemvector);

            //--------------------------------------------------------------------------------------------
          }
          else if (vcs.size() == 1)
          {
            dserror("just one vcs reasonable?! face %d", faceele->Id());
          }
        }  // facet outside or (inside and include_inner)
        else if ((*f)->Position() == GEO::CUT::Point::undecided)
        {
          dserror("the position of this facet is undecided, how to stabilize???");
        }
        else if ((*f)->Position() == GEO::CUT::Point::oncutsurface)
        {
#ifdef DEBUG
          std::cout << "the position of this facet of face " << faceele->Id()
                    << " is oncutsurface, we do not stabilize it!!! " << std::endl;
#endif
          // if a facet lies oncutsurface, then there is only one neighbor, we do not stabilize this
          // facet REMARK: in case of one part of the facet is physical and the other part lies on
          // cutsurface,
          //         then the physical part is stabilized via another facet lying on the same fluid
          //         element's side
        }
        else
        {
          // facet is inside!
          face_type = INPAR::XFEM::face_type_ghost;
        }

      }  // loop facets
    }    // if linear elements
    else if (p_master->Shape() == DRT::Element::hex20 or p_master->Shape() == DRT::Element::hex27 or
             p_master->Shape() == DRT::Element::tet10 or p_master->Shape() == DRT::Element::wedge15)
    {
      GEO::CUT::SideHandle* side = GetFace(faceele, wizard);  // the side of the quadratic element
      //-------------------------------- loop facets of this side -----------------------------
      // facet of current side
      std::vector<GEO::CUT::Facet*> facets;

      side->Facets(facets);  // all facets of this quadratic element side
      if (facets.size() == 0) dserror("there is no facet between two elements with elementhandle!");
      // each facet should have 2 volumecells
      std::vector<std::vector<int>> all_used_nds_master;
      std::vector<std::vector<int>> all_used_nds_slave;
      for (std::vector<GEO::CUT::Facet*>::const_iterator f = facets.begin(); f != facets.end(); f++)
      {
        if ((*f)->Position() == GEO::CUT::Point::outside or
            ((*f)->Position() == GEO::CUT::Point::inside and include_inner))
        {
          GEO::CUT::plain_volumecell_set vcs = (*f)->Cells();
          // how many volumecells found?
          if (vcs.size() ==
              2)  // standard XFEM case (facet between two vcs of two neighbouring cut elements
          {
            //------------------------ create nodal dof sets
            {
              TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");

              GEO::CUT::plain_volumecell_set::iterator vc_it = vcs.begin();

              GEO::CUT::VolumeCell* vc1 = *(vc_it);
              vc_it++;
              GEO::CUT::VolumeCell* vc2 = *(vc_it);

              // get the parent element
              int vc_ele1_id = vc1->ParentElement()->GetParentId();
              int vc_ele2_id = vc2->ParentElement()->GetParentId();

              // which element is the parent element
              if (vc_ele1_id == p_master_id)
              {
                nds_master = vc1->NodalDofSet();
                nds_slave = vc2->NodalDofSet();
              }
              else if (vc_ele2_id == p_master_id)
              {  // switch ele 1 <-> ele 2
                nds_master = vc2->NodalDofSet();
                nds_slave = vc1->NodalDofSet();
              }
              else
              {
                dserror("no element (ele1 and ele2) is the parent element!!! WHY?");
              }
            }
            bool new_nds_master = true;
            bool new_nds_slave = true;
            for (std::vector<std::vector<int>>::iterator i = all_used_nds_master.begin();
                 i != all_used_nds_master.end(); ++i)
            {
              std::vector<int> used_nds_master = *i;
              if (used_nds_master == nds_master)
              {
                new_nds_master = false;
              }
            }
            for (std::vector<std::vector<int>>::iterator i = all_used_nds_slave.begin();
                 i != all_used_nds_slave.end(); ++i)
            {
              std::vector<int> used_nds_slave = *i;
              if (used_nds_slave == nds_slave)
              {
                new_nds_slave = false;
              }
            }
            if (new_nds_master == false and new_nds_slave == false)
            {
              continue;
            }
            if (new_nds_master == true)
            {
              all_used_nds_master.push_back(nds_master);
            }
            if (new_nds_slave == true)
            {
              all_used_nds_slave.push_back(nds_slave);
            }
            //------------------------
            num_edgestab++;
            // at least one element has to be cut
            if (p_master_handle->IsIntersected() or p_slave_handle->IsIntersected())
            {
              num_ghostpenalty++;

              face_type = INPAR::XFEM::face_type_ghost_penalty;
            }
            else
              face_type = INPAR::XFEM::face_type_std;

            XFEM::UTILS::GetVolumeCellMaterial(p_master, matptr_m, (*f)->Position());
            XFEM::UTILS::GetVolumeCellMaterial(p_slave, matptr_s, (*f)->Position());

            //--------------------------------------------------------------------------------------------
            // call evaluate and assemble routine
            AssembleEdgeStabGhostPenalty(eleparams, face_type, faceele, matptr_m, matptr_s,
                nds_master, nds_slave, *xdiscret, systemmatrix, systemvector);
            //--------------------------------------------------------------------------------------------
          }
          else if (vcs.size() == 1)
          {
            dserror("just one vcs reasonable?! face %d", faceele->Id());
          }
        }  // facet outside or (inside and include_inner)
        else if ((*f)->Position() == GEO::CUT::Point::undecided)
        {
          dserror("the position of this facet is undecided, how to stabilize???");
        }
        else if ((*f)->Position() == GEO::CUT::Point::oncutsurface)
        {
#ifdef DEBUG
          std::cout << "the position of this facet of face " << faceele->Id()
                    << " is oncutsurface, we do not stabilize it!!! " << std::endl;
#endif
          // if a facet lies oncutsurface, then there is only one neighbor, we do not stabilize this
          // facet REMARK: in case of one part of the facet is physical and the other part lies on
          // cutsurface,
          //         then the physical part is stabilized via another facet lying on the same fluid
          //         element's side
        }
        else
        {
          // facet is inside!
          face_type = INPAR::XFEM::face_type_ghost;
        }
      }  // loop facets
    }
    else
      dserror("not supported for this elements");
  }  // end second case: element handles for both parent elements
  //------------------------------------------------------------------------------
  // third case: element handle only for master element or for slave element available
  // at most one element cut
  //------------------------------------------------------------------------------
  else if ((p_master_handle != NULL and p_slave_handle == NULL) or
           (p_master_handle == NULL and p_slave_handle != NULL))
  {
    // linear elements
    if (p_master->Shape() == DRT::Element::hex8 or p_master->Shape() == DRT::Element::tet4 or
        p_master->Shape() == DRT::Element::wedge6 or p_master->Shape() == DRT::Element::pyramid5 or
        p_master->Shape() == DRT::Element::hex20 or p_master->Shape() == DRT::Element::hex27 or
        p_master->Shape() == DRT::Element::tet10 or p_master->Shape() == DRT::Element::wedge15)
    {
      GEO::CUT::SideHandle* side = GetFace(faceele, wizard);

      // facet of current side
      std::vector<GEO::CUT::Facet*> facets;
      side->Facets(facets);

      if (p_master->Shape() == DRT::Element::hex8 or p_master->Shape() == DRT::Element::tet4 or
          p_master->Shape() == DRT::Element::wedge6 or p_master->Shape() == DRT::Element::pyramid5)
      {
        if (facets.size() != 1) dserror("there has to be 1 facet equal to the side");
      }

      // get the unique single facet
      GEO::CUT::Facet* f = facets[0];
      if (f->Position() == GEO::CUT::Point::outside or
          (f->Position() == GEO::CUT::Point::inside and include_inner))
      {
        GEO::CUT::plain_volumecell_set vcs = f->Cells();

        if (vcs.size() != 1)
          dserror("there has to be 1 volumecell equal to the side");
        else
        {
          //------------------------ create nodal dof sets
          {
            TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");

            GEO::CUT::VolumeCell* vc = *(vcs.begin());

            // get the parent element
            int vc_ele_id = vc->ParentElement()->Id();
            if (vc_ele_id == -1)
            {
              vc_ele_id = vc->ParentElement()->GetParentId();
            }


            // which element is the parent element
            if (p_master_handle != NULL)
            {
              nds_master = vc->NodalDofSet();

              for (size_t i = 0; i < p_slave_numnode; i++) nds_slave.push_back(0);
            }
            else if (p_slave_handle != NULL)
            {
              for (size_t i = 0; i < p_master_numnode; i++) nds_master.push_back(0);

              nds_slave = vc->NodalDofSet();
            }
            else
              dserror("no element (ele1 and ele2) is the parent element!!! WHY?");
          }
          //------------------------

          num_edgestab++;

          // at most one element can be a cut one
          if (p_master_handle != NULL)
          {
            if (p_master_handle->IsIntersected())
            {
              num_ghostpenalty++;
              face_type = INPAR::XFEM::face_type_ghost_penalty;
            }
            else
              face_type = INPAR::XFEM::face_type_std;
          }
          else if (p_slave_handle != NULL)
          {
            if (p_slave_handle->IsIntersected())
            {
              num_ghostpenalty++;
              face_type = INPAR::XFEM::face_type_ghost_penalty;
            }
            else
              face_type = INPAR::XFEM::face_type_std;
          }
          else
            face_type = INPAR::XFEM::face_type_std;


          // Get materials:
          XFEM::UTILS::GetVolumeCellMaterial(p_master, matptr_m, f->Position());
          XFEM::UTILS::GetVolumeCellMaterial(p_slave, matptr_s, f->Position());

          //--------------------------------------------------------------------------------------------

          // call evaluate and assemble routine
          AssembleEdgeStabGhostPenalty(eleparams, face_type, faceele, matptr_m, matptr_s,
              nds_master, nds_slave, *xdiscret, systemmatrix, systemvector);

          //--------------------------------------------------------------------------------------------
        }

      }  // if outside or (inside and include_inner)
    }
  }  // end last case



  if (gmsh_eos_out)
  {
    ghost_penalty_stab_.insert(std::pair<int, int>(faceele->Id(), num_ghostpenalty));
    edge_based_stab_.insert(std::pair<int, int>(faceele->Id(), num_edgestab));
  }

  //--------------------------------------------------------------------------------------------

  return;
}


/*----------------------------------------------------------------------*
 | calls the evaluate and assemble routine for edge based stabilization |
 | and ghost penaly in the XFEM                            schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::AssembleEdgeStabGhostPenalty(
    Teuchos::ParameterList& eleparams,        ///< element parameter list
    const INPAR::XFEM::FaceType& face_type,   ///< which type of face std, ghost, ghost-penalty
    DRT::ELEMENTS::FluidIntFace* intface,     ///< internal face element
    Teuchos::RCP<MAT::Material>& material_m,  ///< material of the master side
    Teuchos::RCP<MAT::Material>& material_s,  ///< material of the slave side
    std::vector<int>& nds_master,             ///< nodal dofset vector w.r.t. master element
    std::vector<int>& nds_slave,              ///< nodal dofset vector w.r.t. slave element
    DRT::DiscretizationFaces& xdiscret,       ///< XFEM discretization
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
)
{
  // If Saftey check is passed, both elements contain the same material and with the same settings
  XFEM::UTILS::SafetyCheckMaterials(material_m, material_s);

  //======================================================================================
  // call the internal faces stabilization routine for the current side/surface
  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: AssembleEdgeStabGhostPenalty");

  // set action and facetype for elements

  // TODO: set here the right stab-type LPS or EOS
  eleparams.set<int>("action", FLD::EOS_and_GhostPenalty_stabilization);


  // call the egde-based assemble and evaluate routine
  DRT::ELEMENTS::FluidIntFaceImplInterface::Impl(intface)->AssembleInternalFacesUsingNeighborData(
      intface, material_m, nds_master, nds_slave, face_type, eleparams, xdiscret, systemmatrix,
      systemvector);


  return;
}


/*----------------------------------------------------------------------*
 | get the cut side for face's element identified using the sorted      |
 | node ids                                                schott 04/12 |
 *----------------------------------------------------------------------*/
GEO::CUT::SideHandle* XFEM::XFEM_EdgeStab::GetFace(
    DRT::Element* faceele, Teuchos::RCP<GEO::CutWizard> wizard)
{
  TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: GetFace");

  // get faceele's nodes
  const int numnode = faceele->NumNode();
  std::vector<int> nodeids(numnode);

  for (int inode = 0; inode < numnode; inode++)
  {
    nodeids[inode] = faceele->NodeIds()[inode];
  }

  std::sort(nodeids.begin(), nodeids.end());

  return wizard->GetSide(nodeids);
}

/*----------------------------------------------------------------------*
 | reset maps for output                                    kruse 04/15 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::Reset()
{
  ghost_penalty_stab_.clear();
  edge_based_stab_.clear();
}



/*----------------------------------------------------------------------*
 |  prepares edge based stabilization for standard fluid   schott 05/12 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::EvaluateEdgeStabStd(
    Teuchos::ParameterList& eleparams,                ///< element parameter list
    Teuchos::RCP<DRT::Discretization> discret,        ///< discretization
    DRT::ELEMENTS::FluidIntFace* faceele,             ///< face element
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
)
{
  Teuchos::RCP<DRT::DiscretizationFaces> xdiscret =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


  // get the parent fluid elements
  DRT::ELEMENTS::Fluid* p_master = faceele->ParentMasterElement();
  DRT::ELEMENTS::Fluid* p_slave = faceele->ParentSlaveElement();

  size_t p_master_numnode = p_master->NumNode();
  size_t p_slave_numnode = p_slave->NumNode();

  //------------------------------------------------------------------------------
  // simplest case: no element handles for both parent elements
  // two uncut elements / standard fluid case
  //------------------------------------------------------------------------------

  std::vector<int> nds_master(p_master_numnode, 0);
  std::vector<int> nds_slave(p_slave_numnode, 0);

  // Provide material at both sides:
  Teuchos::RCP<MAT::Material> matptr_m;
  Teuchos::RCP<MAT::Material> matptr_s;
  matptr_m = p_master->Material();
  matptr_s = p_slave->Material();

  //--------------------------------------------------------------------------------------------

  // call evaluate and assemble routine
  AssembleEdgeStabGhostPenalty(eleparams, INPAR::XFEM::face_type_std, faceele, matptr_m, matptr_s,
      nds_master, nds_slave, *xdiscret, systemmatrix, systemvector);

  //--------------------------------------------------------------------------------------------

  return;
}

/*----------------------------------------------------------------------*
 |  prepares edge based stabilization for fluid-fluid applications      |
 |  where EOS pressure stab. shall be applied to interface-contributing |
 |  embedded fluid elements                               (kruse 10/14) |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::EvaluateEdgeStabBoundaryGP(
    Teuchos::ParameterList& eleparams,          ///< element parameter list
    Teuchos::RCP<DRT::Discretization> discret,  ///< discretization
    Teuchos::RCP<DRT::Discretization>
        boundarydiscret,  ///< auxiliary discretization of interface-contributing elements
    DRT::ELEMENTS::FluidIntFace* faceele,             ///< face element
    Teuchos::RCP<LINALG::SparseMatrix> systemmatrix,  ///< systemmatrix
    Teuchos::RCP<Epetra_Vector> systemvector          ///< systemvector
)
{
  Teuchos::RCP<DRT::DiscretizationFaces> xdiscret =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


  // get the parent fluid elements
  DRT::ELEMENTS::Fluid* p_master = faceele->ParentMasterElement();
  DRT::ELEMENTS::Fluid* p_slave = faceele->ParentSlaveElement();

  size_t p_master_numnode = p_master->NumNode();
  size_t p_slave_numnode = p_slave->NumNode();


  std::vector<int> nds_master;
  nds_master.reserve(p_master_numnode);

  std::vector<int> nds_slave;
  nds_slave.reserve(p_slave_numnode);

  //------------------------------------------------------------------------------
  // simplest case: no element handles for both parent elements
  // two uncut elements / standard fluid case
  //------------------------------------------------------------------------------

  {
    TEUCHOS_FUNC_TIME_MONITOR("XFEM::Edgestab EOS: create nds");

    for (size_t i = 0; i < p_master_numnode; i++) nds_master.push_back(0);

    for (size_t i = 0; i < p_slave_numnode; i++) nds_slave.push_back(0);
  }

  //--------------------------------------------------------------------------------------------
  // leave, if neither slave nor master element of this face contributes to the fluid-fluid
  // interface
  if (!(boundarydiscret->HaveGlobalElement(p_master->Id()) ||
          boundarydiscret->HaveGlobalElement(p_slave->Id())))
    return;

  // Provide material at both sides:
  Teuchos::RCP<MAT::Material> matptr_m;
  Teuchos::RCP<MAT::Material> matptr_s;
  matptr_m = p_master->Material();
  matptr_s = p_slave->Material();


  // call evaluate and assemble routine
  AssembleEdgeStabGhostPenalty(eleparams, INPAR::XFEM::face_type_boundary_ghost_penalty, faceele,
      matptr_m, matptr_s, nds_master, nds_slave, *xdiscret, systemmatrix, systemvector);

  //--------------------------------------------------------------------------------------------

  return;
}
