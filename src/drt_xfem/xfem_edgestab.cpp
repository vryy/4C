/*!
\file xfem_edgestab.cpp

\brief provides the xfem fluid and ghost penalty stabilization based on EOS/CIP (edge-oriented, continuous interior penalty) scheme

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/




#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_assemblestrategy.H"
#include "../drt_lib/drt_discret_faces.H"

#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_volumecell.H"

#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_intfaces_calc.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_fluid_xfluid/xfluid_defines.H"

#include "../drt_inpar/inpar_xfem.H"
#include "../drt_inpar/inpar_fluid.H"

#include "xfem_fluidwizard.H"

#include "xfem_edgestab.H"


/*----------------------------------------------------------------------*
 |  Constructor for XFEM_EdgeStab                          schott 03/12 |
 *----------------------------------------------------------------------*/
XFEM::XFEM_EdgeStab::XFEM_EdgeStab(
  Teuchos::RCP<XFEM::FluidWizard>              wizard,  ///< fluid wizard
  Teuchos::RCP<DRT::Discretization>            discret  ///< discretization
  ) :
  wizard_(wizard),
  discret_(discret)
{
  ghost_penalty_stab_.clear();
  edge_based_stab_.clear();
} // end constructor


/*----------------------------------------------------------------------*
 |  prepares edge based stabilization and ghost penaly in case of XFEM  |
 |  and calls evaluate routine                             schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::EvaluateEdgeStabGhostPenalty(
    Teuchos::ParameterList &               eleparams,        ///< element parameter list
    Teuchos::RCP<DRT::Discretization>      discret,          ///< discretization
    DRT::ELEMENTS::FluidIntFace *          faceele,          ///< face element
    Teuchos::RCP<LINALG::SparseMatrix>     systemmatrix,     ///< systemmatrix
    Teuchos::RCP<Epetra_Vector>            systemvector,     ///< systemvector
    bool                                   gmsh_eos_out      ///< stabilization gmsh output
)
{


  //====================================================================================================
  // implementation of edge-based fluid stabilization and ghost penalty
  //====================================================================================================

  // EDGE-based fluid stabilization and EDGE-based Ghost-penalty stabilization
  // REMARK: the current implementation of edge-based stabilization is based on the DiscretizationXFEM extension
  //         using additional information about faces between two volume elements
  // * fluid stabilization has to be integrated for all internal faces
  // * ghost penalty has to be integrated if there is at least one cut element
  //   (because all faces between two elements for that at least one element is cut by the interface has to be stabilized)
  //
  // we distinguish different stabilization cases
  //  1. the master element and slave element (connected via current side)
  //     do not have an elementhandle (standard fluid case)
  //     -> standard fluid stabilization
  //                               => EOS(fluid): YES         GHOST-PENALTY: NO
  //  2. element handles for both parent elements
  //     -> stabilization for each facet and corresponding volumecells of parent elements
  //                               => EOS(fluid): YES         GHOST-PENALTY: Yes (if at least one parent element is cut)
  //                                                                         NO  (if both parent elements are uncut)
  //  3. just one elementhandle available (at limit of bounding box)
  //     -> stabilization for each facet and corresponding volumecells of parent elements
  //                               => EOS(fluid): YES         GHOST-PENALTY: Yes (if at least one parent element is cut)
  //                                                                         NO  (if both parent elements are uncut)


  RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


  // get the parent fluid elements
  DRT::ELEMENTS::Fluid* p_master = faceele->ParentMasterElement();
  DRT::ELEMENTS::Fluid* p_slave  = faceele->ParentSlaveElement();

  // get corresponding element handles if available
  GEO::CUT::ElementHandle * p_master_handle   = wizard_->GetElement( p_master );
  GEO::CUT::ElementHandle * p_slave_handle    = wizard_->GetElement( p_slave  );

  size_t p_master_numnode = p_master->NumNode();
  size_t p_slave_numnode  = p_slave->NumNode();

  // get the parent element
  int p_master_id = p_master->Id();
  int p_slave_id  = p_slave->Id();

  std::vector<int> nds_master;
  nds_master.reserve(p_master_numnode);

  std::vector<int> nds_slave;
  nds_slave.reserve(p_slave_numnode);

  INPAR::XFEM::FaceType face_type;

  int num_edgestab = 0;      // how often to stabilize this face for edgebased stabilizations
  int num_ghostpenalty = 0;  // how often to stabilize this face for ghost penalty stabilizations

  //------------------------------------------------------------------------------
  // simplest case: no element handles for both parent elements
  // two uncut elements / standard fluid case
  //------------------------------------------------------------------------------
  if( p_master_handle == NULL and p_slave_handle == NULL)
  {
    num_edgestab++;

    face_type = INPAR::XFEM::face_type_std;

    {
      TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: create nds" );

      for(size_t i=0; i< p_master_numnode; i++)  nds_master.push_back(0);

      for(size_t i=0; i< p_slave_numnode; i++)   nds_slave.push_back(0);
    }

    //--------------------------------------------------------------------------------------------

    // call evaluate and assemble routine
    AssembleEdgeStabGhostPenalty( eleparams,
                                  face_type,
                                  faceele,
                                  nds_master,
                                  nds_slave,
                                  *xdiscret,
                                  systemmatrix,
                                  systemvector);

    //--------------------------------------------------------------------------------------------

  }
  //------------------------------------------------------------------------------
  // second case: element handles for both parent elements
  // two elements that are maybe cut
  //------------------------------------------------------------------------------
  else if( p_master_handle != NULL and p_slave_handle != NULL)
  {
    // linear elements
    if(    p_master->Shape() == DRT::Element::hex8
        or p_master->Shape() == DRT::Element::tet4
        or p_master->Shape() == DRT::Element::wedge6
        or p_master->Shape() == DRT::Element::pyramid5 )
    {

      GEO::CUT::Side* side = GetCutSide(faceele);

      //-------------------------------- loop facets of this side -----------------------------
      // facet of current side
      const std::vector<GEO::CUT::Facet*> facets = side->Facets();

      if(facets.size() == 0) dserror("there is no facet between two elements with elementhandle!");

      // each facet should have 2 volumecells
      for(std::vector<GEO::CUT::Facet*>::const_iterator f=facets.begin(); f!=facets.end(); f++)
      {
        if((*f)->Position() == GEO::CUT::Point::outside /*or (*f)->Position() == GEO::CUT::Point::oncutsurface*/)
        {

          GEO::CUT::plain_volumecell_set vcs = (*f)->Cells();

          // how many volumecells found?
          if(vcs.size() == 2) // standard XFEM case (facet between two vcs of two neighbouring cut elements
          {
            //------------------------ create nodal dof sets
            {
              TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: create nds" );

              GEO::CUT::VolumeCell* vc1 = vcs[0];
              GEO::CUT::VolumeCell* vc2 = vcs[1];

#ifdef DOFSETS_NEW

              // get the parent element
              int vc_ele1_id = vc1->ParentElement()->Id();
              int vc_ele2_id = vc2->ParentElement()->Id();

              // which element is the parent element
              if(vc_ele1_id == p_master_id)
              {
                nds_master = vc1->NodalDofSet();
                nds_slave  = vc2->NodalDofSet();
              }
              else if(vc_ele2_id == p_master_id)
              { // switch ele 1 <-> ele 2
                nds_master = vc2->NodalDofSet();
                nds_slave  = vc1->NodalDofSet();
              }
              else dserror("no element (ele1 and ele2) is the parent element!!! WHY?");
#else
              for(size_t i=0; i< p_master_numnode; i++)  nds_master.push_back(0);

              for(size_t i=0; i< p_slave_numnode; i++)   nds_slave.push_back(0);
#endif

            }
            //------------------------

            num_edgestab++;

            // at least one element has to be cut
            if(p_master_handle->IsCut() or p_slave_handle->IsCut())
            {
              num_ghostpenalty++;

              face_type = INPAR::XFEM::face_type_ghost_penalty;
            }
            else face_type = INPAR::XFEM::face_type_std;


            //--------------------------------------------------------------------------------------------

            // call evaluate and assemble routine
            AssembleEdgeStabGhostPenalty( eleparams,
                                          face_type,
                                          faceele,
                                          nds_master,
                                          nds_slave,
                                          *xdiscret,
                                          systemmatrix,
                                          systemvector);

            //--------------------------------------------------------------------------------------------

          }
          else if(vcs.size() == 1)
          {
            dserror("just one vcs reasonable?! face %d", faceele->Id());
          }
        } // facet outside
        else if((*f)->Position() == GEO::CUT::Point::undecided)
        {
          dserror("the position of this facet is undecided, how to stabilize???");
        }
        else if((*f)->Position() == GEO::CUT::Point::oncutsurface)
        {
#ifdef DEBUG
          cout << "the position of this facet of face " << faceele->Id() << " is oncutsurface, we do not stabilize it!!! " << endl;
#endif
          // if a facet lies oncutsurface, then there is only one neighbor, we do not stabilize this facet
          // REMARK: in case of one part of the facet is physical and the other part lies on cutsurface,
          //         then the physical part is stabilized via another facet lying on the same fluid element's side
        }
        else
        {
          // facet is inside!
          face_type = INPAR::XFEM::face_type_ghost;
        }

      } // loop facets
    } // hex 8 elements
    else dserror("not supported for not hex8 elements");

  } // end second case: element handles for both parent elements
  //------------------------------------------------------------------------------
  // third case: element handle only for master element or for slave element available
  // at most one element cut
  //------------------------------------------------------------------------------
  else if (   (p_master_handle != NULL and p_slave_handle == NULL )
           or (p_master_handle == NULL and p_slave_handle != NULL ) )
  {
    // linear elements
    if(    p_master->Shape() == DRT::Element::hex8
        or p_master->Shape() == DRT::Element::tet4
        or p_master->Shape() == DRT::Element::wedge6
        or p_master->Shape() == DRT::Element::pyramid5 )
    {

      GEO::CUT::Side* side = GetCutSide(faceele);

      // facet of current side
      const std::vector<GEO::CUT::Facet*> facets = side->Facets();

      if(facets.size() != 1) dserror("there has to be 1 facet equal to the side");

      // get the unique single facet
      GEO::CUT::Facet* f = facets[0];

      if(f->Position() == GEO::CUT::Point::outside /*or (*f)->Position() == GEO::CUT::Point::oncutsurface*/)
      {

          GEO::CUT::plain_volumecell_set vcs = f->Cells();

          if(vcs.size() != 1) dserror("there has to be 1 volumecell equal to the side");
          else
          {
            //------------------------ create nodal dof sets
            {
              TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: create nds" );

              GEO::CUT::VolumeCell* vc = vcs[0];

              // get the parent element
              int vc_ele_id = vc->ParentElement()->Id();

#ifdef DOFSETS_NEW
              // which element is the parent element
              if(vc_ele_id == p_master_id)
              {
                nds_master = vc->NodalDofSet();

                for(size_t i=0; i< p_slave_numnode; i++)  nds_slave.push_back(0);
              }
              else if(vc_ele_id == p_slave_id)
              {
                for(size_t i=0; i< p_master_numnode; i++)  nds_master.push_back(0);

                nds_slave  = vc->NodalDofSet();
              }
              else dserror("no element (ele1 and ele2) is the parent element!!! WHY?");
#else
              for(size_t i=0; i< p_master_numnode; i++)  nds_master.push_back(0);

              for(size_t i=0; i< p_slave_numnode; i++)   nds_slave.push_back(0);
#endif

            }
            //------------------------

            num_edgestab++;

            // at most one element can be a cut one
            if(p_master_handle != NULL)
            {
              if(p_master_handle->IsCut())
              {
                num_ghostpenalty++;
                face_type = INPAR::XFEM::face_type_ghost_penalty;
              }
              else face_type = INPAR::XFEM::face_type_std;
            }
            else if(p_slave_handle != NULL)
            {
              if(p_slave_handle->IsCut())
              {
                num_ghostpenalty++;
                face_type = INPAR::XFEM::face_type_ghost_penalty;
              }
              else face_type = INPAR::XFEM::face_type_std;
            }
            else face_type = INPAR::XFEM::face_type_std;

            //--------------------------------------------------------------------------------------------

            // call evaluate and assemble routine
            AssembleEdgeStabGhostPenalty( eleparams,
                                          face_type,
                                          faceele,
                                          nds_master,
                                          nds_slave,
                                          *xdiscret,
                                          systemmatrix,
                                          systemvector);

            //--------------------------------------------------------------------------------------------

          }

      } // if outside
    } // if linear elements

  } // end last case



  if(gmsh_eos_out)
  {
    ghost_penalty_stab_.insert(std::pair<int,int>(faceele->Id(),num_ghostpenalty));
    edge_based_stab_.insert(std::pair<int,int>(faceele->Id(),num_edgestab));
  }

  //--------------------------------------------------------------------------------------------

  return;
}


/*----------------------------------------------------------------------*
 | calls the evaluate and assemble routine for edge based stabilization |
 | and ghost penaly in the XFEM                            schott 03/12 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::AssembleEdgeStabGhostPenalty( Teuchos::ParameterList &               eleparams,        ///< element parameter list
                                                        const INPAR::XFEM::FaceType &          face_type,        ///< which type of face std, ghost, ghost-penalty
                                                        DRT::ELEMENTS::FluidIntFace*           intface,          ///< internal face element
                                                        std::vector<int> &                     nds_master,       ///< nodal dofset vector w.r.t. master element
                                                        std::vector<int> &                     nds_slave,        ///< nodal dofset vector w.r.t. slave element
                                                        DRT::DiscretizationFaces &              xdiscret,         ///< XFEM discretization
                                                        Teuchos::RCP<LINALG::SparseMatrix>     systemmatrix,     ///< systemmatrix
                                                        Teuchos::RCP<Epetra_Vector>            systemvector      ///< systemvector
                                                        )
{


  //======================================================================================
  // call the internal faces stabilization routine for the current side/surface
  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: AssembleEdgeStabGhostPenalty" );

  // set action and facetype for elements
  eleparams.set<int>("action", FLD::EOS_and_GhostPenalty_stabilization);
  eleparams.set("facetype", face_type);

  // call the egde-based assemble and evaluate routine
  DRT::ELEMENTS::FluidIntFaceImplInterface::Impl(intface)->AssembleInternalFacesUsingNeighborData(     intface,
                                                                                                       nds_master,
                                                                                                       nds_slave,
                                                                                                       eleparams,
                                                                                                       xdiscret,
                                                                                                       systemmatrix,
                                                                                                       systemvector);


  return;
}


/*----------------------------------------------------------------------*
 | get the cut side for face's element identified using the sorted      |
 | node ids                                                schott 04/12 |
 *----------------------------------------------------------------------*/
GEO::CUT::Side* XFEM::XFEM_EdgeStab::GetCutSide(DRT::Element* faceele)
{

  TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: FindSide" );

  // get faceele's nodes
  const int numnode = faceele->NumNode();
  std::vector<int> nodeids(numnode);

  for(int inode=0;inode<numnode; inode++)
  {
    nodeids[inode] = faceele->NodeIds()[inode];
  }

  return wizard_->GetSide(nodeids);
}



/*----------------------------------------------------------------------*
 |  prepares edge based stabilization for standard fluid   schott 05/12 |
 *----------------------------------------------------------------------*/
void XFEM::XFEM_EdgeStab::EvaluateEdgeStabStd(
    Teuchos::ParameterList &               eleparams,        ///< element parameter list
    Teuchos::RCP<DRT::Discretization>      discret,          ///< discretization
    DRT::ELEMENTS::FluidIntFace *          faceele,          ///< face element
    Teuchos::RCP<LINALG::SparseMatrix>     systemmatrix,     ///< systemmatrix
    Teuchos::RCP<Epetra_Vector>            systemvector      ///< systemvector
)
{
  RCP<DRT::DiscretizationFaces> xdiscret = Teuchos::rcp_dynamic_cast<DRT::DiscretizationFaces>(discret);
  if (xdiscret == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationFaces.");


  // get the parent fluid elements
  DRT::ELEMENTS::Fluid* p_master = faceele->ParentMasterElement();
  DRT::ELEMENTS::Fluid* p_slave  = faceele->ParentSlaveElement();

  size_t p_master_numnode = p_master->NumNode();
  size_t p_slave_numnode  = p_slave->NumNode();


  std::vector<int> nds_master;
  nds_master.reserve(p_master_numnode);

  std::vector<int> nds_slave;
  nds_slave.reserve(p_slave_numnode);

  //------------------------------------------------------------------------------
  // simplest case: no element handles for both parent elements
  // two uncut elements / standard fluid case
  //------------------------------------------------------------------------------

  {
    TEUCHOS_FUNC_TIME_MONITOR( "XFEM::Edgestab EOS: create nds" );

    for(size_t i=0; i< p_master_numnode; i++)  nds_master.push_back(0);

    for(size_t i=0; i< p_slave_numnode; i++)   nds_slave.push_back(0);
  }

  //--------------------------------------------------------------------------------------------

  // call evaluate and assemble routine
   AssembleEdgeStabGhostPenalty( eleparams,
                                 INPAR::XFEM::face_type_std,
                                 faceele,
                                 nds_master,
                                 nds_slave,
                                 *xdiscret,
                                 systemmatrix,
                                 systemvector);

  //--------------------------------------------------------------------------------------------

  return;
}
