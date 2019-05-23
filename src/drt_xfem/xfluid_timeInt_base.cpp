/*----------------------------------------------------------------------*/
/*!

\brief provides the basic class for XFEM-time-integration, e.g. for Semi-Lagrangean methods

\level 2

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

*/
/*----------------------------------------------------------------------*/


#include "../drt_bele3/bele3.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_utils.H"

#include "../linalg/linalg_utils.H"

#include "../drt_geometry/position_array.H"

#include "../drt_inpar/inpar_xfem.H"

#include "../drt_cut/cut_volumecell.H"
#include "../drt_cut/cut_sidehandle.H"
#include "../drt_cut/cut_elementhandle.H"
#include "../drt_cut/cut_intersection.H"
#include "../drt_cut/cut_position.H"
#include "../drt_cut/cut_cutwizard.H"


#include "../drt_io/io.H"
#include "../drt_io/io_gmsh.H"
#include "../drt_io/io_control.H"

#include "../drt_xfem/xfem_dofset.H"

#include "xfluid_timeInt_base.H"


//#define DEBUG_TIMINT_STD

/*------------------------------------------------------------------------------------------------*
 * basic XFEM time-integration constructor                                           schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_TIMEINT_BASE::XFLUID_TIMEINT_BASE(
    const Teuchos::RCP<DRT::Discretization> discret,      /// background discretization
    const Teuchos::RCP<DRT::Discretization> boundarydis,  /// cut discretization
    Teuchos::RCP<GEO::CutWizard> wizard_old,    /// cut wizard w.r.t. old interface position
    Teuchos::RCP<GEO::CutWizard> wizard_new,    /// cut wizard w.r.t. new interface position
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_old,  /// XFEM dofset w.r.t. old interface position
    Teuchos::RCP<XFEM::XFEMDofSet> dofset_new,  /// XFEM dofset w.r.t. new interface position
    std::vector<Teuchos::RCP<Epetra_Vector>>
        oldVectors,                      /// vector of col-vectors w.r.t. old interface position
    Teuchos::RCP<Epetra_Vector> dispn,   /// old col displacement vector
    Teuchos::RCP<Epetra_Vector> dispnp,  /// col displacment n +1
    const Epetra_Map& olddofcolmap,      /// dofcolmap w.r.t. old interface position
    const Epetra_Map& newdofrowmap,      /// dofcolmap w.r.t. new interface position
    const Teuchos::RCP<std::map<int, std::vector<int>>>
        pbcmap  /// map of periodic boundary conditions
    )
    : discret_(discret),
      boundarydis_(boundarydis),
      wizard_old_(wizard_old),
      wizard_new_(wizard_new),
      dofset_old_(dofset_old),
      dofset_new_(dofset_new),
      olddofcolmap_(olddofcolmap),
      newdofrowmap_(newdofrowmap),
      oldVectors_(oldVectors),
      dispn_(dispn),
      dispnp_(dispnp),
      timeIntData_(Teuchos::null),
      pbcmap_(pbcmap),
      myrank_(discret_->Comm().MyPID()),
      numproc_(discret_->Comm().NumProc()),
      newton_max_iter_(10),  /// maximal number of newton iterations for Semi-Lagrangean algorithm
      limits_tol_(1.0e-10),  /// newton tolerance for Semi-Lagrangean algorithm
      TOL_dist_(1.0e-12)  /// tolerance to find the shortest distance of point to its projection on
                          /// the surface dis
{
  return;
}  // end constructor


/*------------------------------------------------------------------------------------------------*
 * set the computation type with help of the iteration counter                       schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::type(int iter, int iterMax)
{
  // TODO: update this for FSI
  if (iter == 1)
    FGIType_ = FRS1FGI1_;
  else if (iterMax == 1 or iter % iterMax == 1)
    FGIType_ = FRS1FGINot1_;
  else
    FGIType_ = FRSNot1_;

  return;
}  // end function type



/*------------------------------------------------------------------------------------------------*
 * assign the Epetra vectors which shall be computed to the                                       *
 * algorithms data structure                                                         schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::handleVectors(
    std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn)
{
  newVectors_ = newRowVectorsn;

  if (oldVectors_.size() != newVectors_.size())
    dserror(
        "Number of state-vectors at new and old discretization are different!");  // but they have
                                                                                  // different maps
                                                                                  // (row vs col)

}  // end function handleVectors



/*------------------------------------------------------------------------------------------------*
 * check if the current point x2 changed the side compared to x1                     schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_TIMEINT_BASE::changedSideSameTime(
    const bool newTimeStep,    /// new/old timestep for both points x1 and x2
    DRT::Element* ele1,        /// first element where x1 lies in
    LINALG::Matrix<3, 1>& x1,  /// global coordinates of point x1
    DRT::Element* ele2,        /// second element where x2 lies in
    LINALG::Matrix<3, 1>& x2   /// global coordinates of point x2
    ) const
{
#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t\t\t check if point changed the side ";
#endif

  //-----------------------------------------------------------------------
  // special case of equal coordinates x1 and x2 -> no line
  LINALG::Matrix<3, 1> diff(true);
  diff.Update(1.0, x1, -1.0, x2);

  if (diff.Norm2() < 1.0e-13) return false;

  //-----------------------------------------------------------------------
  // standard case of a real line between x1 and x2

  Teuchos::RCP<GEO::CutWizard> wizard = newTimeStep ? wizard_new_ : wizard_old_;

  // REMARK:
  // changing the side of a point at two times (newton steps) with coordinates x1 and x2 is done
  // via changing the cut of the linear trace between x1 and x2 with surrounding cutting sides
  // if there is a cut between the ray tracer and a cutting side, then the point changed the side
  // while moving from Position x1 to Position x2

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // first: find involved elements which are possibly passed through by the ray between x1 and x2

#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t\t\t\t\t --------------------------------- ";
  IO::cout << "\n\t\t\t\t\t Find involved background elements ";
  IO::cout << "\n\t\t\t\t\t --------------------------------- ";
#endif

  int mi = 0;  // we assume only cutting sides with mesh 0 (meshintersection)

  bool check_allsides =
      false;  // flag to check all sides in case of one common node is a non-row node

  std::set<int> cut_sides;

  std::set<int> eids;
  std::set<int> common_nodes;

  if (ele1 == NULL or ele2 == NULL)
  {
    check_allsides = true;  // no efficient search possible
  }
  else
  {
    if (ele1->Id() == ele2->Id())  // line within one element
    {
      eids.insert(ele1->Id());  // just one element to check

#ifdef DEBUG_TIMINT_STD
      IO::cout << "\n\t\t\t\t\t check changing side just for the unique element" << ele1->Id()
               << IO::endl;
#endif
    }
    else if (Neighbors(ele1, ele2, common_nodes))
    {
#ifdef DEBUG_TIMINT_STD
      IO::cout << "\n\t\t\t\t\t check changing side for neighboring elements: " << IO::endl;
#endif

      // get all elements adjacent to common nodes in ele1 and ele2

      // REMARK: if at least one common node is not a row node, then there is possibly one adjacent
      // element missing, but then there are also involved sides missing -> check all sides in
      // cutdiscret, also if it is very slow!

      for (std::set<int>::iterator it = common_nodes.begin(); it != common_nodes.end(); it++)
      {
        DRT::Node* n = discret_->gNode(*it);

        if (n->Owner() != myrank_)
        {
          check_allsides =
              true;  // flag to check all sides in case of one common node is a non-row node
          break;
        }

        const int numele = n->NumElement();

        DRT::Element** elements = n->Elements();

        for (int e_it = 0; e_it < numele; e_it++)
        {
          DRT::Element* ele = elements[e_it];

          eids.insert(ele->Id());

#ifdef DEBUG_TIMINT_STD
          IO::cout << "\n\t\t\t\t\t add element " << ele->Id() << IO::endl;
#endif
        }
      }
    }
    else
    {
      // check all sides ghosted in boundarydis
      check_allsides = true;

#ifdef DEBUG_TIMINT_STD
      IO::cout << "\n\t\t\t\t\t all sides on boundary discretization have to be check, ele1 = "
               << ele1->Id() << " and ele2 = " << ele2->Id() << " are not Neighbors" << IO::endl;
#endif
    }
  }
  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // second: find all boundary sides involved in cutting the determined background elements


#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t\t\t\t\t --------------------------------- ";
  IO::cout << "\n\t\t\t\t\t Find involved sides with possible cuts between trace and side ";
  IO::cout << "\n\t\t\t\t\t --------------------------------- " << IO::endl;
#endif


  // collect all cutting sides
  if (check_allsides)
  {
    // add all sides of cut_discret to check
    for (int i = 0; i < boundarydis_->NumMyColElements(); i++)
    {
      cut_sides.insert(boundarydis_->ElementColMap()->GID(i));
    }
  }
  else
  {
    // loop all found element
    for (std::set<int>::iterator e_it = eids.begin(); e_it != eids.end(); e_it++)
    {
      DRT::Element* ele = discret_->gElement(*e_it);

      GEO::CUT::ElementHandle* eh = wizard->GetElement(ele);

      // no cutsides within this element, then no side changing possible
      if (eh == NULL)
      {
        continue;  // next element
      }

      //--------------------------------------------------------
      // get involved side ids
      GEO::CUT::plain_element_set elements;
      eh->CollectElements(elements);

      // get all side-ids
      for (GEO::CUT::plain_element_set::iterator eles = elements.begin(); eles != elements.end();
           eles++)
      {
        GEO::CUT::Element* sub_ele = *eles;

        GEO::CUT::plain_facet_set facets = sub_ele->Facets();

        for (GEO::CUT::plain_facet_set::const_iterator facet_it = facets.begin();
             facet_it != facets.end(); facet_it++)
        {
          GEO::CUT::Facet* facet = *facet_it;

          GEO::CUT::Side* parent_side = facet->ParentSide();

          bool is_elements_side = sub_ele->OwnedSide(parent_side);

          if (!is_elements_side)  // is a cutting side
          {
            cut_sides.insert(parent_side->Id());

#ifdef DEBUG_TIMINT_STD
            IO::cout << "\n\t\t\t\t\t add side with Id=" << parent_side->Id() << IO::endl;
#endif
          }  // !element's side
        }    // facets
      }      // sub elements
    }
  }


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // third: check cut between the line and all cutting sides
  for (std::set<int>::iterator side_it = cut_sides.begin(); side_it != cut_sides.end(); side_it++)
  {
    // get the side via sidehandle
    GEO::CUT::SideHandle* sh = wizard->GetMeshCuttingSide(*side_it, mi);

#ifdef DEBUG_TIMINT_STD
    IO::cout << "\n\t\t\t\t\t SIDE-CHECK with side=" << *side_it << IO::endl;
#endif

    if (callSideEdgeIntersection(sh, *side_it, x1, x2))
    {
#ifdef DEBUG_TIMINT_STD
      IO::cout << "\n\t\t\t\t\t <<< POINT CHANGED THE SIDE >>>" << IO::endl;
#endif

      return true;
    }
  }

  return false;
}


/*------------------------------------------------------------------------------------------------*
 * check if both element are neighbors sharing at least one common node              schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_TIMEINT_BASE::Neighbors(
    DRT::Element* ele1, DRT::Element* ele2, std::set<int>& common_nodes) const
{
  bool is_neighbor = false;

  const int numnode1 = ele1->NumNode();
  const int numnode2 = ele2->NumNode();

  const int* nodeids1 = ele1->NodeIds();
  const int* nodeids2 = ele2->NodeIds();

  std::vector<int> nodeids1_vec;
  std::vector<int> nodeids2_vec;

  for (int i = 0; i < numnode1; i++) nodeids1_vec.push_back(nodeids1[i]);

  for (int i = 0; i < numnode2; i++) nodeids2_vec.push_back(nodeids2[i]);

  sort(nodeids1_vec.begin(), nodeids1_vec.end());
  sort(nodeids2_vec.begin(), nodeids2_vec.end());

  //------------------------------------------------------
  // find common nodes
  int index1 = 0;
  int index2 = 0;
  while (true)
  {
    // if at least one vector of nodes is checked
    if (index1 >= numnode1 or index2 >= numnode2) break;

    if (nodeids1_vec[index1] < nodeids2_vec[index2])
      index1++;
    else if (nodeids1_vec[index1] > nodeids2_vec[index2])
      index2++;
    else  // nodeids1_vec[index1] == nodeids2_vec[index2]
    {
      common_nodes.insert(nodeids1_vec[index1]);

      is_neighbor = true;  // true if at least one node is common, find also further common nodes

      index1++;
      index2++;
    }
  }

  return is_neighbor;
}


/*------------------------------------------------------------------------------------------------*
 * check if edge between x1 and x2 cuts the side                                     schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_TIMEINT_BASE::callSideEdgeIntersection(GEO::CUT::SideHandle* sh,  /// side handle
    int sid,                                                                        /// side id
    LINALG::Matrix<3, 1>& x1,  /// coordinates of edge's start point
    LINALG::Matrix<3, 1>& x2   /// coordinates of edge's end point
    ) const
{
  switch (sh->Shape())
  {
    case DRT::Element::tri3:
    {
      return callSideEdgeIntersectionT<DRT::Element::tri3>(sh, sid, x1, x2);
    }
    case DRT::Element::quad4:
    {
      return callSideEdgeIntersectionT<DRT::Element::quad4>(sh, sid, x1, x2);
    }
    case DRT::Element::quad8:
    {
      return callSideEdgeIntersectionT<DRT::Element::quad8>(sh, sid, x1, x2);
    }
    case DRT::Element::quad9:
    {
      return callSideEdgeIntersectionT<DRT::Element::quad9>(sh, sid, x1, x2);
    }
    default:
    {
      dserror("unknown side shape");
      break;
    }
  }

  return false;
}

/*------------------------------------------------------------------------------------------------*
 * check if edge between x1 and x2 cuts the side (templated)                         schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType sidetype>
bool XFEM::XFLUID_TIMEINT_BASE::callSideEdgeIntersectionT(GEO::CUT::SideHandle* sh,  /// side handle
    int sid,                                                                         /// side id
    LINALG::Matrix<3, 1>& x1,  /// coordinates of edge's start point
    LINALG::Matrix<3, 1>& x2   /// coordinates of edge's end point
    ) const
{
  const int nsd = 3;
  const int numNodesSurface = DRT::UTILS::DisTypeToNumNodePerEle<sidetype>::numNodePerElement;

  LINALG::Matrix<nsd, 2> xyze_lineElement(true);

  for (int i = 0; i < nsd; i++)
  {
    xyze_lineElement(i, 0) = x1(i);
    xyze_lineElement(i, 1) = x2(i);
  }

  Epetra_SerialDenseMatrix xyze_side;
  sh->Coordinates(xyze_side);

  LINALG::Matrix<nsd, numNodesSurface> xyze_surfaceElement(xyze_side);

  LINALG::Matrix<3, 1> xsi(true);


  Teuchos::RCP<GEO::CUT::IntersectionBase> intersect =
      GEO::CUT::IntersectionBase::Create(DRT::Element::line2, sidetype);
  Teuchos::RCP<GEO::CUT::Options> options =
      Teuchos::rcp(new GEO::CUT::Options());  // Create cut options for intersection (specify to use
                                              // double prec.)
  intersect->Init(xyze_lineElement, xyze_surfaceElement, false, false, false, options.getRawPtr());

  // check also limits during the newton scheme and when converged
  double itol;
  return (static_cast<bool>(intersect->ComputeEdgeSideIntersection(itol)));
}


/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for an integration cell                 schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::callXToXiCoords(const DRT::Element* ele,  /// pointer to element
    LINALG::Matrix<3, 1>& x,   /// global coordinates of point
    LINALG::Matrix<3, 1>& xi,  /// determined local coordinates w.r.t ele
    const std::string state,   ///< state dispn or dispnp?
    bool& pointInDomain        /// lies point in element ?
    ) const
{
  LINALG::SerialDenseMatrix xyz(3, ele->NumNode(), true);
  GEO::fillInitialPositionArray(ele, xyz);

  // add ale displacements to initial position
  if (state != "reference" && dispnp_ != Teuchos::null &&
      dispn_ != Teuchos::null)  // add no displacements in case of state == "reference" or //is ale?
  {
    int nen = ele->NumNode();
    int numdof = ele->NumDofPerNode(*(ele->Nodes()[0]));

    std::vector<int> nds(nen, 0);

    DRT::Element::LocationArray la(1);
    ele->LocationVector(*discret_, nds, la, false);

    // extract local values of the global vectors
    std::vector<double> mydispnp(la[0].lm_.size());

    if (state == "dispnp")
      DRT::UTILS::ExtractMyValues(*dispnp_, mydispnp, la[0].lm_);
    else if (state == "dispn")
      DRT::UTILS::ExtractMyValues(*dispn_, mydispnp, la[0].lm_);
    else
      dserror("XFEM::XFLUID_TIMEINT_BASE::callXToXiCoords: Undefined state!");

    for (int inode = 0; inode < nen; ++inode)  // number of nodes
    {
      for (int idim = 0; idim < 3; ++idim)  // number of dimensions
      {
        (xyz)(idim, inode) +=
            mydispnp[idim + (inode * numdof)];  // attention! dispnp vector has 3+1 dofs for
                                                // displacement (the same as for (u,p))
      }
    }
  }

  callXToXiCoords(xyz, ele->Shape(), x, xi, pointInDomain);
}  // end function callXToXiCoords



/*------------------------------------------------------------------------------------------------*
 * call the computation of local coordinates for a polytop                                        *
 * with corners given by the coordinates                                             schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::callXToXiCoords(
    LINALG::SerialDenseMatrix& nodecoords,     /// node coordinates of element
    DRT::Element::DiscretizationType DISTYPE,  /// discretization type
    LINALG::Matrix<3, 1>& x,                   /// global coordinates of point
    LINALG::Matrix<3, 1>& xi,                  /// determined local coordinates w.r.t ele
    bool& pointInDomain                        /// lies point in element ?
    ) const
{
  switch (DISTYPE)
  {
    case DRT::Element::hex8:
      XToXiCoords<DRT::Element::hex8>(nodecoords, x, xi, pointInDomain);
      break;
    case DRT::Element::hex20:
      XToXiCoords<DRT::Element::hex20>(nodecoords, x, xi, pointInDomain);
      break;
    case DRT::Element::tet4:
      XToXiCoords<DRT::Element::tet4>(nodecoords, x, xi, pointInDomain);
      break;
    default:
      dserror("add your 3D distype and the according transformation!");
      break;
  }  // end switch

  return;
}  // end function callXToXiCoords


//! compute local element coordinates and check whether the according point is inside the element
template <DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_TIMEINT_BASE::XToXiCoords(
    LINALG::SerialDenseMatrix& xyz,  /// node coordinates of element
    LINALG::Matrix<3, 1>& x,         /// global coordinates of point
    LINALG::Matrix<3, 1>& xi,        /// determined local coordinates w.r.t ele
    bool& pointInCell                /// lies point in element?
    ) const
{
  const int nsd = 3;  // dimension
  const int numnode =
      DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;  // number of nodes of element

  LINALG::Matrix<nsd, numnode> xyze(xyz);

  Teuchos::RCP<GEO::CUT::Position> pos =
      GEO::CUT::PositionFactory::BuildPosition<3, DISTYPE>(xyze, x);
  pos->Compute();
  pos->LocalCoordinates(xi);  // local coordinates

  pointInCell = pos->WithinLimitsTol(limits_tol_);  // check if point is in element
};


//! data at an arbitrary point lying in an element
template <const int numnode, DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_TIMEINT_BASE::evalShapeAndDeriv(DRT::Element* element,  /// pointer to element
    LINALG::Matrix<3, 1>& xi,              /// local coordinates of point w.r.t element
    LINALG::Matrix<3, 3>& xji,             /// inverse of jacobian
    LINALG::Matrix<numnode, 1>& shapeFcn,  /// shape functions at point
    LINALG::Matrix<3, numnode>&
        shapeFcnDerivXY,  /// derivatives of shape function w.r.t global coordiantes xyz
    bool compute_deriv    /// shall derivatives and jacobian be computed
    ) const
{
  const int* elenodeids = element->NodeIds();  // nodeids of element
  const int nsd = 3;                           // dimension

  // clear data that shall be filled
  shapeFcn.Clear();
  shapeFcnDerivXY.Clear();

  //-------------------------------------------------------
  DRT::UTILS::shape_function_3D(
      shapeFcn, xi(0), xi(1), xi(2), DISTYPE);  // evaluate shape functions at xi

  if (compute_deriv)
  {
    LINALG::Matrix<nsd, numnode> nodecoords(true);       // node coordinates of the element
    for (size_t nodeid = 0; nodeid < numnode; nodeid++)  // fill node coordinates
    {
      DRT::Node* currnode = discret_->gNode(elenodeids[nodeid]);
      for (int i = 0; i < nsd; i++) nodecoords(i, nodeid) = currnode->X()[i];
    }

    // add ale displacements to initial position for state np
    if (dispnp_ != Teuchos::null)  // is ale?
    {
      int nen = element->NumNode();
      int numdof = element->NumDofPerNode(*(element->Nodes()[0]));

      std::vector<int> nds(nen, 0);
      DRT::Element::LocationArray la(1);
      element->LocationVector(*discret_, nds, la, false);

      // extract local values of the global vectors
      std::vector<double> mydispnp(la[0].lm_.size());
      DRT::UTILS::ExtractMyValues(*dispnp_, mydispnp, la[0].lm_);

      for (int inode = 0; inode < nen; ++inode)  // number of nodes
      {
        for (int idim = 0; idim < nsd; ++idim)  // number of dimensions
        {
          (nodecoords)(idim, inode) +=
              mydispnp[idim + (inode * numdof)];  // attention! dispnp vector has 3+1 dofs for
                                                  // displacement (the same as for (u,p))
        }
      }
    }

    // shape function derivatives w.r.t local coordinates
    LINALG::Matrix<3, numnode> shapeFcnDeriv;
    DRT::UTILS::shape_function_3D_deriv1(shapeFcnDeriv, xi(0), xi(1), xi(2), DISTYPE);

    LINALG::Matrix<nsd, nsd> xjm(true);         // jacobi matrix
    xjm.MultiplyNT(shapeFcnDeriv, nodecoords);  // jacobian J = (dx/dxi)^T
    xji.Clear();
    xji.Invert(xjm);  // jacobian inverted J^(-1) = dxi/dx

    shapeFcnDerivXY.Multiply(xji, shapeFcnDeriv);  // (dN/dx)^T = (dN/dxi)^T * J^(-T)
  }

  return;
};  // end template evalShapeAndDeriv



// explicit instantiations of this function...
template void XFEM::XFLUID_TIMEINT_BASE::evalShapeAndDeriv<8, DRT::Element::hex8>(DRT::Element*,
    LINALG::Matrix<3, 1>&, LINALG::Matrix<3, 3>&, LINALG::Matrix<8, 1>&, LINALG::Matrix<3, 8>&,
    bool) const;
template void XFEM::XFLUID_TIMEINT_BASE::evalShapeAndDeriv<20, DRT::Element::hex20>(DRT::Element*,
    LINALG::Matrix<3, 1>&, LINALG::Matrix<3, 3>&, LINALG::Matrix<20, 1>&, LINALG::Matrix<3, 20>&,
    bool) const;
template void XFEM::XFLUID_TIMEINT_BASE::evalShapeAndDeriv<4, DRT::Element::tet4>(DRT::Element*,
    LINALG::Matrix<3, 1>&, LINALG::Matrix<3, 3>&, LINALG::Matrix<4, 1>&, LINALG::Matrix<3, 4>&,
    bool) const;



/*------------------------------------------------------------------------------------------------*
 * add adjacent elements for a periodic boundary node                            winklmaier 05/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::addPBCelements(
    const DRT::Node* node, std::vector<const DRT::Element*>& eles) const
{
  dserror("what to do in addPBCelements?");

  const DRT::Element* const* elements = node->Elements();  // element around current node

  for (int iele = 0; iele < node->NumElement(); iele++) eles.push_back(elements[iele]);

  // get pbcnode
  bool pbcnodefound = false;  // boolean indicating whether this node is a pbc node
  DRT::Node* pbcnode = NULL;
  findPBCNode(node, pbcnode, pbcnodefound);

  // add elements located around the coupled pbc node
  if (pbcnodefound)
  {
    // get adjacent elements of this node
    const DRT::Element* const* pbcelements = pbcnode->Elements();
    // add elements to list
    for (int iele = 0; iele < pbcnode->NumElement(); iele++)  // = ptToNode->Elements();
    {
      eles.push_back(pbcelements[iele]);
    }
  }  // end if pbcnode true
}


/*------------------------------------------------------------------------------------------------*
 * find the PBC node                                                             winklmaier 05/11 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::findPBCNode(
    const DRT::Node* node, DRT::Node*& pbcnode, bool& pbcnodefound) const
{
  const int nodegid = node->Id();  // global node id
  pbcnodefound = false;            // boolean indicating whether this node is a pbc node
  int coupnodegid = -1;
  // loop all nodes with periodic boundary conditions (master nodes)
  for (std::map<int, std::vector<int>>::const_iterator pbciter = (*pbcmap_).begin();
       pbciter != (*pbcmap_).end(); ++pbciter)
  {
    if (pbciter->first == nodegid)  // node is a pbc master node
    {
      pbcnodefound = true;
      // coupled node is the (first) slave node
      coupnodegid = pbciter->second[0];
      if (pbciter->second.size() != 1)
        dserror("this might need to be modified for more than 1 slave per master");
    }
    else
    {
      // loop all slave nodes
      for (size_t islave = 0; islave < pbciter->second.size(); islave++)
      {
        if (pbciter->second[islave] == nodegid)  // node is a pbc slave node
        {
          pbcnodefound = true;
          coupnodegid = pbciter->first;  // coupled node is the master node
        }
      }  // end loop over slave nodes
    }    // end if
  }      // end loop over pbc map

  // get pbcnode
  if (pbcnodefound) pbcnode = discret_->gNode(coupnodegid);
}  // end if pbcnode true}



/*------------------------------------------------------------------------------------------------*
 * reset a special state with another state in the data class                        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::resetState(
    TimeIntData::state oldState, TimeIntData::state newState) const
{
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->state_ == oldState) data->state_ = newState;
  }

  return;
}  // end function resetState



/*------------------------------------------------------------------------------------------------*
 * clear the data of all nodes having a special state                                schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::clearState(TimeIntData::state state  /// state of time int to clear
    ) const
{
  std::vector<TimeIntData>::iterator data;
  while (true)  // while loop over data to be cleared
  {
    for (data = timeIntData_->begin(); data != timeIntData_->end();
         data++)  // for loop over all data
    {
      if (data->state_ == state)
      {
        timeIntData_->erase(data);
        break;
      }
    }
    if (data == timeIntData_->end()) break;
  }  // end while loop over data to be cleared

  return;
}  // end function clear state



/*------------------------------------------------------------------------------------------------*
 * basic function sending data to dest and receiving data from source                schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::sendData(
    DRT::PackBuffer& dataSend, int& dest, int& source, std::vector<char>& dataRecv) const
{
  std::vector<int> lengthSend(1, 0);
  lengthSend[0] = dataSend().size();
  int size_one = 1;

#ifdef DEBUG
  std::cout << "--- sending " << lengthSend[0] << " bytes: from proc " << myrank_ << " to proc "
            << dest << std::endl;
#endif

  // exporter for sending
  DRT::Exporter exporter(discret_->Comm());

  // send length of the data to be received ...
  MPI_Request req_length_data;
  int length_tag = 0;
  exporter.ISend(myrank_, dest, &(lengthSend[0]), size_one, length_tag, req_length_data);
  // ... and receive length
  std::vector<int> lengthRecv(1, 0);
  exporter.Receive(source, length_tag, lengthRecv, size_one);
  exporter.Wait(req_length_data);

  // send actual data ...
  int data_tag = 4;
  MPI_Request req_data;
  exporter.ISend(myrank_, dest, &(dataSend()[0]), lengthSend[0], data_tag, req_data);

  // ... and receive data
  dataRecv.clear();
  dataRecv.resize(lengthRecv[0]);
  exporter.ReceiveAny(source, data_tag, dataRecv, lengthRecv[0]);
  exporter.Wait(req_data);

#ifdef DEBUG
  std::cout << "--- receiving " << lengthRecv[0] << " bytes: to proc " << myrank_ << " from proc "
            << source << std::endl;
#endif
}  // end sendData



/*------------------------------------------------------------------------------------------------*
 * packing a node for parallel communication only with the basic nodal data                       *
 * without an underlying discretization fitting to the node's new prozessor          schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::packNode(DRT::PackBuffer& dataSend, DRT::Node& node) const
{
  const int nsd = 3;
  DRT::ParObject::AddtoPack(dataSend, node.Id());
  DRT::ParObject::AddtoPack(dataSend, LINALG::Matrix<nsd, 1>(node.X()));
  DRT::ParObject::AddtoPack(dataSend, node.Owner());
}  // end packNode



/*------------------------------------------------------------------------------------------------*
 * unpacking a node after parallel communication only with the basic nodal data                   *
 * without an underlying discretization fitting to the node's new prozessor          schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_TIMEINT_BASE::unpackNode(
    std::vector<char>::size_type& posinData, std::vector<char>& dataRecv, DRT::Node& node) const
{
  const int nsd = 3;              // dimension
  int id;                         // global id
  LINALG::Matrix<nsd, 1> coords;  // coordinates
  int owner;                      // processor

  DRT::ParObject::ExtractfromPack(posinData, dataRecv, id);
  DRT::ParObject::ExtractfromPack(posinData, dataRecv, coords);
  DRT::ParObject::ExtractfromPack(posinData, dataRecv, owner);

  if (owner == myrank_)  // real node with all data
  {
    node = *discret_->gNode(id);
  }
  else  // just id, coords and owner
  {
    double coordinates[nsd];
    for (int dim = 0; dim < nsd; dim++) coordinates[dim] = coords(dim);

    DRT::Node newNode(id, coordinates, owner);
    node = newNode;
  }  // end if correct processor
}  // end function unpackNode



/*------------------------------------------------------------------------------------------------*
 * basic XFEM time-integration constructor for standard degrees of freedom           schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
XFEM::XFLUID_STD::XFLUID_STD(
    XFEM::XFLUID_TIMEINT_BASE& timeInt,  ///< time integration base class object
    const std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>&
        reconstr_method,                      ///< reconstruction map for nodes and its dofsets
    INPAR::XFEM::XFluidTimeInt& timeIntType,  ///< type of time integration
    const Teuchos::RCP<Epetra_Vector> veln,   ///< velocity at time t^n in col map
    const double& dt,                         ///< time step size
    const bool initialize                     ///< is initialization?
    )
    : XFEM::XFLUID_TIMEINT_BASE::XFLUID_TIMEINT_BASE(timeInt),
      timeIntType_(timeIntType),
      veln_(veln),
      dt_(dt)
{
  if (initialize)
  {
    const int nsd = 3;  // dimension
    timeIntData_ = Teuchos::rcp(
        new std::vector<TimeIntData>);  // vector containing all data used for computation

    LINALG::Matrix<nsd, 1> dummyStartpoint;  // dummy startpoint for comparison
    for (int i = 0; i < nsd; i++) dummyStartpoint(i) = 777.777;

      //--------------------------------------------------------------------------------------
      // fill timeIntData_ structure with the data for the nodes which are marked for SEMILAGRANGEAN
      // reconstruction

#ifdef DEBUG_TIMINT_STD
    IO::cout << "\n\t ---------------------------------------------------------------------- ";
    IO::cout << "\n\t PREPARE SEMI-LAGRANGEAN time integration: fill timeIntData for nodes";
    IO::cout << "\n\t ---------------------------------------------------------------------- "
             << IO::endl;
#endif

    // loop over processor nodes
    for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
    {
      DRT::Node* node = discret_->lRowNode(lnodeid);

      std::map<int, std::vector<INPAR::XFEM::XFluidTimeInt>>::const_iterator it =
          reconstr_method.find(node->Id());
      if (it != reconstr_method.end())
      {
        // reconstruction methods just for marked dofsets
        for (size_t i = 0; i < (it->second).size(); i++)
        {
          if ((it->second)[i] == INPAR::XFEM::Xf_TimeInt_STD_by_SL)
          {
#ifdef DEBUG_TIMINT_STD
            IO::cout << "\t * fill timeIntData for node " << node->Id() << " and dofset " << i
                     << IO::endl;
#endif

            LINALG::Matrix<3, 1> nodedispnp(true);
            if (dispnp_ != Teuchos::null)  // is alefluid
            {
              //------------------------------------------------------- add ale disp
              // get node location vector, dirichlet flags and ownerships (discret, nds, la,
              // doDirichlet)
              std::vector<int> lm;
              std::vector<int> dofs;
              dofset_new_->Dof(dofs, node, 0);  // dofs for standard dofset
              for (int j = 0; j < 4; ++j) lm.push_back(dofs[j]);

              LINALG::Matrix<1, 1> nodepredummy(true);
              extractNodalValuesFromVector<1>(nodedispnp, nodepredummy, dispnp_, lm);
            }

            // constructor for standard computation
            timeIntData_->push_back(TimeIntData(*node,  //! node for which SL-algorithm is called
                i,  //! nds (nodal dofset) number w.r.t new interface position, for which SL-algo is
                    //! called
                LINALG::Matrix<nsd, 1>(true),  //!  velocity at point x (=x_Lagr(t^n+1))
                std::vector<LINALG::Matrix<nsd, nsd>>(oldVectors_.size(),
                    LINALG::Matrix<nsd, nsd>(
                        true)),  //! velocity gradient at point x (=x_Lagr(t^n+1))
                std::vector<LINALG::Matrix<1, nsd>>(oldVectors_.size(),
                    LINALG::Matrix<1, nsd>(
                        true)),             //! pressure gradient at point x (=x_Lagr(t^n+1))
                nodedispnp,                 //!  displacement at point x (=x_Lagr(t^n+1))
                dummyStartpoint,            // dummy-startpoint
                1,                          // searchedProcs
                0,                          // counter
                INFINITY,                   // minimal distance
                TimeIntData::predictor_));  // data created for the node
          }
        }  // nodaldofsets
      }    // some dofsets have to be reconstructed
    }      // end loop over processor nodes

    //--------------------------------------------------------------------------------------
    // compute initial points on the right side of the interface to start finding the Lagrangean
    // origin
    startpoints();

    //--------------------------------------------------------------------------------------

    // test loop if all initial startpoints have been computed
    for (std::vector<TimeIntData>::iterator data = timeIntData_->begin();
         data != timeIntData_->end(); data++)
    {
      if (data->startpoint_ == dummyStartpoint)  // startpoint unchanged
        dserror(
            "WARNING! No enriched node on one interface side found!\nThis "
            "indicates that the whole area is at one side of the interface!");
    }  // end loop over nodes
  }

  return;
}  // end constructor



/*------------------------------------------------------------------------------------------------*
 * out of order!                                                                     schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::compute(std::vector<Teuchos::RCP<Epetra_Vector>>& newRowVectorsn)
{
  dserror("Unused function! Use a function of the derived classes");
}



///*------------------------------------------------------------------------------------------------*
// * initialize data when called in a new FGI                                      winklmaier 10/11
// *
// *------------------------------------------------------------------------------------------------*/
// void XFEM::XFLUID_STD::importNewFGIData(
//    const Teuchos::RCP<DRT::Discretization> discret,
//    const Teuchos::RCP<XFEM::DofManager> newdofman,
//    const Teuchos::RCP<COMBUST::FlameFront> flamefront,
//    const Epetra_Map& newdofrowmap,
//    const std::map<DofKey, DofGID>& newNodalDofRowDistrib)
//{
//  discret_ = discret;
//  newdofman_ = newdofman;
//  phinpi_ = phinp_;
//  phinp_ = flamefront->Phinp();
//  flamefront_ = flamefront;
//  newinterfacehandle_ = flamefront->InterfaceHandle();
//  newdofrowmap_ = newdofrowmap;
//  newNodalDofRowDistrib_ = newNodalDofRowDistrib;
//  return;
//}



/*------------------------------------------------------------------------------------------------*
 * identify an element containing a point and additional data                        schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::elementSearch(
    DRT::Element*& ele,        /// pointer to element if point lies in a found element
    LINALG::Matrix<3, 1>& x,   /// global coordiantes of point
    LINALG::Matrix<3, 1>& xi,  /// determined local coordinates w.r.t ele
    bool& found                /// is element found?
    ) const
{
  // REMARK: if ele!= NULL, then check that element first, before loop all row elements

  int startid;  // local row element id
  if (ele == NULL)
    startid = 0;  // start with first local row element
  else
    startid = -1;  // pseudo-id so that id+1 will be 0 (for additional element check, if ele!=NULL)

  DRT::Element* currele = NULL;  // current element

  // loop over elements
  for (int ieleid = startid; ieleid < discret_->NumMyRowElements(); ieleid++)
  {
    // if ele != NULL additional check
    // first it should be checked if it is fitting
    if (ieleid == -1)
    {
      currele = ele;
      ele = NULL;  // reset ele, ele will be set if an element is found finally
    }
    else
      currele = discret_->lRowElement(ieleid);

    callXToXiCoords(currele, x, xi, "dispn", found);

    if (found)
    {
      ele = currele;
      break;
    }
  }  // end loop over processor elements

  return;
}  // end elementSearch



/*------------------------------------------------------------------------------------------------*
 * interpolate velocity and derivatives for a point in an element                    schott 06/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::getGPValues(DRT::Element* ele,  ///< pointer to element
    LINALG::Matrix<3, 1>& xi,                          ///< local coordinates of point w.r.t element
    std::vector<int>& nds,                      ///< nodal dofset of point for elemental nodes
    XFEM::XFEMDofSet& dofset,                   ///< XFEM dofset
    LINALG::Matrix<3, 1>& vel,                  ///< determine velocity at point
    LINALG::Matrix<3, 3>& vel_deriv,            ///< determine velocity derivatives at point
    double& pres,                               ///< pressure
    LINALG::Matrix<1, 3>& pres_deriv,           ///< pressure gradient
    Teuchos::RCP<const Epetra_Vector> vel_vec,  ///< vector used for interpolating at gp
    bool compute_deriv                          ///< shall derivatives be computed?
    ) const
{
  switch (ele->Shape())
  {
    case DRT::Element::hex8:
      getGPValuesT<DRT::Element::hex8>(
          ele, xi, nds, dofset, vel, vel_deriv, pres, pres_deriv, vel_vec, compute_deriv);
      break;
    case DRT::Element::hex20:
      getGPValuesT<DRT::Element::hex20>(
          ele, xi, nds, dofset, vel, vel_deriv, pres, pres_deriv, vel_vec, compute_deriv);
      break;
    case DRT::Element::tet4:
      getGPValuesT<DRT::Element::tet4>(
          ele, xi, nds, dofset, vel, vel_deriv, pres, pres_deriv, vel_vec, compute_deriv);
      break;
    default:
      dserror("add your 3D distype here!");
      break;
  }  // end switch

  return;
}  // end function getGPValues


//! interpolate velocity and derivatives for a point in an element
template <DRT::Element::DiscretizationType DISTYPE>
void XFEM::XFLUID_STD::getGPValuesT(DRT::Element* ele,  ///< pointer to element
    LINALG::Matrix<3, 1>& xi,                   ///< local coordinates of point w.r.t element
    std::vector<int>& nds,                      ///< nodal dofset of point for elemental nodes
    XFEM::XFEMDofSet& dofset,                   ///< XFEM dofset
    LINALG::Matrix<3, 1>& vel,                  ///< determine velocity at point
    LINALG::Matrix<3, 3>& vel_deriv,            ///< determine velocity derivatives at point
    double& pres,                               ///< pressure
    LINALG::Matrix<1, 3>& pres_deriv,           ///< pressure gradient
    Teuchos::RCP<const Epetra_Vector> vel_vec,  ///< vector used for interpolating at gp
    bool compute_deriv                          ///< shall derivatives be computed?
    ) const
{
  if (vel_vec == Teuchos::null) dserror("vector is not filled");

  const int nsd = 3;  // dimension

  const int numdofpernode = nsd + 1;
  const int numnode =
      DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;  // number of element nodes

  //-------------------------------------------------------
  // initialization
  LINALG::Matrix<numnode, 1> shapeFcn(true);         /// shape function at point
  LINALG::Matrix<nsd, numnode> shapeFcnDeriv(true);  /// xyz shape derivatives at point
  LINALG::Matrix<nsd, nsd> xji(true);                /// inverse of jacobian

  //-------------------------------------------------------
  // get element location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
  if ((int)(nds.size()) != numnode)
    dserror("size of nds-vector (%d) != numnode(%d)", nds.size(), numnode);

  std::vector<int> lm;

  for (int inode = 0; inode < numnode; inode++)
  {
    DRT::Node* node = ele->Nodes()[inode];
    std::vector<int> dofs;
    dofset.Dof(dofs, node, nds[inode]);

    for (int j = 0; j < 4; ++j)
    {
      lm.push_back(dofs[j]);
    }
  }

  //-------------------------------------------------------
  // get element-wise velocity/pressure field
  LINALG::Matrix<nsd, numnode> evel(true);
  LINALG::Matrix<numnode, 1> epre(true);


  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*vel_vec, mymatrix, lm);

  for (int inode = 0; inode < numnode; ++inode)  // number of nodes
  {
    for (int idim = 0; idim < nsd; ++idim)  // number of dimensions
    {
      (evel)(idim, inode) = mymatrix[idim + (inode * numdofpernode)];
    }
    (epre)(inode, 0) = mymatrix[nsd + (inode * numdofpernode)];
  }


  //-------------------------------------------------------
  evalShapeAndDeriv<numnode, DISTYPE>(ele, xi, xji, shapeFcn, shapeFcnDeriv, compute_deriv);

  //-------------------------------------------------------
  // interpolate velocity and pressure values at starting point
  vel.Clear();  // set to zero
  vel.Multiply(evel, shapeFcn);

  pres = epre.Dot(shapeFcn);

  if (compute_deriv)
  {
    vel_deriv.MultiplyNT(evel, shapeFcnDeriv);
    pres_deriv.MultiplyTT(epre, shapeFcnDeriv);
  }

  return;
};  // end template getGPValuesT



/*------------------------------------------------------------------------------------------------*
 * Compute starting values for nodes marked for semi-lagrangian algo                 schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::startpoints()
{
#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t ---------------------------------------------------------------------- ";
  IO::cout << "\n\t compute initial start points for finding the Lagrangean origin";
  IO::cout << "\n\t ---------------------------------------------------------------------- "
           << IO::endl;
#endif

  // REMARK: we do not need parallel communication for start point values because
  // structural surface is ghosted on all procs

  // loop over nodes which changed interface side
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->state_ == TimeIntData::basicStd_)  // correct state
    {
      // find a start approximation for Semi-Lagrange backtracking
      // this start approximation can be found by tracking back the projection of current node
      // along structural movement
      ProjectAndTrackback(*data);

    }  // end if correct state
  }    // end loop over nodes which changed interface side


  return;
}  // end startpoints()



/*------------------------------------------------------------------------------------------------*
 * Project the current point onto the structural surface, track it back along structural         *
 * movement and project it back into the fluid domain                                schott 06/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::ProjectAndTrackback(TimeIntData& data)
{
#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t * ProjectAndTrackback for node " << data.node_.Id() << IO::endl;
#endif

  const std::string state = "idispnp";

  const int nsd = 3;
  LINALG::Matrix<nsd, 1> newNodeCoords(
      data.node_.X());  // coords of endpoint of Lagrangian characteristics
  for (int i = 0; i < nsd; ++i) newNodeCoords(i) += data.dispnp_(i);

  // determine the smallest distance of node at t^(n+1) to the current side elements

  GEO::CUT::Node* n_new = wizard_new_->GetNode(data.node_.Id());

  //------------------------------------
  // find all involved nodes and sides (edges) for computing distance to node

  std::set<int> points;
  std::set<int> sides;

  // is there a nodehandle to get information about cutting sides near this node?
  if (n_new != NULL)
  {
    //--------------------------------------------------------
    // get involved side ids for projection and distance computation
    //--------------------------------------------------------

    // set of side-ids involved in cutting the current connection of volumecells at t^(n+1)
    std::map<int, std::vector<GEO::CUT::BoundaryCell*>> bcells_new;

    // the std-set
    const std::set<GEO::CUT::plain_volumecell_set, GEO::CUT::Cmp>& cell_set_new =
        n_new->GetNodalDofSet(data.nds_np_)->VolumeCellComposite();  // always the standard dofset

    // get all side-ids w.r.t to all volumecells contained in current new set around the current
    // node
    for (std::set<GEO::CUT::plain_volumecell_set>::const_iterator adj_eles = cell_set_new.begin();
         adj_eles != cell_set_new.end(); adj_eles++)
    {
      const GEO::CUT::plain_volumecell_set ele_vc = *adj_eles;

      for (GEO::CUT::plain_volumecell_set::const_iterator vcs = ele_vc.begin(); vcs != ele_vc.end();
           vcs++)
      {
        GEO::CUT::VolumeCell* vc = *vcs;

        // get sides involved in creation boundary cells (std::map<sideId,bcells>)
        vc->GetBoundaryCells(bcells_new);
      }
    }

    //--------------------------------------------------------
    // distance computation w.r.t surface involved in cutting the adjacent elements(n_new to
    // sides/edges/nodes)
    //--------------------------------------------------------

    // loop bcs and extract sides and nodes
    for (std::map<int, std::vector<GEO::CUT::BoundaryCell*>>::iterator bcells_it =
             bcells_new.begin();
         bcells_it != bcells_new.end(); bcells_it++)
    {
      int sid = bcells_it->first;
      sides.insert(sid);

      // get the side
      DRT::Element* side = boundarydis_->gElement(sid);
      int numnode = side->NumNode();
      const int* side_nodes = side->NodeIds();

      // get the nodes
      for (int i = 0; i < numnode; i++)
      {
        points.insert(side_nodes[i]);
      }
    }
  }

  if ((points.size() == 0 and sides.size() == 0) or n_new == NULL)
  {
    // node does not carry a node handle -> Semilagrangean algorithm for nodes with too much
    // structural movement
    // TODO: can this improved?!
    // -> compute projection for all sides and nodes, that's really expensive!

    // loop all column sides and extract sides and nodes
    for (int i = 0; i < boundarydis_->NumMyColElements(); i++)
    {
      // get the side
      DRT::Element* side = boundarydis_->lColElement(i);
      int numnode = side->NumNode();
      const int* side_nodes = side->NodeIds();

      sides.insert(side->Id());

      // get the nodes
      for (int i = 0; i < numnode; i++)
      {
        points.insert(side_nodes[i]);
      }
    }
  }

  if (points.size() == 0 and sides.size() == 0)
    dserror("there are no cutting sides for node %d, that cannot be anymore", data.node_.Id());

  //-------------------------------------
  // initialize data holding information about minimal distance

  // smallest distance
  LINALG::Matrix<3, 1> proj_x_np(true);  ///< projected point at t^(n+1)
  LINALG::Matrix<3, 1> proj_x_n(
      true);  ///< projected point at t^n (tracked back along structural movement)
  LINALG::Matrix<3, 1> start_point(true);  ///< final start point for SemiLagrange algo
  int proj_sid = -1;                       ///< id of side that contains the projected point
  std::map<std::vector<int>, std::vector<int>>
      proj_lineid;  ///< std::map< sorted nids, global side IDs >
  // smallest distance w.r.t side
  LINALG::Matrix<2, 1> proj_xi_side(
      true);  ///< local coordinates of projected point if projection w.r.t side
  // smallest distance w.r.t line
  std::map<std::vector<int>, std::vector<double>>
      proj_xi_line;  ///< std::map<sorted nids,local line coordinates w.r.t lines of different
                     ///< sides>
  std::map<std::vector<int>, std::vector<int>> proj_nid_line;  ///< std::map<sorted nids, sideids>
  // smallest distance w.r.t point
  int proj_nid_np = -1;  ///< nid of projected point if projection w.r.t point (structural node)

  double min_dist = INFINITY;  ///< minimal distance

  data.proj_ = TimeIntData::failed_;  ///< projection method yielding minimal distance

  //----------------------------------------------------------------------
  // check distance to sides and its lines
  //----------------------------------------------------------------------
  for (std::set<int>::iterator it = sides.begin(); it != sides.end(); it++)
  {
    //-------------------------------------------------

    DRT::Element* side = boundarydis_->gElement(*it);

    CallProjectOnSide(
        side, state, newNodeCoords, proj_sid, min_dist, proj_x_np, proj_xi_side, data.proj_);

    //------------------------------------
    // edges are checked twice (w.r.t both sides)
    //------------------------------------

    // loop all edges of current side
    std::vector<Teuchos::RCP<DRT::Element>> lines = side->Lines();

    int line_count = 0;

    for (std::vector<Teuchos::RCP<DRT::Element>>::iterator line_it = lines.begin();
         line_it != lines.end(); line_it++)
    {
      CallProjectOnLine(side, &(**line_it),
          line_count,  ///< local line id w.r.t side element
          state,
          newNodeCoords,  ///< node coordinates of point that has to be projected
          min_dist,       ///< minimal distance, potentially updated
          proj_x_np,      ///< projection of point on this side
          proj_xi_line,   ///< std::map<side ID, local coordinates of projection of point w.r.t to
                          ///< this line>
          proj_lineid,    ///< std::map<side ID, local line id>
          proj_nid_line,  ///< std::map<side ID, vec<line Ids>>
          proj_sid,       ///< id of side that contains the projected point
          data.proj_      ///< reference to data
      );

      line_count++;
    }  // loop lines of side
  }    // loop sides


  //------------------------------------
  // check minimal distance w.r.t points (structural nodes)
  //------------------------------------
  // loop all points
  for (std::set<int>::iterator it = points.begin(); it != points.end(); it++)
  {
    // compute distance to following point
    DRT::Node* node = boundarydis_->gNode(*it);

    CallProjectOnPoint(node,  ///< pointer to node
        state,
        newNodeCoords,  ///< node coordinates of point that has to be projected
        min_dist,       ///< minimal distance, potentially updated
        proj_nid_np,    ///< nid id of projected point on surface
        proj_x_np,      ///< projection of point on this side
        proj_sid,       ///< id of side that contains the projected point
        proj_xi_line,  ///< std::map<side ID, local coordinates of projection of point w.r.t to this
                       ///< line>
        proj_lineid,   ///< std::map<side ID, local line id>
        proj_nid_line,  ///< std::map<side ID, vec<line Ids>>
        data.proj_      ///< reference to data
    );


  }  // loop points


  if (data.proj_ == TimeIntData::failed_)
    dserror("projection of node %d not successful!", n_new->Id());

#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t => Projection of node lies on (side=0, line=1, point=2, failed=3): "
           << data.proj_ << " with distance " << min_dist << IO::endl;
#endif


  // track back the projected points on structural surface from t^(n+1)->t^n
  if (data.proj_ == TimeIntData::onSide_)
  {
    DRT::Element* side = boundarydis_->gElement(proj_sid);

    if (side == NULL) dserror("side with id %d not found ", proj_sid);

    // side geometry at initial state t^0
    const int numnodes = side->NumNode();
    DRT::Node** nodes = side->Nodes();
    Epetra_SerialDenseMatrix side_xyze(3, numnodes);
    for (int i = 0; i < numnodes; ++i)
    {
      const double* x = nodes[i]->X();
      std::copy(x, x + 3, &side_xyze(0, i));
    }

    DRT::Element::LocationArray cutla(1);
    side->LocationVector(*boundarydis_, cutla, false);

    ComputeStartPoint_Side(
        side, side_xyze, cutla[0].lm_, proj_xi_side, min_dist, proj_x_n, start_point);
  }
  else if (data.proj_ == TimeIntData::onLine_)
  {
    // check if both lines are identical
    for (std::map<std::vector<int>, std::vector<int>>::iterator line = proj_nid_line.begin();
         line != proj_nid_line.end(); line++)
    {
      std::vector<int>& sides = line->second;

      std::vector<int>& local_lineIds = proj_lineid.find(line->first)->second;
      std::vector<double>& local_xi_line = proj_xi_line.find(line->first)->second;


      int sid_1 = -1;
      int sid_2 = -1;

      if (sides.size() <= 0)
      {
        dserror(
            "there must be at least one found side adjacent to this line, but there are %d sides",
            sides.size());
      }
      else if (sides.size() == 1)
      {
        sid_1 = sides[0];

        // find the second neighboring side
        DRT::Element* side_1 = boundarydis_->gElement(sid_1);

        DRT::Node** nodes = side_1->Nodes();

        std::set<int> neighbor_sides;

        for (int i = 0; i < side_1->NumNode(); i++)
        {
          DRT::Element** surr_sides = nodes[i]->Elements();
          for (int s = 0; s < nodes[i]->NumElement(); s++)
          {
            neighbor_sides.insert(surr_sides[s]->Id());
          }
        }

        // check which neighbor side is the right one
        for (std::set<int>::iterator s = neighbor_sides.begin(); s != neighbor_sides.end(); s++)
        {
          DRT::Element* side_tmp = boundarydis_->gElement(*s);

          std::vector<Teuchos::RCP<DRT::Element>> lines_tmp = side_tmp->Lines();

          for (std::vector<Teuchos::RCP<DRT::Element>>::iterator line_it = lines_tmp.begin();
               line_it != lines_tmp.end(); line_it++)
          {
            DRT::Node** line_nodes = (*line_it)->Nodes();
            std::vector<int> line_nids;
            for (int i = 0; i < (*line_it)->NumNode(); ++i)
            {
              line_nids.push_back(line_nodes[i]->Id());
            }

            sort(line_nids.begin(), line_nids.end());

            if (line_nids == line->first)
            {
              // side found
              if (sid_2 == -1 and side_tmp->Id() != sid_1) sid_2 = side_tmp->Id();
              break;
            }
          }
        }
        if (sid_2 == -1)
        {
#ifdef DEBUG_TIMINT_STD
          std::cout
              << "\t\t\t\tneighbor side not found, current line seems to be a line at the boundary"
              << std::endl;
          std::cout << "\t\t\t\t -> the same distance should have been found for a side ? "
                    << std::endl;
          std::cout << "line nodeids: " << std::endl;
          for (int i = 0; i < (int)(line->first).size(); i++)
            std::cout << "\t nid: " << (line->first)[i] << std::endl;
          std::cout << "\t\t\t\In case your simulation is a pseudo 2d case, every thing is fine "
                       "and we will just use"
                    << std::endl;
          std::cout << "\t\t\t\tthe normal on the first side!" << std::endl;
#endif
        }
      }
      else if (sides.size() == 2)
      {
        sid_1 = sides[0];
        sid_2 = sides[1];
      }
      else
        dserror("non valid number of sides adjacent to line");

      //---------------------------------------------------------
      // side 1
      DRT::Element* side_1 = boundarydis_->gElement(sid_1);

      // side geometry at initial state t^0
      const int numnodes_1 = side_1->NumNode();
      DRT::Node** nodes_1 = side_1->Nodes();
      Epetra_SerialDenseMatrix side_xyze_1(3, numnodes_1);
      for (int i = 0; i < numnodes_1; ++i)
      {
        const double* x = nodes_1[i]->X();
        std::copy(x, x + 3, &side_xyze_1(0, i));
      }

      DRT::Element::LocationArray cutla_1(1);
      side_1->LocationVector(*boundarydis_, cutla_1, false);

      //---------------------------------------------------------
      // side 2
      DRT::Element* side_2;
      Epetra_SerialDenseMatrix side_xyze_2;
      DRT::Element::LocationArray cutla_2(1);
      if (sid_2 != -1)
      {
        side_2 = boundarydis_->gElement(sid_2);
        // side geometry at initial state t^0
        const int numnodes_2 = side_2->NumNode();
        DRT::Node** nodes_2 = side_2->Nodes();
        side_xyze_2.Shape(3, numnodes_2);
        for (int i = 0; i < numnodes_2; ++i)
        {
          const double* x = nodes_2[i]->X();
          std::copy(x, x + 3, &side_xyze_2(0, i));
        }
        side_2->LocationVector(*boundarydis_, cutla_2, false);
      }
      else
      {
        side_2 = NULL;
      }
      //---------------------------------------------------------
      // line geometry at initial state t^0
      std::vector<Teuchos::RCP<DRT::Element>> lines = side_1->Lines();

      Teuchos::RCP<DRT::Element> line_ele = lines[local_lineIds[0]];

      call_get_projxn_Line(side_1, &*line_ele, proj_x_n, local_xi_line[0]);


      ComputeStartPoint_Line(side_1,  ///< pointer to side element
          side_xyze_1,                ///< side's node coordinates
          side_2,                     ///< pointer to side element
          side_xyze_2,                ///< side's node coordinates
          cutla_1[0].lm_,             ///< local map
          cutla_2[0].lm_,             ///< local map
          min_dist,                   ///< distance from point to its projection
          proj_x_n,                   ///< projected point at t^n
          start_point                 ///< final start point
      );
    }
  }
  else if (data.proj_ == TimeIntData::onPoint_)
  {
    std::vector<DRT::Element*> surr_sides;                  ///< pointer to side element
    std::vector<Epetra_SerialDenseMatrix> surr_sides_xyze;  ///< side's node coordinates
    std::vector<std::vector<int>> surr_sides_lm;            ///< local map


    // get the node
    DRT::Node* node = boundarydis_->gNode(proj_nid_np);

    // get all surrounding sides
    DRT::Element** node_sides = node->Elements();

    for (int s = 0; s < node->NumElement(); s++)
    {
      // side geometry at initial state t^0
      DRT::Element* side = node_sides[s];
      const int numnodes = side->NumNode();
      DRT::Node** nodes = side->Nodes();
      Epetra_SerialDenseMatrix side_xyze(3, numnodes);
      for (int i = 0; i < numnodes; ++i)
      {
        const double* x = nodes[i]->X();
        std::copy(x, x + 3, &side_xyze(0, i));
      }

      DRT::Element::LocationArray cutla(1);
      side->LocationVector(*boundarydis_, cutla, false);

      surr_sides.push_back(side);
      surr_sides_xyze.push_back(side_xyze);
      surr_sides_lm.push_back(cutla[0].lm_);
    }

    //----------------------------------------------------

    // its point geometry
    Epetra_SerialDenseMatrix point_xyze(3, 1);

    const double* x = node->X();
    std::copy(x, x + 3, &point_xyze(0, 0));

    std::vector<int> lm = boundarydis_->Dof(0, node);

    // get the coordinates of the projection point at time tn
    const std::string state = "idispn";

    Teuchos::RCP<const Epetra_Vector> matrix_state = boundarydis_->GetState(state);
    if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

    // extract local values of the global vectors
    std::vector<double> mymatrix(lm.size());
    DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

    // add the displacement of the interface
    for (int idim = 0; idim < 3; ++idim)  // number of dimensions
    {
      point_xyze(idim, 0) += mymatrix[idim];  // attention! disp state vector has 3+1 dofs for
                                              // displacement (the same as for (u,p))
    }

    // fill LINALG matrix
    for (int idim = 0; idim < 3; ++idim)  // number of dimensions
    {
      proj_x_n(idim) = point_xyze(idim,
          0);  // attention! disp state vector has 3+1 dofs for displacement (the same as for (u,p))
    }

    //----------------------------------------------------


    ComputeStartPoint_AVG(surr_sides,  ///< pointer to side element
        surr_sides_xyze,               ///< side's node coordinates
        surr_sides_lm,                 ///< local map
        min_dist,                      ///< distance from point to its projection
        proj_x_n,                      ///< projected point at t^n
        start_point                    ///< final start point
    );
  }
  else
    dserror("unknown projection method for tracking back the structural surface point");


  /// set the data
  if (data.proj_ != TimeIntData::failed_)
  {
    // fill data
    //  data.startpoint_.Update(1.0,newNodeCoords,-dt_,nodevel);
    data.startpoint_.Update(1.0, start_point, 0.0);
    data.initialpoint_.Update(1.0, start_point, 0.0);
    data.initial_eid_ = -1;
    data.initial_ele_owner_ = -1;
    data.dMin_ = min_dist;
  }


#if (0)
  // GMSH debug output

  int myrank = 0;  // discret.Comm().MyPID();

  if (myrank == 0) IO::cout << "\n\t ... writing Gmsh output for startvalues ...\n" << IO::flush;

  int step = 0;

  const int step_diff = 100;
  bool screen_out = true;

  // output for Element and Node IDs
  std::ostringstream filename_base_startvalues;
  filename_base_startvalues << "TIMINT"
                            << "_startvalues_" << n_new->Id() << "_";
  const std::string filename_startvalues = IO::GMSH::GetNewFileNameAndDeleteOldFiles(
      filename_base_startvalues.str(), step, step_diff, screen_out, myrank);
  if (true) std::cout << std::endl;
  std::ofstream gmshfilecontent_startvalues(filename_startvalues.c_str());
  gmshfilecontent_startvalues.setf(std::ios::scientific, std::ios::floatfield);
  gmshfilecontent_startvalues.precision(16);

  gmshfilecontent_startvalues << "View \""
                              << "SL "
                              << "node " << n_new->Id() << "\" {\n";

  {
    // 1. line ( x_np -> x_np projected on side )
    gmshfilecontent_startvalues << "SL(" << newNodeCoords(0) << "," << newNodeCoords(1) << ","
                                << newNodeCoords(2) << "," << proj_x_np(0) << "," << proj_x_np(1)
                                << "," << proj_x_np(2) << ")" << std::endl;
    gmshfilecontent_startvalues << "{" << 0 << "," << 1 << "};" << std::endl;
    // 2. line ( x_np projected on side -> x_n projected on side )
    gmshfilecontent_startvalues << "SL(" << proj_x_np(0) << "," << proj_x_np(1) << ","
                                << proj_x_np(2) << "," << proj_x_n(0) << "," << proj_x_n(1) << ","
                                << proj_x_n(2) << ")" << std::endl;
    gmshfilecontent_startvalues << "{" << 1 << "," << 2 << "};" << std::endl;
    // 3. line ( x_n projected on side -> start_point)
    gmshfilecontent_startvalues << "SL(" << proj_x_n(0) << "," << proj_x_n(1) << "," << proj_x_n(2)
                                << "," << start_point(0) << "," << start_point(1) << ","
                                << start_point(2) << ")" << std::endl;
    gmshfilecontent_startvalues << "{" << 2 << "," << 3 << "};" << std::endl;
  }

  gmshfilecontent_startvalues << "};\n";

#endif
}


bool XFEM::XFLUID_STD::FindNearestSurfPoint(
    LINALG::Matrix<3, 1>& x,       ///< coords of point to be projected
    LINALG::Matrix<3, 1>& proj_x,  ///< coords of projected point
    GEO::CUT::VolumeCell* vc,      ///< volumcell on that's cut-surface we want to project
    const std::string state        ///< state n or np?
)
{
  if (vc == NULL) dserror("do not call FindNearestSurfPoint with Null-Pointer for Volumecell");

  //--------------------------------------------------------
  // get involved side ids for projection and distance computation
  //--------------------------------------------------------

  // set of side-ids involved in cutting the current connection of volumecells at t^(n+1)
  std::map<int, std::vector<GEO::CUT::BoundaryCell*>> bcells_new;

  // get sides involved in creation boundary cells (std::map<sideId,bcells>)
  vc->GetBoundaryCells(bcells_new);



  std::set<int> points;
  std::set<int> sides;

  // loop bcs and extract sides and nodes
  for (std::map<int, std::vector<GEO::CUT::BoundaryCell*>>::iterator bcells_it = bcells_new.begin();
       bcells_it != bcells_new.end(); bcells_it++)
  {
    int sid = bcells_it->first;
    sides.insert(sid);

    // get the side
    DRT::Element* side = boundarydis_->gElement(sid);
    int numnode = side->NumNode();
    const int* side_nodes = side->NodeIds();

    // get the nodes
    for (int i = 0; i < numnode; i++)
    {
      points.insert(side_nodes[i]);
    }
  }


  return ProjectToSurface(x, proj_x, state, points, sides);
}


/*------------------------------------------------------------------------------------------------*
 * Project the current point onto the structural surface given by point and side Ids schott 06/12 *
 *------------------------------------------------------------------------------------------------*/
bool XFEM::XFLUID_STD::ProjectToSurface(
    LINALG::Matrix<3, 1>& x,       ///< coords of point to be projected
    LINALG::Matrix<3, 1>& proj_x,  ///< coords of projected point
    const std::string state,       ///< state n or np?
    std::set<int>& points,         ///< node Ids of surface for that distance has to be computed
    std::set<int>&
        sides  ///< side Ids of surface for that and it's lines the distances have to be computed
)
{
#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t * FindNearestSurfPoint for point " << x << IO::endl;
#endif

  if (points.size() == 0 and sides.size() == 0)
    dserror(
        "there are no cutting sides around current point (%d,%d,%d), Projection on surface not "
        "possible here",
        x(0), x(1), x(2));

  //-------------------------------------
  // initialize data holding information about minimal distance

  // smallest distance
  int proj_sid = -1;  ///< id of side that contains the projected point
  std::map<std::vector<int>, std::vector<int>>
      proj_lineid;  ///< std::map< sorted nids, global side IDs >
  // smallest distance w.r.t side
  LINALG::Matrix<2, 1> proj_xi_side(
      true);  ///< local coordinates of projected point if projection w.r.t side
  // smallest distance w.r.t line
  std::map<std::vector<int>, std::vector<double>>
      proj_xi_line;  ///< std::map<sorted nids,local line coordinates w.r.t lines of different
                     ///< sides>
  std::map<std::vector<int>, std::vector<int>> proj_nid_line;  ///< std::map<sorted nids, sideids>
  // smallest distance w.r.t point
  int proj_nid_np = -1;  ///< nid of projected point if projection w.r.t point (structural node)

  double min_dist = INFINITY;  ///< minimal distance

  TimeIntData::projection proj =
      TimeIntData::failed_;  ///< projection method yielding minimal distance

  //----------------------------------------------------------------------
  // check distance to sides and its lines
  //----------------------------------------------------------------------
  for (std::set<int>::iterator it = sides.begin(); it != sides.end(); it++)
  {
    //-------------------------------------------------

    DRT::Element* side = boundarydis_->gElement(*it);

    CallProjectOnSide(side, state, x, proj_sid, min_dist, proj_x, proj_xi_side, proj);

    //------------------------------------
    // edges are checked twice (w.r.t both sides)
    //------------------------------------

    // loop all edges of current side
    std::vector<Teuchos::RCP<DRT::Element>> lines = side->Lines();

    int line_count = 0;

    for (std::vector<Teuchos::RCP<DRT::Element>>::iterator line_it = lines.begin();
         line_it != lines.end(); line_it++)
    {
      CallProjectOnLine(side, &(**line_it),
          line_count,  ///< local line id w.r.t side element
          state,
          x,              ///< node coordinates of point that has to be projected
          min_dist,       ///< minimal distance, potentially updated
          proj_x,         ///< projection of point on this side
          proj_xi_line,   ///< std::map<side ID, local coordinates of projection of point w.r.t to
                          ///< this line>
          proj_lineid,    ///< std::map<side ID, local line id>
          proj_nid_line,  ///< std::map<side ID, vec<line Ids>>
          proj_sid,       ///< id of side that contains the projected point
          proj            ///< reference to data
      );

      line_count++;
    }  // loop lines of side
  }    // loop sides


  //------------------------------------
  // check minimal distance w.r.t points (structural nodes)
  //------------------------------------
  // loop all points
  for (std::set<int>::iterator it = points.begin(); it != points.end(); it++)
  {
    // compute distance to following point
    DRT::Node* node = boundarydis_->gNode(*it);

    CallProjectOnPoint(node,  ///< pointer to node on which we want to project
        state,
        x,             ///< node coordinates of point that has to be projected
        min_dist,      ///< minimal distance, potentially updated
        proj_nid_np,   ///< nid id of projected point on surface
        proj_x,        ///< projection of point on this side
        proj_sid,      ///< id of side that contains the projected point
        proj_xi_line,  ///< std::map<side ID, local coordinates of projection of point w.r.t to this
                       ///< line>
        proj_lineid,   ///< std::map<side ID, local line id>
        proj_nid_line,  ///< std::map<side ID, vec<line Ids>>
        proj            ///< reference to data
    );


  }  // loop points


  if (proj == TimeIntData::failed_)
    dserror("projection of point (%d,%d,%d) not successful!", x(0), x(1), x(2));

#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t => Projection of node lies on (side=0, line=1, point=2, failed=3): " << proj
           << " with distance " << min_dist << " coords: " << proj_x << IO::endl;
#endif

  return (proj != TimeIntData::failed_);
}



void XFEM::XFLUID_STD::ComputeStartPoint_Side(DRT::Element* side,  ///< pointer to side element
    Epetra_SerialDenseMatrix& side_xyze,                           ///< side's node coordinates
    const std::vector<int>& lm,                                    ///< local map
    LINALG::Matrix<2, 1>& xi_side,     ///< local coordinates of projected point w.r.t side
    double& dist,                      ///< distance from point to its projection
    LINALG::Matrix<3, 1>& proj_x_n,    ///< projected point at t^n
    LINALG::Matrix<3, 1>& start_point  ///< final start point
)
{
  LINALG::Matrix<3, 1> normal(true);

  callgetNormalSide_tn(side, normal, side_xyze, lm, proj_x_n, xi_side);

  // map point into fluid domain along normal vector
  start_point.Update(1.0, proj_x_n, dist, normal);

  return;
}

void XFEM::XFLUID_STD::ComputeStartPoint_Line(DRT::Element* side1,  ///< pointer to side element
    Epetra_SerialDenseMatrix& side1_xyze,                           ///< side's node coordinates
    DRT::Element* side2,                                            ///< pointer to side element
    Epetra_SerialDenseMatrix& side2_xyze,                           ///< side's node coordinates
    const std::vector<int>& lm1,                                    ///< local map
    const std::vector<int>& lm2,                                    ///< local map
    double& dist,                      ///< distance from point to its projection
    LINALG::Matrix<3, 1>& proj_x_n,    ///< projected point at t^n
    LINALG::Matrix<3, 1>& start_point  ///< final start point
)
{
  LINALG::Matrix<3, 1> normal_avg(true);  // averaged normal vector
  LINALG::Matrix<3, 1> normal1(true);
  LINALG::Matrix<3, 1> normal2(true);

  LINALG::Matrix<2, 1> xi_side1(true);
  LINALG::Matrix<2, 1> xi_side2(true);


  LINALG::Matrix<3, 1> xi_1_avg(true);
  LINALG::Matrix<3, 1> xi_2_avg(true);

  LINALG::Matrix<3, 1> proj_x_n_dummy1(true);


  for (int i = 0; i < side1->NumNode(); i++)
    xi_1_avg.Update(1.0, DRT::UTILS::getNodeCoordinates(i, side1->Shape()), 1.0);

  xi_1_avg.Scale(1.0 / side1->NumNode());

  xi_side1(0, 0) = xi_1_avg(0, 0);
  xi_side1(1, 0) = xi_1_avg(1, 0);

  callgetNormalSide_tn(side1, normal1, side1_xyze, lm1, proj_x_n_dummy1, xi_side1);

  normal_avg.Update(1.0, normal1, 1.0);

  if (side2 != NULL)  // in case we have side2, use averaged normal
  {
    for (int i = 0; i < side2->NumNode(); i++)
      xi_2_avg.Update(1.0, DRT::UTILS::getNodeCoordinates(i, side2->Shape()), 1.0);

    xi_2_avg.Scale(1.0 / side2->NumNode());

    xi_side2(0, 0) = xi_2_avg(0, 0);
    xi_side2(1, 0) = xi_2_avg(1, 0);

    callgetNormalSide_tn(side2, normal2, side2_xyze, lm2, proj_x_n_dummy1, xi_side2);

    normal_avg.Update(1.0, normal2, 1.0);
  }

  normal_avg.Scale(1.0 / normal_avg.Norm2());

  // map point into fluid domain along normal vector
  start_point.Update(1.0, proj_x_n, dist, normal_avg);

  return;
}


void XFEM::XFLUID_STD::ComputeStartPoint_AVG(
    const std::vector<DRT::Element*>& sides,            ///< pointer to side element
    std::vector<Epetra_SerialDenseMatrix>& sides_xyze,  ///< side's node coordinates
    const std::vector<std::vector<int>>& sides_lm,      ///< local map
    double& dist,                                       ///< distance from point to its projection
    LINALG::Matrix<3, 1>& proj_x_n,                     ///< projected point at t^n
    LINALG::Matrix<3, 1>& start_point                   ///< final start point
)
{
  if (sides.size() != sides_xyze.size() or sides.size() != sides_lm.size())
    dserror("not equal number of sides, xyze-coordinates or lm-vectors ");

  LINALG::Matrix<3, 1> normal_avg(true);  // averaged normal vector

  ///
  for (std::vector<DRT::Element*>::const_iterator it = sides.begin(); it != sides.end(); it++)
  {
    DRT::Element* side = *it;

    LINALG::Matrix<2, 1> side_center(true);

    LINALG::Matrix<3, 1> local_node_coord(true);

    // get the side-center
    for (int i = 0; i < side->NumNode(); i++)
    {
      local_node_coord = DRT::UTILS::getNodeCoordinates(i, side->Shape());
      side_center(0) += local_node_coord(0);
      side_center(1) += local_node_coord(1);
    }

    side_center.Scale(1.0 / side->NumNode());

    // get the normal at the side center
    LINALG::Matrix<3, 1> proj_x_n_dummy1(true);
    LINALG::Matrix<3, 1> side_normal(true);

    callgetNormalSide_tn(side, side_normal, sides_xyze[it - sides.begin()],
        sides_lm[it - sides.begin()], proj_x_n_dummy1, side_center);

    normal_avg += side_normal;
  }

  normal_avg.Scale(1.0 / normal_avg.Norm2());

  // map point into fluid domain along normal vector
  start_point.Update(1.0, proj_x_n, dist, normal_avg);

  return;
}



void XFEM::XFLUID_STD::callgetNormalSide_tn(DRT::Element* side,  ///< pointer to side element
    LINALG::Matrix<3, 1>& normal,                                ///< normal vector w.r.t side
    Epetra_SerialDenseMatrix& side_xyze,                         ///< side's node coordinates
    const std::vector<int>& lm,                                  ///< local map
    LINALG::Matrix<3, 1>& proj_x_n,                              ///< projected point on side
    LINALG::Matrix<2, 1>& xi_side  ///< local coordinates of projected point w.r.t side
)
{
  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3)  // three dofs per node, for standard Dirichlet coupling
  {
    switch (side->Shape())
    {
        //      case DRT::Element::tri3:
        //      {
        //        getNormalSide_tn<DRT::Element::tri3, 3>(normal, side_xyze, lm, proj_x_n, xi_side);
        //        break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        getNormalSide_tn<DRT::Element::tri6, 3>(normal, side_xyze, lm, proj_x_n, xi_side);
        //        break;
        //      }
      case DRT::Element::quad4:
      {
        getNormalSide_tn<DRT::Element::quad4, 3>(normal, side_xyze, lm, proj_x_n, xi_side);
        break;
      }
      case DRT::Element::quad8:
      {
        getNormalSide_tn<DRT::Element::quad8, 3>(normal, side_xyze, lm, proj_x_n, xi_side);
        break;
      }
        //      case DRT::Element::quad9:
        //      {
        //        getNormalSide_tn<DRT::Element::quad9, 3>(normal, side_xyze, lm, proj_x_n,
        //        xi_side); break;
        //      }
      default:
        dserror("unsupported side shape %d", side->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)  // four dofs per node, for standard Dirichlet coupling
  {
    switch (side->Shape())
    {
        //      case DRT::Element::tri3:
        //      {
        //        getNormalSide_tn<DRT::Element::tri3, 4>(normal, side_xyze, lm, proj_x_n, xi_side);
        //        break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        getNormalSide_tn<DRT::Element::tri6, 4>(normal, side_xyze, lm, proj_x_n, xi_side);
        //        break;
        //      }
      case DRT::Element::quad4:
      {
        getNormalSide_tn<DRT::Element::quad4, 4>(normal, side_xyze, lm, proj_x_n, xi_side);
        break;
      }
      case DRT::Element::quad8:
      {
        getNormalSide_tn<DRT::Element::quad8, 4>(normal, side_xyze, lm, proj_x_n, xi_side);
        break;
      }
        //      case DRT::Element::quad9:
        //      {
        //        getNormalSide_tn<DRT::Element::quad9, 4>(normal, side_xyze, lm, proj_x_n,
        //        xi_side); break;
        //      }
      default:
        dserror("unsupported side shape %d", side->Shape());
        break;
    }
  }
}

template <DRT::Element::DiscretizationType side_distype, const int numdof>
void XFEM::XFLUID_STD::getNormalSide_tn(LINALG::Matrix<3, 1>& normal,  ///< normal vector w.r.t side
    Epetra_SerialDenseMatrix& side_xyze,                               ///< side's node coordinates
    const std::vector<int>& lm,                                        ///< local map
    LINALG::Matrix<3, 1>& proj_x_n,                                    ///< projected point on side
    LINALG::Matrix<2, 1>& xi_side  ///< local coordinates of projected point w.r.t side
)
{
  const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;

  // add displacements
  addeidisp<side_distype, numdof>(side_xyze, *boundarydis_, "idispn", lm);

  // fill Linalg matrix with current node coordinates including displacements
  LINALG::Matrix<3, side_nen_> xyze_(side_xyze);


  // Initialization
  LINALG::Matrix<side_nen_, 1> funct(true);  // shape functions
  LINALG::Matrix<2, side_nen_> deriv(true);  // derivatives dr, ds


  LINALG::Matrix<3, 1> x(true);

  LINALG::Matrix<3, 2> derxy(true);
  LINALG::Matrix<3, 1> dx_dr(true);
  LINALG::Matrix<3, 1> dx_ds(true);

  // get current values
  DRT::UTILS::shape_function_2D(funct, xi_side(0), xi_side(1), side_distype);
  DRT::UTILS::shape_function_2D_deriv1(deriv, xi_side(0), xi_side(1), side_distype);

  proj_x_n.Multiply(xyze_, funct);

  derxy.MultiplyNT(xyze_, deriv);


  // set dx_dr and dx_ds
  for (int i = 0; i < 3; i++)
  {
    dx_dr(i) = derxy(i, 0);
    dx_ds(i) = derxy(i, 1);
  }


  normal(0) = dx_dr(1) * dx_ds(2) - dx_ds(1) * dx_dr(2);
  normal(1) = dx_dr(2) * dx_ds(0) - dx_ds(2) * dx_dr(0);
  normal(2) = dx_dr(0) * dx_ds(1) - dx_ds(0) * dx_dr(1);

  normal.Scale(1.0 / normal.Norm2());

  return;
}



/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to side                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::call_get_projxn_Line(
    DRT::Element* side,  ///< pointer to structural side element
    DRT::Element* line,  ///< pointer to structural line of side element
    LINALG::Matrix<3, 1>& proj_x_n, double& xi_line)
{
  //-------------------------------------------------

  // line geometry at initial state t^0

  // line geometry
  const int numnodes = line->NumNode();
  DRT::Node** nodes = line->Nodes();
  Epetra_SerialDenseMatrix line_xyze(3, numnodes);
  for (int i = 0; i < numnodes; ++i)
  {
    const double* x = nodes[i]->X();
    std::copy(x, x + 3, &line_xyze(0, i));
  }


  DRT::Element::LocationArray cutla(1);
  line->LocationVector(*boundarydis_, cutla, false);

  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3)
  {
    switch (line->Shape())
    {
      case DRT::Element::line2:
      {
        get_projxn_Line<DRT::Element::line2, 3>(line_xyze, cutla[0].lm_, proj_x_n, xi_line);
        break;
      }
      default:
        dserror("unsupported line shape %d", line->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)  // four dofs per node, for standard Dirichlet coupling
  {
    switch (line->Shape())
    {
      case DRT::Element::line2:
      {
        get_projxn_Line<DRT::Element::line2, 4>(line_xyze, cutla[0].lm_, proj_x_n, xi_line);
        break;
      }
      default:
        dserror("unsupported line shape %d", line->Shape());
        break;
    }
  }
}


template <DRT::Element::DiscretizationType line_distype, const int numdof>
void XFEM::XFLUID_STD::get_projxn_Line(
    Epetra_SerialDenseMatrix& line_xyze,  ///< line's node coordinates
    const std::vector<int>& lm,           ///< local map
    LINALG::Matrix<3, 1>& proj_x_n,       ///< projected point on side
    double& xi_line                       ///< local coordinates of projected point w.r.t line
)
{
  const int line_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<line_distype>::numNodePerElement;

  // add displacements
  addeidisp<line_distype, numdof>(line_xyze, *boundarydis_, "idispn", lm);

  // fill LINALG matrix
  LINALG::Matrix<3, line_nen_> xyze_(line_xyze);


  // Initialization
  LINALG::Matrix<line_nen_, 1> funct(true);  // shape functions
  LINALG::Matrix<1, line_nen_> deriv(true);  // derivatives dr

  LINALG::Matrix<3, 1> x(true);


  // get current values
  DRT::UTILS::shape_function_1D(funct, xi_line, line_distype);

  // projected point tracked back at t^n
  proj_x_n.Multiply(xyze_, funct);


  return;
}


/*--------------------------------------------------------------------------------
 * add side's or line's interface displacements and set current node coordinates
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, const int numdof>
void XFEM::XFLUID_STD::addeidisp(
    Epetra_SerialDenseMatrix& xyze,     ///< node coordinates of side or line
    const DRT::Discretization& cutdis,  ///< cut discretization
    const std::string state,            ///< state
    const std::vector<int>& lm          ///< local map
)
{
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  LINALG::Matrix<3, nen> eidisp(true);


  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = cutdis.GetState(state);
  if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

  for (int inode = 0; inode < nen; ++inode)  // number of nodes
  {
    for (int idim = 0; idim < 3; ++idim)  // number of dimensions
    {
      (eidisp)(idim, inode) =
          mymatrix[idim + (inode * numdof)];  // attention! disp state vector has 3+1 dofs for
                                              // displacement (the same as for (u,p))
    }
  }

  // add the displacement of the interface
  for (int inode = 0; inode < nen; ++inode)
  {
    xyze(0, inode) += eidisp(0, inode);
    xyze(1, inode) += eidisp(1, inode);
    xyze(2, inode) += eidisp(2, inode);
  }

  return;
}  // addeidisp


/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to side                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::CallProjectOnSide(
    DRT::Element* side,                   ///< pointer to structural side element
    const std::string state,              ///< state n or np?
    LINALG::Matrix<3, 1>& newNodeCoords,  ///< node coordinates of point that has to be projected
    int& proj_sid,                        ///< id of side when projection lies on side
    double& min_dist,                     ///< minimal distance, potentially updated
    LINALG::Matrix<3, 1>& proj_x_np,      ///< projection of point on this side
    LINALG::Matrix<2, 1>&
        proj_xi_side,              ///< local coordinates of projection of point w.r.t to this side
    TimeIntData::projection& proj  ///< reference to data
)
{
  bool on_side = false;                ///< lies projection on side?
  double curr_dist = INFINITY;         ///< resulting distance
  LINALG::Matrix<3, 1> x_side(true);   ///< resulting projected point
  LINALG::Matrix<2, 1> xi_side(true);  ///< local coordinates resulting projected point


  // side geometry at initial state t^0
  const int numnodes = side->NumNode();
  DRT::Node** nodes = side->Nodes();
  Epetra_SerialDenseMatrix side_xyze(3, numnodes);
  for (int i = 0; i < numnodes; ++i)
  {
    const double* x = nodes[i]->X();
    std::copy(x, x + 3, &side_xyze(0, i));
  }


  DRT::Element::LocationArray cutla(1);
  side->LocationVector(*boundarydis_, cutla, false);

  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3)  // three dofs per node, for standard Dirichlet coupling
  {
    switch (side->Shape())
    {
        //      case DRT::Element::tri3:
        //      {
        //        on_side = ProjectOnSide<DRT::Element::tri3,3>(side_xyze,
        //        cutla[0].lm_,state,newNodeCoords,x_side,xi_side); break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        on_side = ProjectOnSide<DRT::Element::tri6,3>(side_xyze,
        //        cutla[0].lm_,state,newNodeCoords,x_side,xi_side); break;
        //      }
      case DRT::Element::quad4:
      {
        on_side = ProjectOnSide<DRT::Element::quad4, 3>(
            side_xyze, cutla[0].lm_, state, newNodeCoords, x_side, xi_side, curr_dist);
        break;
      }
      case DRT::Element::quad8:
      {
        on_side = ProjectOnSide<DRT::Element::quad8, 3>(
            side_xyze, cutla[0].lm_, state, newNodeCoords, x_side, xi_side, curr_dist);
        break;
      }
      case DRT::Element::quad9:
      {
        on_side = ProjectOnSide<DRT::Element::quad9, 3>(
            side_xyze, cutla[0].lm_, state, newNodeCoords, x_side, xi_side, curr_dist);
        break;
      }
      default:
        dserror("unsupported side shape %d", side->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)  // four dofs per node, for standard Dirichlet coupling
  {
    switch (side->Shape())
    {
        //      case DRT::Element::tri3:
        //      {
        //        on_side = ProjectOnSide<DRT::Element::tri3,4>(side_xyze,
        //        cutla[0].lm_,state,newNodeCoords,x_side,xi_side, curr_dist); break;
        //      }
        //      case DRT::Element::tri6:
        //      {
        //        on_side = ProjectOnSide<DRT::Element::tri6,4>(side_xyze,
        //        cutla[0].lm_,state,newNodeCoords,x_side,xi_side, curr_dist); break;
        //      }
      case DRT::Element::quad4:
      {
        on_side = ProjectOnSide<DRT::Element::quad4, 4>(
            side_xyze, cutla[0].lm_, state, newNodeCoords, x_side, xi_side, curr_dist);
        break;
      }
      case DRT::Element::quad8:
      {
        on_side = ProjectOnSide<DRT::Element::quad8, 4>(
            side_xyze, cutla[0].lm_, state, newNodeCoords, x_side, xi_side, curr_dist);
        break;
      }
        //      case DRT::Element::quad9:
        //      {
        //        //      on_side = ProjectOnSide<DRT::Element::quad9,4>(side_xyze,
        //        cutla[0].lm_,state,newNodeCoords,x_side,xi_side, curr_dist); break;
        //      }
      default:
        dserror("unsupported side shape %d", side->Shape());
        break;
    }
  }

  // update minimal distance if possible
  if (on_side)
  {
#ifdef DEBUG_TIMINT_STD
    std::cout << "\t\tprojection of current node lies on side " << side->Id() << " with distance "
              << curr_dist << std::endl;
#endif

    if (curr_dist < min_dist - TOL_dist_ && curr_dist >= 0)
    {
#ifdef DEBUG_TIMINT_STD
      std::cout << "\t\t\t>>> updated smallest distance! " << curr_dist << std::endl;
#endif


      //--------------------
      // set current minimal distance w.r.t side
      proj = TimeIntData::onSide_;
      // set side that contains the projected point
      proj_sid = side->Id();
      // update minimal distance
      min_dist = curr_dist;
      // update projection point
      proj_x_np.Update(1.0, x_side, 0.0);
      // update local coordinates w.r.t side
      proj_xi_side.Update(1.0, xi_side, 0.0);
      //--------------------
    }
    else if (curr_dist < 0)
    {
#ifdef DEBUG_TIMINT_STD
      std::cout << "\t\t\tnegative distance -> do not update" << std::endl;
#endif
    }
  }

  return;
}



/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to side                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::CallProjectOnLine(
    DRT::Element* side,                   ///< pointer to structural side element
    DRT::Element* line,                   ///< pointer to structural line of side element
    int line_count,                       ///< local line id w.r.t side element
    const std::string state,              ///< state n or np?
    LINALG::Matrix<3, 1>& newNodeCoords,  ///< node coordinates of point that has to be projected
    double& min_dist,                     ///< minimal distance, potentially updated
    LINALG::Matrix<3, 1>& proj_x_np,      ///< projection of point on this side
    std::map<std::vector<int>, std::vector<double>>&
        proj_xi_line,  ///< std::map<sorted nids, local line coordinates of projection of point
                       ///< w.r.t sides >
    std::map<std::vector<int>, std::vector<int>>&
        proj_lineid,  ///< std::map<sorted nids, local line id w.r.t sides>
    std::map<std::vector<int>, std::vector<int>>&
        proj_nid_line,             ///< std::map<sorted nids, side Ids>
    int& proj_sid,                 ///< id of side that contains the projected point
    TimeIntData::projection& proj  ///< reference to data
)
{
  bool on_line = false;               ///< is projection on line
  double curr_dist = INFINITY;        ///< resulting distance
  LINALG::Matrix<3, 1> x_line(true);  ///< resulting projected point
  double xi_line = INFINITY;          ///< local coordinates resulting projected point

  //-------------------------------------------------

  // line geometry at initial state t^0

  // line geometry
  const int numnodes = line->NumNode();
  DRT::Node** nodes = line->Nodes();
  Epetra_SerialDenseMatrix line_xyze(3, numnodes);
  for (int i = 0; i < numnodes; ++i)
  {
    const double* x = nodes[i]->X();
    std::copy(x, x + 3, &line_xyze(0, i));
  }


  DRT::Element::LocationArray cutla(1);
  line->LocationVector(*boundarydis_, cutla, false);

  // get number of dofs for this side element
  const int numdofpernode = side->NumDofPerNode(*side->Nodes()[0]);

  if (numdofpernode == 3)
  {
    switch (line->Shape())
    {
      case DRT::Element::line2:
      {
        on_line = ProjectOnLine<DRT::Element::line2, 3>(
            line_xyze, cutla[0].lm_, state, newNodeCoords, x_line, xi_line, curr_dist);
        break;
      }
      default:
        dserror("unsupported line shape %d", line->Shape());
        break;
    }
  }
  else if (numdofpernode == 4)  // four dofs per node, for standard Dirichlet coupling
  {
    switch (line->Shape())
    {
      case DRT::Element::line2:
      {
        on_line = ProjectOnLine<DRT::Element::line2, 4>(
            line_xyze, cutla[0].lm_, state, newNodeCoords, x_line, xi_line, curr_dist);
        break;
      }
      case DRT::Element::line3:
      {
        on_line = ProjectOnLine<DRT::Element::line3, 4>(
            line_xyze, cutla[0].lm_, state, newNodeCoords, x_line, xi_line, curr_dist);
        break;
      }
      default:
        dserror("unsupported line shape %d", line->Shape());
        break;
    }
  }

  // update minimal distance if possible
  if (on_line)
  {
#ifdef DEBUG_TIMINT_STD
    std::cout << "\t\tline-proj: projection lies on line " << line_count << " of side "
              << side->Id() << std::endl;
#endif

    // distance has to be smaller than the last one
    if (curr_dist < (min_dist - TOL_dist_) && curr_dist >= 0)
    {
#ifdef DEBUG_TIMINT_STD
      std::cout.precision(15);
      std::cout << "\t\t\t>>> updated smallest distance!" << curr_dist << std::endl;
#endif

      if (curr_dist < (min_dist - TOL_dist_))  // smaller distance found
      {
        // another line found, reset already found lines
        proj_lineid.clear();
        proj_xi_line.clear();
        proj_nid_line.clear();
      }
      else
      {
        // add lines, that have the same distance up to TOL_dist
      }

      std::vector<int> line_nids;
      for (int i = 0; i < numnodes; ++i)
      {
        line_nids.push_back(nodes[i]->Id());
      }

      sort(line_nids.begin(), line_nids.end());

      // reset side id for projection on side
      proj_sid = -1;

      //--------------------
      // set current minimal distance w.r.t line
      proj = TimeIntData::onLine_;
      // update minimal distance
      min_dist = curr_dist;
      // update projection point
      proj_x_np.Update(1.0, x_line, 0.0);



      std::map<std::vector<int>, std::vector<int>>::iterator lines = proj_nid_line.find(line_nids);

      // line already inserted via other side
      if (lines != proj_nid_line.end())
      {
        // set line id w.r.t side->Id() that contains the projected point
        (lines->second).push_back(side->Id());
        // update local line id w.r.t side
        (proj_lineid.find(line_nids))->second.push_back(line_count);
        // update local coordinates w.r.t line
        (proj_xi_line.find(line_nids))->second.push_back(xi_line);
      }
      else
      {
        std::vector<int> sideids;
        sideids.push_back(side->Id());
        std::vector<int> locallineids;
        locallineids.push_back(line_count);
        std::vector<double> locallineXiCoords;
        locallineXiCoords.push_back(xi_line);

        // set line id w.r.t side->Id() that contains the projected point
        proj_nid_line.insert(std::pair<std::vector<int>, std::vector<int>>(line_nids, sideids));
        // update local line id w.r.t side
        proj_lineid.insert(std::pair<std::vector<int>, std::vector<int>>(line_nids, locallineids));
        // update local coordinates w.r.t line
        proj_xi_line.insert(
            std::pair<std::vector<int>, std::vector<double>>(line_nids, locallineXiCoords));
      }

      //--------------------
    }
    else if (curr_dist < 0)
    {
      dserror(" no negative distance to line possible");
    }
  }  // on_line

  return;
}

/*------------------------------------------------------------------------------------------------*
 * call and prepare the projection of point to point (distance computation)          schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::CallProjectOnPoint(DRT::Node* node,  ///< pointer to node
    const std::string state,                                ///< state n or np?
    LINALG::Matrix<3, 1>& newNodeCoords,  ///< node coordinates of point that has to be projected
    double& min_dist,                     ///< minimal distance, potentially updated
    int& proj_nid_np,                     ///< nid id of projected point on surface
    LINALG::Matrix<3, 1>& proj_x_np,      ///< projection of point on this side
    int& proj_sid,                        ///< id of side that contains the projected point
    std::map<std::vector<int>, std::vector<double>>
        proj_xi_line,  ///< std::map<side ID, local coordinates of projection of point w.r.t to this
                       ///< line>
    std::map<std::vector<int>, std::vector<int>> proj_lineid,  ///< std::map<side ID, local line id>
    std::map<std::vector<int>, std::vector<int>>
        proj_nid_line,             ///< std::map<side ID, vec<line Ids>>
    TimeIntData::projection& proj  ///< projection type
)
{
  double curr_dist = INFINITY;         ///< resulting distance
  LINALG::Matrix<3, 1> x_point(true);  ///< resulting projected point


  // its point geometry
  Epetra_SerialDenseMatrix point_xyze(3, 1);

  const double* x = node->X();
  std::copy(x, x + 3, &point_xyze(0, 0));

  std::vector<int> lm = boundarydis_->Dof(0, node);


  // compute distance between two points
  ProjectOnPoint(point_xyze, lm, state, newNodeCoords, x_point, curr_dist);


  // update minimal distance if possible
  if (curr_dist < min_dist - TOL_dist_ and curr_dist >= 0)
  {
#ifdef DEBUG_TIMINT_STD
    std::cout << "\t\t\t>>> updated smallest distance w.r.t point-proj!" << curr_dist << std::endl;
#endif
    //--------------------
    // set current minimal distance w.r.t point
    proj = TimeIntData::onPoint_;
    // reset side id
    proj_sid = -1;
    // reset line id w.r.t side->Id() that contains the projected point
    proj_lineid.clear();
    proj_xi_line.clear();
    proj_nid_line.clear();
    // update minimal distance
    min_dist = curr_dist;
    // update point with minimal distance
    proj_x_np.Update(1.0, x_point, 0.0);
    // update node id of point with minimal distance
    proj_nid_np = node->Id();
    //--------------------
  }

  return;
}


/*--------------------------------------------------------------------------------*
 * project point from in normal direction onto corresponding side    schott 07/12 *
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType side_distype, const int numdof>
bool XFEM::XFLUID_STD::ProjectOnSide(
    Epetra_SerialDenseMatrix& side_xyze,  ///< side's node coordinates
    const std::vector<int>& lm,           ///< local map
    const std::string state,              ///< state n or np
    LINALG::Matrix<3, 1>& x_gp_lin,       ///< global coordinates of point that has to be projected
    LINALG::Matrix<3, 1>& x_side,         ///< projected point on side
    LINALG::Matrix<2, 1>& xi_side,        ///< local coordinates of projected point w.r.t side
    double& dist                          ///< distance from point to its projection
)
{
  const int side_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<side_distype>::numNodePerElement;

  // add displacements
  addeidisp<side_distype, numdof>(side_xyze, *boundarydis_, state, lm);

  // fill Linalg matrix with current node coordinates including displacements
  LINALG::Matrix<3, side_nen_> xyze_(side_xyze);


  // Initialization
  LINALG::Matrix<side_nen_, 1> funct(true);   // shape functions
  LINALG::Matrix<2, side_nen_> deriv(true);   // derivatives dr, ds
  LINALG::Matrix<3, side_nen_> deriv2(true);  // 2nd derivatives drdr, dsds, drds


  LINALG::Matrix<3, 1> x(true);

  LINALG::Matrix<3, 2> derxy(true);
  LINALG::Matrix<3, 1> dx_dr(true);
  LINALG::Matrix<3, 1> dx_ds(true);

  LINALG::Matrix<3, 3> derxy2(true);
  LINALG::Matrix<3, 1> dx_drdr(true);
  LINALG::Matrix<3, 1> dx_dsds(true);
  LINALG::Matrix<3, 1> dx_drds(true);

  LINALG::Matrix<3, 1> dx_drdr_times_dx_ds(true);
  LINALG::Matrix<3, 1> dx_dr_times_dx_drds(true);
  LINALG::Matrix<3, 1> dx_drds_times_dx_ds(true);
  LINALG::Matrix<3, 1> dx_dr_times_dx_dsds(true);
  LINALG::Matrix<3, 1> dx_dr_times_dx_ds(true);

  LINALG::Matrix<3, 1> residuum(true);  // residuum of the newton iteration
  LINALG::Matrix<3, 3> sysmat(true);    // matrix for the newton system
  LINALG::Matrix<3, 1> incr(true);      // increment of the newton system

  LINALG::Matrix<3, 1> sol(true);  // sol carries xi_1, xi_2, d (distance)

  if (side_distype == DRT::Element::tri3 or side_distype == DRT::Element::tri6)
  {
    sol(0) = 0.333333333333333;
    sol(1) = 0.333333333333333;
  }
  else if (side_distype == DRT::Element::quad4 or side_distype == DRT::Element::quad8 or
           side_distype == DRT::Element::quad9)
  {
    sol(0) = 0.0;
    sol(1) = 0.0;
  }
  else
  {
    dserror("define start side xi-coordinates for unsupported cell type");
  }

  const double absTolIncr = 1.0e-9;  // absolute tolerance for the local coordinates increment
  const double absTolRes = 1.0e-9;   // absolute tolerance for the whole residual
  const double absTOLdist = 1.0e-9;  // absolute tolerance for distance

  int iter = 0;
  const int maxiter = 7;

  bool converged = false;

  while (iter < maxiter && !converged)
  {
    iter++;


    // get current values
    DRT::UTILS::shape_function_2D(funct, sol(0), sol(1), side_distype);
    DRT::UTILS::shape_function_2D_deriv1(deriv, sol(0), sol(1), side_distype);
    DRT::UTILS::shape_function_2D_deriv2(deriv2, sol(0), sol(1), side_distype);

    x.Multiply(xyze_, funct);

    derxy.MultiplyNT(xyze_, deriv);

    derxy2.MultiplyNT(xyze_, deriv2);

    // set dx_dr and dx_ds
    for (int i = 0; i < 3; i++)
    {
      dx_dr(i) = derxy(i, 0);
      dx_ds(i) = derxy(i, 1);

      dx_drdr(i) = derxy2(i, 0);
      dx_dsds(i) = derxy2(i, 1);
      dx_drds(i) = derxy2(i, 2);
    }

    // get vector products

    dx_drdr_times_dx_ds(0) = dx_drdr(1) * dx_ds(2) - dx_ds(1) * dx_drdr(2);
    dx_drdr_times_dx_ds(1) = dx_drdr(2) * dx_ds(0) - dx_ds(2) * dx_drdr(0);
    dx_drdr_times_dx_ds(2) = dx_drdr(0) * dx_ds(1) - dx_ds(0) * dx_drdr(1);

    dx_dr_times_dx_drds(0) = dx_dr(1) * dx_drds(2) - dx_drds(1) * dx_dr(2);
    dx_dr_times_dx_drds(1) = dx_dr(2) * dx_drds(0) - dx_drds(2) * dx_dr(0);
    dx_dr_times_dx_drds(2) = dx_dr(0) * dx_drds(1) - dx_drds(0) * dx_dr(1);

    dx_drds_times_dx_ds(0) = dx_drds(1) * dx_ds(2) - dx_ds(1) * dx_drds(2);
    dx_drds_times_dx_ds(1) = dx_drds(2) * dx_ds(0) - dx_ds(2) * dx_drds(0);
    dx_drds_times_dx_ds(2) = dx_drds(0) * dx_ds(1) - dx_ds(0) * dx_drds(1);

    dx_dr_times_dx_dsds(0) = dx_dr(1) * dx_dsds(2) - dx_dsds(1) * dx_dr(2);
    dx_dr_times_dx_dsds(1) = dx_dr(2) * dx_dsds(0) - dx_dsds(2) * dx_dr(0);
    dx_dr_times_dx_dsds(2) = dx_dr(0) * dx_dsds(1) - dx_dsds(0) * dx_dr(1);

    dx_dr_times_dx_ds(0) = dx_dr(1) * dx_ds(2) - dx_ds(1) * dx_dr(2);
    dx_dr_times_dx_ds(1) = dx_dr(2) * dx_ds(0) - dx_ds(2) * dx_dr(0);
    dx_dr_times_dx_ds(2) = dx_dr(0) * dx_ds(1) - dx_ds(0) * dx_dr(1);

    // define sysmat
    for (int i = 0; i < 3; i++)
    {
      // d/dr
      sysmat(i, 0) = dx_dr(i) - sol(2) * (dx_drdr_times_dx_ds(i) + dx_dr_times_dx_drds(i));

      // d/ds
      sysmat(i, 1) = dx_ds(i) - sol(2) * (dx_drds_times_dx_ds(i) + dx_dr_times_dx_dsds(i));

      // d/d(dist)
      sysmat(i, 2) = -dx_dr_times_dx_ds(i);


      // residual
      residuum(i) = x(i) - sol(2) * dx_dr_times_dx_ds(i) - x_gp_lin(i);
    }



    sysmat.Invert();

    // solve Newton iteration
    incr.Clear();
    incr.Multiply(-1.0, sysmat, residuum);  // incr = -Systemmatrix^-1 * residuum

    // update solution
    sol.Update(1.0, incr, 1.0);

    if ((incr.Norm2() < absTolIncr) && (residuum.Norm2() < absTolRes))
    {
      converged = true;
    }

    // check ° relative criterion for local coordinates (between [-1,1]^2)
    //       ° absolute criterion for distance (-> 0)
    //       ° relative criterion for whole residuum
    if (  // sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1)) <  relTolIncr
        sqrt(incr(0) * incr(0) + incr(1) * incr(1)) < absTolIncr && incr(2) < absTOLdist &&
        residuum.Norm2() < absTolRes)
    {
      converged = true;
    }
  }

  bool on_side = false;

  if (!converged)
  {
#ifdef DEBUG_TIMINT_STD
    std::cout.precision(15);

    std::cout << "increment criterion loc coord "
              //<< sqrt(incr(0)*incr(0)+incr(1)*incr(1))/sqrt(sol(0)*sol(0)+sol(1)*sol(1))
              << sqrt(incr(0) * incr(0) + incr(1) * incr(1)) << " \tabsTOL: " << absTolIncr
              << std::endl;
    std::cout << "absolute criterion for distance " << incr(2) << " \tabsTOL: " << absTOLdist
              << std::endl;
    std::cout << "absolute criterion whole residuum " << residuum.Norm2()
              << " \tabsTOL: " << absTolRes << std::endl;


    std::cout << "sysmat.Invert" << sysmat << std::endl;
    std::cout << "sol-norm " << sol.Norm2() << std::endl;
    std::cout << "sol " << sol << std::endl;
    std::cout << "x_gp_lin" << x_gp_lin << std::endl;
    std::cout << "side " << xyze_ << std::endl;
#endif

    // dserror( "newton scheme in ProjectOnSide not converged! " );

    xi_side(0) = INFINITY;
    xi_side(1) = INFINITY;

    dist = INFINITY;

    on_side = false;
  }
  else
  {
    LINALG::Matrix<3, 1> xsi(true);
    xsi(0) = sol(0);
    xsi(1) = sol(1);
    xsi(2) = 0;

    bool within_limits = WithinLimits<side_distype>(xsi, 1e-10);

    if (within_limits)
    {
      on_side = true;


      // get length of current normal vector
      double normal_length = dx_dr_times_dx_ds.Norm2();


      dist = -sol(2) * normal_length;  // negative sol(2)!!! and scaling with normal length

      // evaluate shape function at solution
      DRT::UTILS::shape_function_2D(funct, sol(0), sol(1), side_distype);

      // get projected gauss point
      x_side.Multiply(xyze_, funct);

      // set local coordinates w.r.t side
      xi_side(0) = sol(0);
      xi_side(1) = sol(1);

#ifdef DEBUG_TIMINT_STD
      std::cout << "\t\tside-proj: converged with local coordinates " << sol(0) << " " << sol(1)
                << " "
                << " dist: " << sol(2) * normal_length << std::endl;
#endif
    }
    else
    {
      on_side = false;

      xi_side(0) = INFINITY;
      xi_side(1) = INFINITY;

      dist = INFINITY;

      x_side(0) = INFINITY;
      x_side(1) = INFINITY;
      x_side(2) = INFINITY;

#ifdef DEBUG_TIMINT_STD
      std::cout << "\t\tside-proj: converged with local coordinates " << sol(0) << " " << sol(1)
                << ". Check if is on side: false" << std::endl;
#endif
    }
  }

  return on_side;
}  // ProjectOnSide

/*--------------------------------------------------------------------------------*
 * project point in normal direction onto corresponding line         schott 07/12 *
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType line_distype, const int numdof>
bool XFEM::XFLUID_STD::ProjectOnLine(
    Epetra_SerialDenseMatrix& line_xyze,  ///< line's node coordinates
    const std::vector<int>& lm,           ///< local map
    const std::string state,              ///< state n or np?
    LINALG::Matrix<3, 1>& x_point_np,     ///< global coordinates of point that has to be projected
    LINALG::Matrix<3, 1>& x_line,         ///< projected point on line
    double& xi_line,                      ///< local coordinates of projected point w.r.t line
    double& dist                          ///< distance from point to its projection
)
{
  bool on_line = false;

  const int line_nen_ = DRT::UTILS::DisTypeToNumNodePerEle<line_distype>::numNodePerElement;

  // add displacements
  addeidisp<line_distype, numdof>(line_xyze, *boundarydis_, state, lm);

  // fill LINALG matrix
  LINALG::Matrix<3, line_nen_> xyze_(line_xyze);

  // linear line
  if (line_nen_ == 2)
  {
    // get line direction
    LINALG::Matrix<3, 1> line_dir(true);  // x2-x1
    LINALG::Matrix<3, 1> dir_tmp(true);   // x - 0.5(x1 + x2)
    for (int isd = 0; isd < 3; isd++)
    {
      line_dir(isd) = xyze_(isd, 1) - xyze_(isd, 0);
      dir_tmp(isd) = x_point_np(isd, 0) - 0.5 * (xyze_(isd, 1) + xyze_(isd, 0));
    }

    double line_length = line_dir.Norm2();  // || (x2-x1) ||

    if (line_length < 1e-12) dserror("line has lenght smaller than 1e-12");

    // 2.0/|| (x2-x1) ||^2 * < x-0.5(x1+x2),x2-x1 >
    xi_line = 2.0 / (line_length * line_length) * line_dir.Dot(dir_tmp);

    // check if projection within line segement between x1 and x2
    LINALG::Matrix<3, 1> xsi(true);
    xsi(0) = xi_line;
    xsi(1) = 0;
    xsi(2) = 0;

    bool within_limits = WithinLimits<line_distype>(xsi, 1e-10);

    if (within_limits)
    {
      on_line = true;

      // set projection
      for (int i = 0; i < 3; i++)
      {
        x_line(i) = 0.5 * (1 + xi_line) * xyze_(i, 1) + 0.5 * (1 - xi_line) * xyze_(i, 0);
      }

      // set distance vector between x and projection
      LINALG::Matrix<3, 1> dist_vec(true);
      for (int i = 0; i < 3; i++)
      {
        dist_vec(i) = x_point_np(i) - x_line(i);
      }

      // compute distance
      dist = dist_vec.Norm2();
    }
    else
    {
      on_line = false;
      xi_line = INFINITY;
      dist = INFINITY;

      x_line(0) = INFINITY;
      x_line(1) = INFINITY;
      x_line(2) = INFINITY;
    }
  }
  else
  {
    dserror("projection on line just for linear lines implemented!");
  }

  return on_line;
}  // ProjectOnLine

/*--------------------------------------------------------------------------------*
 * compute distance (project) between two points                     schott 07/12 *
 *--------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::ProjectOnPoint(
    Epetra_SerialDenseMatrix& point_xyze,  ///< point's node coordinates
    const std::vector<int>& lm,            ///< local map
    const std::string state,               ///< state n or np?
    LINALG::Matrix<3, 1>& x_point_np,      ///< global coordinates of point that has to be projected
    LINALG::Matrix<3, 1>& x_point,         ///< global coordinates of point on surface
    double& dist                           ///< distance from point to its projection
)
{
  Teuchos::RCP<const Epetra_Vector> matrix_state = boundarydis_->GetState(state);
  if (matrix_state == Teuchos::null) dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state, mymatrix, lm);

  // add the displacement of the interface
  for (int idim = 0; idim < 3; ++idim)  // number of dimensions
  {
    (point_xyze)(idim, 0) += mymatrix[idim];  // attention! disp state vector has 3+1 dofs for
                                              // displacement (the same as for (u,p))
  }

  // fill LINALG matrix
  LINALG::Matrix<3, 1> p(point_xyze);

  // compute direction vector between two points
  LINALG::Matrix<3, 1> direction(true);
  direction.Update(1.0, x_point_np, -1.0, p);

  // compute distance
  dist = direction.Norm2();

  x_point.Update(1.0, p, 0.0);

  return;
}  // ProjectOnPoint


/*--------------------------------------------------------------------------------*
 * check if local coordinates are within limits                      schott 07/12 *
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType elementtype>
bool XFEM::XFLUID_STD::WithinLimits(LINALG::Matrix<3, 1>& xsi_, const double TOL)
{
  switch (elementtype)
  {
    case DRT::Element::line2:
    case DRT::Element::line3:
      return (xsi_(0) >= -1 - TOL and xsi_(0) <= 1 + TOL);
    case DRT::Element::tri3:
    case DRT::Element::tri6:
      return (xsi_(0) >= -TOL and xsi_(1) >= -TOL and xsi_(0) <= 1 + TOL and xsi_(1) <= 1 + TOL and
              (xsi_(0) + xsi_(1)) <= 1 + TOL);
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      return (xsi_(0) >= -1 - TOL and xsi_(1) >= -1 - TOL and xsi_(0) <= 1 + TOL and
              xsi_(1) <= 1 + TOL);
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
      return (xsi_(0) >= -1 - TOL and xsi_(1) >= -1 - TOL and xsi_(2) >= -1 - TOL and
              xsi_(0) <= 1 + TOL and xsi_(1) <= 1 + TOL and xsi_(2) <= 1 + TOL);
    case DRT::Element::tet4:
    case DRT::Element::tet10:
      return (xsi_(0) >= -TOL and xsi_(1) >= -TOL and xsi_(2) >= -TOL and
              xsi_(0) + xsi_(1) + xsi_(2) <= 1 + TOL);
    case DRT::Element::pyramid5:
      return (xsi_(0) >= -1 - TOL and xsi_(1) >= -1 - TOL and xsi_(2) >= -TOL and
              xsi_(0) <= 1 + TOL and xsi_(1) <= 1 + TOL and
              (fabs(xsi_(0)) > fabs(xsi_(1)) ? (xsi_(2) + fabs(xsi_(0)) <= 1 + TOL)
                                             : (xsi_(2) + fabs(xsi_(1)) <= 1 + TOL)));
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
      return (xsi_(0) >= -TOL and xsi_(1) >= -TOL and xsi_(2) >= -1 - TOL and xsi_(2) <= 1 + TOL and
              xsi_(0) + xsi_(1) <= 1 + TOL);
    default:
      throw std::runtime_error("unsupported element type in XFEM::XFLUID_STD::WithinLimits");
  }
  return false;
}

/*------------------------------------------------------------------------------------------------*
 * setting the computed data for the standard degrees of freedom into the according               *
 * Epetra Vectors for all handled nodes                                              schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::setFinalData()
{
  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // loop over data
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->state_ != TimeIntData::doneStd_)
      dserror("when data is set, all computation has to be done");

    std::vector<LINALG::Matrix<nsd, 1>>& velValues(data->velValues_);  // velocities of the node
    std::vector<double>& presValues(data->presValues_);                // pressures of the node

    const int gnodeid = data->node_.Id();  // global node id

    //-------------------------------------------------------
    DRT::Node* node = discret_->gNode(gnodeid);

#ifdef DEBUG_TIMINT_STD
    IO::cout << "dofset at new timestep " << data->nds_np_ << IO::endl;
#endif

    if (data->nds_np_ == -1) dserror("cannot get dofs for dofset with number %d", data->nds_np_);

    std::vector<int> dofs;
    dofset_new_->Dof(dofs, node, data->nds_np_);

    if (dofs.size() != (nsd + 1))
      dserror("not the right number of dofs %d for this node ", dofs.size());


    // set velocity dofs
    for (int i = 0; i < nsd; i++)
    {
      for (size_t index = 0; index < vectorSize(data->type_); index++)
        (*newVectors_[index])[newdofrowmap_.LID(dofs[i])] =
            velValues[index](i, 0);  // set the value
    }
    for (size_t index = 0; index < vectorSize(data->type_); index++)
      (*newVectors_[index])[newdofrowmap_.LID(dofs[nsd])] = presValues[index];  // set the value

    data->type_ = TimeIntData::standard_;  // predictor is done, so next time standard
  }                                        // end loop over nodes
}  // end setFinalData



/*------------------------------------------------------------------------------------------------*
 * export start data to neighbour proc                                           winklmaier 06/10 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::exportStartData()
{
  dserror("this function has still to be adapted for XFSI");

  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // destination proc (the "next" one)
  int dest = myrank_ + 1;
  if (myrank_ == (numproc_ - 1)) dest = 0;

  // source proc (the "last" one)
  int source = myrank_ - 1;
  if (myrank_ == 0) source = numproc_ - 1;

  DRT::PackBuffer dataSend;  // data to be sent

  // packing the data
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    packNode(dataSend, data->node_);
    DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
    DRT::ParObject::AddtoPack(dataSend, data->vel_);
    DRT::ParObject::AddtoPack(dataSend, data->velDeriv_);
    DRT::ParObject::AddtoPack(dataSend, data->presDeriv_);
    DRT::ParObject::AddtoPack(dataSend, data->dispnp_);
    DRT::ParObject::AddtoPack(dataSend, data->startpoint_);
    DRT::ParObject::AddtoPack(dataSend, data->searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend, data->counter_);
    DRT::ParObject::AddtoPack(dataSend, data->dMin_);
    DRT::ParObject::AddtoPack(dataSend, (int)data->type_);
  }

  dataSend.StartPacking();

  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    packNode(dataSend, data->node_);
    DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
    DRT::ParObject::AddtoPack(dataSend, data->vel_);
    DRT::ParObject::AddtoPack(dataSend, data->velDeriv_);
    DRT::ParObject::AddtoPack(dataSend, data->presDeriv_);
    DRT::ParObject::AddtoPack(dataSend, data->dispnp_);
    DRT::ParObject::AddtoPack(dataSend, data->startpoint_);
    DRT::ParObject::AddtoPack(dataSend, data->searchedProcs_);
    DRT::ParObject::AddtoPack(dataSend, data->counter_);
    DRT::ParObject::AddtoPack(dataSend, data->dMin_);
    DRT::ParObject::AddtoPack(dataSend, (int)data->type_);
  }

  std::vector<char> dataRecv;
  sendData(dataSend, dest, source, dataRecv);

  // pointer to current position of group of cells in global std::string (counts bytes)
  std::vector<char>::size_type posinData = 0;

  // clear vector that should be filled
  timeIntData_->clear();

  // unpack received data
  while (posinData < dataRecv.size())
  {
    double coords[nsd] = {0.0};
    DRT::Node node(0, (double*)coords, 0);  // initialize node
    int nds_np = -1;
    LINALG::Matrix<nsd, 1> vel;                      // velocity at point x
    std::vector<LINALG::Matrix<nsd, nsd>> velDeriv;  // derivation of velocity at point x
    std::vector<LINALG::Matrix<1, nsd>> presDeriv;   // derivation of pressure at point x
    LINALG::Matrix<nsd, 1> dispnp;                   // displacement at point x
    LINALG::Matrix<nsd, 1> startpoint;               // startpoint
    int searchedProcs;                               // number of searched processors
    int counter;                                     // iteration counter
    double dMin;                                     // minimal distance
    int newtype;                                     // type of the data

    unpackNode(posinData, dataRecv, node);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, nds_np);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, vel);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, velDeriv);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, presDeriv);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, dispnp);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, startpoint);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, searchedProcs);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, counter);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, dMin);
    DRT::ParObject::ExtractfromPack(posinData, dataRecv, newtype);

    timeIntData_->push_back(TimeIntData(node, nds_np, vel, velDeriv, presDeriv, dispnp, startpoint,
        searchedProcs, counter, dMin, (TimeIntData::type)newtype));
  }

  discret_->Comm().Barrier();  // processors wait for each other
}  // end exportStartData



/*------------------------------------------------------------------------------------------------*
 * export final data to node's proc                                                  schott 07/12 *
 *------------------------------------------------------------------------------------------------*/
void XFEM::XFLUID_STD::exportFinalData()
{
#ifdef DEBUG_TIMINT_STD
  IO::cout << "\n\t=============================";
  IO::cout << "\n\t  export Final Data  ";
  IO::cout << "\n\t=============================" << IO::endl;
#endif


  const int nsd = 3;  // 3 dimensions for a 3d fluid element

  // array of vectors which stores data for
  // every processor in one vector
  std::vector<std::vector<TimeIntData>> dataVec(numproc_);

  // fill vectors with the data
  // fill data into a vector where the index corresponds to the destination which is equal to the
  // node's owner
  for (std::vector<TimeIntData>::iterator data = timeIntData_->begin(); data != timeIntData_->end();
       data++)
  {
    if (data->state_ != TimeIntData::doneStd_)
      dserror("All data should be set here, having status 'done'. Thus something is wrong!");
    dataVec[data->node_.Owner()].push_back(*data);
  }

  timeIntData_->clear();
  *timeIntData_ = dataVec[myrank_];  // set final data of own processor
  dataVec[myrank_].clear();          // clear data about current proc

  // send data to the processor where the point lies (1. nearest higher neighbour 2. 2nd nearest
  // higher neighbour...)
  for (int dest = (myrank_ + 1) % numproc_; dest != myrank_;
       dest = (dest + 1) % numproc_)  // dest is the target processor
  {
    // Initialization
    int source = myrank_ - (dest - myrank_);  // source proc (sends (dest-myrank_) far and gets from
                                              // (dest-myrank_) earlier)
    if (source < 0)
      source += numproc_;
    else if (source >= numproc_)
      source -= numproc_;

    DRT::PackBuffer dataSend;

    // pack data to be sent
    for (std::vector<TimeIntData>::iterator data = dataVec[dest].begin();
         data != dataVec[dest].end(); data++)
    {
      DRT::ParObject::AddtoPack(dataSend, data->node_.Id());
      DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
      DRT::ParObject::AddtoPack(dataSend, data->startpoint_);
      DRT::ParObject::AddtoPack(dataSend, data->velValues_);
      DRT::ParObject::AddtoPack(dataSend, data->presValues_);
      DRT::ParObject::AddtoPack(dataSend, data->type_);
    }

    dataSend.StartPacking();

    for (std::vector<TimeIntData>::iterator data = dataVec[dest].begin();
         data != dataVec[dest].end(); data++)
    {
      DRT::ParObject::AddtoPack(dataSend, data->node_.Id());
      DRT::ParObject::AddtoPack(dataSend, data->nds_np_);
      DRT::ParObject::AddtoPack(dataSend, data->startpoint_);
      DRT::ParObject::AddtoPack(dataSend, data->velValues_);
      DRT::ParObject::AddtoPack(dataSend, data->presValues_);
      DRT::ParObject::AddtoPack(dataSend, data->type_);
    }

    // clear the no more needed data
    dataVec[dest].clear();

    std::vector<char> dataRecv;
    sendData(dataSend, dest, source, dataRecv);

    // pointer to current position of group of cells in global std::string (counts bytes)
    std::vector<char>::size_type posinData = 0;

    // unpack received data
    while (posinData < dataRecv.size())
    {
      int gid;                                        // global id of node
      int nds_np;                                     // dofset number of node at new timestep
      LINALG::Matrix<nsd, 1> startpoint;              // startpoint
      std::vector<LINALG::Matrix<nsd, 1>> velValues;  // velocity values
      std::vector<double> presValues;                 // pressure values
      int newtype;                                    // type of the data

      DRT::ParObject::ExtractfromPack(posinData, dataRecv, gid);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, nds_np);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, startpoint);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, velValues);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, presValues);
      DRT::ParObject::ExtractfromPack(posinData, dataRecv, newtype);

      LINALG::Matrix<3, 1> nodedispnp(true);
      if (dispnp_ != Teuchos::null)  // is alefluid
      {
        //------------------------------------------------------- add ale disp
        // get node location vector, dirichlet flags and ownerships (discret, nds, la, doDirichlet)
        std::vector<int> lm;
        std::vector<int> dofs;
        dofset_new_->Dof(dofs, discret_->gNode(gid), 0);  // dofs for standard dofset
        for (int j = 0; j < 4; ++j) lm.push_back(dofs[j]);

        LINALG::Matrix<1, 1> nodepredummy(true);
        extractNodalValuesFromVector<1>(nodedispnp, nodepredummy, dispnp_, lm);
      }

      timeIntData_->push_back(TimeIntData(*discret_->gNode(gid), nds_np, nodedispnp, startpoint,
          velValues, presValues, (TimeIntData::type)newtype));
    }  // end loop over number of nodes to get

    // processors wait for each other
    discret_->Comm().Barrier();
  }  // end loop over processors
}  // end exportfinalData
