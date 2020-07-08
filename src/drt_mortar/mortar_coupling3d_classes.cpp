/*-----------------------------------------------------------------------*/
/*! \file
\brief A class for mortar coupling of ONE slave element and ONE master
       element of a mortar interface in 3D (definition of sub-classes).

\level 1

*/
/*-----------------------------------------------------------------------*/

#include "mortar_coupling3d.H"
#include "mortar_coupling3d_classes.H"
#include "mortar_projector.H"
#include "mortar_integrator.H"
#include "mortar_defines.H"
#include "mortar_element.H"

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_node.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/09|
 *----------------------------------------------------------------------*/
MORTAR::IntElement::IntElement(int lid, int id, int owner, MORTAR::MortarElement* parele,
    const DRT::Element::DiscretizationType& shape, const int numnode, const int* nodeids,
    std::vector<DRT::Node*> nodes, const bool isslave, const bool rewind)
    : MORTAR::MortarElement(id, owner, shape, numnode, nodeids, isslave),
      lid_(lid),
      rewind_(rewind),
      parele_(parele)
{
  if ((int)nodes.size() != numnode) dserror("some inconsistency");

  // check for consistency of nodeids and nodes
  // for nurbs, the nodes are not actual nodes in the
  // discretization, so just skip that part.
  if (ParShape() != DRT::Element::nurbs9)
    for (int i = 0; i < numnode; ++i)
      if (nodes[i]->Id() != nodeids[i])
        dserror("ERROR: IntElement: Inconsistency Nodes and NodeIds!");

  nodes_.clear();
  nodes_ptr_.clear();
  std::vector<int> empty_dofs(3, -2);

  for (int i = 0; i < numnode; ++i)
    nodes_.push_back(
        MortarNode(nodeids[i], nodes[i]->X(), nodes[i]->Owner(), 3, empty_dofs, isslave));
  for (int i = 0; i < numnode; ++i) nodes_ptr_.push_back(&(nodes_[i]));

  if (numnode > 0) BuildNodalPointers(&nodes[0]);

  // as discretization is already evaluated, compute area
  // (data container has to be initialized first)
  InitializeDataContainer();
  MoData().Area() = ComputeArea();

  return;
}

/*----------------------------------------------------------------------*
 |  map IntElement coords to Element coords (public)          popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::IntElement::MapToParent(const double* xi, double* parxi)
{
  // outdated (popp 05/2016)
  // - affine mapping is only correct for undistorted planar elements
  // - in general we need nonlinear projection procedure
  dserror("ERROR: MapToParent() function is outdated");

  // *********************************************************************
  // do mapping for given IntElement and Element
  // *********************************************************** quad9 ***
  if (ParShape() == DRT::Element::quad9)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        parxi[0] = 0.5 * xi[0] - 0.5;
        parxi[1] = 0.5 * xi[1] - 0.5;
        break;
      }
      case 1:
      {
        parxi[0] = 0.5 * xi[0] + 0.5;
        parxi[1] = 0.5 * xi[1] - 0.5;
        break;
      }
      case 2:
      {
        parxi[0] = 0.5 * xi[0] + 0.5;
        parxi[1] = 0.5 * xi[1] + 0.5;
        break;
      }
      case 3:
      {
        parxi[0] = 0.5 * xi[0] - 0.5;
        parxi[1] = 0.5 * xi[1] + 0.5;
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // *********************************************************** quad8 ***
  else if (ParShape() == DRT::Element::quad8)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        parxi[0] = xi[0] - 1.0;
        parxi[1] = xi[1] - 1.0;
        break;
      }
      case 1:
      {
        parxi[0] = -xi[1] + 1.0;
        parxi[1] = xi[0] - 1.0;
        break;
      }
      case 2:
      {
        parxi[0] = -xi[0] + 1.0;
        parxi[1] = -xi[1] + 1.0;
        break;
      }
      case 3:
      {
        parxi[0] = xi[1] - 1.0;
        parxi[1] = -xi[0] + 1.0;
        break;
      }
      case 4:
      {
        parxi[0] = 0.5 * xi[0] - 0.5 * xi[1];
        parxi[1] = 0.5 * xi[0] + 0.5 * xi[1];
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // ************************************************************ tri6 ***
  else if (ParShape() == DRT::Element::tri6)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        parxi[0] = 0.5 * xi[0];
        parxi[1] = 0.5 * xi[1];
        break;
      }
      case 1:
      {
        parxi[0] = 0.5 * xi[0] + 0.5;
        parxi[1] = 0.5 * xi[1];
        break;
      }
      case 2:
      {
        parxi[0] = 0.5 * xi[0];
        parxi[1] = 0.5 * xi[1] + 0.5;
        break;
      }
      case 3:
      {
        parxi[0] = -0.5 * xi[0] + 0.5;
        parxi[1] = -0.5 * xi[1] + 0.5;
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // *********************************************************** quad4 ***
  else if (ParShape() == DRT::Element::quad4)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        parxi[0] = xi[0];
        parxi[1] = xi[1];
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // ************************************************************ tri3 ***
  else if (ParShape() == DRT::Element::tri3)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        parxi[0] = xi[0];
        parxi[1] = xi[1];
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // ************************************************************ nurbs9 ***
  else if (ParShape() == DRT::Element::nurbs9)
  {
    if (Lid() != 0) dserror("nurbs9 should only have one integration element");
    // TODO: There is not necessarily a constant mapping from the IntEle
    // to the parent ele. Actually, we want to have the GP at a certain
    // spatial location, which is determined by the integration element.
    // However, for higher order Elements, the projection from the (bi-)
    // linear IntEle to a distorted higher order (e.g. NURBS) ele
    // might be more complicated. It is still to be seen, if this has
    // a notable effect.
    if (!rewind_)
    {
      parxi[0] = xi[0];
      parxi[1] = xi[1];
    }
    else
    {
      parxi[0] = xi[1];
      parxi[1] = xi[0];
    }
  }
  // ************************************************************ nurbs9 ***

  // ********************************************************* invalid ***
  else
    dserror("ERROR: MapToParent called for invalid parent element type!");
  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  map IntElement coord derivatives to Element (public)      popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::IntElement::MapToParent(const std::vector<GEN::pairedvector<int, double>>& dxi,
    std::vector<GEN::pairedvector<int, double>>& dparxi)
{
  // outdated (popp 05/2016)
  // - affine mapping is only correct for undistorted planar elements
  // - in general we need nonlinear projection procedure
  dserror("ERROR: MapToParent() function is outdated");

  // map iterator
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // *********************************************************************
  // do mapping for given IntElement and Element
  // *********************************************************** quad9 ***
  if (ParShape() == DRT::Element::quad9)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      case 1:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      case 2:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      case 3:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // *********************************************************** quad8 ***
  else if (ParShape() == DRT::Element::quad8)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[0][p->first] += (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[1][p->first] += (p->second);
        break;
      }
      case 1:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[1][p->first] += (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[0][p->first] -= (p->second);
        break;
      }
      case 2:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[0][p->first] -= (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[1][p->first] -= (p->second);
        break;
      }
      case 3:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[1][p->first] -= (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[0][p->first] += (p->second);
        break;
      }
      case 4:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
        {
          dparxi[0][p->first] += 0.5 * (p->second);
          dparxi[1][p->first] += 0.5 * (p->second);
        }
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
        {
          dparxi[0][p->first] -= 0.5 * (p->second);
          dparxi[1][p->first] += 0.5 * (p->second);
        }
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // ************************************************************ tri6 ***
  else if (ParShape() == DRT::Element::tri6)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      case 1:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      case 2:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] += 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] += 0.5 * (p->second);
        break;
      }
      case 3:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p)
          dparxi[0][p->first] -= 0.5 * (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p)
          dparxi[1][p->first] -= 0.5 * (p->second);
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // *********************************************************** quad4 ***
  else if (ParShape() == DRT::Element::quad4)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[0][p->first] = (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[1][p->first] = (p->second);
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // ************************************************************ tri3 ***
  else if (ParShape() == DRT::Element::tri3)
  {
    // do mapping according to sub-element id
    switch (Lid())
    {
      case 0:
      {
        for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[0][p->first] = (p->second);
        for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[1][p->first] = (p->second);
        break;
      }
      default:
      {
        dserror("ERROR: MapToParent: Invalid local IntElement Id!");
        break;
      }
    }
  }
  // ************************************************************ nurbs9 ***
  else if (ParShape() == DRT::Element::nurbs9)
  {
    if (Lid() != 0) dserror("nurbs9 should only have one integration element");
    if (!rewind_)
    {
      for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[0][p->first] = (p->second);
      for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[1][p->first] = (p->second);
    }
    else
    {
      for (CI p = dxi[1].begin(); p != dxi[1].end(); ++p) dparxi[0][p->first] = (p->second);
      for (CI p = dxi[0].begin(); p != dxi[0].end(); ++p) dparxi[1][p->first] = (p->second);
    }
  }
  // ************************************************************ nurbs9 ***

  // ********************************************************* invalid ***
  else
    dserror("ERROR: MapToParent called for invalid parent element type!");
  // *********************************************************************

  return true;
}

void MORTAR::IntElement::NodeLinearization(
    std::vector<std::vector<GEN::pairedvector<int, double>>>& nodelin)
{
  switch (parele_->Shape())
  {
    // for all Lagrange Finite elements we can associate them directly with
    // the interpolatory nodes of the parent element
    case quad4:
    case quad8:
    case quad9:
    case tri3:
    case tri6:
    {
      // resize the linearizations
      nodelin.resize(NumNode(), std::vector<GEN::pairedvector<int, double>>(3, 1));

      // loop over all intEle nodes
      for (int in = 0; in < NumNode(); ++in)
      {
        MORTAR::MortarNode* mrtrnode = dynamic_cast<MORTAR::MortarNode*>(Nodes()[in]);
        for (int dim = 0; dim < 3; ++dim) nodelin[in][dim][mrtrnode->Dofs()[dim]] += 1.;
      }
      break;
    }
    case nurbs9:
    {
      // resize the linearizations
      nodelin.resize(
          NumNode(), std::vector<GEN::pairedvector<int, double>>(3, 3 * (parele_->NumNode())));

      // parameter space coords of pseudo nodes
      double pseudo_nodes_param_coords[4][2];

      if (rewind_)
      {
        pseudo_nodes_param_coords[0][0] = -1.;
        pseudo_nodes_param_coords[0][1] = -1.;
        pseudo_nodes_param_coords[1][0] = -1.;
        pseudo_nodes_param_coords[1][1] = +1.;
        pseudo_nodes_param_coords[2][0] = +1.;
        pseudo_nodes_param_coords[2][1] = +1.;
        pseudo_nodes_param_coords[3][0] = +1.;
        pseudo_nodes_param_coords[3][1] = -1.;
      }
      else
      {
        pseudo_nodes_param_coords[0][0] = -1.;
        pseudo_nodes_param_coords[0][1] = -1.;
        pseudo_nodes_param_coords[1][0] = +1.;
        pseudo_nodes_param_coords[1][1] = -1.;
        pseudo_nodes_param_coords[2][0] = +1.;
        pseudo_nodes_param_coords[2][1] = +1.;
        pseudo_nodes_param_coords[3][0] = -1.;
        pseudo_nodes_param_coords[3][1] = +1.;
      }

      // loop over all pseudo-nodes
      for (int pn = 0; pn < NumNode(); ++pn)
      {
        double xi[2] = {pseudo_nodes_param_coords[pn][0], pseudo_nodes_param_coords[pn][1]};

        // evaluate shape functions at pseudo node param coords
        LINALG::SerialDenseVector sval(9);
        LINALG::SerialDenseMatrix sderiv(9, 2);
        parele_->EvaluateShape(xi, sval, sderiv, 9, true);

        // loop over all parent element control points
        for (int cp = 0; cp < parele_->NumNode(); ++cp)
        {
          MORTAR::MortarNode* mrtrcp = dynamic_cast<MORTAR::MortarNode*>(parele_->Nodes()[cp]);

          // loop over all dimensions
          for (int dim = 0; dim < 3; ++dim) nodelin.at(pn).at(dim)[mrtrcp->Dofs()[dim]] += sval(cp);
        }
      }
      break;
    }
    default:
    {
      dserror("unknown type of parent element shape");
      break;
    }
  }
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::IntCell::IntCell(int id, int nvertices, LINALG::Matrix<3, 3>& coords, double* auxn,
    const DRT::Element::DiscretizationType& shape,
    std::vector<GEN::pairedvector<int, double>>& linv1,
    std::vector<GEN::pairedvector<int, double>>& linv2,
    std::vector<GEN::pairedvector<int, double>>& linv3,
    std::vector<GEN::pairedvector<int, double>>& linauxn)
    : id_(id), slaveId_(-1), masterId_(-1), nvertices_(nvertices), coords_(coords), shape_(shape)
{
  // store auxiliary plane normal
  for (int k = 0; k < 3; ++k) Auxn()[k] = auxn[k];

  if (shape == DRT::Element::tri3)
  {
    // compute area of IntCell
    double t1[3] = {0.0, 0.0, 0.0};
    double t2[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; ++k)
    {
      t1[k] = Coords()(k, 1) - Coords()(k, 0);
      t2[k] = Coords()(k, 2) - Coords()(k, 0);
    }

    double t1xt2[3] = {0.0, 0.0, 0.0};
    t1xt2[0] = t1[1] * t2[2] - t1[2] * t2[1];
    t1xt2[1] = t1[2] * t2[0] - t1[0] * t2[2];
    t1xt2[2] = t1[0] * t2[1] - t1[1] * t2[0];
    area_ = 0.5 * sqrt(t1xt2[0] * t1xt2[0] + t1xt2[1] * t1xt2[1] + t1xt2[2] * t1xt2[2]);
  }
  else if (shape == DRT::Element::line2)
  {
    // compute length of IntLine
    double v[3] = {0.0, 0.0, 0.0};
    v[0] = Coords()(0, 0) - Coords()(0, 1);
    v[1] = Coords()(1, 0) - Coords()(1, 1);
    v[2] = Coords()(2, 0) - Coords()(2, 1);

    area_ = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (area_ < 1e-12)
    {
      std::cout << "v0 = " << Coords()(0, 0) << "  " << Coords()(1, 0) << "  " << Coords()(2, 0)
                << std::endl;
      std::cout << "v1 = " << Coords()(0, 1) << "  " << Coords()(1, 1) << "  " << Coords()(2, 1)
                << std::endl;
      dserror("ERROR: INTCELL has no length!");
    }
  }

  // store vertex linearizations and auxn linearization
  linvertex_.resize(3);
  linvertex_[0] = linv1;
  linvertex_[1] = linv2;
  linvertex_[2] = linv3;  // dummy for line2
  linauxn_ = linauxn;

  return;
}


/*----------------------------------------------------------------------*
 |  Get global coords for given local coords (IntCell)        popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::IntCell::LocalToGlobal(const double* xi, double* globcoord, int inttype)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  if (Shape() == DRT::Element::tri3 or Shape() == DRT::Element::line2)
  {
    // collect fundamental data
    LINALG::Matrix<3, 1> val;
    LINALG::Matrix<3, 2> deriv;

    // Evaluate shape, get nodal coords and interpolate global coords
    EvaluateShape(xi, val, deriv);
    for (int i = 0; i < 3; ++i) globcoord[i] = 0.0;

    for (int i = 0; i < NumVertices(); ++i)
    {
      if (inttype == 0)
      {
        // use shape function values for interpolation
        globcoord[0] += val(i) * Coords()(0, i);
        globcoord[1] += val(i) * Coords()(1, i);
        globcoord[2] += val(i) * Coords()(2, i);
      }
      else if (inttype == 1)
      {
        // use shape function derivatives xi for interpolation
        globcoord[0] += deriv(i, 0) * Coords()(0, i);
        globcoord[1] += deriv(i, 0) * Coords()(1, i);
        globcoord[2] += deriv(i, 0) * Coords()(2, i);
      }
      else if (inttype == 2)
      {
        if (Shape() == DRT::Element::line2)
          dserror("ERROR: for line2 elements only 1 parameter space coordinate");

        // use shape function derivatives eta for interpolation
        globcoord[0] += deriv(i, 1) * Coords()(0, i);
        globcoord[1] += deriv(i, 1) * Coords()(1, i);
        globcoord[2] += deriv(i, 1) * Coords()(2, i);
      }
      else
        dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");
    }
  }


  return true;
}

/*----------------------------------------------------------------------*
 |  output for integration cell                              farah 01/16|
 *----------------------------------------------------------------------*/
void MORTAR::IntCell::Print()
{
  std::cout << "Slave  ID= " << GetSlaveId() << std::endl;
  std::cout << "Master ID= " << GetMasterId() << std::endl;
  std::cout << "Coordinates for vertex 0 = " << Coords()(0, 0) << " " << Coords()(1, 0) << " "
            << Coords()(2, 0) << std::endl;
  std::cout << "Coordinates for vertex 1 = " << Coords()(0, 1) << " " << Coords()(1, 1) << " "
            << Coords()(2, 1) << std::endl;
  std::cout << "Coordinates for vertex 2 = " << Coords()(0, 2) << " " << Coords()(1, 2) << " "
            << Coords()(2, 2) << std::endl;

  return;
}


/*----------------------------------------------------------------------*
 |  Evaluate shape functions (IntCell)                        popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::IntCell::EvaluateShape(
    const double* xi, LINALG::Matrix<3, 1>& val, LINALG::Matrix<3, 2>& deriv)
{
  if (!xi) dserror("ERROR: EvaluateShape (IntCell) called with xi=NULL");

  // 3noded triangular element
  if (Shape() == DRT::Element::tri3)
  {
    val(0) = 1.0 - xi[0] - xi[1];
    val(1) = xi[0];
    val(2) = xi[1];
    deriv(0, 0) = -1.0;
    deriv(0, 1) = -1.0;
    deriv(1, 0) = 1.0;
    deriv(1, 1) = 0.0;
    deriv(2, 0) = 0.0;
    deriv(2, 1) = 1.0;
  }
  else if (Shape() == DRT::Element::line2)
  {
    val(0) = 0.5 * (1 - xi[0]);
    val(1) = 0.5 * (1 + xi[0]);
    deriv(0, 0) = -0.5;
    deriv(1, 0) = 0.5;
  }

  // unknown case
  else
    dserror("ERROR: EvaluateShape (IntCell) called for type != tri3/line2");

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant (IntCell)                   popp 11/08|
 *----------------------------------------------------------------------*/
double MORTAR::IntCell::Jacobian()
{
  double jac = 0.0;

  // 2D linear case (2noded line element)
  if (Shape() == DRT::Element::tri3)
    jac = Area() * 2.0;
  else if (Shape() == DRT::Element::line2)
    jac = Area() * 0.5;
  // unknown case
  else
    dserror("ERROR: Jacobian (IntCell) called for unknown ele type!");

  return jac;
}

/*----------------------------------------------------------------------*
 |  Evaluate directional deriv. of Jacobian det. AuxPlane     popp 03/09|
 *----------------------------------------------------------------------*/
void MORTAR::IntCell::DerivJacobian(GEN::pairedvector<int, double>& derivjac)
{
  // define iterator
  typedef GEN::pairedvector<int, double>::const_iterator CI;

  // 1d line element
  if (Shape() == DRT::Element::line2)
  {
    // compute length of IntLine
    double v[3] = {0.0, 0.0, 0.0};
    v[0] = Coords()(0, 0) - Coords()(0, 1);
    v[1] = Coords()(1, 0) - Coords()(1, 1);
    v[2] = Coords()(2, 0) - Coords()(2, 1);

    double l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    double linv = 1.0 / l;
    double fac = 0.25 * linv;

    // linearizarion of v
    std::vector<GEN::pairedvector<int, double>> vg(3, 1000);

    // first entry (x component lin)
    for (CI p = GetDerivVertex(0)[0].begin(); p != GetDerivVertex(0)[0].end(); ++p)
      vg[0][p->first] += (p->second);
    for (CI p = GetDerivVertex(1)[0].begin(); p != GetDerivVertex(1)[0].end(); ++p)
      vg[0][p->first] -= (p->second);

    // first entry (y component lin)
    for (CI p = GetDerivVertex(0)[1].begin(); p != GetDerivVertex(0)[1].end(); ++p)
      vg[1][p->first] += (p->second);
    for (CI p = GetDerivVertex(1)[1].begin(); p != GetDerivVertex(1)[1].end(); ++p)
      vg[1][p->first] -= (p->second);

    // first entry (z component lin)
    for (CI p = GetDerivVertex(0)[2].begin(); p != GetDerivVertex(0)[2].end(); ++p)
      vg[2][p->first] += (p->second);
    for (CI p = GetDerivVertex(1)[2].begin(); p != GetDerivVertex(1)[2].end(); ++p)
      vg[2][p->first] -= (p->second);

    // linearizarion of v^t * v
    GEN::pairedvector<int, double> vv(1000);

    // delta v^T * v
    for (CI p = vg[0].begin(); p != vg[0].end(); ++p) vv[p->first] += v[0] * (p->second);
    for (CI p = vg[1].begin(); p != vg[1].end(); ++p) vv[p->first] += v[1] * (p->second);
    for (CI p = vg[2].begin(); p != vg[2].end(); ++p) vv[p->first] += v[2] * (p->second);

    // v^T * delta v
    for (CI p = vg[0].begin(); p != vg[0].end(); ++p) vv[p->first] += v[0] * (p->second);
    for (CI p = vg[1].begin(); p != vg[1].end(); ++p) vv[p->first] += v[1] * (p->second);
    for (CI p = vg[2].begin(); p != vg[2].end(); ++p) vv[p->first] += v[2] * (p->second);

    // fac * vv
    for (CI p = vv.begin(); p != vv.end(); ++p) derivjac[p->first] += fac * (p->second);
  }
  // 2D linear case (2noded line element)
  else if (Shape() == DRT::Element::tri3)
  {
    // metrics routine gives local basis vectors
    static std::vector<double> gxi(3);
    static std::vector<double> geta(3);

    for (int k = 0; k < 3; ++k)
    {
      gxi[k] = Coords()(k, 1) - Coords()(k, 0);
      geta[k] = Coords()(k, 2) - Coords()(k, 0);
    }

    // cross product of gxi and geta
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = gxi[1] * geta[2] - gxi[2] * geta[1];
    cross[1] = gxi[2] * geta[0] - gxi[0] * geta[2];
    cross[2] = gxi[0] * geta[1] - gxi[1] * geta[0];

    // inverse jacobian
    const double jacinv =
        1.0 / sqrt(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2]);


    // *********************************************************************
    // compute Jacobian derivative
    // *********************************************************************
    // first vertex (Coords(k,0)) is part of gxi and geta
    for (CI p = GetDerivVertex(0)[0].begin(); p != GetDerivVertex(0)[0].end(); ++p)
    {
      derivjac[p->first] -= jacinv * cross[1] * gxi[2] * (p->second);
      derivjac[p->first] += jacinv * cross[1] * geta[2] * (p->second);
      derivjac[p->first] += jacinv * cross[2] * gxi[1] * (p->second);
      derivjac[p->first] -= jacinv * cross[2] * geta[1] * (p->second);
    }
    for (CI p = GetDerivVertex(0)[1].begin(); p != GetDerivVertex(0)[1].end(); ++p)
    {
      derivjac[p->first] += jacinv * cross[0] * gxi[2] * (p->second);
      derivjac[p->first] -= jacinv * cross[0] * geta[2] * (p->second);
      derivjac[p->first] -= jacinv * cross[2] * gxi[0] * (p->second);
      derivjac[p->first] += jacinv * cross[2] * geta[0] * (p->second);
    }
    for (CI p = GetDerivVertex(0)[2].begin(); p != GetDerivVertex(0)[2].end(); ++p)
    {
      derivjac[p->first] -= jacinv * cross[0] * gxi[1] * (p->second);
      derivjac[p->first] += jacinv * cross[0] * geta[1] * (p->second);
      derivjac[p->first] += jacinv * cross[1] * gxi[0] * (p->second);
      derivjac[p->first] -= jacinv * cross[1] * geta[0] * (p->second);
    }

    // second vertex (Coords(k,1)) is part of gxi
    for (CI p = GetDerivVertex(1)[0].begin(); p != GetDerivVertex(1)[0].end(); ++p)
    {
      derivjac[p->first] -= jacinv * cross[1] * geta[2] * (p->second);
      derivjac[p->first] += jacinv * cross[2] * geta[1] * (p->second);
    }
    for (CI p = GetDerivVertex(1)[1].begin(); p != GetDerivVertex(1)[1].end(); ++p)
    {
      derivjac[p->first] += jacinv * cross[0] * geta[2] * (p->second);
      derivjac[p->first] -= jacinv * cross[2] * geta[0] * (p->second);
    }
    for (CI p = GetDerivVertex(1)[2].begin(); p != GetDerivVertex(1)[2].end(); ++p)
    {
      derivjac[p->first] -= jacinv * cross[0] * geta[1] * (p->second);
      derivjac[p->first] += jacinv * cross[1] * geta[0] * (p->second);
    }

    // third vertex (Coords(k,2)) is part of geta
    for (CI p = GetDerivVertex(2)[0].begin(); p != GetDerivVertex(2)[0].end(); ++p)
    {
      derivjac[p->first] += jacinv * cross[1] * gxi[2] * (p->second);
      derivjac[p->first] -= jacinv * cross[2] * gxi[1] * (p->second);
    }
    for (CI p = GetDerivVertex(2)[1].begin(); p != GetDerivVertex(2)[1].end(); ++p)
    {
      derivjac[p->first] -= jacinv * cross[0] * gxi[2] * (p->second);
      derivjac[p->first] += jacinv * cross[2] * gxi[0] * (p->second);
    }
    for (CI p = GetDerivVertex(2)[2].begin(); p != GetDerivVertex(2)[2].end(); ++p)
    {
      derivjac[p->first] += jacinv * cross[0] * gxi[1] * (p->second);
      derivjac[p->first] -= jacinv * cross[1] * gxi[0] * (p->second);
    }
  }

  // unknown case
  else
    dserror("ERROR: DerivJacobian (IntCell) called for unknown ele type!");

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Vertex::Vertex(std::vector<double> coord, Vertex::vType type, std::vector<int> nodeids,
    Vertex* next, Vertex* prev, bool intersect, bool entryexit, Vertex* neighbor, double alpha)
    : coord_(coord),
      type_(type),
      nodeids_(nodeids),
      next_(next),
      prev_(prev),
      intersect_(intersect),
      entryexit_(entryexit),
      neighbor_(neighbor),
      alpha_(alpha)
{
  // empty constructor body
  return;
}

/*----------------------------------------------------------------------*
 |  cctor (public)                                            popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Vertex::Vertex(const Vertex& old)
    : coord_(old.coord_),
      type_(old.type_),
      nodeids_(old.nodeids_),
      next_(old.next_),
      prev_(old.prev_),
      intersect_(old.intersect_),
      entryexit_(old.entryexit_),
      neighbor_(old.neighbor_),
      alpha_(old.alpha_)
{
  // empty copy constructor body
  return;
}
