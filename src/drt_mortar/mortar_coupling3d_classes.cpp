/*!----------------------------------------------------------------------
\file mortar_coupling3d_classes.cpp
\brief A class for mortar coupling of ONE slave element and ONE master
       element of a mortar interface in 3D (definition of sub-classes).

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/

#include "mortar_coupling3d.H"
#include "mortar_coupling3d_classes.H"
#include "mortar_projector.H"
#include "mortar_integrator.H"
#include "mortar_defines.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../drt_lib/drt_node.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 03/09|
 *----------------------------------------------------------------------*/
MORTAR::IntElement::IntElement(int lid, int id, int owner,
                               const DRT::Element::DiscretizationType& parshape,
                               const DRT::Element::DiscretizationType& shape,
                               const int numnode,
                               const int* nodeids,
                               std::vector<DRT::Node*> nodes,
                               const bool isslave) :
MORTAR::MortarElement(id,owner,shape,numnode,nodeids,isslave),
lid_(lid),
parshape_(parshape)
{
  // check for consistency of nodeids and nodes
  for (int i=0;i<NumNode();++i)
    if (nodes[i]->Id()!=nodeids[i])
      dserror("ERROR: IntElement: Inconsistency Nodes and NodeIds!");

  // store given nodes into class variable
  quadnode_.resize((int)nodes.size());
  for (int i=0;i<(int)nodes.size();++i)
    quadnode_[i] = nodes[i];

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
  // *********************************************************************
  // do mapping for given IntElement and Element
  // *********************************************************** quad9 ***
  if (ParShape()==DRT::Element::quad9)
  {
    // do mapping according to sub-element id
    switch(Lid())
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
  else if (ParShape()==DRT::Element::quad8)
  {
    // do mapping according to sub-element id
    switch(Lid())
    {
    case 0:
    {
      parxi[0] =  xi[0] - 1;
      parxi[1] =  xi[1] - 1;
      break;
    }
    case 1:
    {
      parxi[0] = -xi[1] + 1;
      parxi[1] =  xi[0] - 1;
      break;
    }
    case 2:
    {
      parxi[0] = -xi[0] + 1;
      parxi[1] = -xi[1] + 1;
      break;
    }
    case 3:
    {
      parxi[0] =  xi[1] - 1;
      parxi[1] = -xi[0] + 1;
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
  else if (ParShape()==DRT::Element::tri6)
  {
    // do mapping according to sub-element id
    switch(Lid())
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
  else if (ParShape()==DRT::Element::quad4)
  {
    // do mapping according to sub-element id
    switch(Lid())
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
  else if (ParShape()==DRT::Element::tri3)
  {
    // do mapping according to sub-element id
    switch(Lid())
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
  // ********************************************************* invalid ***
  else
    dserror("ERROR: MapToParent called for invalid parent element type!");
  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  map IntElement coord derivatives to Element (public)      popp 03/09|
 *----------------------------------------------------------------------*/
bool MORTAR::IntElement::MapToParent(const std::vector<std::map<int,double> >& dxi,
                                     std::vector<std::map<int,double> >& dparxi)
{
  // map iterator
  typedef std::map<int,double>::const_iterator CI;

  // *********************************************************************
  // do mapping for given IntElement and Element
  // *********************************************************** quad9 ***
  if (ParShape()==DRT::Element::quad9)
  {
    // do mapping according to sub-element id
    switch(Lid())
    {
    case 0:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 0.5 * (p->second);
      break;
    }
    case 1:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 0.5 * (p->second);
      break;
    }
    case 2:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 0.5 * (p->second);
      break;
    }
    case 3:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
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
  else if (ParShape()==DRT::Element::quad8)
  {
    // do mapping according to sub-element id
    switch(Lid())
    {
    case 0:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 1.0 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 1.0 * (p->second);
      break;
    }
    case 1:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[1][p->first] += 1.0 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[0][p->first] -= 1.0 * (p->second);
      break;
    }
    case 2:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] -= 1.0 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] -= 1.0 * (p->second);
      break;
    }
    case 3:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[1][p->first] -= 1.0 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[0][p->first] += 1.0 * (p->second);
      break;
    }
    case 4:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
      {
        dparxi[0][p->first] += 0.5 * (p->second);
        dparxi[1][p->first] += 0.5 * (p->second);
      }
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
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
  else if (ParShape()==DRT::Element::tri6)
  {
    // do mapping according to sub-element id
    switch(Lid())
    {
    case 0:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 0.5 * (p->second);
      break;
    }
    case 1:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 0.5 * (p->second);
      break;
    }
    case 2:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] += 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] += 0.5 * (p->second);
      break;
    }
    case 3:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] -= 0.5 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
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
  else if (ParShape()==DRT::Element::quad4)
  {
    // do mapping according to sub-element id
    switch(Lid())
    {
    case 0:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] = 1.0 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] = 1.0 * (p->second);
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
  else if (ParShape()==DRT::Element::tri3)
  {
    // do mapping according to sub-element id
    switch(Lid())
    {
    case 0:
    {
      for (CI p=dxi[0].begin();p!=dxi[0].end();++p)
        dparxi[0][p->first] = 1.0 * (p->second);
      for (CI p=dxi[1].begin();p!=dxi[1].end();++p)
        dparxi[1][p->first] = 1.0 * (p->second);
      break;
    }
    default:
    {
      dserror("ERROR: MapToParent: Invalid local IntElement Id!");
      break;
    }
    }
  }
  // ********************************************************* invalid ***
  else
    dserror("ERROR: MapToParent called for invalid parent element type!");
  // *********************************************************************

  return true;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::IntCell::IntCell(int id, int nvertices, Epetra_SerialDenseMatrix& coords,
                         double* auxn, const DRT::Element::DiscretizationType& shape,
                         std::vector<std::map<int,double> >& linv1,
                         std::vector<std::map<int,double> >& linv2,
                         std::vector<std::map<int,double> >& linv3,
                         std::vector<std::map<int,double> >& linauxn) :
id_(id),
nvertices_(nvertices),
coords_(coords),
shape_(shape)
{
   // check nvertices_ and shape_
  if (nvertices_!=3) dserror("ERROR: Integration cell must have 3 vertices");
  if (shape_!=DRT::Element::tri3) dserror("ERROR: Integration cell must be tri3");

  // check dimensions of coords_
  if (coords_.M() != 3) dserror("ERROR: Inconsistent coord matrix");
  if (coords_.N() != nvertices_) dserror("ERROR: Inconsistent coord matrix");

  // store auxiliary plane normal
  for (int k=0;k<3;++k) Auxn()[k] = auxn[k];

  // compute area of IntCell
  double t1[3] = {0.0, 0.0, 0.0};
  double t2[3] = {0.0, 0.0, 0.0};
  for (int k=0;k<3;++k)
  {
    t1[k]=Coords()(k,1)-Coords()(k,0);
    t2[k]=Coords()(k,2)-Coords()(k,0);
  }

  double t1xt2[3] = {0.0, 0.0, 0.0};
  t1xt2[0] = t1[1]*t2[2]-t1[2]*t2[1];
  t1xt2[1] = t1[2]*t2[0]-t1[0]*t2[2];
  t1xt2[2] = t1[0]*t2[1]-t1[1]*t2[0];
  area_ = 0.5*sqrt(t1xt2[0]*t1xt2[0]+t1xt2[1]*t1xt2[1]+t1xt2[2]*t1xt2[2]);

  // store vertex linearizations and auxn linearization
  linvertex_.resize(3);
  linvertex_[0] = linv1;
  linvertex_[1] = linv2;
  linvertex_[2] = linv3;
  linauxn_ = linauxn;

  return;
}

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords (IntCell)        popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::IntCell::LocalToGlobal(const double* xi,
                                    double* globcoord,
                                    int inttype)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  // collect fundamental data
  int nnodes = NumVertices();
  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes,2,true);

  // Evaluate shape, get nodal coords and interpolate global coords
  EvaluateShape(xi, val, deriv);
  for (int i=0;i<3;++i) globcoord[i]=0.0;

  for (int i=0;i<nnodes;++i)
  {
    if (inttype==0)
    {
      // use shape function values for interpolation
      globcoord[0]+=val[i]*Coords()(0,i);
      globcoord[1]+=val[i]*Coords()(1,i);
      globcoord[2]+=val[i]*Coords()(2,i);
    }
    else if (inttype==1)
    {
      // use shape function derivatives xi for interpolation
      globcoord[0]+=deriv(i,0)*Coords()(0,i);
      globcoord[1]+=deriv(i,0)*Coords()(1,i);
      globcoord[2]+=deriv(i,0)*Coords()(2,i);
    }
    else if (inttype==2)
    {
      // use shape function derivatives eta for interpolation
      globcoord[0]+=deriv(i,1)*Coords()(0,i);
      globcoord[1]+=deriv(i,1)*Coords()(1,i);
      globcoord[2]+=deriv(i,1)*Coords()(2,i);
    }
    else
      dserror("ERROR: Invalid interpolation type requested, only 0,1,2!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate shape functions (IntCell)                        popp 11/08|
 *----------------------------------------------------------------------*/
bool MORTAR::IntCell::EvaluateShape(const double* xi,
                                    LINALG::SerialDenseVector& val,
                                    LINALG::SerialDenseMatrix& deriv)
{
  if (!xi)
    dserror("ERROR: EvaluateShape (IntCell) called with xi=NULL");

  // 3noded triangular element
  if(Shape()==DRT::Element::tri3)
  {
    val[0] = 1-xi[0]-xi[1];
    val[1] = xi[0];
    val[2] = xi[1];
    deriv(0,0) = -1.0; deriv(0,1) = -1.0;
    deriv(1,0) =  1.0; deriv(1,1) =  0.0;
    deriv(2,0) =  0.0; deriv(2,1) =  1.0;
  }

  // unknown case
  else dserror("ERROR: EvaluateShape (IntCell) called for type != tri3");

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant (IntCell)                   popp 11/08|
 *----------------------------------------------------------------------*/
double MORTAR::IntCell::Jacobian(double* xi)
{
  double jac = 0.0;
  std::vector<double> gxi(3);
  std::vector<double> geta(3);

  // 2D linear case (2noded line element)
  if (Shape()==DRT::Element::tri3)
    jac = Area()*2;

  // unknown case
  else dserror("ERROR: Jacobian (IntCell) called for unknown ele type!");

  return jac;
}

/*----------------------------------------------------------------------*
 |  Evaluate directional deriv. of Jacobian det. AuxPlane     popp 03/09|
 *----------------------------------------------------------------------*/
void MORTAR::IntCell::DerivJacobian(double* xi, std::map<int,double>& derivjac)
{
  // metrics routine gives local basis vectors
  std::vector<double> gxi(3);
  std::vector<double> geta(3);

  for (int k=0;k<3;++k)
  {
    gxi[k]=Coords()(k,1)-Coords()(k,0);
    geta[k]=Coords()(k,2)-Coords()(k,0);
  }

  // cross product of gxi and geta
  double cross[3] = {0.0, 0.0, 0.0};
  cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
  cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
  cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];

  double jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  typedef std::map<int,double>::const_iterator CI;

  // 2D linear case (2noded line element)
  if (Shape()==DRT::Element::tri3)
  {
    // *********************************************************************
    // compute Jacobian derivative
    // *********************************************************************
    // first vertex (Coords(k,0)) is part of gxi and geta
    for (CI p=GetDerivVertex(0)[0].begin();p!=GetDerivVertex(0)[0].end();++p)
    {
      derivjac[p->first] -= 1/jac*cross[1]*gxi[2]*(p->second);
      derivjac[p->first] += 1/jac*cross[1]*geta[2]*(p->second);
      derivjac[p->first] += 1/jac*cross[2]*gxi[1]*(p->second);
      derivjac[p->first] -= 1/jac*cross[2]*geta[1]*(p->second);
    }
    for (CI p=GetDerivVertex(0)[1].begin();p!=GetDerivVertex(0)[1].end();++p)
    {
      derivjac[p->first] += 1/jac*cross[0]*gxi[2]*(p->second);
      derivjac[p->first] -= 1/jac*cross[0]*geta[2]*(p->second);
      derivjac[p->first] -= 1/jac*cross[2]*gxi[0]*(p->second);
      derivjac[p->first] += 1/jac*cross[2]*geta[0]*(p->second);
    }
    for (CI p=GetDerivVertex(0)[2].begin();p!=GetDerivVertex(0)[2].end();++p)
    {
      derivjac[p->first] -= 1/jac*cross[0]*gxi[1]*(p->second);
      derivjac[p->first] += 1/jac*cross[0]*geta[1]*(p->second);
      derivjac[p->first] += 1/jac*cross[1]*gxi[0]*(p->second);
      derivjac[p->first] -= 1/jac*cross[1]*geta[0]*(p->second);
    }

    // second vertex (Coords(k,1)) is part of gxi
    for (CI p=GetDerivVertex(1)[0].begin();p!=GetDerivVertex(1)[0].end();++p)
    {
      derivjac[p->first] -= 1/jac*cross[1]*geta[2]*(p->second);
      derivjac[p->first] += 1/jac*cross[2]*geta[1]*(p->second);
    }
    for (CI p=GetDerivVertex(1)[1].begin();p!=GetDerivVertex(1)[1].end();++p)
    {
      derivjac[p->first] += 1/jac*cross[0]*geta[2]*(p->second);
      derivjac[p->first] -= 1/jac*cross[2]*geta[0]*(p->second);
    }
    for (CI p=GetDerivVertex(1)[2].begin();p!=GetDerivVertex(1)[2].end();++p)
    {
      derivjac[p->first] -= 1/jac*cross[0]*geta[1]*(p->second);
      derivjac[p->first] += 1/jac*cross[1]*geta[0]*(p->second);
    }

    // third vertex (Coords(k,2)) is part of geta
    for (CI p=GetDerivVertex(2)[0].begin();p!=GetDerivVertex(2)[0].end();++p)
    {
      derivjac[p->first] += 1/jac*cross[1]*gxi[2]*(p->second);
      derivjac[p->first] -= 1/jac*cross[2]*gxi[1]*(p->second);
    }
    for (CI p=GetDerivVertex(2)[1].begin();p!=GetDerivVertex(2)[1].end();++p)
    {
      derivjac[p->first] -= 1/jac*cross[0]*gxi[2]*(p->second);
      derivjac[p->first] += 1/jac*cross[2]*gxi[0]*(p->second);
    }
    for (CI p=GetDerivVertex(2)[2].begin();p!=GetDerivVertex(2)[2].end();++p)
    {
      derivjac[p->first] += 1/jac*cross[0]*gxi[1]*(p->second);
      derivjac[p->first] -= 1/jac*cross[1]*gxi[0]*(p->second);
    }
  }

  // unknown case
  else dserror("ERROR: DerivJacobian (IntCell) called for unknown ele type!");

  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                             popp 11/08|
 *----------------------------------------------------------------------*/
MORTAR::Vertex::Vertex(std::vector<double> coord, Vertex::vType type, std::vector<int> nodeids,
                        Vertex* next, Vertex* prev, bool intersect,
                        bool entryexit, Vertex* neighbor, double alpha) :
coord_(coord),
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
MORTAR::Vertex::Vertex(const Vertex& old) :
coord_(old.coord_),
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

