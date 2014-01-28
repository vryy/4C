/*!----------------------------------------------------------------------
\file mortar_shape_utils.H

</pre>

<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              farah 01/14 |
 *----------------------------------------------------------------------*/
#ifndef VOLMORTAR_SHAPE_H_
#define VOLMORTAR_SHAPE_H_

/*----------------------------------------------------------------------*
 | Header                                                   farah 01/14 |
 *----------------------------------------------------------------------*/
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_fem_general/drt_utils_integration.H"


/*----------------------------------------------------------------------*
 | Utils                                                    farah 01/14 |
 *----------------------------------------------------------------------*/
namespace VOLMORTAR
{

namespace UTILS
{

template<DRT::Element::DiscretizationType distype>
bool LocalToGlobal(DRT::Element& ele,
                   const double* xi,
                   double* globcoord);

template<DRT::Element::DiscretizationType distype>
double Jacobian(const double* xi, DRT::Element& ele);

template<class V>
void volmortar_shape_function_2D(
  V&                                 funct,   ///< to be filled with shape function values
  double&                            xi0,       ///< xi0 coordinate
  double&                            xi1,       ///< xi1 coordinate
  DRT::Element::DiscretizationType  shape    ///< distuinguish between mortar shape
  );

template<class V>
void volmortar_shape_function_3D(
  V&                                 funct,     ///< to be filled with shape function values
  double&                            r,       ///< xi0 coordinate
  double&                            s,       ///< xi1 coordinate
  double&                            t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType  shape       ///< distuinguish between mortar shape
  );

template<class V>
void volmortar_dualshape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  DRT::Element&                            ele,
  const double&                            xi0,       ///< xi0 coordinate
  const double&                            xi1,       ///< xi1 coordinate
  const DRT::Element::DiscretizationType&  shape    ///< distuinguish between mortar shape
  );

template<class V>
void volmortar_dualshape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  DRT::Element&                            ele,
  const double&                            xi0,       ///< xi0 coordinate
  const double&                            xi1,       ///< xi1 coordinate
  const double&                            xi2,       ///< xi1 coordinate
  const DRT::Element::DiscretizationType&  shape    ///< distuinguish between mortar shape
  );

/*----------------------------------------------------------------------*
 |  Get global coords for given local coords (ref pos)       farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool LocalToGlobal(DRT::Element& ele,
                   const double* xi,
                   double* globcoord)
{
  // check input
  if (!xi) dserror("ERROR: LocalToGlobal called with xi=NULL");
  if (!globcoord) dserror("ERROR: LocalToGlobal called with globcoord=NULL");

  static const int n    = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

  for (int i=0;i<ndim;++i)
    globcoord[i]=0.0;

  LINALG::Matrix<n,1>      val;
  if(ndim==2)      DRT::UTILS::shape_function_2D(val,xi[0],xi[1],distype);
  else if(ndim==3) DRT::UTILS::shape_function_3D(val,xi[0],xi[1],xi[2],distype);
  else dserror("Wrong Dimension");

  LINALG::Matrix<ndim,n>   coord;

  for (int i=0;i<n;++i)
  {
    for(int j=0;j<ndim;++j)
    {
      coord(j,i)   = mynodes[i]->X()[j];

      // use shape function values for interpolation
      globcoord[j] += val(i)*coord(j,i);
    }
  }


  return true;
};

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant                            farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double Jacobian(const double* xi, DRT::Element& ele)
{
  //! nn_: number of master element nodes
  static const int nn = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions ("+1" due to considering only interface elements)
  static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

  double jac = 0.0;
  std::vector<double> gxi(3);
  std::vector<double> geta(3);

  // 2D bilinear case (4noded quadrilateral element)
  // 2D quadratic case (6noded triangular element)
  // 2D serendipity case (8noded quadrilateral element)
  // 2D biquadratic case (9noded quadrilateral element)
  if (distype==DRT::Element::quad4 || distype==DRT::Element::tri6 || distype==DRT::Element::quad8 || distype==DRT::Element::quad9)
  {
    // metrics routine gives local basis vectors
    LINALG::Matrix<ndim,nn> deriv;

    // get shape function values and derivatives at xi
    DRT::UTILS::shape_function_2D_deriv1(deriv,xi[0],xi[1],distype);

    // build basis vectors gxi and geta
    for (int i=0;i<nn;++i)
    {
      // first local basis vector
      gxi[0] += deriv(0,i)*ele.Nodes()[i]->X()[0];
      gxi[1] += deriv(0,i)*ele.Nodes()[i]->X()[1];
      gxi[2] += deriv(0,i)*ele.Nodes()[i]->X()[2];

      // second local basis vector
      geta[0] += deriv(1,i)*ele.Nodes()[i]->X()[0];
      geta[1] += deriv(1,i)*ele.Nodes()[i]->X()[1];
      geta[2] += deriv(1,i)*ele.Nodes()[i]->X()[2];
    }

    // cross product of gxi and geta
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
    cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
    cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
    jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);
  }

  // 3D Hex8
  else if (distype == DRT::Element::hex8)
  {
    LINALG::Matrix<3,8> derivs;
    const double r = xi[0];
    const double s = xi[1];
    const double t = xi[2];

    DRT::UTILS::shape_function_3D_deriv1(derivs,r,s,t,DRT::Element::hex8);


    LINALG::Matrix<8,3> xrefe;
    for (int i=0; i<8; ++i)
    {
      DRT::Node** nodes=ele.Nodes();
      if(!nodes) dserror("Nodes() returned null pointer");

      xrefe(i,0) = nodes[i]->X()[0];
      xrefe(i,1) = nodes[i]->X()[1];
      xrefe(i,2) = nodes[i]->X()[2];
    }

    LINALG::Matrix<3,3> invJ;
    invJ.Clear();

    invJ.Multiply(derivs,xrefe);
    jac = invJ.Invert();
    if (jac <= 0.0)
      dserror("Element Jacobian mapping %10.5e <= 0.0",jac);
  }

  // unknown case
  else
    dserror("ERROR: Jacobian called for unknown element type!");

  return jac;
}
/*----------------------------------------------------------------------*
 |  evaluate std. shape 2D                                   farah 01/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_shape_function_2D(
  V&                                 funct,   ///< to be filled with shape function values
  double&                            xi0,       ///< xi0 coordinate
  double&                            xi1,       ///< xi1 coordinate
  DRT::Element::DiscretizationType  shape    ///< distuinguish between mortar shape
  )
{
  switch (shape)
  {
  // *********************************************************************
  // 2D standard linear shape functions (tri3)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case DRT::Element::tri3:
  {
    funct(0) = 1-xi0-xi1;
    funct(1) = xi0;
    funct(2) = xi1;
    break;
  }
  // *********************************************************************
  // 2D standard bilinear shape functions (quad4)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case DRT::Element::quad4:
  {
    funct(0) = 0.25*(1-xi0)*(1-xi1);
    funct(1) = 0.25*(1+xi0)*(1-xi1);
    funct(2) = 0.25*(1+xi0)*(1+xi1);
    funct(3) = 0.25*(1-xi0)*(1+xi1);
    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate std. shape 3D                                   farah 01/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_shape_function_3D(
  V&                                 funct,     ///< to be filled with shape function values
  double&                            r,       ///< xi0 coordinate
  double&                            s,       ///< xi1 coordinate
  double&                            t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType  shape       ///< distuinguish between mortar shape
  )
{
  const double Q18 = 1.0/8.0;
//  const double Q12 = 1.0/2.0;
//  const double Q14 = 1.0/4.0;

  switch (shape)
  {
  case DRT::Element::hex8:
  {
      const double rp=1.0+r;
      const double rm=1.0-r;
      const double sp=1.0+s;
      const double sm=1.0-s;
      const double tp=1.0+t;
      const double tm=1.0-t;

      funct(0)=Q18*rm*sm*tm;
      funct(1)=Q18*rp*sm*tm;
      funct(2)=Q18*rp*sp*tm;
      funct(3)=Q18*rm*sp*tm;
      funct(4)=Q18*rm*sm*tp;
      funct(5)=Q18*rp*sm*tp;
      funct(6)=Q18*rp*sp*tp;
      funct(7)=Q18*rm*sp*tp;

      break;
  }
  case DRT::Element::tet4:
  {
      const double t1=1.0-r-s-t;
      const double t2=r;
      const double t3=s;
      const double t4=t;

      funct(0)= t1;
      funct(1)= t2;
      funct(2)= t3;
      funct(3)= t4;
      break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate dual shape 2D                                   farah 01/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_dualshape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  DRT::Element&                            ele,
  const double&                            xi0,       ///< xi0 coordinate
  const double&                            xi1,       ///< xi1 coordinate
  const DRT::Element::DiscretizationType&  shape    ///< distuinguish between mortar shape
  )
{
  switch (shape)
  {
  // *********************************************************************
  // 2D dual linear shape functions (tri3)
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case DRT::Element::tri3:
  {
    funct(0) = 3-4*xi0-4*xi1;
    funct(1) = 4*xi0-1;
    funct(2) = 4*xi1-1;

    break;
  }
  // *********************************************************************
  // 2D dual bilinear shape functions (quad4)
  // (used for interpolation of Lagrange mutliplier field)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case DRT::Element::quad4:
  {
    int nnodes = ele.NumNode();
    // get solution matrix with dual parameters
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);
    LINALG::Matrix<4,1>      stdval;

    // establish fundamental data
    double detg = 0.0;

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    MORTAR::ElementIntegrator integrator(ele.Shape());

    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);

    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i,0), integrator.Coordinate(i,1)};
      volmortar_shape_function_2D(stdval,gpc[0],gpc[1],DRT::Element::quad4);
      detg = Jacobian<DRT::Element::quad4>(gpc,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=integrator.Weight(i)*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*stdval(j)*detg;
        }
    }

    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    double xi[2] = {xi0,xi1};
    volmortar_shape_function_2D(stdval,xi[0],xi[1],DRT::Element::quad4);

    // evaluate dual shape functions
    for (int i=0;i<nnodes;++i)
      funct(i)=0;

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
        funct(i)+=ae(i,j)*stdval(j);

    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate dual shape 3D                                   farah 01/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_dualshape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  DRT::Element&                            ele,
  const double&                            xi0,       ///< xi0 coordinate
  const double&                            xi1,       ///< xi1 coordinate
  const double&                            xi2,       ///< xi1 coordinate
  const DRT::Element::DiscretizationType&  shape    ///< distuinguish between mortar shape
  )
{
  switch (shape)
  {
  // *********************************************************************
  // 3D dual trilinear shape functions (hex8)
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case DRT::Element::hex8:
  {
    // initialize GaußRule
    Epetra_SerialDenseMatrix coords;
    std::vector<double> weights;
    int ngp;

    DRT::UTILS::GaussRule3D mygaussrule=DRT::UTILS::intrule_hex_27point;

    const DRT::UTILS::IntegrationPoints3D intpoints(mygaussrule);
    ngp = intpoints.nquad;
    coords.Reshape(ngp,3);
    weights.resize(ngp);
    for (int i=0;i<ngp;++i)
    {
      coords(i,0)=intpoints.qxg[i][0];
      coords(i,1)=intpoints.qxg[i][1];
      coords(i,2)=intpoints.qxg[i][2];
      weights[i]=intpoints.qwgt[i];
    }

    int nnodes = ele.NumNode();
    // get solution matrix with dual parameters
    LINALG::Matrix<8,1>      stdval;

    // establish fundamental data
    double detg = 0.0;

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    LINALG::SerialDenseMatrix me(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix de(nnodes,nnodes,true);
    LINALG::SerialDenseMatrix ae(nnodes,nnodes);

    for (int i=0;i<ngp;++i)
    {
      double gpc[3] = {coords(i,0), coords(i,1) , coords(i,2)};
      volmortar_shape_function_3D(stdval,gpc[0],gpc[1],gpc[2],DRT::Element::hex8);
      detg = Jacobian<DRT::Element::hex8>(gpc,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=weights[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*weights[i]*stdval(j)*detg;
        }
    }

    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    double xi[3] = {xi0,xi1,xi2};
    volmortar_shape_function_3D(stdval,xi[0],xi[1],xi[2],DRT::Element::hex8);

    // evaluate dual shape functions
    for (int i=0;i<nnodes;++i)
      funct(i)=0;

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
        funct(i)+=ae(i,j)*stdval(j);

    break;
  }
  // *********************************************************************
  // 3D dual linear shape functions (tet4)
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case DRT::Element::tet4:
  {
    funct(0) = -5.0 + 3*xi0 + 3*xi1 + 3*xi2;
    funct(1) =  2.0 + 3*xi0 - 3*xi1 - 3*xi2;
    funct(2) =  2.0 - 3*xi0 + 3*xi1 - 3*xi2;
    funct(3) =  2.0 - 3*xi0 - 3*xi1 + 3*xi2;

    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

}//namespace UTILS
}//namespace MORTAR



#endif /* MORTAR_SHAPE_UTILS_H_ */
