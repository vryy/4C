/*!----------------------------------------------------------------------
\file volmortar_shape.H

\brief shape functions for the volmortar framework

\level 1

<pre>
\maintainer Philipp Farah
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
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_mortar/mortar_element.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_inpar/inpar_volmortar.H"

/*----------------------------------------------------------------------*
 | Utils                                                    farah 01/14 |
 *----------------------------------------------------------------------*/
namespace VOLMORTAR
{

namespace UTILS
{
//=====================================================================================
//=====================================================================================
template<class V>
void volmortar_shape_function_3D_deriv(
  V&                                       deriv1,  ///< to be filled with shape function values
  double&                                  r,       ///< xi0 coordinate
  double&                                  s,       ///< xi1 coordinate
  double&                                  t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType         shape    ///< distinguish between mortar shape
  );

template<DRT::Element::DiscretizationType distype>
bool LocalToGlobal(DRT::Element&           ele,     ///< element which is considered
                   const double*           xi,      ///< para. coordinates
                   double*                 globcoord);

template<DRT::Element::DiscretizationType distype>
double Jacobian(
    const double*                          xi,      ///< para. coordinates
    DRT::Element&                          ele);    ///< element which is considered

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
double nurbs_Jacobian(
  W&                                       deriv,   ///< to be filled with shape function deriv
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  );

template<class V>
void volmortar_shape_function_1D(
  V&                                       funct,   ///< to be filled with shape function values
  double&                                  xi0,     ///< xi0 coordinate
  DRT::Element::DiscretizationType         shape    ///< distinguish between mortar shape
  );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_shape_function_1D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function values
  const U*                                 xi,      ///< xi0 coordinate
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  );

template<class V>
void volmortar_shape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  const double&                            xi0,     ///< xi0 coordinate
  const double&                            xi1,     ///< xi1 coordinate
  DRT::Element::DiscretizationType         shape    ///< distinguish between mortar shape
  );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_shape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function values
  const U*                                 xi,      ///< xi0 coordinate
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  );

template<class V>
void volmortar_shape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  const double&                            r,       ///< xi0 coordinate
  const double&                            s,       ///< xi1 coordinate
  const double&                            t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType         shape    ///< distinguish between mortar shape
  );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_shape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function values
  const U*                                 xi,      ///< xi0 coordinate
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  );

template<class V>
void volmortar_shape_function_3D_modified(
  V&                                       funct,   ///< to be filled with shape function values
  const double&                            r,       ///< xi0 coordinate
  const double&                            s,       ///< xi1 coordinate
  const double&                            t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType         shape    ///< distinguish between mortar shape
  );

template<DRT::Element::DiscretizationType distype, class V, class T>
void volmortar_dualshape_function_1D(
    V&                                       funct,   ///< to be filled with shape function values
    DRT::Element&                            ele,     ///< element which is considered
    const T*                                 xi,      ///< para. coordinates
    INPAR::VOLMORTAR::DualQuad               quadtype ///< type of quadratic element modification
    );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_dualshape_function_1D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  );

template<DRT::Element::DiscretizationType distype, class V, class T>
void volmortar_dualshape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  DRT::Element&                            ele,     ///< element which is considered
  const T*                                 xi,      ///< para. coordinates
  INPAR::VOLMORTAR::DualQuad               quadtype ///< type of quadratic element modification
  );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_dualshape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  );

template<DRT::Element::DiscretizationType distype, class V, class T>
void volmortar_dualshape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  DRT::Element&                            ele,     ///< element which is considered
  const T*                                 xi,      ///< para. coordinates
  INPAR::VOLMORTAR::DualQuad               quadtype ///< type of quadratic element modification
  );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_dualshape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  );

// general evaluation routine for std shape functions
template<DRT::Element::DiscretizationType distype, class V, class T>
void shape_function(
    V&                                     f,       ///< to be filled with shape function values
    const T*                               xi,      ///< para. coordinates
    INPAR::VOLMORTAR::DualQuad             dualquad = INPAR::VOLMORTAR::dualquad_no_mod);

// general evaluation routine for dual shape functions
template<DRT::Element::DiscretizationType distype, class V, class T>
void dual_shape_function(
    V&                                     f,       ///< to be filled with shape function values
    const T&                               xi,      ///< para. coordinates
    DRT::Element&                          ele,     ///< element which is considered
    INPAR::VOLMORTAR::DualQuad             dualquad = INPAR::VOLMORTAR::dualquad_no_mod);

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void nurbs_shape_function(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  );

template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void nurbs_dualshape_function(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  );
//=====================================================================================
//=====================================================================================

/*----------------------------------------------------------------------*
 |  get optimal gauss rule for volmortar integration         farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
  struct DisTypeToOptGaussRule {};
  template<> struct DisTypeToOptGaussRule<DRT::Element::hex8>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_8point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_64point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_64point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_4point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::tet10>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_45point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::wedge6>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_wedge_6point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::wedge15> {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_wedge_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::pyramid5>{static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_pyramid_8point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs8>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_8point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs27> {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_4point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_3point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_6point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs4>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_4point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs9>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_3point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs2>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs3>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_3point; };

  /*----------------------------------------------------------------------*
   | NURBS info                                                farah 09/14|
   *----------------------------------------------------------------------*/
  template <DRT::Element::DiscretizationType distype>
  bool IsNurbs()
  {
    bool isnurbs = true;
    switch (distype)
    {
    case DRT::Element::nurbs2:
    case DRT::Element::nurbs3:
    case DRT::Element::nurbs4:
    case DRT::Element::nurbs8:
    case DRT::Element::nurbs9:
    case DRT::Element::nurbs27:
      isnurbs = true;
      break;
    case DRT::Element::line2:
    case DRT::Element::line3:
    case DRT::Element::line4:
    case DRT::Element::line5:
    case DRT::Element::line6:
    case DRT::Element::quad4:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet4:
    case DRT::Element::tet10:
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
    case DRT::Element::pyramid5:
      isnurbs = false;
      break;
    default:
      dserror("Distype unknown!");
    }
    return isnurbs;
  }

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

  switch(distype)
  {
  case DRT::Element::line2:
  case DRT::Element::line3:
  {
    // metrics routine gives local basis vectors
    LINALG::Matrix<ndim,nn> deriv;

    // get shape function values and derivatives at xi
    DRT::UTILS::shape_function_1D_deriv1(deriv,xi[0],distype);

    // build basis vectors gxi and geta
    for (int i=0;i<nn;++i)
    {
      // first local basis vector
      gxi[0] += deriv(0,i)*ele.Nodes()[i]->X()[0];
      gxi[1] += deriv(0,i)*ele.Nodes()[i]->X()[1];
      gxi[2] += deriv(0,i)*ele.Nodes()[i]->X()[2];
    }

    // second local basis vector
    geta[0] = 0.0;
    geta[1] = 0.0;
    geta[2] = 1.0;

    // cross product of gxi and geta
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
    cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
    cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
    jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

    break;
  }
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
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

    break;
  }
  case DRT::Element::tet4:
  {
    LINALG::Matrix<4,4> jacob;
    for (int i=0; i<4; i++)  jacob(0,i)=1;
    for (int row=0;row<3;row++)
      for (int col=0;col<4;col++)
        jacob(row+1,col)= ele.Nodes()[col]->X()[row];

    jac=jacob.Determinant()/6.0;

    break;
  }
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::pyramid5:
  {
    LINALG::Matrix<ndim,nn> derivs;
    const double r = xi[0];
    const double s = xi[1];
    const double t = xi[2];

    DRT::UTILS::shape_function_3D_deriv1(derivs,r,s,t,distype);

    LINALG::Matrix<nn,ndim> xrefe;
    for (int i=0; i<nn; ++i)
    {
      DRT::Node** nodes=ele.Nodes();
      if(!nodes) dserror("Nodes() returned null pointer");

      xrefe(i,0) = nodes[i]->X()[0];
      xrefe(i,1) = nodes[i]->X()[1];
      xrefe(i,2) = nodes[i]->X()[2];
    }

    LINALG::Matrix<ndim,ndim> invJ;
    invJ.Clear();

    invJ.Multiply(derivs,xrefe);
    jac = invJ.Invert();
    if (jac <= 0.0)
      dserror("Element Jacobian mapping %10.5e <= 0.0",jac);

    break;
  }

  default:
      dserror("shape unknown\n");
  }

  return jac;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian determinant for NURBS                  farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
double nurbs_Jacobian(
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  )
{
  //! nn_: number of master element nodes
  static const int nn = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  //! number of space dimensions ("+1" due to considering only interface elements)
  static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

  double jac = 0.0;
  std::vector<double> gxi(3);
  std::vector<double> geta(3);

  switch(distype)
  {
  case DRT::Element::nurbs2:
  case DRT::Element::nurbs3:
  {
    // build basis vectors gxi and geta
    for (int i=0;i<nn;++i)
    {
      // first local basis vector
      gxi[0] += deriv(0,i)*ele.Nodes()[i]->X()[0];
      gxi[1] += deriv(0,i)*ele.Nodes()[i]->X()[1];
      gxi[2] += deriv(0,i)*ele.Nodes()[i]->X()[2];
    }

    // second local basis vector
    geta[0] = 0.0;
    geta[1] = 0.0;
    geta[2] = 1.0;

    // cross product of gxi and geta
    double cross[3] = {0.0, 0.0, 0.0};
    cross[0] = gxi[1]*geta[2]-gxi[2]*geta[1];
    cross[1] = gxi[2]*geta[0]-gxi[0]*geta[2];
    cross[2] = gxi[0]*geta[1]-gxi[1]*geta[0];
    jac = sqrt(cross[0]*cross[0]+cross[1]*cross[1]+cross[2]*cross[2]);

    break;
  }
  case DRT::Element::nurbs9:
  case DRT::Element::nurbs4:
  {
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

    break;
  }
  case DRT::Element::nurbs8:
  case DRT::Element::nurbs27:
  {
    LINALG::Matrix<nn,ndim> xrefe;
    for (int i=0; i<nn; ++i)
    {
      DRT::Node** nodes=ele.Nodes();
      if(!nodes) dserror("Nodes() returned null pointer");

      xrefe(i,0) = nodes[i]->X()[0];
      xrefe(i,1) = nodes[i]->X()[1];
      xrefe(i,2) = nodes[i]->X()[2];
    }

    LINALG::Matrix<ndim,ndim> invJ;
    invJ.Clear();

    invJ.Multiply(deriv,xrefe);
    jac = invJ.Invert();
    if (jac <= 0.0)
      dserror("Element Jacobian mapping %10.5e <= 0.0",jac);

    break;
  }

  default:
      dserror("shape unknown\n");
  }

  return jac;
}

/*----------------------------------------------------------------------*
 |  evaluate std. shape 1D                                   farah 09/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_shape_function_1D(
  V&                                 funct,  ///< to be filled with shape function values
  const double&                            xi0,    ///< xi0 coordinate
  DRT::Element::DiscretizationType   shape   ///< distinguish between mortar shape
  )
{
  switch (shape)
  {
  // *********************************************************************
  // 1D standard linear shape functions (line2)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case DRT::Element::line2:
  {
    funct(0) = 0.5*(1.0-xi0);
    funct(1) = 0.5*(1.0+xi0);
    break;
  }
  // *********************************************************************
  // 1D standard quadratic shape functions (line3)
  // (used for interpolation of displacement field)
  // *********************************************************************
  case DRT::Element::line3:
  {
    funct(0) = 0.5*xi0*(xi0-1.0);
    funct(1) = 0.5*xi0*(xi0+1.0);
    funct(2) = (1.0-xi0)*(1.0+xi0);
    break;
  }
  default:
      dserror("shape unknown\n");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate std. nurbs shape 1D                             farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_shape_function_1D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  )
{
  switch (distype)
  {
  case DRT::Element::nurbs2:
  case DRT::Element::nurbs3:
  {
    DRT::NURBS::UTILS::nurbs_get_1D_funct_deriv
    (funct,
     deriv,
     xi,
     knots,
     weights,
     distype);

    break;
  }
  default:
      dserror("shape unknown\n");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate std. shape 2D                                   farah 01/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_shape_function_2D(
  V&                                 funct,  ///< to be filled with shape function values
  const double&                      xi0,    ///< xi0 coordinate
  const double&                      xi1,    ///< xi1 coordinate
  DRT::Element::DiscretizationType   shape   ///< distinguish between mortar shape
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
    funct(0) = 1.0 - xi0 - xi1;
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
    funct(0) = 0.25*(1.0-xi0)*(1.0-xi1);
    funct(1) = 0.25*(1.0+xi0)*(1.0-xi1);
    funct(2) = 0.25*(1.0+xi0)*(1.0+xi1);
    funct(3) = 0.25*(1.0-xi0)*(1.0+xi1);
    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate std. nurbs shape 2D                             farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_shape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  )
{
  switch (distype)
  {
  case DRT::Element::nurbs4:
  case DRT::Element::nurbs9:
  {
    DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
    (funct,
     deriv,
     xi,
     knots,
     weights,
     distype);

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
  V&                                 funct,   ///< to be filled with shape function values
  const double&                      r,       ///< xi0 coordinate
  const double&                      s,       ///< xi1 coordinate
  const double&                      t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType  shape     ///< distinguish between mortar shape
  )
{
  switch (shape)
  {
  case DRT::Element::hex8:
  {
    const double rp  = 1.0+r;
    const double rm  = 1.0-r;
    const double sp  = 1.0+s;
    const double sm  = 1.0-s;
    const double tp  = 1.0+t;
    const double tm  = 1.0-t;
    const double Q18 = 1.0/8.0;

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
  case DRT::Element::hex20:
  {
    const double Q18 = 1.0/8.0;
    const double rp  = 1.0+r;
    const double rm  = 1.0-r;
    const double sp  = 1.0+s;
    const double sm  = 1.0-s;
    const double tp  = 1.0+t;
    const double tm  = 1.0-t;
    const double rrm = 1.0-r*r;
    const double ssm = 1.0-s*s;
    const double ttm = 1.0-t*t;

    // corner nodes
    funct(0) =Q18*rm*sm*tm*(rm+sm+tm-5.0);
    funct(1) =Q18*rp*sm*tm*(rp+sm+tm-5.0);
    funct(2) =Q18*rp*sp*tm*(rp+sp+tm-5.0);
    funct(3) =Q18*rm*sp*tm*(rm+sp+tm-5.0);
    funct(4) =Q18*rm*sm*tp*(rm+sm+tp-5.0);
    funct(5) =Q18*rp*sm*tp*(rp+sm+tp-5.0);
    funct(6) =Q18*rp*sp*tp*(rp+sp+tp-5.0);
    funct(7) =Q18*rm*sp*tp*(rm+sp+tp-5.0);

    // centernodes, bottom surface
    funct( 8)=0.25*rrm*sm*tm;
    funct( 9)=0.25*rp*ssm*tm;
    funct(10)=0.25*rrm*sp*tm;
    funct(11)=0.25*rm*ssm*tm;

    // centernodes, rs-plane
    funct(12)=0.25*rm*sm*ttm;
    funct(13)=0.25*rp*sm*ttm;
    funct(14)=0.25*rp*sp*ttm;
    funct(15)=0.25*rm*sp*ttm;

    // centernodes, top surface
    funct(16)=0.25*rrm*sm*tp;
    funct(17)=0.25*rp*ssm*tp;
    funct(18)=0.25*rrm*sp*tp;
    funct(19)=0.25*rm*ssm*tp;

    break;
  }
  case DRT::Element::hex27:
  {
    const double rm1 = 0.5*r*(r - 1.0);
    const double r00 = (1.0 - r*r);
    const double rp1 = 0.5*r*(r + 1.0);
    const double sm1 = 0.5*s*(s - 1.0);
    const double s00 = (1.0 - s*s);
    const double sp1 = 0.5*s*(s + 1.0);
    const double tm1 = 0.5*t*(t - 1.0);
    const double t00 = (1.0 - t*t);
    const double tp1 = 0.5*t*(t + 1.0);

    funct( 0) = rm1*sm1*tm1;
    funct( 1) = rp1*sm1*tm1;
    funct( 2) = rp1*sp1*tm1;
    funct( 3) = rm1*sp1*tm1;
    funct( 4) = rm1*sm1*tp1;
    funct( 5) = rp1*sm1*tp1;
    funct( 6) = rp1*sp1*tp1;
    funct( 7) = rm1*sp1*tp1;
    funct( 8) = r00*sm1*tm1;
    funct( 9) = s00*tm1*rp1;
    funct(10) = r00*tm1*sp1;
    funct(11) = s00*rm1*tm1;
    funct(12) = t00*rm1*sm1;
    funct(13) = t00*sm1*rp1;
    funct(14) = t00*rp1*sp1;
    funct(15) = t00*rm1*sp1;
    funct(16) = r00*sm1*tp1;
    funct(17) = s00*rp1*tp1;
    funct(18) = r00*sp1*tp1;
    funct(19) = s00*rm1*tp1;
    funct(20) = r00*s00*tm1;
    funct(21) = r00*t00*sm1;
    funct(22) = s00*t00*rp1;
    funct(23) = r00*t00*sp1;
    funct(24) = s00*t00*rm1;
    funct(25) = r00*s00*tp1;
    funct(26) = r00*s00*t00;
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
  case DRT::Element::tet10:
  {
    const double u=1.0-r-s-t;

    funct(0) =u*(2*u -1.0);
    funct(1) =r*(2*r -1.0);
    funct(2) =s*(2*s -1.0);
    funct(3) =t*(2*t -1.0);
    funct(4) =4*r*u;
    funct(5) =4*r*s;
    funct(6) =4*s*u;
    funct(7) =4*t*u;
    funct(8) =4*r*t;
    funct(9) =4*s*t;
    break;
  }
  case DRT::Element::pyramid5:
  {
    DRT::UTILS::shape_function_3D(funct,r,s,t,shape);

    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate std. nurbs shape 3D                             farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_shape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  )
{
  switch (distype)
  {
  case DRT::Element::nurbs8:
  case DRT::Element::nurbs27:
  {
    DRT::NURBS::UTILS::nurbs_get_3D_funct_deriv
    (funct,
     deriv,
     xi,
     knots,
     weights,
     distype);

    break;
  }
  default:
      dserror("shape unknown\n");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate modified shape 3D (for quadr. dual elements)    farah 05/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_shape_function_3D_modified(
  V&                                 funct,   ///< to be filled with shape function values
  const double&                      r,       ///< xi0 coordinate
  const double&                      s,       ///< xi1 coordinate
  const double&                      t,       ///< xi2 coordinate
  DRT::Element::DiscretizationType  shape     ///< distinguish between mortar shape
  )
{
  switch (shape)
  {
  case DRT::Element::hex20:
  {
    const double Q18 = 1.0/8.0;
    const double rp  = 1.0+r;
    const double rm  = 1.0-r;
    const double sp  = 1.0+s;
    const double sm  = 1.0-s;
    const double tp  = 1.0+t;
    const double tm  = 1.0-t;
    const double rrm = 1.0-r*r;
    const double ssm = 1.0-s*s;
    const double ttm = 1.0-t*t;

    LINALG::Matrix<20,1> valtmp;

    // corner nodes
    valtmp(0) =Q18*rm*sm*tm*(rm+sm+tm-5.0);
    valtmp(1) =Q18*rp*sm*tm*(rp+sm+tm-5.0);
    valtmp(2) =Q18*rp*sp*tm*(rp+sp+tm-5.0);
    valtmp(3) =Q18*rm*sp*tm*(rm+sp+tm-5.0);
    valtmp(4) =Q18*rm*sm*tp*(rm+sm+tp-5.0);
    valtmp(5) =Q18*rp*sm*tp*(rp+sm+tp-5.0);
    valtmp(6) =Q18*rp*sp*tp*(rp+sp+tp-5.0);
    valtmp(7) =Q18*rm*sp*tp*(rm+sp+tp-5.0);

    // centernodes, bottom surface
    valtmp( 8)=0.25*rrm*sm*tm;
    valtmp( 9)=0.25*rp*ssm*tm;
    valtmp(10)=0.25*rrm*sp*tm;
    valtmp(11)=0.25*rm*ssm*tm;

    // centernodes, rs-plane
    valtmp(12)=0.25*rm*sm*ttm;
    valtmp(13)=0.25*rp*sm*ttm;
    valtmp(14)=0.25*rp*sp*ttm;
    valtmp(15)=0.25*rm*sp*ttm;

    // centernodes, top surface
    valtmp(16)=0.25*rrm*sm*tp;
    valtmp(17)=0.25*rp*ssm*tp;
    valtmp(18)=0.25*rrm*sp*tp;
    valtmp(19)=0.25*rm*ssm*tp;

    // *******************************
    // Basis-Trafo:

    const double alpha = 0.3;

    // corner nodes
    funct(0) = valtmp(0) + alpha * (valtmp(8)  + valtmp(11) + valtmp(12));
    funct(1) = valtmp(1) + alpha * (valtmp(8)  + valtmp(9)  + valtmp(13));
    funct(2) = valtmp(2) + alpha * (valtmp(9)  + valtmp(10) + valtmp(14));
    funct(3) = valtmp(3) + alpha * (valtmp(10) + valtmp(11) + valtmp(15));
    funct(4) = valtmp(4) + alpha * (valtmp(12) + valtmp(16) + valtmp(19));
    funct(5) = valtmp(5) + alpha * (valtmp(13) + valtmp(16) + valtmp(17));
    funct(6) = valtmp(6) + alpha * (valtmp(14) + valtmp(17) + valtmp(18));
    funct(7) = valtmp(7) + alpha * (valtmp(15) + valtmp(18) + valtmp(19));

    //edge nodes
    funct(8)  = valtmp(8)  * (1.0-2.0*alpha);
    funct(9)  = valtmp(9)  * (1.0-2.0*alpha);
    funct(10) = valtmp(10) * (1.0-2.0*alpha);
    funct(11) = valtmp(11) * (1.0-2.0*alpha);
    funct(12) = valtmp(12) * (1.0-2.0*alpha);
    funct(13) = valtmp(13) * (1.0-2.0*alpha);
    funct(14) = valtmp(14) * (1.0-2.0*alpha);
    funct(15) = valtmp(15) * (1.0-2.0*alpha);
    funct(16) = valtmp(16) * (1.0-2.0*alpha);
    funct(17) = valtmp(17) * (1.0-2.0*alpha);
    funct(18) = valtmp(18) * (1.0-2.0*alpha);
    funct(19) = valtmp(19) * (1.0-2.0*alpha);

    break;
  }
  case DRT::Element::tet10:
  {
    const double u=1.0-r-s-t;
    LINALG::Matrix<10,1> valtmp;

    valtmp(0) = u*(2.0*u -1.0);
    valtmp(1) = r*(2.0*r -1.0);
    valtmp(2) = s*(2.0*s -1.0);
    valtmp(3) = t*(2.0*t -1.0);
    valtmp(4) = 4.0*r*u;
    valtmp(5) = 4.0*r*s;
    valtmp(6) = 4.0*s*u;
    valtmp(7) = 4.0*t*u;
    valtmp(8) = 4.0*r*t;
    valtmp(9) = 4.0*s*t;

    // *******************************
    // Basis-Trafo:

    const double alpha = 0.3;

    // corner nodes
    funct(0) = valtmp(0) + alpha * (valtmp(4)  + valtmp(6) + valtmp(7));
    funct(1) = valtmp(1) + alpha * (valtmp(4)  + valtmp(5) + valtmp(8));
    funct(2) = valtmp(2) + alpha * (valtmp(5)  + valtmp(6) + valtmp(9));
    funct(3) = valtmp(3) + alpha * (valtmp(7)  + valtmp(8) + valtmp(9));

    // edge nodes
    funct(4)  = valtmp(4)  * (1.0-2.0*alpha);
    funct(5)  = valtmp(5)  * (1.0-2.0*alpha);
    funct(6)  = valtmp(6)  * (1.0-2.0*alpha);
    funct(7)  = valtmp(7)  * (1.0-2.0*alpha);
    funct(8)  = valtmp(8)  * (1.0-2.0*alpha);
    funct(9)  = valtmp(9)  * (1.0-2.0*alpha);

    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate dual shape 1D                                   farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class T>
void volmortar_dualshape_function_1D(
    V&                                       funct,   ///< to be filled with shape function values
    DRT::Element&                            ele,     ///< element which is considered
    const T*                                 xi,      ///< para. coordinates
    INPAR::VOLMORTAR::DualQuad               quadtype ///< type of quadratic element modification
    )
{
  switch (distype)
  {
  // *********************************************************************
  // 1D dual linear shape functions (line2)
  // (used for interpolation of Lagrange multiplier field)
  // *********************************************************************
  case DRT::Element::line2:
  {
    funct(0) = 0.5*(1-3.0*xi[0]);
    funct(1) = 0.5*(1+3.0*xi[0]);

    break;
  }
  // *********************************************************************
  // *********************************************************************
  case DRT::Element::line3:
  {
    // establish fundamental data
    //! number of space dimensions
    static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of element nodes
    static const int nnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get gauss rule
    const DRT::UTILS::IntPointsAndWeights<ndim> intpoints(DisTypeToOptGaussRule<distype>::rule);

    double detg = 0.0;

    LINALG::Matrix<nnodes,nnodes> me(true);
    LINALG::Matrix<nnodes,nnodes> de(true);
    LINALG::Matrix<nnodes,nnodes> ae;

    LINALG::Matrix<nnodes,1> stdval;

    for(int i=0;i<nnodes;++i)
      funct(i)=0.0;

    for (int i=0;i<intpoints.IP().nquad;++i)
    {
      double gpc[1] = {intpoints.IP().qxg[i][0]};

      shape_function<distype>(stdval,gpc);
      detg = Jacobian<distype>(gpc,ele);

      for (int j=0;j<nnodes;++j)
      {
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=intpoints.IP().qwgt[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*intpoints.IP().qwgt[i]*stdval(j)*detg;
        }
      }
    }

    // get solution matrix with dual parameters
    LINALG::InvertAndMultiplyByCholesky<nnodes>(me,de,ae);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    volmortar_shape_function_1D(stdval,xi[0],DRT::Element::line3);

    // evaluate dual shape functions
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
 |  evaluate dual NURBS shapes 1D                            farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_dualshape_function_1D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  )
{
  switch (distype)
  {
  // *********************************************************************
  // *********************************************************************
  case DRT::Element::nurbs2:
  case DRT::Element::nurbs3:
  {
    //! number of space dimensions
    static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of element nodes
    static const int nnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get gauss rule
    const DRT::UTILS::IntPointsAndWeights<ndim> intpoints(DisTypeToOptGaussRule<distype>::rule);

    // get solution matrix with dual parameters
    LINALG::Matrix<nnodes,1> stdval;
    W refderiv = deriv;

    // establish fundamental data
    double detg = 0.0;

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    LINALG::Matrix<nnodes,nnodes> me(true);
    LINALG::Matrix<nnodes,nnodes> de(true);
    LINALG::Matrix<nnodes,nnodes> ae;

    for (int i=0;i<intpoints.IP().nquad;++i)
    {
      double gpc[1] = {intpoints.IP().qxg[i][0]};
      nurbs_shape_function<distype>(stdval,refderiv,gpc,weights,knots);
      detg = nurbs_Jacobian<distype>(refderiv,gpc,weights,knots,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=intpoints.IP().qwgt[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*intpoints.IP().qwgt[i]*stdval(j)*detg;
        }
    }

    // calcute coefficient matrix
    LINALG::InvertAndMultiplyByCholesky<nnodes>(me,de,ae);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    nurbs_shape_function<distype>(stdval,refderiv,xi,weights,knots);

    // evaluate dual shape functions
    for (int i=0;i<nnodes;++i)
    {
      funct(i)   = 0.0;
      deriv(0,i) = 0.0;
    }

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        funct(i)   += ae(i,j)*stdval(j);
        deriv(0,i) += ae(i,j)*refderiv(0,j);
      }

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
template<DRT::Element::DiscretizationType distype, class V, class T>
void volmortar_dualshape_function_2D(
    V&                                       funct,   ///< to be filled with shape function values
    DRT::Element&                            ele,     ///< element which is considered
    const T*                                 xi,      ///< para. coordinates
    INPAR::VOLMORTAR::DualQuad               quadtype ///< type of quadratic element modification
    )
{
  switch (distype)
  {
  // *********************************************************************
  // 2D dual linear shape functions (tri3)
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case DRT::Element::tri3:
  {
    funct(0) = 3.0 - 4.0*xi[0] - 4.0*xi[1];
    funct(1) = 4.0 * xi[0]     - 1.0;
    funct(2) = 4.0 * xi[1]     - 1.0;

    break;
  }
  // *********************************************************************
  // *********************************************************************
  case DRT::Element::tri6:
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  case DRT::Element::nurbs9:
  {
    //! number of space dimensions
    static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of element nodes
    static const int nnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get gauss rule
    const DRT::UTILS::IntPointsAndWeights<ndim> intpoints(DisTypeToOptGaussRule<distype>::rule);

    // get solution matrix with dual parameters
    LINALG::Matrix<nnodes,1>      stdval;

    // establish fundamental data
    double detg = 0.0;
    LINALG::Matrix<nnodes,nnodes> me(true);
    LINALG::Matrix<nnodes,nnodes> de(true);
    LINALG::Matrix<nnodes,nnodes> ae;

    for (int i=0;i<intpoints.IP().nquad;++i)
    {
      double gpc[2] = {intpoints.IP().qxg[i][0], intpoints.IP().qxg[i][1]};
      shape_function<distype>(stdval,gpc);

      detg = Jacobian<distype>(gpc,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=intpoints.IP().qwgt[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*intpoints.IP().qwgt[i]*stdval(j)*detg;
        }
    }
    // invert bi-ortho matrix me
    LINALG::InvertAndMultiplyByCholesky<nnodes>(me,de,ae);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    shape_function<distype>(stdval,xi);

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
 |  evaluate dual NURBS shapes 2D                            farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_dualshape_function_2D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  )
{
  switch (distype)
  {
  // *********************************************************************
  // *********************************************************************
  case DRT::Element::nurbs4:
  case DRT::Element::nurbs9:
  {
    //! number of space dimensions
    static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of element nodes
    static const int nnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get gauss rule
    const DRT::UTILS::IntPointsAndWeights<ndim> intpoints(DisTypeToOptGaussRule<distype>::rule);

    // get solution matrix with dual parameters
    LINALG::Matrix<nnodes,1> stdval;
    W refderiv = deriv;

    // establish fundamental data
    double detg = 0.0;

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    LINALG::Matrix<nnodes,nnodes> me(true);
    LINALG::Matrix<nnodes,nnodes> de(true);
    LINALG::Matrix<nnodes,nnodes> ae;

    for (int i=0;i<intpoints.IP().nquad;++i)
    {
      double gpc[2] = {intpoints.IP().qxg[i][0], intpoints.IP().qxg[i][1]};
      nurbs_shape_function<distype>(stdval,refderiv,gpc,weights,knots);
      detg = nurbs_Jacobian<distype>(refderiv,gpc,weights,knots,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=intpoints.IP().qwgt[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*intpoints.IP().qwgt[i]*stdval(j)*detg;
        }
    }

    // calcute coefficient matrix
    LINALG::InvertAndMultiplyByCholesky<nnodes>(me,de,ae);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    nurbs_shape_function<distype>(stdval,refderiv,xi,weights,knots);

    // evaluate dual shape functions
    for (int i=0;i<nnodes;++i)
    {
      funct(i)   = 0.0;
      deriv(0,i) = 0.0;
      deriv(1,i) = 0.0;
    }

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        funct(i)   += ae(i,j)*stdval(j);
        deriv(0,i) += ae(i,j)*refderiv(0,j);
        deriv(1,i) += ae(i,j)*refderiv(1,j);
      }

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
template<DRT::Element::DiscretizationType distype, class V, class T>
void volmortar_dualshape_function_3D(
    V&                                       funct,   ///< to be filled with shape function values
    DRT::Element&                            ele,     ///< element which is considered
    const T*                                 xi,      ///< para. coordinates
    INPAR::VOLMORTAR::DualQuad               quadtype ///< type of quadratic element modification
    )
{
  switch (distype)
  {
  // *********************************************************************
  // *********************************************************************
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  case DRT::Element::tet10:
  case DRT::Element::pyramid5:
  {
    //! number of space dimensions
    static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of element nodes
    static const int nnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get gauss rule
    const DRT::UTILS::IntPointsAndWeights<ndim> intpoints(DisTypeToOptGaussRule<distype>::rule);

    // get solution matrix with dual parameters
    LINALG::Matrix<nnodes,1>      stdval;

    // establish fundamental data
    double detg = 0.0;

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    LINALG::Matrix<nnodes,nnodes> me(true);
    LINALG::Matrix<nnodes,nnodes> de(true);
    LINALG::Matrix<nnodes,nnodes> ae;

    for (int i=0;i<intpoints.IP().nquad;++i)
    {
      double gpc[3] = {intpoints.IP().qxg[i][0], intpoints.IP().qxg[i][1] , intpoints.IP().qxg[i][2]};
      shape_function<distype>(stdval,gpc,quadtype);
      detg = Jacobian<distype>(gpc,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=intpoints.IP().qwgt[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*intpoints.IP().qwgt[i]*stdval(j)*detg;
        }
    }

    // calcute coefficient matrix
    LINALG::InvertAndMultiplyByCholesky<nnodes>(me,de,ae);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    shape_function<distype>(stdval,xi,quadtype);

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
  // *********************************************************************
  case DRT::Element::tet4:
  {
    funct(0) =  4.0 - 5.0*xi[0] - 5.0*xi[1] - 5.0*xi[2];
    funct(1) = -1.0 + 5.0*xi[0];
    funct(2) = -1.0 + 5.0*xi[1];
    funct(3) = -1.0 + 5.0*xi[2];

    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate dual NURBS shapes 3D                            farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void volmortar_nurbs_dualshape_function_3D(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  )
{
  switch (distype)
  {
  // *********************************************************************
  // *********************************************************************
  case DRT::Element::nurbs8:
  case DRT::Element::nurbs27:
  {
    //! number of space dimensions
    static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

    //! number of element nodes
    static const int nnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get gauss rule
    const DRT::UTILS::IntPointsAndWeights<ndim> intpoints(DisTypeToOptGaussRule<distype>::rule);

    // get solution matrix with dual parameters
    LINALG::Matrix<nnodes,1> stdval;
    W refderiv = deriv;

    // establish fundamental data
    double detg = 0.0;

    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    LINALG::Matrix<nnodes,nnodes> me(true);
    LINALG::Matrix<nnodes,nnodes> de(true);
    LINALG::Matrix<nnodes,nnodes> ae;

    for (int i=0;i<intpoints.IP().nquad;++i)
    {
      double gpc[3] = {intpoints.IP().qxg[i][0], intpoints.IP().qxg[i][1] , intpoints.IP().qxg[i][2]};
      nurbs_shape_function<distype>(stdval,refderiv,gpc,weights,knots);
      detg = nurbs_Jacobian<distype>(refderiv,gpc,weights,knots,ele);

      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=intpoints.IP().qwgt[i]*stdval(j)*stdval(k)*detg;
          de(j,k)+=(j==k)*intpoints.IP().qwgt[i]*stdval(j)*detg;
        }
    }

    // calcute coefficient matrix
    LINALG::InvertAndMultiplyByCholesky<nnodes>(me,de,ae);

    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    nurbs_shape_function<distype>(stdval,refderiv,xi,weights,knots);

    // evaluate dual shape functions
    for (int i=0;i<nnodes;++i)
    {
      funct(i)   = 0.0;
      deriv(0,i) = 0.0;
      deriv(1,i) = 0.0;
      deriv(2,i) = 0.0;
    }

    for (int i=0;i<nnodes;++i)
      for (int j=0;j<nnodes;++j)
      {
        funct(i)   += ae(i,j)*stdval(j);
        deriv(0,i) += ae(i,j)*refderiv(0,j);
        deriv(1,i) += ae(i,j)*refderiv(1,j);
        deriv(2,i) += ae(i,j)*refderiv(2,j);
      }

    break;
  }
  default:
      dserror("shape unknown\n");
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate modified deriv 3D                               farah 05/14|
 *----------------------------------------------------------------------*/
template<class V>
void volmortar_shape_function_3D_deriv(
  V&                                 deriv1,   ///< to be filled with shape function values
  double&                            r,        ///< xi0 coordinate
  double&                            s,        ///< xi1 coordinate
  double&                            t,        ///< xi2 coordinate
  DRT::Element::DiscretizationType  shape      ///< distinguish between mortar shape
  )
{
  switch (shape)
  {
  case DRT::Element::hex20:
  {
    // form basic values
    const double rp  = 1.0+r;
    const double rm  = 1.0-r;
    const double sp  = 1.0+s;
    const double sm  = 1.0-s;
    const double tp  = 1.0+t;
    const double tm  = 1.0-t;
    const double rrm = 1.0-r*r;
    const double ssm = 1.0-s*s;
    const double ttm = 1.0-t*t;
    const double Q18 = 1.0/8.0;

    LINALG::Matrix<3,20> valtmp;

    // corner nodes
    valtmp(0, 0) =-Q18*sm*tm*(2.0*rm+sm+tm-5.0);
    valtmp(1, 0) =-Q18*tm*rm*(2.0*sm+tm+rm-5.0);
    valtmp(2, 0) =-Q18*rm*sm*(2.0*tm+rm+sm-5.0);

    valtmp(0, 1) = Q18*sm*tm*(2.0*rp+sm+tm-5.0);
    valtmp(1, 1) =-Q18*tm*rp*(2.0*sm+tm+rp-5.0);
    valtmp(2, 1) =-Q18*rp*sm*(2.0*tm+rp+sm-5.0);

    valtmp(0, 2) = Q18*sp*tm*(2.0*rp+sp+tm-5.0);
    valtmp(1, 2) = Q18*tm*rp*(2.0*sp+tm+rp-5.0);
    valtmp(2, 2) =-Q18*rp*sp*(2.0*tm+rp+sp-5.0);

    valtmp(0, 3) =-Q18*sp*tm*(2.0*rm+sp+tm-5.0);
    valtmp(1, 3) = Q18*tm*rm*(2.0*sp+tm+rm-5.0);
    valtmp(2, 3) =-Q18*rm*sp*(2.0*tm+rm+sp-5.0);

    valtmp(0, 4) =-Q18*sm*tp*(2.0*rm+sm+tp-5.0);
    valtmp(1, 4) =-Q18*tp*rm*(2.0*sm+tp+rm-5.0);
    valtmp(2, 4) = Q18*rm*sm*(2.0*tp+rm+sm-5.0);

    valtmp(0, 5) = Q18*sm*tp*(2.0*rp+sm+tp-5.0);
    valtmp(1, 5) =-Q18*tp*rp*(2.0*sm+tp+rp-5.0);
    valtmp(2, 5) = Q18*rp*sm*(2.0*tp+rp+sm-5.0);

    valtmp(0, 6) = Q18*sp*tp*(2.0*rp+sp+tp-5.0);
    valtmp(1, 6) = Q18*tp*rp*(2.0*sp+tp+rp-5.0);
    valtmp(2, 6) = Q18*rp*sp*(2.0*tp+rp+sp-5.0);

    valtmp(0, 7) =-Q18*sp*tp*(2.0*rm+sp+tp-5.0);
    valtmp(1, 7) = Q18*tp*rm*(2.0*sp+tp+rm-5.0);
    valtmp(2, 7) = Q18*rm*sp*(2.0*tp+rm+sp-5.0);

    // centernodes, bottom surface
    valtmp(0, 8) =-0.5*r*sm*tm;
    valtmp(1, 8) =-0.25*rrm*tm;
    valtmp(2, 8) =-0.25*rrm*sm;

    valtmp(0, 9) = 0.25*ssm*tm;
    valtmp(1, 9) =-0.5*s*tm*rp;
    valtmp(2, 9) =-0.25*ssm*rp;

    valtmp(0, 10) =-0.5*r*sp*tm;
    valtmp(1, 10) = 0.25*rrm*tm;
    valtmp(2, 10) =-0.25*rrm*sp;

    valtmp(0, 11)=-0.25*ssm*tm;
    valtmp(1, 11)=-0.5*s*tm*rm;
    valtmp(2, 11)=-0.25*ssm*rm;

    // centernodes, rs-plane
    valtmp(0, 12)=-0.25*sm*ttm;
    valtmp(1, 12)=-0.25*ttm*rm;
    valtmp(2, 12)=-0.5*t*rm*sm;

    valtmp(0, 13)= 0.25*sm*ttm;
    valtmp(1, 13)=-0.25*ttm*rp;
    valtmp(2, 13)=-0.5*t*rp*sm;

    valtmp(0, 14)= 0.25*sp*ttm;
    valtmp(1, 14)= 0.25*ttm*rp;
    valtmp(2, 14)=-0.5*t*rp*sp;

    valtmp(0, 15)=-0.25*sp*ttm;
    valtmp(1, 15)= 0.25*ttm*rm;
    valtmp(2, 15)=-0.5*t*rm*sp;

    // centernodes, top surface
    valtmp(0, 16)=-0.5*r*sm*tp;
    valtmp(1, 16)=-0.25*rrm*tp;
    valtmp(2, 16)= 0.25*rrm*sm;

    valtmp(0, 17)= 0.25*ssm*tp;
    valtmp(1, 17)=-0.5*s*tp*rp;
    valtmp(2, 17)= 0.25*ssm*rp;

    valtmp(0, 18)=-0.5*r*sp*tp;
    valtmp(1, 18)= 0.25*rrm*tp;
    valtmp(2, 18)= 0.25*rrm*sp;

    valtmp(0, 19)=-0.25*ssm*tp;
    valtmp(1, 19)=-0.5*s*tp*rm;
    valtmp(2, 19)= 0.25*ssm*rm;

    // *******************************
    // Basis-Trafo:

    const double alpha = 0.1;

    // corner nodes
    deriv1(0,0) = valtmp(0,0) + alpha * (valtmp(0,8)  + valtmp(0,11) + valtmp(0,12));
    deriv1(1,0) = valtmp(1,0) + alpha * (valtmp(1,8)  + valtmp(1,11) + valtmp(1,12));
    deriv1(2,0) = valtmp(2,0) + alpha * (valtmp(2,8)  + valtmp(2,11) + valtmp(2,12));

    deriv1(0,1) = valtmp(0,1) + alpha * (valtmp(0,8)  + valtmp(0,9)  + valtmp(0,13));
    deriv1(1,1) = valtmp(1,1) + alpha * (valtmp(1,8)  + valtmp(1,9)  + valtmp(1,13));
    deriv1(2,1) = valtmp(2,1) + alpha * (valtmp(2,8)  + valtmp(2,9)  + valtmp(2,13));

    deriv1(0,2) = valtmp(0,2) + alpha * (valtmp(0,9)  + valtmp(0,10) + valtmp(0,14));
    deriv1(1,2) = valtmp(1,2) + alpha * (valtmp(1,9)  + valtmp(1,10) + valtmp(1,14));
    deriv1(2,2) = valtmp(2,2) + alpha * (valtmp(2,9)  + valtmp(2,10) + valtmp(2,14));

    deriv1(0,3) = valtmp(0,3) + alpha * (valtmp(0,10) + valtmp(0,11) + valtmp(0,15));
    deriv1(1,3) = valtmp(1,3) + alpha * (valtmp(1,10) + valtmp(1,11) + valtmp(1,15));
    deriv1(2,3) = valtmp(2,3) + alpha * (valtmp(2,10) + valtmp(2,11) + valtmp(2,15));

    deriv1(0,4) = valtmp(0,4) + alpha * (valtmp(0,12) + valtmp(0,16) + valtmp(0,19));
    deriv1(1,4) = valtmp(1,4) + alpha * (valtmp(1,12) + valtmp(1,16) + valtmp(1,19));
    deriv1(2,4) = valtmp(2,4) + alpha * (valtmp(2,12) + valtmp(2,16) + valtmp(2,19));

    deriv1(0,5) = valtmp(0,5) + alpha * (valtmp(0,13) + valtmp(0,16) + valtmp(0,17));
    deriv1(1,5) = valtmp(1,5) + alpha * (valtmp(1,13) + valtmp(1,16) + valtmp(1,17));
    deriv1(2,5) = valtmp(2,5) + alpha * (valtmp(2,13) + valtmp(2,16) + valtmp(2,17));

    deriv1(0,6) = valtmp(0,6) + alpha * (valtmp(0,14) + valtmp(0,17) + valtmp(0,18));
    deriv1(1,6) = valtmp(1,6) + alpha * (valtmp(1,14) + valtmp(1,17) + valtmp(1,18));
    deriv1(2,6) = valtmp(2,6) + alpha * (valtmp(2,14) + valtmp(2,17) + valtmp(2,18));

    deriv1(0,7) = valtmp(0,7) + alpha * (valtmp(0,15) + valtmp(0,18) + valtmp(0,19));
    deriv1(1,7) = valtmp(1,7) + alpha * (valtmp(1,15) + valtmp(1,18) + valtmp(1,19));
    deriv1(2,7) = valtmp(2,7) + alpha * (valtmp(2,15) + valtmp(2,18) + valtmp(2,19));

    //edge nodes
    deriv1(0,8)  = valtmp(0,8)  * (1.0-3.0*alpha);
    deriv1(1,8)  = valtmp(1,8)  * (1.0-3.0*alpha);
    deriv1(2,8)  = valtmp(2,8)  * (1.0-3.0*alpha);

    deriv1(0,9)  = valtmp(0,9)  * (1.0-3.0*alpha);
    deriv1(1,9)  = valtmp(1,9)  * (1.0-3.0*alpha);
    deriv1(2,9)  = valtmp(2,9)  * (1.0-3.0*alpha);

    deriv1(0,10) = valtmp(0,10) * (1.0-3.0*alpha);
    deriv1(1,10) = valtmp(1,10) * (1.0-3.0*alpha);
    deriv1(2,10) = valtmp(2,10) * (1.0-3.0*alpha);

    deriv1(0,11) = valtmp(0,11) * (1.0-3.0*alpha);
    deriv1(1,11) = valtmp(1,11) * (1.0-3.0*alpha);
    deriv1(2,11) = valtmp(2,11) * (1.0-3.0*alpha);

    deriv1(0,12) = valtmp(0,12) * (1.0-3.0*alpha);
    deriv1(1,12) = valtmp(1,12) * (1.0-3.0*alpha);
    deriv1(2,12) = valtmp(2,12) * (1.0-3.0*alpha);

    deriv1(0,13) = valtmp(0,13) * (1.0-3.0*alpha);
    deriv1(1,13) = valtmp(1,13) * (1.0-3.0*alpha);
    deriv1(2,13) = valtmp(2,13) * (1.0-3.0*alpha);

    deriv1(0,14) = valtmp(0,14) * (1.0-3.0*alpha);
    deriv1(1,14) = valtmp(1,14) * (1.0-3.0*alpha);
    deriv1(2,14) = valtmp(2,14) * (1.0-3.0*alpha);

    deriv1(0,15) = valtmp(0,15) * (1.0-3.0*alpha);
    deriv1(1,15) = valtmp(1,15) * (1.0-3.0*alpha);
    deriv1(2,15) = valtmp(2,15) * (1.0-3.0*alpha);

    deriv1(0,16) = valtmp(0,16) * (1.0-3.0*alpha);
    deriv1(1,16) = valtmp(1,16) * (1.0-3.0*alpha);
    deriv1(2,16) = valtmp(2,16) * (1.0-3.0*alpha);

    deriv1(0,17) = valtmp(0,17) * (1.0-3.0*alpha);
    deriv1(1,17) = valtmp(1,17) * (1.0-3.0*alpha);
    deriv1(2,17) = valtmp(2,17) * (1.0-3.0*alpha);

    deriv1(0,18) = valtmp(0,18) * (1.0-3.0*alpha);
    deriv1(1,18) = valtmp(1,18) * (1.0-3.0*alpha);
    deriv1(2,18) = valtmp(2,18) * (1.0-3.0*alpha);

    deriv1(0,19) = valtmp(0,19) * (1.0-3.0*alpha);
    deriv1(1,19) = valtmp(1,19) * (1.0-3.0*alpha);
    deriv1(2,19) = valtmp(2,19) * (1.0-3.0*alpha);

    break;
  }
  default:
      dserror("shape unknown\n");
  }

  return;
}


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
  if (IsNurbs<distype>()) dserror("ERROR: Lagr. LocalToGlobal called for NURBS!");

  static const int n    = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
  static const int ndim = DRT::UTILS::DisTypeToDim<distype>::dim;

  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("ERROR: LocalToGlobal: Null pointer!");

  for (int i=0;i<ndim;++i)
    globcoord[i]=0.0;

  LINALG::Matrix<n,1>    val;
  LINALG::Matrix<ndim,n> coord;

  shape_function<distype>(val,xi);

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
 |  evaluate shapes                                          farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class T>
void shape_function(
    V&                             f,
    const T*                       xi,
    INPAR::VOLMORTAR::DualQuad     dualquad)
{
  switch (DRT::UTILS::DisTypeToDim<distype>::dim)
  {
  case 1:
  {
    volmortar_shape_function_1D(f, xi[0], distype);
    break;
  }
  case 2:
  {
    volmortar_shape_function_2D(f, xi[0],xi[1], distype);
    break;
  }
  case 3:
  {
    // modified shape function
    if(dualquad == INPAR::VOLMORTAR::dualquad_quad_mod)
      volmortar_shape_function_3D_modified(f, xi[0],xi[1], xi[2],distype);

    // non-modified shape function
    else if (dualquad == INPAR::VOLMORTAR::dualquad_no_mod)
      volmortar_shape_function_3D(f, xi[0],xi[1], xi[2],distype);

    // not implemented
    else
      dserror("No lin modification for quadratic elements possible!");

    break;
  }
  default:
    dserror("dimension of the element is not correct");
    break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate NURBS shapes                                    farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void nurbs_shape_function(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots    ///< knot vectors
  )
{
  switch (DRT::UTILS::DisTypeToDim<distype>::dim)
  {
  case 1:
  {
    volmortar_nurbs_shape_function_1D<distype>(funct,deriv,xi,weights,knots);
    break;
  }
  case 2:
  {
    volmortar_nurbs_shape_function_2D<distype>(funct,deriv,xi,weights,knots);
    break;
  }
  case 3:
  {
    volmortar_nurbs_shape_function_3D<distype>(funct,deriv,xi,weights,knots);
    break;
  }
  default:
    dserror("dimension of the element is not correct");
    break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate dual shapes                                     farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class T>
void dual_shape_function(
    V&                             f,
    const T&                       xi,
    DRT::Element&                  ele,
    INPAR::VOLMORTAR::DualQuad     dualquad)
{
  switch (DRT::UTILS::DisTypeToDim<distype>::dim)
  {
  case 1:
  {
    volmortar_dualshape_function_1D<distype>(f,ele, xi,dualquad);
    break;
  }
  case 2:
  {
    volmortar_dualshape_function_2D<distype>(f,ele, xi,dualquad);
    break;
  }
  case 3:
  {
    volmortar_dualshape_function_3D<distype>(f,ele, xi,dualquad);
    break;
  }
  default:
    dserror("dimension of the element is not correct");
    break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate dual NURBS shapes                               farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype, class V, class W, class U, class T>
void nurbs_dualshape_function(
  V&                                       funct,   ///< to be filled with shape function values
  W&                                       deriv,   ///< to be filled with shape function derivs
  const U*                                 xi,      ///< xi coordinates
  T&                                       weights, ///< control point weights
  std::vector<Epetra_SerialDenseVector>&   knots,   ///< knot vectors
  DRT::Element&                            ele      ///< element which is considered
  )
{
  switch (DRT::UTILS::DisTypeToDim<distype>::dim)
  {
  case 1:
  {
    volmortar_nurbs_dualshape_function_1D<distype>(funct,deriv,xi,weights,knots,ele);
    break;
  }
  case 2:
  {
    volmortar_nurbs_dualshape_function_2D<distype>(funct,deriv,xi,weights,knots,ele);
    break;
  }
  case 3:
  {
    volmortar_nurbs_dualshape_function_3D<distype>(funct,deriv,xi,weights,knots,ele);
    break;
  }
  default:
    dserror("dimension of the element is not correct");
    break;
  }
  return;
}

}//namespace UTILS
}//namespace VOLMORTAR

#endif /* VOLMORTAR_SHAPE_H_ */
