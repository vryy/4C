/*!----------------------------------------------------------------------
\file mortar_projector.cpp

\brief A class to perform projections of nodes onto opposing MortarElements

\level 1

\maintainer Philipp Farah, Alexander Seitz

*-----------------------------------------------------------------------*/

#include "mortar_projector.H"
#include "mortar_interface.H"
#include "mortar_element.H"
#include "mortar_node.H"
#include "mortar_defines.H"
#include "mortar_calc_utils.H"

#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

/*----------------------------------------------------------------------*
 |  impl. for aux.-plane based projection                    farah 01/14|
 *----------------------------------------------------------------------*/
MORTAR::MortarProjector* MORTAR::MortarProjector::Impl(MortarElement& ele)
{
  switch (ele.Shape())
  {
  case DRT::Element::quad4:
  {
    return MortarProjectorCalc<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return MortarProjectorCalc<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return MortarProjectorCalc<DRT::Element::quad9>::Instance();
  }
  case DRT::Element::tri3:
  {
    return MortarProjectorCalc<DRT::Element::tri3>::Instance();
  }
  case DRT::Element::tri6:
  {
    return MortarProjectorCalc<DRT::Element::tri6>::Instance();
  }
  case DRT::Element::line2:
  {
    return MortarProjectorCalc<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return MortarProjectorCalc<DRT::Element::line3>::Instance();
  }
    //==================================================
    //                     NURBS
    //==================================================
  case DRT::Element::nurbs2:
  {
    return MortarProjectorCalc<DRT::Element::nurbs2>::Instance();
  }
  case DRT::Element::nurbs3:
  {
    return MortarProjectorCalc<DRT::Element::nurbs3>::Instance();
  }
  case DRT::Element::nurbs4:
  {
    return MortarProjectorCalc<DRT::Element::nurbs4>::Instance();
  }
  case DRT::Element::nurbs8:
  {
    return MortarProjectorCalc<DRT::Element::nurbs8>::Instance();
  }
  case DRT::Element::nurbs9:
  {
    return MortarProjectorCalc<DRT::Element::nurbs9>::Instance();
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.",
        ele.Shape(), ele.NumNode());
    break;
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 |  impl. for element based projection                       farah 04/14|
 *----------------------------------------------------------------------*/
MORTAR::MortarProjector* MORTAR::MortarProjector::Impl(MortarElement& sele,
    MortarElement& mele)
{
  switch (sele.Shape())
  {
  case DRT::Element::quad4:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad4,
          DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad4,
          DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad4,
          DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad4,
          DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad4,
          DRT::Element::tri6>::Instance();
    }
    case DRT::Element::nurbs9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad4,
          DRT::Element::nurbs9>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::quad8:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad8,
          DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad8,
          DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad8,
          DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad8,
          DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad8,
          DRT::Element::tri6>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::quad9:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad9,
          DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad9,
          DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad9,
          DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad9,
          DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::quad9,
          DRT::Element::tri6>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::tri3:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri3,
          DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri3,
          DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri3,
          DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri3, DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri3, DRT::Element::tri6>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::tri6:
  {
    switch (mele.Shape())
    {
    case DRT::Element::quad4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri6,
          DRT::Element::quad4>::Instance();
    }
    case DRT::Element::quad8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri6,
          DRT::Element::quad8>::Instance();
    }
    case DRT::Element::quad9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri6,
          DRT::Element::quad9>::Instance();
    }
    case DRT::Element::tri3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri6, DRT::Element::tri3>::Instance();
    }
    case DRT::Element::tri6:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::tri6, DRT::Element::tri6>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::line2:
  {
    switch (mele.Shape())
    {
    case DRT::Element::line2:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::line2,
          DRT::Element::line2>::Instance();
    }
    case DRT::Element::line3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::line2,
          DRT::Element::line3>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::line3:
  {
    switch (mele.Shape())
    {
    case DRT::Element::line2:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::line3,
          DRT::Element::line2>::Instance();
    }
    case DRT::Element::line3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::line3,
          DRT::Element::line3>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
    //==================================================
    //                     NURBS
    //==================================================
  case DRT::Element::nurbs2:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs2:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs2,
          DRT::Element::nurbs2>::Instance();
    }
    case DRT::Element::nurbs3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs2,
          DRT::Element::nurbs3>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::nurbs3:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs2:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs3,
          DRT::Element::nurbs2>::Instance();
    }
    case DRT::Element::nurbs3:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs3,
          DRT::Element::nurbs3>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::nurbs4:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs4,
          DRT::Element::nurbs4>::Instance();
    }
    case DRT::Element::nurbs8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs4,
          DRT::Element::nurbs8>::Instance();
    }
    case DRT::Element::nurbs9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs4,
          DRT::Element::nurbs9>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::nurbs8:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs8,
          DRT::Element::nurbs4>::Instance();
    }
    case DRT::Element::nurbs8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs8,
          DRT::Element::nurbs8>::Instance();
    }
    case DRT::Element::nurbs9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs8,
          DRT::Element::nurbs9>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  case DRT::Element::nurbs9:
  {
    switch (mele.Shape())
    {
    case DRT::Element::nurbs4:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs9,
          DRT::Element::nurbs4>::Instance();
    }
    case DRT::Element::nurbs8:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs9,
          DRT::Element::nurbs8>::Instance();
    }
    case DRT::Element::nurbs9:
    {
      return MortarProjectorCalc_EleBased<DRT::Element::nurbs9,
          DRT::Element::nurbs9>::Instance();
    }
    default:
      dserror("Element shape not supported!");
      break;
    }
    break;
  }
  default:
    dserror("Element shape %d (%d nodes) not activated. Just do it.",
        sele.Shape(), sele.NumNode());
    break;
  }
  return NULL;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
MORTAR::MortarProjectorCalc<distype>::MortarProjectorCalc()
{
  //nothing
}

/*----------------------------------------------------------------------*
 |  ctor ele-based (public)                                  farah 04/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::MortarProjectorCalc_EleBased()
{
  //nothing
}

/*----------------------------------------------------------------------*
 |  Instance (public)                                        farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
MORTAR::MortarProjectorCalc<distype> * MORTAR::MortarProjectorCalc<distype>::Instance(
    bool create)
{
  static MortarProjectorCalc<distype> * instance;
  if (create)
  {
    if (instance == NULL)
      instance = new MortarProjectorCalc<distype>();
  }
  else
  {
    if (instance != NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 |  Instance ele-based (public)                              farah 04/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM> * MORTAR::MortarProjectorCalc_EleBased<
    distypeS, distypeM>::Instance(bool create)
{
  static MortarProjectorCalc_EleBased<distypeS, distypeM> * instance;
  if (create)
  {
    if (instance == NULL)
      instance = new MortarProjectorCalc_EleBased<distypeS, distypeM>();
  }
  else
  {
    if (instance != NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 |  Done (public)                                             farah 01/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void MORTAR::MortarProjectorCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 |  Done ele-based (public)                                  farah 04/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
void MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}


/*----------------------------------------------------------------------*
 |  Project a node along its nodal normal (public)            popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectNodalNormal(
    MORTAR::MortarNode& node, MORTAR::MortarElement& ele, double* xi)
{
  bool ok = true;
  if (ndim_ == 2)
  {
    // local Newton iteration for xi, start in the element middle
    double eta[2] =
    { 0.0, 0.0 };
    double f = 0.0;
    double df = 0.0;
    int k = 0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      f = EvaluateFNodalNormal(node, ele, eta);
      if (abs(f) < MORTARCONVTOL)
        break;
      df = EvaluateGradFNodalNormal(node, ele, eta);
      if (abs(df) < 1.0e-12)
      {
        xi[0] = 1.0e12;
        return false;
        dserror("ERROR: Singular Jacobian for projection");
      }
      eta[0] += (-f) / df;
    }

    // get the result
    xi[0] = eta[0];

    // Newton iteration unconverged
    if (abs(f) > MORTARCONVTOL)
    {
      ok = false;
      xi[0] = 1.0e12;

      // Here (S->M projection) we only give a warning, no error!!!
      // This iteration sometimes diverges, when the projection is far off.
      // These cases are harmless, as these nodes then do not participate in
      // the overlap detection anyway!
      //std::cout << "***WARNING*** ProjectNodalNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and MortarElementID " << ele.Id() << std::endl;
    }

    /*
     // Newton iteration converged
     else
     {
     std::cout << "Newton iteration converged in " << k << " step(s)!" << std::endl
     << "The result is: " << xi[0] << std::endl;
     }*/
  }

  else
    dserror("ERROR: ProjectNodalNormal: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a node along element's normal field (public)      popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectElementNormal(
    MORTAR::MortarNode& node,
    MORTAR::MortarElement& ele,
    double* xi)
{
  bool ok = true;
  if (ndim_ == 2)
  {
    // local Newton iteration for xi, start in the element middle
    double eta[2] =
    { 0.0, 0.0 };
    double f = 0.0;
    double df = 0.0;
    int k = 0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      f = EvaluateFElementNormal(node, ele, eta);
      if (abs(f) < MORTARCONVTOL)
        break;
      df = EvaluateGradFElementNormal(node, ele, eta);
      if (abs(df) < 1.0e-12)
      {
        ok = false;
        xi[0] = 1.0e12;
        return ok;
        dserror("ERROR: Singular Jacobian for projection");
      }
      eta[0] += (-f) / df;
    }

    // get the result
    xi[0] = eta[0];

    // Newton iteration unconverged
    if (abs(f) > MORTARCONVTOL)
    {
      ok = false;
      xi[0] = 1.0e12;

      // Here (M->S projection) we only give a warning, no error!!!
      // This iteration sometimes diverges, when the projection is far off.
      // These cases are harmless, as these nodes then do not participate in
      // the overlap detection anyway!
      //std::cout << "***WARNING*** ProjectElementNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and MortarElementID " << ele.Id() << std::endl;
    }

    /*
     // Newton iteration converged
     else
     {
     std::cout << "Newton iteration converged in " << k << " step(s)!" << std::endl
     << "The result is: " << xi[0] << std::endl;
     }*/
  }

  else
    dserror("ERROR: ProjectElementNormal: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (public)           popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
bool MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::ProjectGaussPoint2D(
    MORTAR::MortarElement& gpele, const double* gpeta,
    MORTAR::MortarElement& ele, double* xi)
{
  bool ok = true;
  if (ndim_ == 2)
  {
    LINALG::Matrix<ns_, 1> val;
    LINALG::Matrix<ndim_, ns_> coord;

    DRT::Node** mynodes = gpele.Nodes();
    if (!mynodes)
      dserror("ERROR: ProjectGaussPoint: Null pointer!");

    // get shape function values and derivatives at gpeta
    if (distypeS == DRT::Element::nurbs2 || distypeS == DRT::Element::nurbs3)
    {
      LINALG::SerialDenseVector auxval(ns_);
      LINALG::SerialDenseMatrix deriv(ns_, 1);
      gpele.EvaluateShape(gpeta, auxval, deriv, gpele.NumNode());

      for (int i = 0; i < ns_; ++i)
        val(i) = auxval(i);
    }
    else
      DRT::UTILS::shape_function_1D(val, gpeta[0], distypeS);

    // get interpolated GP normal and GP coordinates
    double gpn[ndim_];
    double gpx[ndim_];
    for (int i = 0; i < ndim_; ++i)
    {
      gpn[i] = 0.0;
      gpx[i] = 0.0;
    }

    for (int i = 0; i < ns_; ++i)
    {
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(mynodes[i]);

      for (int j = 0; j < ndim_; ++j)
      {
        gpn[j] += val(i) * mymrtrnode->MoData().n()[j];

        coord(j, i) = mymrtrnode->xspatial()[j];

        gpx[j] += val(i) * coord(j, i);
      }
    }

    // local Newton iteration for xi, start in the element middle
    double eta[2] =
    { 0.0, 0.0 };
    double f = 0.0;
    double df = 0.0;
    int k = 0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      f = EvaluateFGaussPoint2D(gpx, gpn, ele, eta);
      if (abs(f) < MORTARCONVTOL)
        break;
      df = EvaluateGradFGaussPoint2D(gpn, ele, eta);
      if (abs(df) < 1.0e-12)
        dserror("ERROR: Singular Jacobian for projection");
      eta[0] += (-f) / df;
    }

    // get the result
    xi[0] = eta[0];

    // Newton iteration unconverged
    if (abs(f) > MORTARCONVTOL)
    {
      ok = false;
      xi[0] = 1.0e12;

      return ok;

//      dserror("ERROR: ProjectGaussPoint: Newton unconverged for GP at xi=%d"
//          " from MortarElementID %i", gpeta[0], gpele.Id());
    }
  }

  else
    dserror("ERROR: ProjectGaussPoint: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal for hermit(public) farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
bool MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::ProjectGaussPointHermit(MORTAR::MortarElement& gpele,
                                                const double* gpeta, MORTAR::MortarElement& ele,
                                                double* xi)
{
  bool ok = true;

  DRT::Node** mynodes = gpele.Nodes();
  double gpn[3] = {0.0, 0.0, 0.0};
  double gpx[3] = {0.0, 0.0, 0.0};
  int nnodes = 0;
  int status = 0;
  int features[2] = {0,0};

  gpele.AdjEleStatus(features);
  status = features[0];
  nnodes = features[1];

  LINALG::SerialDenseVector val(nnodes);
  LINALG::SerialDenseMatrix deriv(nnodes ,2, true);
  LINALG::SerialDenseMatrix coord(3,nnodes);
  if(!mynodes) dserror("ERROR: ProjectGaussPoint: Null pointer!");

  DRT::Node* snodes[4] = {0,0,0,0};
  gpele.HermitEleNodes(snodes, status);

  // get shape function values and derivatives at gpeta

  gpele.EvaluateShape(gpeta,val,deriv,nnodes,false);
  gpele.AdjNodeCoords(coord, status);

  for (int j=0; j<nnodes; ++j)
  {
    gpx[0]+=val[j]*coord(0,j);
    gpx[1]+=val[j]*coord(1,j);
    gpx[2]+=val[j]*coord(2,j);

    MortarNode* cnode = dynamic_cast<MortarNode*>(snodes[j]);
    gpn[0]+=val[j]*cnode->MoData().n()[0];
    gpn[1]+=val[j]*cnode->MoData().n()[1];
    gpn[2]+=val[j]*cnode->MoData().n()[2];
  }

  double gpn_l2 = sqrt(pow(gpn[0],2.0) + pow(gpn[1],2.0) + pow(gpn[2],2.0));

  gpn[0] /= gpn_l2;
  gpn[1] /= gpn_l2;
  gpn[2] /= gpn_l2;

  // local Newton iteration for xi, start in the element middle
  double eta[2] = {0.0, 0.0};
  double f = 0.0;
  double df = 0.0;
  int k=0;

  for (k=0;k<MORTARMAXITER;++k)
  {
    f=EvaluateFGaussPointHermit(gpx,gpn,ele,eta);
    if (abs(f) < MORTARCONVTOL) break;
    df=EvaluateGradFGaussPointHermit(gpn,ele,eta);
    if (abs(df)<1.0e-12) dserror("ERROR: Singular Jacobian for projection");
    eta[0]+=(-f)/df;
  }
  // get the result
  xi[0]=eta[0];

  // Newton iteration unconverged
  if (abs(f) > MORTARCONVTOL)
  {
    ok = false;
    xi[0] = 1.0e12;
  }

  return ok;
}


/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case -- Hermit (public)       farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
double MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::EvaluateFGaussPointHermit(const double* gpx,
                                                    const double* gpn,
                                                    MORTAR::MortarElement& ele,
                                                    const double* eta)
{

  double fval = 0.0;

   // build interpolation of master node coordinates for current eta
   double nx[3] = {0.0, 0.0, 0.0};
   ele.LocalToGlobal(eta,nx,0);

   // subtract GP coordinates
   nx[0]-=gpx[0];
   nx[1]-=gpx[1];
   nx[2]-=gpx[2];

   // calculate F
   fval = nx[0]*gpn[1]-nx[1]*gpn[0];
   return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate grad F for Gauss point case -- Hermit (public)  farah 09/14|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
double MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::EvaluateGradFGaussPointHermit(const double* gpn,
                                                        MORTAR::MortarElement& ele,
                                                        const double* eta)
{
  double fgrad = 0.0;

  // build interpolation of master node coordinates for current eta
  // use shape function derivatives for interpolation (hence "1")
  double nxeta[3] = {0.0, 0.0, 0.0};
  ele.LocalToGlobal(eta,nxeta,1);

  // calculate GradF
  fgrad = nxeta[0]*gpn[1]-nxeta[1]*gpn[0];

  return fgrad;
}

/*----------------------------------------------------------------------*
 | Check projection for warped elements quad4 elements       farah 01/13|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
bool MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::CheckProjection4AUXPLANE(
    MORTAR::MortarElement& ele, double* ngp, double* globgp)
{
  if (ele.Shape() == DRT::Element::tri3)
    dserror("ELEMENT SHAPE TRI3 -- NO WARPING");

    if (ele.Shape()!=DRT::Element::quad4)
    {
      return true;
    }

    int nnode = ele.NumNode();
    DRT::Node** mynodes = ele.Nodes();
    if (!mynodes) dserror("ERROR: Project: Null pointer!");

    //compute base-vectors
    std::vector<double> t0(3);
    std::vector<double> t1(3);
    std::vector<double> auxn(3);
    std::vector<double> auxc(3);
    std::vector<double> proj_gp(3);
    LINALG::Matrix<3, 3> P;
    LINALG::Matrix<3,3> T;
    double n[3]=
    { 0.0 , 0.0 , 0.0};
    double length_t=0.0;
    double length_n=0.0;
    double a1=0.0;
    bool all_negative=true;
    MortarNode* mycnode_0=0;
    MortarNode* mycnode_1=0;
    MortarNode* mycnode_2=0;

    //project gp onto auxn
    for (int i=0;i<nnode;++i)//loop over edges
    {
      if(i==0)
      {
        mycnode_1 = dynamic_cast<MortarNode*> (mynodes[0]);
        if (!mycnode_1) dserror("ERROR: Project: Null pointer!");

        mycnode_0 = dynamic_cast<MortarNode*> (mynodes[3]);
        if (!mycnode_0) dserror("ERROR: Project: Null pointer!");

        mycnode_2 = dynamic_cast<MortarNode*> (mynodes[1]);
        if (!mycnode_2) dserror("ERROR: Project: Null pointer!");
      }
      if(i==3)
      {
        mycnode_1 = dynamic_cast<MortarNode*> (mynodes[3]);
        if (!mycnode_1) dserror("ERROR: Project: Null pointer!");

        mycnode_0 = dynamic_cast<MortarNode*> (mynodes[2]);
        if (!mycnode_0) dserror("ERROR: Project: Null pointer!");

        mycnode_2 = dynamic_cast<MortarNode*> (mynodes[0]);
        if (!mycnode_2) dserror("ERROR: Project: Null pointer!");
      }
      if (i==1 || i==2)
      {
        mycnode_1 = dynamic_cast<MortarNode*> (mynodes[i]);
        if (!mycnode_1) dserror("ERROR: Project: Null pointer!");

        mycnode_0 = dynamic_cast<MortarNode*> (mynodes[i-1]);
        if (!mycnode_0) dserror("ERROR: Project: Null pointer!");

        mycnode_2 = dynamic_cast<MortarNode*> (mynodes[i+1]);
        if (!mycnode_2) dserror("ERROR: Project: Null pointer!");
      }

      //span vectors -- edges
      t0[0]=mycnode_0->xspatial()[0] - mycnode_1->xspatial()[0];
      t0[1]=mycnode_0->xspatial()[1] - mycnode_1->xspatial()[1];
      t0[2]=mycnode_0->xspatial()[2] - mycnode_1->xspatial()[2];

      t1[0]=mycnode_2->xspatial()[0] - mycnode_1->xspatial()[0];
      t1[1]=mycnode_2->xspatial()[1] - mycnode_1->xspatial()[1];
      t1[2]=mycnode_2->xspatial()[2] - mycnode_1->xspatial()[2];

      auxc[0]=mycnode_1->xspatial()[0];
      auxc[1]=mycnode_1->xspatial()[1];
      auxc[2]=mycnode_1->xspatial()[2];

      auxn[0]=t1[1]*t0[2] - t1[2]*t0[1];
      auxn[1]=t1[2]*t0[0] - t1[0]*t0[2];
      auxn[2]=t1[0]*t0[1] - t1[1]*t0[0];

      //fill Projection matrix P
      for (int j=0;j<3;++j)
      P(j,2)=-ngp[j];

      for (int j=0;j<3;++j)
      P(j,0)=t0[j];

      for (int j=0;j<3;++j)
      P(j,1)=t1[j];

      P.Invert();
      double lambda1= P(0,0)*(globgp[0]-auxc[0]) + P(0,1)*(globgp[1]-auxc[1]) + P(0,2)*(globgp[2]-auxc[2]);
      double lambda2= P(1,0)*(globgp[0]-auxc[0]) + P(1,1)*(globgp[1]-auxc[1]) + P(1,2)*(globgp[2]-auxc[2]);
      //double alph= P(2,0)*(globgp[0]-auxc[0]) + P(2,1)*(globgp[1]-auxc[1]) + P(2,2)*(globgp[2]-auxc[2]);

      proj_gp[0]=lambda1*t0[0] + lambda2*t1[0]+auxc[0];
      proj_gp[1]=lambda1*t0[1] + lambda2*t1[1]+auxc[1];
      proj_gp[2]=lambda1*t0[2] + lambda2*t1[2]+auxc[2];
      //std::cout << "proj_gp_AUX4PLANE= " << proj_gp[0] << " "  << proj_gp[1] << " " << proj_gp[2] << std::endl;

      //check
      //edge 1
      n[0]=-(t0[1]*auxn[2] - t0[2]*auxn[1]);
      n[1]=-(t0[2]*auxn[0] - t0[0]*auxn[2]);
      n[2]=-(t0[0]*auxn[1] - t0[1]*auxn[0]);

      length_t=sqrt(t0[0]*t0[0] + t0[1]*t0[1] +t0[2]*t0[2]);
      length_n=sqrt(n[0]*n[0] + n[1]*n[1] +n[2]*n[2]);
      //fill Projection matrix T
      for (int j=0;j<3;++j)
      T(j,0)=n[j]/length_n;

      for (int j=0;j<3;++j)
      T(j,1)=t0[j]/length_t;

      for (int j=0;j<3;++j)
      T(j,2)=auxc[j];

      T.Invert();
      a1=T(0,0)*proj_gp[0] + T(0,1)*proj_gp[1] + T(0,2)*proj_gp[2];

      if (a1>0.0)
      all_negative=false;

      //edge 2
      n[0]=(t1[1]*auxn[2] - t1[2]*auxn[1]);
      n[1]=(t1[2]*auxn[0] - t1[0]*auxn[2]);
      n[2]=(t1[0]*auxn[1] - t1[1]*auxn[0]);

      length_t=sqrt(t1[0]*t1[0] + t1[1]*t1[1] +t1[2]*t1[2]);
      length_n=sqrt(n[0]*n[0] + n[1]*n[1] +n[2]*n[2]);
      //fill Projection matrix T
      for (int j=0;j<3;++j)
      T(j,0)=n[j]/length_n;

      for (int j=0;j<3;++j)
      T(j,1)=t1[j]/length_t;

      for (int j=0;j<3;++j)
      T(j,2)=auxc[j];

      T.Invert();

      a1=T(0,0)*proj_gp[0] + T(0,1)*proj_gp[1] + T(0,2)*proj_gp[2];
      //a2=T(1,0)*proj_gp[0] + T(1,1)*proj_gp[1] + T(1,2)*proj_gp[2];

      if (a1>0.0)
      all_negative=false;
    }

    if(all_negative==false)
    return true;

    return false; // creates error in ProjectGaussPoint3D()
  }

    /*----------------------------------------------------------------------*
     |  Project a Gauss point along its normal (3D)               popp 11/08|
     *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
bool MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::ProjectGaussPoint3D(
    MORTAR::MortarElement& gpele, const double* gpeta,
    MORTAR::MortarElement& ele, double* xi, double& par)
{
  if (ndim_ == 3)
  {
    LINALG::Matrix<ns_, 1> val;
    LINALG::Matrix<ndim_, ns_> coord;
    coord.Clear();

    DRT::Node** mypoints = gpele.Points();
    DRT::Node** mynodes = gpele.Nodes();
    if (!mypoints)
      dserror("ERROR: ProjectGaussPoint: Null pointer!");

    // get shape function values and derivatives at gpeta
    if(distypeS==DRT::Element::nurbs4 || distypeS==DRT::Element::nurbs8 ||
       distypeS==DRT::Element::nurbs9)
    {
      LINALG::SerialDenseVector auxval(ns_);
      LINALG::SerialDenseMatrix deriv(ns_,2);
      gpele.EvaluateShape(gpeta, auxval, deriv, gpele.NumNode());

      for(int i=0;i<ns_;++i)
      val(i)=auxval(i);
    }
    else
    {
      DRT::UTILS::shape_function_2D (val,gpeta[0],gpeta[1],distypeS);
    }

    // get interpolated GP normal and GP coordinates
    double gpn[3];
    double gpx[3];
    for (int i=0;i<ndim_;++i)
    {
      gpn[i] = 0.0;
      gpx[i] = 0.0;
    }

    for (int i=0;i<ns_;++i)
    {
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (mypoints[i]);

      for (int j=0;j<ndim_;++j)
      {
        coord(j,i) = mymrtrnode->xspatial()[j];

        gpx[j] += val(i)*coord(j,i);
      }
    }

    for (int i=0;i<gpele.NumNode();++i)
    {
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (mynodes[i]);
      for (int j=0;j<ndim_;++j)
      {
        gpn[j] += val(i)*mymrtrnode->MoData().n()[j];
      }
    }

    // start in the element center
    DRT::Element::DiscretizationType dt = ele.Shape();
    double eta[2] =
    { 0.0, 0.0};
    if (dt==DRT::Element::tri3 || dt==DRT::Element::tri6)
    {
      eta[0] = 1.0/3.0;
      eta[1] = 1.0/3.0;
    }

    // auxiliary variable
    double alpha = 0.0;

    // function f (vector-valued)
    double f[3] =
    { 0.0, 0.0, 0.0};

    // gradient of f (df/deta[0], df/deta[1], df/dalpha)
    LINALG::Matrix<3, 3> df;
    // start iteration
    int k= 0;
    double conv=0.0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      EvaluateFGaussPoint3D(f, gpx, gpn, ele, eta, alpha);
      conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
      if (conv <= MORTARCONVTOL)
        break;
      EvaluateGradFGaussPoint3D(df, gpx, gpn, ele, eta, alpha);

      // safety check: if projection normal is parallel to the master element --> det can be zero
      double det = df.Determinant();
      if (det > -1e-12 and det < 1e-12)
      {
        std::cout
            << "WARNING: GPProjection parallel to master element --> GP skipped for this master element!"
            << std::endl;
        // leave here
        return false;
      }

      // solve deta = - inv(df) * f
      double jacdet = df.Invert();
      if (abs(jacdet) < 1.0e-12)
        dserror("ERROR: Singular Jacobian for projection");

      // update eta and alpha
      eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
      eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
      alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];

      //Projection Check
//      if (k==MORTARMAXITER-1)
//      {
//        bool check = CheckProjection4AUXPLANE(ele, gpn,gpx);
//        if (check==false)
//          dserror("!!! STOP !!!   -->   Projection Error: Newton unconverged but GP on mele !!!");
//      }
    }

    // Newton iteration unconverged
    if (conv > MORTARCONVTOL)
    {
      xi[0] = 1e12;
      xi[1] = 1e12;
    }
    else
    {
      xi[0] = eta[0];
      xi[1] = eta[1];
    }

    par = alpha;
  }

  else
    dserror("ERROR: ProjectGaussPoint: Called 3D version for 2D problem!");

  return true;
}
/*----------------------------------------------------------------------*
 |  Project a Gauss point along AuxPlane normal (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectGaussPointAuxn3D(
    const double* globgp, const double* auxn, MORTAR::MortarElement& ele,
    double* xi, double& par)
{
  if (ndim_ == 3)
  {
    // start in the element center
    DRT::Element::DiscretizationType dt = ele.Shape();
    double eta[2] =
    { 0.0, 0.0 };
    if (dt == DRT::Element::tri3 || dt == DRT::Element::tri6)
    {
      eta[0] = 1.0 / 3.0;
      eta[1] = 1.0 / 3.0;
    }

    // auxiliary variable
    double alpha = 0.0;

    // function f (vector-valued)
    double f[3] =
    { 0.0, 0.0, 0.0 };

    // gradient of f (df/deta[0], df/deta[1], df/dalpha)
    LINALG::Matrix<3, 3> df;

    // start iteration
    int k = 0;
    double conv = 0.0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      EvaluateFGaussPointAuxn3D(f, globgp, auxn, ele, eta, alpha);
      conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
      //std::cout << "Iteration " << k << ": -> |f|=" << conv << std::endl;
      if (conv <= MORTARCONVTOL)
        break;
      EvaluateGradFGaussPointAuxn3D(df, globgp, auxn, ele, eta, alpha);

      // solve deta = - inv(df) * f
      double jacdet = df.Invert();
      if (abs(jacdet) < 1.0e-12)
      {
        dserror("ERROR: Singular Jacobian for projection");
      }

      // update eta and alpha
      eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
      eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
      alpha  += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
    }

    // Newton iteration unconverged
    if (conv > MORTARCONVTOL)
    {
      std::cout << "res= " << conv << std::endl;
      dserror("ERROR: ProjectGaussPointAuxn3D: Newton unconverged for GP"
          "at xi = (%f,%f,%f) onto MortarElementID %i", globgp[0], globgp[1],
          globgp[2], ele.Id());
    }


    // Newton iteration converged
    xi[0] = eta[0];
    xi[1] = eta[1];
    par   = alpha;
    //std::cout << "Newton iteration converged in " << k << " steps!" << std::endl;
  }

  else
    dserror("ERROR: ProjectGaussPoint: Called 3D version for 2D problem!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNormal3D(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist)
{
  if(ndim_!=3)
    dserror("ERROR: ProjectSNodeByMNormal3D is only for 3D problems!");

  // start in the element center
  double eta[2] =  { 0.0, 0.0 };
  if (distype == DRT::Element::tri3 || distype == DRT::Element::tri6)
  {
    eta[0] = 1.0 / 3.0;
    eta[1] = 1.0 / 3.0;
  }

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  double f[3] =  { 0.0, 0.0, 0.0 };

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3]         = { 0.0, 0.0, 0.0 };
    double xs[3]         = { 0.0, 0.0, 0.0 };
    double unormal[3]    = { 0.0, 0.0, 0.0 };
    double normal[3]     = { 0.0, 0.0, 0.0 };
    double normalpart[3] = { 0.0, 0.0, 0.0 };

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta,unormal);
    for (int i = 0; i < 3; ++i)
      normalpart[i] = unormal[i]*alpha;

    for (int i = 0; i < 3; ++i)
      normal[i] = unormal[i]*length;

    // calc xslave
    for (int i = 0; i < 3; ++i)
      xs[i] = snode.xspatial()[i];

    //calculate F
    for (int i = 0; i < 3; ++i)
      f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL)
      break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // master coordinate grad
    double meta0[3] = { 0.0, 0.0, 0.0 }; // x,xi_0
    double meta1[3] = { 0.0, 0.0, 0.0 }; // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta1, 2);

    // normal grad
    LINALG::Matrix<3,n_> secderiv;
    double meta00[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_0
    double meta11[3]=  { 0.0, 0.0, 0.0 }; // x,xi_1 xi_1
    double meta01[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_1

    DRT::UTILS::shape_function_2D_deriv2(secderiv,eta[0],eta[1],distype);

    for (int i = 0;i<n_;++i)
    {
      MortarNode* mymnode = dynamic_cast<MortarNode*> (mele.Nodes()[i]);
      if (!mymnode)
        dserror("ERROR: Null pointer!");
      for (int d = 0;d<3;++d)
      {
        meta00[d] += secderiv(0,i) * mymnode->xspatial()[d];
        meta11[d] += secderiv(1,i) * mymnode->xspatial()[d];
        meta01[d] += secderiv(2,i) * mymnode->xspatial()[d];
      }
    }

    double naux_0[3] = { 0.0, 0.0, 0.0 };
    double naux_1[3] = { 0.0, 0.0, 0.0 };

    // normal grad xi_0
    naux_0[0] = (meta00[1]*meta1[2]-meta00[2]*meta1[1]);
    naux_0[1] = (meta00[2]*meta1[0]-meta00[0]*meta1[2]);
    naux_0[2] = (meta00[0]*meta1[1]-meta00[1]*meta1[0]);

    naux_0[0] += (meta0[1]*meta01[2]-meta0[2]*meta01[1]);
    naux_0[1] += (meta0[2]*meta01[0]-meta0[0]*meta01[2]);
    naux_0[2] += (meta0[0]*meta01[1]-meta0[1]*meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1]*meta1[2]-meta01[2]*meta1[1]);
    naux_1[1] = (meta01[2]*meta1[0]-meta01[0]*meta1[2]);
    naux_1[2] = (meta01[0]*meta1[1]-meta01[1]*meta1[0]);

    naux_1[0] += (meta0[1]*meta11[2]-meta0[2]*meta11[1]);
    naux_1[1] += (meta0[2]*meta11[0]-meta0[0]*meta11[2]);
    naux_1[2] += (meta0[0]*meta11[1]-meta0[1]*meta11[0]);

    double n_0[3] = { 0.0, 0.0, 0.0 };
    double n_1[3] = { 0.0, 0.0, 0.0 };

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i=0;i<3; ++i)
      fac0 += naux_0[i]* normal[i];

    for (int i = 0;i<3;++i)
      n_0[i] = naux_0[i]/length - fac0*normal[i]/(length*length*length);

    for (int i=0;i<3; ++i)
      fac1 += naux_1[i]* normal[i];

    for (int i = 0;i<3;++i)
      n_1[i] = naux_1[i]/length - fac1*normal[i]/(length*length*length);

    // evaluate function f gradient
    for (int i = 0; i < 3; ++i)
    {
      df(i, 0) = meta0[i] + n_0[i];
      df(i, 1) = meta1[i] + n_1[i];
      df(i, 2) = unormal[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
      dserror("ERROR: Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha  += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  } // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
    dserror("ERROR: Projector not converged!");

  // Newton iteration converged
  xi[0] = eta[0];
  xi[1] = eta[1];
  mele.ComputeUnitNormalAtXi(eta,normal);
  dist  = alpha;

  // bye bye
  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNormal3DLin(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  if(ndim_!=3)
    dserror("ERROR: ProjectSNodeByMNormal3D is only for 3D problems!");

  // start in the element center
  double eta[2] =  { 0.0, 0.0 };
  if (distype == DRT::Element::tri3 || distype == DRT::Element::tri6)
  {
    eta[0] = 1.0 / 3.0;
    eta[1] = 1.0 / 3.0;
  }

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  double f[3] =  { 0.0, 0.0, 0.0 };

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3]         = { 0.0, 0.0, 0.0 };
    double xs[3]         = { 0.0, 0.0, 0.0 };
    double unormal[3]    = { 0.0, 0.0, 0.0 };
    double normal[3]     = { 0.0, 0.0, 0.0 };
    double normalpart[3] = { 0.0, 0.0, 0.0 };

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta,unormal);
    for (int i = 0; i < 3; ++i)
      normalpart[i] = unormal[i]*alpha;

    for (int i = 0; i < 3; ++i)
      normal[i] = unormal[i]*length;

    // calc xslave
    for (int i = 0; i < 3; ++i)
      xs[i] = snode.xspatial()[i];

    //calculate F
    for (int i = 0; i < 3; ++i)
      f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL)
      break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // master coordinate grad
    double meta0[3] = { 0.0, 0.0, 0.0 }; // x,xi_0
    double meta1[3] = { 0.0, 0.0, 0.0 }; // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta1, 2);

    // normal grad
    LINALG::Matrix<3,n_> secderiv;
    double meta00[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_0
    double meta11[3]=  { 0.0, 0.0, 0.0 }; // x,xi_1 xi_1
    double meta01[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_1

    DRT::UTILS::shape_function_2D_deriv2(secderiv,eta[0],eta[1],distype);

    for (int i = 0;i<n_;++i)
    {
      MortarNode* mymnode = dynamic_cast<MortarNode*> (mele.Nodes()[i]);
      if (!mymnode)
        dserror("ERROR: Null pointer!");
      for (int d = 0;d<3;++d)
      {
        meta00[d] += secderiv(0,i) * mymnode->xspatial()[d];
        meta11[d] += secderiv(1,i) * mymnode->xspatial()[d];
        meta01[d] += secderiv(2,i) * mymnode->xspatial()[d];
      }
    }

    double naux_0[3] = { 0.0, 0.0, 0.0 };
    double naux_1[3] = { 0.0, 0.0, 0.0 };

    // normal grad xi_0
    naux_0[0] = (meta00[1]*meta1[2]-meta00[2]*meta1[1]);
    naux_0[1] = (meta00[2]*meta1[0]-meta00[0]*meta1[2]);
    naux_0[2] = (meta00[0]*meta1[1]-meta00[1]*meta1[0]);

    naux_0[0] += (meta0[1]*meta01[2]-meta0[2]*meta01[1]);
    naux_0[1] += (meta0[2]*meta01[0]-meta0[0]*meta01[2]);
    naux_0[2] += (meta0[0]*meta01[1]-meta0[1]*meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1]*meta1[2]-meta01[2]*meta1[1]);
    naux_1[1] = (meta01[2]*meta1[0]-meta01[0]*meta1[2]);
    naux_1[2] = (meta01[0]*meta1[1]-meta01[1]*meta1[0]);

    naux_1[0] += (meta0[1]*meta11[2]-meta0[2]*meta11[1]);
    naux_1[1] += (meta0[2]*meta11[0]-meta0[0]*meta11[2]);
    naux_1[2] += (meta0[0]*meta11[1]-meta0[1]*meta11[0]);

    double n_0[3] = { 0.0, 0.0, 0.0 };
    double n_1[3] = { 0.0, 0.0, 0.0 };

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i=0;i<3; ++i)
      fac0 += naux_0[i]* normal[i];

    for (int i = 0;i<3;++i)
      n_0[i] = naux_0[i]/length - fac0*normal[i]/(length*length*length);

    for (int i=0;i<3; ++i)
      fac1 += naux_1[i]* normal[i];

    for (int i = 0;i<3;++i)
      n_1[i] = naux_1[i]/length - fac1*normal[i]/(length*length*length);

    // evaluate function f gradient
    for (int i = 0; i < 3; ++i)
    {
      df(i, 0) = meta0[i] + n_0[i];
      df(i, 1) = meta1[i] + n_1[i];
      df(i, 2) = unormal[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
      dserror("ERROR: Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha  += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  } // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
    dserror("ERROR: Projector not converged!");

  // Newton iteration converged
  xi[0] = eta[0];
  xi[1] = eta[1];
  mele.ComputeUnitNormalAtXi(eta,normal);
  dist  = alpha;

  // bye bye
  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 05/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNodalNormal2DLin(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  if(ndim_!=2)
    dserror("ERROR: ProjectSNodeByMNormal2DLin is only for 2D problems!");

  // start in the element center
  double eta[2] =  { 0.0, 0.0 };

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  double f[3] =  { 0.0, 0.0, 0.0 };

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3]         = { 0.0, 0.0, 0.0 };
    double xs[3]         = { 0.0, 0.0, 0.0 };
    double unormal[3]    = { 0.0, 0.0, 0.0 };
    double normal_k[3]   = { 0.0, 0.0, 0.0 };
    double normalpart[3] = { 0.0, 0.0, 0.0 };

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeAveragedUnitNormalAtXi(eta,unormal);

    for (int i = 0; i < 3; ++i)
      normalpart[i] = unormal[i]*alpha;

    for (int i = 0; i < 3; ++i)
      normal_k[i] = unormal[i]*length;

    // calc xslave
    for (int i = 0; i < 3; ++i)
      xs[i] = snode.xspatial()[i];

    //calculate F
    for (int i = 0; i < 3; ++i)
      f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL)
      break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    double meta0[3] = { 0.0, 0.0, 0.0 }; // x,xi_0
    double meta1[3] = { 0.0, 0.0, 1.0 }; // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);

    // normal grad
    LINALG::Matrix<1,n_> secderiv;
    double meta00[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_0
    double meta11[3]=  { 0.0, 0.0, 0.0 }; // x,xi_1 xi_1
    double meta01[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_1

    DRT::UTILS::shape_function_1D_deriv2(secderiv,eta[0],distype);

    for (int i = 0;i<n_;++i)
    {
      MortarNode* mymnode = dynamic_cast<MortarNode*> (mele.Nodes()[i]);
      if (!mymnode)
        dserror("ERROR: Null pointer!");
      for (int d = 0;d<3;++d)
        meta00[d] += secderiv(0,i) * mymnode->xspatial()[d];
    }

    double naux_0[3] = { 0.0, 0.0, 0.0 };
    double naux_1[3] = { 0.0, 0.0, 0.0 };

    // normal grad xi_0
    naux_0[0] = (meta00[1]*meta1[2]-meta00[2]*meta1[1]);
    naux_0[1] = (meta00[2]*meta1[0]-meta00[0]*meta1[2]);
    naux_0[2] = (meta00[0]*meta1[1]-meta00[1]*meta1[0]);

    naux_0[0] += (meta0[1]*meta01[2]-meta0[2]*meta01[1]);
    naux_0[1] += (meta0[2]*meta01[0]-meta0[0]*meta01[2]);
    naux_0[2] += (meta0[0]*meta01[1]-meta0[1]*meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1]*meta1[2]-meta01[2]*meta1[1]);
    naux_1[1] = (meta01[2]*meta1[0]-meta01[0]*meta1[2]);
    naux_1[2] = (meta01[0]*meta1[1]-meta01[1]*meta1[0]);

    naux_1[0] += (meta0[1]*meta11[2]-meta0[2]*meta11[1]);
    naux_1[1] += (meta0[2]*meta11[0]-meta0[0]*meta11[2]);
    naux_1[2] += (meta0[0]*meta11[1]-meta0[1]*meta11[0]);

    double n_0[3] = { 0.0, 0.0, 0.0 };
    double n_1[3] = { 0.0, 0.0, 0.0 };

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i=0;i<3; ++i)
      fac0 += naux_0[i]* normal_k[i];

    for (int i = 0;i<3;++i)
      n_0[i] = naux_0[i]/length - fac0*normal_k[i]/(length*length*length);

    for (int i=0;i<3; ++i)
      fac1 += naux_1[i]* normal_k[i];

    for (int i = 0;i<3;++i)
      n_1[i] = naux_1[i]/length - fac1*normal_k[i]/(length*length*length);

    // evaluate function f gradient
    for (int i = 0; i < 3; ++i)
    {
      df(i, 0) = meta0[i] + n_0[i];
      df(i, 1) = meta1[i] + n_1[i];
      df(i, 2) = unormal[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
      dserror("ERROR: Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha  += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  } // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
    dserror("ERROR: Projector not converged!");

  //**********************************************
  //   Get stuff                                //
  //**********************************************
  // Newton iteration converged
  xi[0] = eta[0];
  double normallength = mele.ComputeUnitNormalAtXi(eta,normal);
  dist  = alpha;


  //**********************************************
  //   Lin deta = - inv(dF) * Lin F             //
  //**********************************************
  // prepare linearizations
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  std::vector<GEN::pairedvector<int,double> > etaLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > fLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > xmLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > normalpartLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > xsLin(3,1000);

  //--------------------------
  // master part:
  LINALG::Matrix<n_,1> val;
  DRT::UTILS::shape_function_1D(val,eta[0],distype);

  for (int i = 0;i<n_;++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("ERROR: Cannot find master node");
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    for(int k = 0; k<2; ++k)
      (xmLin[k])[mnode->Dofs()[k]] += val(i) ;
  }

  //--------------------------
  // normal part:
  std::vector<GEN::pairedvector<int,double> > x_0Lin(3,1000);
  std::vector<GEN::pairedvector<int,double> > auxnormalLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > auxnormalunitLin(3,1000);

  LINALG::Matrix<1,n_>   deriv1;
  DRT::UTILS::shape_function_1D_deriv1(deriv1,eta[0],distype);
  for (int i = 0;i<n_;++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("ERROR: Cannot find master node");
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    for(int k = 0; k<2; ++k)
      (x_0Lin[k])[mnode->Dofs()[k]] += deriv1(i) ;
  }

  // cross product linearization
  for (_CI p=x_0Lin[1].begin();p!=x_0Lin[1].end();++p)
    (auxnormalLin[0])[p->first] += (p->second);
  for (_CI p=x_0Lin[0].begin();p!=x_0Lin[0].end();++p)
    (auxnormalLin[1])[p->first] -= (p->second);

  // calc normalpart without alpha
  double linnormalaux[3] = {0.0, 0.0, 0.0};
  linnormalaux[0] = normal[0]*normallength;
  linnormalaux[1] = normal[1]*normallength;
  linnormalaux[2] = normal[2]*normallength;

  // derivative weighting matrix for current element
  LINALG::Matrix<3, 3>W;
  const double lcubeinv = 1.0 / (normallength * normallength * normallength);

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      W(j, k) = -lcubeinv * linnormalaux[j] * linnormalaux[k];
      if (j == k)
        W(j, k) += 1 / normallength;
    }
  }

  // row loop
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=auxnormalLin[0].begin();p!=auxnormalLin[0].end();++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j,0);

    for (_CI p=auxnormalLin[1].begin();p!=auxnormalLin[1].end();++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j,1);

    for (_CI p=auxnormalLin[2].begin();p!=auxnormalLin[2].end();++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j,2);
  }

  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=auxnormalunitLin[j].begin();p!=auxnormalunitLin[j].end();++p)
      (normalpartLin[j])[p->first] += (p->second) * alpha;
  }


  //--------------------------
  // slave part:
  for(int k = 0; k<2; ++k)
    (xsLin[k])[snode.Dofs()[k]] += 1.0;

  // All terms:
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=xmLin[j].begin();p!=xmLin[j].end();++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p=normalpartLin[j].begin();p!=normalpartLin[j].end();++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p=xsLin[j].begin();p!=xsLin[j].end();++p)
      (fLin[j])[p->first] -= (p->second);
  }


  for(int i = 0; i<3 ; ++i)
  {
    for(int k = 0; k<3 ; ++k)
    {
      for (_CI p=fLin[i].begin();p!=fLin[i].end();++p)
        (etaLin[k])[p->first] -= (p->second) * df(k,i);
    }
  }

  //**********************************************
  //   Lin N                                    //
  //**********************************************
  std::vector<GEN::pairedvector<int,double> > x_0Linnew(3,1000);
  std::vector<GEN::pairedvector<int,double> > normaltolineLinaux(3,1000);

  LINALG::Matrix<1,n_>   deriv;
  DRT::UTILS::shape_function_1D_deriv1(deriv,eta[0],distype);

  LINALG::Matrix<1,n_>   deriv2;
  DRT::UTILS::shape_function_1D_deriv2(deriv2,eta[0],distype);
  for (int i = 0;i<n_;++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("ERROR: Cannot find master node");
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    for(int k = 0; k<2; ++k)
      (x_0Linnew[k])[mnode->Dofs()[k]] += deriv(i);

    for(int k = 0; k<2; ++k)
    {
      for (_CI p=etaLin[0].begin();p!=etaLin[0].end();++p)
        (x_0Linnew[k])[p->first] += (p->second) * deriv2(i) * mnode->xspatial()[k];
    }
  }

  // cross product linearization
  for (_CI p=x_0Linnew[1].begin();p!=x_0Linnew[1].end();++p)
    (normaltolineLinaux[0])[p->first] += (p->second);
  for (_CI p=x_0Linnew[0].begin();p!=x_0Linnew[0].end();++p)
    (normaltolineLinaux[1])[p->first] -= (p->second);

  //normalize lin
  LINALG::Matrix<3, 3>Wfinal;
//  const double lcubeinv = 1.0 / (normallength * normallength * normallength);

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      Wfinal(j, k) = -lcubeinv * linnormalaux[j] * linnormalaux[k];
      if (j == k)
        Wfinal(j, k) += 1 / normallength;
    }
  }

  // row loop
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=normaltolineLinaux[0].begin();p!=normaltolineLinaux[0].end();++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j,0);

    for (_CI p=normaltolineLinaux[1].begin();p!=normaltolineLinaux[1].end();++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j,1);

    for (_CI p=normaltolineLinaux[2].begin();p!=normaltolineLinaux[2].end();++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j,2);
  }

  // bye bye
  return true;
}



/*----------------------------------------------------------------------*
 |  Project snode onto melement with master element normal   farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNormal2D(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist)
{
  if(ndim_!=2)
    dserror("ERROR: ProjectSNodeByMNormal2D is only for 2D problems!");

  // start in the element center
  double eta[2] =  { 0.0, 0.0 };

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  double f[3] =  { 0.0, 0.0, 0.0 };

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3]         = { 0.0, 0.0, 0.0 };
    double xs[3]         = { 0.0, 0.0, 0.0 };
    double unormal[3]    = { 0.0, 0.0, 0.0 };
    double normal[3]     = { 0.0, 0.0, 0.0 };
    double normalpart[3] = { 0.0, 0.0, 0.0 };

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta,unormal);
    for (int i = 0; i < 3; ++i)
      normalpart[i] = unormal[i]*alpha;

    for (int i = 0; i < 3; ++i)
      normal[i] = unormal[i]*length;

    // calc xslave
    for (int i = 0; i < 3; ++i)
      xs[i] = snode.xspatial()[i];

    //calculate F
    for (int i = 0; i < 3; ++i)
      f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL)
      break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    double meta0[3] = { 0.0, 0.0, 0.0 }; // x,xi_0
    double meta1[3] = { 0.0, 0.0, 1.0 }; // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);

    // normal grad
    LINALG::Matrix<1,n_> secderiv;
    double meta00[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_0
    double meta11[3]=  { 0.0, 0.0, 0.0 }; // x,xi_1 xi_1
    double meta01[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_1

    DRT::UTILS::shape_function_1D_deriv2(secderiv,eta[0],distype);

    for (int i = 0;i<n_;++i)
    {
      MortarNode* mymnode = dynamic_cast<MortarNode*> (mele.Nodes()[i]);
      if (!mymnode)
        dserror("ERROR: Null pointer!");
      for (int d = 0;d<3;++d)
        meta00[d] += secderiv(0,i) * mymnode->xspatial()[d];
    }

    double naux_0[3] = { 0.0, 0.0, 0.0 };
    double naux_1[3] = { 0.0, 0.0, 0.0 };

    // normal grad xi_0
    naux_0[0] = (meta00[1]*meta1[2]-meta00[2]*meta1[1]);
    naux_0[1] = (meta00[2]*meta1[0]-meta00[0]*meta1[2]);
    naux_0[2] = (meta00[0]*meta1[1]-meta00[1]*meta1[0]);

    naux_0[0] += (meta0[1]*meta01[2]-meta0[2]*meta01[1]);
    naux_0[1] += (meta0[2]*meta01[0]-meta0[0]*meta01[2]);
    naux_0[2] += (meta0[0]*meta01[1]-meta0[1]*meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1]*meta1[2]-meta01[2]*meta1[1]);
    naux_1[1] = (meta01[2]*meta1[0]-meta01[0]*meta1[2]);
    naux_1[2] = (meta01[0]*meta1[1]-meta01[1]*meta1[0]);

    naux_1[0] += (meta0[1]*meta11[2]-meta0[2]*meta11[1]);
    naux_1[1] += (meta0[2]*meta11[0]-meta0[0]*meta11[2]);
    naux_1[2] += (meta0[0]*meta11[1]-meta0[1]*meta11[0]);

    double n_0[3] = { 0.0, 0.0, 0.0 };
    double n_1[3] = { 0.0, 0.0, 0.0 };

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i=0;i<3; ++i)
      fac0 += naux_0[i]* normal[i];

    for (int i = 0;i<3;++i)
      n_0[i] = naux_0[i]/length - fac0*normal[i]/(length*length*length);

    for (int i=0;i<3; ++i)
      fac1 += naux_1[i]* normal[i];

    for (int i = 0;i<3;++i)
      n_1[i] = naux_1[i]/length - fac1*normal[i]/(length*length*length);

    // evaluate function f gradient
    for (int i = 0; i < 3; ++i)
    {
      df(i, 0) = meta0[i] + n_0[i];
      df(i, 1) = meta1[i] + n_1[i];
      df(i, 2) = unormal[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
      dserror("ERROR: Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha  += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  } // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
    dserror("ERROR: Projector not converged!");

  // Newton iteration converged
  xi[0] = eta[0];
  mele.ComputeUnitNormalAtXi(eta,normal);
  dist  = alpha;

  // bye bye
  return true;
}

/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal + Lin     farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNormal2DLin(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  if(ndim_!=2)
    dserror("ERROR: ProjectSNodeByMNormal2D is only for 2D problems!");

  // start in the element center
  double eta[2] =  { 0.0, 0.0 };

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  double f[3] =  { 0.0, 0.0, 0.0 };

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3]         = { 0.0, 0.0, 0.0 };
    double xs[3]         = { 0.0, 0.0, 0.0 };
    double unormal[3]    = { 0.0, 0.0, 0.0 };
    double normal_k[3]   = { 0.0, 0.0, 0.0 };
    double normalpart[3] = { 0.0, 0.0, 0.0 };

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta,unormal);
    for (int i = 0; i < 3; ++i)
      normalpart[i] = unormal[i]*alpha;

    for (int i = 0; i < 3; ++i)
      normal_k[i] = unormal[i]*length;

    // calc xslave
    for (int i = 0; i < 3; ++i)
      xs[i] = snode.xspatial()[i];

    //calculate F
    for (int i = 0; i < 3; ++i)
      f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL)
      break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    double meta0[3] = { 0.0, 0.0, 0.0 }; // x,xi_0
    double meta1[3] = { 0.0, 0.0, 1.0 }; // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);

    // normal grad
    LINALG::Matrix<1,n_> secderiv;
    double meta00[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_0
    double meta11[3]=  { 0.0, 0.0, 0.0 }; // x,xi_1 xi_1
    double meta01[3]=  { 0.0, 0.0, 0.0 }; // x,xi_0 xi_1

    DRT::UTILS::shape_function_1D_deriv2(secderiv,eta[0],distype);

    for (int i = 0;i<n_;++i)
    {
      MortarNode* mymnode = dynamic_cast<MortarNode*> (mele.Nodes()[i]);
      if (!mymnode)
        dserror("ERROR: Null pointer!");
      for (int d = 0;d<3;++d)
        meta00[d] += secderiv(0,i) * mymnode->xspatial()[d];
    }

    double naux_0[3] = { 0.0, 0.0, 0.0 };
    double naux_1[3] = { 0.0, 0.0, 0.0 };

    // normal grad xi_0
    naux_0[0] = (meta00[1]*meta1[2]-meta00[2]*meta1[1]);
    naux_0[1] = (meta00[2]*meta1[0]-meta00[0]*meta1[2]);
    naux_0[2] = (meta00[0]*meta1[1]-meta00[1]*meta1[0]);

    naux_0[0] += (meta0[1]*meta01[2]-meta0[2]*meta01[1]);
    naux_0[1] += (meta0[2]*meta01[0]-meta0[0]*meta01[2]);
    naux_0[2] += (meta0[0]*meta01[1]-meta0[1]*meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1]*meta1[2]-meta01[2]*meta1[1]);
    naux_1[1] = (meta01[2]*meta1[0]-meta01[0]*meta1[2]);
    naux_1[2] = (meta01[0]*meta1[1]-meta01[1]*meta1[0]);

    naux_1[0] += (meta0[1]*meta11[2]-meta0[2]*meta11[1]);
    naux_1[1] += (meta0[2]*meta11[0]-meta0[0]*meta11[2]);
    naux_1[2] += (meta0[0]*meta11[1]-meta0[1]*meta11[0]);

    double n_0[3] = { 0.0, 0.0, 0.0 };
    double n_1[3] = { 0.0, 0.0, 0.0 };

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i=0;i<3; ++i)
      fac0 += naux_0[i]* normal_k[i];

    for (int i = 0;i<3;++i)
      n_0[i] = naux_0[i]/length - fac0*normal_k[i]/(length*length*length);

    for (int i=0;i<3; ++i)
      fac1 += naux_1[i]* normal_k[i];

    for (int i = 0;i<3;++i)
      n_1[i] = naux_1[i]/length - fac1*normal_k[i]/(length*length*length);

    // evaluate function f gradient
    for (int i = 0; i < 3; ++i)
    {
      df(i, 0) = meta0[i] + n_0[i];
      df(i, 1) = meta1[i] + n_1[i];
      df(i, 2) = unormal[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12)
      dserror("ERROR: Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha  += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  } // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
    dserror("ERROR: Projector not converged!");

  //**********************************************
  //   Get stuff                                //
  //**********************************************
  // Newton iteration converged
  xi[0] = eta[0];
  double normallength = mele.ComputeUnitNormalAtXi(eta,normal);
  dist  = alpha;


  //**********************************************
  //   Lin deta = - inv(dF) * Lin F             //
  //**********************************************
  // prepare linearizations
  typedef GEN::pairedvector<int,double>::const_iterator _CI;

  std::vector<GEN::pairedvector<int,double> > etaLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > fLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > xmLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > normalpartLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > xsLin(3,1000);

  //--------------------------
  // master part:
  LINALG::Matrix<n_,1> val;
  DRT::UTILS::shape_function_1D(val,eta[0],distype);

  for (int i = 0;i<n_;++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("ERROR: Cannot find master node");
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    for(int k = 0; k<2; ++k)
      (xmLin[k])[mnode->Dofs()[k]] += val(i) ;
  }

  //--------------------------
  // normal part:
  std::vector<GEN::pairedvector<int,double> > x_0Lin(3,1000);
  std::vector<GEN::pairedvector<int,double> > auxnormalLin(3,1000);
  std::vector<GEN::pairedvector<int,double> > auxnormalunitLin(3,1000);

  LINALG::Matrix<1,n_>   deriv1;
  DRT::UTILS::shape_function_1D_deriv1(deriv1,eta[0],distype);
  for (int i = 0;i<n_;++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("ERROR: Cannot find master node");
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    for(int k = 0; k<2; ++k)
      (x_0Lin[k])[mnode->Dofs()[k]] += deriv1(i) ;
  }

  // cross product linearization
  for (_CI p=x_0Lin[1].begin();p!=x_0Lin[1].end();++p)
    (auxnormalLin[0])[p->first] += (p->second);
  for (_CI p=x_0Lin[0].begin();p!=x_0Lin[0].end();++p)
    (auxnormalLin[1])[p->first] -= (p->second);

  // calc normalpart without alpha
  double linnormalaux[3] = {0.0, 0.0, 0.0};
  linnormalaux[0] = normal[0]*normallength;
  linnormalaux[1] = normal[1]*normallength;
  linnormalaux[2] = normal[2]*normallength;

  // derivative weighting matrix for current element
  LINALG::Matrix<3, 3>W;
  const double lcubeinv = 1.0 / (normallength * normallength * normallength);

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      W(j, k) = -lcubeinv * linnormalaux[j] * linnormalaux[k];
      if (j == k)
        W(j, k) += 1 / normallength;
    }
  }

  // row loop
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=auxnormalLin[0].begin();p!=auxnormalLin[0].end();++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j,0);

    for (_CI p=auxnormalLin[1].begin();p!=auxnormalLin[1].end();++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j,1);

    for (_CI p=auxnormalLin[2].begin();p!=auxnormalLin[2].end();++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j,2);
  }

  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=auxnormalunitLin[j].begin();p!=auxnormalunitLin[j].end();++p)
      (normalpartLin[j])[p->first] += (p->second) * alpha;
  }


  //--------------------------
  // slave part:
  for(int k = 0; k<2; ++k)
    (xsLin[k])[snode.Dofs()[k]] += 1.0;

  // All terms:
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=xmLin[j].begin();p!=xmLin[j].end();++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p=normalpartLin[j].begin();p!=normalpartLin[j].end();++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p=xsLin[j].begin();p!=xsLin[j].end();++p)
      (fLin[j])[p->first] -= (p->second);
  }


  for(int i = 0; i<3 ; ++i)
  {
    for(int k = 0; k<3 ; ++k)
    {
      for (_CI p=fLin[i].begin();p!=fLin[i].end();++p)
        (etaLin[k])[p->first] -= (p->second) * df(k,i);
    }
  }

  //**********************************************
  //   Lin N                                    //
  //**********************************************
  std::vector<GEN::pairedvector<int,double> > x_0Linnew(3,1000);
  std::vector<GEN::pairedvector<int,double> > normaltolineLinaux(3,1000);

  LINALG::Matrix<1,n_>   deriv;
  DRT::UTILS::shape_function_1D_deriv1(deriv,eta[0],distype);

  LINALG::Matrix<1,n_>   deriv2;
  DRT::UTILS::shape_function_1D_deriv2(deriv2,eta[0],distype);
  for (int i = 0;i<n_;++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("ERROR: Cannot find master node");
    MortarNode* mnode = dynamic_cast<MortarNode*>(node);

    for(int k = 0; k<2; ++k)
      (x_0Linnew[k])[mnode->Dofs()[k]] += deriv(i);

    for(int k = 0; k<2; ++k)
    {
      for (_CI p=etaLin[0].begin();p!=etaLin[0].end();++p)
        (x_0Linnew[k])[p->first] += (p->second) * deriv2(i) * mnode->xspatial()[k];
    }
  }

  // cross product linearization
  for (_CI p=x_0Linnew[1].begin();p!=x_0Linnew[1].end();++p)
    (normaltolineLinaux[0])[p->first] += (p->second);
  for (_CI p=x_0Linnew[0].begin();p!=x_0Linnew[0].end();++p)
    (normaltolineLinaux[1])[p->first] -= (p->second);

  //normalize lin
  LINALG::Matrix<3, 3>Wfinal;
//  const double lcubeinv = 1.0 / (normallength * normallength * normallength);

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      Wfinal(j, k) = -lcubeinv * linnormalaux[j] * linnormalaux[k];
      if (j == k)
        Wfinal(j, k) += 1 / normallength;
    }
  }

  // row loop
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p=normaltolineLinaux[0].begin();p!=normaltolineLinaux[0].end();++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j,0);

    for (_CI p=normaltolineLinaux[1].begin();p!=normaltolineLinaux[1].end();++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j,1);

    for (_CI p=normaltolineLinaux[2].begin();p!=normaltolineLinaux[2].end();++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j,2);
  }

  // bye bye
  return true;
}

/*----------------------------------------------------------------------*
 |  Project snode onto melement with master element normal   farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNormal(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist)
{
  if(ndim_==2)
  {
    ProjectSNodeByMNormal2D(
        snode,
        mele,
        xi,
        normal,
        dist);
  }
  else if(ndim_==3)
  {
    ProjectSNodeByMNormal3D(
        snode,
        mele,
        xi,
        normal,
        dist);
  }
  else
  {
    dserror("ERROR: wrong dimension!");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master element normal   farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNodalNormalLin(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  if(ndim_==2)
  {
    ProjectSNodeByMNodalNormal2DLin(
        snode,
        mele,
        xi,
        normal,
        dist,
        normaltolineLin);
  }
  else if(ndim_==3)
  {
    dserror("not yet implemented!");
  }
  else
  {
    dserror("ERROR: wrong dimension!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 01/16|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::ProjectSNodeByMNormalLin(
    MORTAR::MortarNode&    snode,
    MORTAR::MortarElement& mele,
    double* xi,
    double* normal,
    double& dist,
    std::vector<GEN::pairedvector<int,double> >& normaltolineLin)
{
  if(ndim_==2)
  {
    ProjectSNodeByMNormal2DLin(
        snode,
        mele,
        xi,
        normal,
        dist,
        normaltolineLin);
  }
  else if(ndim_==3)
  {
    dserror("ERROR: Not yet implemented!");
  }
  else
  {
    dserror("ERROR: wrong dimension!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for nodal normal case (public)                 popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateFNodalNormal(
    MORTAR::MortarNode& node, MORTAR::MortarElement& ele, const double* eta)
{
  /* Evaluate the function F(eta) = ( Ni * xim - xs ) x ns,
   or to be more precise the third component of this vector function!

   Ni  shape functions of element
   xim coords of element nodes (master side)
   xs  coords of node to be projected (slave side)
   ns   outward normal of node to be projected (slave side)          */
  double fval = 0.0;

  // build interpolation of master node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nx, 0);

  // subtract slave node coordinates
  for (int i = 0; i < ndim_; ++i)
    nx[i] -= node.xspatial()[i];

  //calculate F
  fval = nx[0] * node.MoData().n()[1] - nx[1] * node.MoData().n()[0];

  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for nodal normal case (public)             popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateGradFNodalNormal(
    MORTAR::MortarNode& node, MORTAR::MortarElement& ele, const double* eta)
{
  /* Evaluate the function GradF(eta)
   = Ni,eta * xim * nys - Ni,eta * yim * nxs,

   Ni,eta    shape function derivatives of element
   xim, yim  coords of element nodes (master side)
   nxs, nys   outward normal of node to be projected (slave side)   */

  double fgrad = 0.0;

  // build interpolation of master node coordinates for current eta
  // use shape function derivatives for interpolation (hence "1")
  double nxeta[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nxeta, 1);

  // calculate GradF
  fgrad = nxeta[0] * node.MoData().n()[1] - nxeta[1] * node.MoData().n()[0];

  return fgrad;

}

/*----------------------------------------------------------------------*
 |  Evaluate F for element normal case (public)               popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateFElementNormal(
    MORTAR::MortarNode& node, MORTAR::MortarElement& ele, const double* eta)
{
  /* Evaluate the function F(eta) = ( Ni * xis - xm ) x ( Nj * njs),
   or to be more precise the third component of this vector function!

   Ni  shape functions of element
   xis coords of element nodes (slave side)
   xm  coords of node to be projected (master side)
   nis outward normals of element nodes (slave side)                */

  double fval = 0.0;

  // collect necessary data (slave side)
  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes)
    dserror("ERROR: EvaluateFElementNormal: Null pointer!");

  LINALG::Matrix<n_, 1>      val;
  LINALG::Matrix<ndim_, n_>  coord;

  // get shape function values and derivatives at gpeta
  if (distype == DRT::Element::nurbs2 || distype == DRT::Element::nurbs3)
  {
    LINALG::SerialDenseVector auxval(n_);
    LINALG::SerialDenseMatrix deriv(n_, 1);
    ele.EvaluateShape(eta, auxval, deriv, ele.NumNode());

    for (int i = 0; i < n_; ++i)
      val(i) = auxval(i);
  }
  else
    DRT::UTILS::shape_function_1D(val, eta[0], distype);

  // get interpolated normal and proj. coordinates for current eta
  double nn[ndim_];
  double nx[ndim_];
  for (int j = 0; j < ndim_; ++j)
  {
    nn[j] = 0.0;
    nx[j] = 0.0;
  }

  for (int i = 0; i < n_; ++i)
  {
    MortarNode* mymrtrnode = dynamic_cast<MortarNode*>(mynodes[i]);

    for (int j = 0; j < ndim_; ++j)
    {
      nn[j] += val(i) * mymrtrnode->MoData().n()[j];

      coord(j, i) = mymrtrnode->xspatial()[j];

      nx[j] += val(i) * coord(j, i);
    }
  }

  // subtract master node coordinates
  for (int j = 0; j < ndim_; ++j)
    nx[j] -= node.xspatial()[j];

  // calculate F
  fval = nx[0] * nn[1] - nx[1] * nn[0];

  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for element normal case (public)           popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
double MORTAR::MortarProjectorCalc<distype>::EvaluateGradFElementNormal(
    MORTAR::MortarNode& node, MORTAR::MortarElement& ele, const double* eta)
{
  if (ndim_ == 3)
    dserror("This Projector is only for 2D Problems!");

    /* Evaluate the function GradF(eta)
     = ( Ni,eta * xis ) * ( Nj * nyjs )
     + ( Ni * xis - xm ) * ( Nj,eta * nyjs )
     - ( Ni,eta * yis ) * ( Nj * nxjs )
     - ( Ni * yis - ym ) * ( Nj,eta * nxjs )

     Ni,eta      shape function derivatives of element
     xis, yis   coords of element nodes (slave side)
     xm, ym     coords of node to be projected (master side)
     nxjs, nyjs outward normals of element nodes (slave side)         */

    double fgrad = 0.0;

    // collect necessary data (slave side)
    LINALG::Matrix<n_, 1> val;
    LINALG::Matrix<ndim_-1,n_> deriv;
    LINALG::Matrix<ndim_,n_> coord;

    DRT::Node** mynodes = ele.Nodes();
    if(!mynodes) dserror("ERROR: EvaluateGradFElementNormal: Null pointer!");

    // get shape function values and derivatives at gpeta
    if(distype==DRT::Element::nurbs2 || distype==DRT::Element::nurbs3)
    {
      LINALG::SerialDenseVector auxval(n_);
      LINALG::SerialDenseMatrix auxderiv(n_,1);
      ele.EvaluateShape(eta, auxval, auxderiv, ele.NumNode());

      for(int i=0;i<n_;++i)
      {
        val(i)=auxval(i);
        deriv(0,i)=auxderiv(i,0);
      }
    }
    else
    {
      DRT::UTILS::shape_function_1D (val,eta[0],distype);
      DRT::UTILS::shape_function_1D_deriv1 (deriv,eta[0],distype);
    }

    // get interpolated normal and proj. coordinates for current eta
    double nn[ndim_];
    double nneta[ndim_];
    double nx[ndim_];
    double nxeta[ndim_];
    for (int j=0;j<ndim_;++j)
    {
      nn[j]=0.0;
      nneta[j]=0.0;
      nx[j]=0.0;
      nxeta[j]=0.0;
    }

    for (int i=0;i<n_;++i)
    {
      MortarNode* mymrtrnode = dynamic_cast<MortarNode*> (mynodes[i]);

      for (int j=0;j<ndim_;++j)
      {
        nn[j] += val(i)*mymrtrnode->MoData().n()[j];
        nneta[j] += deriv(0,i)*mymrtrnode->MoData().n()[j];
        coord(j,i) = mymrtrnode->xspatial()[j];
        nx[j] += val(i)*coord(j,i);
        nxeta[j] += deriv(0,i)*coord(j,i);
      }
    }

    // subtract master node coordinates
    for(int j=0;j<ndim_;++j)
    nx[j]-=node.xspatial()[j];

    // calculate GradF
    fgrad = nxeta[0]*nn[1] + nx[0]*nneta[1]
    - nxeta[1]*nn[0] - nx[1]*nneta[0];

    return fgrad;
  }

    /*----------------------------------------------------------------------*
     |  Evaluate F for Gauss point case (public)                  popp 01/08|
     *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
double MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::EvaluateFGaussPoint2D(
    const double* gpx, const double* gpn, MORTAR::MortarElement& ele,
    const double* eta)
{
  /* Evaluate the function F(eta) = ( Ni * xim - gpx ) x gpn,
   or to be more precise the third component of this vector function!

   Ni  shape functions of element (master side)
   xim coords of element nodes (master side)
   gpx coords of GP to be projected (slave side)
   gpn outward normal of GP to be projected (slave side)          */

  double fval = 0.0;

  // build interpolation of master node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nx, 0);

  // subtract GP coordinates
  nx[0] -= gpx[0];
  nx[1] -= gpx[1];

  // calculate F
  fval = nx[0] * gpn[1] - nx[1] * gpn[0];

  return fval;

}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (public)              popp 01/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
double MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::EvaluateGradFGaussPoint2D(
    const double* gpn, MORTAR::MortarElement& ele, const double* eta)
{
  /* Evaluate the function GradF(eta)
   = Ni,eta * xim * gpny - Ni,eta * yim * gpnx,

   Ni,eta     shape function derivatives of element (master side)
   xim, yim   coords of element nodes (master side)
   gpnx, gpny outward normal of GP to be projected (slave side)   */

  double fgrad = 0.0;

  // build interpolation of master node coordinates for current eta
  // use shape function derivatives for interpolation (hence "1")
  double nxeta[ndim_];
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nxeta, 1);

  // calculate GradF
  fgrad = nxeta[0] * gpn[1] - nxeta[1] * gpn[0];

  return fgrad;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (3D)                      popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
bool MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::EvaluateFGaussPoint3D(
    double* f, const double* gpx, const double* gpn, MORTAR::MortarElement& ele,
    const double* eta, const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi - alpha * gpn - gpx
   which is a vector-valued function with 3 components!

   Ni      shape functions of element to project on
   xi      coords of nodes of element to project on
   gpx     coords of GP to be projected
   gpn     normal of GP along which to project                  */

  // build interpolation of ele node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nx, 0);

  // evaluate function f
  for (int i = 0; i < ndim_; ++i)
    f[i] = nx[i] - alpha * gpn[i] - gpx[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (3D)                  popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distypeS,
    DRT::Element::DiscretizationType distypeM>
bool MORTAR::MortarProjectorCalc_EleBased<distypeS, distypeM>::EvaluateGradFGaussPoint3D(
    LINALG::Matrix<3, 3>& fgrad,
    const double* gpx,
    const double* gpn,
    MORTAR::MortarElement& ele,
    const double* eta,
    const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
   - alpha * gpn - gpx, which is a (3x3)-matrix!

   Ni      shape functions of element to project on
   xi      coords of nodes of element to project on
   gpx     coords of GP to be projected
   gpn     normal of GP along which to project                  */

  // build interpolation of ele node coordinates for current eta
//  double nxeta1[3] = {0.0, 0.0, 0.0};
//  double nxeta2[3] = {0.0, 0.0, 0.0};
//  ele.LocalToGlobal(eta,nxeta1,1);
//  ele.LocalToGlobal(eta,nxeta2,2);
  double nxeta1[ndim_];
  double nxeta2[ndim_];
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nxeta1, 1);
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nxeta2, 2);

  //evaluate function f gradient
  for (int i = 0; i < ndim_; ++i)
  {
    fgrad(i, 0) = nxeta1[i];
    fgrad(i, 1) = nxeta2[i];
    fgrad(i, 2) = -gpn[i];
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for AuxPlane Gauss point case (3D)             popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::EvaluateFGaussPointAuxn3D(double* f,
    const double* globgp, const double* auxn, MORTAR::MortarElement& ele,
    const double* eta, const double& alpha)
{
  /* Evaluate the function F(eta,alpha) = Ni * xi - alpha * auxn - globgp
   which is a vector-valued function with 3 components!

   Ni      shape functions of element to project on
   xi      coords of nodes of element to project on
   globgp  coords of AuxPlaneGP to be projected
   auxn    normal of AuxPlane along which to project            */

  // build interpolation of ele node coordinates for current eta
  double nx[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nx, 0);

  // evaluate function f
  for (int i = 0; i < ndim_; ++i)
    f[i] = nx[i] - alpha * auxn[i] - globgp[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for AuxPlane Gauss point case (3D)         popp 11/08|
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool MORTAR::MortarProjectorCalc<distype>::EvaluateGradFGaussPointAuxn3D(
    LINALG::Matrix<3, 3>& fgrad,
    const double* globgp,
    const double* auxn,
    MORTAR::MortarElement& ele,
    const double* eta,
    const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
   - alpha * auxn - globgp, which is a (3x3)-matrix!

   Ni      shape functions of element to project on
   xi      coords of nodes of element to project on
   globgp  coords of AuxPlaneGP to be projected
   auxn    normal of AuxPlane along which to project            */

  // build interpolation of ele node coordinates for current eta
  double nxeta1[ndim_];
  double nxeta2[ndim_];
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nxeta1, 1);
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nxeta2, 2);

  // evaluate function f gradient
  for (int i = 0; i < ndim_; ++i)
  {
    fgrad(i, 0) = nxeta1[i];
    fgrad(i, 1) = nxeta2[i];
    fgrad(i, 2) = -auxn[i];
  }

  return true;
}

