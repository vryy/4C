/*-----------------------------------------------------------------------*/
/*! \file
\brief A class to perform projections of nodes onto opposing MORTAR::Elements

\level 1

*/
/*-----------------------------------------------------------------------*/

#include "baci_mortar_projector.hpp"

#include "baci_contact_element.hpp"
#include "baci_contact_node.hpp"
#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_linalg_serialdensematrix.hpp"
#include "baci_linalg_serialdensevector.hpp"
#include "baci_linalg_utils_sparse_algebra_math.hpp"
#include "baci_mortar_calc_utils.hpp"
#include "baci_mortar_defines.hpp"
#include "baci_mortar_element.hpp"
#include "baci_mortar_interface.hpp"
#include "baci_mortar_node.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  impl. for aux.-plane based projection                    farah 01/14|
 *----------------------------------------------------------------------*/
MORTAR::Projector* MORTAR::Projector::Impl(MORTAR::Element& ele)
{
  switch (ele.Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      return ProjectorCalc<CORE::FE::CellType::quad4>::Instance();
    }
    case CORE::FE::CellType::quad8:
    {
      return ProjectorCalc<CORE::FE::CellType::quad8>::Instance();
    }
    case CORE::FE::CellType::quad9:
    {
      return ProjectorCalc<CORE::FE::CellType::quad9>::Instance();
    }
    case CORE::FE::CellType::tri3:
    {
      return ProjectorCalc<CORE::FE::CellType::tri3>::Instance();
    }
    case CORE::FE::CellType::tri6:
    {
      return ProjectorCalc<CORE::FE::CellType::tri6>::Instance();
    }
    case CORE::FE::CellType::line2:
    {
      return ProjectorCalc<CORE::FE::CellType::line2>::Instance();
    }
    case CORE::FE::CellType::line3:
    {
      return ProjectorCalc<CORE::FE::CellType::line3>::Instance();
    }
      //==================================================
      //                     NURBS
      //==================================================
    case CORE::FE::CellType::nurbs2:
    {
      return ProjectorCalc<CORE::FE::CellType::nurbs2>::Instance();
    }
    case CORE::FE::CellType::nurbs3:
    {
      return ProjectorCalc<CORE::FE::CellType::nurbs3>::Instance();
    }
    case CORE::FE::CellType::nurbs4:
    {
      return ProjectorCalc<CORE::FE::CellType::nurbs4>::Instance();
    }
    case CORE::FE::CellType::nurbs8:
    {
      return ProjectorCalc<CORE::FE::CellType::nurbs8>::Instance();
    }
    case CORE::FE::CellType::nurbs9:
    {
      return ProjectorCalc<CORE::FE::CellType::nurbs9>::Instance();
    }
    default:
      dserror("Element shape %d (%d nodes) not activated. Just do it.", ele.Shape(), ele.NumNode());
      break;
  }
  return nullptr;
}

/*----------------------------------------------------------------------*
 |  impl. for element based projection                       farah 04/14|
 *----------------------------------------------------------------------*/
MORTAR::Projector* MORTAR::Projector::Impl(MORTAR::Element& sele, MORTAR::Element& mele)
{
  switch (sele.Shape())
  {
    case CORE::FE::CellType::quad4:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad4,
              CORE::FE::CellType::quad4>::Instance();
        }
        case CORE::FE::CellType::quad8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad4,
              CORE::FE::CellType::quad8>::Instance();
        }
        case CORE::FE::CellType::quad9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad4,
              CORE::FE::CellType::quad9>::Instance();
        }
        case CORE::FE::CellType::tri3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad4,
              CORE::FE::CellType::tri3>::Instance();
        }
        case CORE::FE::CellType::tri6:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad4,
              CORE::FE::CellType::tri6>::Instance();
        }
        case CORE::FE::CellType::nurbs9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad4,
              CORE::FE::CellType::nurbs9>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::quad8:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad8,
              CORE::FE::CellType::quad4>::Instance();
        }
        case CORE::FE::CellType::quad8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad8,
              CORE::FE::CellType::quad8>::Instance();
        }
        case CORE::FE::CellType::quad9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad8,
              CORE::FE::CellType::quad9>::Instance();
        }
        case CORE::FE::CellType::tri3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad8,
              CORE::FE::CellType::tri3>::Instance();
        }
        case CORE::FE::CellType::tri6:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad8,
              CORE::FE::CellType::tri6>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::quad9:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad9,
              CORE::FE::CellType::quad4>::Instance();
        }
        case CORE::FE::CellType::quad8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad9,
              CORE::FE::CellType::quad8>::Instance();
        }
        case CORE::FE::CellType::quad9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad9,
              CORE::FE::CellType::quad9>::Instance();
        }
        case CORE::FE::CellType::tri3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad9,
              CORE::FE::CellType::tri3>::Instance();
        }
        case CORE::FE::CellType::tri6:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::quad9,
              CORE::FE::CellType::tri6>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::tri3:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri3,
              CORE::FE::CellType::quad4>::Instance();
        }
        case CORE::FE::CellType::quad8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri3,
              CORE::FE::CellType::quad8>::Instance();
        }
        case CORE::FE::CellType::quad9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri3,
              CORE::FE::CellType::quad9>::Instance();
        }
        case CORE::FE::CellType::tri3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri3,
              CORE::FE::CellType::tri3>::Instance();
        }
        case CORE::FE::CellType::tri6:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri3,
              CORE::FE::CellType::tri6>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::tri6:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::quad4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri6,
              CORE::FE::CellType::quad4>::Instance();
        }
        case CORE::FE::CellType::quad8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri6,
              CORE::FE::CellType::quad8>::Instance();
        }
        case CORE::FE::CellType::quad9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri6,
              CORE::FE::CellType::quad9>::Instance();
        }
        case CORE::FE::CellType::tri3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri6,
              CORE::FE::CellType::tri3>::Instance();
        }
        case CORE::FE::CellType::tri6:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::tri6,
              CORE::FE::CellType::tri6>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::line2:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::line2:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::line2,
              CORE::FE::CellType::line2>::Instance();
        }
        case CORE::FE::CellType::line3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::line2,
              CORE::FE::CellType::line3>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::line3:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::line2:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::line3,
              CORE::FE::CellType::line2>::Instance();
        }
        case CORE::FE::CellType::line3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::line3,
              CORE::FE::CellType::line3>::Instance();
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
    case CORE::FE::CellType::nurbs2:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs2:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs2,
              CORE::FE::CellType::nurbs2>::Instance();
        }
        case CORE::FE::CellType::nurbs3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs2,
              CORE::FE::CellType::nurbs3>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::nurbs3:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs2:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs3,
              CORE::FE::CellType::nurbs2>::Instance();
        }
        case CORE::FE::CellType::nurbs3:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs3,
              CORE::FE::CellType::nurbs3>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::nurbs4:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs4,
              CORE::FE::CellType::nurbs4>::Instance();
        }
        case CORE::FE::CellType::nurbs8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs4,
              CORE::FE::CellType::nurbs8>::Instance();
        }
        case CORE::FE::CellType::nurbs9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs4,
              CORE::FE::CellType::nurbs9>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::nurbs8:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs8,
              CORE::FE::CellType::nurbs4>::Instance();
        }
        case CORE::FE::CellType::nurbs8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs8,
              CORE::FE::CellType::nurbs8>::Instance();
        }
        case CORE::FE::CellType::nurbs9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs8,
              CORE::FE::CellType::nurbs9>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    case CORE::FE::CellType::nurbs9:
    {
      switch (mele.Shape())
      {
        case CORE::FE::CellType::nurbs4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs9,
              CORE::FE::CellType::nurbs4>::Instance();
        }
        case CORE::FE::CellType::nurbs8:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs9,
              CORE::FE::CellType::nurbs8>::Instance();
        }
        case CORE::FE::CellType::nurbs9:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs9,
              CORE::FE::CellType::nurbs9>::Instance();
        }
        case CORE::FE::CellType::quad4:
        {
          return ProjectorCalcEleBased<CORE::FE::CellType::nurbs9,
              CORE::FE::CellType::quad4>::Instance();
        }
        default:
          dserror("Element shape not supported!");
          break;
      }
      break;
    }
    default:
      dserror(
          "Element shape %d (%d nodes) not activated. Just do it.", sele.Shape(), sele.NumNode());
      break;
  }
  return nullptr;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            farah 01/14|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
MORTAR::ProjectorCalc<distype>::ProjectorCalc()
{
  // nothing
}

/*----------------------------------------------------------------------*
 |  ctor ele-based (public)                                  farah 04/14|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::ProjectorCalcEleBased()
{
  // nothing
}

template <CORE::FE::CellType distype>
MORTAR::ProjectorCalc<distype>* MORTAR::ProjectorCalc<distype>::Instance(
    CORE::UTILS::SingletonAction action)
{
  static auto singleton_owner = CORE::UTILS::MakeSingletonOwner(
      []() {
        return std::unique_ptr<MORTAR::ProjectorCalc<distype>>(
            new MORTAR::ProjectorCalc<distype>());
      });

  return singleton_owner.Instance(action);
}

template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
MORTAR::ProjectorCalcEleBased<distypeS, distypeM>*
MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::Instance(CORE::UTILS::SingletonAction action)
{
  static CORE::UTILS::SingletonOwner<MORTAR::ProjectorCalcEleBased<distypeS, distypeM>>
      singleton_owner(
          []()
          {
            return std::unique_ptr<MORTAR::ProjectorCalcEleBased<distypeS, distypeM>>(
                new MORTAR::ProjectorCalcEleBased<distypeS, distypeM>());
          });

  return singleton_owner.Instance(action);
}


/*----------------------------------------------------------------------*
 |  Project a node along its nodal normal (public)            popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectNodalNormal(
    MORTAR::Node& node, MORTAR::Element& ele, double* xi)
{
  bool ok = true;
  if (ndim_ == 2)
  {
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double f = 0.0;
    double df = 0.0;
    int k = 0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      f = EvaluateFNodalNormal(node, ele, eta);
      if (abs(f) < MORTARCONVTOL) break;
      df = EvaluateGradFNodalNormal(node, ele, eta);
      if (abs(df) < 1.0e-12)
      {
        xi[0] = 1.0e12;
        return false;
        dserror("Singular Jacobian for projection");
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
      // std::cout << "***WARNING*** ProjectNodalNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and MORTAR::ElementID " << ele.Id() << std::endl;
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
    dserror("ProjectNodalNormal: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a node along element's normal field (public)      popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectElementNormal(
    MORTAR::Node& node, MORTAR::Element& ele, double* xi)
{
  bool ok = true;
  if (ndim_ == 2)
  {
    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double f = 0.0;
    double df = 0.0;
    int k = 0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      f = EvaluateFElementNormal(node, ele, eta);
      if (abs(f) < MORTARCONVTOL) break;
      df = EvaluateGradFElementNormal(node, ele, eta);
      if (abs(df) < 1.0e-12)
      {
        ok = false;
        xi[0] = 1.0e12;
        return ok;
        dserror("Singular Jacobian for projection");
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
      // std::cout << "***WARNING*** ProjectElementNormal:" << " Newton unconverged for NodeID "
      //     << node.Id() << " and MORTAR::ElementID " << ele.Id() << std::endl;
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
    dserror("ProjectElementNormal: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (public)           popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
bool MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::ProjectGaussPoint2D(
    MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele, double* xi)
{
  bool ok = true;
  if (ndim_ == 2)
  {
    CORE::LINALG::Matrix<ns_, 1> val;
    CORE::LINALG::Matrix<ndim_, ns_> coord;

    DRT::Node** mynodes = gpele.Nodes();
    if (!mynodes) dserror("ProjectGaussPoint: Null pointer!");

    // get shape function values and derivatives at gpeta
    if (distypeS == CORE::FE::CellType::nurbs2 || distypeS == CORE::FE::CellType::nurbs3)
    {
      CORE::LINALG::SerialDenseVector auxval(ns_);
      CORE::LINALG::SerialDenseMatrix deriv(ns_, 1);
      gpele.EvaluateShape(gpeta, auxval, deriv, gpele.NumNode());

      for (int i = 0; i < ns_; ++i) val(i) = auxval(i);
    }
    else
      CORE::FE::shape_function_1D(val, gpeta[0], distypeS);

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
      Node* mymrtrnode = dynamic_cast<Node*>(mynodes[i]);

      for (int j = 0; j < ndim_; ++j)
      {
        gpn[j] += val(i) * mymrtrnode->MoData().n()[j];

        coord(j, i) = mymrtrnode->xspatial()[j];

        gpx[j] += val(i) * coord(j, i);
      }
    }

    // local Newton iteration for xi, start in the element middle
    double eta[2] = {0.0, 0.0};
    double f = 0.0;
    double df = 0.0;
    int k = 0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      f = EvaluateFGaussPoint2D(gpx, gpn, ele, eta);
      if (abs(f) < MORTARCONVTOL) break;
      df = EvaluateGradFGaussPoint2D(gpn, ele, eta);
      if (abs(df) < 1.0e-12) dserror("Singular Jacobian for projection");
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

      //      dserror("ProjectGaussPoint: Newton unconverged for GP at xi=%d"
      //          " from MORTAR::ElementID %i", gpeta[0], gpele.Id());
    }
  }

  else
    dserror("ProjectGaussPoint: Called 2D version for 3D problem!");

  return ok;
}

/*----------------------------------------------------------------------*
 | Check projection for warped elements quad4 elements       farah 01/13|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
bool MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::CheckProjection4AUXPLANE(
    MORTAR::Element& ele, double* ngp, double* globgp)
{
  if (ele.Shape() == CORE::FE::CellType::tri3) dserror("ELEMENT SHAPE TRI3 -- NO WARPING");

  if (ele.Shape() != CORE::FE::CellType::quad4)
  {
    return true;
  }

  int nnode = ele.NumNode();
  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("Project: Null pointer!");

  // compute base-vectors
  std::vector<double> t0(3);
  std::vector<double> t1(3);
  std::vector<double> auxn(3);
  std::vector<double> auxc(3);
  std::vector<double> proj_gp(3);
  CORE::LINALG::Matrix<3, 3> P;
  CORE::LINALG::Matrix<3, 3> T;
  std::array<double, 3> n = {0.0, 0.0, 0.0};
  double length_t = 0.0;
  double length_n = 0.0;
  double a1 = 0.0;
  bool all_negative = true;
  Node* mycnode_0 = nullptr;
  Node* mycnode_1 = nullptr;
  Node* mycnode_2 = nullptr;

  // project gp onto auxn
  for (int i = 0; i < nnode; ++i)  // loop over edges
  {
    if (i == 0)
    {
      mycnode_1 = dynamic_cast<Node*>(mynodes[0]);
      if (!mycnode_1) dserror("Project: Null pointer!");

      mycnode_0 = dynamic_cast<Node*>(mynodes[3]);
      if (!mycnode_0) dserror("Project: Null pointer!");

      mycnode_2 = dynamic_cast<Node*>(mynodes[1]);
      if (!mycnode_2) dserror("Project: Null pointer!");
    }
    if (i == 3)
    {
      mycnode_1 = dynamic_cast<Node*>(mynodes[3]);
      if (!mycnode_1) dserror("Project: Null pointer!");

      mycnode_0 = dynamic_cast<Node*>(mynodes[2]);
      if (!mycnode_0) dserror("Project: Null pointer!");

      mycnode_2 = dynamic_cast<Node*>(mynodes[0]);
      if (!mycnode_2) dserror("Project: Null pointer!");
    }
    if (i == 1 || i == 2)
    {
      mycnode_1 = dynamic_cast<Node*>(mynodes[i]);
      if (!mycnode_1) dserror("Project: Null pointer!");

      mycnode_0 = dynamic_cast<Node*>(mynodes[i - 1]);
      if (!mycnode_0) dserror("Project: Null pointer!");

      mycnode_2 = dynamic_cast<Node*>(mynodes[i + 1]);
      if (!mycnode_2) dserror("Project: Null pointer!");
    }

    // span vectors -- edges
    t0[0] = mycnode_0->xspatial()[0] - mycnode_1->xspatial()[0];
    t0[1] = mycnode_0->xspatial()[1] - mycnode_1->xspatial()[1];
    t0[2] = mycnode_0->xspatial()[2] - mycnode_1->xspatial()[2];

    t1[0] = mycnode_2->xspatial()[0] - mycnode_1->xspatial()[0];
    t1[1] = mycnode_2->xspatial()[1] - mycnode_1->xspatial()[1];
    t1[2] = mycnode_2->xspatial()[2] - mycnode_1->xspatial()[2];

    auxc[0] = mycnode_1->xspatial()[0];
    auxc[1] = mycnode_1->xspatial()[1];
    auxc[2] = mycnode_1->xspatial()[2];

    auxn[0] = t1[1] * t0[2] - t1[2] * t0[1];
    auxn[1] = t1[2] * t0[0] - t1[0] * t0[2];
    auxn[2] = t1[0] * t0[1] - t1[1] * t0[0];

    // fill Projection matrix P
    for (int j = 0; j < 3; ++j) P(j, 2) = -ngp[j];

    for (int j = 0; j < 3; ++j) P(j, 0) = t0[j];

    for (int j = 0; j < 3; ++j) P(j, 1) = t1[j];

    P.Invert();
    double lambda1 = P(0, 0) * (globgp[0] - auxc[0]) + P(0, 1) * (globgp[1] - auxc[1]) +
                     P(0, 2) * (globgp[2] - auxc[2]);
    double lambda2 = P(1, 0) * (globgp[0] - auxc[0]) + P(1, 1) * (globgp[1] - auxc[1]) +
                     P(1, 2) * (globgp[2] - auxc[2]);
    // double alph= P(2,0)*(globgp[0]-auxc[0]) + P(2,1)*(globgp[1]-auxc[1]) +
    // P(2,2)*(globgp[2]-auxc[2]);

    proj_gp[0] = lambda1 * t0[0] + lambda2 * t1[0] + auxc[0];
    proj_gp[1] = lambda1 * t0[1] + lambda2 * t1[1] + auxc[1];
    proj_gp[2] = lambda1 * t0[2] + lambda2 * t1[2] + auxc[2];
    // std::cout << "proj_gp_AUX4PLANE= " << proj_gp[0] << " "  << proj_gp[1] << " " << proj_gp[2]
    // << std::endl;

    // check
    // edge 1
    n[0] = -(t0[1] * auxn[2] - t0[2] * auxn[1]);
    n[1] = -(t0[2] * auxn[0] - t0[0] * auxn[2]);
    n[2] = -(t0[0] * auxn[1] - t0[1] * auxn[0]);

    length_t = sqrt(t0[0] * t0[0] + t0[1] * t0[1] + t0[2] * t0[2]);
    length_n = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    // fill Projection matrix T
    for (int j = 0; j < 3; ++j) T(j, 0) = n[j] / length_n;

    for (int j = 0; j < 3; ++j) T(j, 1) = t0[j] / length_t;

    for (int j = 0; j < 3; ++j) T(j, 2) = auxc[j];

    T.Invert();
    a1 = T(0, 0) * proj_gp[0] + T(0, 1) * proj_gp[1] + T(0, 2) * proj_gp[2];

    if (a1 > 0.0) all_negative = false;

    // edge 2
    n[0] = (t1[1] * auxn[2] - t1[2] * auxn[1]);
    n[1] = (t1[2] * auxn[0] - t1[0] * auxn[2]);
    n[2] = (t1[0] * auxn[1] - t1[1] * auxn[0]);

    length_t = sqrt(t1[0] * t1[0] + t1[1] * t1[1] + t1[2] * t1[2]);
    length_n = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    // fill Projection matrix T
    for (int j = 0; j < 3; ++j) T(j, 0) = n[j] / length_n;

    for (int j = 0; j < 3; ++j) T(j, 1) = t1[j] / length_t;

    for (int j = 0; j < 3; ++j) T(j, 2) = auxc[j];

    T.Invert();

    a1 = T(0, 0) * proj_gp[0] + T(0, 1) * proj_gp[1] + T(0, 2) * proj_gp[2];
    // a2=T(1,0)*proj_gp[0] + T(1,1)*proj_gp[1] + T(1,2)*proj_gp[2];

    if (a1 > 0.0) all_negative = false;
  }

  if (all_negative == false) return true;

  return false;  // creates error in ProjectGaussPoint3D()
}

/*----------------------------------------------------------------------*
 |  Project a Gauss point along its normal (3D)               popp 11/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
bool MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::ProjectGaussPoint3D(
    MORTAR::Element& gpele, const double* gpeta, MORTAR::Element& ele, double* xi, double& par)
{
  if (ndim_ == 3)
  {
    CORE::LINALG::Matrix<ns_, 1> val;
    CORE::LINALG::Matrix<ndim_, ns_> coord;
    coord.Clear();

    DRT::Node** mypoints = gpele.Points();
    DRT::Node** mynodes = gpele.Nodes();
    if (!mypoints) dserror("ProjectGaussPoint: Null pointer!");

    // get shape function values and derivatives at gpeta
    if (distypeS == CORE::FE::CellType::nurbs4 || distypeS == CORE::FE::CellType::nurbs8 ||
        distypeS == CORE::FE::CellType::nurbs9)
    {
      CORE::LINALG::SerialDenseVector auxval(ns_);
      CORE::LINALG::SerialDenseMatrix deriv(ns_, 2);
      gpele.EvaluateShape(gpeta, auxval, deriv, gpele.NumNode());

      for (int i = 0; i < ns_; ++i) val(i) = auxval(i);
    }
    else
    {
      CORE::FE::shape_function_2D(val, gpeta[0], gpeta[1], distypeS);
    }

    // get interpolated GP normal and GP coordinates
    double gpn[3];
    double gpx[3];
    for (int i = 0; i < ndim_; ++i)
    {
      gpn[i] = 0.0;
      gpx[i] = 0.0;
    }

    for (int i = 0; i < ns_; ++i)
    {
      Node* mymrtrnode = dynamic_cast<Node*>(mypoints[i]);

      for (int j = 0; j < ndim_; ++j)
      {
        coord(j, i) = mymrtrnode->xspatial()[j];

        gpx[j] += val(i) * coord(j, i);
      }
    }

    for (int i = 0; i < gpele.NumNode(); ++i)
    {
      Node* mymrtrnode = dynamic_cast<Node*>(mynodes[i]);
      for (int j = 0; j < ndim_; ++j)
      {
        gpn[j] += val(i) * mymrtrnode->MoData().n()[j];
      }
    }

    // start in the element center
    CORE::FE::CellType dt = ele.Shape();
    double eta[2] = {0.0, 0.0};

    if (dt == CORE::FE::CellType::tri3 || dt == CORE::FE::CellType::tri6)
    {
      eta[0] = 1.0 / 3.0;
      eta[1] = 1.0 / 3.0;
    }

    // auxiliary variable
    double alpha = 0.0;

    // function f (vector-valued)
    double f[3] = {0.0, 0.0, 0.0};

    // gradient of f (df/deta[0], df/deta[1], df/dalpha)
    CORE::LINALG::Matrix<3, 3> df;
    // start iteration
    int k = 0;
    double conv = 0.0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      EvaluateFGaussPoint3D(f, gpx, gpn, ele, eta, alpha);
      conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
      if (conv <= MORTARCONVTOL) break;
      EvaluateGradFGaussPoint3D(df, gpx, gpn, ele, eta, alpha);

      // safety check: if projection normal is parallel to the master element --> det can be zero
      double det = df.Determinant();
      if (det > -1e-12 and det < 1e-12)
      {
        std::cout << "WARNING: GPProjection parallel to master element --> GP skipped for this "
                     "master element!"
                  << std::endl;
        // leave here
        return false;
      }

      // solve deta = - inv(df) * f
      double jacdet = df.Invert();
      if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

      // update eta and alpha
      eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
      eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
      alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];

      // Projection Check
      //      if (k==MORTARMAXITER-1)
      //      {
      //        bool check = CheckProjection4AUXPLANE(ele, gpn,gpx);
      //        if (check==false)
      //          dserror("!!! STOP !!!   -->   Projection Error: Newton unconverged but GP on mele
      //          !!!");
      //      }
    }

    // Newton iteration unconverged
    if (conv > MORTARCONVTOL)
    {
      xi[0] = 1e12;
      xi[1] = 1e12;
      par = alpha;

      return false;
    }
    else
    {
      xi[0] = eta[0];
      xi[1] = eta[1];
      par = alpha;
    }
  }

  else
    dserror("ProjectGaussPoint: Called 3D version for 2D problem!");

  return true;
}
/*----------------------------------------------------------------------*
 |  Project a Gauss point along AuxPlane normal (3D)          popp 11/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectGaussPointAuxn3D(
    const double* globgp, const double* auxn, MORTAR::Element& ele, double* xi, double& par)
{
  if (ndim_ == 3)
  {
    // start in the element center
    CORE::FE::CellType dt = ele.Shape();
    double eta[2] = {0.0, 0.0};

    switch (dt)
    {
      case CORE::FE::CellType::tri3:
      case CORE::FE::CellType::tri6:
      {
        eta[0] = 1.0 / 3.0;
        eta[1] = 1.0 / 3.0;

        break;
      }
      default:
        break;
    }

    // auxiliary variable
    double& alpha = par;
    alpha = 0.0;

    // function f (vector-valued)
    double f[3] = {0.0, 0.0, 0.0};

    // gradient of f (df/deta[0], df/deta[1], df/dalpha)
    CORE::LINALG::Matrix<3, 3> df;

    // start iteration
    int k = 0;
    double conv = 0.0;

    for (k = 0; k < MORTARMAXITER; ++k)
    {
      EvaluateFGaussPointAuxn3D(f, globgp, auxn, ele, eta, alpha);
      conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
      // std::cout << "Iteration " << k << ": -> |f|=" << conv << std::endl;
      if (conv <= MORTARCONVTOL) break;
      EvaluateGradFGaussPointAuxn3D(df, globgp, auxn, ele, eta, alpha);

      // solve deta = - inv(df) * f
      double jacdet = df.Invert();
      if (abs(jacdet) < 1.0e-12)
      {
        dserror("Singular Jacobian for projection");
      }

      // update eta and alpha
      eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
      eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
      alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
    }

    // Newton iteration unconverged
    if (conv > MORTARCONVTOL)
    {
      //      std::cout << "res= " << conv << std::endl;
      //      dserror("ProjectGaussPointAuxn3D: Newton unconverged for GP"
      //          "at xi = (%f,%f,%f) onto MORTAR::ElementID %i", globgp[0], globgp[1],
      //          globgp[2], ele.Id());
      xi[0] = 1e12;
      xi[1] = 1e12;
      return false;
    }

    // Newton iteration converged
    else
    {
      xi[0] = eta[0];
      xi[1] = eta[1];
    }

    // std::cout << "Newton iteration converged in " << k << " steps!" << std::endl;
  }

  else
    dserror("ProjectGaussPoint: Called 3D version for 2D problem!");

  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNormal3D(
    MORTAR::Node& snode, MORTAR::Element& mele, double* xi, double* normal, double& dist)
{
  if (ndim_ != 3) dserror("ProjectSNodeByMNormal3D is only for 3D problems!");

  // start in the element center
  double eta[2] = {0.0, 0.0};
  switch (distype)
  {
    case CORE::FE::CellType::tri3:
    case CORE::FE::CellType::tri6:
    {
      eta[0] = 1.0 / 3.0;
      eta[1] = 1.0 / 3.0;
      break;
    }
    default:
      break;
  }

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  std::array<double, 3> f = {0.0, 0.0, 0.0};

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  CORE::LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};
    double unormal[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 0.0};
    double normalpart[3] = {0.0, 0.0, 0.0};

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta, unormal);
    for (int i = 0; i < 3; ++i) normalpart[i] = unormal[i] * alpha;

    for (int i = 0; i < 3; ++i) normal[i] = unormal[i] * length;

    // calc xslave
    for (int i = 0; i < 3; ++i) xs[i] = snode.xspatial()[i];

    // calculate F
    for (int i = 0; i < 3; ++i) f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL) break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // master coordinate grad
    double meta0[3] = {0.0, 0.0, 0.0};  // x,xi_0
    double meta1[3] = {0.0, 0.0, 0.0};  // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta1, 2);

    // normal grad
    CORE::LINALG::Matrix<3, n_> secderiv;
    std::array<double, 3> meta00 = {0.0, 0.0, 0.0};  // x , xi_0 xi_0
    std::array<double, 3> meta11 = {0.0, 0.0, 0.0};  // x , xi_1 xi_1
    std::array<double, 3> meta01 = {0.0, 0.0, 0.0};  // x , xi_0 xi_1

    CORE::FE::shape_function_2D_deriv2(secderiv, eta[0], eta[1], distype);

    for (int i = 0; i < n_; ++i)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[i]);
      if (!mymnode) dserror("Null pointer!");
      for (int d = 0; d < 3; ++d)
      {
        meta00[d] += secderiv(0, i) * mymnode->xspatial()[d];
        meta11[d] += secderiv(1, i) * mymnode->xspatial()[d];
        meta01[d] += secderiv(2, i) * mymnode->xspatial()[d];
      }
    }

    std::array<double, 3> naux_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> naux_1 = {0.0, 0.0, 0.0};

    // normal grad xi_0
    naux_0[0] = (meta00[1] * meta1[2] - meta00[2] * meta1[1]);
    naux_0[1] = (meta00[2] * meta1[0] - meta00[0] * meta1[2]);
    naux_0[2] = (meta00[0] * meta1[1] - meta00[1] * meta1[0]);

    naux_0[0] += (meta0[1] * meta01[2] - meta0[2] * meta01[1]);
    naux_0[1] += (meta0[2] * meta01[0] - meta0[0] * meta01[2]);
    naux_0[2] += (meta0[0] * meta01[1] - meta0[1] * meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1] * meta1[2] - meta01[2] * meta1[1]);
    naux_1[1] = (meta01[2] * meta1[0] - meta01[0] * meta1[2]);
    naux_1[2] = (meta01[0] * meta1[1] - meta01[1] * meta1[0]);

    naux_1[0] += (meta0[1] * meta11[2] - meta0[2] * meta11[1]);
    naux_1[1] += (meta0[2] * meta11[0] - meta0[0] * meta11[2]);
    naux_1[2] += (meta0[0] * meta11[1] - meta0[1] * meta11[0]);

    std::array<double, 3> n_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> n_1 = {0.0, 0.0, 0.0};

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i = 0; i < 3; ++i) fac0 += naux_0[i] * normal[i];

    for (int i = 0; i < 3; ++i)
      n_0[i] = naux_0[i] / length - fac0 * normal[i] / (length * length * length);

    for (int i = 0; i < 3; ++i) fac1 += naux_1[i] * normal[i];

    for (int i = 0; i < 3; ++i)
      n_1[i] = naux_1[i] / length - fac1 * normal[i] / (length * length * length);

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
    if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  }  // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL) dserror("Projector not converged!");

  // Newton iteration converged
  xi[0] = eta[0];
  xi[1] = eta[1];
  mele.ComputeUnitNormalAtXi(eta, normal);
  dist = alpha;

  // bye bye
  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNormal3DLin(MORTAR::Node& snode,
    MORTAR::Element& mele, double* xi, double* normal, double& dist,
    std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin)
{
  if (ndim_ != 3) dserror("ProjectSNodeByMNormal3D is only for 3D problems!");

  // start in the element center
  double eta[2] = {0.0, 0.0};
  if (distype == CORE::FE::CellType::tri3 || distype == CORE::FE::CellType::tri6)
  {
    eta[0] = 1.0 / 3.0;
    eta[1] = 1.0 / 3.0;
  }

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  std::array<double, 3> f = {0.0, 0.0, 0.0};

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  CORE::LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};
    double unormal[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 0.0};
    double normalpart[3] = {0.0, 0.0, 0.0};

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta, unormal);
    for (int i = 0; i < 3; ++i) normalpart[i] = unormal[i] * alpha;

    for (int i = 0; i < 3; ++i) normal[i] = unormal[i] * length;

    // calc xslave
    for (int i = 0; i < 3; ++i) xs[i] = snode.xspatial()[i];

    // calculate F
    for (int i = 0; i < 3; ++i) f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL) break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************
    // master coordinate grad
    std::array<double, 3> meta0 = {0.0, 0.0, 0.0};  // x , xi_0
    std::array<double, 3> meta1 = {0.0, 0.0, 0.0};  // x , xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta1, 2);

    // normal grad
    CORE::LINALG::Matrix<3, n_> secderiv;
    std::array<double, 3> meta00 = {0.0, 0.0, 0.0};  // x , xi_0 xi_0
    std::array<double, 3> meta11 = {0.0, 0.0, 0.0};  // x , xi_1 xi_1
    std::array<double, 3> meta01 = {0.0, 0.0, 0.0};  // x , xi_0 xi_1

    CORE::FE::shape_function_2D_deriv2(secderiv, eta[0], eta[1], distype);

    for (int i = 0; i < n_; ++i)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[i]);
      if (!mymnode) dserror("Null pointer!");
      for (int d = 0; d < 3; ++d)
      {
        meta00[d] += secderiv(0, i) * mymnode->xspatial()[d];
        meta11[d] += secderiv(1, i) * mymnode->xspatial()[d];
        meta01[d] += secderiv(2, i) * mymnode->xspatial()[d];
      }
    }

    std::array<double, 3> naux_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> naux_1 = {0.0, 0.0, 0.0};

    // normal grad xi_0
    naux_0[0] = (meta00[1] * meta1[2] - meta00[2] * meta1[1]);
    naux_0[1] = (meta00[2] * meta1[0] - meta00[0] * meta1[2]);
    naux_0[2] = (meta00[0] * meta1[1] - meta00[1] * meta1[0]);

    naux_0[0] += (meta0[1] * meta01[2] - meta0[2] * meta01[1]);
    naux_0[1] += (meta0[2] * meta01[0] - meta0[0] * meta01[2]);
    naux_0[2] += (meta0[0] * meta01[1] - meta0[1] * meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1] * meta1[2] - meta01[2] * meta1[1]);
    naux_1[1] = (meta01[2] * meta1[0] - meta01[0] * meta1[2]);
    naux_1[2] = (meta01[0] * meta1[1] - meta01[1] * meta1[0]);

    naux_1[0] += (meta0[1] * meta11[2] - meta0[2] * meta11[1]);
    naux_1[1] += (meta0[2] * meta11[0] - meta0[0] * meta11[2]);
    naux_1[2] += (meta0[0] * meta11[1] - meta0[1] * meta11[0]);

    std::array<double, 3> n_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> n_1 = {0.0, 0.0, 0.0};

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i = 0; i < 3; ++i) fac0 += naux_0[i] * normal[i];

    for (int i = 0; i < 3; ++i)
      n_0[i] = naux_0[i] / length - fac0 * normal[i] / (length * length * length);

    for (int i = 0; i < 3; ++i) fac1 += naux_1[i] * normal[i];

    for (int i = 0; i < 3; ++i)
      n_1[i] = naux_1[i] / length - fac1 * normal[i] / (length * length * length);

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
    if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  }  // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL) dserror("Projector not converged!");

  // Newton iteration converged
  xi[0] = eta[0];
  xi[1] = eta[1];
  mele.ComputeUnitNormalAtXi(eta, normal);
  dist = alpha;

  // bye bye
  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 08/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNodalNormal3DLin(MORTAR::Node& snode,
    MORTAR::Element& mele, double* xi, double* normal, double& dist,
    std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin)
{
  if (ndim_ != 3) dserror("ProjectSNodeByMNormal3DLin is only for 3D problems!");

  // start in the element center
  std::array<double, 2> eta = {0.0, 0.0};

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  std::array<double, 3> f = {0.0, 0.0, 0.0};

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  CORE::LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************
    // get shape function value
    CORE::LINALG::Matrix<n_, 1> mval;
    CORE::LINALG::Matrix<2, n_> mderiv;

    CORE::FE::shape_function_2D(mval, eta[0], eta[1], distype);
    CORE::FE::shape_function_2D_deriv1(mderiv, eta[0], eta[1], distype);

    // build interpolation of master node coordinates for current eta
    std::array<double, 3> xm = {0.0, 0.0, 0.0};
    std::array<double, 3> xs = {0.0, 0.0, 0.0};
    std::array<double, 3> normalnewton = {0.0, 0.0, 0.0};
    std::array<double, 3> normalpart = {0.0, 0.0, 0.0};

    // calc xmaster
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");
      for (int i = 0; i < 3; ++i)
      {
        xm[i] += mval(j) * mymnode->xspatial()[i];
      }
    }

    // calc xslave
    for (int i = 0; i < 3; ++i) xs[i] = snode.xspatial()[i];

    // calc normal part
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");
      for (int i = 0; i < 3; ++i)
      {
        normalnewton[i] += mval(j) * mymnode->MoData().n()[i];
      }
    }

    for (int i = 0; i < 3; ++i) normalpart[i] = normalnewton[i] * alpha;

    // calculate F
    for (int i = 0; i < 3; ++i) f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    //    std::cout << "k = " << k << "  tol= " << conv << std::endl;
    if (conv <= MORTARCONVTOL) break;
    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    std::array<double, 3> meta0 = {0.0, 0.0, 0.0};  // x , xi_0
    std::array<double, 3> meta1 = {0.0, 0.0, 0.0};  // x , xi_1

    // calc xmaster
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");
      for (int i = 0; i < 3; ++i)
      {
        meta0[i] += mderiv(0, j) * mymnode->xspatial()[i];
        meta1[i] += mderiv(1, j) * mymnode->xspatial()[i];
      }
    }

    std::array<double, 3> n_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> n_1 = {0.0, 0.0, 0.0};

    // calc normal part
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");

      for (int i = 0; i < 3; ++i)
      {
        n_0[i] += mderiv(0, j) * mymnode->MoData().n()[i];
        n_1[i] += mderiv(1, j) * mymnode->MoData().n()[i];
      }
    }

    // evaluate function f gradient
    for (int i = 0; i < 3; ++i)
    {
      df(i, 0) = meta0[i] + alpha * n_0[i];
      df(i, 1) = meta1[i] + alpha * n_1[i];
      df(i, 2) = normalnewton[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  }  // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
  {
    return false;
    dserror("Projector not converged!");
  }

  //**********************************************
  //   Get stuff                                //
  //**********************************************
  // Newton iteration converged
  xi[0] = eta[0];
  xi[1] = eta[1];
  dist = alpha;

  // get shape function value
  CORE::LINALG::Matrix<n_, 1> mval;
  CORE::LINALG::Matrix<2, n_> mderiv;

  CORE::FE::shape_function_2D(mval, eta[0], eta[1], distype);
  CORE::FE::shape_function_2D_deriv1(mderiv, eta[0], eta[1], distype);

  // calc normal part
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;
  for (int j = 0; j < n_; ++j)
  {
    Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
    if (!mymnode) dserror("Null pointer!");

    for (int i = 0; i < 3; ++i)
    {
      normal[i] += mval(j) * mymnode->MoData().n()[i];
    }
  }

  //**********************************************
  //   Lin deta = - inv(dF) * Lin F             //
  //**********************************************
  // prepare linearizations
  typedef CORE::GEN::Pairedvector<int, double>::const_iterator _CI;

  // get linsize
  int linsize = 0;
  for (int i = 0; i < n_; ++i)
  {
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);
    linsize += mnode->GetLinsize();
  }

  std::vector<CORE::GEN::Pairedvector<int, double>> xmLin(3, n_);  // nnode entry per dimension
  std::vector<CORE::GEN::Pairedvector<int, double>> xsLin(3, 1);   // one entry per dimension

  std::vector<CORE::GEN::Pairedvector<int, double>> normalpartLin(
      3, linsize);  // linsize of all mnodes
  std::vector<CORE::GEN::Pairedvector<int, double>> auxnormalLin(
      3, linsize);  // linsize of all mnodes

  std::vector<CORE::GEN::Pairedvector<int, double>> etaLin(3, linsize + n_ + 1);  // added all sizes
  std::vector<CORE::GEN::Pairedvector<int, double>> fLin(3, linsize + n_ + 1);    // added all sizes


  //--------------------------
  // master part:
  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    Node* mnode = dynamic_cast<Node*>(node);

    for (int k = 0; k < 3; ++k) (xmLin[k])[mnode->Dofs()[k]] += mval(i);
  }

  //--------------------------
  // normal part:
  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);

    for (int k = 0; k < 3; ++k)
    {
      for (_CI p = mnode->Data().GetDerivN()[k].begin(); p != mnode->Data().GetDerivN()[k].end();
           ++p)
      {
        (auxnormalLin[k])[p->first] += mval(i) * (p->second);
      }
    }
  }

  // *alpha
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = auxnormalLin[j].begin(); p != auxnormalLin[j].end(); ++p)
      (normalpartLin[j])[p->first] += (p->second) * alpha;
  }

  //--------------------------
  // slave part:
  for (int k = 0; k < 3; ++k) (xsLin[k])[snode.Dofs()[k]] += 1.0;

  // All terms:
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = xmLin[j].begin(); p != xmLin[j].end(); ++p) (fLin[j])[p->first] += (p->second);

    for (_CI p = normalpartLin[j].begin(); p != normalpartLin[j].end(); ++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p = xsLin[j].begin(); p != xsLin[j].end(); ++p) (fLin[j])[p->first] -= (p->second);
  }

  for (int i = 0; i < 3; ++i)
  {
    for (int k = 0; k < 3; ++k)
    {
      for (_CI p = fLin[i].begin(); p != fLin[i].end(); ++p)
        (etaLin[k])[p->first] -= (p->second) * df(k, i);
    }
  }

  //**********************************************
  //   Lin N                                    //
  //**********************************************
  std::vector<CORE::GEN::Pairedvector<int, double>> n_eta0_deriv(
      3, linsize + n_ + 1);  // added all sizes
  std::vector<CORE::GEN::Pairedvector<int, double>> n_eta1_deriv(
      3, linsize + n_ + 1);                                                 // added all sizes
  std::vector<CORE::GEN::Pairedvector<int, double>> n_n_deriv(3, linsize);  // linsize

  for (int k = 0; k < n_; ++k)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[k];
    if (!node) dserror("Cannot find master node");
    CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);

    for (_CI p = etaLin[0].begin(); p != etaLin[0].end(); ++p)
    {
      (n_eta0_deriv[0])[p->first] += mderiv(0, k) * (p->second) * mnode->MoData().n()[0];
      (n_eta0_deriv[1])[p->first] += mderiv(0, k) * (p->second) * mnode->MoData().n()[1];
      (n_eta0_deriv[2])[p->first] += mderiv(0, k) * (p->second) * mnode->MoData().n()[2];
    }

    for (_CI p = etaLin[1].begin(); p != etaLin[1].end(); ++p)
    {
      (n_eta1_deriv[0])[p->first] += mderiv(1, k) * (p->second) * mnode->MoData().n()[0];
      (n_eta1_deriv[1])[p->first] += mderiv(1, k) * (p->second) * mnode->MoData().n()[1];
      (n_eta1_deriv[2])[p->first] += mderiv(1, k) * (p->second) * mnode->MoData().n()[2];
    }

    for (_CI p = mnode->Data().GetDerivN()[0].begin(); p != mnode->Data().GetDerivN()[0].end(); ++p)
      (n_n_deriv[0])[p->first] += mval(k) * (p->second);
    for (_CI p = mnode->Data().GetDerivN()[1].begin(); p != mnode->Data().GetDerivN()[1].end(); ++p)
      (n_n_deriv[1])[p->first] += mval(k) * (p->second);
    for (_CI p = mnode->Data().GetDerivN()[2].begin(); p != mnode->Data().GetDerivN()[2].end(); ++p)
      (n_n_deriv[2])[p->first] += mval(k) * (p->second);
  }


  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = n_eta1_deriv[j].begin(); p != n_eta1_deriv[j].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second);
    for (_CI p = n_eta0_deriv[j].begin(); p != n_eta0_deriv[j].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second);
    for (_CI p = n_n_deriv[j].begin(); p != n_n_deriv[j].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second);
  }

  // bye bye
  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 05/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNodalNormal2DLin(MORTAR::Node& snode,
    MORTAR::Element& mele, double* xi, double* normal, double& dist,
    std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin)
{
  if (ndim_ != 2) dserror("ProjectSNodeByMNormal2DLin is only for 2D problems!");

  // start in the element center
  double eta[2] = {0.0, 0.0};

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  std::array<double, 2> f = {0.0, 0.0};

  // gradient of f (df/deta[0], df/dalpha)
  CORE::LINALG::Matrix<2, 2> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************
    // get shape function value
    CORE::LINALG::Matrix<n_, 1> mval;
    CORE::LINALG::Matrix<1, n_> mderiv;

    CORE::FE::shape_function_1D(mval, eta[0], distype);
    CORE::FE::shape_function_1D_deriv1(mderiv, eta[0], distype);

    // build interpolation of master node coordinates for current eta
    std::array<double, 3> xm = {0.0, 0.0, 0.0};
    std::array<double, 3> xs = {0.0, 0.0, 0.0};
    std::array<double, 3> normalnewton = {0.0, 0.0, 0.0};
    std::array<double, 3> normalpart = {0.0, 0.0, 0.0};

    // calc xmaster
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");
      for (int i = 0; i < 2; ++i)
      {
        xm[i] += mval(j) * mymnode->xspatial()[i];
      }
    }

    // calc xslave
    for (int i = 0; i < 2; ++i) xs[i] = snode.xspatial()[i];

    // calc normal part
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");
      for (int i = 0; i < 2; ++i)
      {
        normalnewton[i] += mval(j) * mymnode->MoData().n()[i];
      }
    }

    for (int i = 0; i < 3; ++i) normalpart[i] = normalnewton[i] * alpha;

    // calculate F
    for (int i = 0; i < 2; ++i) f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1]);
    //    std::cout << "k = " << k << "  tol= " << conv << std::endl;
    if (conv <= MORTARCONVTOL) break;
    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    std::array<double, 3> meta0 = {0.0, 0.0, 0.0};  // x , xi_0

    // calc xmaster
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");
      for (int i = 0; i < 2; ++i)
      {
        meta0[i] += mderiv(0, j) * mymnode->xspatial()[i];
      }
    }

    std::array<double, 3> n_0 = {0.0, 0.0, 0.0};

    // calc normal part
    for (int j = 0; j < n_; ++j)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
      if (!mymnode) dserror("Null pointer!");

      for (int i = 0; i < 2; ++i)
      {
        n_0[i] += mderiv(0, j) * mymnode->MoData().n()[i];
      }
    }

    // evaluate function f gradient
    for (int i = 0; i < 2; ++i)
    {
      df(i, 0) = meta0[i] + alpha * n_0[i];
      df(i, 1) = normalnewton[i];
    }

    //**********************************************
    //   solve deta = - inv(dF) * F               //
    //**********************************************
    double jacdet = df.Invert();
    if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1];
    alpha += -df(1, 0) * f[0] - df(1, 1) * f[1];
  }  // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL)
  {
    return false;
    dserror("Projector not converged!");
  }

  //**********************************************
  //   Get stuff                                //
  //**********************************************
  // Newton iteration converged
  xi[0] = eta[0];
  xi[1] = eta[1];
  dist = alpha;

  // get shape function value
  CORE::LINALG::Matrix<n_, 1> mval;
  CORE::LINALG::Matrix<1, n_> mderiv;

  CORE::FE::shape_function_1D(mval, eta[0], distype);
  CORE::FE::shape_function_1D_deriv1(mderiv, eta[0], distype);

  // calc normal part
  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;
  for (int j = 0; j < n_; ++j)
  {
    Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[j]);
    if (!mymnode) dserror("Null pointer!");

    for (int i = 0; i < 2; ++i)
    {
      normal[i] += mval(j) * mymnode->MoData().n()[i];
    }
  }

  //**********************************************
  //   Lin deta = - inv(dF) * Lin F             //
  //**********************************************
  // prepare linearizations
  typedef CORE::GEN::Pairedvector<int, double>::const_iterator _CI;

  std::vector<CORE::GEN::Pairedvector<int, double>> etaLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> fLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> xmLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> normalpartLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> xsLin(3, 1000);

  //--------------------------
  // master part:
  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    Node* mnode = dynamic_cast<Node*>(node);

    for (int k = 0; k < 2; ++k) (xmLin[k])[mnode->Dofs()[k]] += mval(i);
  }

  //--------------------------
  // normal part:
  std::vector<CORE::GEN::Pairedvector<int, double>> x_0Lin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> auxnormalLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> auxnormalunitLin(3, 1000);

  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);

    for (int k = 0; k < 3; ++k)
    {
      for (_CI p = mnode->Data().GetDerivN()[k].begin(); p != mnode->Data().GetDerivN()[k].end();
           ++p)
      {
        (auxnormalLin[k])[p->first] += mval(i) * (p->second);
      }
    }
  }

  // *alpha
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = auxnormalLin[j].begin(); p != auxnormalLin[j].end(); ++p)
      (normalpartLin[j])[p->first] += (p->second) * alpha;
  }

  //--------------------------
  // slave part:
  for (int k = 0; k < 2; ++k) (xsLin[k])[snode.Dofs()[k]] += 1.0;

  // All terms:
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = xmLin[j].begin(); p != xmLin[j].end(); ++p) (fLin[j])[p->first] += (p->second);

    for (_CI p = normalpartLin[j].begin(); p != normalpartLin[j].end(); ++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p = xsLin[j].begin(); p != xsLin[j].end(); ++p) (fLin[j])[p->first] -= (p->second);
  }

  for (int i = 0; i < 2; ++i)
  {
    for (int k = 0; k < 2; ++k)
    {
      for (_CI p = fLin[i].begin(); p != fLin[i].end(); ++p)
        (etaLin[k])[p->first] -= (p->second) * df(k, i);
    }
  }

  //**********************************************
  //   Lin N                                    //
  //**********************************************
  std::vector<CORE::GEN::Pairedvector<int, double>> n_eta0_deriv(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> n_eta1_deriv(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> n_n_deriv(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> normaltolineLinaux(3, 1000);

  for (int k = 0; k < n_; ++k)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[k];
    if (!node) dserror("Cannot find master node");
    CONTACT::Node* mnode = dynamic_cast<CONTACT::Node*>(node);

    for (_CI p = etaLin[0].begin(); p != etaLin[0].end(); ++p)
    {
      (n_eta0_deriv[0])[p->first] += mderiv(0, k) * (p->second) * mnode->MoData().n()[0];
      (n_eta0_deriv[1])[p->first] += mderiv(0, k) * (p->second) * mnode->MoData().n()[1];
      (n_eta0_deriv[2])[p->first] += mderiv(0, k) * (p->second) * mnode->MoData().n()[2];
    }

    for (_CI p = mnode->Data().GetDerivN()[0].begin(); p != mnode->Data().GetDerivN()[0].end(); ++p)
      (n_n_deriv[0])[p->first] += mval(k) * (p->second);
    for (_CI p = mnode->Data().GetDerivN()[1].begin(); p != mnode->Data().GetDerivN()[1].end(); ++p)
      (n_n_deriv[1])[p->first] += mval(k) * (p->second);
    for (_CI p = mnode->Data().GetDerivN()[2].begin(); p != mnode->Data().GetDerivN()[2].end(); ++p)
      (n_n_deriv[2])[p->first] += mval(k) * (p->second);
  }


  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = n_eta0_deriv[j].begin(); p != n_eta0_deriv[j].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second);
    for (_CI p = n_n_deriv[j].begin(); p != n_n_deriv[j].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second);
  }

  // bye bye
  return true;
}



/*----------------------------------------------------------------------*
 |  Project snode onto melement with master element normal   farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNormal2D(
    MORTAR::Node& snode, MORTAR::Element& mele, double* xi, double* normal, double& dist)
{
  if (ndim_ != 2) dserror("ProjectSNodeByMNormal2D is only for 2D problems!");

  // start in the element center
  double eta[2] = {0.0, 0.0};

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  std::array<double, 3> f = {0.0, 0.0, 0.0};

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  CORE::LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};
    double unormal[3] = {0.0, 0.0, 0.0};
    double normal[3] = {0.0, 0.0, 0.0};
    double normalpart[3] = {0.0, 0.0, 0.0};

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta, unormal);
    for (int i = 0; i < 3; ++i) normalpart[i] = unormal[i] * alpha;

    for (int i = 0; i < 3; ++i) normal[i] = unormal[i] * length;

    // calc xslave
    for (int i = 0; i < 3; ++i) xs[i] = snode.xspatial()[i];

    // calculate F
    for (int i = 0; i < 3; ++i) f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL) break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    double meta0[3] = {0.0, 0.0, 0.0};  // x,xi_0
    double meta1[3] = {0.0, 0.0, 1.0};  // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);

    // normal grad
    CORE::LINALG::Matrix<1, n_> secderiv;
    std::array<double, 3> meta00 = {0.0, 0.0, 0.0};  // x , xi_0 xi_0
    std::array<double, 3> meta11 = {0.0, 0.0, 0.0};  // x , xi_1 xi_1
    std::array<double, 3> meta01 = {0.0, 0.0, 0.0};  // x , xi_0 xi_1

    CORE::FE::shape_function_1D_deriv2(secderiv, eta[0], distype);

    for (int i = 0; i < n_; ++i)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[i]);
      if (!mymnode) dserror("Null pointer!");
      for (int d = 0; d < 3; ++d) meta00[d] += secderiv(0, i) * mymnode->xspatial()[d];
    }

    std::array<double, 3> naux_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> naux_1 = {0.0, 0.0, 0.0};

    // normal grad xi_0
    naux_0[0] = (meta00[1] * meta1[2] - meta00[2] * meta1[1]);
    naux_0[1] = (meta00[2] * meta1[0] - meta00[0] * meta1[2]);
    naux_0[2] = (meta00[0] * meta1[1] - meta00[1] * meta1[0]);

    naux_0[0] += (meta0[1] * meta01[2] - meta0[2] * meta01[1]);
    naux_0[1] += (meta0[2] * meta01[0] - meta0[0] * meta01[2]);
    naux_0[2] += (meta0[0] * meta01[1] - meta0[1] * meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1] * meta1[2] - meta01[2] * meta1[1]);
    naux_1[1] = (meta01[2] * meta1[0] - meta01[0] * meta1[2]);
    naux_1[2] = (meta01[0] * meta1[1] - meta01[1] * meta1[0]);

    naux_1[0] += (meta0[1] * meta11[2] - meta0[2] * meta11[1]);
    naux_1[1] += (meta0[2] * meta11[0] - meta0[0] * meta11[2]);
    naux_1[2] += (meta0[0] * meta11[1] - meta0[1] * meta11[0]);

    std::array<double, 3> n_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> n_1 = {0.0, 0.0, 0.0};

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i = 0; i < 3; ++i) fac0 += naux_0[i] * normal[i];

    for (int i = 0; i < 3; ++i)
      n_0[i] = naux_0[i] / length - fac0 * normal[i] / (length * length * length);

    for (int i = 0; i < 3; ++i) fac1 += naux_1[i] * normal[i];

    for (int i = 0; i < 3; ++i)
      n_1[i] = naux_1[i] / length - fac1 * normal[i] / (length * length * length);

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
    if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  }  // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL) dserror("Projector not converged!");

  // Newton iteration converged
  xi[0] = eta[0];
  mele.ComputeUnitNormalAtXi(eta, normal);
  dist = alpha;

  // bye bye
  return true;
}

/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal + Lin     farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNormal2DLin(MORTAR::Node& snode,
    MORTAR::Element& mele, double* xi, double* normal, double& dist,
    std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin)
{
  if (ndim_ != 2) dserror("ProjectSNodeByMNormal2D is only for 2D problems!");

  // start in the element center
  double eta[2] = {0.0, 0.0};

  // auxiliary variable for distance
  double alpha = 0.0;

  // function f (vector-valued)
  std::array<double, 3> f = {0.0, 0.0, 0.0};

  // gradient of f (df/deta[0], df/deta[1], df/dalpha)
  CORE::LINALG::Matrix<3, 3> df;

  // start iteration
  int k = 0;
  double conv = 0.0;

  for (k = 0; k < MORTARMAXITER; ++k)
  {
    //**********************************************
    //  F CALCULATION                             //
    //**********************************************

    // build interpolation of master node coordinates for current eta
    double xm[3] = {0.0, 0.0, 0.0};
    double xs[3] = {0.0, 0.0, 0.0};
    double unormal[3] = {0.0, 0.0, 0.0};
    double normal_k[3] = {0.0, 0.0, 0.0};
    double normalpart[3] = {0.0, 0.0, 0.0};

    // calc xmaster
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, xm, 0);

    // calc normal part
    double length = mele.ComputeUnitNormalAtXi(eta, unormal);
    for (int i = 0; i < 3; ++i) normalpart[i] = unormal[i] * alpha;

    for (int i = 0; i < 3; ++i) normal_k[i] = unormal[i] * length;

    // calc xslave
    for (int i = 0; i < 3; ++i) xs[i] = snode.xspatial()[i];

    // calculate F
    for (int i = 0; i < 3; ++i) f[i] = xm[i] + normalpart[i] - xs[i];

    // check for convergence
    conv = sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
    if (conv <= MORTARCONVTOL) break;

    //**********************************************
    //   F GRADIENT CALCULATION                   //
    //**********************************************

    // master coordinate grad
    double meta0[3] = {0.0, 0.0, 0.0};  // x,xi_0
    double meta1[3] = {0.0, 0.0, 1.0};  // x,xi_1
    MORTAR::UTILS::LocalToGlobal<distype>(mele, eta, meta0, 1);

    // normal grad
    CORE::LINALG::Matrix<1, n_> secderiv;
    std::array<double, 3> meta00 = {0.0, 0.0, 0.0};  // x , xi_0 xi_0
    std::array<double, 3> meta11 = {0.0, 0.0, 0.0};  // x , xi_1 xi_1
    std::array<double, 3> meta01 = {0.0, 0.0, 0.0};  // x , xi_0 xi_1

    CORE::FE::shape_function_1D_deriv2(secderiv, eta[0], distype);

    for (int i = 0; i < n_; ++i)
    {
      Node* mymnode = dynamic_cast<Node*>(mele.Nodes()[i]);
      if (!mymnode) dserror("Null pointer!");
      for (int d = 0; d < 3; ++d) meta00[d] += secderiv(0, i) * mymnode->xspatial()[d];
    }

    std::array<double, 3> naux_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> naux_1 = {0.0, 0.0, 0.0};

    // normal grad xi_0
    naux_0[0] = (meta00[1] * meta1[2] - meta00[2] * meta1[1]);
    naux_0[1] = (meta00[2] * meta1[0] - meta00[0] * meta1[2]);
    naux_0[2] = (meta00[0] * meta1[1] - meta00[1] * meta1[0]);

    naux_0[0] += (meta0[1] * meta01[2] - meta0[2] * meta01[1]);
    naux_0[1] += (meta0[2] * meta01[0] - meta0[0] * meta01[2]);
    naux_0[2] += (meta0[0] * meta01[1] - meta0[1] * meta01[0]);

    // normal grad xi_1
    naux_1[0] = (meta01[1] * meta1[2] - meta01[2] * meta1[1]);
    naux_1[1] = (meta01[2] * meta1[0] - meta01[0] * meta1[2]);
    naux_1[2] = (meta01[0] * meta1[1] - meta01[1] * meta1[0]);

    naux_1[0] += (meta0[1] * meta11[2] - meta0[2] * meta11[1]);
    naux_1[1] += (meta0[2] * meta11[0] - meta0[0] * meta11[2]);
    naux_1[2] += (meta0[0] * meta11[1] - meta0[1] * meta11[0]);

    std::array<double, 3> n_0 = {0.0, 0.0, 0.0};
    std::array<double, 3> n_1 = {0.0, 0.0, 0.0};

    double fac0 = 0.0;
    double fac1 = 0.0;

    for (int i = 0; i < 3; ++i) fac0 += naux_0[i] * normal_k[i];

    for (int i = 0; i < 3; ++i)
      n_0[i] = naux_0[i] / length - fac0 * normal_k[i] / (length * length * length);

    for (int i = 0; i < 3; ++i) fac1 += naux_1[i] * normal_k[i];

    for (int i = 0; i < 3; ++i)
      n_1[i] = naux_1[i] / length - fac1 * normal_k[i] / (length * length * length);

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
    if (abs(jacdet) < 1.0e-12) dserror("Singular Jacobian for projection");

    // update eta and alpha
    eta[0] += -df(0, 0) * f[0] - df(0, 1) * f[1] - df(0, 2) * f[2];
    eta[1] += -df(1, 0) * f[0] - df(1, 1) * f[1] - df(1, 2) * f[2];
    alpha += -df(2, 0) * f[0] - df(2, 1) * f[1] - df(2, 2) * f[2];
  }  // end newton loop

  // Newton iteration unconverged
  if (conv > MORTARCONVTOL) dserror("Projector not converged!");

  //**********************************************
  //   Get stuff                                //
  //**********************************************
  // Newton iteration converged
  xi[0] = eta[0];
  double normallength = mele.ComputeUnitNormalAtXi(eta, normal);
  dist = alpha;


  //**********************************************
  //   Lin deta = - inv(dF) * Lin F             //
  //**********************************************
  // prepare linearizations
  typedef CORE::GEN::Pairedvector<int, double>::const_iterator _CI;

  std::vector<CORE::GEN::Pairedvector<int, double>> etaLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> fLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> xmLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> normalpartLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> xsLin(3, 1000);

  //--------------------------
  // master part:
  CORE::LINALG::Matrix<n_, 1> val;
  CORE::FE::shape_function_1D(val, eta[0], distype);

  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    Node* mnode = dynamic_cast<Node*>(node);

    for (int k = 0; k < 2; ++k) (xmLin[k])[mnode->Dofs()[k]] += val(i);
  }

  //--------------------------
  // normal part:
  std::vector<CORE::GEN::Pairedvector<int, double>> x_0Lin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> auxnormalLin(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> auxnormalunitLin(3, 1000);

  CORE::LINALG::Matrix<1, n_> deriv1;
  CORE::FE::shape_function_1D_deriv1(deriv1, eta[0], distype);
  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    Node* mnode = dynamic_cast<Node*>(node);

    for (int k = 0; k < 2; ++k) (x_0Lin[k])[mnode->Dofs()[k]] += deriv1(i);
  }

  // cross product linearization
  for (_CI p = x_0Lin[1].begin(); p != x_0Lin[1].end(); ++p)
    (auxnormalLin[0])[p->first] += (p->second);
  for (_CI p = x_0Lin[0].begin(); p != x_0Lin[0].end(); ++p)
    (auxnormalLin[1])[p->first] -= (p->second);

  // calc normalpart without alpha
  std::array<double, 3> linnormalaux = {0.0, 0.0, 0.0};
  linnormalaux[0] = normal[0] * normallength;
  linnormalaux[1] = normal[1] * normallength;
  linnormalaux[2] = normal[2] * normallength;

  // derivative weighting matrix for current element
  CORE::LINALG::Matrix<3, 3> W;
  const double lcubeinv = 1.0 / (normallength * normallength * normallength);

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      W(j, k) = -lcubeinv * linnormalaux[j] * linnormalaux[k];
      if (j == k) W(j, k) += 1 / normallength;
    }
  }

  // row loop
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = auxnormalLin[0].begin(); p != auxnormalLin[0].end(); ++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j, 0);

    for (_CI p = auxnormalLin[1].begin(); p != auxnormalLin[1].end(); ++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j, 1);

    for (_CI p = auxnormalLin[2].begin(); p != auxnormalLin[2].end(); ++p)
      (auxnormalunitLin[j])[p->first] += (p->second) * W(j, 2);
  }

  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = auxnormalunitLin[j].begin(); p != auxnormalunitLin[j].end(); ++p)
      (normalpartLin[j])[p->first] += (p->second) * alpha;
  }


  //--------------------------
  // slave part:
  for (int k = 0; k < 2; ++k) (xsLin[k])[snode.Dofs()[k]] += 1.0;

  // All terms:
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = xmLin[j].begin(); p != xmLin[j].end(); ++p) (fLin[j])[p->first] += (p->second);

    for (_CI p = normalpartLin[j].begin(); p != normalpartLin[j].end(); ++p)
      (fLin[j])[p->first] += (p->second);

    for (_CI p = xsLin[j].begin(); p != xsLin[j].end(); ++p) (fLin[j])[p->first] -= (p->second);
  }


  for (int i = 0; i < 3; ++i)
  {
    for (int k = 0; k < 3; ++k)
    {
      for (_CI p = fLin[i].begin(); p != fLin[i].end(); ++p)
        (etaLin[k])[p->first] -= (p->second) * df(k, i);
    }
  }

  //**********************************************
  //   Lin N                                    //
  //**********************************************
  std::vector<CORE::GEN::Pairedvector<int, double>> x_0Linnew(3, 1000);
  std::vector<CORE::GEN::Pairedvector<int, double>> normaltolineLinaux(3, 1000);

  CORE::LINALG::Matrix<1, n_> deriv;
  CORE::FE::shape_function_1D_deriv1(deriv, eta[0], distype);

  CORE::LINALG::Matrix<1, n_> deriv2;
  CORE::FE::shape_function_1D_deriv2(deriv2, eta[0], distype);
  for (int i = 0; i < n_; ++i)
  {
    // get master node
    DRT::Node* node = mele.Nodes()[i];
    if (!node) dserror("Cannot find master node");
    Node* mnode = dynamic_cast<Node*>(node);

    for (int k = 0; k < 2; ++k) (x_0Linnew[k])[mnode->Dofs()[k]] += deriv(i);

    for (int k = 0; k < 2; ++k)
    {
      for (_CI p = etaLin[0].begin(); p != etaLin[0].end(); ++p)
        (x_0Linnew[k])[p->first] += (p->second) * deriv2(i) * mnode->xspatial()[k];
    }
  }

  // cross product linearization
  for (_CI p = x_0Linnew[1].begin(); p != x_0Linnew[1].end(); ++p)
    (normaltolineLinaux[0])[p->first] += (p->second);
  for (_CI p = x_0Linnew[0].begin(); p != x_0Linnew[0].end(); ++p)
    (normaltolineLinaux[1])[p->first] -= (p->second);

  // normalize lin
  CORE::LINALG::Matrix<3, 3> Wfinal;
  //  const double lcubeinv = 1.0 / (normallength * normallength * normallength);

  for (int j = 0; j < 3; ++j)
  {
    for (int k = 0; k < 3; ++k)
    {
      Wfinal(j, k) = -lcubeinv * linnormalaux[j] * linnormalaux[k];
      if (j == k) Wfinal(j, k) += 1 / normallength;
    }
  }

  // row loop
  for (int j = 0; j < 3; ++j)
  {
    for (_CI p = normaltolineLinaux[0].begin(); p != normaltolineLinaux[0].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j, 0);

    for (_CI p = normaltolineLinaux[1].begin(); p != normaltolineLinaux[1].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j, 1);

    for (_CI p = normaltolineLinaux[2].begin(); p != normaltolineLinaux[2].end(); ++p)
      (normaltolineLin[j])[p->first] += (p->second) * Wfinal(j, 2);
  }

  // bye bye
  return true;
}

/*----------------------------------------------------------------------*
 |  Project snode onto melement with master element normal   farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNormal(
    MORTAR::Node& snode, MORTAR::Element& mele, double* xi, double* normal, double& dist)
{
  if (ndim_ == 2)
  {
    ProjectSNodeByMNormal2D(snode, mele, xi, normal, dist);
  }
  else if (ndim_ == 3)
  {
    ProjectSNodeByMNormal3D(snode, mele, xi, normal, dist);
  }
  else
  {
    dserror("wrong dimension!");
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Project snode onto melement with master element normal   farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNodalNormalLin(MORTAR::Node& snode,
    MORTAR::Element& mele, double* xi, double* normal, double& dist,
    std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin)
{
  bool success = false;

  if (ndim_ == 2)
  {
    success = ProjectSNodeByMNodalNormal2DLin(snode, mele, xi, normal, dist, normaltolineLin);
  }
  else if (ndim_ == 3)
  {
    success = ProjectSNodeByMNodalNormal3DLin(snode, mele, xi, normal, dist, normaltolineLin);
  }
  else
  {
    dserror("wrong dimension!");
  }

  return success;
}

/*----------------------------------------------------------------------*
 |  Project snode onto melement with master normal           farah 01/16|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::ProjectSNodeByMNormalLin(MORTAR::Node& snode,
    MORTAR::Element& mele, double* xi, double* normal, double& dist,
    std::vector<CORE::GEN::Pairedvector<int, double>>& normaltolineLin)
{
  if (ndim_ == 2)
  {
    ProjectSNodeByMNormal2DLin(snode, mele, xi, normal, dist, normaltolineLin);
  }
  else if (ndim_ == 3)
  {
    dserror("Not yet implemented!");
  }
  else
  {
    dserror("wrong dimension!");
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for nodal normal case (public)                 popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double MORTAR::ProjectorCalc<distype>::EvaluateFNodalNormal(
    MORTAR::Node& node, MORTAR::Element& ele, const double* eta)
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
  for (int i = 0; i < ndim_; ++i) nx[i] -= node.xspatial()[i];

  // calculate F
  fval = nx[0] * node.MoData().n()[1] - nx[1] * node.MoData().n()[0];

  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for nodal normal case (public)             popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double MORTAR::ProjectorCalc<distype>::EvaluateGradFNodalNormal(
    MORTAR::Node& node, MORTAR::Element& ele, const double* eta)
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
template <CORE::FE::CellType distype>
double MORTAR::ProjectorCalc<distype>::EvaluateFElementNormal(
    MORTAR::Node& node, MORTAR::Element& ele, const double* eta)
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
  if (!mynodes) dserror("EvaluateFElementNormal: Null pointer!");

  CORE::LINALG::Matrix<n_, 1> val;
  CORE::LINALG::Matrix<ndim_, n_> coord;

  // get shape function values and derivatives at gpeta
  if (distype == CORE::FE::CellType::nurbs2 || distype == CORE::FE::CellType::nurbs3)
  {
    CORE::LINALG::SerialDenseVector auxval(n_);
    CORE::LINALG::SerialDenseMatrix deriv(n_, 1);
    ele.EvaluateShape(eta, auxval, deriv, ele.NumNode());

    for (int i = 0; i < n_; ++i) val(i) = auxval(i);
  }
  else
    CORE::FE::shape_function_1D(val, eta[0], distype);

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
    Node* mymrtrnode = dynamic_cast<Node*>(mynodes[i]);

    for (int j = 0; j < ndim_; ++j)
    {
      nn[j] += val(i) * mymrtrnode->MoData().n()[j];

      coord(j, i) = mymrtrnode->xspatial()[j];

      nx[j] += val(i) * coord(j, i);
    }
  }

  // subtract master node coordinates
  for (int j = 0; j < ndim_; ++j) nx[j] -= node.xspatial()[j];

  // calculate F
  fval = nx[0] * nn[1] - nx[1] * nn[0];

  return fval;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for element normal case (public)           popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
double MORTAR::ProjectorCalc<distype>::EvaluateGradFElementNormal(
    MORTAR::Node& node, MORTAR::Element& ele, const double* eta)
{
  if (ndim_ == 3) dserror("This Projector is only for 2D Problems!");

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
  CORE::LINALG::Matrix<n_, 1> val;
  CORE::LINALG::Matrix<ndim_ - 1, n_> deriv;
  CORE::LINALG::Matrix<ndim_, n_> coord;

  DRT::Node** mynodes = ele.Nodes();
  if (!mynodes) dserror("EvaluateGradFElementNormal: Null pointer!");

  // get shape function values and derivatives at gpeta
  if (distype == CORE::FE::CellType::nurbs2 || distype == CORE::FE::CellType::nurbs3)
  {
    CORE::LINALG::SerialDenseVector auxval(n_);
    CORE::LINALG::SerialDenseMatrix auxderiv(n_, 1);
    ele.EvaluateShape(eta, auxval, auxderiv, ele.NumNode());

    for (int i = 0; i < n_; ++i)
    {
      val(i) = auxval(i);
      deriv(0, i) = auxderiv(i, 0);
    }
  }
  else
  {
    CORE::FE::shape_function_1D(val, eta[0], distype);
    CORE::FE::shape_function_1D_deriv1(deriv, eta[0], distype);
  }

  // get interpolated normal and proj. coordinates for current eta
  double nn[ndim_];
  double nneta[ndim_];
  double nx[ndim_];
  double nxeta[ndim_];
  for (int j = 0; j < ndim_; ++j)
  {
    nn[j] = 0.0;
    nneta[j] = 0.0;
    nx[j] = 0.0;
    nxeta[j] = 0.0;
  }

  for (int i = 0; i < n_; ++i)
  {
    Node* mymrtrnode = dynamic_cast<Node*>(mynodes[i]);

    for (int j = 0; j < ndim_; ++j)
    {
      nn[j] += val(i) * mymrtrnode->MoData().n()[j];
      nneta[j] += deriv(0, i) * mymrtrnode->MoData().n()[j];
      coord(j, i) = mymrtrnode->xspatial()[j];
      nx[j] += val(i) * coord(j, i);
      nxeta[j] += deriv(0, i) * coord(j, i);
    }
  }

  // subtract master node coordinates
  for (int j = 0; j < ndim_; ++j) nx[j] -= node.xspatial()[j];

  // calculate GradF
  fgrad = nxeta[0] * nn[1] + nx[0] * nneta[1] - nxeta[1] * nn[0] - nx[1] * nneta[0];

  return fgrad;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for Gauss point case (public)                  popp 01/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
double MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::EvaluateFGaussPoint2D(
    const double* gpx, const double* gpn, MORTAR::Element& ele, const double* eta)
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
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
double MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::EvaluateGradFGaussPoint2D(
    const double* gpn, MORTAR::Element& ele, const double* eta)
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
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
bool MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::EvaluateFGaussPoint3D(double* f,
    const double* gpx, const double* gpn, MORTAR::Element& ele, const double* eta,
    const double& alpha)
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
  for (int i = 0; i < ndim_; ++i) f[i] = nx[i] - alpha * gpn[i] - gpx[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for Gauss point case (3D)                  popp 11/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distypeS, CORE::FE::CellType distypeM>
bool MORTAR::ProjectorCalcEleBased<distypeS, distypeM>::EvaluateGradFGaussPoint3D(
    CORE::LINALG::Matrix<3, 3>& fgrad, const double* gpx, const double* gpn, MORTAR::Element& ele,
    const double* eta, const double& alpha)
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
  double* nxeta1 = &fgrad(0, 0);
  double* nxeta2 = &fgrad(0, 1);
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nxeta1, 1);
  MORTAR::UTILS::LocalToGlobal<distypeM>(ele, eta, nxeta2, 2);

  // evaluate function f gradient
  for (int i = 0; i < ndim_; ++i) fgrad(i, 2) = -gpn[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate F for AuxPlane Gauss point case (3D)             popp 11/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::EvaluateFGaussPointAuxn3D(double* f, const double* globgp,
    const double* auxn, MORTAR::Element& ele, const double* eta, const double& alpha)
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
  for (int i = 0; i < ndim_; ++i) f[i] = nx[i] - alpha * auxn[i] - globgp[i];

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate GradF for AuxPlane Gauss point case (3D)         popp 11/08|
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype>
bool MORTAR::ProjectorCalc<distype>::EvaluateGradFGaussPointAuxn3D(
    CORE::LINALG::Matrix<3, 3>& fgrad, const double* globgp, const double* auxn,
    MORTAR::Element& ele, const double* eta, const double& alpha)
{
  /* Evaluate the gradient of the function F(eta,alpha) = Ni * xi -
   - alpha * auxn - globgp, which is a (3x3)-matrix!

   Ni      shape functions of element to project on
   xi      coords of nodes of element to project on
   globgp  coords of AuxPlaneGP to be projected
   auxn    normal of AuxPlane along which to project            */

  // build interpolation of ele node coordinates for current eta
  //  double nxeta1[ndim_];
  //  double nxeta2[ndim_];
  double* nxeta1 = &fgrad(0, 0);
  double* nxeta2 = &fgrad(0, 1);
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nxeta1, 1);
  MORTAR::UTILS::LocalToGlobal<distype>(ele, eta, nxeta2, 2);

  // evaluate function f gradient
  for (int i = 0; i < ndim_; ++i) fgrad(i, 2) = -auxn[i];

  return true;
}

FOUR_C_NAMESPACE_CLOSE
