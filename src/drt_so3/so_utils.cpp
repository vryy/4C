/*----------------------------------------------------------------------*/
/*! \file

\brief A collection of helper methods for solid elements

\level 1
\maintainer Amadeus Gebauer
 *-----------------------------------------------------------------------*/

#include "so_utils.H"

#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fiber/drt_fiber_node.H"

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::CalcR(const DRT::Element* ele, const std::vector<double>& disp,
    LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, DRT::UTILS::DisTypeToDim<distype>::dim>&
        R)
{
  // number of nodes per element
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // spatial dimension
  const int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;

  if (disp.size() != nsd * nen) dserror("mismatch in dimensions");

  LINALG::Matrix<nsd, 1> xi_ele_center =
      DRT::UTILS::getLocalCenterPosition<nsd>(distype);  // depending on distype

  LINALG::Matrix<nen, nsd> xrefe;  // X, material coord. of element
  LINALG::Matrix<nen, nsd> xcurr;  // x, current  coord. of element
  for (int i = 0; i < nen; ++i)
  {
    for (int d = 0; d < nsd; ++d)
    {
      xrefe(i, d) = ele->Nodes()[i]->X()[d];
      xcurr(i, d) = ele->Nodes()[i]->X()[d] + disp[i * nsd + d];
    }
  }
  LINALG::Matrix<nsd, nen> deriv;
  DRT::UTILS::shape_function_deriv1<distype>(xi_ele_center, deriv);

  LINALG::Matrix<nsd, nsd> jac;
  LINALG::Matrix<nsd, nsd> defgrd;
  LINALG::Matrix<nsd, nen> deriv_xyz;
  jac.Multiply(deriv, xrefe);
  jac.Invert();
  deriv_xyz.Multiply(jac, deriv);
  defgrd.MultiplyTT(xcurr, deriv_xyz);

  // Calculate rotcurr from defgrd
  LINALG::Matrix<nsd, nsd> Q(true);
  LINALG::Matrix<nsd, nsd> S(true);
  LINALG::Matrix<nsd, nsd> VT(true);
  LINALG::SVD<nsd, nsd>(defgrd, Q, S, VT);
  R.MultiplyNN(Q, VT);
}


/*----------------------------------------------------------------------*
 | interpolate nodal fibers to gauss point                              |
 |                                pfaller (adapted from hoermann) 07/18 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::UTILS::NodalFiber(DRT::Node** nodes,
    const std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement,
        1>>& shapefcts,
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1>>& gpfiber1,
    std::vector<LINALG::Matrix<DRT::UTILS::DisTypeToDim<distype>::dim, 1>>& gpfiber2)
{
  // number of nodes per element
  const int nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

  // spatial dimension
  const int nsd = DRT::UTILS::DisTypeToDim<distype>::dim;

  // number of gauss points
  const int ngp = shapefcts.size();

  std::vector<Epetra_SerialDenseVector> fibernodes1(nen, Epetra_SerialDenseVector(3));
  std::vector<Epetra_SerialDenseVector> fibernodes2(nen, Epetra_SerialDenseVector(3));
  std::vector<Epetra_SerialDenseVector> cirnodes(nen, Epetra_SerialDenseVector(3));
  std::vector<Epetra_SerialDenseVector> tannodes(nen, Epetra_SerialDenseVector(3));
  std::vector<double> helixnodes(nen);
  std::vector<double> transversenodes(nen);

  // read cosy information from node
  for (int inode = 0; inode < nen; ++inode)
  {
    DRT::FIBER::FiberNode* fnode = dynamic_cast<DRT::FIBER::FiberNode*>(nodes[inode]);
    if (!fnode) dserror("No fiber direction defined on nodes or elements");
    for (int j = 0; j < nsd; ++j)
    {
      fibernodes1[inode](j) = fnode->Fiber1()[j];
      fibernodes2[inode](j) = fnode->Fiber2()[j];
      cirnodes[inode](j) = fnode->Cir()[j];
      tannodes[inode](j) = fnode->Tan()[j];
    }
    helixnodes[inode] = fnode->Helix();
    transversenodes[inode] = fnode->Transverse();
  }

  // case fibers 1 are interpolated
  if ((fibernodes1[0]).Norm2() > 1e-13)
  {
    // interpolate fibers to integration points
    for (int i = 0; i < nsd; ++i)
      for (int q = 0; q < ngp; ++q)
        for (int j = 0; j < nen; ++j) gpfiber1[q](i, 0) += shapefcts[q](j) * fibernodes1[j](i);
    for (int q = 0; q < ngp; ++q) gpfiber1[q].Scale(1.0 / gpfiber1[q].Norm2());

    // fiber 2
    if ((fibernodes2[0]).Norm2() > 1e-13)
    {
      // interpolate fibers to integration points
      for (int i = 0; i < nsd; ++i)
        for (int q = 0; q < ngp; ++q)
          for (int j = 0; j < nen; ++j) gpfiber2[q](i, 0) += shapefcts[q](j) * fibernodes2[j](i);
      for (int q = 0; q < ngp; ++q) gpfiber2[q].Scale(1.0 / gpfiber2[q].Norm2());
    }
  }
  // case coordinate system is interpolated
  else
  {
    std::vector<LINALG::Matrix<nsd, 1>> cirgp(ngp);
    std::vector<LINALG::Matrix<nsd, 1>> tangp(ngp);
    std::vector<LINALG::Matrix<nsd, 1>> radgp(ngp);
    std::vector<double> helixgp(ngp);
    std::vector<double> transversegp(ngp);
    Epetra_SerialDenseVector tmpvector(3);
    Epetra_SerialDenseVector ftmp(3);
    Epetra_SerialDenseVector stmp(3);

    // interpolate circumferential and tangential directions to gauss points
    for (int q = 0; q < ngp; ++q)
    {
      LINALG::Matrix<nsd, 1> tmptan;
      for (int j = 0; j < nen; ++j)
      {
        for (int i = 0; i < nsd; ++i)
        {
          cirgp[q](i, 0) += shapefcts[q](j) * cirnodes[j](i);
          tmptan(i, 0) += shapefcts[q](j) * tannodes[j](i);
        }
        helixgp[q] += shapefcts[q](j) * helixnodes[j];
        transversegp[q] += shapefcts[q](j) * transversenodes[j];
      }
      cirgp[q].Scale(1.0 / cirgp[q].Norm2());
      tmptan.Scale(1.0 / tmptan.Norm2());

      // ensure that the circumferential direction is orthogonal to the tangential direction
      for (int i = 0; i < nsd; ++i)
        tangp[q](i, 0) = tmptan(i, 0) - tmptan.Dot(cirgp[q]) * cirgp[q](i, 0);
      tangp[q].Scale(1.0 / tangp[q].Norm2());
    }

    // get rad direction at gp
    for (int q = 0; q < ngp; ++q)
    {
      // Cross product and normalize
      tmpvector(0) = cirgp[q](1, 0) * tangp[q](2, 0) - cirgp[q](2, 0) * tangp[q](1, 0);
      tmpvector(1) = cirgp[q](2, 0) * tangp[q](0, 0) - cirgp[q](0, 0) * tangp[q](2, 0);
      tmpvector(2) = cirgp[q](0, 0) * tangp[q](1, 0) - cirgp[q](1, 0) * tangp[q](0, 0);

      for (int i = 0; i < nsd; ++i) radgp[q] = tmpvector(i) / tmpvector.Norm2();
    }

    double scale = PI / 180.;
    for (int q = 0; q < ngp; ++q)
    {
      // fiber 1 vector in cir,tan,rad cosy
      double fc = cos(helixgp[q] * scale) * cos(transversegp[q] * scale);
      double ft = sin(helixgp[q] * scale) * cos(transversegp[q] * scale);
      double fr = sin(transversegp[q] * scale);

      // fiber 2 vector in cir,tan,rad cosy
      double sc = cos((helixgp[q] + 90.0) * scale) * cos(transversegp[q] * scale);
      double st = sin((helixgp[q] + 90.0) * scale) * cos(transversegp[q] * scale);
      double sr = sin(transversegp[q] * scale);

      // fiber 1/2 vector in global cosy
      for (int i = 0; i < nsd; ++i)
      {
        ftmp(i) = fc * cirgp[q](i, 0) + ft * tangp[q](i, 0) + fr * radgp[q](i, 0);
        stmp(i) = sc * cirgp[q](i, 0) + st * tangp[q](i, 0) + sr * radgp[q](i, 0);
      }
      for (int i = 0; i < nsd; ++i)
      {
        gpfiber1[q](i, 0) = ftmp(i) / ftmp.Norm2();
        gpfiber2[q](i, 0) = stmp(i) / stmp.Norm2();
      }
    }
  }
  return;
}

template void DRT::ELEMENTS::UTILS::CalcR<DRT::Element::tet10>(
    const DRT::Element*, const std::vector<double>&, LINALG::Matrix<3, 3>&);
template void DRT::ELEMENTS::UTILS::NodalFiber<DRT::Element::tet10>(DRT::Node**,
    const std::vector<LINALG::Matrix<10, 1>>&, std::vector<LINALG::Matrix<3, 1>>&,
    std::vector<LINALG::Matrix<3, 1>>&);
