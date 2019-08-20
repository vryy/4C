/*----------------------------------------------------------------------*/
/*! \file

\brief routines doing coordinate transformation between various coordinate systems
       a integration cell is lying in
--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE BACI SOON

\level 3

\maintainer Martin Kronbichler
*----------------------------------------------------------------------*/


#ifndef INTEGRATIONCELL_COORDTRAFO_H
#define INTEGRATIONCELL_COORDTRAFO_H


#include "element_coordtrafo.H"
#include "integrationcell.H"


namespace GEO
{
  /////////////////////  Domain Integration Cell ////////////////////////////////

  //! map position from eta to xi space
  inline void mapEtaToXi3D(const GEO::DomainIntCell& cell,
      const LINALG::Matrix<3, 1>& pos_eta_domain, LINALG::Matrix<3, 1>& pos_xsi_domain)
  {
    // get cell node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(cell.CellNodalPosXiDomain());
    GEO::elementToCurrentCoordinates(cell.Shape(), xyze_cell, pos_eta_domain, pos_xsi_domain);
    return;
  }


  //! compute determinant from eta to xi space
  template <DRT::Element::DiscretizationType distype>
  static inline double detEtaToXi3DT(
      const LINALG::SerialDenseMatrix& xyze_cell, const LINALG::Matrix<3, 1>& pos_eta_domain)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    // get shape function derivative
    static LINALG::Matrix<3, numNodes> deriv;
    DRT::UTILS::shape_function_3D_deriv1(
        deriv, pos_eta_domain(0), pos_eta_domain(1), pos_eta_domain(2), distype);

    // get jacobian matrix ...
    static LINALG::Matrix<3, 3> xjm;
    xjm.Clear();

    for (int inode = 0; inode < numNodes; ++inode)
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) xjm(i, j) += deriv(i, inode) * xyze_cell(j, inode);

    // ...and return its determinant
    return xjm.Determinant();
  }


  //! get determinant between eta and xi space
  inline double detEtaToXi3D(
      const GEO::DomainIntCell& cell, const LINALG::Matrix<3, 1>& pos_eta_domain)
  {
    // get node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(cell.CellNodalPosXiDomain());
    switch (cell.Shape())
    {
      case DRT::Element::tet4:
        return detEtaToXi3DT<DRT::Element::tet4>(xyze_cell, pos_eta_domain);
      case DRT::Element::tet10:
        return detEtaToXi3DT<DRT::Element::tet10>(xyze_cell, pos_eta_domain);
      case DRT::Element::hex8:
        return detEtaToXi3DT<DRT::Element::hex8>(xyze_cell, pos_eta_domain);
      case DRT::Element::hex20:
        return detEtaToXi3DT<DRT::Element::hex20>(xyze_cell, pos_eta_domain);
      case DRT::Element::hex27:
        return detEtaToXi3DT<DRT::Element::hex27>(xyze_cell, pos_eta_domain);
      case DRT::Element::wedge6:
        return detEtaToXi3DT<DRT::Element::wedge6>(xyze_cell, pos_eta_domain);
      case DRT::Element::wedge15:
        return detEtaToXi3DT<DRT::Element::wedge15>(xyze_cell, pos_eta_domain);
      case DRT::Element::pyramid5:
        return detEtaToXi3DT<DRT::Element::pyramid5>(xyze_cell, pos_eta_domain);
      default:
        std::cout << DistypeToString(cell.Shape()) << std::endl;
        dserror("add your 3D distype to this switch!");
    }
    return -1.0;
  }

  //! get determinant between eta and xi space
  inline double detXiD_over_EtaD(
      const GEO::DomainIntCell& cell, const LINALG::Matrix<3, 1>& pos_eta_domain)
  {
    return GEO::detEtaToXi3D(cell, pos_eta_domain);
  }

  /////////////////////  Boundary Integration Cell ////////////////////////////////

  //! map position from eta^boundary to xi^boundary space
  template <DRT::Element::DiscretizationType distype, class M>
  inline void mapEtaBToXiBT(const M& xyze_cell, const LINALG::Matrix<2, 1>& pos_eta_boundary,
      LINALG::Matrix<2, 1>& pos_xsi_boundary)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // create shape function vectors
    static LINALG::Matrix<numNodes, 1> funct;
    DRT::UTILS::shape_function_2D(funct, pos_eta_boundary(0), pos_eta_boundary(1), distype);

    // does this make sense
    // translate position into from cell coordinates to element coordinates
    pos_xsi_boundary.Clear();
    for (int j = 0; j < numNodes; ++j)
      for (int i = 0; i < 2; ++i) pos_xsi_boundary(i) += xyze_cell(i, j) * funct(j);

    return;
  }


  //! map position from eta^boundary to xi^boundary space
  inline void mapEtaBToXiB(const GEO::BoundaryIntCell& cell,
      const LINALG::Matrix<2, 1>& pos_eta_boundary, LINALG::Matrix<2, 1>& pos_xsi_boundary)
  {
    // get node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(
        cell.CellNodalPosXiBoundary());  // TODO check was xi domain before
    switch (cell.Shape())
    {
      case DRT::Element::tri3:
        mapEtaBToXiBT<DRT::Element::tri3>(xyze_cell, pos_eta_boundary, pos_xsi_boundary);
        break;
      case DRT::Element::tri6:
        mapEtaBToXiBT<DRT::Element::tri6>(xyze_cell, pos_eta_boundary, pos_xsi_boundary);
        break;
      case DRT::Element::quad4:
        mapEtaBToXiBT<DRT::Element::quad4>(xyze_cell, pos_eta_boundary, pos_xsi_boundary);
        break;
      case DRT::Element::quad8:
        mapEtaBToXiBT<DRT::Element::quad8>(xyze_cell, pos_eta_boundary, pos_xsi_boundary);
        break;
      case DRT::Element::quad9:
        mapEtaBToXiBT<DRT::Element::quad9>(xyze_cell, pos_eta_boundary, pos_xsi_boundary);
        break;
      default:
        std::cout << DistypeToString(cell.Shape()) << std::endl;
        dserror("add your 2D distype to this switch!");
    }
    return;
  }


  //! map position from eta^boundary to xi^domain space
  inline void mapEtaBToXiD(const GEO::BoundaryIntCell& cell,
      const LINALG::Matrix<2, 1>& pos_eta_boundary, LINALG::Matrix<3, 1>& pos_xsi_domain)
  {
    // get cell node coordinates in xi_domain
    const LINALG::SerialDenseMatrix& xyze_cell(cell.CellNodalPosXiDomain());
    GEO::elementToCurrentCoordinates(cell.Shape(), xyze_cell, pos_eta_boundary, pos_xsi_domain);
    return;
  }


  //! get determinant between eta^boundary and xi^boundary space
  template <DRT::Element::DiscretizationType distype>
  static inline double detEtaBToXiBT(
      const LINALG::SerialDenseMatrix& xyze_cell, const LINALG::Matrix<2, 1>& pos_eta_boundary)
  {
    const int numNodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    // get shape function derivative
    static LINALG::Matrix<2, numNodes> deriv;
    DRT::UTILS::shape_function_2D_deriv1(deriv, pos_eta_boundary(0), pos_eta_boundary(1), distype);

    // get jacobian matrix ...
    static LINALG::Matrix<2, 2> xjm;
    xjm.Clear();

    for (int inode = 0; inode < numNodes; ++inode)
      for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j) xjm(i, j) += deriv(i, inode) * xyze_cell(j, inode);

    // ...and return its determinant
    return xjm.Determinant();
  }


  //! get determinant between eta^boundary and xi^boundary space
  inline double detEtaBToXiB(
      const GEO::BoundaryIntCell& cell, const LINALG::Matrix<2, 1>& pos_eta_boundary)
  {
    // get node coordinates
    const LINALG::SerialDenseMatrix& xyze_cell(cell.CellNodalPosXiBoundary());
    switch (cell.Shape())
    {
      case DRT::Element::tri3:
        return detEtaBToXiBT<DRT::Element::tri3>(xyze_cell, pos_eta_boundary);
      case DRT::Element::tri6:
        return detEtaBToXiBT<DRT::Element::tri6>(xyze_cell, pos_eta_boundary);
      case DRT::Element::quad4:
        return detEtaBToXiBT<DRT::Element::quad4>(xyze_cell, pos_eta_boundary);
      case DRT::Element::quad8:
        return detEtaBToXiBT<DRT::Element::quad8>(xyze_cell, pos_eta_boundary);
      case DRT::Element::quad9:
        return detEtaBToXiBT<DRT::Element::quad9>(xyze_cell, pos_eta_boundary);
      default:
        std::cout << DistypeToString(cell.Shape()) << std::endl;
        dserror("add your 2D distype to this switch!");
    }
    return -1.0;
  }

  //! get determinant: det( xi^boundary / eta^boundary)
  inline double detXiB_over_EtaB(
      const GEO::BoundaryIntCell& cell, const LINALG::Matrix<2, 1>& pos_eta_boundary)
  {
    return detEtaBToXiB(cell, pos_eta_boundary);
  }

}  // namespace GEO

#endif
