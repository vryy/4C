/*----------------------------------------------------------------------*/
/*! \file

\brief routines doing coordinate transformation between various coordinate systems
       a integration cell is lying in
--> THIS FUNCTIONALITY IS JUST USED IN COMBUST AND WILL LEAVE 4C SOON

\level 3

*----------------------------------------------------------------------*/


#ifndef FOUR_C_DISCRETIZATION_GEOMETRY_INTEGRATIONCELL_COORDTRAFO_HPP
#define FOUR_C_DISCRETIZATION_GEOMETRY_INTEGRATIONCELL_COORDTRAFO_HPP


#include "4C_config.hpp"

#include "4C_discretization_geometry_element_coordtrafo.hpp"
#include "4C_discretization_geometry_integrationcell.hpp"

FOUR_C_NAMESPACE_OPEN


namespace CORE::GEO
{
  //! map position from eta^boundary to xi^domain space
  inline void mapEtaBToXiD(const CORE::GEO::BoundaryIntCell& cell,
      const CORE::LINALG::Matrix<2, 1>& pos_eta_boundary,
      CORE::LINALG::Matrix<3, 1>& pos_xsi_domain)
  {
    // get cell node coordinates in xi_domain
    const CORE::LINALG::SerialDenseMatrix& xyze_cell(cell.cell_nodal_pos_xi_domain());
    CORE::GEO::elementToCurrentCoordinates(
        cell.Shape(), xyze_cell, pos_eta_boundary, pos_xsi_domain);
    return;
  }
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif
