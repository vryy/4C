/*----------------------------------------------------------------------*/
/*! \file

\brief Scalar types for different kind of geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef FOUR_C_GEOMETRY_PAIR_SCALAR_TYPES_HPP
#define FOUR_C_GEOMETRY_PAIR_SCALAR_TYPES_HPP


#include "4C_config.hpp"

#include "4C_utils_fad.hpp"

FOUR_C_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  //! Scalar type to be used for line to volume pairs.
  template <typename Line, typename Volume>
  using line_to_volume_scalar_type = typename Core::FADUtils::HigherOrderFadType<1,
      Sacado::Fad::SLFad<double, Line::n_dof_ + Volume::n_dof_>>::type;

  //! Scalar type to be used for line to surface pairs without averaged current normals.
  template <typename Line, typename Surface>
  using line_to_surface_scalar_type = typename Core::FADUtils::HigherOrderFadType<1,
      Sacado::ELRFad::SLFad<double, Line::n_dof_ + Surface::n_dof_>>::type;

  //! First order FAD scalar type to be used for line to surface patch pairs with averaged current
  //! normals.
  using line_to_surface_patch_scalar_type_1st_order =
      typename Core::FADUtils::HigherOrderFadType<1, Sacado::ELRFad::DFad<double>>::type;

  //! First order FAD scalar type to be used for line to surface patch pairs with nurbs
  //! discretization.
  template <typename Line, typename Surface>
  using line_to_surface_patch_scalar_type_fixed_size_1st_order =
      typename Core::FADUtils::HigherOrderFadType<1,
          Sacado::ELRFad::SLFad<double, Line::n_dof_ + Surface::n_dof_>>::type;

  //! Second order FAD scalar type to be used for line to surface patch pairs with averaged current
  //! normals.
  using line_to_surface_patch_scalar_type =
      typename Core::FADUtils::HigherOrderFadType<2, Sacado::ELRFad::DFad<double>>::type;

  //! Scalar type to be used for line to surface patch pairs with nurbs discretization.
  template <typename Line, typename Surface>
  using line_to_surface_patch_scalar_type_fixed_size =
      typename Core::FADUtils::HigherOrderFadType<2,
          Sacado::ELRFad::SLFad<double, Line::n_dof_ + Surface::n_dof_>>::type;
}  // namespace GEOMETRYPAIR

FOUR_C_NAMESPACE_CLOSE

#endif
