/*----------------------------------------------------------------------*/
/*! \file

\brief Scalar types for different kind of geometry pairs.

\level 1
*/
// End doxygen header.


#ifndef BACI_GEOMETRY_PAIR_SCALAR_TYPES_H
#define BACI_GEOMETRY_PAIR_SCALAR_TYPES_H


#include "baci_config.H"

#include "baci_utils_fad.H"

BACI_NAMESPACE_OPEN

namespace GEOMETRYPAIR
{
  //! Scalar type to be used for line to volume pairs.
  template <typename line, typename volume>
  using line_to_volume_scalar_type = typename CORE::FADUTILS::HigherOrderFadType<1,
      Sacado::Fad::SLFad<double, line::n_dof_ + volume::n_dof_>>::type;

  //! Scalar type to be used for line to surface pairs without averaged current normals.
  template <typename line, typename surface>
  using line_to_surface_scalar_type = typename CORE::FADUTILS::HigherOrderFadType<1,
      Sacado::ELRFad::SLFad<double, line::n_dof_ + surface::n_dof_>>::type;

  //! First order FAD scalar type to be used for line to surface patch pairs with averaged current
  //! normals.
  using line_to_surface_patch_scalar_type_1st_order =
      typename CORE::FADUTILS::HigherOrderFadType<1, Sacado::ELRFad::DFad<double>>::type;

  //! First order FAD scalar type to be used for line to surface patch pairs with nurbs
  //! discretization.
  template <typename line, typename surface>
  using line_to_surface_patch_scalar_type_fixed_size_1st_order =
      typename CORE::FADUTILS::HigherOrderFadType<1,
          Sacado::ELRFad::SLFad<double, line::n_dof_ + surface::n_dof_>>::type;

  //! Second order FAD scalar type to be used for line to surface patch pairs with averaged current
  //! normals.
  using line_to_surface_patch_scalar_type =
      typename CORE::FADUTILS::HigherOrderFadType<2, Sacado::ELRFad::DFad<double>>::type;

  //! Scalar type to be used for line to surface patch pairs with nurbs discretization.
  template <typename line, typename surface>
  using line_to_surface_patch_scalar_type_fixed_size =
      typename CORE::FADUTILS::HigherOrderFadType<2,
          Sacado::ELRFad::SLFad<double, line::n_dof_ + surface::n_dof_>>::type;
}  // namespace GEOMETRYPAIR

BACI_NAMESPACE_CLOSE

#endif
