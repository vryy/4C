/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_aniso.cpp

\level 3

 *----------------------------------------------------------------------*/


#include "4C_scatra_ele_calc_aniso.hpp"

#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_fem_shapefunctions.hpp"
#include "4C_fem_general_utils_gder2.hpp"
#include "4C_fem_general_utils_nurbs_shapefunctions.hpp"
#include "4C_fem_geometry_position_array.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_newtonianfluid.hpp"
#include "4C_nurbs_discret_nurbs_utils.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_parameter_std.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_utils_singleton_owner.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>*
Discret::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = Core::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcAniso<distype, probdim>>(
            new ScaTraEleCalcAniso<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      Core::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
Discret::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::ScaTraEleCalcAniso(
    const int numdofpernode, const int numscal, const std::string& disname)
    : Discret::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(
          numdofpernode, numscal, disname)
{
  // get diffusion manager for anisotropic diffusivity / diffusivities (in case of systems)
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerAniso<nsd_>(my::numscal_));
}


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::calc_rhs_diff(
    Core::LinAlg::SerialDenseVector& erhs, const int k, const double rhsfac)
{
  const Core::LinAlg::Matrix<nsd_, 1>& gradphi = my::scatravarmanager_->GradPhi(k);

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    double laplawf(0.0);
    get_laplacian_weak_form_rhs(laplawf, diff_manager()->GetAnisotropicDiff(k), gradphi, vi);
    erhs[fvi] -= rhsfac * laplawf;
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
void Discret::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::calc_mat_diff(
    Core::LinAlg::SerialDenseMatrix& emat, const int k, const double timefacfac)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + k;
      double laplawf(0.0);
      get_laplacian_weak_form(laplawf, diff_manager()->GetAnisotropicDiff(k), ui, vi);
      emat(fvi, fui) += timefacfac * laplawf;
    }
  }
  return;
}


// template classes

// 1D elements
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::line2, 1>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::line2, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::line2, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::line3, 1>;

// 2D elements
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::tri3, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::tri3, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::tri6, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::quad4, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::quad4, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::quad8>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::quad9, 2>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::nurbs9, 2>;

// 3D elements
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::hex8, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::hex20>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::hex27, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::tet4, 3>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::tet10, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::wedge6>;
template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::pyramid5, 3>;
// template class Discret::ELEMENTS::ScaTraEleCalcAniso<Core::FE::CellType::nurbs27>;

FOUR_C_NAMESPACE_CLOSE
