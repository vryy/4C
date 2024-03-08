/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_aniso.cpp

\level 3

 *----------------------------------------------------------------------*/


#include "baci_scatra_ele_calc_aniso.hpp"

#include "baci_discretization_fem_general_utils_fem_shapefunctions.hpp"
#include "baci_discretization_fem_general_utils_gder2.hpp"
#include "baci_discretization_fem_general_utils_nurbs_shapefunctions.hpp"
#include "baci_discretization_geometry_position_array.hpp"
#include "baci_global_data.hpp"
#include "baci_lib_condition_utils.hpp"
#include "baci_lib_discret.hpp"
#include "baci_lib_element.hpp"
#include "baci_lib_utils.hpp"
#include "baci_mat_list.hpp"
#include "baci_mat_newtonianfluid.hpp"
#include "baci_mat_scatra_aniso.hpp"
#include "baci_nurbs_discret_nurbs_utils.hpp"
#include "baci_scatra_ele.hpp"
#include "baci_scatra_ele_parameter_std.hpp"
#include "baci_scatra_ele_parameter_timint.hpp"
#include "baci_utils_singleton_owner.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::Instance(
    const int numdofpernode, const int numscal, const std::string& disname)
{
  static auto singleton_map = CORE::UTILS::MakeSingletonMap<std::string>(
      [](const int numdofpernode, const int numscal, const std::string& disname)
      {
        return std::unique_ptr<ScaTraEleCalcAniso<distype, probdim>>(
            new ScaTraEleCalcAniso<distype, probdim>(numdofpernode, numscal, disname));
      });

  return singleton_map[disname].Instance(
      CORE::UTILS::SingletonAction::create, numdofpernode, numscal, disname);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::ScaTraEleCalcAniso(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // get diffusion manager for anisotropic diffusivity / diffusivities (in case of systems)
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerAniso<nsd_>(my::numscal_));
}



/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::Materials(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point
)
{
  if (material->MaterialType() == INPAR::MAT::m_scatra_aniso)
    MatScaTraAniso(material, k, densn, densnp, densam, visc, iquad);
  else
    dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::MatScaTraAniso(
    const Teuchos::RCP<const MAT::Material> material,  //!< pointer to current material
    const int k,                                       //!< id of current scalar
    double& densn,                                     //!< density at t_(n)
    double& densnp,                                    //!< density at t_(n+1) or t_(n+alpha_F)
    double& densam,                                    //!< density at t_(n+alpha_M)
    double& visc,                                      //!< fluid viscosity
    const int iquad                                    //!< id of current gauss point (default = -1)
)
{
  const Teuchos::RCP<const MAT::ScatraMatAniso>& actmat =
      Teuchos::rcp_dynamic_cast<const MAT::ScatraMatAniso>(material);

  // get constant diffusivity
  CORE::LINALG::Matrix<nsd_, nsd_> difftensor(true);
  CORE::LINALG::Matrix<3, 1> diff = actmat->Diffusivity();

  for (unsigned i = 0; i < nsd_; i++) difftensor(i, i) = diff(i);

  DiffManager()->SetAnisotropicDiff(difftensor, k);

  return;
}  // ScaTraEleCalcAniso<distype>::MatScaTra


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::CalcRHSDiff(
    CORE::LINALG::SerialDenseVector& erhs, const int k, const double rhsfac)
{
  const CORE::LINALG::Matrix<nsd_, 1>& gradphi = my::scatravarmanager_->GradPhi(k);

  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf, DiffManager()->GetAnisotropicDiff(k), gradphi, vi);
    erhs[fvi] -= rhsfac * laplawf;
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::CalcMatDiff(
    CORE::LINALG::SerialDenseMatrix& emat, const int k, const double timefacfac)
{
  for (unsigned vi = 0; vi < nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < nen_; ++ui)
    {
      const int fui = ui * my::numdofpernode_ + k;
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf, DiffManager()->GetAnisotropicDiff(k), ui, vi);
      emat(fvi, fui) += timefacfac * laplawf;
    }
  }
  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<CORE::FE::CellType::nurbs27>;

BACI_NAMESPACE_CLOSE
