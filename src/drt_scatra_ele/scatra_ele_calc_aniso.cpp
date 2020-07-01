/*----------------------------------------------------------------------*/
/*! \file

\brief scatra_ele_calc_aniso.cpp

\level 3

 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_aniso.H"

#include "scatra_ele.H"

#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_lib/drt_globalproblem.H"  // for time curve in body force
#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_mat/matlist.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/scatra_mat_aniso.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>*
DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::Instance(const int numdofpernode,
    const int numscal, const std::string& disname, const ScaTraEleCalcAniso* delete_me)
{
  static std::map<std::string, ScaTraEleCalcAniso<distype, probdim>*> instances;

  if (delete_me == NULL)
  {
    if (instances.find(disname) == instances.end())
      instances[disname] =
          new ScaTraEleCalcAniso<distype, probdim>(numdofpernode, numscal, disname);
  }

  else
  {
    for (typename std::map<std::string, ScaTraEleCalcAniso<distype, probdim>*>::iterator i =
             instances.begin();
         i != instances.end(); ++i)
      if (i->second == delete_me)
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[disname];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(0, 0, "", this);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::ScaTraEleCalcAniso(
    const int numdofpernode, const int numscal, const std::string& disname)
    : DRT::ELEMENTS::ScaTraEleCalc<distype, probdim>::ScaTraEleCalc(numdofpernode, numscal, disname)
{
  // get diffusion manager for anisotropic diffusivity / diffusivities (in case of systems)
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerAniso<my::nsd_>(my::numscal_));
}



/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
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
template <DRT::Element::DiscretizationType distype, int probdim>
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
  LINALG::Matrix<my::nsd_, my::nsd_> difftensor(true);
  LINALG::Matrix<3, 1> diff = actmat->Diffusivity();

  for (unsigned i = 0; i < my::nsd_; i++) difftensor(i, i) = diff(i);

  DiffManager()->SetAnisotropicDiff(difftensor, k);

  return;
}  // ScaTraEleCalcAniso<distype>::MatScaTra


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::CalcRHSDiff(
    Epetra_SerialDenseVector& erhs, const int k, const double rhsfac)
{
  const LINALG::Matrix<my::nsd_, 1>& gradphi = my::scatravarmanager_->GradPhi(k);

  for (unsigned vi = 0; vi < my::nen_; ++vi)
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
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype, probdim>::CalcMatDiff(
    Epetra_SerialDenseMatrix& emat, const int k, const double timefacfac)
{
  for (unsigned vi = 0; vi < my::nen_; ++vi)
  {
    const int fvi = vi * my::numdofpernode_ + k;

    for (unsigned ui = 0; ui < my::nen_; ++ui)
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
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::nurbs27>;
