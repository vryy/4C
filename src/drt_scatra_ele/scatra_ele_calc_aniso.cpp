/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_aniso.cpp

 \brief

 <pre>
   Maintainer: Lasse Jagschies
               jagschies@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10365
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_aniso.H"

#include "scatra_ele.H"

#include "scatra_ele_parameter.H"
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
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcAniso<distype> * DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcAniso<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcAniso<distype>(numdofpernode,numscal);
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::ScaTraEleCalcAniso(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{
  // get diffusion manager for anisotropic diffusivity / diffusivities (in case of systems)
  my::diffmanager_ = Teuchos::rcp(new ScaTraEleDiffManagerAniso<my::nsd_>(my::numscal_));
}



/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  if (material->MaterialType() == INPAR::MAT::m_scatra_aniso)
    MatScaTraAniso(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
  else dserror("Material type is not supported");

  return;
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::MatScaTraAniso(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad   //!< id of current gauss point (default = -1)
  )
{
  // dynamic cast to anisotropic diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerAniso<my::nsd_> > diffmanageraniso = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerAniso<my::nsd_> >(diffmanager);

  int leleid = -1;
  if(DRT::Problem::Instance()->ProblemType()==prb_acou) leleid = DRT::Problem::Instance()->GetDis("scatra")->ElementColMap()->LID(my::eid_);

  const Teuchos::RCP<const MAT::ScatraMatAniso>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMatAniso>(material);

  // get constant diffusivity
  LINALG::Matrix<my::nsd_,my::nsd_> difftensor(true);
  LINALG::Matrix<3,1> diff = actmat->Diffusivity(leleid);

  for (int i=0; i<my::nsd_; i++)
    difftensor(i,i) = diff(i);

  diffmanageraniso->SetAnisotropicDiff(difftensor,k);

  // get reaction coefficient
  reamanager->SetReaCoeff(actmat->ReaCoeff(leleid),k);

  return;
} // ScaTraEleCalcAniso<distype>::MatScaTra


/*-------------------------------------------------------------------- *
 |  standard Galerkin diffusive term on right hand side     ehrl 11/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::CalcRHSDiff(
  Epetra_SerialDenseVector&                erhs,
  const int                                k,
  const double                             rhsfac,
  Teuchos::RCP<ScaTraEleDiffManager>       diffmanager,
  const LINALG::Matrix<my::nsd_,1>&        gradphi
  )
{
  // dynamic cast to anisotropic diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerAniso<my::nsd_> > diffmanageraniso = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerAniso<my::nsd_> >(diffmanager);

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    double laplawf(0.0);
    GetLaplacianWeakFormRHS(laplawf,diffmanageraniso->GetAnisotropicDiff(k), gradphi,vi);
    erhs[fvi] -= rhsfac*laplawf;
  }

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of diffusive element matrix                ehrl 11/13 |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAniso<distype>::CalcMatDiff(
  Epetra_SerialDenseMatrix&                emat,
  const int                                k,
  const double                             timefacfac,
  Teuchos::RCP<ScaTraEleDiffManager>       diffmanager
  )
{
  // dynamic cast to anisotropic diffusion manager
  Teuchos::RCP<ScaTraEleDiffManagerAniso<my::nsd_> > diffmanageraniso = Teuchos::rcp_dynamic_cast<ScaTraEleDiffManagerAniso<my::nsd_> >(diffmanager);

  for (int vi=0; vi<my::nen_; ++vi)
  {
    const int fvi = vi*my::numdofpernode_+k;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;
      double laplawf(0.0);
      GetLaplacianWeakForm(laplawf,diffmanageraniso->GetAnisotropicDiff(k),ui,vi);
      emat(fvi,fui) += timefacfac*laplawf;
    }
  }
  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcAniso<DRT::Element::nurbs27>;

