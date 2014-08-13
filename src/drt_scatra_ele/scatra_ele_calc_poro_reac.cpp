/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_poro_reac.cpp

 \brief

 <pre>
   Maintainer: Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_poro_reac.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/scatra_mat.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal),
    DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode,numscal),
    DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal)
{

}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype> * DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcPoroReac<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcPoroReac<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}

/*-------------------------------------------------------------------------------*
 |  Evaluate including check for ScaTra in Porous medium with reaction   vuong 08/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::Evaluate(
  DRT::ELEMENTS::Transport*  ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  const std::vector<int>&    lm,
  Epetra_SerialDenseMatrix&  elemat1_epetra,
  Epetra_SerialDenseMatrix&  elemat2_epetra,
  Epetra_SerialDenseVector&  elevec1_epetra,
  Epetra_SerialDenseVector&  elevec2_epetra,
  Epetra_SerialDenseVector&  elevec3_epetra
  )
{

  if (not advreac::isinit_)
  {
    advreac::iscoupled_=advreac::IsCoupledAndRead(discretization);
    advreac::isinit_=true;
  }
  return my::Evaluate(
      ele,
      params,
      discretization,
      lm,
      elemat1_epetra,
      elemat2_epetra,
      elevec1_epetra,
      elevec2_epetra,
      elevec3_epetra
      );

  //Todo: extend for evaluation of coupling terms
}


/*----------------------------------------------------------------------*
 |  Material ScaTra                                          thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::MatScaTra(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  Teuchos::RCP<ScaTraEleDiffManager>      diffmanager,  //!< diffusion manager handling diffusivity / diffusivities (in case of systems) or (thermal conductivity/specific heat) in case of loma
  Teuchos::RCP<ScaTraEleReaManager>       reamanager,   //!< reaction manager
  double&                                 visc,     //!< fluid viscosity
  const int                               iquad     //!< id of current gauss point
  )
{
  if(iquad==-1)
    dserror("no gauss point given for evaluation of scatra material. Check your input file.");

  //get porosity from structure material
  const double porosity = poro::GetPorosityAtGP(iquad);

  const Teuchos::RCP<const MAT::ScatraMat>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

  // set diffusivity (scaled with porosity)
  poro::SetDiffusivity(actmat,k,diffmanager,porosity);

  // set/calculate reaction coefficient (scaled with porosity)
  if (not advreac::iscoupled_)
    poro::SetReaCoefficient(actmat,k,reamanager,porosity);
  else
    advreac::SetAdvancedReactionTerms(reamanager,k,porosity);

  // set densities (scaled with porosity)
  poro::SetDensities(porosity,densn,densnp,densam);

  return;
} // ScaTraEleCalcPoroReac<distype>::MatScaTra





// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs27>;

