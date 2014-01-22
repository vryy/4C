/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.cpp

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
   Maintainer: Moritz Thon
               thon@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_advanced_reaction.H"

#include "scatra_ele_parameter.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/biofilm.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype> * DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  bool create )
{
  static ScaTraEleCalcAdvReac<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new ScaTraEleCalcAdvReac<distype>(numdofpernode,numscal);
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
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal)
  : DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal)
{

}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::Materials(
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
  if (material->MaterialType() == INPAR::MAT::m_scatra)
    my::MatScaTra(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
  else if (material->MaterialType() == INPAR::MAT::m_biofilm)
    MatBioFilm(material,k,densn,densnp,densam,diffmanager,reamanager,visc,iquad);
  else dserror("Material type %i is not supported",material->MaterialType());

  return;
}

/*----------------------------------------------------------------------*
 |  Material BioFilm                                         ehrl 11/13 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::MatBioFilm(
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
  const Teuchos::RCP<const MAT::Biofilm>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Biofilm>(material);

  // get constant diffusivity
  diffmanager->SetIsotropicDiff(actmat->Diffusivity(),k);

  // get substrate concentration at n+1 or n+alpha_F at integration point
  const double csnp = my::funct_.Dot(my::ephinp_[k]);

  // set reaction coefficient
  reamanager->SetReaCoeff(actmat->ComputeReactionCoeff(csnp),k);
  // set derivative of reaction coefficient
  reamanager->SetReaCoeffDeriv(actmat->ComputeReactionCoeffDeriv(csnp),k);

  // set density at various time steps and density gradient factor to 1.0/0.0
  densn      = 1.0;
  densnp     = 1.0;
  densam     = 1.0;

  return;
}

/*------------------------------------------------------------------- *
 |  calculation of reactive element matrix                ehrl 11/13  |
 *--------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::CalcMatReact(
  Epetra_SerialDenseMatrix&          emat,
  const int                          k,
  const double                       timefacfac,
  const double                       timetaufac,
  const double                       taufac,
  const double                       densnp,
  const double                       phinp,
  Teuchos::RCP<ScaTraEleReaManager>  reamanager,
  const LINALG::Matrix<my::nen_,1>&      conv,
  const LINALG::Matrix<my::nen_,1>&      sgconv,
  const LINALG::Matrix<my::nen_,1>&      diff
  )
{
  my::CalcMatReact(emat,k,timefacfac,timetaufac,taufac,densnp,phinp,reamanager,conv,sgconv,diff);

  const double fac_reac        = timefacfac*densnp*reamanager->GetReaCoeffDeriv(k)*phinp;
  const double timetaufac_reac = timetaufac*densnp*reamanager->GetReaCoeffDeriv(k)*phinp;

  //----------------------------------------------------------------
  // standard Galerkin reactive term
  //----------------------------------------------------------------
  for (int vi=0; vi<my::nen_; ++vi)
  {
    const double v = fac_reac*my::funct_(vi);
    const int fvi = vi*my::numdofpernode_+k;

    for (int ui=0; ui<my::nen_; ++ui)
    {
      const int fui = ui*my::numdofpernode_+k;

      emat(fvi,fui) += v*my::funct_(ui);
    }
  }

  //----------------------------------------------------------------
  // stabilization of reactive term
  //----------------------------------------------------------------
  if(my::scatrapara_->StabType()!=INPAR::SCATRA::stabtype_no_stabilization)
  {
    double densreataufac = timetaufac_reac*densnp;
    // convective stabilization of reactive term (in convective form)
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = densreataufac*(conv(vi)+sgconv(vi)+my::scatrapara_->USFEMGLSFac()*1.0/my::scatraparatimint_->TimeFac()*my::funct_(vi));
      const int fvi = vi*my::numdofpernode_+k;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;

        emat(fvi,fui) += v*my::funct_(ui);
      }
    }

    if (my::use2ndderiv_)
    {
      // diffusive stabilization of reactive term
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double v = my::scatrapara_->USFEMGLSFac()*timetaufac_reac*diff(vi);
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+k;

          emat(fvi,fui) -= v*my::funct_(ui);
        }
      }
    }

    //----------------------------------------------------------------
    // reactive stabilization
    //----------------------------------------------------------------
    densreataufac = my::scatrapara_->USFEMGLSFac()*timetaufac_reac*densnp;

    // reactive stabilization of convective (in convective form) and reactive term
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = densreataufac*my::funct_(vi);
      const int fvi = vi*my::numdofpernode_+k;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+k;

        emat(fvi,fui) += v*(conv(ui)+reamanager->GetReaCoeff(k)*my::funct_(ui));
      }
    }

    if (my::use2ndderiv_)
    {
      // reactive stabilization of diffusive term
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double v = my::scatrapara_->USFEMGLSFac()*timetaufac_reac*my::funct_(vi);
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+k;

          emat(fvi,fui) -= v*diff(ui);
        }
      }
    }


    if (not my::scatraparatimint_->IsStationary())
    {
      // reactive stabilization of transient term
      for (int vi=0; vi<my::nen_; ++vi)
      {
        const double v = my::scatrapara_->USFEMGLSFac()*taufac*densnp*reamanager->GetReaCoeff(k)*densnp*my::funct_(vi);
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+k;

          emat(fvi,fui) += v*my::funct_(ui);
        }
      }

      if (my::use2ndderiv_ and reamanager->GetReaCoeff(k)!=0.0)
        dserror("Second order reactive stabilization is not fully implemented!! ");
    }
  }

  return;
}


// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs27>;
