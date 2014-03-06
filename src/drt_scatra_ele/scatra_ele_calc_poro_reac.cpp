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

// TODO: Do we need this part?? found in ScaTraEleCalc<distype>::MatScaTra
//  // in case of multifractal subgrid-scales, read Schmidt number
//  if (scatrapara_->TurbModel() == INPAR::FLUID::multifractal_subgrid_scales or scatrapara_->RBSubGrVel()
//      or scatrapara_->TurbModel() == INPAR::FLUID::dynamic_smagorinsky)
//  {
//    //access fluid discretization
//    Teuchos::RCP<DRT::Discretization> fluiddis = Teuchos::null;
//    fluiddis = DRT::Problem::Instance()->GetDis("fluid");
//    //get corresponding fluid element (it has the same global ID as the scatra element)
//    DRT::Element* fluidele = fluiddis->gElement(eid_);
//    if (fluidele == NULL)
//      dserror("Fluid element %i not on local processor", eid_);
//
//    // get fluid material
//    Teuchos::RCP<MAT::Material> fluidmat = fluidele->Material();
//    if(fluidmat->MaterialType() != INPAR::MAT::m_fluid)
//      dserror("Invalid fluid material for passive scalar transport in turbulent flow!");
//
//    const Teuchos::RCP<const MAT::NewtonianFluid>& actfluidmat = Teuchos::rcp_dynamic_cast<const MAT::NewtonianFluid>(fluidmat);
//
//    // get constant dynamic viscosity
//    visc = actfluidmat->Viscosity();
//    densn = actfluidmat->Density();
//    densnp = actfluidmat->Density();
//    densam = actfluidmat->Density();
//
//    if (densam != 1.0 or densnp != 1.0 or densn != 1.0)
//       dserror("Check your diffusivity! Dynamic diffusivity required!");
//   }

    if(iquad==-1)
      dserror("no gauss point given for evaluation of scatra material. Check your input file.");

    //access structure discretization
    Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
    structdis = DRT::Problem::Instance()->GetDis("structure");
    //get corresponding fluid element (it has the same global ID as the scatra element)
    DRT::Element* structele = structdis->gElement(my::eid_);
    if (structele == NULL)
      dserror("Structure element %i not on local processor", my::eid_);

    const Teuchos::RCP<const MAT::StructPoro>& structmat
              = Teuchos::rcp_dynamic_cast<const MAT::StructPoro>(structele->Material());
    if(structmat == Teuchos::null)
      dserror("invalid structure material for poroelasticity");

    const double porosity   = structmat->GetPorosityAtGP(iquad);

    const Teuchos::RCP<const MAT::ScatraMat>& actmat
      = Teuchos::rcp_dynamic_cast<const MAT::ScatraMat>(material);

    dsassert(my::numdofpernode_==1,"more than 1 dof per node for SCATRA material");

    // set diffusivity (scaled with porosity)
    diffmanager->SetIsotropicDiff(actmat->Diffusivity()*porosity,k);

    // set/calculate reaction coefficient (scaled with porosity)
    if (not advreac::iscoupled_)
      reamanager->SetReaCoeff(actmat->ReaCoeff()*porosity,k); // set reaction coefficient (scaled with porosity)
    else
    {
      // dynamic cast to Advanced_Reaction-specific reaction manager
      Teuchos::RCP<ScaTraEleReaManagerAdvReac> reamanageradvreac = Teuchos::rcp_dynamic_cast<ScaTraEleReaManagerAdvReac>(reamanager);

      reamanageradvreac->SetReaBodyForce( advreac::CalcReaBodyForceTerm(k)*porosity ,k);
      reamanageradvreac->SetReaCoeff( advreac::CalcReaCoeff(k)*porosity ,k);
      for (int j=0; j<my::numscal_ ;j++)
      {
        reamanageradvreac->SetReaBodyForceDerivMatrix( advreac::CalcReaBodyForceDerivMatrix(k,j)*porosity ,k,j );
        reamanager->SetReaCoeffDerivMatrix( advreac::CalcReaCoeffDerivMatrix(k,j)*porosity ,k,j );
      }
    }
    // set densities (scaled with porosity)
    densn = porosity;
    densnp = porosity;
    densam = porosity;

    return;
  } // ScaTraEleCalcPoroReac<distype>::MatScaTra





// template classes

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::line2>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex20>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::wedge6>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReac<DRT::Element::nurbs27>;

