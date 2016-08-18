/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_poro_reac_ECM.cpp

\brief scatra_ele_calc_poro_reac_ECM.cpp

\level 2

<pre>
  \maintainer Anh-Tu Vuong
              vuong@lnm.mw.tum.de
              http://www.lnm.mw.tum.de
              089 - 289-15251
</pre>
 *----------------------------------------------------------------------*/

#include "scatra_ele_calc_poro_reac_ECM.H"
#include "scatra_ele_parameter_std.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

#include "../drt_mat/structporo.H"
#include "../drt_mat/structporo_reaction_ecm.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/scatra_mat_poro_ecm.H"
#include "../drt_mat/matlist_reactions.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::ScaTraEleCalcPoroReacECM(const int numdofpernode,const int numscal,const std::string& disname)
: DRT::ELEMENTS::ScaTraEleCalc<distype>::ScaTraEleCalc(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcPoro<distype>::ScaTraEleCalcPoro(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype>::ScaTraEleCalcAdvReac(numdofpernode,numscal,disname),
  DRT::ELEMENTS::ScaTraEleCalcPoroReac<distype>::ScaTraEleCalcPoroReac(numdofpernode,numscal,disname)
{
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerPoroReacECM(my::numscal_));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype> * DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcPoroReacECM* delete_me )
{
  static std::map<std::string,ScaTraEleCalcPoroReacECM<distype>* >  instances;

  if(delete_me == NULL)
  {
    if(instances.find(disname) == instances.end())
      instances[disname] = new ScaTraEleCalcPoroReacECM<distype>(numdofpernode,numscal,disname);
  }

  else
  {
    for( typename std::map<std::string,ScaTraEleCalcPoroReacECM<distype>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
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
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::Materials(
  const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
  const int                               k,        //!< id of current scalar
  double&                                 densn,    //!< density at t_(n)
  double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
  double&                                 densam,   //!< density at t_(n+alpha_M)
  double&                                 visc,         //!< fluid viscosity
  const int                               iquad         //!< id of current gauss point
  )
{
  switch(material->MaterialType())
  {
  case INPAR::MAT::m_scatra:
    pororeac::MatScaTra(material,k,densn,densnp,densam,visc,iquad);
    break;
  default:
    dserror("Material type %i is not supported",material->MaterialType());
   break;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 10/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
  // call poro base class to compute porosity
  poro::ComputePorosity(ele);

  // call base class
  advreac::GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

  return;
}

/*-----------------------------------------------------------------------------------------*
 |  get numcond, stoich list, reaction coefficient, couplingtpye from material  vuong 09/15 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::GetAdvancedReactionCoefficients(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int           iquad
  )
{
  const Teuchos::RCP<const MAT::MatListReactions>& actmat = Teuchos::rcp_dynamic_cast<const MAT::MatListReactions>(material);

  if(actmat==Teuchos::null)
    dserror("cast to MatListReactions failed");

  advreac::numcond_= actmat->NumReac();

  //We always have to reinitialize these vectors since our elements are singleton
  advreac::stoich_.resize(advreac::numcond_);
  advreac::couplingtype_.resize(advreac::numcond_);
  advreac::reaccoeff_.resize(advreac::numcond_);
  advreac::couprole_.resize(advreac::numcond_);
  advreac::reacstart_.resize(advreac::numcond_);

  {
    const Teuchos::RCP<const MAT::MatListReactions>& actmat = Teuchos::rcp_dynamic_cast<const MAT::MatListReactions>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<advreac::numcond_;++k)
    {
      int matid = actmat->ReacID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Teuchos::RCP<MAT::ScatraMatPoroECM> scatramat = Teuchos::rcp_dynamic_cast<MAT::ScatraMatPoroECM>(singlemat);

      if(scatramat != Teuchos::null)
      {
        Teuchos::RCP<MAT::StructPoroReactionECM> structmat = Teuchos::rcp_dynamic_cast<MAT::StructPoroReactionECM>(my::ele_->Material(1));
        if(structmat == Teuchos::null)
          dserror("cast to MAT::StructPoroReactionECM failed!");
        double structpot = ComputeStructChemPotential(structmat,iquad);

        scatramat->ComputeReacCoeff(structpot);
      }
    }
  }

  for (int i=0;i<advreac::numcond_;i++)
  {
    const int reacid = actmat->ReacID(i);
    const Teuchos::RCP<const MAT::ScatraReactionMat>& reacmat = Teuchos::rcp_dynamic_cast<const MAT::ScatraReactionMat>(actmat->MaterialById(reacid));

    advreac::stoich_[i] = *(reacmat->Stoich()); //get stoichometrie
    advreac::couplingtype_[i] = reacmat->Coupling(); //get coupling type
    advreac::reaccoeff_[i] = reacmat->ReacCoeff(); //get reaction coefficient
    advreac::couprole_[i] = *(reacmat->Couprole());
    advreac::reacstart_[i] = reacmat->ReacStart(); //get reaction start coefficient
  }
}

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    vuong 19/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<distype>::ComputeStructChemPotential(
    Teuchos::RCP<MAT::StructPoroReactionECM>& structmat,
    const int gp
  )
{

  //gauss point displacements
  LINALG::Matrix<my::nsd_,1> dispint(false);
  dispint.Multiply(my::edispnp_,my::funct_);

  // transposed jacobian "dX/ds"
  LINALG::Matrix<my::nsd_,my::nsd_> xjm0;
  xjm0.MultiplyNT(my::deriv_,poro::xyze0_);

  // inverse of transposed jacobian "ds/dX"
  LINALG::Matrix<my::nsd_,my::nsd_> xji0(true);
  xji0.Invert(xjm0);

  // inverse of transposed jacobian "ds/dX"
  const double det0= xjm0.Determinant();

  my::xjm_.MultiplyNT(my::deriv_,my::xyze_);
  const double det = my::xjm_.Determinant();

  // determinant of deformationgradient det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  const double J = det/det0;

  // ----------------------compute derivatives N_XYZ_ at gp w.r.t. material coordinates
  ///first derivatives of shape functions w.r.t. material coordinates
  LINALG::Matrix<my::nsd_,my::nen_> N_XYZ;
  N_XYZ.Multiply(xji0,my::deriv_);

  // -------------------------(material) deformation gradient F = d xyze_ / d XYZE = xyze_ * N_XYZ_^T
  static LINALG::Matrix<my::nsd_,my::nsd_>          defgrd(false);
  defgrd.MultiplyNT(my::xyze_,N_XYZ);

  // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
  static LINALG::Matrix<6,1> glstrain(true);
  glstrain.Clear();
  //if (kinemtype_ == INPAR::STR::kinem_nonlinearTotLag)
  {
    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<my::nsd_,my::nsd_> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);
    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    if(my::nsd_==3)
    {
      glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
      glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
      glstrain(3) = cauchygreen(0,1);
      glstrain(4) = cauchygreen(1,2);
      glstrain(5) = cauchygreen(2,0);
    }
    else if(my::nsd_==2)
    {
      glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
      glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
      glstrain(2) = 0.0;
      glstrain(3) = cauchygreen(0,1);
      glstrain(4) = 0.0;
      glstrain(5) = 0.0;
    }
  }

  //fluid pressure at gauss point
  const double pres = my::eprenp_.Dot(my::funct_);

  double pot = 0.0;

  structmat->ChemPotential(glstrain,poro::DiffManager()->GetPorosity(0),pres,J,my::eid_,pot,gp);

  return pot;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::line3>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::quad9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::nurbs9>;
//template class DRT::ELEMENTS::ScaTraEleCalcPoroReacECM<DRT::Element::nurbs27>;
