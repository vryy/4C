/*----------------------------------------------------------------------*/
/*!
 \file scatra_ele_calc_advanced_reaction.cpp

 \brief main file containing routines for calculation of scatra element with advanced reaction terms


 <pre>
 \level 2

   \maintainer Moritz Thon
               thon@mhpc.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-10364
 </pre>
 *----------------------------------------------------------------------*/


#include "scatra_ele_calc_advanced_reaction.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_element.H"

//MATERIALS
#include "../drt_mat/biofilm.H"
#include "../drt_mat/scatra_growth_scd.H"
#include "../drt_mat/growth_scd.H"
#include "../drt_mat/growth_law.H"
#include "../drt_mat/matlist_reactions.H"
#include "../drt_mat/scatra_mat.H"
#include "../drt_mat/matlist.H"

  //! note for the implementation of the homogenous scatra coupling:
  //! assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0
  //!
  //! if we assume the reaction is depending on the product of all
  //! reactants (this corresponds to couplingtype "simple_multiplicative"),
  //! the corresponding equations are: \partial_t A = -(4*1*B)*A  (negative since reactant)
  //!                              \partial_t B = -(4*2*A)*B  (negative since reactant)
  //!                              \partial_t C = + 4*3*A*B   (positive since product)
  //!
  //! this equation is in BACI achieved by the MAT_scatra_reaction material:
  //! ----------------------------------------------------------MATERIALS
  //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 4 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 5 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -2 3 REACOEFF 4.0 COUPLING simple_multiplicative
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
  //! hence we have to calculate and set K(c)=(4*B;8*A;0) and f(c)=(0;0;12*A*B) and corresponding derivatives.


  //! note for the implementation of the michaelis-menten scatra coupling feature:
  //! reactant A promotes C and reactant B influences only until certain limit
  //!
  //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 4 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 5 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -1 1 REACOEFF 5.0 COUPLING michaelis_menten ROLE -1 3 0
  //!
  //! the corresponding equations are
  //!            \partial_t A = -5*A*(B/(3+B))
  //!            \partial_t B = -5*A*(B/(3+B))
  //!            \partial_t C = 5*A*(B/(3+B))
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
  //! hence we have to calculate and set K(c)=(5*(B/(3+B));5*A*(1/(3+B));0) and f(c)=(0;0;5*A*(B/(3+B))) and corresponding derivatives.
  //!
  //! Thereby ROLE does describe of how to build the reaction term (negative value -1.3: simply multiply by A,
  //! positive value 3.2: multiply by B/(B+3.2) ) and STOICH does describe on which scalar the reaction term should be applied.
  //! Here another example:
  //!   //! MAT 3 MAT_scatra_reaction NUMSCAL 4 STOICH +2 -1 0 0  REACOEFF 2.1 COUPLING michaelis_menten ROLE -1 0 1.2 3.0
  //! the corresponding equations are
  //!            \partial_t A = 2*2.1*A*C/(C+1.2)*D/(D+3.0)
  //!            \partial_t B = -1*2.1*A*C/(C+1.2)*D/(D+3.0)
  //!            \partial_t C = 0
  //!            \partial_t D = 0

  //! note for the implementation of the constant scatra coupling feature:
  //! Product A is constantly produced
  //!
  //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 1 MATIDS 2 NUMREAC 1 REACIDS 3 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 1 STOICH 2 REACOEFF 5.0 COUPLING constant
  //!
  //! the corresponding equations are
  //!            \partial_t A = 5*2
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to depend linearly on c_i
  //! hence we have to calculate and set K(c)=(0) and f(c)=(5*2) and zero derivatives.

  //! note for the implementation of the reacstart feature:
  //! Assume concentration A is reproducing with reaction coefficient 1.0 and if the concentration
  //! exceeds some threshold 2.0 if starts to react A->3*B with reacion coefficient 4.0.
  //!
  //! the corresponding equations are:
  //!            \partial_t A = -(-1.0)*A - 4.0*(A - 2.0)_{+} (first termn postive, since equivalent as reactant with negative reaction coefficient)
  //!            \partial_t B = 3.0*4.0 (A - 2.0)_{+}   (positive since product)
  //!
  //! this equation is in BACI achieved by the boundary condition:
  //!   //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 1 MATIDS 2 NUMREAC 2 REACIDS 3 4 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 2 STOICH -1 0 REACOEFF -1.0 COUPLING simple_multiplicative
  //! MAT 4 MAT_scatra_reaction NUMSCAL 2 STOICH -1 3 REACOEFF 4.0 COUPLING simple_multiplicative REACSTART 2.0
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
  //! hence we have to calculate and set K(c)=(-A + 4*(A-2)_{+};0) and f(c)=(0;12*(A-2)_{+}) and corresponding derivatives.

  //! note for the implementation of the power law feature:
  //! assume the following reaction: 1*A + 2*B  --> 3*C with reaction coefficient 4.0
  //!
  //! if we assume the reaction is depending on the product of all
  //! reactants via a power law (this corresponds to couplingtype "power_multiplicative"),
  //! the corresponding equations are: \partial_t A = -(4*1*B^2)*A^3  (negative since reactant)
  //!                              \partial_t B = -(4*2*A^3)*B^2  (negative since reactant)
  //!                              \partial_t C = + 4*3*A^3*B^2   (positive since product)
  //!
  //! this equation is in BACI achieved by the MAT_scatra_reaction material:
  //! ----------------------------------------------------------MATERIALS
  //! MAT 1 MAT_matlist_reactions LOCAL No NUMMAT 3 MATIDS 2 4 5 NUMREAC 1 REACIDS 3 END //collect Concentrations
  //! MAT 2 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 4 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 5 MAT_scatra DIFFUSIVITY 0.0
  //! MAT 3 MAT_scatra_reaction NUMSCAL 3 STOICH -1 -2 3 REACOEFF 4.0 COUPLING power_multiplicative ROLE 3 2 0
  //!
  //! implementation is of form: \partial_t c_i + K_i(c)*c_i = f_i(c), were f_i(c) is supposed not to linear depend on c_i
  //! hence we have to calculate and set K(c)=(4*B;8*A;0) and f(c)=(0;0;12*A*B) and corresponding derivatives.

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype,int probdim>
DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim> * DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::Instance(
  const int numdofpernode,
  const int numscal,
  const std::string& disname,
  const ScaTraEleCalcAdvReac *delete_me )
{
  static std::map<std::pair<std::string,int>,ScaTraEleCalcAdvReac<distype,probdim>* > instances;

  std::pair<std::string,int> key(disname,numdofpernode);

  if(delete_me == NULL)
  {
    if(instances.find(key) == instances.end())
      instances[key] = new ScaTraEleCalcAdvReac<distype,probdim>(numdofpernode,numscal,disname);
  }

  else
  {
    // since we keep several instances around in the general case, we need to
    // find which of the instances to delete with this call. This is done by
    // letting the object to be deleted hand over the 'this' pointer, which is
    // located in the map and deleted
    for( typename std::map<std::pair<std::string,int>,ScaTraEleCalcAdvReac<distype,probdim>* >::iterator i=instances.begin(); i!=instances.end(); ++i )
      if ( i->second == delete_me )
      {
        delete i->second;
        instances.erase(i);
        return NULL;
      }
    dserror("Could not locate the desired instance. Internal error.");
  }

  return instances[key];
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( 0, 0, "", this );
}


/*----------------------------------------------------------------------*
 *  constructor---------------------------                   thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::ScaTraEleCalcAdvReac(const int numdofpernode,const int numscal,const std::string& disname)
  : DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::ScaTraEleCalc(numdofpernode,numscal,disname)
{
  my::reamanager_ = Teuchos::rcp(new ScaTraEleReaManagerAdvReac(my::numscal_));

  // safety check
  if(not my::scatrapara_->TauGP())
    dserror("For advanced reactions, tau needs to be evaluated by integration-point evaluations!");
}

/*----------------------------------------------------------------------*
 |  get the material constants  (private)                      thon 09/14|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::GetMaterialParams(
  const DRT::Element* ele,       //!< the element we are dealing with
  double&             densn,     //!< density at t_(n)
  double&             densnp,    //!< density at t_(n+1) or t_(n+alpha_F)
  double&             densam,    //!< density at t_(n+alpha_M)
  double&             visc,      //!< fluid viscosity
  const int           iquad      //!< id of current gauss point
  )
{
  // get the material
  Teuchos::RCP<MAT::Material> material = ele->Material();

  // We may have some reactive and some non-reactive elements in one discretisation.
  // But since the calculation classes are singleton, we have to reset all reactive stuff in case
  // of non-reactive elements:
  ReaManager()->Clear(my::numscal_);

  if (material->MaterialType() == INPAR::MAT::m_matlist)
  {
    const Teuchos::RCP<const MAT::MatList> actmat = Teuchos::rcp_dynamic_cast<const MAT::MatList>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      Materials(singlemat,k,densn,densnp,densam,visc,iquad);
    }
  }

  else if (material->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    const Teuchos::RCP<MAT::MatListReactions> actmat = Teuchos::rcp_dynamic_cast<MAT::MatListReactions>(material);
    if (actmat->NumMat() != my::numscal_) dserror("Not enough materials in MatList.");

    for (int k = 0;k<my::numscal_;++k)
    {
      int matid = actmat->MatID(k);
      Teuchos::RCP< MAT::Material> singlemat = actmat->MaterialById(matid);

      //Note: order is important here!!
      Materials(singlemat,k,densn,densnp,densam,visc,iquad);

      SetAdvancedReactionTerms(k,actmat); //every reaction calculation stuff happens in here!!
    }
  }

  else
  {
    Materials(material,0,densn,densnp,densam,visc,iquad);
  }

  return;
} //ScaTraEleCalc::GetMaterialParams

/*----------------------------------------------------------------------*
 |  evaluate single material  (protected)                    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::Materials(
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
    my::MatScaTra(material,k,densn,densnp,densam,visc,iquad);
    break;
  case INPAR::MAT::m_biofilm:
    MatBioFilm(material,k,densn,densnp,densam,visc,iquad);
    break;
  case INPAR::MAT::m_scatra_growth_scd:
    MatGrowthScd(material,k,densn,densnp,densam,visc,iquad);
    break;
  default:
    dserror("Material type %i is not supported",material->MaterialType());
    break;
  }
  return;
}


/*----------------------------------------------------------------------*
 |  Material BioFilm                                         thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::MatBioFilm(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,         //!< fluid viscosity
    const int                               iquad         //!< id of current gauss point
  )
{
  const Teuchos::RCP<const MAT::Biofilm>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::Biofilm>(material);

  // get constant diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->Diffusivity(),k);

  // get substrate concentration at n+1 or n+alpha_F at integration point
  const double csnp = my::scatravarmanager_->Phinp(k);
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = ReaManager();

  // set reaction coefficient
  remanager->SetReaCoeff(actmat->ComputeReactionCoeff(csnp),k);

  // set derivative of reaction coefficient
  remanager->SetReaCoeffDerivMatrix(actmat->ComputeReactionCoeffDeriv(csnp),k,k);

  // set density at various time steps and density gradient factor to 1.0/0.0
  densn      = 1.0;
  densnp     = 1.0;
  densam     = 1.0;

  return;
}

/*----------------------------------------------------------------------*
 |  Material GrowthScd                                       vuong 01/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::MatGrowthScd(
    const Teuchos::RCP<const MAT::Material> material, //!< pointer to current material
    const int                               k,        //!< id of current scalar
    double&                                 densn,    //!< density at t_(n)
    double&                                 densnp,   //!< density at t_(n+1) or t_(n+alpha_F)
    double&                                 densam,   //!< density at t_(n+alpha_M)
    double&                                 visc,         //!< fluid viscosity
    const int                               iquad         //!< id of current gauss point
  )
{
  dsassert(my::numdofpernode_==1,"more than 1 dof per node for ScatraGrowthScd material");

  if(iquad < 0) dserror("ScatraGrowthScd material has to be evaluated at gauss point!");

  const Teuchos::RCP<const MAT::ScatraGrowthScd>& actmat
    = Teuchos::rcp_dynamic_cast<const MAT::ScatraGrowthScd>(material);

  // get and save constant diffusivity
  my::diffmanager_->SetIsotropicDiff(actmat->Diffusivity(),k);

  //strategy to obtain theta from the structure at equivalent gauss-point
  //access structure discretization
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  //get corresponding structure element (it has the same global ID as the scatra element)
  DRT::Element* structele = structdis->gElement(my::eid_);
  if (structele == NULL)
    dserror("Structure element %i not on local processor", my::eid_);

  const Teuchos::RCP<const MAT::GrowthScd>& structmat
          = Teuchos::rcp_dynamic_cast<const MAT::GrowthScd>(structele->Material());
  if (structmat == Teuchos::null)
    dserror("dynamic cast of structure material GrowthScd failed.");
  if(structmat->MaterialType() != INPAR::MAT::m_growth_volumetric_scd)
    dserror("invalid structure material for scalar dependent growth");

  if (structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_linear or
      structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_exponential)
  {
    const double theta    = structmat->Gettheta_atgp(iquad);
    const double dtheta   = structmat->Getdtheta_atgp(iquad);
    const double thetaold = structmat->Getthetaold_atgp(iquad);
    const double detFe    = structmat->GetdetFe_atgp(iquad);

    // get substrate concentration at n+1 or n+alpha_F at integration point
    const double csnp = my::scatravarmanager_->Phinp(k);
    const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = ReaManager();

    // set reaction coefficient
    remanager->SetReaCoeff(actmat->ComputeReactionCoeff(csnp,theta,dtheta,detFe),k);
    // set derivative of reaction coefficient
    remanager->SetReaCoeffDerivMatrix(actmat->ComputeReactionCoeffDeriv(csnp,theta,thetaold,1.0),k,k);

    // set density at various time steps and density gradient factor to 1.0/0.0
    densn      = 1.0;
    densnp     = 1.0;
    densam     = 1.0;
  }
  else if (structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_ac or
           structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_ac_radial or
           structmat->Parameter()->growthlaw_->MaterialType() == INPAR::MAT::m_growth_ac_radial_refconc )

    {
    dserror("In the case of MAT_GrowthAC or MAT_GrowthACNormal one should not end up in here, "
        "since the growth does only change the scalars field size/volume. And this is already"
        " cared due to the conservative formulation you hopefully use!");
    }
  else
  {
    dserror("Your growth law is not a valid one!");
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Get right hand side including reaction bodyforce term    thon 02/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::GetRhsInt(
  double&      rhsint,   //!< rhs containing bodyforce at Gauss point
  const double densnp,  //!< density at t_(n+1)
  const int    k        //!< index of current scalar
  )
{
                                       //... + reaction terms not depending on phi(k) -> source term
  rhsint = my::bodyforce_[k].Dot(my::funct_) + densnp*ReaManager()->GetReaBodyForce(k);

  return;
}

/*--------------------------------------------------------------------------- *
 |  calculation of reactive element matrix for coupled reactions  thon 02/14  |
 *----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::CalcMatReact(
  Epetra_SerialDenseMatrix&          emat,
  const int                          k,
  const double                       timefacfac,
  const double                       timetaufac,
  const double                       taufac,
  const double                       densnp,
  const LINALG::Matrix<my::nen_,1>&      sgconv,
  const LINALG::Matrix<my::nen_,1>&      diff
  )
{
  // -----------------first care for Term K(c)*(\partial_c c)=Id*K(c)--------------------------------------

  my::CalcMatReact(emat,k,timefacfac,timetaufac,taufac,densnp,sgconv,diff);

  const double&                       phinp = my::scatravarmanager_->Phinp(k);
  const LINALG::Matrix<my::nen_,1>&   conv  = my::scatravarmanager_->Conv();

  // -----------------second care for Term (\partial_c K(c)) .* c - (\partial_c f_{reabody}(c))------------

  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = ReaManager();

  LINALG::Matrix<my::nen_,1> functint = my::funct_;
  if (not my::scatrapara_->MatGP())
    functint = funct_elementcenter_;

  for (int j=0; j<my::numscal_ ;j++)
  {
    const double fac_reac        = timefacfac*densnp*( remanager->GetReaCoeffDerivMatrix(k,j)*phinp - remanager->GetReaBodyForceDerivMatrix(k,j) );
    const double timetaufac_reac = timetaufac*densnp*( remanager->GetReaCoeffDerivMatrix(k,j)*phinp - remanager->GetReaBodyForceDerivMatrix(k,j) );

    //----------------------------------------------------------------
    // standard Galerkin reactive term
    //----------------------------------------------------------------
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac_reac*functint(vi);
      const int fvi = vi*my::numdofpernode_+k;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+j;

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
        const double v = densreataufac*(conv(vi)+sgconv(vi)+my::scatrapara_->USFEMGLSFac()*1.0/my::scatraparatimint_->TimeFac()*functint(vi));
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+j;

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
            const int fui = ui*my::numdofpernode_+j;

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
        const double v = densreataufac*functint(vi);
        const int fvi = vi*my::numdofpernode_+k;

        for (int ui=0; ui<my::nen_; ++ui)
        {
          const int fui = ui*my::numdofpernode_+j;

          emat(fvi,fui) += v*(conv(ui)+remanager->GetReaCoeff(k)*my::funct_(ui));
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
            const int fui = ui*my::numdofpernode_+j;

            emat(fvi,fui) -= v*diff(ui);
          }
        }
      }


      if (not my::scatraparatimint_->IsStationary())
      {
        // reactive stabilization of transient term
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = my::scatrapara_->USFEMGLSFac()*taufac*densnp*remanager->GetReaCoeff(k)*densnp*functint(vi);
          const int fvi = vi*my::numdofpernode_+k;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+j;

            emat(fvi,fui) += v*my::funct_(ui);
          }
        }

        if (my::use2ndderiv_ and remanager->GetReaCoeff(k)!=0.0)
          dserror("Second order reactive stabilization is not fully implemented!! ");
      }
    }
  } //end for
  return;
}


/*-------------------------------------------------------------------------------*
 |  set reac. body force, reaction coefficient and derivatives        thon 09/14 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::SetAdvancedReactionTerms(
    const int                                 k,          //!< index of current scalar
    const Teuchos::RCP<MAT::MatListReactions> matreaclist //!< index of current scalar
    )
{
  const Teuchos::RCP<ScaTraEleReaManagerAdvReac> remanager = ReaManager();

  remanager->AddToReaBodyForce( matreaclist->CalcReaBodyForceTerm(k,my::scatravarmanager_->Phinp()) ,k );

  remanager->SetReaCoeff( matreaclist->CalcReaCoeff(k,my::scatravarmanager_->Phinp()) ,k);

  for (int j=0; j<my::numscal_ ;j++)
  {
    remanager->AddToReaBodyForceDerivMatrix( matreaclist->CalcReaBodyForceDerivMatrix(k,j,my::scatravarmanager_->Phinp()) ,k,j );

    remanager->SetReaCoeffDerivMatrix( matreaclist->CalcReaCoeffDerivMatrix(k,j,my::scatravarmanager_->Phinp()) ,k,j );
  }

}

/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at ele. center   jhoer 11/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
double DRT::ELEMENTS::ScaTraEleCalcAdvReac<distype,probdim>::EvalShapeFuncAndDerivsAtEleCenter()
{
  const double vol = my::EvalShapeFuncAndDerivsAtEleCenter();

  //shape function at element center
  funct_elementcenter_ = my::funct_;

  return vol;

} //ScaTraImpl::EvalShapeFuncAndDerivsAtEleCenter


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line2,1>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line2,2>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line2,3>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::line3,1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri3,2>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri3,3>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tri6,2>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad4,2>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad4,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::quad9,2>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs9,2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex8,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::hex27,3>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tet4,3>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::tet10,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::pyramid5,3>;
//template class DRT::ELEMENTS::ScaTraEleCalcAdvReac<DRT::Element::nurbs27>;
