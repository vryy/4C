/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch_diffcond.cpp

\brief evaluation of scatra elements for elch

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15252
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "scatra_ele.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/elchphase.H"

#include "scatra_ele_calc_elch_diffcond.H"


/*-----------------------------------------------------------------------*
  |  Set scatra element parameter                             ehrl 01/14 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CheckElchElementParameter(
    DRT::Element*              ele
  )
{
  // get the material
  Teuchos::RCP<const MAT::Material> material = ele->Material();

  // 1) Check material specific options
  // 2) Check if numdofpernode, numscal is set correctly
  if (material->MaterialType() == INPAR::MAT::m_elchmat)
  {
    const Teuchos::RCP<const MAT::ElchMat>& actmat
          = Teuchos::rcp_dynamic_cast<const MAT::ElchMat>(material);

    int numphase = actmat->NumPhase();

    // access mat_elchmat: container material for porous structures in elch
    if (numphase != 1) dserror("In the moment a single phase is only allowed.");

    // 1) loop over single phases
    for (int iphase=0; iphase < actmat->NumPhase();++iphase)
    {
      // access phase material
      const int phaseid = actmat->PhaseID(iphase);
      Teuchos::RCP<const MAT::Material> singlephase = actmat->PhaseById(phaseid);

      // dynmic cast: get access to mat_phase
      const Teuchos::RCP<const MAT::ElchPhase>& actphase
                = Teuchos::rcp_dynamic_cast<const MAT::ElchPhase>(singlephase);

      // Check if numdofpernode, numscal is set correctly
      int nummat = actphase->NumMat();
      // enough materials defined
      if (nummat != my::numscal_)
        dserror("The number of scalars defined in the material ElchMat does not correspond with "
                "the number of materials defined in the material MatPhase.");

      int numdofpernode = 0;
      if (ElchPara()->CurSolVar())
        numdofpernode = nummat+DRT::Problem::Instance()->NDim()+numphase;
      else
        numdofpernode = nummat+numphase;

      if(numdofpernode != my::numdofpernode_)
        dserror("The chosen element formulation (e.g. current as solution variable) "
                "does not correspond with the number of dof's defined in your material");

      // 2) loop over materials of the single phase
      for (int imat=0; imat < actphase->NumMat();++imat)
      {
        const int matid = actphase->MatID(imat);
        Teuchos::RCP<const MAT::Material> singlemat = actphase->MatById(matid);

        if(singlemat->MaterialType() == INPAR::MAT::m_newman)
        {
          // Material Newman is derived for a binary electrolyte utilizing the ENC to condense the non-reacting species
          if(my::numscal_>1)
            dserror("Material Newman is only valid for one scalar (binary electrolyte utilizing the ENC)");
        }
      }
    }
  }

  return;
}


/*---------------------------------------------------------------------------------------------*
 | calculate element mass matrix and element residual for initial time derivative   fang 03/15 |
 *---------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalcInitialTimeDerivative(
    DRT::Element*               ele,              //!< current element
    Epetra_SerialDenseMatrix&   emat,             //!< element matrix
    Epetra_SerialDenseVector&   erhs,             //!< element residual
    Teuchos::ParameterList&     params,           //!< parameter list
    DRT::Discretization&        discretization,   //!< discretization
    const std::vector<int>&     lm                //!< location vector
    )
{
  // call base class routine
  myelch::CalcInitialTimeDerivative(
      ele,
      emat,
      erhs,
      params,
      discretization,
      lm
      );

  // dummy mass matrix and zero residual entries for the electric current dofs
  if(ElchPara()->CurSolVar())
  {
    // integration points and weights
    const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      for(int idim=0;idim<my::nsd_;++idim)
      {
        for (int vi=0; vi<my::nen_; ++vi)
        {
          const double v = fac*my::funct_(vi); // no density required here
          const int fvi = vi*my::numdofpernode_+my::numscal_+1+idim;

          for (int ui=0; ui<my::nen_; ++ui)
          {
            const int fui = ui*my::numdofpernode_+my::numscal_+1+idim;

            emat(fvi,fui) += v*my::funct_(ui);
          }

          erhs[vi*my::numdofpernode_+my::numscal_+1+idim]=0.0;
        }
      }
    }
  }

  // In the moment the diffusion manager contains the porosity at the last Gauss point (previous call my::CalcInitialTimeDerivative())
  // Since the whole approach is valid only for constant porosities
  // we do not fill the diffusion manager again at the element center

  // The solution variable is the initial time derivative
  // Therefore, we have to correct rhs by the initial porosity
  // attention: this procedure is only valid for a constant porosity in the beginning
  erhs.Scale(1.0/DiffManager()->GetPhasePoro(0));

  return;
}


/*----------------------------------------------------------------------*
 |  CorrectRHSFromCalcRHSLinMass                             ehrl 06/14 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CorrectRHSFromCalcRHSLinMass(
    Epetra_SerialDenseVector&     erhs,
    const int                     k,
    const double                  fac,
    const double                  densnp,
    const double                  phinp
  )
{
  // fac->-fac to change sign of rhs
  if (my::scatraparatimint_->IsIncremental())
     my::CalcRHSLinMass(erhs,k,0.0,-fac,0.0,DiffManager()->GetPhasePoro(0));
  else
    dserror("Must be incremental!");

  return;
}


/*-----------------------------------------------------------------------------------------*
 | extract element based or nodal values and return extracted values of phinp   fang 01/15 |
 *-----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
const std::vector<double> DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::ExtractElementAndNodeValues(
    DRT::Element*              ele,
    Teuchos::ParameterList&    params,
    DRT::Discretization&       discretization,
    const std::vector<int>&    lm
    )
{
  // call base class routine
  const std::vector<double> myphinp = myelch::ExtractElementAndNodeValues(ele,params,discretization,lm);

  // get current density at element nodes if required
  if(ElchPara()->CurSolVar())
  {
    // get values for current at element nodes
    for (int ien=0;ien<my::nen_;++ien)
    {
      for(int idim=0; idim<my::nsd_; ++idim)
      {
        //current is stored after potential
        ecurnp_(idim,ien) = myphinp[ien*my::numdofpernode_+(my::numscal_+1)+idim];
      }
    }
  }
  else
    ecurnp_.Clear();

  return myphinp;
}


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     ae 05/15|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalculateFlux(
    LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const int                       k,          //!< index of current scalar
    const double                    fac         //!< integration factor
)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
  case INPAR::SCATRA::flux_total_domain:
    // call base class routine
    myelectrode::CalculateFlux(q,fluxtype,k,fac);

    // no break statement here!
  case INPAR::SCATRA::flux_diffusive_domain:
    // diffusive flux contribution
    q.Update(-DiffManager()->GetIsotropicDiff(k),VarManager()->GradPhi(k),1.0);
    // flux due to ohmic overpotential
    q.Update(-DiffManager()->GetTransNum(k)*VarManager()->InvFVal(k)*DiffManager()->GetCond(),VarManager()->GradPot(),1.0);
    // flux due to concentration overpotential
    q.Update(-DiffManager()->GetTransNum(k)*VarManager()->RTFFCVal(k)*DiffManager()->GetCond()*DiffManager()->GetThermFac()*(ElchPara()->NewmanConstA()+(ElchPara()->NewmanConstB()*DiffManager()->GetTransNum(k)))*VarManager()->ConIntInv(k),VarManager()->GradPhi(k),1.0);
    break;
  default:
    dserror("received illegal flag inside flux evaluation for whole domain"); break;
  };

  return;
} // ScaTraCalc::CalculateFlux


/*----------------------------------------------------------------------*
  |  calculate weighted mass flux (no reactive flux so far)     ae 05/15|
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::CalculateCurrent(
    LINALG::Matrix<my::nsd_,1>&     q,          //!< flux of species k
    const INPAR::SCATRA::FluxType   fluxtype,   //!< type fo flux
    const double                    fac         //!< integration factor
)
{
  /*
  * Actually, we compute here a weighted (and integrated) form of the fluxes!
  * On time integration level, these contributions are then used to calculate
  * an L2-projected representation of fluxes.
  * Thus, this method here DOES NOT YET provide flux values that are ready to use!!
  /                                                         \
  |                /   \                               /   \  |
  | w, -D * nabla | phi | + u*phi - frt*z_k*c_k*nabla | pot | |
  |                \   /                               \   /  |
  \                      [optional]      [ELCH]               /
  */

  // add different flux contributions as specified by user input
  switch (fluxtype)
  {
  case INPAR::SCATRA::flux_total_domain:
  case INPAR::SCATRA::flux_diffusive_domain:
    // ohmic flux contribution
    q.Update(-DiffManager()->GetCond(),VarManager()->GradPot(),1.0);
    // diffusion overpotential flux contribution
    for (int k = 0; k<my::numscal_; ++k)
      q.Update(-VarManager()->RTF()/ElchPara()->NewmanConstC()*DiffManager()->GetCond()*DiffManager()->GetThermFac()*(ElchPara()->NewmanConstA()+(ElchPara()->NewmanConstB()*DiffManager()->GetTransNum(k)))*VarManager()->ConIntInv(k),VarManager()->GradPhi(k),1.0);

    break;
  default:
    dserror("received illegal flag inside flux evaluation for whole domain"); break;
  };


  return;
} // ScaTraCalc::CalculateCurrent


/*------------------------------------------------------------------------------*
 | set internal variables for diffusion-conduction formulation       fang 02/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<distype>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchDiffCond(my::funct_,my::derxy_,my::ephinp_,my::ephin_,myelch::epotnp_,my::econvelnp_,my::ehist_,DiffManager(),ecurnp_);

  return;
}


// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElchDiffCond<DRT::Element::nurbs27>;
