/*--------------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_elch.cpp

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
#include "scatra_ele_action.H"
#include "scatra_ele_parameter_elch.H"

#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_discret.H"  // for time curve in body force
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_globalproblem.H"  // consistency check of formulation and material

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_nurbs_discret/drt_nurbs_utils.H"

#include "../drt_mat/elchmat.H"
#include "../drt_mat/newman.H"
#include "../drt_mat/elchphase.H"

#include "scatra_ele_calc_elch.H"


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleCalcElch<distype>::EvaluateAction(
    DRT::Element*                 ele,
    Teuchos::ParameterList&       params,
    DRT::Discretization&          discretization,
    const SCATRA::Action&         action,
    DRT::Element::LocationArray&  la,
    Epetra_SerialDenseMatrix&     elemat1_epetra,
    Epetra_SerialDenseMatrix&     elemat2_epetra,
    Epetra_SerialDenseVector&     elevec1_epetra,
    Epetra_SerialDenseVector&     elevec2_epetra,
    Epetra_SerialDenseVector&     elevec3_epetra
    )
{
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::check_scatra_element_parameter:
  {
    CheckElchElementParameter(ele);
    break;
  }

  case SCATRA::integrate_shape_functions:
  {
    // calculate integral of shape functions
    const Epetra_IntSerialDenseVector dofids = params.get<Epetra_IntSerialDenseVector>("dofids");
    my::IntegrateShapeFunctions(ele,elevec1_epetra,dofids);

    break;
  }

  case SCATRA::calc_flux_domain:
  {
    //--------------------------------------------------------------------------------
    // extract element based or nodal values
    //--------------------------------------------------------------------------------

    // get velocity values at the nodes
    const Teuchos::RCP<Epetra_MultiVector> velocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,my::evelnp_,velocity,my::nsd_);
    const Teuchos::RCP<Epetra_MultiVector> convelocity = params.get< Teuchos::RCP<Epetra_MultiVector> >("convective velocity field");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,my::econvelnp_,convelocity,my::nsd_);

    // need current values of transported scalar
    // -> extract local values from global vectors
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");
    std::vector<double> myphinp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,la[0].lm_);

    // fill all element arrays
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0;k<my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        my::ephinp_[k](i,0) = myphinp[k+(i*my::numdofpernode_)];
      }
      epotnp_(i) = myphinp[i*my::numdofpernode_+my::numscal_];
    } // for i

    //----------------------------------------------------------------------
    // calculation of element volume both for tau at ele. cent. and int. pt.
    //----------------------------------------------------------------------
    my::EvalShapeFuncAndDerivsAtEleCenter();

    //----------------------------------------------------------------------
    // get material and stabilization parameters (evaluation at element center)
    //----------------------------------------------------------------------
    // density at t_(n)
    double densn(1.0);
    // density at t_(n+1) or t_(n+alpha_F)
    double densnp(1.0);
    // density at t_(n+alpha_M)
    double densam(1.0);

    // fluid viscosity
    double visc(0.0);

    // material parameter at the element center
    if (not my::scatrapara_->MatGP())
      this->GetMaterialParams(ele,densn,densnp,densam,visc);

    //----------------------------------------------------------------------
    // integration loop for one element
    //----------------------------------------------------------------------
    // integration points and weights
    const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

    for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
    {
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

      //----------------------------------------------------------------------
      // get material parameters (evaluation at integration point)
      //----------------------------------------------------------------------
      if (my::scatrapara_->MatGP())
        this->GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

      // set internal variables
      SetInternalVariablesForMatAndRHS();

      // access control parameter for flux calculation
      INPAR::SCATRA::FluxType fluxtype = ElchPara()->WriteFlux();
      Teuchos::RCP<std::vector<int> > writefluxids = ElchPara()->WriteFluxIds();

      // do a loop for systems of transported scalars
      for (std::vector<int>::iterator it = writefluxids->begin(); it!=writefluxids->end(); ++it)
      {
        int k=0;
        // Actually, we compute here a weighted (and integrated) form of the fluxes!
        // On time integration level, these contributions are then used to calculate
        // an L2-projected representation of fluxes.
        // Thus, this method here DOES NOT YET provide flux values that are ready to use!!

        // allocate and initialize!
        LINALG::Matrix<my::nsd_,1> q(true);

        if((*it) != my::numdofpernode_)
        {
          k=(*it)-1;
          CalculateFlux(q,fluxtype,k,fac);
        }
        else if ((*it) == my::numdofpernode_)
        {
          k=my::numdofpernode_-1;
          CalculateCurrent(q,fluxtype,fac);
        }
        else
          dserror("Flux id, which should be calculated, does not exit in the dof set.");

        // integrate and assemble everything into the "flux" vector
        for (int vi=0; vi < my::nen_; vi++)
        {
          const int fvi = vi*my::numdofpernode_+k;
          elevec1_epetra[fvi] += fac*my::funct_(vi)*q(0);
          elevec2_epetra[fvi] += fac*my::funct_(vi)*q(1);
          if(my::nsd_<3)
            elevec3_epetra[fvi] = 0.0;
          else
            elevec3_epetra[fvi] += fac*my::funct_(vi)*q(2);
        } // vi
      }
    }

    break;
  }
  case SCATRA::calc_error:
  {
    // check if length suffices
    if (elevec1_epetra.Length() < 1) dserror("Result vector too short");

    // need current solution
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> myphinp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,la[0].lm_);

    // fill element arrays
    for (int i=0;i<my::nen_;++i)
    {
      // split for each transported scalar, insert into element arrays
      for (int k = 0; k< my::numscal_; ++k)
      {
        my::ephinp_[k](i) = myphinp[k+(i*my::numdofpernode_)];
      }
      // get values for el. potential at element nodes
      epotnp_(i) = myphinp[i*my::numdofpernode_+my::numscal_];
    } // for i

    CalErrorComparedToAnalytSolution(
      ele,
      params,
      elevec1_epetra);

    break;
  }

  case SCATRA::calc_elch_conductivity:
  {
    // get flag if effective conductivity should be calculated
    bool effCond = params.get<bool>("effCond");
    // get flag if the inverse of the conductivity should be calculated -> specific resistance
    bool specresist = params.get<bool>("specresist");

    // extract local values from the global vector
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    std::vector<double> myphinp(la[0].lm_.size());
    DRT::UTILS::ExtractMyValues(*phinp,myphinp,la[0].lm_);

    // fill element arrays
    for (int i=0;i<my::nen_;++i)
    {
      for (int k = 0; k< my::numscal_; ++k)
      {
        // split for each transported scalar, insert into element arrays
        my::ephinp_[k](i,0) = myphinp[k+(i*my::numdofpernode_)];
      }
    } // for i

    // elevec1_epetra(0):          conductivity of ionic species 0
    // elevec1_epetra(numscal_-1): conductivity of ionic species (numscal_-1)
    // elevec1_epetra(numscal_):   conductivity of the electrolyte solution (sum_k sigma(k))
    // elevec1_epetra(numscal_+1): domain integral
    CalculateConductivity(ele,ElchPara()->EquPot(), elevec1_epetra, effCond, specresist);
    break;
  }

  default:
  {
    my::EvaluateAction(ele,
                       params,
                       discretization,
                       action,
                       la,
                       elemat1_epetra,
                       elemat2_epetra,
                       elevec1_epetra,
                       elevec2_epetra,
                       elevec3_epetra);
    break;
  }
  } // switch(action)

  return 0;
}


/*----------------------------------------------------------------------*
  |  Calculate conductivity (ELCH) (private)                   gjb 07/09 |
  *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalculateConductivity(
  const DRT::Element*               ele,
  const enum INPAR::ELCH::EquPot    equpot,
  Epetra_SerialDenseVector&         sigma_domint,
  bool                              effCond,
  bool                              specresist
  )
{
  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    //----------------------------------------------------------------------
    // get material and stabilization parameters (evaluation at element center)
    //----------------------------------------------------------------------
    // density at t_(n)
    double densn(1.0);
    // density at t_(n+1) or t_(n+alpha_F)
    double densnp(1.0);
    // density at t_(n+alpha_M)
    double densam(1.0);

    // fluid viscosity
    double visc(0.0);

    // material parameter at the element center
    this->GetMaterialParams(ele,densn,densnp,densam,visc,iquad);

    // calculate integrals of (inverted) scalar(s) and domain
    for (int i=0; i<my::nen_; i++)
    {
      double sigma_all(0.0);
      std::vector<double> sigma(my::numscal_, 0.0);
      // compute the conductivity (1/(\Omega m) = 1 Siemens / m)
      GetConductivity(equpot, sigma_all, sigma, effCond);

      const double fac_funct_i = fac*my::funct_(i);

      // sigma_domint(0):          conductivity of ionic species 0
      // sigma_domint(numscal_-1): conductivity of ionic species (numscal_-1)
      // sigma_domint(numscal_):   conductivity of the electrolyte solution (sum_k sigma(k))
      // sigma_domint(numscal_+1): domain integral
      for (int k = 0; k < my::numscal_; k++)
      {
        sigma_domint[k]        += sigma[k]*fac_funct_i;
      }

      // calculation of conductivity or specific resistance of electrolyte solution
      if(specresist == false)
        sigma_domint[my::numscal_]   += sigma_all*fac_funct_i;
      else
        sigma_domint[my::numscal_]   += 1/sigma_all*fac_funct_i;

      // domain integral
      sigma_domint[my::numscal_+1] += fac_funct_i;
    }
  } // loop over integration points

  return;

} //ScaTraEleCalcElch()


/*---------------------------------------------------------------------------------------------*
 | calculate element mass matrix and element residual for initial time derivative   fang 03/15 |
 *---------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::CalcInitialTimeDerivative(
    DRT::Element*                 ele,              //!< current element
    Epetra_SerialDenseMatrix&     emat,             //!< element matrix
    Epetra_SerialDenseVector&     erhs,             //!< element residual
    Teuchos::ParameterList&       params,           //!< parameter list
    DRT::Discretization&          discretization,   //!< discretization
    DRT::Element::LocationArray&  la                //!< location array
    )
{
  // call base class routine
  my::CalcInitialTimeDerivative(
      ele,
      emat,
      erhs,
      params,
      discretization,
      la
      );

  // do another loop over the integration points for the potential:
  // At this point time is not so critical since the CalcInitialTimeDerivative()
  // is called once in the beginning of the simulation!

  // we put a dummy mass matrix for the electric potential dofs
  // here in order to have a regular matrix in the lower right block of the whole system-matrix
  // An identity matrix would cause problems with ML solver in the SIMPLE
  // schemes since ML needs to have off-diagonal entries for the aggregation!

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // element integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad);

    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double v = fac*my::funct_(vi); // no density required here
      const int fvi = vi*my::numdofpernode_+my::numscal_;

      for (int ui=0; ui<my::nen_; ++ui)
      {
        const int fui = ui*my::numdofpernode_+my::numscal_;

        emat(fvi,fui) += v*my::funct_(ui);
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------------------*
 | finite difference check on element level (for debugging only) (protected)   fang 10/14 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleCalcElch<distype>::FDCheck(
  DRT::Element*                              ele,
  Epetra_SerialDenseMatrix&                  emat,
  Epetra_SerialDenseVector&                  erhs,
  Epetra_SerialDenseVector&                  subgrdiff
  )
{
  // screen output
  std::cout << "FINITE DIFFERENCE CHECK FOR ELEMENT " << ele->Id();

  // make a copy of state variables to undo perturbations later
  std::vector<LINALG::Matrix<my::nen_,1> > ephinp_original(my::numscal_);
  for(int k=0; k<my::numscal_; ++k)
    for(int i=0; i<my::nen_; ++i)
      ephinp_original[k](i,0) = my::ephinp_[k](i,0);
  LINALG::Matrix<my::nen_,1> epotnp_original(true);
  for(int i=0; i<my::nen_; ++i)
    epotnp_original(i) = epotnp_(i);

  // generalized-alpha time integration requires a copy of history variables as well
  std::vector<LINALG::Matrix<my::nen_,1> > ehist_original(my::numscal_);
  if(my::scatraparatimint_->IsGenAlpha())
  {
    for(int k=0; k<my::numscal_; ++k)
      for(int i=0; i<my::nen_; ++i)
        ehist_original[k](i,0)  = my::ehist_[k](i,0);
  }

  // initialize element matrix and vectors for perturbed state
  Epetra_SerialDenseMatrix emat_dummy(emat);
  Epetra_SerialDenseVector erhs_perturbed(erhs);
  Epetra_SerialDenseVector subgrdiff_dummy(subgrdiff);

  // initialize counter for failed finite difference checks
  unsigned counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  // loop over columns of element matrix by first looping over nodes and then over dofs at each node
  for(int inode=0; inode<my::nen_; ++inode)
  {
    for(int idof=0; idof<my::numdofpernode_; ++idof)
    {
      // number of current column of element matrix
      unsigned col = inode*my::numdofpernode_+idof;

      // clear element matrix and vectors for perturbed state
      emat_dummy.Scale(0.0);
      erhs_perturbed.Scale(0.0);
      subgrdiff_dummy.Scale(0.0);

      // fill state vectors with original state variables
      for(int k=0; k<my::numscal_; ++k)
        for(int i=0; i<my::nen_; ++i)
          my::ephinp_[k](i,0) = ephinp_original[k](i,0);
      for(int i=0; i<my::nen_; ++i)
        epotnp_(i) = epotnp_original(i);
      if(my::scatraparatimint_->IsGenAlpha())
        for(int k=0; k<my::numscal_; ++k)
          for(int i=0; i<my::nen_; ++i)
            my::ehist_[k](i,0)  = ehist_original[k](i,0);

      // impose perturbation
      if(my::scatraparatimint_->IsGenAlpha())
      {
        // perturbation on concentration
        if(idof < my::numscal_)
        {
          // perturbation of phi(n+alphaF), not of phi(n+1) => scale epsilon by factor alphaF
          my::ephinp_[idof](inode,0) += my::scatraparatimint_->AlphaF() * my::scatrapara_->FDCheckEps();

          // perturbation of phi(n+alphaF) by alphaF*epsilon corresponds to perturbation of phidtam (stored in ehist_)
          // by alphaM*epsilon/(gamma*dt); note: alphaF/timefac = alphaM/(gamma*dt)
          my::ehist_[idof](inode,0) += my::scatraparatimint_->AlphaF() / my::scatraparatimint_->TimeFac() * my::scatrapara_->FDCheckEps();
        }

        // perturbation on electric potential
        else
          epotnp_(inode) += my::scatraparatimint_->AlphaF() * my::scatrapara_->FDCheckEps();
      }
      else
      {
        // perturbation on concentration
        if(idof < my::numscal_)
          my::ephinp_[idof](inode,0) += my::scatrapara_->FDCheckEps();
        else
          epotnp_(inode) += my::scatrapara_->FDCheckEps();
      }

      // calculate element right-hand side vector for perturbed state
      Sysmat(ele,emat_dummy,erhs_perturbed,subgrdiff_dummy);

      // Now we compare the difference between the current entries in the element matrix
      // and their finite difference approximations according to
      // entries ?= (-erhs_perturbed + erhs_original) / epsilon

      // Note that the element right-hand side equals the negative element residual.
      // To account for errors due to numerical cancellation, we additionally consider
      // entries - erhs_original / epsilon ?= -erhs_perturbed / epsilon

      // Note that we still need to evaluate the first comparison as well. For small entries in the element
      // matrix, the second comparison might yield good agreement in spite of the entries being wrong!
      for(int row=0; row<my::numdofpernode_*my::nen_; ++row)
      {
        // get current entry in original element matrix
        const double entry = emat(row,col);

        // finite difference suggestion (first divide by epsilon and then subtract for better conditioning)
        const double fdval = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps() + erhs(row) / my::scatrapara_->FDCheckEps();

        // confirm accuracy of first comparison
        if(abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
          dserror("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in first comparison
        const double abserr1 = entry - fdval;
        if(abs(abserr1) > abs(maxabserr))
          maxabserr = abserr1;
        double relerr1(0.);
        if(abs(entry) > 1.e-17)
          relerr1 = abserr1 / abs(entry);
        else if(abs(fdval) > 1.e-17)
          relerr1 = abserr1 / abs(fdval);
        if(abs(relerr1) > abs(maxrelerr))
          maxrelerr = relerr1;

        // evaluate first comparison
        if(abs(relerr1) > my::scatrapara_->FDCheckTol())
        {
          if(!counter)
            std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
          std::cout << "emat[" << row << "," << col << "]:  " << entry << "   ";
          std::cout << "finite difference suggestion:  " << fdval << "   ";
          std::cout << "absolute error:  " << abserr1 << "   ";
          std::cout << "relative error:  " << relerr1 << std::endl;

          counter++;
        }

        // first comparison OK
        else
        {
          // left-hand side in second comparison
          const double left  = entry - erhs(row) / my::scatrapara_->FDCheckEps();

          // right-hand side in second comparison
          const double right = -erhs_perturbed(row) / my::scatrapara_->FDCheckEps();

          // confirm accuracy of second comparison
          if(abs(right) > 1.e-17 and abs(right) < 1.e-15)
            dserror("Finite difference check involves values too close to numerical zero!");

          // absolute and relative errors in second comparison
          const double abserr2 = left - right;
          if(abs(abserr2) > abs(maxabserr))
            maxabserr = abserr2;
          double relerr2(0.);
          if(abs(left) > 1.e-17)
            relerr2 = abserr2 / abs(left);
          else if(abs(right) > 1.e-17)
            relerr2 = abserr2 / abs(right);
          if(abs(relerr2) > abs(maxrelerr))
            maxrelerr = relerr2;

          // evaluate second comparison
          if(abs(relerr2) > my::scatrapara_->FDCheckTol())
          {
            if(!counter)
              std::cout << " --> FAILED AS FOLLOWS:" << std::endl;
            std::cout << "emat[" << row << "," << col << "]-erhs[" << row << "]/eps:  " << left << "   ";
            std::cout << "-erhs_perturbed[" << row << "]/eps:  " << right << "   ";
            std::cout << "absolute error:  " << abserr2 << "   ";
            std::cout << "relative error:  " << relerr2 << std::endl;

            counter++;
          }
        }
      }
    }
  }

  // screen output in case finite difference check is passed
  if(!counter)
    std::cout << " --> PASSED WITH MAXIMUM ABSOLUTE ERROR " << maxabserr << " AND MAXIMUM RELATIVE ERROR " << maxrelerr << std::endl;

  // undo perturbations of state variables
  for(int k=0; k<my::numscal_; ++k)
    for(int i=0; i<my::nen_; ++i)
      my::ephinp_[k](i,0) = ephinp_original[k](i,0);
  for(int i=0; i<my::nen_; ++i)
    epotnp_(i) = epotnp_original(i);
  if(my::scatraparatimint_->IsGenAlpha())
    for(int k=0; k<my::numscal_; ++k)
      for(int i=0; i<my::nen_; ++i)
        my::ehist_[k](i,0)  = ehist_original[k](i,0);

  return;
}


// template classes

// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::line3>;

// 2D elements
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri3>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad4>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs9>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex8>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::hex27>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet4>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::tet10>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::pyramid5>;
//template class DRT::ELEMENTS::ScaTraEleCalcElch<DRT::Element::nurbs27>;
