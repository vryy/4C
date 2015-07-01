/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_boundary_calc_elch.cpp

\brief evaluation of scatra boundary terms at integration points

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_ele_boundary_calc_elch.H"
#include "scatra_ele_parameter_elch.H"
#include "scatra_ele_utils_elch.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 01/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::ScaTraEleBoundaryCalcElch(const int numdofpernode, const int numscal) :
  // constructor of base class
  my::ScaTraEleBoundaryCalc(numdofpernode,numscal),

  // instance of utility class supporting element evaluation
  utils_(ScaTraEleUtilsElch<distype>::Instance(numdofpernode,numscal))
{
  // replace parameter class for standard scalar transport by parameter class for electrochemistry problems
  my::scatraparams_ = DRT::ELEMENTS::ScaTraEleParameterElch::Instance();

  return;
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::EvaluateAction(
    DRT::FaceElement*                   ele,
    Teuchos::ParameterList&             params,
    DRT::Discretization&                discretization,
    SCATRA::BoundaryAction              action,
    std::vector<int>&                   lm,
    Epetra_SerialDenseMatrix&           elemat1_epetra,
    Epetra_SerialDenseMatrix&           elemat2_epetra,
    Epetra_SerialDenseVector&           elevec1_epetra,
    Epetra_SerialDenseVector&           elevec2_epetra,
    Epetra_SerialDenseVector&           elevec3_epetra
)
{
  // determine and evaluate action
  switch(action)
  {
  case SCATRA::bd_calc_elch_boundary_kinetics:
  {
    CalcElchBoundaryKinetics(
        ele,
        params,
        discretization,
        lm,
        elemat1_epetra,
        elevec1_epetra,
        1.
        );

    break;
  }

  case SCATRA::bd_calc_elch_linearize_nernst:
  {
    CalcNernstLinearization(ele,
        params,
        discretization,
        lm,
        elemat1_epetra,
        elevec1_epetra);

    break;
  }

  case SCATRA::bd_calc_elch_cell_voltage:
  {
    CalcCellVoltage(ele,params,discretization,lm,elevec1_epetra);

    break;
  }

  default:
  {
    my::EvaluateAction(
        ele,
        params,
        discretization,
        action,
        lm,
        elemat1_epetra,
        elemat2_epetra,
        elevec1_epetra,
        elevec2_epetra,
        elevec3_epetra
        );

   break;
  }
  } // switch action

  return 0;
}


/*----------------------------------------------------------------------*
 | evaluate electrode kinetics boundary condition                 ehrl  |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::CalcElchBoundaryKinetics(
    DRT::FaceElement*                 ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization&              discretization,
    std::vector<int>&                 lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra,
    const double                      scalar
    )
{
  // get actual values of transported scalars
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if(phinp == Teuchos::null)
    dserror("Cannot get state vector 'phinp'");

  // extract local values from the global vector
  std::vector<double> ephinp(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

  // get history variable (needed for double layer modeling)
  Teuchos::RCP<const Epetra_Vector> hist = discretization.GetState("hist");
  if(phinp == Teuchos::null)
    dserror("Cannot get state vector 'hist'");

  // extract local values from the global vector
  std::vector<double> ehist(lm.size());
  DRT::UTILS::ExtractMyValues(*hist,ehist,lm);

  // get current condition
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
  if(cond == Teuchos::null)
    dserror("Cannot access condition 'ElchBoundaryKinetics'");

  // access parameters of the condition
  const int                 kinetics = cond->GetInt("kinetic model");
  double                    pot0 = cond->GetDouble("pot");
  const int                 curvenum = cond->GetInt("curve");
  const int                 nume = cond->GetInt("e-");
  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neuman condition)
  // but the electrode status is evaluated
  const int                 zerocur = cond->GetInt("zero_cur");
  if(nume < 0)
    dserror("The convention for electrochemical reactions at the electrodes does not allow \n"
        "a negative number of transferred electrons");

  // convention for stoichiometric coefficients s_i:
  // Sum_i (s_i  M_i^(z_i)) -> n e- (n needs to be positive)
  const std::vector<int>*   stoich = cond->GetMutable<std::vector<int> >("stoich");
  if((unsigned int)my::numscal_ != (*stoich).size())
    dserror("Electrode kinetics: number of stoichiometry coefficients %u does not match"
            " the number of ionic species %d", (*stoich).size(), my::numscal_);

  // the classical implementations of kinetic electrode models does not support
  // more than one reagent or product!! There are alternative formulations
  // as e.g. Newman (2004), pp. 205, eq. 8.6 with 8.10
  {
    int reactspecies = 0;
    for(int kk=0; kk<my::numscal_; ++kk)
      reactspecies += abs((*stoich)[kk]);

    if(reactspecies>1 and (kinetics==INPAR::SCATRA::butler_volmer or kinetics == INPAR::SCATRA::butler_volmer_yang1997 or
        kinetics == INPAR::SCATRA::tafel or kinetics == INPAR::SCATRA::linear))
      dserror("Kinetic model Butler-Volmer / Butler-Volmer-Yang / Tafel and Linear: \n"
          "Only one educt and no product is allowed in the implemented version");
  }

  // access input parameter
  const double frt = ElchParams()->FRT();
  if (frt<=0.0)
    dserror("A negative factor frt is not possible by definition");

  // get control parameter from parameter list
  const bool   is_stationary = my::scatraparamstimint_->IsStationary();
  const double time = my::scatraparamstimint_->Time();
  double       timefac = 1.0;
  double       rhsfac  = 1.0;
  // find out whether we shell use a time curve and get the factor
  // this feature can be also used for stationary "pseudo time loops"
  if (curvenum>=0)
  {
    const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
    // adjust potential at metal side accordingly
    pot0 *= curvefac;
  }

  //TODO: Ist es besser über parameter?
  if(!(params.get<bool>("calc_status",false)))
  {
    if(not is_stationary)
    {
      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
      timefac = my::scatraparamstimint_->TimeFac();
      if (timefac < 0.0) dserror("time factor is negative.");
      // for correct scaling of rhs contribution (see below)
      rhsfac =  1/my::scatraparamstimint_->AlphaF();
    }

    if(zerocur == 0)
    {
      EvaluateElchBoundaryKinetics(
          ele,
          elemat1_epetra,
          elevec1_epetra,
          ephinp,
          ehist,
          timefac,
          ele->ParentElement()->Material(),
          cond,
          nume,
          *stoich,
          kinetics,
          pot0,
          frt,
          scalar
      );
    }

    // realize correct scaling of rhs contribution for gen.alpha case
    // with dt*(gamma/alpha_M) = timefac/alpha_F
    // matrix contributions are already scaled correctly with
    // timefac=dt*(gamma*alpha_F/alpha_M)
    elevec1_epetra.Scale(rhsfac);
  }
  else
  {
    // get actual values of transported scalars
    Teuchos::RCP<const Epetra_Vector> phidtnp = discretization.GetState("phidtnp");
    if(phidtnp == Teuchos::null)
      dserror("Cannot get state vector 'ephidtnp'");
    // extract local values from the global vector
    std::vector<double> ephidtnp(lm.size());
    DRT::UTILS::ExtractMyValues(*phidtnp,ephidtnp,lm);

    if(not is_stationary)
    {
      // One-step-Theta:    timefacrhs = theta*dt
      // BDF2:              timefacrhs = 2/3 * dt
      // generalized-alpha: timefacrhs = (gamma/alpha_M) * dt
      timefac = my::scatraparamstimint_->TimeFacRhs();
      if(timefac < 0.)
        dserror("time factor is negative.");
    }

    ElectrodeStatus(
        ele,
        elevec1_epetra,
        params,
        cond,
        ephinp,
        ephidtnp,
        kinetics,
        *stoich,
        nume,
        pot0,
        frt,
        timefac,
        scalar
        );
  }
}


/*----------------------------------------------------------------------*
 | calculate linearization of nernst equation                     ehrl  |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::CalcNernstLinearization(
    DRT::FaceElement*                 ele,
    Teuchos::ParameterList&           params,
    DRT::Discretization&              discretization,
    std::vector<int>&                 lm,
    Epetra_SerialDenseMatrix&         elemat1_epetra,
    Epetra_SerialDenseVector&         elevec1_epetra
    )
{
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition> >("condition");
    if (cond == Teuchos::null) dserror("Cannot access condition 'ElchBoundaryKinetics'");

  const int  kinetics = cond->GetInt("kinetic model");

  // Nernst-BC
  if(kinetics == INPAR::SCATRA::nernst)
  {
    // get actual values of transported scalars
    Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
    if (phinp==Teuchos::null) dserror("Cannot get state vector 'phinp'");

    // extract local values from the global vector
    std::vector<double> ephinp(lm.size());
    DRT::UTILS::ExtractMyValues(*phinp,ephinp,lm);

    // access parameters of the condition
    double       pot0 = cond->GetDouble("pot");
    const int    curvenum = cond->GetInt("curve");
    const int    nume = cond->GetInt("e-");
    const double e0 = cond->GetDouble("e0");
    const double c0 = cond->GetDouble("c0");

    if(nume < 0)
      dserror("The convention for electrochemical reactions at the electrodes does not allow \n"
          "a negative number of transferred electrons");

    const std::vector<int>*   stoich = cond->GetMutable<std::vector<int> >("stoich");
    if((unsigned int)my::numscal_ != (*stoich).size())
      dserror("Electrode kinetics: number of stoichiometry coefficients %u does not match"
              " the number of ionic species %d", (*stoich).size(), my::numscal_);

    // access input parameter
    const double frt = ElchParams()->FRT();
    if (frt<0.0)
      dserror("A negative factor frt is not possible by definition");

    const double time = my::scatraparamstimint_->Time();

    if (curvenum>=0)
    {
      const double curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
      // adjust potential at metal side accordingly
      pot0 *= curvefac;
    }

    // concentration values of reactive species at element nodes
    std::vector<LINALG::Matrix<my::nen_,1> > conreact(my::numscal_);

    // el. potential values at element nodes
    LINALG::Matrix<my::nen_,1> pot(true);

    for (int inode=0; inode< my::nen_;++inode)
    {
      for(int kk=0; kk<my::numscal_; kk++)
        conreact[kk](inode) = ephinp[inode*my::numdofpernode_+kk];

      pot(inode) = ephinp[inode*my::numdofpernode_+my::numscal_];
    }

    // concentration of active species at integration point
    std::vector<double> conint(my::numscal_,0.0);

    // index of reactive species (starting from zero)
    // loop over all scalars
    for(int k=0; k<my::numscal_; ++k)
    {
      if((*stoich)[k]==0)
      continue;

     /*----------------------------------------------------------------------*
      |               start loop over integration points                     |
      *----------------------------------------------------------------------*/
      const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);
      for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
      {
        const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

        // elch-specific values at integration point:
        // concentration is evaluated at all GP since some reaction models depend on all concentrations
        for(int kk=0;kk<my::numscal_;++kk)
          conint[kk] = my::funct_.Dot(conreact[kk]);

        // el. potential at integration point
        const double potint = my::funct_.Dot(pot);

        if (c0 < EPS12) dserror("reference concentration is too small (c0 < 1.0E-12) : %f",c0);

         for (int vi=0; vi<my::nen_; ++vi)
         {
           for (int ui=0; ui<my::nen_; ++ui)
           {
             elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+k) += fac*my::funct_(vi)/(frt*conint[k]*nume)*my::funct_(ui);
             elemat1_epetra(vi*my::numdofpernode_+my::numscal_,ui*my::numdofpernode_+my::numscal_) += fac*my::funct_(vi)*my::funct_(ui);
           }

           // -----right-hand-side
           elevec1_epetra[vi*my::numdofpernode_+my::numscal_] += fac*my::funct_(vi)*(pot0-e0-potint-log(conint[k]/c0)/(frt*nume));
         }
      } // end of loop over integration points gpid
    } // end loop over scalars
  }  // end if(kinetics == INPAR::SCATRA::nernst)
}


/*----------------------------------------------------------------------*
 | calculate cell voltage                                    fang 01/15 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::CalcCellVoltage(
    const DRT::Element*               ele,              //!< the element we are dealing with
    Teuchos::ParameterList&           params,           //!< parameter list
    DRT::Discretization&              discretization,   //!< discretization
    const std::vector<int>&           lm,               //!< location vector
    Epetra_SerialDenseVector&         scalars           //!< result vector for scalar integrals to be computed
    )
{
  // get global state vector
  Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if(phinp == Teuchos::null)
    dserror("Cannot get state vector \"phinp\"!");

  // extract local nodal values of electric potential from global state vector
  std::vector<double> ephinpvec(lm.size());
  DRT::UTILS::ExtractMyValues(*phinp,ephinpvec,lm);
  LINALG::Matrix<my::nen_,1> epotnp(true);
  for(int inode=0; inode<my::nen_; ++inode)
    epotnp(inode,0) = ephinpvec[inode*my::numdofpernode_+my::numscal_];

  // initialize variables for electric potential and domain integrals
  double intpotential(0.);
  double intdomain(0.);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints,iquad);

    // calculate potential and domain integrals
    for (int vi=0; vi<my::nen_; ++vi)
    {
      const double funct_vi_fac = my::funct_(vi)*fac;

      // potential integral
      intpotential += funct_vi_fac*epotnp(vi,0);

      // domain integral
      intdomain += funct_vi_fac;
    }
  } // loop over integration points

  // safety check
  if(scalars.Length() != 2)
    dserror("Result vector for cell voltage computation has invalid length!");

  // write results for electric potential and domain integrals into result vector
  scalars(0) = intpotential;
  scalars(1) = intdomain;

  return;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition (private) gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::EvaluateElchBoundaryKinetics(
    const DRT::Element*                 ele,        ///< the actual boundary element
    Epetra_SerialDenseMatrix&           emat,       ///< element-matrix
    Epetra_SerialDenseVector&           erhs,       ///< element-rhs
    const std::vector<double>&          ephinp,     ///< actual conc. and pot. values
    const std::vector<double>&          ehist,      ///< element history vector
    double                              timefac,    ///< time factor
    Teuchos::RCP<const MAT::Material>   material,   ///< the material
    Teuchos::RCP<DRT::Condition>        cond,       ///< the condition
    const int                           nume,       ///< number of transferred electrons
    const std::vector<int>              stoich,     ///< stoichiometry of the reaction
    const int                           kinetics,   ///< desired electrode kinetics model
    const double                        pot0,       ///< actual electrode potential on metal side
    const double                        frt,        ///< factor F/RT
    const double                        scalar      ///< scaling factor for element matrix and residual contributions
)
{
  // for pre-multiplication of i0 with 1/(F z_k)
  const double faraday = INPAR::ELCH::faraday_const;    // unit of F: C/mol or mC/mmol or µC/µmol

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  std::vector<LINALG::Matrix<my::nen_,1> > conreact(my::numscal_);

  // el. potential values at element nodes
  LINALG::Matrix<my::nen_,1> pot(true);

  // element history vector for potential at electrode
  LINALG::Matrix<my::nen_,1> phihist(true);

  for (int inode=0; inode< my::nen_;++inode)
  {
    for(int kk=0; kk<my::numscal_; kk++)
      conreact[kk](inode) = ephinp[inode*my::numdofpernode_+kk];

    pot(inode) = ephinp[inode*my::numdofpernode_+my::numscal_];
    phihist(inode) = ehist[inode*my::numdofpernode_+my::numscal_];
  }

  // loop over all scalars
  for(int k=0; k<my::numscal_; ++k)
  {
    if(stoich[k]==0)
      continue;

    //(- N^(d+m)*n) = j = s_k / (nume * faraday * z_e-) * i
    //                  = s_k / (nume * faraday * (-1)) * i
    //                    |_______fns_________________|
    // see, e.g. in Ehrl et al., "A computational approach for the simulation of natural convection in
    // electrochemical cells", JCP, 2012
    double fns = -1.0/faraday/nume;
    // stoichiometry as a consequence of the reaction convention
    fns*=stoich[k];

    // get valence of the single reactant
    const double valence_k = GetValence(material,k);

   /*----------------------------------------------------------------------*
    |               start loop over integration points                     |
    *----------------------------------------------------------------------*/
    for(int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

      // get boundary porosity from condition if available, or set equal to volume porosity otherwise
      double epsilon = cond->GetDouble("epsilon");
      if(epsilon == -1)
        epsilon = scalar;
      else if(epsilon <= 0 or epsilon > 1)
        dserror("Boundary porosity has to be between 0 and 1, or -1 by default!");

      // call utility class for element evaluation
      utils_->EvaluateElchKineticsAtIntegrationPoint(
          ele,
          emat,
          erhs,
          conreact,
          pot,
          phihist,
          timefac,
          fac,
          my::funct_,
          material,
          cond,
          nume,
          stoich,
          valence_k,
          kinetics,
          pot0,
          frt,
          fns,
          epsilon,
          k
      );
    } // loop over integration points
  } // loop over all scalars

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::EvaluateElchBoundaryKinetics


/*----------------------------------------------------------------------*
 | calculate electrode kinetics status information             gjb 01/09|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::ElectrodeStatus(
    const DRT::Element*           ele,        ///< the actual boundary element
    Epetra_SerialDenseVector&     scalars,    ///< scalars to be computed
    Teuchos::ParameterList&       params,     ///< the parameter list
    Teuchos::RCP<DRT::Condition>  cond,       ///< the condition
    const std::vector<double>&    ephinp,     ///< current conc. and potential values
    const std::vector<double>&    ephidtnp,   ///< time derivative vector evaluated at t_{n+1}
    const int                     kinetics,   ///< desired electrode kinetics model
    const std::vector<int>        stoich,     ///< stoichiometry of the reaction
    const int                     nume,       ///<  number of transferred electrons
    const double                  pot0,       ///< actual electrode potential on metal side at t_{n+1}
    const double                  frt,        ///< factor F/RT
    const double                  timefac,    ///< factor due to time discretization
    const double                  scalar      ///< scaling factor for current related quantities
)
{
  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values at the electrode.
  // A different approach is not possible (without major hacks) since the time-integration scheme is
  // necessary to perform galvanostatic simulations, for instance.
  // Think about: double layer effects for genalpha time-integration scheme

  double faraday = INPAR::ELCH::faraday_const;    // unit of F: C/mol or mC/mmol or µC/µmol

  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neuman condition)
  // but the electrode status is evaluated
  const int zerocur = cond->GetInt("zero_cur");

  // get boundary porosity from condition if available, or set equal to volume porosity otherwise
  double epsilon = cond->GetDouble("epsilon");
  if(epsilon == -1)
    epsilon = scalar;
  else if(epsilon <= 0 or epsilon > 1)
    dserror("Boundary porosity has to be between 0 and 1, or -1 by default!");

  // get variables with their current values
  // TODO
  // current integrals: (i = epsilon i^E ) is calculated in case of porous media
  double currentintegral(0.);
  double currentdlintegral(0.);
  double boundaryint(0.);
  double electpotentialint(0.);
  double overpotentialint(0.);
  double electdiffpotint(0.);
  double opencircuitpotint(0.);
  double concentrationint(0.);
  double currderiv(0.);
  double currentresidual(0.);
  double boundaryint_porous(0.);   // only implemented for Butler-Volmer and Butler-Volmer-Yang so far!

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_> intpoints(SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // concentration values of reactive species at element nodes
  std::vector<LINALG::Matrix<my::nen_,1> > conreact(my::numscal_);

  // el. potential values at element nodes
  LINALG::Matrix<my::nen_,1> pot(true);
  LINALG::Matrix<my::nen_,1> potdtnp(true);
  for (int inode=0; inode< my::nen_;++inode)
  {
    for (int kk=0; kk<my::numscal_; ++kk)
      conreact[kk](inode) = ephinp[inode*my::numdofpernode_+kk];

    pot(inode) = ephinp[inode*my::numdofpernode_+my::numscal_];
    potdtnp(inode) = ephidtnp[inode*my::numdofpernode_+my::numscal_];
  }

  // concentration of active species at integration point
  std::vector<double> conint(my::numscal_,0.0);

  bool statistics = false;
  // index of reactive species (starting from zero)
  for(int k=0; k<my::numscal_; ++k)
  {
    // only the first oxidized species O is considered for statistics
    // statistics of other species result directly from the oxidized species (current density, ...)
    // or need to be implemented (surface concentration, OCV, ...)
    if(stoich[k]>=0)
      continue;

    statistics = true;

    // loop over integration points
    for (int gpid=0; gpid<intpoints.IP().nquad; gpid++)
    {
      const double fac = my::EvalShapeFuncAndIntFac(intpoints,gpid);

      // elch-specific values at integration point:
      for (int kk=0; kk<my::numscal_; ++kk)
        conint[kk] = my::funct_.Dot(conreact[kk]);

      // el. potential at integration point
      const double potint = my::funct_.Dot(pot);

      // history term of el. potential at integration point
      const double potdtnpint = my::funct_.Dot(potdtnp);

      // electrode potential difference epd at integration point
      double epd = (pot0 - potint);

      // linearization of current w.r.t applied electrode potential "pot0"
      double linea(0.0);

      // concentration-dependent Butler-Volmer law(s)
      switch(kinetics)
      {
      case INPAR::SCATRA::butler_volmer:
      case INPAR::SCATRA::butler_volmer_yang1997:
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha_a");
        const double alphac = cond->GetDouble("alpha_c");
        double       i0 = cond->GetDouble("i0");
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);

        const double dlcap = cond->GetDouble("dl_spec_cap");
        double pot0dtnp = 0.0;
        double pot0hist = 0.0;
        if(dlcap!=0.0)
        {
          pot0dtnp = cond->GetDouble("pot0dtnp");
          pot0hist = cond->GetDouble("pot0hist");
        }

        // opencircuit potential is assumed to be zero here
        double ocp = 0.0;
        // surface overpotential based on opencircuit potential
        double eta = 0.0;
        // electrode potential
        double elepot = 0.0;

        if(zerocur==0)
        {
          eta = epd - ocp;
          elepot = pot0;
        }
        else if(zerocur==1)
        {
          elepot = potint+ocp;
          epd  = ocp;
        }
        else
          dserror("The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

        double expterm(0.0);
        if (kinetics==INPAR::SCATRA::butler_volmer)
        {
          // general Butler-Volmer
          expterm = std::pow(conint[k]/refcon,gamma) * (exp(alphaa*frt*eta)-exp((-alphac)*frt*eta));
          linea = std::pow(conint[k]/refcon,gamma) * frt*((alphaa*exp(alphaa*frt*eta)) + (alphac*exp((-alphac)*frt*eta)));
        }
        if (kinetics==INPAR::SCATRA::butler_volmer_yang1997)
        {
          if (((conint[k]/refcon)<EPS13) && (gamma < 1.0))
          {// prevents NaN's in the current density evaluation
            expterm = (exp(alphaa*frt*eta)-(pow(EPS13/refcon,gamma)*exp((-alphac)*frt*eta)));
            linea = ((alphaa)*frt*exp(alphaa*frt*eta))+(pow(EPS13/refcon,gamma)*alphac*frt*exp((-alphac)*frt*eta));
          }
          else
          {
            expterm = (exp(alphaa*frt*eta)-(pow(conint[k]/refcon,gamma)*exp((-alphac)*frt*eta)));
            linea = ((alphaa)*frt*exp(alphaa*frt*eta))+(pow(conint[k]/refcon,gamma)*alphac*frt*exp((-alphac)*frt*eta));
          }
        }

        // scan for NaNs due to negative concentrations under exponent gamma
        if (std::isnan(expterm) or std::isnan(linea))
          dserror("NaN detected in electrode status calculation");

        // compute integrals
        electpotentialint += elepot*fac;
        overpotentialint += eta*fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        currentintegral += epsilon*i0*expterm*fac; // the negative(!) normal flux density
        boundaryint += fac;
        boundaryint_porous += fac*epsilon;
        concentrationint += epsilon*conint[k]*fac;

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += epsilon*i0*linea*timefac*fac;
        currentresidual += epsilon*i0 * expterm * timefac *fac;

        if (dlcap != 0.0)
        {
          currentdlintegral+=fac*dlcap*(pot0dtnp-potdtnpint);

          // add contributions due to double-layer capacitance
          // positive due to redefinition of the exchange current density
          currderiv += epsilon*fac*dlcap;
          currentresidual += epsilon*fac*dlcap*(pot0-pot0hist-(timefac*potdtnpint));
        }
        break;
      }
      // Tafel law:
      // implementation of cathodic path: i_n = i_0 * (-exp(-alpha * frt* eta)
      // -> cathodic reaction path: i_0 > 0 and alpha > 0
      // -> anodic reacton path:    i_0 < 0 and alpha < 0
      case INPAR::SCATRA::tafel:
      {
        // read model-specific parameter
        const double alpha = cond->GetDouble("alpha");
        double       i0 = cond->GetDouble("i0");

        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
        const double dlcap = cond->GetDouble("dl_spec_cap");
        if(dlcap!=0.0) dserror("double layer charging is not implemented for Tafel electrode kinetics");

        // opencircuit potential is assumed to be zero here
        double ocp = 0.0;
        // surface overpotential based on opencircuit potential
        double eta = 0.0;
        // electrode potential
        double elepot = 0.0;

        if(zerocur==0)
        {
          eta = epd - ocp;
          elepot = pot0;
        }
        else if(zerocur==1)
        {
          elepot = potint+ocp;
          epd  = ocp;
        }
        else
          dserror("The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

        const double expterm = std::pow(conint[k]/refcon,gamma) * (-exp((-alpha)*frt*eta));
        linea = std::pow(conint[k]/refcon,gamma) * frt*(alpha*exp((-alpha)*frt*eta));
        // compute integrals
        electpotentialint += elepot*fac;
        overpotentialint += eta * fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        currentintegral += i0 * expterm * fac; // the negative(!) normal flux density
        boundaryint += fac;
        boundaryint_porous += fac*epsilon;
        concentrationint += conint[k]*fac;

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += scalar*i0*linea*timefac*fac;
        currentresidual += scalar*i0*expterm*timefac*fac;

        break;
      }
      // linear law:  i_n = frt*i_0*((alpha_a+alpha_c)*(V_M - phi)) -> changed 10/13
      // previously implemented: i_n = i_0*(alphaa*frt*(V_M - phi) + 1.0)
      //                         -> linearization in respect to anodic branch!!
      //                         this is not the classical version of a linear electrode kinetics law
      case INPAR::SCATRA::linear:
      {
        // read model-specific parameter
        const double alphaa = cond->GetDouble("alpha");
        double       i0 = cond->GetDouble("i0");
        if (i0 < -EPS14) dserror("i0 is negative, \n"
                                 "a positive definition is necessary due to generalized reaction models: %f",i0);
        const double gamma = cond->GetDouble("gamma");
        const double refcon = cond->GetDouble("refcon");
        if (refcon < EPS12) dserror("reference concentration is too small: %f",refcon);
        const double dlcap = cond->GetDouble("dl_spec_cap");
        double pot0dtnp = 0.0;
        double pot0hist = 0.0;
        if(dlcap!=0.0)
        {
          pot0dtnp = cond->GetDouble("pot0dtnp");
          pot0hist = cond->GetDouble("pot0hist");
        }

        // opencircuit potential is assumed to be zero here
        double ocp = 0.0;
        // surface overpotential based on opencircuit potential
        double eta = 0.0;
        // electrode potential
        double elepot = 0.0;

        if(zerocur==0)
        {
          eta = epd - ocp;
          elepot = pot0;
        }
        else if(zerocur==1)
        {
          elepot = potint+ocp;
          epd  = ocp;
        }
        else
          dserror("The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

        // compute integrals
        electpotentialint += elepot*fac;
        overpotentialint += eta * fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        currentintegral += i0 * pow(conint[k]/refcon,gamma)*(alphaa*frt*eta) * fac; // the negative(!) normal flux density
        boundaryint += fac;
        boundaryint_porous += fac*epsilon;
        concentrationint += conint[k]*fac;

        // tangent and rhs (= negative residual) for galvanostatic equation
        linea = std::pow(conint[k]/refcon,gamma)*(alphaa*frt);
        currderiv += scalar*i0*linea*timefac*fac;
        currentresidual += scalar*i0*pow(conint[k]/refcon,gamma)*(alphaa*frt*eta)*timefac*fac;

        if (dlcap != 0.0)
        {
          currentdlintegral+=fac*dlcap*(pot0dtnp-potdtnpint);

          // add contributions due to double-layer capacitance
          // positive due to redefinition of the exchange current density
          currderiv += fac*dlcap;
          currentresidual += fac*dlcap*(pot0-pot0hist-(timefac*potdtnpint));
        }

        break;
      }
      case INPAR::SCATRA::butler_volmer_newman:
      {
        // "Electrochemical systems"
        // Newman ad Thomas-Alyea, 2004
        // General stoichiometry: pp. 212-213, e.q. 8.26
        // consideration of a elementary step of the form:
        // Sum_i s_i M_i ->  ne-
        // n is one if charge transfer is involved, multiple electron transfers "being unlikely in
        // an elementary step

        const double k_a = cond->GetDouble("k_a");
        const double k_c = cond->GetDouble("k_c");
        const double beta = cond->GetDouble("beta");
        const double dlcap = cond->GetDouble("dl_spec_cap");
        if(dlcap!=0.0) dserror("double layer charging is not implemented for Butler-Volmer-Newman electrode kinetics");
        if(zerocur!=0) dserror("The electrode kinetics flag zero_cur is not implemented for this specific kinetic model.");

        //reaction order of the cathodic and anodic reactants of ionic species k
        std::vector<int> q(my::numscal_,0);
        std::vector<int> p(my::numscal_,0);

        for(int kk=0; kk<my::numscal_; kk++)
        {
          //according to the convention: anodic reactant is positiv
          if(stoich[kk] > 0)
          {
            q[kk] = 0;
            p[kk] = stoich[kk];
          }
          //according to the convention: cathodic reactant is negative
          else
          {
            q[kk]= -stoich[kk];
            p[kk] = 0;
          }
        }

        // linearization of current w.r.t applied electrode potential "pot0"
        double linea(0.0);
        double expterma(0.0);
        double exptermc(0.0);
        double expterm(0.0);
        double pow_conint_p = 1.0;      //product over i (c_i)^(p_i)
        double pow_conint_q = 1.0;      //product over i (c_i)^(q_i)

        // overpotential based on opencircuit potential
        double eta = 0.0;
        // electrode potential
        double elepot = 0.0;

        // open circuit potential (ocp): time dependent electrode surface concentrations
        // defined in Newman, 2004, p. 211, eq. 8.20
        double ocp = 1/frt/nume*log(k_c/k_a);
        for(int kk=0;kk<my::numscal_;++kk)
        {
          ocp +=  1/frt/nume*(q[kk]-p[kk])*log(conint[kk]);
          //safety check
          if((q[kk]-p[kk])!=-stoich[kk])
            dserror("stoichiometry factors and the factors q,p do not correlate!!");
        }

        if(zerocur==0)
        {
          // overpotential based on open circuit potential
          eta = epd - ocp;
          elepot = pot0;
        }
        else if(zerocur==1)
        {
          elepot = potint+ocp;
          epd  = ocp;
        }
        else
          dserror("The electrode kinetics flag zero_cur has only two options: false (=0) or true (=1).");

        for(int kk=0; kk<my::numscal_; ++kk)
        {
          if ((conint[kk]) < EPS13)
          {
            pow_conint_p *= std::pow(EPS13,p[kk]);
            pow_conint_q *= std::pow(EPS13,q[kk]);
#ifdef DEBUG
            std::cout<<"WARNING: Rel. Conc. of species"<<kk<<" in Butler-Volmer formula is zero/negative: "<<(conint[kk])<<std::endl;
            std::cout<<"-> Replacement value: pow(EPS,p[ispec]) = "<< pow(EPS13,p[kk]) << " pow(1.0E-16,q[i]) = "<< pow(EPS13,q[kk]) <<std::endl;
#endif
          }
          else
          {
            pow_conint_p *= std::pow((conint[kk]),p[kk]);
            pow_conint_q *= std::pow((conint[kk]),q[kk]);
          }
        }
        expterma = exp((1-beta)*nume*frt*epd);
        exptermc = exp(-beta*nume*frt*epd);
        linea =  nume*faraday*(frt*nume*((k_a*(1-beta)*expterma*pow_conint_p)-(k_c*(-1)*beta*exptermc*pow_conint_q)));

        // scan for NaNs due to negative concentrations under exponent gamma
        if (std::isnan(expterm) or std::isnan(linea))
          dserror("NaN detected in electrode status calculation");

        // compute integrals
        currentintegral += scalar* nume*faraday*((k_a*expterma*pow_conint_p)-(k_c*exptermc*pow_conint_q))*fac;
        boundaryint += fac;
        boundaryint_porous += fac*epsilon;
        electpotentialint += elepot * fac;
        overpotentialint += eta * fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        concentrationint += conint[k]*fac;

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += scalar*linea*fac*timefac;
        currentresidual += scalar*nume*faraday*((k_a*expterma*pow_conint_p)-(k_c*exptermc*pow_conint_q))*timefac*fac;

        break;
      }
      case INPAR::SCATRA::butler_volmer_bard:
      {
        // "Electrochemical Methods Fundamentals and Applications"
        // Bard and Faulkner, 2001, pp. 94 ff; pp. 99 eq. 3.4.10
        // reaction model for a one-step, one-electron process
        // O + e -> R
        const double e0 = cond->GetDouble("e0");
        const double k0 = cond->GetDouble("k0");
        const double beta = cond->GetDouble("beta");
        const double c_c0 = cond->GetDouble("c_c0");
        const double c_a0 = cond->GetDouble("c_a0");
        const double dlcap = cond->GetDouble("dl_spec_cap");
        if(dlcap!=0.0) dserror("double layer charging is not implemented for Butler-Volmer-Bard electrode kinetics");
        if(zerocur!=0) dserror("The electrode kinetics flag zero_cur is not implemented for this specific kinetic model.");

        if(nume!=1)
          dserror("electron != 1; \n "
              "this Butler-Volmer-equation (Bard/Faulkner) works for elementary steps (one electron) only!");

        // only one reactant and product are supported by the basic model
        // only stoichiometry of 1
        {
          int check1 = 0;
          int check2 = 0;
          for(int kk=0; kk<my::numscal_;kk++)
          {
            if(abs(stoich[kk])>1)
              dserror("Stoichiometry is larger than 1!! \n"
                      "This is not supported by the reaction model based on Bard");

            check1 += abs(stoich[kk]);
            check2 += stoich[kk];
          }
          if (check1>2 or check1==0)
            dserror("More than one reactant or product defined!! \n"
                "This is not supported by the reaction model based on Bard");

          // In the moment it is not checked if two products (and no reactants) and vis versa are defined
        }

        // reactant or product not a species in the electrolyte
        // -> concentration = 1.0
        double conctermc = 1.0;
        double concterma = 1.0;

        // concentration terms for anodic and cathodic reaction
        // only one reactant and product are supported by the basic model
        // only stoichiometry of 1
        for(int kk=0; kk<my::numscal_;kk++)
        {
          if(stoich[kk]==1)
            concterma = conint[kk]/c_a0;
          else if(stoich[kk]==-1)
            conctermc = conint[kk]/c_c0;
        }

        // equilibrium potential (equilpot):
        // defined in Bard, 2001, p.98, eq. 3.4.3
        const double equilpot = e0 + (log(c_c0/c_a0))/(frt*nume);
        // overpotential based on equilibrium potential
        const double eta_equilpot = epd - equilpot;
        // difference between equilibrium potential and open circuit potential:
        // -> equilpot: depends on initial electrode surface concentration
        // -> ocp:      depends on actual electrode surface concentration

        // open circuit potential (ocp): time dependent electrode surface concentrations
        // defined in Newman, 2004, p. 211, eq. 8.20
        const double ocp = e0 + 1/frt/nume*log(conctermc/concterma);
        // overpotential based on open circuit potential
        const double eta = epd - ocp;

        const double expterma = exp((1-beta)*(frt*nume)*eta_equilpot);
        const double exptermc = exp(-beta*(frt*nume)*eta_equilpot);
        const double linea = concterma*(1-beta)*(frt*nume)*expterma+conctermc*beta*(frt*nume)*exptermc;

        // scan for NaNs due to negative concentrations under exponent gamma
        if (std::isnan(expterma) or std::isnan(exptermc) or std::isnan(linea))
          dserror("NaN detected in electrode status calculation");

        const double i0 = faraday*k0*pow(c_c0,1-beta)*pow(c_a0,beta);

        // compute integrals
        overpotentialint += eta*fac;
        electdiffpotint += epd*fac;
        opencircuitpotint += ocp*fac;
        currentintegral += scalar*i0*(concterma*expterma-conctermc*exptermc)*fac; // the negative(!) normal flux density
        boundaryint += fac;
        boundaryint_porous += fac*epsilon;
        concentrationint += conint[k]*fac;  //concentration-output for the first species only

        // tangent and rhs (= negative residual) for galvanostatic equation
        currderiv += scalar*i0*linea*timefac*fac;
        currentresidual += scalar*i0*(concterma*expterma-conctermc*exptermc)*timefac*fac;

        break;
      } // end Butler-Volmer-Bard
      case INPAR::SCATRA::nernst:
      {
        const double e0 = cond->GetDouble("e0");
        const double c0 = cond->GetDouble("c0");
        if(zerocur!=0) dserror("The electrode kinetics flag zero_cur is not implemented for this specific kinetic model.");

        // compute integrals
        overpotentialint += potint * fac;
        boundaryint += fac;
        boundaryint_porous += fac*epsilon;
        concentrationint += conint[k]*fac;

        opencircuitpotint+= e0 + (log(concentrationint/boundaryint/c0))/(frt*nume);
        opencircuitpotint*=boundaryint;
        break;
      }
      default:
        dserror("Kinetic model not implemented");
        break;
      }
    }  // loop over integration points
    // stop loop over ionic species after one evaluation
    break;
  }  // loop over scalars

  if(statistics == false)
    dserror("There is no oxidized species O (stoich<0) defined in your input file!! \n"
            " Statistics could not be evaluated");

  // write results into result vector
  scalars(0) = currentintegral;
  scalars(1) = currentdlintegral;
  scalars(2) = boundaryint;
  scalars(3) = electpotentialint;
  scalars(4) = overpotentialint;
  scalars(5) = electdiffpotint;
  scalars(6) = opencircuitpotint;
  scalars(7) = concentrationint;
  scalars(8) = currderiv;
  scalars(9) = currentresidual;
  scalars(10) = boundaryint_porous;

  return;
} // DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype>::ElectrodeStatus


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::quad4>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::quad9>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::tri3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::tri6>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::line2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::line3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::nurbs3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<DRT::Element::nurbs9>;
