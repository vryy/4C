/*----------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra boundary terms at integration points

\level 2

 */
/*----------------------------------------------------------------------*/
#include "4C_scatra_ele_boundary_calc_elch.hpp"

#include "4C_global_data.hpp"
#include "4C_lib_discret.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_ele_utils_elch.hpp"
#include "4C_utils_function_of_time.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | protected constructor for singletons                      fang 01/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::ScaTraEleBoundaryCalcElch(
    const int numdofpernode, const int numscal, const std::string& disname)
    :  // constructor of base class
      my::ScaTraEleBoundaryCalc(numdofpernode, numscal, disname),

      // instance of parameter class for electrochemistry problems
      elchparams_(DRT::ELEMENTS::ScaTraEleParameterElch::Instance(disname)),

      // instance of utility class supporting element evaluation
      utils_(ScaTraEleUtilsElch<distype>::Instance(numdofpernode, numscal, disname))
{
}


/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::EvaluateAction(
    DRT::FaceElement* ele,                            //!< boundary element
    Teuchos::ParameterList& params,                   //!< parameter list
    DRT::Discretization& discretization,              //!< discretization
    SCATRA::BoundaryAction action,                    //!< action
    DRT::Element::LocationArray& la,                  //!< location array
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  //!< element matrix 1
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,  //!< element matrix 2
    CORE::LINALG::SerialDenseVector& elevec1_epetra,  //!< element right-hand side vector 1
    CORE::LINALG::SerialDenseVector& elevec2_epetra,  //!< element right-hand side vector 2
    CORE::LINALG::SerialDenseVector& elevec3_epetra   //!< element right-hand side vector 3
)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_elch_linearize_nernst:
    {
      CalcNernstLinearization(ele, params, discretization, la, elemat1_epetra, elevec1_epetra);

      break;
    }

    case SCATRA::BoundaryAction::calc_elch_cell_voltage:
    {
      CalcCellVoltage(ele, params, discretization, la, elevec1_epetra);

      break;
    }

    default:
    {
      my::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra, elemat2_epetra,
          elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch action

  return 0;
}


/*----------------------------------------------------------------------*
 | process an electrode kinetics boundary condition          fang 08/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::CalcElchBoundaryKinetics(
    DRT::FaceElement* ele,                            ///< current element
    Teuchos::ParameterList& params,                   ///< parameter list
    DRT::Discretization& discretization,              ///< discretization
    DRT::Element::LocationArray& la,                  ///< location array
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,  ///< element matrix
    CORE::LINALG::SerialDenseVector& elevec1_epetra,  ///< element right-hand side vector
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // state and history variables at element nodes
  my::ExtractNodeValues(discretization, la);
  std::vector<CORE::LINALG::Matrix<nen_, 1>> ehist(
      my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
  my::ExtractNodeValues(ehist, discretization, la, "hist");

  // get current condition
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null) FOUR_C_THROW("Cannot access condition 'ElchBoundaryKinetics'");

  // access parameters of the condition
  const auto kinetics = cond->parameters().Get<int>("kinetic model");
  auto pot0 = cond->parameters().Get<double>("pot");
  const auto curvenum = cond->parameters().Get<int>("funct");
  const auto nume = cond->parameters().Get<int>("e-");
  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neuman
  // condition) but the electrode status is evaluated
  const auto zerocur = cond->parameters().Get<int>("zero_cur");
  if (nume < 0)
    FOUR_C_THROW(
        "The convention for electrochemical reactions at the electrodes does not allow \n"
        "a negative number of transferred electrons");

  // convention for stoichiometric coefficients s_i:
  // Sum_i (s_i  M_i^(z_i)) -> n e- (n needs to be positive)
  const auto* stoich = &cond->parameters().Get<std::vector<int>>("stoich");
  if ((unsigned int)my::numscal_ != (*stoich).size())
    FOUR_C_THROW(
        "Electrode kinetics: number of stoichiometry coefficients %u does not match"
        " the number of ionic species %d",
        (*stoich).size(), my::numscal_);

  // the classical implementations of kinetic electrode models does not support
  // more than one reagent or product!! There are alternative formulations
  // as e.g. Newman (2004), pp. 205, eq. 8.6 with 8.10
  {
    int reactspecies = 0;
    for (int kk = 0; kk < my::numscal_; ++kk) reactspecies += abs((*stoich)[kk]);

    if (reactspecies > 1 and (kinetics == INPAR::ELCH::butler_volmer or
                                 kinetics == INPAR::ELCH::butler_volmer_yang1997 or
                                 kinetics == INPAR::ELCH::tafel or kinetics == INPAR::ELCH::linear))
      FOUR_C_THROW(
          "Kinetic model Butler-Volmer / Butler-Volmer-Yang / Tafel and Linear: \n"
          "Only one educt and no product is allowed in the implemented version");
  }

  // access input parameter
  const double frt = elchparams_->FRT();

  // get control parameter from parameter list
  const bool is_stationary = my::scatraparamstimint_->IsStationary();
  const double time = my::scatraparamstimint_->Time();
  double timefac = 1.0;
  double rhsfac = 1.0;
  // find out whether we shell use a time curve and get the factor
  // this feature can be also used for stationary "pseudo time loops"
  if (curvenum >= 0)
  {
    const double curvefac =
        GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfTime>(curvenum).Evaluate(
            time);
    // adjust potential at metal side accordingly
    pot0 *= curvefac;
  }

  // TODO: Ist es besser ueber parameter?
  if (!(params.get<bool>("calc_status", false)))
  {
    if (not is_stationary)
    {
      // One-step-Theta:    timefac = theta*dt
      // BDF2:              timefac = 2/3 * dt
      // generalized-alpha: timefac = (gamma*alpha_F/alpha_M) * dt
      timefac = my::scatraparamstimint_->TimeFac();
      if (timefac < 0.0) FOUR_C_THROW("time factor is negative.");
      // for correct scaling of rhs contribution (see below)
      rhsfac = 1 / my::scatraparamstimint_->AlphaF();
    }

    if (zerocur == 0)
    {
      EvaluateElchBoundaryKinetics(ele, elemat1_epetra, elevec1_epetra, my::ephinp_, ehist, timefac,
          ele->ParentElement()->Material(), cond, nume, *stoich, kinetics, pot0, frt, scalar);
    }

    // realize correct scaling of rhs contribution for gen.alpha case
    // with dt*(gamma/alpha_M) = timefac/alpha_F
    // matrix contributions are already scaled correctly with
    // timefac=dt*(gamma*alpha_F/alpha_M)
    elevec1_epetra.scale(rhsfac);
  }
  else
  {
    // extract local values from the global vector
    std::vector<CORE::LINALG::Matrix<nen_, 1>> ephidtnp(
        my::numdofpernode_, CORE::LINALG::Matrix<nen_, 1>(true));
    my::ExtractNodeValues(ephidtnp, discretization, la, "phidtnp");

    if (not is_stationary)
    {
      // One-step-Theta:    timefacrhs = theta*dt
      // BDF2:              timefacrhs = 2/3 * dt
      // generalized-alpha: timefacrhs = (gamma/alpha_M) * dt
      timefac = my::scatraparamstimint_->TimeFacRhs();
      if (timefac < 0.) FOUR_C_THROW("time factor is negative.");
    }

    EvaluateElectrodeStatus(ele, elevec1_epetra, params, cond, my::ephinp_, ephidtnp, kinetics,
        *stoich, nume, pot0, frt, timefac, scalar);
  }

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::CalcElchBoundaryKinetics


/*----------------------------------------------------------------------*
 | calculate linearization of nernst equation                     ehrl  |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::CalcNernstLinearization(
    DRT::FaceElement* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra)
{
  Teuchos::RCP<DRT::Condition> cond = params.get<Teuchos::RCP<DRT::Condition>>("condition");
  if (cond == Teuchos::null) FOUR_C_THROW("Cannot access condition 'ElchBoundaryKinetics'");

  const auto kinetics = cond->parameters().Get<int>("kinetic model");

  // Nernst-BC
  if (kinetics == INPAR::ELCH::nernst)
  {
    // extract local values from the global vector
    my::ExtractNodeValues(discretization, la);

    // access parameters of the condition
    auto pot0 = cond->parameters().Get<double>("pot");
    const auto curvenum = cond->parameters().Get<int>("funct");
    const auto nume = cond->parameters().Get<int>("e-");
    const auto e0 = cond->parameters().Get<double>("e0");
    const auto c0 = cond->parameters().Get<double>("c0");

    if (nume < 0)
      FOUR_C_THROW(
          "The convention for electrochemical reactions at the electrodes does not allow \n"
          "a negative number of transferred electrons");

    const auto* stoich = &cond->parameters().Get<std::vector<int>>("stoich");
    if ((unsigned int)my::numscal_ != (*stoich).size())
      FOUR_C_THROW(
          "Electrode kinetics: number of stoichiometry coefficients %u does not match"
          " the number of ionic species %d",
          (*stoich).size(), my::numscal_);

    // access input parameter
    const double frt = elchparams_->FRT();

    const double time = my::scatraparamstimint_->Time();

    if (curvenum >= 0)
    {
      const double curvefac =
          GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfTime>(curvenum).Evaluate(
              time);
      // adjust potential at metal side accordingly
      pot0 *= curvefac;
    }

    // concentration of active species at integration point
    std::vector<double> conint(my::numscal_, 0.0);

    // index of reactive species (starting from zero)
    // loop over all scalars
    for (int k = 0; k < my::numscal_; ++k)
    {
      if ((*stoich)[k] == 0) continue;

      /*----------------------------------------------------------------------*
       |               start loop over integration points                     |
       *----------------------------------------------------------------------*/
      const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
          SCATRA::DisTypeToOptGaussRule<distype>::rule);
      for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
      {
        const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

        // elch-specific values at integration point:
        // concentration is evaluated at all GP since some reaction models depend on all
        // concentrations
        for (int kk = 0; kk < my::numscal_; ++kk) conint[kk] = my::funct_.Dot(my::ephinp_[kk]);

        // el. potential at integration point
        const double potint = my::funct_.Dot(my::ephinp_[my::numscal_]);

        if (c0 < 1e-12)
          FOUR_C_THROW("reference concentration is too small (c0 < 1.0E-12) : %f", c0);

        for (int vi = 0; vi < nen_; ++vi)
        {
          for (int ui = 0; ui < nen_; ++ui)
          {
            elemat1_epetra(vi * my::numdofpernode_ + my::numscal_, ui * my::numdofpernode_ + k) +=
                fac * my::funct_(vi) / (frt * conint[k] * nume) * my::funct_(ui);
            elemat1_epetra(vi * my::numdofpernode_ + my::numscal_,
                ui * my::numdofpernode_ + my::numscal_) += fac * my::funct_(vi) * my::funct_(ui);
          }

          // -----right-hand-side
          elevec1_epetra[vi * my::numdofpernode_ + my::numscal_] +=
              fac * my::funct_(vi) * (pot0 - e0 - potint - log(conint[k] / c0) / (frt * nume));
        }
      }  // end of loop over integration points gpid
    }    // end loop over scalars
  }      // end if(kinetics == INPAR::ELCH::nernst)
}


/*----------------------------------------------------------------------*
 | calculate cell voltage                                    fang 01/15 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::CalcCellVoltage(
    const DRT::Element* ele,                  //!< the element we are dealing with
    Teuchos::ParameterList& params,           //!< parameter list
    DRT::Discretization& discretization,      //!< discretization
    DRT::Element::LocationArray& la,          //!< location array
    CORE::LINALG::SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
)
{
  // extract local nodal values of electric potential from global state vector
  my::ExtractNodeValues(discretization, la);

  // initialize variables for electric potential and domain integrals
  double intpotential(0.);
  double intdomain(0.);

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndIntFac(intpoints, iquad);

    // calculate potential and domain integrals
    for (int vi = 0; vi < nen_; ++vi)
      // potential integral
      intpotential += my::ephinp_[my::numscal_](vi, 0) * my::funct_(vi) * fac;

    // domain integral
    intdomain += fac;
  }  // loop over integration points

  // safety check
  if (scalars.length() != 2)
    FOUR_C_THROW("Result vector for cell voltage computation has invalid length!");

  // write results for electric potential and domain integrals into result vector
  scalars(0) = intpotential;
  scalars(1) = intdomain;

  return;
}


/*----------------------------------------------------------------------*
 | evaluate an electrode kinetics boundary condition          gjb 01/09 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::EvaluateElchBoundaryKinetics(
    const DRT::Element* ele,                ///< current element
    CORE::LINALG::SerialDenseMatrix& emat,  ///< element matrix
    CORE::LINALG::SerialDenseVector& erhs,  ///< element right-hand side vector
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ehist,  ///< nodal history vector
    double timefac,                                           ///< time factor
    Teuchos::RCP<const CORE::MAT::Material> material,         ///< material
    Teuchos::RCP<DRT::Condition> cond,  ///< electrode kinetics boundary condition
    const int nume,                     ///< number of transferred electrons
    const std::vector<int> stoich,      ///< stoichiometry of the reaction
    const int kinetics,                 ///< desired electrode kinetics model
    const double pot0,                  ///< electrode potential on metal side
    const double frt,                   ///< factor F/RT
    const double scalar  ///< scaling factor for element matrix and right-hand side contributions
)
{
  // for pre-multiplication of i0 with 1/(F z_k)
  const double faraday = elchparams_->Faraday();

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over all scalars
  for (int k = 0; k < my::numscal_; ++k)
  {
    if (stoich[k] == 0) continue;

    //(- N^(d+m)*n) = j = s_k / (nume * faraday * z_e-) * i
    //                  = s_k / (nume * faraday * (-1)) * i
    //                    |_______fns_________________|
    // see, e.g. in Ehrl et al., "A computational approach for the simulation of natural convection
    // in electrochemical cells", JCP, 2012
    double fns = -1.0 / faraday / nume;
    // stoichiometry as a consequence of the reaction convention
    fns *= stoich[k];

    // get valence of the single reactant
    const double valence_k = GetValence(material, k);

    /*----------------------------------------------------------------------*
     |               start loop over integration points                     |
     *----------------------------------------------------------------------*/
    for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
    {
      const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

      // get boundary porosity from condition if available, or set equal to volume porosity
      // otherwise
      auto epsilon = cond->parameters().Get<double>("epsilon");
      if (epsilon == -1)
        epsilon = scalar;
      else if (epsilon <= 0 or epsilon > 1)
        FOUR_C_THROW("Boundary porosity has to be between 0 and 1, or -1 by default!");

      // call utility class for element evaluation
      utils_->EvaluateElchKineticsAtIntegrationPoint(ele, emat, erhs, ephinp, ehist, timefac, fac,
          my::funct_, cond, nume, stoich, valence_k, kinetics, pot0, frt, fns, epsilon, k);
    }  // loop over integration points
  }    // loop over all scalars

  return;
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::EvaluateElchBoundaryKinetics


/*----------------------------------------------------------------------*
 | evaluate electrode kinetics status information             gjb 01/09 |
 *----------------------------------------------------------------------*/
template <CORE::FE::CellType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::EvaluateElectrodeStatus(
    const DRT::Element* ele,                   ///< current element
    CORE::LINALG::SerialDenseVector& scalars,  ///< scalars to be integrated
    Teuchos::ParameterList& params,            ///< parameter list
    Teuchos::RCP<DRT::Condition> cond,         ///< condition
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>&
        ephinp,  ///< nodal values of concentration and electric potential
    const std::vector<CORE::LINALG::Matrix<nen_, 1>>& ephidtnp,  ///< nodal time derivative vector
    const int kinetics,             ///< desired electrode kinetics model
    const std::vector<int> stoich,  ///< stoichiometry of the reaction
    const int nume,                 ///< number of transferred electrons
    const double pot0,              ///< electrode potential on metal side
    const double frt,               ///< factor F/RT
    const double timefac,           ///< time factor
    const double scalar             ///< scaling factor for current related quantities
)
{
  // Warning:
  // Specific time integration parameter are set in the following function.
  // In the case of a genalpha-time integration scheme the solution vector phiaf_ at time n+af
  // is passed to the element evaluation routine. Therefore, the electrode status is evaluate at a
  // different time (n+af) than our output routine (n+1), resulting in slightly different values at
  // the electrode. A different approach is not possible (without major hacks) since the
  // time-integration scheme is necessary to perform galvanostatic simulations, for instance. Think
  // about: double layer effects for genalpha time-integration scheme

  // if zero=1=true, the current flow across the electrode is zero (comparable to do-nothing Neumann
  // condition) but the electrode status is evaluated
  const int zerocur = cond->parameters().Get<int>("zero_cur");

  // get boundary porosity from condition if available, or set equal to volume porosity otherwise
  auto epsilon = cond->parameters().Get<double>("epsilon");
  if (epsilon == -1)
    epsilon = scalar;
  else if (epsilon <= 0 or epsilon > 1)
    FOUR_C_THROW("Boundary porosity has to be between 0 and 1, or -1 by default!");

  // integration points and weights
  const CORE::FE::IntPointsAndWeights<nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  bool statistics = false;

  // index of reactive species (starting from zero)
  for (int k = 0; k < my::numscal_; ++k)
  {
    // only the first oxidized species O is considered for statistics
    // statistics of other species result directly from the oxidized species (current density, ...)
    // or need to be implemented (surface concentration, OCV, ...)
    if (stoich[k] >= 0) continue;

    statistics = true;

    // loop over integration points
    for (int gpid = 0; gpid < intpoints.IP().nquad; gpid++)
    {
      const double fac = my::EvalShapeFuncAndIntFac(intpoints, gpid);

      // call utility class for element evaluation
      utils_->EvaluateElectrodeStatusAtIntegrationPoint(ele, scalars, params, cond, ephinp,
          ephidtnp, my::funct_, zerocur, kinetics, stoich, nume, pot0, frt, timefac, fac, epsilon,
          k);
    }  // loop over integration points

    // stop loop over ionic species after one evaluation (see also comment above)
    break;
  }  // loop over scalars

  // safety check
  if (statistics == false)
    FOUR_C_THROW(
        "There is no oxidized species O (stoich<0) defined in your input file!! \n"
        " Statistics could not be evaluated");
}  // DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<distype, probdim>::EvaluateElectrodeStatus


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::quad4, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::quad8, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::quad9, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::tri6, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::line3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::nurbs3, 2>;
template class DRT::ELEMENTS::ScaTraEleBoundaryCalcElch<CORE::FE::CellType::nurbs9, 3>;

FOUR_C_NAMESPACE_CLOSE
