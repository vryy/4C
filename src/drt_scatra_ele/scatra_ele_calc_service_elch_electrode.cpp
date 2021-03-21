/*--------------------------------------------------------------------------*/
/*! \file

\brief evaluation of scatra elements for conservation of mass concentration and electronic charge
within electrodes

\level 2

*/
/*--------------------------------------------------------------------------*/
#include "scatra_ele_calc_elch_electrode.H"
#include "scatra_ele_parameter_elch_manifold.H"
#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/material.H"

/*----------------------------------------------------------------------*
 | evaluate action                                           fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::calc_elch_electrode_soc_and_c_rate:
    {
      CalculateElectrodeSOCAndCRate(ele, discretization, la, elevec1_epetra);

      break;
    }
    default:
    {
      myelch::EvaluateAction(ele, params, discretization, action, la, elemat1_epetra,
          elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);

      break;
    }
  }  // switch(action)

  return 0;
}


/*----------------------------------------------------------------------------------------------------------*
 | validity check with respect to input parameters, degrees of freedom, number of scalars etc. fang
 02/15 |
 *----------------------------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CheckElchElementParameter(
    DRT::Element* ele  //!< current element
)
{
  // safety checks
  if (ele->Material()->MaterialType() != INPAR::MAT::m_electrode) dserror("Invalid material type!");

  if (my::numscal_ != 1) dserror("Invalid number of transported scalars!");
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CheckElchElementParameter


/*----------------------------------------------------------------------*
 | get conductivity                                          fang 02/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::GetConductivity(
    const enum INPAR::ELCH::EquPot equpot,  //!< type of closing equation for electric potential
    double& sigma_all,                      //!< conductivity of electrolyte solution
    std::vector<double>& sigma,  //!< conductivity or a single ion + overall electrolyte solution
    bool effCond)
{
  // use precomputed conductivity
  sigma_all = DiffManager()->GetCond();
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::GetConductivity


/*----------------------------------------------------------------------*
 | calculate weighted current density                        fang 07/16 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalculateCurrent(
    LINALG::Matrix<my::nsd_, 1>& q,          //!< flux of species k
    const INPAR::SCATRA::FluxType fluxtype,  //!< type fo flux
    const double fac                         //!< integration factor
)
{
  /*
  Actually, we compute here a weighted (and integrated) form of the current density!
  On time integration level, these contributions are then used to calculate
  an L2-projected representation of the current density.
  Thus, this method here DOES NOT YET provide current density values that are ready to use!

  /                           \
  |                    /   \  |
  | w, -sigma * nabla | phi | |
  |                    \   /  |
  \                           /
  */

  switch (fluxtype)
  {
    case INPAR::SCATRA::flux_diffusive:
    case INPAR::SCATRA::flux_total:
    {
      // ohmic contribution to current density
      q.Update(-DiffManager()->GetCond(), VarManager()->GradPot());
      break;
    }

    default:
    {
      dserror("Invalid flux type!");
      break;
    }
  }
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateCurrent


/*----------------------------------------------------------------------*
 | calculate electrode state of charge and C rate            fang 01/15 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalculateElectrodeSOCAndCRate(
    const DRT::Element* const& ele,             //!< the element we are dealing with
    const DRT::Discretization& discretization,  //!< discretization
    DRT::Element::LocationArray& la,            //!< location array
    Epetra_SerialDenseVector& scalars  //!< result vector for scalar integrals to be computed
)
{
  // safety check
  if (my::numscal_ != 1)
    dserror("Electrode state of charge can only be computed for one transported scalar!");

  // get global state vectors
  const Teuchos::RCP<const Epetra_Vector> phinp = discretization.GetState("phinp");
  if (phinp == Teuchos::null) dserror("Cannot get state vector \"phinp\"!");
  const Teuchos::RCP<const Epetra_Vector> phidtnp = discretization.GetState("phidtnp");
  if (phidtnp == Teuchos::null) dserror("Cannot get state vector \"phidtnp\"!");

  // extract local nodal values from global state vectors
  DRT::UTILS::ExtractMyValues(*phinp, my::ephinp_, la[0].lm_);
  static std::vector<LINALG::Matrix<my::nen_, 1>> ephidtnp(2);
  DRT::UTILS::ExtractMyValues(*phidtnp, ephidtnp, la[0].lm_);

  // initialize variables for integrals of concentration, its time derivative, and domain
  double intconcentration(0.);
  double intconcentrationtimederiv(0.);
  double intdomain(0.);

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  // loop over integration points
  for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

    // calculate integrals of concentration and its time derivative
    for (unsigned vi = 0; vi < my::nen_; ++vi)
    {
      const double vi_fac = my::funct_(vi) * fac;

      // integral of concentration
      intconcentration += my::ephinp_[0](vi) * vi_fac;

      // integral of time derivative of concentration
      intconcentrationtimederiv += ephidtnp[0](vi) * vi_fac;
    }

    // domain integral
    intdomain += fac;
  }  // loop over integration points

  // safety check
  if (not my::scatrapara_->IsAle() and scalars.Length() != 3)
    dserror(
        "Result vector for electrode state of charge and C rate computation has invalid length!");

  // write results for integrals of concentration, its time derivative, and domain into result
  // vector
  scalars(0) = intconcentration;
  scalars(1) = intconcentrationtimederiv;
  scalars(2) = intdomain;

  // additional computations in case of ALE
  if (my::scatrapara_->IsAle())
  {
    // extract velocities
    const Teuchos::RCP<const Epetra_Vector> vel = discretization.GetState(1, "velocity field");
    if (vel == Teuchos::null) dserror("Cannot get state vector \"velocity field\"!");
    DRT::UTILS::ExtractMyValues(*vel, my::evelnp_, la[1].lm_);

    // initialize additional variables for integrals related to velocity divergence
    double intdivv(0.);
    double intcdivv(0.);
    double intvgradc(0.);

    // loop over integration points
    for (int iquad = 0; iquad < intpoints.IP().nquad; ++iquad)
    {
      // evaluate values of shape functions and domain integration factor at current integration
      // point
      const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, iquad);

      // compute internal variables at current integration point
      VarManager()->SetInternalVariablesElchElectrodeSOCAndCRate(
          my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_);

      // compute velocity and its divergence
      static LINALG::Matrix<my::nsd_, 1> v;
      v.Multiply(my::evelnp_, my::funct_);
      double divv(0.);
      my::GetDivergence(divv, my::evelnp_);

      // integral of velocity divergence
      const double divv_fac = divv * fac;
      intdivv += divv_fac;

      // integral of concentration times velocity divergence
      intcdivv += my::scatravarmanager_->Phinp(0) * divv_fac;

      // integral of velocity times concentration gradient
      intvgradc += v.Dot(my::scatravarmanager_->GradPhi(0)) * fac;
    }  // loop over integration points

    // safety check
    if (scalars.Length() != 6)
      dserror(
          "Result vector for electrode state of charge and C rate computation has invalid length!");

    // write results for integrals related to velocity divergence into result vector
    scalars(3) = intdivv;
    scalars(4) = intcdivv;
    scalars(5) = intvgradc;
  }
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateElectrodeSOCAndCRate


/*---------------------------------------------------------------------*
 | calculate weighted mass flux (no reactive flux so far)   fang 02/15 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalculateFlux(
    LINALG::Matrix<my::nsd_, 1>& q,          //!< flux of species k
    const INPAR::SCATRA::FluxType fluxtype,  //!< type fo flux
    const int k                              //!< index of current scalar
)
{
  /*
    Actually, we compute here a weighted (and integrated) form of the fluxes!
    On time integration level, these contributions are then used to calculate
    an L2-projected representation of fluxes.
    Thus, this method here DOES NOT YET provide flux values that are ready to use!
  */

  // add convective flux contribution
  switch (fluxtype)
  {
    case INPAR::SCATRA::flux_diffusive:
    case INPAR::SCATRA::flux_total:
    {
      // diffusive flux contribution
      q.Update(-DiffManager()->GetIsotropicDiff(k), VarManager()->GradPhi(k));
      break;
    }

    default:
    {
      dserror("received illegal flag inside flux evaluation for whole domain");
      break;
    }
  }
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalculateFlux


/*----------------------------------------------------------------------------------------*
 | calculate error of numerical solution with respect to analytical solution   fang 10/16 |
 *----------------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalErrorComparedToAnalytSolution(
    const DRT::Element* ele,          //!< element
    Teuchos::ParameterList& params,   //!< parameter list
    Epetra_SerialDenseVector& errors  //!< vector containing L2 and H1 error norms
)
{
  // call base class routine
  myelch::CalErrorComparedToAnalytSolution(ele, params, errors);
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::CalErrorComparedToAnalytSolution


/*------------------------------------------------------------------------------*
 | set internal variables for electrodes                             fang 02/15 |
 *------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::SetInternalVariablesForMatAndRHS()
{
  // set internal variables
  VarManager()->SetInternalVariablesElchElectrode(
      my::funct_, my::derxy_, my::ephinp_, my::ephin_, my::econvelnp_, my::ehist_);
}  // DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype>::SetInternalVariablesForMatAndRHS()

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
void DRT::ELEMENTS::ScaTraEleCalcElchElectrode<distype, probdim>::CalcScaTraScaTraManifoldFlux(
    const DRT::Element* ele, Teuchos::ParameterList& params, DRT::Discretization& discretization,
    DRT::Element::LocationArray& la, Epetra_SerialDenseMatrix& eslavematrix,
    Epetra_SerialDenseMatrix& emastermatrix, Epetra_SerialDenseVector& eslaveresidual)
{
  const int nsd = static_cast<int>(my::nsd_);
  const int nen = static_cast<int>(my::nen_);

  // get states from this field and respective coupled field
  auto manifoldfield = discretization.GetState("phinp");
  if (manifoldfield == Teuchos::null) dserror("Cannot get state vector phinp!");
  auto coupledfield = elchmanifoldparams_->EvaluateMasterSide()
                          ? discretization.GetState(2, "imasterscatra")
                          : discretization.GetState(2, "scalarfield");
  if (coupledfield == Teuchos::null) dserror("Cannot get state vector scalarfield!");

  const int differentiationtype =
      params.get<int>("differentiationtype", static_cast<int>(SCATRA::DifferentiationType::none));

  // extract local nodal values from global state vectors
  std::vector<LINALG::Matrix<nen, 1>> ecoupledfield(
      my::numdofpernode_, LINALG::Matrix<nen, 1>(true));
  DRT::UTILS::ExtractMyValues(*manifoldfield, my::ephinp_, la[0].lm_);
  DRT::UTILS::ExtractMyValues(*coupledfield, ecoupledfield, la[2].lm_);

  const int kinetic_model = elchmanifoldparams_->KineticModel();

  // integration points and weights
  const DRT::UTILS::IntPointsAndWeights<my::nsd_ele_> intpoints(
      SCATRA::DisTypeToOptGaussRule<distype>::rule);

  for (int gpid = 0; gpid < intpoints.IP().nquad; ++gpid)
  {
    // evaluate values of shape functions and domain integration factor at current integration point
    const double fac = my::EvalShapeFuncAndDerivsAtIntPoint(intpoints, gpid);

    // evaluate overall integration factors
    const double timefacfac = my::scatraparatimint_->TimeFac() * fac;
    const double timefacrhsfac = my::scatraparatimint_->TimeFacRhs() * fac;

#ifdef DEBUG
    if (timefacfac < 0.0 or timefacrhsfac < 0.0) dserror("Integration factor is negative!");
#endif

    const double manifoldpotint = my::funct_.Dot(my::ephinp_[1]);
    const double coupledpotint = my::funct_.Dot(ecoupledfield[1]);

    switch (kinetic_model)
    {
      case INPAR::SSI::kinetics_constantinterfaceresistance:
      {
        const double resistance = elchmanifoldparams_->Resistance();
        const int num_electrons = elchmanifoldparams_->NumElectrons();

        const double inv_fluxresistance = 1.0 / resistance;

        const double j = (manifoldpotint - coupledpotint) * inv_fluxresistance;

        switch (differentiationtype)
        {
          case static_cast<int>(SCATRA::DifferentiationType::elch):
          {
            // residual
            const double jtimefacrhsfac = timefacrhsfac * j;

            // linearizations
            const double dj_dpot_slave_timefacfac = timefacfac * inv_fluxresistance;
            const double dj_dpot_master_timefacfac = -dj_dpot_slave_timefacfac;

            for (int vi = 0; vi < nen; ++vi)
            {
              const int row_conc = vi * 2;
              const int row_pot = vi * 2 + 1;

              for (int ui = 0; ui < nen; ++ui)
              {
                const int col_conc = ui * 2;
                const int col_pot = ui * 2 + 1;

                eslavematrix(row_conc, col_conc) +=
                    my::funct_(vi) * dj_dpot_slave_timefacfac * my::funct_(ui);
                eslavematrix(row_pot, col_pot) +=
                    my::funct_(vi) * num_electrons * dj_dpot_slave_timefacfac * my::funct_(ui);
              }

              const int nen_master = nen;
              for (int ui = 0; ui < nen_master; ++ui)
              {
                const int col_conc = ui * 2;
                const int col_pot = ui * 2 + 1;

                emastermatrix(row_conc, col_conc) +=
                    my::funct_(vi) * dj_dpot_master_timefacfac * my::funct_(ui);
                emastermatrix(row_pot, col_pot) +=
                    my::funct_(vi) * num_electrons * dj_dpot_master_timefacfac * my::funct_(ui);
              }

              eslaveresidual[row_conc] -= my::funct_(vi) * jtimefacrhsfac;
              eslaveresidual[row_pot] -= my::funct_(vi) * num_electrons * jtimefacrhsfac;
            }
            break;
          }
          case static_cast<int>(SCATRA::DifferentiationType::disp):
          {
            static LINALG::Matrix<1, nsd * nen> dJ_dmesh(false);
            const double J = my::xjm_.Determinant();
            for (unsigned node = 0; node < nen; node++)
              for (int dim = 0; dim < nsd; dim++)
                dJ_dmesh(dim + node * nsd) = J * my::derxy_(dim, node);

            const double timefacwgt = my::scatraparatimint_->TimeFac() * intpoints.IP().qwgt[gpid];
            const double dj_dd_slave_timefacwgt = timefacwgt * j;

            // loop over matrix rows
            for (int vi = 0; vi < nen; ++vi)
            {
              const int row_conc = vi * my::numdofpernode_;
              const int row_pot = row_conc + 1;
              const double vi_dj_dd_slave = my::funct_(vi) * dj_dd_slave_timefacwgt;

              // loop over matrix columns
              for (int ui = 0; ui < nen; ++ui)
              {
                // loop over spatial dimensions
                for (int dim = 0; dim < nsd; ++dim)
                {
                  const int col_disp = ui * nsd + dim;

                  // compute linearizations w.r.t. slave-side structural displacements
                  eslavematrix(row_conc, col_disp) += vi_dj_dd_slave * dJ_dmesh(col_disp);
                  eslavematrix(row_pot, col_disp) +=
                      num_electrons * vi_dj_dd_slave * dJ_dmesh(col_disp);
                }
              }
            }
            break;
          }
          default:
          {
            dserror("Unknown differentiation type");
            break;
          }
        }
        break;
      }
      case INPAR::SSI::kinetics_noflux:
      {
        // do nothing
        break;
      }
      default:
      {
        dserror("kinetic model not implemented.");
        break;
      }
    }
  }
}

// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcElchElectrode<DRT::Element::nurbs27>;