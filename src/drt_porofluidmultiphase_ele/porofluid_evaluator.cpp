/*----------------------------------------------------------------------*/
/*! \file
 \brief helper class for evaluation of the governing equation of multiphase porou flow

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/


#include "porofluid_evaluator.H"

#include "porofluid_variablemanager.H"
#include "porofluid_phasemanager.H"

#include "porofluidmultiphase_ele_parameter.H"

#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"

// necessary for function
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | factory method                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
Teuchos::RCP<DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<nsd, nen>>
DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<nsd, nen>::CreateEvaluator(
    const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
    const POROFLUIDMULTIPHASE::Action& action, int numdofpernode, int numfluidphases,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager)
{
  // the evaluator
  Teuchos::RCP<EvaluatorInterface<nsd, nen>> evaluator = Teuchos::null;

  bool inittimederiv = false;
  if (action == POROFLUIDMULTIPHASE::calc_initial_time_deriv) inittimederiv = true;

  // check if we also have to evaluate additional volume fraction terms
  const bool hasvolfracs = (numdofpernode - numfluidphases > 0);

  // determine action
  switch (action)
  {
    // calculate true pressures and saturation
    case POROFLUIDMULTIPHASE::calc_initial_time_deriv:
    case POROFLUIDMULTIPHASE::calc_mat_and_rhs:
    case POROFLUIDMULTIPHASE::calc_fluid_struct_coupl_mat:
    case POROFLUIDMULTIPHASE::calc_fluid_scatra_coupl_mat:
    {
      // initialize the evaluator for the multi phase element
      Teuchos::RCP<MultiEvaluator<nsd, nen>> evaluator_multiphase =
          Teuchos::rcp(new MultiEvaluator<nsd, nen>());

      // build evaluators for all but last fluid phase
      for (int curphase = 0; curphase < numfluidphases - 1; curphase++)
      {
        // initialize the evaluator for the current phase
        Teuchos::RCP<MultiEvaluator<nsd, nen>> evaluator_phase =
            Teuchos::rcp(new MultiEvaluator<nsd, nen>());

        // temporary interfaces
        Teuchos::RCP<EvaluatorInterface<nsd, nen>> tmpevaluator = Teuchos::null;
        Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

        // Note: this term cancels because of the formulation w.r.t. the material formulation of the
        // solid add evaluator for the conservative term (w, v \nabla \cdot S )
        // assembler = Teuchos::rcp(new AssembleStandard(curphase,inittimederiv));
        // tmpevaluator = Teuchos::rcp(new EvaluatorConv<nsd, nen>(assembler,curphase));
        // evaluator_phase->AddEvaluator(tmpevaluator);

        // add evaluator for the convective conservative term (w, S \nabla \cdot v )
        if (para.IsAle())
        {
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorSatDivVel<nsd, nen>(assembler, curphase));
          evaluator_phase->AddEvaluator(tmpevaluator);
        }
        // add evaluator for Biot stabilization
        if (para.BiotStab())
        {
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorBiotStab<nsd, nen>(assembler, curphase));
          evaluator_phase->AddEvaluator(tmpevaluator);
        }

        // add evaluator for the diffusive term (\nabla w, K \nabla p)
        // the diffusive term is also assembled into the last phase
        assembler = Teuchos::rcp(
            new AssembleAlsoIntoOtherPhase(curphase, numfluidphases - 1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorDiff<nsd, nen>(assembler, curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);

        // add evaluator for the reactive term
        if (phasemanager.IsReactive(curphase))
        {
          // the reactive term is also assembled into the last phase
          assembler = Teuchos::rcp(
              new AssembleAlsoIntoOtherPhase(curphase, numfluidphases - 1, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorReac<nsd, nen>(assembler, curphase));
          evaluator_phase->AddEvaluator(tmpevaluator);
        }

        // add evaluators for the instationary terms
        if (not para.IsStationary())
        {
          // add evaluator for the instationary pressure term
          // the term is also assembled into the last phase
          assembler = Teuchos::rcp(
              new AssembleAlsoIntoOtherPhase(curphase, numfluidphases - 1, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorMassPressure<nsd, nen>(assembler, curphase));
          evaluator_phase->AddEvaluator(tmpevaluator);

          // add evaluator for the instationary solid pressure term
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator =
              Teuchos::rcp(new EvaluatorMassSolidPressureSat<nsd, nen>(assembler, curphase));
          evaluator_phase->AddEvaluator(tmpevaluator);

          // add evaluator for the instationary saturation term
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorMassSaturation<nsd, nen>(assembler, curphase));
          evaluator_phase->AddEvaluator(tmpevaluator);
        }

        // add evaluators for the additional terms in fluid equations introduced by volume fractions
        if (hasvolfracs)
        {
          // add evaluators for the instationary terms
          if (not para.IsStationary())
          {
            // add evaluator for the instationary solid pressure term
            assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
            tmpevaluator =
                Teuchos::rcp(new EvaluatorVolFracAddInstatTermsSat<nsd, nen>(assembler, curphase));
            evaluator_phase->AddEvaluator(tmpevaluator);
          }

          if (para.IsAle())
          {
            assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
            tmpevaluator =
                Teuchos::rcp(new EvaluatorVolFracAddDivVelTermSat<nsd, nen>(assembler, curphase));
            evaluator_phase->AddEvaluator(tmpevaluator);
          }
        }

        // add the evaluator of the phase to the multiphase evaluator
        evaluator_multiphase->AddEvaluator(evaluator_phase);
      }

      // build evaluators for the last fluid phase
      {
        const int curphase = numfluidphases - 1;

        // initialize the evaluator for the last phase
        Teuchos::RCP<MultiEvaluator<nsd, nen>> evaluator_lastphase =
            Teuchos::rcp(new MultiEvaluator<nsd, nen>());

        // temporary interfaces
        Teuchos::RCP<EvaluatorInterface<nsd, nen>> tmpevaluator = Teuchos::null;
        Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

        // add evaluator for the convective conservative term (w, \nabla \cdot v )
        if (para.IsAle())
        {
          assembler = Teuchos::rcp(new AssembleStandard(curphase, false));
          tmpevaluator = Teuchos::rcp(new EvaluatorDivVel<nsd, nen>(assembler, curphase));
          evaluator_lastphase->AddEvaluator(tmpevaluator);
        }
        // add evaluator for Biot stabilization
        if (para.BiotStab())
        {
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorBiotStab<nsd, nen>(assembler, curphase));
          evaluator_lastphase->AddEvaluator(tmpevaluator);
        }

        // add evaluator for the diffusive term (\nabla w, K \nabla p)
        assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorDiff<nsd, nen>(assembler, curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);

        // add evaluator for the reactive term
        if (phasemanager.IsReactive(curphase))
        {
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorReac<nsd, nen>(assembler, curphase));
          evaluator_lastphase->AddEvaluator(tmpevaluator);
        }

        // add evaluators for the instationary terms
        if (not para.IsStationary())
        {
          // add evaluator for the instationary pressure term
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorMassPressure<nsd, nen>(assembler, curphase));
          evaluator_lastphase->AddEvaluator(tmpevaluator);

          // add evaluator for the instationary solid pressure term
          assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
          tmpevaluator =
              Teuchos::rcp(new EvaluatorMassSolidPressure<nsd, nen>(assembler, curphase));
          evaluator_lastphase->AddEvaluator(tmpevaluator);
        }

        // add evaluators for the additional terms in fluid equations introduced by volume fractions
        if (hasvolfracs)
        {
          // add evaluators for the instationary terms
          if (not para.IsStationary())
          {
            // add evaluator for the instationary solid pressure term
            assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
            tmpevaluator =
                Teuchos::rcp(new EvaluatorVolFracAddInstatTerms<nsd, nen>(assembler, curphase));
            evaluator_lastphase->AddEvaluator(tmpevaluator);
          }

          if (para.IsAle())
          {
            assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
            tmpevaluator =
                Teuchos::rcp(new EvaluatorVolFracAddDivVelTerm<nsd, nen>(assembler, curphase));
            evaluator_lastphase->AddEvaluator(tmpevaluator);
          }
        }

        // add the evaluator of the phase to the multiphase evaluator
        evaluator_multiphase->AddEvaluator(evaluator_lastphase);
      }

      // evaluate the additional volume fraction terms in the volume fraction equations
      if (hasvolfracs)
      {
        // initialize the evaluator for the volume fractions
        Teuchos::RCP<MultiEvaluator<nsd, nen>> evaluator_volfrac =
            Teuchos::rcp(new MultiEvaluator<nsd, nen>());

        // temporary interfaces
        Teuchos::RCP<EvaluatorInterface<nsd, nen>> tmpevaluator = Teuchos::null;
        Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

        // 1) volume fraction terms
        // ----------------------------------------------------------------- add evaluators for the
        // instationary terms
        if (not para.IsStationary())
        {
          assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorVolFracInstat<nsd, nen>(assembler, -1));
          evaluator_volfrac->AddEvaluator(tmpevaluator);
        }

        // add evaluators for the mesh-divergence term
        if (para.IsAle())
        {
          assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
          tmpevaluator = Teuchos::rcp(new EvaluatorVolFracDivVel<nsd, nen>(assembler, -1));
          evaluator_volfrac->AddEvaluator(tmpevaluator);
        }

        // diffusive term
        assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorVolFracDiff<nsd, nen>(assembler, -1));
        evaluator_volfrac->AddEvaluator(tmpevaluator);

        // reactive term
        assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorVolFracReac<nsd, nen>(assembler, -1));
        evaluator_volfrac->AddEvaluator(tmpevaluator);

        // additional flux term
        assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorVolFracAddFlux<nsd, nen>(assembler, -1));
        evaluator_volfrac->AddEvaluator(tmpevaluator);

        // 2) volume fraction pressure terms
        // -------------------------------------------------------- diffusive term
        assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorVolFracPressureDiff<nsd, nen>(assembler, -1));
        evaluator_volfrac->AddEvaluator(tmpevaluator);

        // reactive term
        assembler = Teuchos::rcp(new AssembleStandard(-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorVolFracPressureReac<nsd, nen>(assembler, -1));
        evaluator_volfrac->AddEvaluator(tmpevaluator);

        // add the evaluator of the volfractions to the multiphase evaluator
        evaluator_multiphase->AddEvaluator(evaluator_volfrac);
      }

      evaluator = evaluator_multiphase;
      break;
    }
    case POROFLUIDMULTIPHASE::calc_pres_and_sat:
    {
      // initialize the evaluator for the multi phase element
      Teuchos::RCP<MultiEvaluator<nsd, nen>> evaluator_multiphase =
          Teuchos::rcp(new MultiEvaluator<nsd, nen>());

      // temporary interfaces
      Teuchos::RCP<EvaluatorInterface<nsd, nen>> tmpevaluator = Teuchos::null;

      // initialize temporary assembler
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

      // build evaluators for all phases (fluid and volfrac)
      // volfrac does not actually need pressures and saturations --> set to -1 in evaluator
      for (int iphase = 0; iphase < numdofpernode; iphase++)
      {
        assembler = Teuchos::rcp(new AssembleStandard(iphase, false));
        tmpevaluator =
            Teuchos::rcp(new EvaluatorPressureAndSaturation<nsd, nen>(assembler, iphase));
        evaluator_multiphase->AddEvaluator(tmpevaluator);
      }
      evaluator = evaluator_multiphase;

      break;
    }
    case POROFLUIDMULTIPHASE::calc_solidpressure:
    {
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1, false));
      evaluator = Teuchos::rcp(new EvaluatorSolidPressure<nsd, nen>(assembler, -1));

      break;
    }
    case POROFLUIDMULTIPHASE::calc_porosity:
    {
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1, false));
      evaluator = Teuchos::rcp(new EvaluatorPorosity<nsd, nen>(assembler, -1));

      break;
    }
    case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
    {
      // initialize the evaluator for the multi phase element
      Teuchos::RCP<MultiEvaluator<nsd, nen>> evaluator_multiphase =
          Teuchos::rcp(new MultiEvaluator<nsd, nen>());

      // temporary interfaces
      Teuchos::RCP<EvaluatorInterface<nsd, nen>> tmpevaluator = Teuchos::null;

      // initialize temporary assembler
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

      assembler = Teuchos::rcp(new AssembleStandard(-1, false));
      tmpevaluator = Teuchos::rcp(new ReconstructFluxLinearization<nsd, nen>(assembler, -1));
      evaluator_multiphase->AddEvaluator(tmpevaluator);

      // build evaluators for all fluid phases
      for (int iphase = 0; iphase < numfluidphases; iphase++)
      {
        assembler = Teuchos::rcp(new AssembleStandard(iphase, false));
        tmpevaluator = Teuchos::rcp(new ReconstructFluxRHS<nsd, nen>(assembler, iphase));
        evaluator_multiphase->AddEvaluator(tmpevaluator);
      }
      evaluator = evaluator_multiphase;

      break;
    }
    case POROFLUIDMULTIPHASE::calc_valid_dofs:
    {
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1, false));
      evaluator = Teuchos::rcp(new EvaluatorValidVolFracPressures<nsd, nen>(assembler, -1));

      break;
    }
    case POROFLUIDMULTIPHASE::calc_domain_integrals:
    {
      int numscal = 0;
      if (para.HasScalar()) numscal = phasemanager.NumScal();
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1, false));
      evaluator = Teuchos::rcp(new EvaluatorDomainIntegrals<nsd, nen>(
          assembler, -1, para.DomainIntFunctions(), numscal));
      break;
    }
    default:
    {
      dserror("unknown action for evaluation class!");
      break;
    }
  }  // switch(action)

  // done
  return evaluator;
}

/*-----------------------------------------------------------------------------------*
 | linearization of a term scaled with saturation after fluid dofs  kremheller 03/18 |
 *-----------------------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(
    Epetra_SerialDenseMatrix& mymat, const LINALG::Matrix<nen, 1>& funct, const double prefac,
    const int numdofpernode, const int numfluidphases, const int curphase, const int phasetoadd,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager)
{
  for (int vi = 0; vi < nen; ++vi)
  {
    const double v = prefac * funct(vi);
    const int fvi = vi * numdofpernode + phasetoadd;

    for (int ui = 0; ui < nen; ++ui)
    {
      const double vfunct = v * funct(ui);
      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        const int fui = ui * numdofpernode + idof;

        mymat(fvi, fui) += vfunct * phasemanager.SaturationDeriv(curphase, idof);
      }
    }
  }

  return;
}
/*-----------------------------------------------------------------------------------*
 | linearization of a term scaled with porosity after fluid dofs    kremheller 03/18 |
 *-----------------------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBase<nsd, nen>::PorosityLinearizationFluid(
    Epetra_SerialDenseMatrix& mymat, const LINALG::Matrix<nen, 1>& funct, const double prefac,
    const int numdofpernode, const int phasetoadd,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager)
{
  for (int vi = 0; vi < nen; ++vi)
  {
    const double v = prefac * funct(vi);
    const int fvi = vi * numdofpernode + phasetoadd;

    for (int ui = 0; ui < nen; ++ui)
    {
      const double vfunct = v * funct(ui);
      for (int idof = 0; idof < numdofpernode; ++idof)
      {
        const int fui = ui * numdofpernode + idof;

        mymat(fvi, fui) += vfunct * phasemanager.PorosityDeriv(idof);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure                 |
 | of a term scaled with div (v^s)                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBase<nsd, nen>::CalcDivVelODMesh(
    Epetra_SerialDenseMatrix& mymat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& deriv, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nsd>& xjm, const LINALG::Matrix<nsd, nsd>& gridvelderiv,
    const double timefacfac, const double fac, const double det, const int numdofpernode,
    const int phasetoadd)
{
  // d (div v_s)/d d_n+1 = derxy * 1.0/theta/dt * d_n+1
  // prefactor is fac since timefacfac/theta/dt = fac
  CalcLinFacODMesh(mymat, funct, derxy, fac, numdofpernode, phasetoadd);

  // shapederivatives see fluid_ele_calc_poro.cpp
  if (nsd == 3)
  {
    const double gridvelderiv_0_0 = gridvelderiv(0, 0);
    const double gridvelderiv_0_1 = gridvelderiv(0, 1);
    const double gridvelderiv_0_2 = gridvelderiv(0, 2);
    const double gridvelderiv_1_0 = gridvelderiv(1, 0);
    const double gridvelderiv_1_1 = gridvelderiv(1, 1);
    const double gridvelderiv_1_2 = gridvelderiv(1, 2);
    const double gridvelderiv_2_0 = gridvelderiv(2, 0);
    const double gridvelderiv_2_1 = gridvelderiv(2, 1);
    const double gridvelderiv_2_2 = gridvelderiv(2, 2);

    const double xjm_0_0 = xjm(0, 0);
    const double xjm_0_1 = xjm(0, 1);
    const double xjm_0_2 = xjm(0, 2);
    const double xjm_1_0 = xjm(1, 0);
    const double xjm_1_1 = xjm(1, 1);
    const double xjm_1_2 = xjm(1, 2);
    const double xjm_2_0 = xjm(2, 0);
    const double xjm_2_1 = xjm(2, 1);
    const double xjm_2_2 = xjm(2, 2);

#define derxjm_(r, c, d, i) derxjm_##r##c##d(i)

#define derxjm_001(ui) (deriv(2, ui) * xjm_1_2 - deriv(1, ui) * xjm_2_2)
#define derxjm_002(ui) (deriv(1, ui) * xjm_2_1 - deriv(2, ui) * xjm_1_1)

#define derxjm_100(ui) (deriv(1, ui) * xjm_2_2 - deriv(2, ui) * xjm_1_2)
#define derxjm_102(ui) (deriv(2, ui) * xjm_1_0 - deriv(1, ui) * xjm_2_0)

#define derxjm_200(ui) (deriv(2, ui) * xjm_1_1 - deriv(1, ui) * xjm_2_1)
#define derxjm_201(ui) (deriv(1, ui) * xjm_2_0 - deriv(2, ui) * xjm_1_0)

#define derxjm_011(ui) (deriv(0, ui) * xjm_2_2 - deriv(2, ui) * xjm_0_2)
#define derxjm_012(ui) (deriv(2, ui) * xjm_0_1 - deriv(0, ui) * xjm_2_1)

#define derxjm_110(ui) (deriv(2, ui) * xjm_0_2 - deriv(0, ui) * xjm_2_2)
#define derxjm_112(ui) (deriv(0, ui) * xjm_2_0 - deriv(2, ui) * xjm_0_0)

#define derxjm_210(ui) (deriv(0, ui) * xjm_2_1 - deriv(2, ui) * xjm_0_1)
#define derxjm_211(ui) (deriv(2, ui) * xjm_0_0 - deriv(0, ui) * xjm_2_0)

#define derxjm_021(ui) (deriv(1, ui) * xjm_0_2 - deriv(0, ui) * xjm_1_2)
#define derxjm_022(ui) (deriv(0, ui) * xjm_1_1 - deriv(1, ui) * xjm_0_1)

#define derxjm_120(ui) (deriv(0, ui) * xjm_1_2 - deriv(1, ui) * xjm_0_2)
#define derxjm_122(ui) (deriv(1, ui) * xjm_0_0 - deriv(0, ui) * xjm_1_0)

#define derxjm_220(ui) (deriv(1, ui) * xjm_0_1 - deriv(0, ui) * xjm_1_1)
#define derxjm_221(ui) (deriv(0, ui) * xjm_1_0 - deriv(1, ui) * xjm_0_0)

    for (int ui = 0; ui < nen; ++ui)
    {
      const double v0 =
          gridvelderiv_1_0 * derxjm_(0, 0, 1, ui) + gridvelderiv_1_1 * derxjm_(0, 1, 1, ui) +
          gridvelderiv_1_2 * derxjm_(0, 2, 1, ui) + gridvelderiv_2_0 * derxjm_(0, 0, 2, ui) +
          gridvelderiv_2_1 * derxjm_(0, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(0, 2, 2, ui);

      const double v1 =
          gridvelderiv_0_0 * derxjm_(1, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(1, 1, 0, ui) +
          gridvelderiv_0_2 * derxjm_(1, 2, 0, ui) + gridvelderiv_2_0 * derxjm_(1, 0, 2, ui) +
          gridvelderiv_2_1 * derxjm_(1, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(1, 2, 2, ui);

      const double v2 =
          gridvelderiv_0_0 * derxjm_(2, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(2, 1, 0, ui) +
          gridvelderiv_0_2 * derxjm_(2, 2, 0, ui) + gridvelderiv_1_0 * derxjm_(2, 0, 1, ui) +
          gridvelderiv_1_1 * derxjm_(2, 1, 1, ui) + gridvelderiv_1_2 * derxjm_(2, 2, 1, ui);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + phasetoadd;
        const double v = timefacfac / det * funct(vi);

        mymat(fvi, ui * 3 + 0) += v * v0;

        mymat(fvi, ui * 3 + 1) += v * v1;

        mymat(fvi, ui * 3 + 2) += v * v2;
      }
    }
  }
  else if (nsd == 2)
  {
    const double gridvelderiv_0_0 = gridvelderiv(0, 0);
    const double gridvelderiv_0_1 = gridvelderiv(0, 1);
    const double gridvelderiv_1_0 = gridvelderiv(1, 0);
    const double gridvelderiv_1_1 = gridvelderiv(1, 1);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + phasetoadd;
      const double v = timefacfac / det * funct(vi);
      for (int ui = 0; ui < nen; ++ui)
      {
        mymat(fvi, ui * 2) +=
            v * (-gridvelderiv_1_0 * deriv(1, ui) + gridvelderiv_1_1 * deriv(0, ui));

        mymat(fvi, ui * 2 + 1) +=
            v * (+gridvelderiv_0_0 * deriv(1, ui) - gridvelderiv_0_1 * deriv(0, ui));
      }
    }
  }
  else
    dserror("shapederivatives not implemented for 1D!");

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure                 |
 | (Fac = Jacobian determinant)                        kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBase<nsd, nen>::CalcLinFacODMesh(
    Epetra_SerialDenseMatrix& mymat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const double vrhs, const int numdofpernode,
    const int phasetoadd)
{
  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J
  //* N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space,
  // i.e. det(dx/ds) in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac *
  // N_x
  //              fac        = J              --> d(fac)/dd        = fac * N_x

  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;
    const double v = vrhs * funct(vi);

    for (int ui = 0; ui < nen; ++ui)
    {
      for (int idim = 0; idim < nsd; ++idim)
      {
        const int fui = ui * nsd + idim;
        mymat(fvi, fui) += v * derxy(idim, ui);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure                 |
 | (diffusive term)                                    kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBase<nsd, nen>::CalcDiffODMesh(
    Epetra_SerialDenseMatrix& mymat, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm,
    const LINALG::Matrix<nsd, 1>& diffflux, const LINALG::Matrix<nsd, 1>& refgrad,
    const LINALG::Matrix<nsd, 1>& grad, const double timefacfac, const double difffac,
    const int numdofpernode, const int phasetoadd)
{
  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J
  //* N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space,
  // i.e. det(dx/ds) in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac *
  // N_x
  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;

    // laplacian in weak form
    double laplawf(0.0);
    for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j);
    double v = -laplawf * timefacfac;

    for (int ui = 0; ui < nen; ++ui)
    {
      for (int idim = 0; idim < nsd; ++idim)
      {
        const int fui = ui * nsd + idim;
        mymat(fvi, fui) += v * derxy(idim, ui);
      }
    }
  }

  //----------------------------------------------------------------
  // standard Galerkin terms  -- "shapederivatives" diffusive term
  //----------------------------------------------------------------
  // see scatra_ele_calc_OD.cpp

  if (nsd == 3)
  {
    const double xjm_0_0 = xjm(0, 0);
    const double xjm_0_1 = xjm(0, 1);
    const double xjm_0_2 = xjm(0, 2);
    const double xjm_1_0 = xjm(1, 0);
    const double xjm_1_1 = xjm(1, 1);
    const double xjm_1_2 = xjm(1, 2);
    const double xjm_2_0 = xjm(2, 0);
    const double xjm_2_1 = xjm(2, 1);
    const double xjm_2_2 = xjm(2, 2);

    const double grad_0 = grad(0);
    const double grad_1 = grad(1);
    const double grad_2 = grad(2);

    for (unsigned vi = 0; vi < nen; ++vi)
    {
      const double deriv_vi_0 = deriv(0, vi);
      const double deriv_vi_1 = deriv(1, vi);
      const double deriv_vi_2 = deriv(2, vi);

      const int fvi = vi * numdofpernode + phasetoadd;

      for (unsigned ui = 0; ui < nen; ++ui)
      {
        const double v00 =
            +grad_1 * (deriv_vi_0 * (deriv(2, ui) * xjm_1_2 - deriv(1, ui) * xjm_2_2) +
                          deriv_vi_1 * (deriv(0, ui) * xjm_2_2 - deriv(2, ui) * xjm_0_2) +
                          deriv_vi_2 * (deriv(1, ui) * xjm_0_2 - deriv(0, ui) * xjm_1_2)) +
            grad_2 * (deriv_vi_0 * (deriv(1, ui) * xjm_2_1 - deriv(2, ui) * xjm_1_1) +
                         deriv_vi_1 * (deriv(2, ui) * xjm_0_1 - deriv(0, ui) * xjm_2_1) +
                         deriv_vi_2 * (deriv(0, ui) * xjm_1_1 - deriv(1, ui) * xjm_0_1));
        const double v01 =
            +grad_0 * (deriv_vi_0 * (deriv(1, ui) * xjm_2_2 - deriv(2, ui) * xjm_1_2) +
                          deriv_vi_1 * (deriv(2, ui) * xjm_0_2 - deriv(0, ui) * xjm_2_2) +
                          deriv_vi_2 * (deriv(0, ui) * xjm_1_2 - deriv(1, ui) * xjm_0_2)) +
            grad_2 * (deriv_vi_0 * (deriv(2, ui) * xjm_1_0 - deriv(1, ui) * xjm_2_0) +
                         deriv_vi_1 * (deriv(0, ui) * xjm_2_0 - deriv(2, ui) * xjm_0_0) +
                         deriv_vi_2 * (deriv(1, ui) * xjm_0_0 - deriv(0, ui) * xjm_1_0));
        const double v02 =
            +grad_0 * (deriv_vi_0 * (deriv(2, ui) * xjm_1_1 - deriv(1, ui) * xjm_2_1) +
                          deriv_vi_1 * (deriv(0, ui) * xjm_2_1 - deriv(2, ui) * xjm_0_1) +
                          deriv_vi_2 * (deriv(1, ui) * xjm_0_1 - deriv(0, ui) * xjm_1_1)) +
            grad_1 * (deriv_vi_0 * (deriv(1, ui) * xjm_2_0 - deriv(2, ui) * xjm_1_0) +
                         deriv_vi_1 * (deriv(2, ui) * xjm_0_0 - deriv(0, ui) * xjm_2_0) +
                         deriv_vi_2 * (deriv(0, ui) * xjm_1_0 - deriv(1, ui) * xjm_0_0));

        mymat(fvi, ui * nsd + 0) += difffac * v00;
        mymat(fvi, ui * nsd + 1) += difffac * v01;
        mymat(fvi, ui * nsd + 2) += difffac * v02;
      }
    }

    const double refgrad_0 = refgrad(0);
    const double refgrad_1 = refgrad(1);
    const double refgrad_2 = refgrad(2);

    for (unsigned vi = 0; vi < nen; ++vi)
    {
      const double derxy_vi_0 = derxy(0, vi);
      const double derxy_vi_1 = derxy(1, vi);
      const double derxy_vi_2 = derxy(2, vi);

      const int fvi = vi * numdofpernode + phasetoadd;

      for (unsigned ui = 0; ui < nen; ++ui)
      {
        const double v00 =
            +derxy_vi_1 * (refgrad_0 * (deriv(2, ui) * xjm_1_2 - deriv(1, ui) * xjm_2_2) +
                              refgrad_1 * (deriv(0, ui) * xjm_2_2 - deriv(2, ui) * xjm_0_2) +
                              refgrad_2 * (deriv(1, ui) * xjm_0_2 - deriv(0, ui) * xjm_1_2)) +
            derxy_vi_2 * (refgrad_0 * (deriv(1, ui) * xjm_2_1 - deriv(2, ui) * xjm_1_1) +
                             refgrad_1 * (deriv(2, ui) * xjm_0_1 - deriv(0, ui) * xjm_2_1) +
                             refgrad_2 * (deriv(0, ui) * xjm_1_1 - deriv(1, ui) * xjm_0_1));
        const double v01 =
            +derxy_vi_0 * (refgrad_0 * (deriv(1, ui) * xjm_2_2 - deriv(2, ui) * xjm_1_2) +
                              refgrad_1 * (deriv(2, ui) * xjm_0_2 - deriv(0, ui) * xjm_2_2) +
                              refgrad_2 * (deriv(0, ui) * xjm_1_2 - deriv(1, ui) * xjm_0_2)) +
            derxy_vi_2 * (refgrad_0 * (deriv(2, ui) * xjm_1_0 - deriv(1, ui) * xjm_2_0) +
                             refgrad_1 * (deriv(0, ui) * xjm_2_0 - deriv(2, ui) * xjm_0_0) +
                             refgrad_2 * (deriv(1, ui) * xjm_0_0 - deriv(0, ui) * xjm_1_0));
        const double v02 =
            +derxy_vi_0 * (refgrad_0 * (deriv(2, ui) * xjm_1_1 - deriv(1, ui) * xjm_2_1) +
                              refgrad_1 * (deriv(0, ui) * xjm_2_1 - deriv(2, ui) * xjm_0_1) +
                              refgrad_2 * (deriv(1, ui) * xjm_0_1 - deriv(0, ui) * xjm_1_1)) +
            derxy_vi_1 * (refgrad_0 * (deriv(1, ui) * xjm_2_0 - deriv(2, ui) * xjm_1_0) +
                             refgrad_1 * (deriv(2, ui) * xjm_0_0 - deriv(0, ui) * xjm_2_0) +
                             refgrad_2 * (deriv(0, ui) * xjm_1_0 - deriv(1, ui) * xjm_0_0));

        mymat(fvi, ui * nsd + 0) += difffac * v00;
        mymat(fvi, ui * nsd + 1) += difffac * v01;
        mymat(fvi, ui * nsd + 2) += difffac * v02;
      }
    }
  }
  else if (nsd == 2)
  {
    {
      const double grad_0 = grad(0);
      const double grad_1 = grad(1);

      for (int vi = 0; vi < nen; ++vi)
      {
        const double deriv_vi_0 = deriv(0, vi);
        const double deriv_vi_1 = deriv(1, vi);

        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double v00 = +grad_1 * (-deriv_vi_0 * deriv(1, ui) + deriv_vi_1 * deriv(0, ui));
          const double v01 = +grad_0 * (deriv_vi_0 * deriv(1, ui) - deriv_vi_1 * deriv(0, ui));

          mymat(fvi, ui * nsd + 0) += difffac * v00;
          mymat(fvi, ui * nsd + 1) += difffac * v01;
        }
      }
    }

    const double refgrad_0 = refgrad(0);
    const double refgrad_1 = refgrad(1);

    for (unsigned vi = 0; vi < nen; ++vi)
    {
      const double derxy_vi_0 = derxy(0, vi);
      const double derxy_vi_1 = derxy(1, vi);

      const int fvi = vi * numdofpernode + phasetoadd;

      for (unsigned ui = 0; ui < nen; ++ui)
      {
        const double v00 = +derxy_vi_1 * (-refgrad_0 * deriv(1, ui) + refgrad_1 * deriv(0, ui));
        const double v01 = +derxy_vi_0 * (refgrad_0 * deriv(1, ui) - refgrad_1 * deriv(0, ui));

        mymat(fvi, ui * nsd + 0) += difffac * v00;
        mymat(fvi, ui * nsd + 1) += difffac * v01;
      }
    }
  }
  else
    dserror("shapederivatives not implemented for 1D!");

  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();

  // convective term in convective form
  /*
       /                               \
      |                                 |
      | prefac * v * nabla * Dphi  , q  |
      |                                 |
       \                               /
  */
  const double prefac = timefacfac;

  static LINALG::Matrix<nen, 1> conv;
  // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
  conv.MultiplyTN(derxy, *variablemanager.ConVelnp());

  for (int vi = 0; vi < nen; ++vi)
  {
    const double v = prefac * funct(vi);
    const int fvi = vi * numdofpernode + phasetoadd;

    for (int ui = 0; ui < nen; ++ui)
      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        const int fui = ui * numdofpernode + idof;

        mymat(fvi, fui) += v * conv(ui) * phasemanager.SaturationDeriv(curphase, idof);
      }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  double conv_sat = 0.0;
  for (int idof = 0; idof < numfluidphases; ++idof)
  {
    // convective term
    const double conv_phi = variablemanager.ConVelnp()->Dot((*variablemanager.GradPhinp())[idof]);
    conv_sat += rhsfac * phasemanager.SaturationDeriv(curphase, idof) * conv_phi;
  }
  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;

    myvec[fvi] -= conv_sat * funct(vi);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd, nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& deriv, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nsd>& xjm, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd, nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  double vrhs = rhsfac * variablemanager.DivConVelnp();

  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;

    myvec[fvi] -= vrhs * funct(vi);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  static LINALG::Matrix<nsd, nsd> gridvelderiv(true);
  gridvelderiv.MultiplyNT(*(variablemanager.EConVelnp()), deriv);

  // OD mesh - div vel term
  EvaluatorBase<nsd, nen>::CalcDivVelODMesh(mymat, funct, deriv, derxy, xjm, gridvelderiv,
      timefacfac, fac, det, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // call base class
  EvaluatorDivVel<nsd, nen>::EvaluateMatrixAndAssemble(elemat, funct, derxy, curphase, phasetoadd,
      numdofpernode, phasemanager, variablemanager, timefacfac, fac, inittimederiv);

  // no linearization needed in case of initial time derivative calculation
  if (!inittimederiv)
  {
    const int numfluidphases = phasemanager.NumFluidPhases();
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const double consfac = timefacfac * variablemanager.DivConVelnp();

    // call base class for saturation linearization
    EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(
        mymat, funct, consfac, numdofpernode, numfluidphases, curphase, phasetoadd, phasemanager);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorDivVel<nsd, nen>::EvaluateVectorAndAssemble(elevec, funct, derxy, xyze, curphase,
      phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * rhsfac, phasemanager.Saturation(curphase) * fac,
      inittimederiv);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // call base class with scaled factors
  EvaluatorDivVel<nsd, nen>::EvaluateMatrixODStructAndAssemble(elemat, funct, deriv, derxy, xjm,
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * timefacfac, phasemanager.Saturation(curphase) * fac, det);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  dserror("Biot stabilization is still missing");
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  dserror("Biot stabilization is still missing");

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  dserror("Biot stabilization is still missing");
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const int numfluidphases = phasemanager.NumFluidPhases();

    const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

    // current pressure gradient
    static LINALG::Matrix<nsd, 1> gradpres(true);
    gradpres.Clear();

    // compute the pressure gradient from the phi gradients
    for (int idof = 0; idof < numfluidphases; ++idof)
      gradpres.Update(phasemanager.PressureDeriv(curphase, idof), gradphi[idof], 1.0);

    double abspressgrad = 0.0;
    for (int i = 0; i < nsd; i++) abspressgrad += gradpres(i) * gradpres(i);
    abspressgrad = sqrt(abspressgrad);

    // permeability tensor
    static LINALG::Matrix<nsd, nsd> permeabilitytensor(true);
    phasemanager.PermeabilityTensor(curphase, permeabilitytensor);

    static LINALG::Matrix<nsd, nen> diffflux(true);
    diffflux.Multiply(permeabilitytensor, derxy);
    diffflux.Scale(
        phasemanager.RelPermeability(curphase) / phasemanager.DynViscosity(curphase, abspressgrad));

    // helper variable for linearization
    static LINALG::Matrix<nsd, 1> diffflux_relpermeability(true);

    if (not phasemanager.HasConstantRelPermeability(curphase))
    {
      diffflux_relpermeability.Multiply(permeabilitytensor, gradpres);
      diffflux_relpermeability.Scale(phasemanager.RelPermeabilityDeriv(curphase) /
                                     phasemanager.DynViscosity(curphase, abspressgrad));
    }
    else
      diffflux_relpermeability.Scale(0.0);

    //----------------------------------------------------------------
    // diffusive term and linearization of relative permeability w.r.t. dof
    //----------------------------------------------------------------
    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + phasetoadd;
      double laplawf_relpermeability(0.0);

      // helper variable for linearization
      for (int j = 0; j < nsd; j++)
        laplawf_relpermeability += derxy(j, vi) * diffflux_relpermeability(j);

      for (int ui = 0; ui < nen; ++ui)
      {
        double laplawf(0.0);
        for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j, ui);

        for (int idof = 0; idof < numfluidphases; ++idof)
        {
          const int fui = ui * numdofpernode + idof;
          mymat(fvi, fui) += timefacfac * (laplawf * phasemanager.PressureDeriv(curphase, idof) +
                                              funct(ui) * laplawf_relpermeability *
                                                  phasemanager.SaturationDeriv(curphase, idof));
        }
      }
    }

    //----------------------------------------------------------------
    // linearization of dynamic viscosity w.r.t. dof
    //----------------------------------------------------------------
    if (not phasemanager.HasConstantDynViscosity(curphase))
    {
      // derivative of abspressgrad w.r.t. pressure gradient
      static LINALG::Matrix<nsd, 1> dabspressgraddpresgradp(true);
      dabspressgraddpresgradp.Scale(0.0);
      // avoid division by zero
      if (abspressgrad > 1.0e-12)
        for (int i = 0; i < nsd; i++) dabspressgraddpresgradp(i) = gradpres(i) / abspressgrad;

      static LINALG::Matrix<nsd, 1> diffflux2(true);
      diffflux2.Multiply(permeabilitytensor, gradpres);
      // d (1/visc) / d abspressgrad = -1.0 * visc^(-2) * d visc / d abspressgrad
      diffflux2.Scale(-1.0 * phasemanager.RelPermeability(curphase) /
                      phasemanager.DynViscosity(curphase, abspressgrad) /
                      phasemanager.DynViscosity(curphase, abspressgrad) *
                      phasemanager.DynViscosityDeriv(curphase, abspressgrad));

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + phasetoadd;
        double laplawf = 0.0;
        for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux2(j);
        for (int ui = 0; ui < nen; ++ui)
        {
          double gradpderxy = 0.0;
          for (int j = 0; j < nsd; j++) gradpderxy += derxy(j, ui) * dabspressgraddpresgradp(j);
          for (int idof = 0; idof < numfluidphases; ++idof)
          {
            const int fui = ui * numdofpernode + idof;
            // d abspressgrad / d phi = d abspressgrad / d gradp * d gradp / d phi =
            //                        = d abspressgrad / d gradp * d / d phi( d p / d phi * d phi /
            //                        d x) = = d abspressgrad / d gradp * derxy * d p / d phi
            // Note: FD-Check might fail here due to kink in formulation of cell-adherence model-law
            mymat(fvi, fui) +=
                timefacfac * laplawf * phasemanager.PressureDeriv(curphase, idof) * gradpderxy;
          }
        }
      }
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // current pressure gradient
  static LINALG::Matrix<nsd, 1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof = 0; idof < numfluidphases; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase, idof), gradphi[idof], 1.0);

  double abspressgrad = 0.0;
  for (int i = 0; i < nsd; i++) abspressgrad += gradpres(i) * gradpres(i);
  abspressgrad = sqrt(abspressgrad);

  // diffusion tensor
  static LINALG::Matrix<nsd, nsd> difftensor(true);
  phasemanager.PermeabilityTensor(curphase, difftensor);
  difftensor.Scale(
      phasemanager.RelPermeability(curphase) / phasemanager.DynViscosity(curphase, abspressgrad));

  static LINALG::Matrix<nsd, 1> diffflux(true);
  diffflux.Multiply(difftensor, gradpres);

  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;

    // laplacian in weak form
    double laplawf(0.0);
    for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j);
    myvec[fvi] -= rhsfac * laplawf;
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd, nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& deriv, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nsd>& xjm, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // current pressure gradient
  static LINALG::Matrix<nsd, 1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof = 0; idof < numfluidphases; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase, idof), gradphi[idof], 1.0);

  double abspressgrad = 0.0;
  for (int i = 0; i < nsd; i++) abspressgrad += gradpres(i) * gradpres(i);
  abspressgrad = sqrt(abspressgrad);

  // diffusion tensor
  static LINALG::Matrix<nsd, nsd> difftensor(true);
  phasemanager.PermeabilityTensor(curphase, difftensor);
  difftensor.Scale(
      phasemanager.RelPermeability(curphase) / phasemanager.DynViscosity(curphase, abspressgrad));

  // TODO: anisotropic difftensor and
  //       non-constant viscosity (because of pressure gradient, probably not really necessary)
  static LINALG::Matrix<nsd, 1> diffflux(true);
  diffflux.Multiply(difftensor, gradpres);

  // diffusive pre-factor for linearization
  const double v = difftensor(0, 0) * timefacfac / det;

  // gradient of pressure w.r.t. reference coordinates
  static LINALG::Matrix<nsd, 1> refgradpres(true);
  refgradpres.Clear();

  // gradient of phi w.r.t. reference coordinates
  std::vector<LINALG::Matrix<nsd, 1>> refgradphi(
      numfluidphases, LINALG::Matrix<nsd, 1>(true));  // static LINALG::Matrix<nsd,1> refgradphi;
  for (int idof = 0; idof < numfluidphases; ++idof) refgradphi[idof].Multiply(xjm, gradphi[idof]);

  // compute the pressure gradient from the phi gradients
  for (int idof = 0; idof < numfluidphases; ++idof)
    refgradpres.Update(phasemanager.PressureDeriv(curphase, idof), refgradphi[idof], 1.0);

  // OD mesh - diffusive term
  EvaluatorBase<nsd, nen>::CalcDiffODMesh(mymat, deriv, derxy, xjm, diffflux, refgradpres, gradpres,
      timefacfac, v, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd, nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    // TODO a constant density is assumed here
    double scaledtimefacfac = timefacfac / phasemanager.Density(curphase);

    //----------------------------------------------------------------
    // reaction terms
    //----------------------------------------------------------------
    for (int vi = 0; vi < nen; ++vi)
    {
      const double v = scaledtimefacfac * funct(vi);
      const int fvi = vi * numdofpernode + phasetoadd;

      for (int ui = 0; ui < nen; ++ui)
      {
        const double vfunct = v * funct(ui);
        for (int idof = 0; idof < numdofpernode; ++idof)
        {
          const int fui = ui * numdofpernode + idof;

          // rhs ---> -
          mymat(fvi, fui) -= vfunct * phasemanager.ReacDeriv(curphase, idof);
        }
      }
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  if (not phasemanager.IsReactive(curphase)) return;

  // TODO a constant density is assumed here
  double scale = 1.0 / phasemanager.Density(curphase);

  double vrhs = scale * rhsfac * phasemanager.ReacTerm(curphase);

  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;
    // rhs ---> +
    myvec[fvi] += vrhs * funct(vi);
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd, nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& deriv, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nsd>& xjm, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  if (not phasemanager.IsReactive(curphase)) return;

  // TODO a constant density is assumed here
  double scale = 1.0 / phasemanager.Density(curphase);
  double vrhs = scale * timefacfac * phasemanager.ReacTerm(curphase);

  // linearization of porosity (may appear in reaction term)
  //-------------dreac/dd = dreac/dporosity * dporosity/dd = dreac/dporosity * dporosity/dJ * dJ/dd
  //= dreac/dporosity * dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
  // (dx/ds) * ( det(dX/ds) )^-1

  if (phasemanager.PorosityDependsOnStruct())
  {
    vrhs += timefacfac * scale * phasemanager.ReacDerivPorosity(curphase) *
            phasemanager.JacobianDefGrad() * phasemanager.PorosityDerivWrtJacobianDefGrad();
  }

  // linearization of mesh motion (Jacobian)
  // 1) linearization of fac +
  // 2) possible linearization w.r.t porosity
  // rhs ---> -
  EvaluatorBase<nsd, nen>::CalcLinFacODMesh(
      mymat, funct, derxy, -1.0 * vrhs, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd, nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  if (not phasemanager.IsReactive(curphase)) return;

  const int numscal = phasemanager.NumScal();

  double vrhs = 1.0 / phasemanager.Density(curphase) * timefacfac;

  // linearization of reaction term w.r.t scalars
  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;
    const double v = vrhs * funct(vi);

    for (int ui = 0; ui < nen; ++ui)
    {
      const double vfunct = v * funct(ui);
      for (int iscal = 0; iscal < numscal; ++iscal)
      {
        const int fui = ui * numscal + iscal;
        // rhs ---> -
        mymat(fvi, fui) -= vfunct * phasemanager.ReacDerivScalar(curphase, iscal);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // in case of an incompressible fluid phase 'curphase' (1/K = 0) the term cancels out
  if (phasemanager.IncompressibleFluidPhase(curphase)) return;

  const int numfluidphases = phasemanager.NumFluidPhases();

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  // saturation
  const double saturation = phasemanager.Saturation(curphase);

  // inverse bulk modulus of phase (compressibility)
  const double invbulkmodulus = phasemanager.InvBulkmodulus(curphase);

  // pre factor
  const double facfacmass = fac * phasemanager.Porosity() * saturation * invbulkmodulus;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (int vi = 0; vi < nen; ++vi)
  {
    const double v = facfacmass * funct(vi);
    const int fvi = vi * numdofpernode + phasetoadd;

    for (int ui = 0; ui < nen; ++ui)
    {
      const double vfunct = v * funct(ui);
      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        const int fui = ui * numdofpernode + idof;

        mymat(fvi, fui) += vfunct * phasemanager.PressureDeriv(curphase, idof);
      }
    }
  }

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if (!inittimederiv)
  {
    //----------------------------------------------------------------
    // linearization of saturation w.r.t. dof
    //----------------------------------------------------------------
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      double hist = 0.0;
      // if(curphase==phasetoadd) // bug fix??
      hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass2 =
          fac * phasemanager.PressureDeriv(curphase, phasetoadd) * (phinp[phasetoadd] - hist);

      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        if (idof != phasetoadd)
          facfacmass2 += timefacfac * phasemanager.PressureDeriv(curphase, idof) *
                         (*variablemanager.Phidtnp())[idof];
      }

      facfacmass2 *= phasemanager.Porosity() * invbulkmodulus;

      // call base class for saturation linearization
      EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(mymat, funct, facfacmass2,
          numdofpernode, numfluidphases, curphase, phasetoadd, phasemanager);
    }

    //----------------------------------------------------------------
    // linearization of porosity w.r.t. dof
    //----------------------------------------------------------------
    if (phasemanager.PorosityDependsOnFluid())
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      double hist = 0.0;
      // if(curphase==phasetoadd)  //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass =
          fac * phasemanager.PressureDeriv(curphase, phasetoadd) * (phinp[phasetoadd] - hist);

      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        if (idof != phasetoadd)
          facfacmass += timefacfac * phasemanager.PressureDeriv(curphase, idof) *
                        (*variablemanager.Phidtnp())[idof];
      }

      facfacmass *= saturation * invbulkmodulus;

      // call base class:
      EvaluatorBase<nsd, nen>::PorosityLinearizationFluid(
          mymat, funct, facfacmass, numdofpernode, phasetoadd, phasemanager);
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // in case of an incompressible fluid phase 'curphase' (1/K = 0) the term cancels out
  if (phasemanager.IncompressibleFluidPhase(curphase)) return;

  // for the initial time derivative calculation no transient terms enter the rhs
  if (!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    double vtrans = GetRhsTrans(
        curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, rhsfac, fac);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + phasetoadd;

      myvec[fvi] -= vtrans * funct(vi);
    }
  }  // !inittimederiv
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // in case of an incompressible fluid phase 'curphase' (1/K = 0) the term cancels out
  if (phasemanager.IncompressibleFluidPhase(curphase)) return;

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  double vtrans = GetRhsTrans(
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, timefacfac, fac);

  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd =
  // dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
  // (dx/ds) * ( det(dX/ds) )^-1 in our case: vtrans is scaled with porosity in GetRhsTrans -->
  // scale it with 1.0/porosity here

  if (phasemanager.PorosityDependsOnStruct())
    vtrans += vtrans * 1.0 / phasemanager.Porosity() * phasemanager.JacobianDefGrad() *
              phasemanager.PorosityDerivWrtJacobianDefGrad();

  // linearization of mesh motion (Jacobian)
  // 1) linearization of fac +
  // 2) possible linearization w.r.t porosity
  EvaluatorBase<nsd, nen>::CalcLinFacODMesh(mymat, funct, derxy, vtrans, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate transient term at GP                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd, nen>::GetRhsTrans(int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac)
{
  const int numfluidphases = phasemanager.NumFluidPhases();
  // read data from managers
  double hist = 0.0;
  // if(curphase==phasetoadd) //bugfix??
  hist = (*variablemanager.Hist())[phasetoadd];
  // std::cout << "hist = " << hist << std::endl;
  const double porosity = phasemanager.Porosity();
  const std::vector<double>& phinp = *variablemanager.Phinp();
  const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

  // saturation
  const double saturation = phasemanager.Saturation(curphase);

  // inverse bulk modulus of phase (compressibility)
  const double invbulkmodulus = phasemanager.InvBulkmodulus(curphase);

  double vtrans = 0.0;

  // TODO check for Genalpha
  // compute scalar at integration point
  vtrans = fac * phasemanager.PressureDeriv(curphase, phasetoadd) * (phinp[phasetoadd] - hist);
  for (int idof = 0; idof < numfluidphases; ++idof)
    if (idof != phasetoadd)
      vtrans += rhsfac * phasemanager.PressureDeriv(curphase, idof) * phidtnp[idof];

  vtrans *= porosity * saturation * invbulkmodulus;

  return vtrans;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // in case of an incompressible solid (1/K_s = 0) the term cancels out
  if (phasemanager.IncompressibleSolid()) return;

  const int numfluidphases = phasemanager.NumFluidPhases();

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  //  get inverse bulkmodulus (=compressiblity)
  // TODO linearization of bulkmodulus
  const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  {
    const double facfacmass = fac * (1.0 - phasemanager.Porosity()) * invsolidbulkmodulus;
    for (int vi = 0; vi < nen; ++vi)
    {
      const double v = facfacmass * funct(vi);
      const int fvi = vi * numdofpernode + phasetoadd;

      for (int ui = 0; ui < nen; ++ui)
      {
        const double vfunct = v * funct(ui);
        for (int idof = 0; idof < numfluidphases; ++idof)
        {
          const int fui = ui * numdofpernode + idof;

          mymat(fvi, fui) += vfunct * phasemanager.SolidPressureDeriv(idof);
        }
      }
    }
  }

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if (!inittimederiv)
  {
    //----------------------------------------------------------------
    // linearization of solid pressure derivative w.r.t. dof
    //----------------------------------------------------------------
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
      double hist = 0.0;
      // if(curphase==phasetoadd) //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass3 = (1.0 - phasemanager.Porosity()) * invsolidbulkmodulus;

      std::vector<double> val(numfluidphases, 0.0);

      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        if (idof == phasetoadd)
          for (int jdof = 0; jdof < numfluidphases; ++jdof)
            val[jdof] +=
                fac * phasemanager.SolidPressureDerivDeriv(idof, jdof) * (phinp[idof] - hist);
        else
          for (int jdof = 0; jdof < numfluidphases; ++jdof)
            val[jdof] +=
                timefacfac * phasemanager.SolidPressureDerivDeriv(idof, jdof) * (phidtnp[idof]);
      }

      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = facfacmass3 * funct(vi);
        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int idof = 0; idof < numfluidphases; ++idof)
          {
            const int fui = ui * numdofpernode + idof;

            mymat(fvi, fui) += vfunct * val[idof];
          }
        }
      }
    }

    //----------------------------------------------------------------
    // linearization of porosity w.r.t. dof
    //----------------------------------------------------------------
    if (phasemanager.PorosityDependsOnFluid())
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
      double hist = 0.0;
      // if(curphase==phasetoadd) //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass3 = -1.0 * invsolidbulkmodulus;

      std::vector<double> val(numdofpernode, 0.0);

      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        const double solidpressurederiv = phasemanager.SolidPressureDeriv(idof);
        if (idof == phasetoadd)
          for (int jdof = 0; jdof < numdofpernode; ++jdof)
            val[jdof] +=
                fac * solidpressurederiv * phasemanager.PorosityDeriv(jdof) * (phinp[idof] - hist);
        else
          for (int jdof = 0; jdof < numdofpernode; ++jdof)
            val[jdof] += timefacfac * solidpressurederiv * phasemanager.PorosityDeriv(jdof) *
                         (phidtnp[idof]);
      }

      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = facfacmass3 * funct(vi);
        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int idof = 0; idof < numdofpernode; ++idof)
          {
            const int fui = ui * numdofpernode + idof;

            mymat(fvi, fui) += vfunct * val[idof];
          }
        }
      }
    }
  }  // !inittimederiv
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // in case of an incompressible solid (1/K_s = 0) the term cancels out
  if (phasemanager.IncompressibleSolid()) return;

  // for the initial time derivative calculation no transient terms enter the rhs
  if (!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    double vtrans = GetRhsTrans(
        curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, rhsfac, fac);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + phasetoadd;

      myvec[fvi] -= vtrans * funct(vi);
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // in case of an incompressible solid (1/K_s = 0) the term cancels out
  if (phasemanager.IncompressibleSolid()) return;

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  double vtrans = GetRhsTrans(
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, timefacfac, fac);

  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd =
  // dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
  // (dx/ds) * ( det(dX/ds) )^-1 in our case: vtrans is scaled with (1.0-porosity) in GetRhsTrans
  // --> scale it with 1.0/(1.0-porosity) here

  if (phasemanager.PorosityDependsOnStruct())
    vtrans += vtrans * (-1.0) / (1.0 - phasemanager.Porosity()) * phasemanager.JacobianDefGrad() *
              phasemanager.PorosityDerivWrtJacobianDefGrad();

  // linearization of mesh motion (Jacobian)
  // 1) linearization of fac +
  // 2) possible linearization w.r.t porosity
  EvaluatorBase<nsd, nen>::CalcLinFacODMesh(mymat, funct, derxy, vtrans, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate transient term at GP                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd, nen>::GetRhsTrans(
    int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac)
{
  const int numfluidphases = phasemanager.NumFluidPhases();
  // read data from managers
  double hist = 0.0;
  // if(curphase==phasetoadd) //bugfix??
  hist = (*variablemanager.Hist())[phasetoadd];
  const double porosity = phasemanager.Porosity();
  const std::vector<double>& phinp = *variablemanager.Phinp();
  const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

  //  get inverse bulkmodulus (=compressiblity)
  const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

  // TODO check genalpha
  // compute scalar at integration point
  double vtrans = fac * phasemanager.SolidPressureDeriv(phasetoadd) * (phinp[phasetoadd] - hist);

  for (int idof = 0; idof < numfluidphases; ++idof)
  {
    if (idof != phasetoadd)
    {
      vtrans += rhsfac * phasemanager.SolidPressureDeriv(idof) * phidtnp[idof];
    }
  }

  vtrans *= (1.0 - porosity) * invsolidbulkmodulus;

  return vtrans;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorMassSolidPressure<nsd, nen>::EvaluateMatrixAndAssemble(elemat, funct, derxy, curphase,
      phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * timefacfac, phasemanager.Saturation(curphase) * fac,
      inittimederiv);

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if (!inittimederiv)
  {
    const int numfluidphases = phasemanager.NumFluidPhases();
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    //  get inverse bulkmodulus (=compressiblity)
    // TODO linearization of bulkmodulus
    const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

    //----------------------------------------------------------------
    // linearization of saturation w.r.t. dof
    //----------------------------------------------------------------
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
      double hist = 0.0;
      // if(curphase==phasetoadd) //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass =
          fac * phasemanager.SolidPressureDeriv(phasetoadd) * (phinp[phasetoadd] - hist);

      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        if (idof != phasetoadd)
          facfacmass += timefacfac * phasemanager.SolidPressureDeriv(idof) * phidtnp[idof];
      }

      facfacmass *= (1.0 - phasemanager.Porosity()) * invsolidbulkmodulus;

      // call base class for saturation linearization
      EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(mymat, funct, facfacmass, numdofpernode,
          numfluidphases, curphase, phasetoadd, phasemanager);
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorMassSolidPressure<nsd, nen>::EvaluateVectorAndAssemble(elevec, funct, derxy, xyze,
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * rhsfac, phasemanager.Saturation(curphase) * fac,
      inittimederiv);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // call base class with scaled factors
  EvaluatorMassSolidPressure<nsd, nen>::EvaluateMatrixODStructAndAssemble(elemat, funct, deriv,
      derxy, xjm, curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * timefacfac, phasemanager.Saturation(curphase) * fac, det);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();

  //----------------------------------------------------------------
  // linearization of saturation w.r.t. dof
  //----------------------------------------------------------------
  {
    const double facfacmass = fac * phasemanager.Porosity();
    // call base class for saturation linearization
    EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(mymat, funct, facfacmass, numdofpernode,
        numfluidphases, curphase, phasetoadd, phasemanager);
  }

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if (!inittimederiv)
  {
    //----------------------------------------------------------------
    // linearization of porosity w.r.t. dof
    //----------------------------------------------------------------
    if (phasemanager.PorosityDependsOnFluid())
    {
      // read data from manager
      double hist = 0.0;
      // if(curphase==phasetoadd) //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

      // TODO genalpha
      // compute scalar at integration point
      double facfacmass =
          fac * phasemanager.SaturationDeriv(curphase, phasetoadd) * (phinp[phasetoadd] - hist);

      for (int idof = 0; idof < numfluidphases; ++idof)
        if (phasetoadd != idof)
          facfacmass += timefacfac * phasemanager.SaturationDeriv(curphase, idof) * phidtnp[idof];

      // call base class:
      EvaluatorBase<nsd, nen>::PorosityLinearizationFluid(
          mymat, funct, facfacmass, numdofpernode, phasetoadd, phasemanager);
    }

    //----------------------------------------------------------------
    // linearization of derivative of saturation w.r.t. dof
    //----------------------------------------------------------------
    {
      // read data from manager
      double hist = 0.0;
      hist = (*variablemanager.Hist())[phasetoadd];
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

      /*for (int iphase=0; iphase < numdofpernode; iphase++)
      {
        std::cout << iphase << " =====================================" << std::endl;
        for (int jphase = 0; jphase < numdofpernode; jphase++)
        {
          for (int kphase = 0; kphase < numdofpernode; kphase++)
          {
            std::cout << std::setprecision(8) <<
      phasemanager.SaturationDerivDeriv(iphase,jphase,kphase) << "  ";
          }
          std::cout << "\n";
        }
      }*/

      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = funct(vi);
        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui) * phasemanager.Porosity();
          for (int idof = 0; idof < numfluidphases; ++idof)
          {
            if (idof == phasetoadd)
            {
              for (int jdof = 0; jdof < numfluidphases; ++jdof)
              {
                const int fui = ui * numdofpernode + jdof;
                mymat(fvi, fui) += fac * vfunct *
                                   phasemanager.SaturationDerivDeriv(curphase, phasetoadd, jdof) *
                                   (phinp[phasetoadd] - hist);
              }
            }
            else
            {
              for (int jdof = 0; jdof < numfluidphases; ++jdof)
              {
                const int fui = ui * numdofpernode + jdof;
                mymat(fvi, fui) += timefacfac * vfunct *
                                   phasemanager.SaturationDerivDeriv(curphase, idof, jdof) *
                                   (phidtnp[idof]);
              }
            }
          }
        }
      }
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // for the initial time derivative calculation no transient terms enter the rhs
  if (!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    double vtrans = GetRhsTrans(
        curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, rhsfac, fac);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + phasetoadd;

      myvec[fvi] -= vtrans * funct(vi);
    }
  }  // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  double vtrans = GetRhsTrans(
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, timefacfac, fac);

  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd =
  // dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
  // (dx/ds) * ( det(dX/ds) )^-1 in our case: vtrans is scaled with porosity in GetRhsTrans -->
  // scale it with 1.0/porosity here

  if (phasemanager.PorosityDependsOnStruct())
    vtrans += vtrans * 1.0 / phasemanager.Porosity() * phasemanager.JacobianDefGrad() *
              phasemanager.PorosityDerivWrtJacobianDefGrad();

  // linearization of mesh motion (Jacobian)
  // 1) linearization of fac +
  // 2) possible linearization w.r.t porosity
  EvaluatorBase<nsd, nen>::CalcLinFacODMesh(mymat, funct, derxy, vtrans, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate transient term at GP                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd, nen>::GetRhsTrans(
    int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac)
{
  const int numfluidphases = phasemanager.NumFluidPhases();
  // read data from manager
  double hist = 0.0;
  // if(curphase==phasetoadd) //bugfix??
  hist = (*variablemanager.Hist())[phasetoadd];

  const double porosity = phasemanager.Porosity();
  const std::vector<double>& phinp = *variablemanager.Phinp();
  const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

  // TODO genalpha
  // compute scalar at integration point
  double vtrans =
      fac * phasemanager.SaturationDeriv(curphase, phasetoadd) * (phinp[phasetoadd] - hist);

  for (int idof = 0; idof < numfluidphases; ++idof)
    if (phasetoadd != idof)
      vtrans += rhsfac * phasemanager.SaturationDeriv(curphase, idof) * phidtnp[idof];
  // note: for one-step theta: rhsfac*phidtnp = theta*dt*(phinp-phin)/theta/dt+(1-theta)*phidtn
  //                                          = phinp - hist
  vtrans *= porosity;

  return vtrans;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get vectors to be filled
  // pressure
  Epetra_SerialDenseVector& pressure = *elevec[0];
  // saturation
  Epetra_SerialDenseVector& saturation = *elevec[1];
  // counter
  Epetra_SerialDenseVector& counter = *elevec[2];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // FLUID Phases:
  if (curphase < numfluidphases)
  {
    for (int inode = 0; inode < nen; inode++)
    {
      // save the pressure value
      pressure[inode * numdofpernode + curphase] +=
          fac * funct(inode) * phasemanager.Pressure(curphase);
      // save the saturation value
      saturation[inode * numdofpernode + curphase] +=
          fac * funct(inode) * phasemanager.Saturation(curphase);
      // mark the evaluated node
      counter[inode * numdofpernode + curphase] += fac * funct(inode);
    }
  }
  // VOLFRAC Phases:
  else if (curphase < numfluidphases + numvolfrac)
  {
    // dummy way: set pressures and saturations to -1
    // TODO: is there a better way to do it ??
    for (int inode = 0; inode < nen; inode++)
    {
      // save the pressure value
      pressure[inode * numdofpernode + curphase] += fac * funct(inode) * (-1.0);
      // save the saturation value
      saturation[inode * numdofpernode + curphase] += fac * funct(inode) * (-1.0);
      // mark the evaluated node
      counter[inode * numdofpernode + curphase] += fac * funct(inode);
    }
  }
  // VOLFRAC PRESSURE Phases:
  else if (curphase < numdofpernode)
  {
    // dummy way: set saturations to -1
    // TODO: is there a better way to do it ??
    for (int inode = 0; inode < nen; inode++)
    {
      // save the pressure value
      pressure[inode * numdofpernode + curphase] +=
          fac * funct(inode) * phasemanager.VolFracPressure(curphase - numfluidphases - numvolfrac);
      // save the saturation value
      saturation[inode * numdofpernode + curphase] += fac * funct(inode) * (-1.0);
      // mark the evaluated node
      counter[inode * numdofpernode + curphase] += fac * funct(inode);
    }
  }
  else
    dserror("wrong value for curphase: %i", curphase);
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get vectors to be filled
  Epetra_SerialDenseVector& solidpressure = *elevec[0];
  Epetra_SerialDenseVector& counter = *elevec[1];

  for (int inode = 0; inode < nen; inode++)
  {
    // save the pressure value
    solidpressure[inode] += fac * funct(inode) * (phasemanager.SolidPressure());
    // mark the evaluated node
    counter[inode] += fac * funct(inode);
  }
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorValidVolFracPressures<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorValidVolFracPressures<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  Epetra_SerialDenseVector& valid_volfracpress = *elevec[1];
  Epetra_SerialDenseVector& valid_volfracspec = *elevec[2];

  for (int inode = 0; inode < nen; inode++)
  {
    for (int idof = numfluidphases + numvolfrac; idof < numdofpernode; idof++)
    {
      const int fvi = inode * numdofpernode + idof;

      const bool evaluatevolfracpress =
          variablemanager.ElementHasValidVolFracPressure(idof - numfluidphases - numvolfrac);

      const bool evaluatevolfracspec =
          variablemanager.ElementHasValidVolFracSpecies(idof - numfluidphases - numvolfrac);

      if (evaluatevolfracpress) valid_volfracpress[fvi] = 1.0;
      if (evaluatevolfracspec) valid_volfracspec[fvi] = 1.0;
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorValidVolFracPressures<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorValidVolFracPressures<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 04/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 04/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get vectors to be filled
  Epetra_SerialDenseVector& porosity = *elevec[0];
  Epetra_SerialDenseVector& counter = *elevec[1];

  for (int inode = 0; inode < nen; inode++)
  {
    // save the porosity value
    porosity[inode] += fac * funct(inode) * phasemanager.Porosity();
    // mark the evaluated node
    counter[inode] += fac * funct(inode);
  }
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 04/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 03/19 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDomainIntegrals<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 03/19 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDomainIntegrals<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)

{
  // get vector to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // get the variables + constants
  std::vector<std::pair<std::string, double>> constants;
  // pressures + saturations + fluiddensities + porosity + volfracs + volfracpressures +
  // volfracdensities + scalars + numdim (x,y and possibly z)
  constants.reserve(3 * numfluidphases + 1 + 3 * numvolfrac + numscal_ + nsd);

  std::vector<std::pair<std::string, double>> variables;
  variables.reserve(0);

  // set pressure, saturation and density values as constants
  for (int k = 0; k < numfluidphases; k++)
  {
    std::ostringstream temp;
    temp << k + 1;
    constants.push_back(std::pair<std::string, double>("p" + temp.str(), phasemanager.Pressure(k)));
    constants.push_back(
        std::pair<std::string, double>("S" + temp.str(), phasemanager.Saturation(k)));
    constants.push_back(
        std::pair<std::string, double>("DENS" + temp.str(), phasemanager.Density(k)));
  }

  // set porosity value as constant
  constants.push_back(std::pair<std::string, double>("porosity", phasemanager.Porosity()));

  // set volfrac, volfrac pressure and volfrac density values as constants
  for (int k = 0; k < numvolfrac; k++)
  {
    std::ostringstream temp;
    temp << k + 1;
    constants.push_back(std::pair<std::string, double>("VF" + temp.str(), phasemanager.VolFrac(k)));
    constants.push_back(
        std::pair<std::string, double>("VFP" + temp.str(), phasemanager.VolFracPressure(k)));
    constants.push_back(
        std::pair<std::string, double>("VFDENS" + temp.str(), phasemanager.VolFracDensity(k)));
  }

  // set scalar values as constants
  for (int k = 0; k < numscal_; k++)
  {
    std::ostringstream temp;
    temp << k + 1;
    constants.push_back(
        std::pair<std::string, double>("phi" + temp.str(), variablemanager.Scalarnp()->at(k)));
  }

  // calculate the coordinates of the gauss point
  std::vector<double> coords(nsd, 0.0);
  for (int idim = 0; idim < nsd; idim++)
  {
    for (int inode = 0; inode < nen; inode++)
    {
      coords[idim] += funct(inode) * xyze(idim, inode);
    }
  }

  // set values as constants in function
  constants.push_back(std::pair<std::string, double>("x", coords[0]));
  constants.push_back(std::pair<std::string, double>("y", coords[1]));
  if (nsd == 3) constants.push_back(std::pair<std::string, double>("z", coords[2]));

  // call the functions and integrate value (multiply with fac)
  for (unsigned int i = 0; i < domainint_funct_.size(); i++)
  {
    myvec[i] += Function(domainint_funct_[i] - 1)
                    .DRT::UTILS::VariableExprFunction::Evaluate(0, variables, constants) *
                fac;
  }
};

/*----------------------------------------------------------------------*
 | cast to VarExp-function                                              |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
inline DRT::UTILS::VariableExprFunction&
DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDomainIntegrals<nsd, nen>::Function(int functnum) const
{
  try
  {
    DRT::UTILS::VariableExprFunction& funct =
        dynamic_cast<DRT::UTILS::VariableExprFunction&>(DRT::Problem::Instance()->Funct(functnum));
    if (funct.NumberComponents() != 1)
      dserror("only one component allowed for domain integral functions");
    return funct;
  }
  catch (std::bad_cast& exp)
  {
    dserror(
        "Cast to VarExp Function failed! For domain integrals only 'VARFUNCTION' functions are "
        "allowed!\n"
        "Check your input file!");
    return dynamic_cast<DRT::UTILS::VariableExprFunction&>(
        DRT::Problem::Instance()->Funct(functnum));
  }
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/19 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDomainIntegrals<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 03/19 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDomainIntegrals<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // get matrixes to fill
  Epetra_SerialDenseMatrix& linearization = *elemat[0];
  // Compute element matrix. For L2-projection
  for (int vi = 0; vi < nen; ++vi)
  {
    const double v = fac * funct(vi);
    for (int ui = 0; ui < nen; ++ui)
    {
      linearization(vi, ui) += v * funct(ui);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // nothing to do
  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // get matrixes to fill
  Epetra_SerialDenseMatrix& rhs = *elemat[1];

  const int numfluidphases = phasemanager.NumFluidPhases();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // current pressure gradient
  LINALG::Matrix<nsd, 1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof = 0; idof < numfluidphases; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase, idof), gradphi[idof], 1.0);

  double abspressgrad = 0.0;
  for (int i = 0; i < nsd; i++) abspressgrad += gradpres(i) * gradpres(i);
  abspressgrad = sqrt(abspressgrad);

  // diffusion tensor
  LINALG::Matrix<nsd, nsd> difftensor(true);
  phasemanager.PermeabilityTensor(curphase, difftensor);
  difftensor.Scale(
      phasemanager.RelPermeability(curphase) / phasemanager.DynViscosity(curphase, abspressgrad));

  // diffusive flux
  static LINALG::Matrix<nsd, 1> diffflux(true);
  diffflux.Multiply(-1.0, difftensor, gradpres);

  // Compute element vectors. For L2-Projection
  for (int node_i = 0; node_i < nen; node_i++)
  {
    for (int j = 0; j < nsd; j++)
    {
      rhs(node_i, nsd * curphase + j) += funct(node_i) * fac * diffflux(j);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // nothing to do
  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTerms<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];
  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  //----------------------------------------------------------------
  // 1) standard Galerkin transient term
  //----------------------------------------------------------------
  {
    const double facfacmass = -fac;
    for (int vi = 0; vi < nen; ++vi)
    {
      const double v = facfacmass * funct(vi);
      const int fvi = vi * numdofpernode + phasetoadd;

      for (int ui = 0; ui < nen; ++ui)
      {
        const double vfunct = v * funct(ui);
        for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ++ivolfrac)
        {
          const int fui = ui * numdofpernode + ivolfrac;

          mymat(fvi, fui) += vfunct;
        }
      }
    }
  }

  //----------------------------------------------------------------
  // 2) - sum_volfrac porosity_volfrac/K_s * d p_s / d t
  //----------------------------------------------------------------
  if (not phasemanager.IncompressibleSolid())
  {
    const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();
    const double sumaddvolfrac = phasemanager.SumAddVolFrac();

    //----------------------------------------------------------------
    // standard Galerkin transient term
    //----------------------------------------------------------------
    {
      const double facfacmass = fac * (-sumaddvolfrac) * invsolidbulkmodulus;
      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = facfacmass * funct(vi);
        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int idof = 0; idof < numfluidphases; ++idof)
          {
            const int fui = ui * numdofpernode + idof;

            mymat(fvi, fui) += vfunct * phasemanager.SolidPressureDeriv(idof);
          }
        }
      }
    }
    // for the initial time derivative calculation no additional derivatives are needed
    if (!inittimederiv)
    {
      //----------------------------------------------------------------
      // linearization of solid pressure derivative w.r.t. dof
      //----------------------------------------------------------------
      {
        const std::vector<double>& phinp = *variablemanager.Phinp();
        const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
        double hist = 0.0;
        // if(curphase==phasetoadd) //bugfix??
        hist = (*variablemanager.Hist())[phasetoadd];

        const double sumaddvolfrac = phasemanager.SumAddVolFrac();

        double facfacmass3 = -sumaddvolfrac * invsolidbulkmodulus;

        std::vector<double> val(numfluidphases, 0.0);

        for (int idof = 0; idof < numfluidphases; ++idof)
        {
          if (idof == phasetoadd)
            for (int jdof = 0; jdof < numfluidphases; ++jdof)
              val[jdof] +=
                  fac * phasemanager.SolidPressureDerivDeriv(idof, jdof) * (phinp[idof] - hist);
          else
            for (int jdof = 0; jdof < numfluidphases; ++jdof)
              val[jdof] +=
                  timefacfac * phasemanager.SolidPressureDerivDeriv(idof, jdof) * (phidtnp[idof]);
        }

        for (int vi = 0; vi < nen; ++vi)
        {
          const double v = facfacmass3 * funct(vi);
          const int fvi = vi * numdofpernode + phasetoadd;

          for (int ui = 0; ui < nen; ++ui)
          {
            const double vfunct = v * funct(ui);
            for (int idof = 0; idof < numfluidphases; ++idof)
            {
              const int fui = ui * numdofpernode + idof;

              mymat(fvi, fui) += vfunct * val[idof];
            }
          }
        }
      }
      //----------------------------------------------------------------
      // linearization of sum_volfrac porosity_volfrac w.r.t. dof
      //----------------------------------------------------------------
      double hist = 0.0;
      // if(curphase==phasetoadd) //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

      // TODO check genalpha
      // compute scalar at integration point
      double vtrans =
          fac * phasemanager.SolidPressureDeriv(phasetoadd) * (phinp[phasetoadd] - hist);

      for (int idof = 0; idof < numfluidphases; ++idof)
      {
        if (idof != phasetoadd)
        {
          vtrans += timefacfac * phasemanager.SolidPressureDeriv(idof) * phidtnp[idof];
        }
      }
      vtrans *= -invsolidbulkmodulus;
      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = vtrans * funct(vi);
        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ++ivolfrac)
          {
            const int fui = ui * numdofpernode + ivolfrac;

            mymat(fvi, fui) += vfunct;
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTerms<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get vector to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  // for the initial time derivative calculation no transient terms enter the rhs
  if (!inittimederiv)
  {
    const double vrhs =
        GetRhs(curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, rhsfac, fac);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + phasetoadd;

      myvec[fvi] -= vrhs * funct(vi);
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTerms<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const double vrhs =
      GetRhs(curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, timefacfac, fac);

  // linearization of mesh motion (Jacobian)
  EvaluatorBase<nsd, nen>::CalcLinFacODMesh(mymat, funct, derxy, vrhs, numdofpernode, phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTerms<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate rhs term at GP                             kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTerms<nsd, nen>::GetRhs(
    int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac)
{
  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  const std::vector<double>& phidtnp = *variablemanager.Phidtnp();

  double vrhs = 0.0;

  // sum^volfrac \frac{\partial phi_volfrac}{\partial t}
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
    vrhs -= rhsfac * phidtnp[ivolfrac];

  // \frac{-\sum^volfrac \phi_volfrac) }{K_s} \frac{\partial p^s}{\partial t}
  if (not phasemanager.IncompressibleSolid())
  {
    double hist = 0.0;
    // if(curphase==phasetoadd) //bugfix??
    hist = (*variablemanager.Hist())[phasetoadd];
    const std::vector<double>& phinp = *variablemanager.Phinp();

    //  get inverse bulkmodulus (=compressiblity)
    const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

    // TODO check genalpha
    // compute scalar at integration point
    double vtrans = fac * phasemanager.SolidPressureDeriv(phasetoadd) * (phinp[phasetoadd] - hist);

    for (int idof = 0; idof < numfluidphases; ++idof)
    {
      if (idof != phasetoadd)
      {
        vtrans += rhsfac * phasemanager.SolidPressureDeriv(idof) * phidtnp[idof];
      }
    }
    vtrans *= -invsolidbulkmodulus;
    const double sumaddvolfrac = phasemanager.SumAddVolFrac();
    vrhs += vtrans * sumaddvolfrac;
  }

  return vrhs;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTerm<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];
    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    //----------------------------------------------------------------
    // - sum_volfrac porosity_volfrac * div v_s
    //----------------------------------------------------------------
    {
      const double prefac = -timefacfac * variablemanager.DivConVelnp();
      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = prefac * funct(vi);
        const int fvi = vi * numdofpernode + phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ++ivolfrac)
          {
            const int fui = ui * numdofpernode + ivolfrac;

            mymat(fvi, fui) += vfunct;
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTerm<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get vector to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const double vrhs = -rhsfac * phasemanager.SumAddVolFrac() * variablemanager.DivConVelnp();

  for (int vi = 0; vi < nen; ++vi)
  {
    const int fvi = vi * numdofpernode + phasetoadd;

    myvec[fvi] -= vrhs * funct(vi);
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTerm<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const double sumaddvolfrac = phasemanager.SumAddVolFrac();

  LINALG::Matrix<nsd, nsd> gridvelderiv(true);
  gridvelderiv.MultiplyNT(*(variablemanager.EConVelnp()), deriv);

  // OD mesh - div vel term
  EvaluatorBase<nsd, nen>::CalcDivVelODMesh(mymat, funct, deriv, derxy, xjm, gridvelderiv,
      timefacfac * sumaddvolfrac * (-1.0), fac * sumaddvolfrac * (-1.0), det, numdofpernode,
      phasetoadd);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTerm<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTermsSat<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorVolFracAddInstatTerms<nsd, nen>::EvaluateMatrixAndAssemble(elemat, funct, derxy,
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      timefacfac * phasemanager.Saturation(curphase), fac * phasemanager.Saturation(curphase),
      inittimederiv);

  // we do not need additional linearizations if we calculate the initial time derivative
  if (!inittimederiv)
  {
    const int numfluidphases = phasemanager.NumFluidPhases();
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    //----------------------------------------------------------------
    // linearization of saturation w.r.t. dof
    //----------------------------------------------------------------

    // first: get rhs
    const double vrhs = EvaluatorVolFracAddInstatTerms<nsd, nen>::GetRhs(
        curphase, phasetoadd, numdofpernode, phasemanager, variablemanager, timefacfac, fac);

    // call base class for saturation linearization
    EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(
        mymat, funct, vrhs, numdofpernode, numfluidphases, curphase, phasetoadd, phasemanager);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTermsSat<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorVolFracAddInstatTerms<nsd, nen>::EvaluateVectorAndAssemble(elevec, funct, derxy, xyze,
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * rhsfac, phasemanager.Saturation(curphase) * fac,
      inittimederiv);

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTermsSat<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // call base class with scaled factors
  EvaluatorVolFracAddInstatTerms<nsd, nen>::EvaluateMatrixODStructAndAssemble(elemat, funct, deriv,
      derxy, xjm, curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * timefacfac, phasemanager.Saturation(curphase) * fac, det);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddInstatTermsSat<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTermSat<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorVolFracAddDivVelTerm<nsd, nen>::EvaluateMatrixAndAssemble(elemat, funct, derxy, curphase,
      phasetoadd, numdofpernode, phasemanager, variablemanager,
      timefacfac * phasemanager.Saturation(curphase), fac * phasemanager.Saturation(curphase),
      inittimederiv);

  // we do not need additional linearizations if we calculate the initial time derivative
  if (!inittimederiv)
  {
    const int numfluidphases = phasemanager.NumFluidPhases();
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const double vrhs = -timefacfac * phasemanager.SumAddVolFrac() * variablemanager.DivConVelnp();

    // call base class for saturation linearization
    EvaluatorBase<nsd, nen>::SaturationLinearizationFluid(
        mymat, funct, vrhs, numdofpernode, numfluidphases, curphase, phasetoadd, phasemanager);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTermSat<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // call base class with scaled factors
  EvaluatorVolFracAddDivVelTerm<nsd, nen>::EvaluateVectorAndAssemble(elevec, funct, derxy, xyze,
      curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * rhsfac, phasemanager.Saturation(curphase) * fac,
      inittimederiv);

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTermSat<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // call base class with scaled factors
  EvaluatorVolFracAddDivVelTerm<nsd, nen>::EvaluateMatrixODStructAndAssemble(elemat, funct, deriv,
      derxy, xjm, curphase, phasetoadd, numdofpernode, phasemanager, variablemanager,
      phasemanager.Saturation(curphase) * timefacfac, phasemanager.Saturation(curphase) * fac, det);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 09/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddDivVelTermSat<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracInstat<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfrac - numfluidphases);

    //----------------------------------------------------------------
    // standard Galerkin transient term
    //----------------------------------------------------------------
    for (int vi = 0; vi < nen; ++vi)
    {
      const double v = funct(vi) * fac;
      const int fvi_volfrac = vi * numdofpernode + ivolfrac;
      const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

      for (int ui = 0; ui < nen; ++ui)
      {
        const int fui = ui * numdofpernode + ivolfrac;

        mymat(fvi_volfrac, fui) += v * funct(ui);
        if (evaluatevolfracpress) mymat(fvi_volfracpress, fui) += v * funct(ui);
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracInstat<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // for the initial time derivative calculation no transient terms enter the rhs
  if (!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    // loop over all volume fractions
    for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
    {
      const double hist = (*variablemanager.Hist())[ivolfrac];
      const double phinp = phasemanager.VolFrac(ivolfrac - numfluidphases);
      const double vtrans = fac * (phinp - hist);

      const bool evaluatevolfracpress =
          variablemanager.ElementHasValidVolFracPressure(ivolfrac - numfluidphases);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi_volfrac = vi * numdofpernode + ivolfrac;
        const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

        myvec[fvi_volfrac] -= vtrans * funct(vi);
        if (evaluatevolfracpress) myvec[fvi_volfracpress] -= vtrans * funct(vi);
      }
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracInstat<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    const double hist = (*variablemanager.Hist())[ivolfrac];
    const double phinp = phasemanager.VolFrac(ivolfrac - numfluidphases);
    const double vtrans = fac * (phinp - hist);

    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfrac - numfluidphases);

    // linearization of mesh motion
    //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} =
    // J * N_x
    // J denotes the determinant of the Jacobian of the mapping between current and parameter space,
    // i.e. det(dx/ds) in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac *
    // N_x
    //              fac        = J              --> d(fac)/dd        = fac * N_x
    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi_volfrac = vi * numdofpernode + ivolfrac;
      const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

      const double v = funct(vi) * vtrans;

      for (int ui = 0; ui < nen; ++ui)
      {
        for (int idim = 0; idim < nsd; ++idim)
        {
          const int fui = ui * nsd + idim;

          mymat(fvi_volfrac, fui) += v * derxy(idim, ui);
          if (evaluatevolfracpress) mymat(fvi_volfracpress, fui) += v * derxy(idim, ui);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracInstat<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDivVel<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // no linearization needed in case of initial time derivative calculation
  if (!inittimederiv)
  {
    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const double consfac = timefacfac * variablemanager.DivConVelnp();

    // loop over all volume fractions
    for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
    {
      const bool evaluatevolfracpress =
          variablemanager.ElementHasValidVolFracPressure(ivolfrac - numfluidphases);

      for (int vi = 0; vi < nen; ++vi)
      {
        const double v = consfac * funct(vi);
        const int fvi_volfrac = vi * numdofpernode + ivolfrac;
        const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

        for (int ui = 0; ui < nen; ++ui)
        {
          const int fui = ui * numdofpernode + ivolfrac;
          mymat(fvi_volfrac, fui) += v * funct(ui);
          if (evaluatevolfracpress) mymat(fvi_volfracpress, fui) += v * funct(ui);
        }
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDivVel<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  double vrhs = rhsfac * variablemanager.DivConVelnp();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    const double v = vrhs * phasemanager.VolFrac(ivolfrac - numfluidphases);
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfrac - numfluidphases);


    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi_volfrac = vi * numdofpernode + ivolfrac;
      const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

      myvec[fvi_volfrac] -= v * funct(vi);
      if (evaluatevolfracpress) myvec[fvi_volfracpress] -= v * funct(vi);
    }
  }
  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDivVel<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfrac - numfluidphases);

    const double vrhs = fac * phasemanager.VolFrac(ivolfrac - numfluidphases);
    // d (div v_s)/d d_n+1 = derxy * 1.0/theta/dt * d_n+1
    // prefactor is fac since timefacfac/theta/dt = fac
    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi_volfrac = vi * numdofpernode + ivolfrac;
      const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;
      const double v = vrhs * funct(vi);

      for (int ui = 0; ui < nen; ++ui)
      {
        for (int idim = 0; idim < nsd; ++idim)
        {
          const int fui = ui * nsd + idim;

          mymat(fvi_volfrac, fui) += v * derxy(idim, ui);
          if (evaluatevolfracpress) mymat(fvi_volfracpress, fui) += v * derxy(idim, ui);
        }
      }
    }

    // shapederivatives see fluid_ele_calc_poro.cpp
    LINALG::Matrix<nsd, nsd> gridvelderiv(true);
    gridvelderiv.MultiplyNT(*(variablemanager.EConVelnp()), deriv);

    if (nsd == 3)
    {
      const double gridvelderiv_0_0 = gridvelderiv(0, 0);
      const double gridvelderiv_0_1 = gridvelderiv(0, 1);
      const double gridvelderiv_0_2 = gridvelderiv(0, 2);
      const double gridvelderiv_1_0 = gridvelderiv(1, 0);
      const double gridvelderiv_1_1 = gridvelderiv(1, 1);
      const double gridvelderiv_1_2 = gridvelderiv(1, 2);
      const double gridvelderiv_2_0 = gridvelderiv(2, 0);
      const double gridvelderiv_2_1 = gridvelderiv(2, 1);
      const double gridvelderiv_2_2 = gridvelderiv(2, 2);

      const double xjm_0_0 = xjm(0, 0);
      const double xjm_0_1 = xjm(0, 1);
      const double xjm_0_2 = xjm(0, 2);
      const double xjm_1_0 = xjm(1, 0);
      const double xjm_1_1 = xjm(1, 1);
      const double xjm_1_2 = xjm(1, 2);
      const double xjm_2_0 = xjm(2, 0);
      const double xjm_2_1 = xjm(2, 1);
      const double xjm_2_2 = xjm(2, 2);

      for (int ui = 0; ui < nen; ++ui)
      {
        const double v0 =
            gridvelderiv_1_0 * derxjm_(0, 0, 1, ui) + gridvelderiv_1_1 * derxjm_(0, 1, 1, ui) +
            gridvelderiv_1_2 * derxjm_(0, 2, 1, ui) + gridvelderiv_2_0 * derxjm_(0, 0, 2, ui) +
            gridvelderiv_2_1 * derxjm_(0, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(0, 2, 2, ui);

        const double v1 =
            gridvelderiv_0_0 * derxjm_(1, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(1, 1, 0, ui) +
            gridvelderiv_0_2 * derxjm_(1, 2, 0, ui) + gridvelderiv_2_0 * derxjm_(1, 0, 2, ui) +
            gridvelderiv_2_1 * derxjm_(1, 1, 2, ui) + gridvelderiv_2_2 * derxjm_(1, 2, 2, ui);

        const double v2 =
            gridvelderiv_0_0 * derxjm_(2, 0, 0, ui) + gridvelderiv_0_1 * derxjm_(2, 1, 0, ui) +
            gridvelderiv_0_2 * derxjm_(2, 2, 0, ui) + gridvelderiv_1_0 * derxjm_(2, 0, 1, ui) +
            gridvelderiv_1_1 * derxjm_(2, 1, 1, ui) + gridvelderiv_1_2 * derxjm_(2, 2, 1, ui);

        for (int vi = 0; vi < nen; ++vi)
        {
          const int fvi_volfrac = vi * numdofpernode + ivolfrac;
          const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

          const double v =
              timefacfac / det * funct(vi) * phasemanager.VolFrac(ivolfrac - numfluidphases);

          mymat(fvi_volfrac, ui * 3 + 0) += v * v0;
          if (evaluatevolfracpress) mymat(fvi_volfracpress, ui * 3 + 0) += v * v0;

          mymat(fvi_volfrac, ui * 3 + 1) += v * v1;
          if (evaluatevolfracpress) mymat(fvi_volfracpress, ui * 3 + 1) += v * v1;

          mymat(fvi_volfrac, ui * 3 + 2) += v * v2;
          if (evaluatevolfracpress) mymat(fvi_volfracpress, ui * 3 + 2) += v * v2;
        }
      }
    }
    else if (nsd == 2)
    {
      const double gridvelderiv_0_0 = gridvelderiv(0, 0);
      const double gridvelderiv_0_1 = gridvelderiv(0, 1);
      const double gridvelderiv_1_0 = gridvelderiv(1, 0);
      const double gridvelderiv_1_1 = gridvelderiv(1, 1);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi_volfrac = vi * numdofpernode + ivolfrac;
        const int fvi_volfracpress = vi * numdofpernode + numvolfrac + ivolfrac;

        const double v =
            timefacfac / det * funct(vi) * phasemanager.VolFrac(ivolfrac - numfluidphases);

        for (int ui = 0; ui < nen; ++ui)
        {
          mymat(fvi_volfrac, ui * 2) +=
              v * (-gridvelderiv_1_0 * deriv(1, ui) + gridvelderiv_1_1 * deriv(0, ui));
          if (evaluatevolfracpress)
            mymat(fvi_volfracpress, ui * 2) +=
                v * (-gridvelderiv_1_0 * deriv(1, ui) + gridvelderiv_1_1 * deriv(0, ui));

          mymat(fvi_volfrac, ui * 2 + 1) +=
              v * (+gridvelderiv_0_0 * deriv(1, ui) - gridvelderiv_0_1 * deriv(0, ui));
          if (evaluatevolfracpress)
            mymat(fvi_volfracpress, ui * 2 + 1) +=
                v * (+gridvelderiv_0_0 * deriv(1, ui) - gridvelderiv_0_1 * deriv(0, ui));
        }
      }
    }
    else
      dserror("shapederivatives not implemented for 1D!");
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDivVel<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDiff<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    // loop over all volume fractions
    for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
    {
      // get difftensor and diffusive flux
      LINALG::Matrix<nsd, nsd> difftensor(true);
      phasemanager.DiffTensorVolFrac(ivolfrac - numfluidphases, difftensor);

      static LINALG::Matrix<nsd, nen> diffflux(true);
      diffflux.Multiply(difftensor, derxy);

      // diffusive term
      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + ivolfrac;

        for (int ui = 0; ui < nen; ++ui)
        {
          double laplawf(0.0);
          for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j, ui);

          const int fui = ui * numdofpernode + ivolfrac;
          mymat(fvi, fui) += timefacfac * laplawf;
        }
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDiff<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    // diffusion tensor
    LINALG::Matrix<nsd, nsd> difftensor(true);
    phasemanager.DiffTensorVolFrac(ivolfrac - numfluidphases, difftensor);

    static LINALG::Matrix<nsd, 1> diffflux(true);
    diffflux.Multiply(difftensor, gradphi[ivolfrac]);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi * numdofpernode + ivolfrac;

      // laplacian in weak form
      double laplawf(0.0);
      for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j);
      myvec[fvi] -= rhsfac * laplawf;
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDiff<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    // diffusion tensor
    LINALG::Matrix<nsd, nsd> difftensor(true);
    phasemanager.DiffTensorVolFrac(ivolfrac - numfluidphases, difftensor);

    static LINALG::Matrix<nsd, 1> diffflux(true);
    diffflux.Multiply(difftensor, gradphi[ivolfrac]);

    // TODO: anisotropic difftensor
    const double v = difftensor(0, 0) * timefacfac / det;

    // gradient of phi w.r.t. reference coordinates
    LINALG::Matrix<nsd, 1> refgradphi(true);
    refgradphi.Multiply(xjm, gradphi[ivolfrac]);

    // OD mesh - diffusive term
    EvaluatorBase<nsd, nen>::CalcDiffODMesh(mymat, deriv, derxy, xjm, diffflux, refgradphi,
        gradphi[ivolfrac], timefacfac, v, numdofpernode, ivolfrac);
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracDiff<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracReac<nsd, nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>& elemat, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    // loop over all volume fractions
    for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
    {
      if (phasemanager.IsReactive(ivolfrac))
      {
        double scaledtimefacfac =
            timefacfac / phasemanager.VolFracDensity(ivolfrac - numfluidphases);
        //----------------------------------------------------------------
        // reaction terms
        //----------------------------------------------------------------
        for (int vi = 0; vi < nen; ++vi)
        {
          const double v = scaledtimefacfac * funct(vi);
          const int fvi = vi * numdofpernode + ivolfrac;

          for (int ui = 0; ui < nen; ++ui)
          {
            const double vfunct = v * funct(ui);
            for (int idof = 0; idof < numdofpernode; ++idof)
            {
              const int fui = ui * numdofpernode + idof;

              // rhs ---> -
              mymat(fvi, fui) -= vfunct * phasemanager.ReacDeriv(ivolfrac, idof);
            }
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracReac<nsd, nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>& elevec, const LINALG::Matrix<nen, 1>& funct,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nen>& xyze, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    if (phasemanager.IsReactive(ivolfrac))
    {
      double scale = 1.0 / phasemanager.VolFracDensity(ivolfrac - numfluidphases);

      double vrhs = scale * rhsfac * phasemanager.ReacTerm(ivolfrac);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + ivolfrac;
        // rhs ---> +
        myvec[fvi] += vrhs * funct(vi);
      }
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracReac<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    if (phasemanager.IsReactive(ivolfrac))
    {
      // TODO a constant density is assumed here
      double scale = 1.0 / phasemanager.VolFracDensity(ivolfrac - numfluidphases);

      double vrhs = scale * timefacfac * phasemanager.ReacTerm(ivolfrac);

      // linearization of porosity (may appear in reaction term)
      //-------------dreac/dd = dreac/dporosity * dporosity/dd = dreac/dporosity * dporosity/dJ *
      // dJ/dd = dreac/dporosity * dporosity/dJ * J * N_x
      // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
      // (dx/ds) * ( det(dX/ds) )^-1

      if (phasemanager.PorosityDependsOnStruct())
      {
        vrhs += timefacfac * scale * phasemanager.ReacDerivPorosity(ivolfrac) *
                phasemanager.JacobianDefGrad() * phasemanager.PorosityDerivWrtJacobianDefGrad();
      }

      // linearization of mesh motion (Jacobian)
      // 1) linearization of fac +
      // 2) possible linearization w.r.t porosity
      // rhs ---> -
      EvaluatorBase<nsd, nen>::CalcLinFacODMesh(
          mymat, funct, derxy, -1.0 * vrhs, numdofpernode, ivolfrac);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracReac<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();
  const int numscal = phasemanager.NumScal();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    if (phasemanager.IsReactive(ivolfrac))
    {
      double vrhs = 1.0 / phasemanager.VolFracDensity(ivolfrac - numfluidphases) * timefacfac;

      // linearization of reaction term w.r.t scalars
      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + ivolfrac;
        const double v = vrhs * funct(vi);

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int iscal = 0; iscal < numscal; ++iscal)
          {
            const int fui = ui * numscal + iscal;
            // rhs ---> -
            mymat(fvi, fui) -= vfunct * phasemanager.ReacDerivScalar(ivolfrac, iscal);
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddFlux<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    // loop over all volume fractions
    for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
    {
      if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases))
      {
        // only in case of additional flux we have access to numscal
        const int numscal = phasemanager.NumScal();
        const std::vector<LINALG::Matrix<nsd, 1>>& gradscalarnp = *variablemanager.GradScalarnp();

        // loop over scalars
        for (int iscal = 0; iscal < numscal; iscal++)
        {
          if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases, iscal))
          {
            // diffusion tensor and diffusive flux
            LINALG::Matrix<nsd, nsd> difftensoraddflux(true);
            for (int i = 0; i < nsd; i++)
              difftensoraddflux(i, i) = phasemanager.ScalarDiff(ivolfrac - numfluidphases, iscal);

            static LINALG::Matrix<nsd, 1> diffflux(true);
            diffflux.Multiply(difftensoraddflux, gradscalarnp[iscal]);

            for (int vi = 0; vi < nen; ++vi)
            {
              const int fvi = vi * numdofpernode + ivolfrac;
              double laplawf = 0.0;
              for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j);
              for (int ui = 0; ui < nen; ++ui)
              {
                const double vfunct = timefacfac * funct(ui) * laplawf;
                // derivative w.r.t. fluid phases
                for (int idof = 0; idof < numfluidphases; ++idof)
                {
                  const int fui = ui * numdofpernode + idof;

                  // chemotaxis
                  if (phasemanager.ScalarToPhase(iscal).species_type ==
                      MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
                  {
                    if (phasemanager.ScalarToPhase(iscal).phaseID > phasemanager.NumFluidPhases())
                      dserror("Wrong PhaseID");
                    // 1) saturation deriv
                    // 2) porosity deriv
                    mymat(fvi, fui) +=
                        vfunct *
                        (phasemanager.SaturationDeriv(
                             phasemanager.ScalarToPhase(iscal).phaseID, idof) *
                                phasemanager.VolFrac(ivolfrac - numfluidphases) *
                                phasemanager.Porosity() +
                            phasemanager.PorosityDeriv(idof) *
                                phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID) *
                                phasemanager.VolFrac(ivolfrac - numfluidphases));
                  }
                  // haptotaxis
                  else if (phasemanager.ScalarToPhase(iscal).species_type ==
                           MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
                    // derivative of solid phase volume fraction w.r.t. all fluid phases = 0
                    break;
                  else
                    dserror("AddScalarDependentFlux only possible for species in fluid or solid!");
                }
                // derivative w.r.t. volfrac phases
                for (int jvolfrac = numfluidphases; jvolfrac < numdofpernode; ++jvolfrac)
                {
                  const int fui = ui * numdofpernode + jvolfrac;
                  // haptotaxis
                  if (phasemanager.ScalarToPhase(iscal).species_type ==
                      MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
                  {
                    // 1) derivative w.r.t. current volume fraction ivolfrac
                    if (ivolfrac == jvolfrac)
                      mymat(fvi, fui) +=
                          vfunct * (1.0 - phasemanager.Porosity() - phasemanager.SumAddVolFrac());
                    // 2) derivative of solid phase volume fraction w.r.t. all volume fractions = 0
                  }
                  // chemotaxis
                  else if (phasemanager.ScalarToPhase(iscal).species_type ==
                           MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
                  {
                    // 1) derivative w.r.t. current volume fraction ivolfrac
                    if (ivolfrac == jvolfrac)
                      mymat(fvi, fui) += vfunct * (phasemanager.Saturation(
                                                       phasemanager.ScalarToPhase(iscal).phaseID) *
                                                      phasemanager.Porosity());
                    // 2) porosity deriv w.r.t. all volume fractions
                    mymat(fvi, fui) +=
                        vfunct *
                        (phasemanager.PorosityDeriv(jvolfrac) *
                            phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID) *
                            phasemanager.VolFrac(ivolfrac - numfluidphases));
                  }
                  else
                    dserror("AddScalarDependentFlux only possible for species in fluid or solid!");
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddFlux<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases))
    {
      // only in case of additional flux we have access to numscal
      const int numscal = phasemanager.NumScal();
      const std::vector<LINALG::Matrix<nsd, 1>>& gradscalarnp = *variablemanager.GradScalarnp();

      // loop over scalars
      for (int iscal = 0; iscal < numscal; iscal++)
      {
        if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases, iscal))
        {
          // diffusion tensor
          LINALG::Matrix<nsd, nsd> difftensoraddflux(true);
          for (int i = 0; i < nsd; i++)
            difftensoraddflux(i, i) = phasemanager.ScalarDiff(ivolfrac - numfluidphases, iscal);

          // haptotaxis: scale with volfrac * (1 - porosity - sumaddvolfrac)
          if (phasemanager.ScalarToPhase(iscal).species_type ==
              MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
          {
            difftensoraddflux.Scale(phasemanager.VolFrac(ivolfrac - numfluidphases) *
                                    (1.0 - phasemanager.Porosity() - phasemanager.SumAddVolFrac()));
          }
          // chemotaxis: scale with volfrac * S * porosity
          else if (phasemanager.ScalarToPhase(iscal).species_type ==
                   MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
          {
            difftensoraddflux.Scale(
                phasemanager.VolFrac(ivolfrac - numfluidphases) * phasemanager.Porosity() *
                phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID));
          }
          else
            dserror("AddScalarDependentFlux only possible for species in fluid or solid!");

          static LINALG::Matrix<nsd, 1> diffflux(true);
          diffflux.Multiply(difftensoraddflux, gradscalarnp[iscal]);
          for (int vi = 0; vi < nen; ++vi)
          {
            const int fvi = vi * numdofpernode + ivolfrac;

            // laplacian in weak form
            double laplawf(0.0);
            for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j);
            myvec[fvi] -= rhsfac * laplawf;
          }
        }
      }
    }
  }
  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddFlux<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases))
    {
      // only in case of additional flux we have access to numscal
      const int numscal = phasemanager.NumScal();
      const std::vector<LINALG::Matrix<nsd, 1>>& gradscalarnp = *variablemanager.GradScalarnp();

      // loop over all scalars
      for (int iscal = 0; iscal < numscal; iscal++)
      {
        if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases, iscal))
        {
          // diffusion tensor
          LINALG::Matrix<nsd, nsd> difftensoraddflux(true);
          for (int i = 0; i < nsd; i++)
            difftensoraddflux(i, i) = phasemanager.ScalarDiff(ivolfrac - numfluidphases, iscal);

          static LINALG::Matrix<nsd, 1> diffflux(true);
          // haptotaxis: scale with volfrac * (1 - porosity - sumaddvolfrac)
          if (phasemanager.ScalarToPhase(iscal).species_type ==
              MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
          {
            diffflux.Multiply(phasemanager.VolFrac(ivolfrac - numfluidphases) *
                                  (1 - phasemanager.Porosity() - phasemanager.SumAddVolFrac()),
                difftensoraddflux, gradscalarnp[iscal]);
          }
          // chemotaxis: scale with volfrac * S * porosity
          else if (phasemanager.ScalarToPhase(iscal).species_type ==
                   MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
          {
            diffflux.Multiply(
                phasemanager.VolFrac(ivolfrac - numfluidphases) * phasemanager.Porosity() *
                    phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID),
                difftensoraddflux, gradscalarnp[iscal]);
          }
          else
            dserror("AddScalarDependentFlux only possible for species in fluid or solid!");

          double v(0.0);
          // TODO: anisotropic difftensor
          // haptotaxis: scale with volfrac * (1 - porosity - sumaddvolfrac)
          if (phasemanager.ScalarToPhase(iscal).species_type ==
              MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
          {
            v = difftensoraddflux(0, 0) * timefacfac / det *
                phasemanager.VolFrac(ivolfrac - numfluidphases) *
                (1 - phasemanager.Porosity() - phasemanager.SumAddVolFrac());
          }
          // chemotaxis: scale with volfrac * S * porosity
          else if (phasemanager.ScalarToPhase(iscal).species_type ==
                   MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
          {
            v = difftensoraddflux(0, 0) * timefacfac / det *
                phasemanager.VolFrac(ivolfrac - numfluidphases) * phasemanager.Porosity() *
                phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID);
          }
          else
            dserror("AddScalarDependentFlux only possible for species in fluid or solid!");

          // gradient of phi w.r.t. reference coordinates
          LINALG::Matrix<nsd, 1> refgradscalarnp(true);
          refgradscalarnp.Multiply(xjm, gradscalarnp[iscal]);

          // 1)
          // -----------------------------------------------------------------------------------------------------------------------------------------
          // OD mesh - diffusive term
          EvaluatorBase<nsd, nen>::CalcDiffODMesh(mymat, deriv, derxy, xjm, diffflux,
              refgradscalarnp, gradscalarnp[iscal], timefacfac, v, numdofpernode, ivolfrac);

          // 2)
          // -----------------------------------------------------------------------------------------------------------------------------------------
          // linearization of porosity
          //-------------dreac/dd = dreac/dporosity * dporosity/dd = dreac/dporosity * dporosity/dJ
          //* dJ/dd = dreac/dporosity * dporosity/dJ * J * N_x
          // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) =
          // det (dx/ds) * ( det(dX/ds) )^-1

          if (phasemanager.PorosityDependsOnStruct())
          {
            static LINALG::Matrix<nsd, 1> diffflux2(true);

            // haptotaxis
            if (phasemanager.ScalarToPhase(iscal).species_type ==
                MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
            {
              diffflux2.Multiply(phasemanager.VolFrac(ivolfrac - numfluidphases) * (-1.0) *
                                     phasemanager.JacobianDefGrad() *
                                     phasemanager.PorosityDerivWrtJacobianDefGrad(),
                  difftensoraddflux, gradscalarnp[iscal]);
            }
            // chemotaxis
            else if (phasemanager.ScalarToPhase(iscal).species_type ==
                     MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
            {
              diffflux2.Multiply(
                  phasemanager.VolFrac(ivolfrac - numfluidphases) *
                      phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID) *
                      phasemanager.JacobianDefGrad() *
                      phasemanager.PorosityDerivWrtJacobianDefGrad(),
                  difftensoraddflux, gradscalarnp[iscal]);
            }
            else
              dserror("AddScalarDependentFlux only possible for species in fluid or solid!");

            for (int vi = 0; vi < nen; ++vi)
            {
              const int fvi = vi * numdofpernode + ivolfrac;
              double laplawf = 0.0;
              for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux2(j);
              const double v = laplawf * timefacfac;

              for (int ui = 0; ui < nen; ++ui)
              {
                for (int idim = 0; idim < nsd; ++idim)
                {
                  const int fui = ui * nsd + idim;
                  mymat(fvi, fui) += v * derxy(idim, ui);
                }
              }
            }
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 08/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracAddFlux<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfrac = numfluidphases; ivolfrac < numfluidphases + numvolfrac; ivolfrac++)
  {
    if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases))
    {
      const int numscal = phasemanager.NumScal();
      // loop over all scalars
      for (int iscal = 0; iscal < numscal; iscal++)
      {
        if (phasemanager.HasAddScalarDependentFlux(ivolfrac - numfluidphases, iscal))
        {
          // diffusion tensor
          LINALG::Matrix<nsd, nsd> difftensoraddflux(true);
          for (int i = 0; i < nsd; i++)
            difftensoraddflux(i, i) = phasemanager.ScalarDiff(ivolfrac - numfluidphases, iscal);

          // haptotaxis: scale with volfrac * (1 - porosity - sumaddvolfrac)
          if (phasemanager.ScalarToPhase(iscal).species_type ==
              MAT::ScatraMatMultiPoro::SpeciesType::species_in_solid)
            difftensoraddflux.Scale(phasemanager.VolFrac(ivolfrac - numfluidphases) *
                                    (1.0 - phasemanager.Porosity() - phasemanager.SumAddVolFrac()));
          // chemotaxis: scale with volfrac * S * porosity
          else if (phasemanager.ScalarToPhase(iscal).species_type ==
                   MAT::ScatraMatMultiPoro::SpeciesType::species_in_fluid)
            difftensoraddflux.Scale(
                phasemanager.VolFrac(ivolfrac - numfluidphases) * phasemanager.Porosity() *
                phasemanager.Saturation(phasemanager.ScalarToPhase(iscal).phaseID));
          else
            dserror("AddScalarDependentFlux only possible for species in fluid or solid!");

          static LINALG::Matrix<nsd, nen> diffflux(true);
          diffflux.Multiply(difftensoraddflux, derxy);

          // diffusive term
          for (int vi = 0; vi < nen; ++vi)
          {
            const int fvi = vi * numdofpernode + ivolfrac;

            for (int ui = 0; ui < nen; ++ui)
            {
              double laplawf(0.0);
              for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j, ui);

              const int fui = ui * numscal + iscal;
              mymat(fvi, fui) += timefacfac * laplawf;
            }
          }

          // additional linearization of receptor kinetic law
          if (phasemanager.HasReceptorKineticLaw(ivolfrac - numfluidphases, iscal))
          {
            // get scalars
            const std::vector<double> scalars = *variablemanager.Scalarnp();
            const std::vector<LINALG::Matrix<nsd, 1>>& gradscalarnp =
                *variablemanager.GradScalarnp();

            // corrrectly scale difftensor
            // d diff / d omega = D_0*(-1.0)*w_half/(w_half+w)^2
            //                    value from above * (-1.0)/(w_half+w)
            difftensoraddflux.Scale(
                -1.0 / (phasemanager.OmegaHalf(ivolfrac - numfluidphases, iscal) + scalars[iscal]));

            static LINALG::Matrix<nsd, 1> diffflux2(true);
            diffflux2.Multiply(difftensoraddflux, gradscalarnp[iscal]);

            for (int vi = 0; vi < nen; ++vi)
            {
              const int fvi = vi * numdofpernode + ivolfrac;

              double laplawf = 0.0;
              for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux2(j);

              for (int ui = 0; ui < nen; ++ui)
              {
                const int fui = ui * numscal + iscal;
                mymat(fvi, fui) += timefacfac * funct(ui) * laplawf;
              }
            }
          }
        }
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureDiff<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    // loop over all volume fraction pressures
    for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
         ivolfracpress++)
    {
      const bool evaluatevolfracpress = variablemanager.ElementHasValidVolFracPressure(
          ivolfracpress - numvolfrac - numfluidphases);

      if (evaluatevolfracpress)
      {
        // get permeability tensor and diffusive flux
        LINALG::Matrix<nsd, nsd> permeabilitytensorvolfracpress(true);
        phasemanager.PermeabilityTensorVolFracPressure(
            ivolfracpress - numfluidphases - numvolfrac, permeabilitytensorvolfracpress);
        permeabilitytensorvolfracpress.Scale(
            1.0 / phasemanager.DynViscosityVolFracPressure(
                      ivolfracpress - numfluidphases - numvolfrac, -1.0));  // TODO: change -1.0

        static LINALG::Matrix<nsd, nen> diffflux(true);
        diffflux.Multiply(permeabilitytensorvolfracpress, derxy);

        // diffusive term
        for (int vi = 0; vi < nen; ++vi)
        {
          const int fvi = vi * numdofpernode + ivolfracpress;

          for (int ui = 0; ui < nen; ++ui)
          {
            double laplawf(0.0);
            for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j, ui);

            const int fui = ui * numdofpernode + ivolfracpress;
            mymat(fvi, fui) += timefacfac * laplawf;
          }
        }

        if (not phasemanager.HasConstantDynViscosityVolFracPressure(
                ivolfracpress - numfluidphases - numvolfrac))
          dserror(
              "only constant dynamic viscosities possible for volume fraction pressures so far");
      }
    }
  }


  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureDiff<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // loop over all volume fractions
  for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
       ivolfracpress++)
  {
    double absgradphi = 0.0;
    for (int idim = 0; idim < nsd; idim++)
    {
      absgradphi += gradphi[ivolfracpress](idim, 0) * gradphi[ivolfracpress](idim, 0);
    }
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfracpress - numvolfrac - numfluidphases);

    if (evaluatevolfracpress)
    {
      // get permeability tensor
      LINALG::Matrix<nsd, nsd> permeabilitytensorvolfracpress(true);
      phasemanager.PermeabilityTensorVolFracPressure(
          ivolfracpress - numfluidphases - numvolfrac, permeabilitytensorvolfracpress);
      permeabilitytensorvolfracpress.Scale(
          1.0 / phasemanager.DynViscosityVolFracPressure(
                    ivolfracpress - numfluidphases - numvolfrac, -1.0));

      static LINALG::Matrix<nsd, 1> diffflux(true);
      diffflux.Multiply(permeabilitytensorvolfracpress, gradphi[ivolfracpress]);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + ivolfracpress;

        // laplacian in weak form
        double laplawf(0.0);
        for (int j = 0; j < nsd; j++) laplawf += derxy(j, vi) * diffflux(j);
        myvec[fvi] -= rhsfac * laplawf;
      }
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureDiff<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  const std::vector<LINALG::Matrix<nsd, 1>>& gradphi = *variablemanager.GradPhinp();

  // loop over all volume fractions
  for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
       ivolfracpress++)
  {
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfracpress - numvolfrac - numfluidphases);

    if (evaluatevolfracpress)
    {
      // get permeability tensor
      LINALG::Matrix<nsd, nsd> permeabilitytensorvolfracpress(true);
      phasemanager.PermeabilityTensorVolFracPressure(
          ivolfracpress - numfluidphases - numvolfrac, permeabilitytensorvolfracpress);
      permeabilitytensorvolfracpress.Scale(
          1.0 / phasemanager.DynViscosityVolFracPressure(
                    ivolfracpress - numfluidphases - numvolfrac, -1.0));

      static LINALG::Matrix<nsd, 1> diffflux(true);
      diffflux.Multiply(permeabilitytensorvolfracpress, gradphi[ivolfracpress]);

      // TODO: anisotropic difftensor
      const double v = permeabilitytensorvolfracpress(0, 0) * timefacfac / det;

      // gradient of phi w.r.t. reference coordinates
      LINALG::Matrix<nsd, 1> refgradphi(true);
      refgradphi.Multiply(xjm, gradphi[ivolfracpress]);

      // OD mesh - diffusive term
      EvaluatorBase<nsd, nen>::CalcDiffODMesh(mymat, deriv, derxy, xjm, diffflux, refgradphi,
          gradphi[ivolfracpress], timefacfac, v, numdofpernode, ivolfracpress);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureDiff<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureReac<nsd,
    nen>::EvaluateMatrixAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, bool inittimederiv)
{
  // we do not need the matrix if we calculate the initial time derivative
  if (!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const int numfluidphases = phasemanager.NumFluidPhases();
    const int numvolfrac = phasemanager.NumVolFrac();

    // loop over all volume fractions
    for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
         ivolfracpress++)
    {
      const bool evaluatevolfracpress = variablemanager.ElementHasValidVolFracPressure(
          ivolfracpress - numvolfrac - numfluidphases);


      if (phasemanager.IsReactive(ivolfracpress) && evaluatevolfracpress)
      {
        double scaledtimefacfac =
            timefacfac / phasemanager.VolFracDensity(ivolfracpress - numfluidphases - numvolfrac);
        //----------------------------------------------------------------
        // reaction terms
        //----------------------------------------------------------------
        for (int vi = 0; vi < nen; ++vi)
        {
          const double v = scaledtimefacfac * funct(vi);
          const int fvi = vi * numdofpernode + ivolfracpress;

          for (int ui = 0; ui < nen; ++ui)
          {
            const double vfunct = v * funct(ui);
            for (int idof = 0; idof < numdofpernode; ++idof)
            {
              const int fui = ui * numdofpernode + idof;

              // rhs ---> -
              mymat(fvi, fui) -= vfunct * phasemanager.ReacDeriv(ivolfracpress, idof);
            }
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureReac<nsd,
    nen>::EvaluateVectorAndAssemble(std::vector<Epetra_SerialDenseVector*>& elevec,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy,
    const LINALG::Matrix<nsd, nen>& xyze, int curphase, int phasetoadd, int numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double rhsfac,
    double fac, bool inittimederiv)
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
       ivolfracpress++)
  {
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfracpress - numvolfrac - numfluidphases);


    if (phasemanager.IsReactive(ivolfracpress) && evaluatevolfracpress)
    {
      double scale = 1.0 / phasemanager.VolFracDensity(ivolfracpress - numfluidphases - numvolfrac);

      double vrhs = scale * rhsfac * phasemanager.ReacTerm(ivolfracpress);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + ivolfracpress;
        // rhs ---> +
        myvec[fvi] += vrhs * funct(vi);
      }
    }
  }

  return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureReac<nsd,
    nen>::EvaluateMatrixODStructAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& deriv,
    const LINALG::Matrix<nsd, nen>& derxy, const LINALG::Matrix<nsd, nsd>& xjm, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac, double det)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();

  // loop over all volume fractions
  for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
       ivolfracpress++)
  {
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfracpress - numvolfrac - numfluidphases);


    if (phasemanager.IsReactive(ivolfracpress) && evaluatevolfracpress)
    {
      // TODO a constant density is assumed here
      double scale = 1.0 / phasemanager.VolFracDensity(ivolfracpress - numfluidphases - numvolfrac);

      double vrhs = scale * timefacfac * phasemanager.ReacTerm(ivolfracpress);

      // linearization of porosity (may appear in reaction term)
      //-------------dreac/dd = dreac/dporosity * dporosity/dd = dreac/dporosity * dporosity/dJ *
      // dJ/dd = dreac/dporosity * dporosity/dJ * J * N_x
      // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det
      // (dx/ds) * ( det(dX/ds) )^-1

      if (phasemanager.PorosityDependsOnStruct())
      {
        vrhs += timefacfac * scale * phasemanager.ReacDerivPorosity(ivolfracpress) *
                phasemanager.JacobianDefGrad() * phasemanager.PorosityDerivWrtJacobianDefGrad();
      }

      // linearization of mesh motion (Jacobian)
      // 1) linearization of fac +
      // 2) possible linearization w.r.t porosity
      // rhs ---> -
      EvaluatorBase<nsd, nen>::CalcLinFacODMesh(
          mymat, funct, derxy, -1.0 * vrhs, numdofpernode, ivolfracpress);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 02/18 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorVolFracPressureReac<nsd,
    nen>::EvaluateMatrixODScatraAndAssemble(std::vector<Epetra_SerialDenseMatrix*>& elemat,
    const LINALG::Matrix<nen, 1>& funct, const LINALG::Matrix<nsd, nen>& derxy, int curphase,
    int phasetoadd, int numdofpernode, const POROFLUIDMANAGER::PhaseManagerInterface& phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd, nen>& variablemanager, double timefacfac,
    double fac)
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const int numfluidphases = phasemanager.NumFluidPhases();
  const int numvolfrac = phasemanager.NumVolFrac();
  const int numscal = phasemanager.NumScal();

  // loop over all volume fractions
  for (int ivolfracpress = numfluidphases + numvolfrac; ivolfracpress < numdofpernode;
       ivolfracpress++)
  {
    const bool evaluatevolfracpress =
        variablemanager.ElementHasValidVolFracPressure(ivolfracpress - numvolfrac - numfluidphases);

    if (phasemanager.IsReactive(ivolfracpress) && evaluatevolfracpress)
    {
      double vrhs = 1.0 / phasemanager.VolFracDensity(ivolfracpress - numfluidphases - numvolfrac) *
                    timefacfac;

      // linearization of reaction term w.r.t scalars
      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi * numdofpernode + ivolfracpress;
        const double v = vrhs * funct(vi);

        for (int ui = 0; ui < nen; ++ui)
        {
          const double vfunct = v * funct(ui);
          for (int iscal = 0; iscal < numscal; ++iscal)
          {
            const int fui = ui * numscal + iscal;
            // rhs ---> -
            mymat(fvi, fui) -= vfunct * phasemanager.ReacDerivScalar(ivolfracpress, iscal);
          }
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
// line 2
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<1, 2>;

// line 3
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<1, 3>;

// 2D elements
// tri3
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<2, 3>;
// quad4
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<2, 4>;

// quad9 and nurbs9
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<2, 9>;

// 3D elements
// hex8
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3, 8>;

// hex27
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3, 27>;
// tet4
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3, 4>;
// tet10
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3, 10>;
// pyramid5
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3, 5>;
