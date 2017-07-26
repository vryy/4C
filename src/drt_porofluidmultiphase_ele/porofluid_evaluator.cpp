/*----------------------------------------------------------------------*/
/*!
 \file porofluid_evaluator.cpp

 \brief helper class for evaluation of the governing equation of multiphase porou flow

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de251
 *----------------------------------------------------------------------*/


#include "porofluid_evaluator.H"

#include "porofluid_variablemanager.H"
#include "porofluid_phasemanager.H"

#include "porofluidmultiphase_ele_parameter.H"

#include "../drt_mat/fluidporo_multiphase.H"
#include "../drt_mat/fluidporo_singlephase.H"

#include "../drt_mat/structporo.H"

/*----------------------------------------------------------------------*
 | factory method                                           vuong 08/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
Teuchos::RCP< DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<nsd,nen> >
DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<nsd,nen>::CreateEvaluator(
    const DRT::ELEMENTS::PoroFluidMultiPhaseEleParameter& para,
    const POROFLUIDMULTIPHASE::Action&                    action,
    int                                                   numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&        phasemanager)
{
  // the evaluator
  Teuchos::RCP<EvaluatorInterface<nsd, nen> > evaluator = Teuchos::null;

  bool inittimederiv = false;
  if(action == POROFLUIDMULTIPHASE::calc_initial_time_deriv)
    inittimederiv = true;
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
    Teuchos::RCP<MultiEvaluator<nsd, nen> > evaluator_multiphase
    = Teuchos::rcp(new MultiEvaluator<nsd, nen>());

    // build evaluators for all but last phase
    for (int curphase = 0; curphase < numdofpernode-1; curphase++)
    {
      // initialize the evaluator for the current phase
      Teuchos::RCP<MultiEvaluator<nsd, nen> > evaluator_phase =
          Teuchos::rcp(new MultiEvaluator<nsd, nen>());

      // temporary interfaces
      Teuchos::RCP<EvaluatorInterface<nsd, nen> > tmpevaluator = Teuchos::null;
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

      // Note: this term cancels because of the formulation w.r.t. the material formulation of the solid
      // add evaluator for the conservative term (w, v \nabla \cdot S )
      //assembler = Teuchos::rcp(new AssembleStandard(curphase,inittimederiv));
      //tmpevaluator = Teuchos::rcp(new EvaluatorConv<nsd, nen>(assembler,curphase));
      //evaluator_phase->AddEvaluator(tmpevaluator);

      // add evaluator for the convective conservative term (w, S \nabla \cdot v )
      if (para.IsAle())
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase,inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorSatDivVel<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);
      }
      // add evaluator for Biot stabilization
      if (para.BiotStab())
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase,inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorBiotStab<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);
      }

      // add evaluator for the diffusive term (\nabla w, K \nabla p)
      // the diffusive term is also assembled into the last phase
      assembler = Teuchos::rcp(new AssembleAlsoIntoOtherPhase(curphase, numdofpernode-1, inittimederiv));
      tmpevaluator = Teuchos::rcp(new EvaluatorDiff<nsd, nen>(assembler,curphase));
      evaluator_phase->AddEvaluator(tmpevaluator);

      // add evaluator for the reactive term
      if (phasemanager.IsReactive(curphase))
      {
        // the reactive term is also assembled into the last phase
        assembler = Teuchos::rcp(
            new AssembleAlsoIntoOtherPhase(curphase, numdofpernode-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorReac<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);
      }

      // add evaluators for the instationary terms
      if (not para.IsStationary())
      {
        // add evaluator for the instationary pressure term
        // the term is also assembled into the last phase
        assembler = Teuchos::rcp(
            new AssembleAlsoIntoOtherPhase(curphase, numdofpernode-1, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassPressure<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);

        // add evaluator for the instationary solid pressure term
        assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassSolidPressureSat<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);

        // add evaluator for the instationary saturation term
        assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassSaturation<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);
      }

      // add the evaluator of the phase to the multiphase evaluator
      evaluator_multiphase->AddEvaluator(evaluator_phase);

    }

    // build evaluators for the last phase
    {
      const int curphase = numdofpernode-1;

      // initialize the evaluator for the last phase
      Teuchos::RCP<MultiEvaluator<nsd, nen> > evaluator_lastphase =
          Teuchos::rcp(new MultiEvaluator<nsd, nen>());

      // temporary interfaces
      Teuchos::RCP<EvaluatorInterface<nsd, nen> > tmpevaluator = Teuchos::null;
      Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

      // add evaluator for the convective conservative term (w, \nabla \cdot v )
      if (para.IsAle())
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase, false));
        tmpevaluator = Teuchos::rcp(new EvaluatorDivVel<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);
      }
      // add evaluator for Biot stabilization
      if (para.BiotStab())
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase,inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorBiotStab<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);
      }

      // add evaluator for the diffusive term (\nabla w, K \nabla p)
      assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
      tmpevaluator = Teuchos::rcp(new EvaluatorDiff<nsd, nen>(assembler,curphase));
      evaluator_lastphase->AddEvaluator(tmpevaluator);

      // add evaluator for the reactive term
      if (phasemanager.IsReactive(curphase))
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorReac<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);
      }

      // add evaluators for the instationary terms
      if (not para.IsStationary())
      {
        // add evaluator for the instationary pressure term
        assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassPressure<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);

        // add evaluator for the instationary solid pressure term
        assembler = Teuchos::rcp(new AssembleStandard(curphase, inittimederiv));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassSolidPressure<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);
      }

      // add the evaluator of the phase to the multiphase evaluator
      evaluator_multiphase->AddEvaluator(evaluator_lastphase);
    }

    evaluator = evaluator_multiphase;
    break;
  }
  case POROFLUIDMULTIPHASE::calc_pres_and_sat:
  {
    // initialize the evaluator for the multi phase element
    Teuchos::RCP<MultiEvaluator<nsd, nen> > evaluator_multiphase
    = Teuchos::rcp(new MultiEvaluator<nsd, nen>());

    // temporary interfaces
    Teuchos::RCP<EvaluatorInterface<nsd, nen> > tmpevaluator = Teuchos::null;

    // initialize temporary assembler
    Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

    // build evaluators for all but last phase
    for (int iphase = 0; iphase < numdofpernode; iphase++)
    {
      assembler = Teuchos::rcp(new AssembleStandard(iphase, false));
      tmpevaluator = Teuchos::rcp(new EvaluatorPressureAndSaturation<nsd, nen>(assembler,iphase));
      evaluator_multiphase->AddEvaluator(tmpevaluator);
    }
    evaluator = evaluator_multiphase;

    break;
  }
  case POROFLUIDMULTIPHASE::calc_solidpressure:
  {
    Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1, false));
    evaluator = Teuchos::rcp(new EvaluatorSolidPressure<nsd, nen>(assembler,-1));

    break;
  }
  case POROFLUIDMULTIPHASE::calc_porosity:
  {
    Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1, false));
    evaluator = Teuchos::rcp(new EvaluatorPorosity<nsd, nen>(assembler,-1));

    break;
  }
  case POROFLUIDMULTIPHASE::recon_flux_at_nodes:
  {
    // initialize the evaluator for the multi phase element
    Teuchos::RCP<MultiEvaluator<nsd, nen> > evaluator_multiphase
    = Teuchos::rcp(new MultiEvaluator<nsd, nen>());

    // temporary interfaces
    Teuchos::RCP<EvaluatorInterface<nsd, nen> > tmpevaluator = Teuchos::null;

    // initialize temporary assembler
    Teuchos::RCP<AssembleInterface> assembler = Teuchos::null;

    assembler = Teuchos::rcp(new AssembleStandard(-1, false));
    tmpevaluator = Teuchos::rcp(new ReconstructFluxLinearization<nsd, nen>(assembler,-1));
    evaluator_multiphase->AddEvaluator(tmpevaluator);

    // build evaluators for all phases
    for (int iphase = 0; iphase < numdofpernode; iphase++)
    {
      assembler = Teuchos::rcp(new AssembleStandard(iphase, false));
      tmpevaluator = Teuchos::rcp(new ReconstructFluxRHS<nsd, nen>(assembler,iphase));
      evaluator_multiphase->AddEvaluator(tmpevaluator);
    }
    evaluator = evaluator_multiphase;

    break;
  }
  default:
  {
    dserror("unknown action for evaluation class!");
    break;
  }
} // switch(action)

// done
  return evaluator;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

   // convective term in convective form
   /*
        /                               \
       |                                 |
       | prefac * v * nabla * Dphi  , q  |
       |                                 |
        \                               /
   */
   const double prefac = timefacfac;

   LINALG::Matrix<nen,1>    conv;
   // convective part in convective form: rho*u_x*N,x+ rho*u_y*N,y
   conv.MultiplyTN(derxy,*variablemanager.ConVelnp());

   for (int vi=0; vi<nen; ++vi)
   {
     const double v = prefac*funct(vi);
     const int fvi = vi*numdofpernode+phasetoadd;

     for (int ui=0; ui<nen; ++ui)
       for (int idof=0; idof<numdofpernode; ++idof)
       {
         const int fui = ui*numdofpernode+idof;

         mymat(fvi,fui) += v*conv(ui)*phasemanager.SaturationDeriv(curphase,idof);
       }
   }
   return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  double conv_sat = 0.0;
  for (int idof=0; idof<numdofpernode; ++idof)
  {
    // convective term
    const double conv_phi = variablemanager.ConVelnp()->Dot((*variablemanager.GradPhinp())[idof]);
    conv_sat += rhsfac*phasemanager.SaturationDeriv(curphase,idof)*conv_phi;
  }
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    myvec[fvi] -= conv_sat*funct(vi);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorConv<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  double vrhs = rhsfac*variablemanager.DivConVelnp();

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    myvec[fvi] -= vrhs*funct(vi);
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  // d (div v_s)/d d_n+1 = derxy * 1.0/theta/dt * d_n+1
  // prefactor is fac since timefacfac/theta/dt = fac
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    const double v = fac*funct(vi);

    for (int ui=0; ui<nen; ++ui)
    {
      for (int idim=0; idim<nsd; ++idim)
      {
        const int fui = ui*nsd+idim;

        mymat(fvi,fui) += v*derxy(idim,ui);
      }
    }
  }

  // shapederivatives see fluid_ele_calc_poro.cpp
  LINALG::Matrix<nsd,nsd> gridvelderiv(true);
  gridvelderiv.MultiplyNT(*(variablemanager.EConVelnp()),deriv);

  if(nsd == 3)
  {
    const double gridvelderiv_0_0   = gridvelderiv(0, 0);
    const double gridvelderiv_0_1   = gridvelderiv(0, 1);
    const double gridvelderiv_0_2   = gridvelderiv(0, 2);
    const double gridvelderiv_1_0   = gridvelderiv(1, 0);
    const double gridvelderiv_1_1   = gridvelderiv(1, 1);
    const double gridvelderiv_1_2   = gridvelderiv(1, 2);
    const double gridvelderiv_2_0   = gridvelderiv(2, 0);
    const double gridvelderiv_2_1   = gridvelderiv(2, 1);
    const double gridvelderiv_2_2   = gridvelderiv(2, 2);

    const double xjm_0_0   = xjm(0, 0);
    const double xjm_0_1   = xjm(0, 1);
    const double xjm_0_2   = xjm(0, 2);
    const double xjm_1_0   = xjm(1, 0);
    const double xjm_1_1   = xjm(1, 1);
    const double xjm_1_2   = xjm(1, 2);
    const double xjm_2_0   = xjm(2, 0);
    const double xjm_2_1   = xjm(2, 1);
    const double xjm_2_2   = xjm(2, 2);

#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

#define derxjm_001(ui) (deriv(2, ui)*xjm_1_2 - deriv(1, ui)*xjm_2_2)
#define derxjm_002(ui) (deriv(1, ui)*xjm_2_1 - deriv(2, ui)*xjm_1_1)

#define derxjm_100(ui) (deriv(1, ui)*xjm_2_2 - deriv(2, ui)*xjm_1_2)
#define derxjm_102(ui) (deriv(2, ui)*xjm_1_0 - deriv(1, ui)*xjm_2_0)

#define derxjm_200(ui) (deriv(2, ui)*xjm_1_1 - deriv(1, ui)*xjm_2_1)
#define derxjm_201(ui) (deriv(1, ui)*xjm_2_0 - deriv(2, ui)*xjm_1_0)

#define derxjm_011(ui) (deriv(0, ui)*xjm_2_2 - deriv(2, ui)*xjm_0_2)
#define derxjm_012(ui) (deriv(2, ui)*xjm_0_1 - deriv(0, ui)*xjm_2_1)

#define derxjm_110(ui) (deriv(2, ui)*xjm_0_2 - deriv(0, ui)*xjm_2_2)
#define derxjm_112(ui) (deriv(0, ui)*xjm_2_0 - deriv(2, ui)*xjm_0_0)

#define derxjm_210(ui) (deriv(0, ui)*xjm_2_1 - deriv(2, ui)*xjm_0_1)
#define derxjm_211(ui) (deriv(2, ui)*xjm_0_0 - deriv(0, ui)*xjm_2_0)

#define derxjm_021(ui) (deriv(1, ui)*xjm_0_2 - deriv(0, ui)*xjm_1_2)
#define derxjm_022(ui) (deriv(0, ui)*xjm_1_1 - deriv(1, ui)*xjm_0_1)

#define derxjm_120(ui) (deriv(0, ui)*xjm_1_2 - deriv(1, ui)*xjm_0_2)
#define derxjm_122(ui) (deriv(1, ui)*xjm_0_0 - deriv(0, ui)*xjm_1_0)

#define derxjm_220(ui) (deriv(1, ui)*xjm_0_1 - deriv(0, ui)*xjm_1_1)
#define derxjm_221(ui) (deriv(0, ui)*xjm_1_0 - deriv(1, ui)*xjm_0_0)

    for (int ui = 0; ui < nen; ++ui)
    {

      const double v0 =    gridvelderiv_1_0 * derxjm_(0,0,1,ui)
                         + gridvelderiv_1_1 * derxjm_(0,1,1,ui)
                         + gridvelderiv_1_2 * derxjm_(0,2,1,ui)
                         + gridvelderiv_2_0 * derxjm_(0,0,2,ui)
                         + gridvelderiv_2_1 * derxjm_(0,1,2,ui)
                         + gridvelderiv_2_2 * derxjm_(0,2,2,ui);

      const double v1 =    gridvelderiv_0_0 * derxjm_(1,0,0,ui)
                         + gridvelderiv_0_1 * derxjm_(1,1,0,ui)
                         + gridvelderiv_0_2 * derxjm_(1,2,0,ui)
                         + gridvelderiv_2_0 * derxjm_(1,0,2,ui)
                         + gridvelderiv_2_1 * derxjm_(1,1,2,ui)
                         + gridvelderiv_2_2 * derxjm_(1,2,2,ui);

      const double v2 =    gridvelderiv_0_0 * derxjm_(2,0,0,ui)
                         + gridvelderiv_0_1 * derxjm_(2,1,0,ui)
                         + gridvelderiv_0_2 * derxjm_(2,2,0,ui)
                         + gridvelderiv_1_0 * derxjm_(2,0,1,ui)
                         + gridvelderiv_1_1 * derxjm_(2,1,1,ui)
                         + gridvelderiv_1_2 * derxjm_(2,2,1,ui);

      for (int vi = 0; vi < nen; ++vi)
      {
        const int fvi = vi*numdofpernode+phasetoadd;
        const double v = timefacfac / det * funct(vi);

        mymat(fvi, ui * 3 + 0) += v * v0;

        mymat(fvi, ui * 3 + 1) += v * v1;

        mymat(fvi, ui * 3 + 2) += v * v2;
      }
    }

  }
  else if(nsd == 2)
  {

    const double gridvelderiv_0_0   = gridvelderiv(0, 0);
    const double gridvelderiv_0_1   = gridvelderiv(0, 1);
    const double gridvelderiv_1_0   = gridvelderiv(1, 0);
    const double gridvelderiv_1_1   = gridvelderiv(1, 1);

    for (int vi = 0; vi < nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;
      const double v = timefacfac / det * funct(vi);
      for (int ui = 0; ui < nen; ++ui)
      {
        mymat(fvi, ui * 2    ) += v * (- gridvelderiv_1_0 * deriv(1,ui)
                                       + gridvelderiv_1_1 * deriv(0,ui));

        mymat(fvi, ui * 2 + 1) += v * (+ gridvelderiv_0_0 * deriv(1,ui)
                                       - gridvelderiv_0_1 * deriv(0,ui));
      }
    }
  }
  else
    dserror("shapederivatives not implemented for 1D!");

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDivVel<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // call base class
  EvaluatorDivVel<nsd,nen>::EvaluateMatrixAndAssemble(
      elemat,
      funct,
      derxy,
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      timefacfac,
      fac,
      inittimederiv
    );

  // no linearization needed in case of initial time derivative calculation
  if(!inittimederiv)
  {
  // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const double consfac = timefacfac*variablemanager.DivConVelnp();
    for (int vi=0; vi<nen; ++vi)
    {
      const double v = consfac*funct(vi);
      const int fvi = vi*numdofpernode+phasetoadd;

      for (int ui=0; ui<nen; ++ui)
        for (int idof=0; idof<numdofpernode; ++idof)
        {
          const int fui = ui*numdofpernode+idof;

          mymat(fvi,fui) += v*phasemanager.SaturationDeriv(curphase,idof)*funct(ui);
        }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // call base class with scaled factors
  EvaluatorDivVel<nsd,nen>::EvaluateVectorAndAssemble(
      elevec,
      funct,
      derxy,
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      phasemanager.Saturation(curphase)*rhsfac,
      phasemanager.Saturation(curphase)*fac,
      inittimederiv);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  // call base class with scaled factors
  EvaluatorDivVel<nsd,nen>::EvaluateMatrixODStructAndAssemble(
      elemat,
      funct,
      deriv,
      derxy,
      xjm,
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      phasemanager.Saturation(curphase)*timefacfac,
      phasemanager.Saturation(curphase)*fac,
      det);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSatDivVel<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{

  dserror("Biot stabilization is still missing");
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  dserror("Biot stabilization is still missing");

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  dserror("Biot stabilization is still missing");
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorBiotStab<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // we do not need the matrix if we calculate the initial time derivative
  if(!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    const std::vector<LINALG::Matrix<nsd,1> >& gradphi = *variablemanager.GradPhinp();

    // current pressure gradient
    LINALG::Matrix<nsd,1> gradpres(true);
    gradpres.Clear();

    // compute the pressure gradient from the phi gradients
    for (int idof=0; idof<numdofpernode; ++idof)
      gradpres.Update(phasemanager.PressureDeriv(curphase,idof),gradphi[idof],1.0);

    double abspressgrad=0.0;
    for (int i=0; i<nsd; i++)
      abspressgrad+=gradpres(i)*gradpres(i);
    abspressgrad=sqrt(abspressgrad);

    // permeability tensor
    LINALG::Matrix<nsd,nsd> permeabilitytensor(true);
    phasemanager.PermeabilityTensor(curphase,permeabilitytensor);

    static LINALG::Matrix<nsd,nen> diffflux(true);
    diffflux.Multiply(permeabilitytensor,derxy);
    diffflux.Scale(phasemanager.RelPermeability(curphase)/phasemanager.DynViscosity(curphase, abspressgrad));

    // diffusive term
    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;

      for (int ui=0; ui<nen; ++ui)
      {
        double laplawf(0.0);
        for (int j = 0; j<nsd; j++)
          laplawf += derxy(j, vi)*diffflux(j, ui);
        for (int idof=0; idof<numdofpernode; ++idof)
        {
          const int fui = ui*numdofpernode+idof;
          mymat(fvi,fui) += timefacfac*laplawf*phasemanager.PressureDeriv(curphase,idof);
          //std::cout << timefacfac*laplawf*phasemanager.PressureDeriv(curphase,idof) << std::endl;
        }
      }
    }
    //----------------------------------------------------------------
    // linearization of relative permeability w.r.t. dof
    //----------------------------------------------------------------
    if(not phasemanager.HasConstantRelPermeability(curphase))
    {
      static LINALG::Matrix<nsd,1> diffflux2(true);
      diffflux2.Multiply(permeabilitytensor,gradpres);
      diffflux2.Scale(phasemanager.RelPermeabilityDeriv(curphase)/phasemanager.DynViscosity(curphase, abspressgrad));

      for (int vi=0; vi<nen; ++vi)
      {
        const int fvi = vi*numdofpernode+phasetoadd;
        double laplawf = 0.0;
        for (int j = 0; j<nsd; j++)
          laplawf += derxy(j, vi)*diffflux2(j);
        for (int ui=0; ui<nen; ++ui)
        {
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;
            mymat(fvi,fui) += timefacfac*funct(ui)*laplawf*phasemanager.SaturationDeriv(curphase,idof);
          }
        }
      }
    }
    //----------------------------------------------------------------
    // linearization of dynamic viscosity w.r.t. dof
    //----------------------------------------------------------------
    if(not phasemanager.HasConstantDynViscosity(curphase))
    {
      // derivative of abspressgrad w.r.t. pressure gradient
      static LINALG::Matrix<nsd,1> dabspressgraddpresgradp(true);
      // avoid division by zero
      if(abspressgrad > 1.0e-12)
        for (int i = 0; i < nsd; i++)
          dabspressgraddpresgradp(i) = gradpres(i)/abspressgrad;

      static LINALG::Matrix<nsd,1> diffflux2(true);
      diffflux2.Multiply(permeabilitytensor,gradpres);
      // d (1/visc) / d abspressgrad = -1.0 * visc^(-2) * d visc / d abspressgrad
      diffflux2.Scale(-1.0*phasemanager.RelPermeability(curphase)/phasemanager.DynViscosity(curphase, abspressgrad)
          /phasemanager.DynViscosity(curphase, abspressgrad)*phasemanager.DynViscosityDeriv(curphase, abspressgrad));

      for (int vi=0; vi<nen; ++vi)
      {
        const int fvi = vi*numdofpernode+phasetoadd;
        double laplawf = 0.0;
        for (int j = 0; j<nsd; j++)
          laplawf += derxy(j, vi)*diffflux2(j);
        for (int ui=0; ui<nen; ++ui)
        {
          double gradpderxy = 0.0;
          for (int j = 0; j<nsd; j++)
            gradpderxy += derxy(j, ui)*dabspressgraddpresgradp(j);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;
            // d abspressgrad / d phi = d abspressgrad / d gradp * d gradp / d phi =
            //                        = d abspressgrad / d gradp * d / d phi( d p / d phi * d phi / d x) =
            //                        = d abspressgrad / d gradp * derxy * d p / d phi
            // Note: FD-Check might fail here due to kink in formulation of cell-adherence model-law
            mymat(fvi,fui) += timefacfac*laplawf*phasemanager.PressureDeriv(curphase,idof)*gradpderxy;
          }
        }
      }
    }
  } // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const std::vector<LINALG::Matrix<nsd,1> >& gradphi = *variablemanager.GradPhinp();

  // current pressure gradient
  LINALG::Matrix<nsd,1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof=0; idof<numdofpernode; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase,idof),gradphi[idof],1.0);

  double abspressgrad=0.0;
  for (int i=0; i<nsd; i++)
    abspressgrad+=gradpres(i)*gradpres(i);
  abspressgrad=sqrt(abspressgrad);

  // diffusion tensor
  LINALG::Matrix<nsd,nsd> difftensor(true);
  phasemanager.PermeabilityTensor(curphase,difftensor);
  difftensor.Scale(phasemanager.RelPermeability(curphase)/phasemanager.DynViscosity(curphase, abspressgrad));

  static LINALG::Matrix<nsd,1> diffflux(true);
  diffflux.Multiply(difftensor,gradpres);

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    // laplacian in weak form
    double laplawf(0.0);
    for (int j = 0; j<nsd; j++)
      laplawf += derxy(j,vi)*diffflux(j);
    myvec[fvi] -= rhsfac*laplawf;
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  const std::vector<LINALG::Matrix<nsd,1> >& gradphi = *variablemanager.GradPhinp();

  // current pressure gradient
  LINALG::Matrix<nsd,1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof=0; idof<numdofpernode; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase,idof),gradphi[idof],1.0);

  double abspressgrad=0.0;
  for (int i=0; i<nsd; i++)
    abspressgrad+=gradpres(i)*gradpres(i);
  abspressgrad=sqrt(abspressgrad);

  // diffusion tensor
  LINALG::Matrix<nsd,nsd> difftensor(true);
  phasemanager.PermeabilityTensor(curphase,difftensor);
  difftensor.Scale(phasemanager.RelPermeability(curphase)/phasemanager.DynViscosity(curphase, abspressgrad));

  static LINALG::Matrix<nsd,1> diffflux(true);
  diffflux.Multiply(difftensor,gradpres);

  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space, i.e. det(dx/ds)
  // in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac * N_x
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    // laplacian in weak form
    double laplawf(0.0);
    for (int j = 0; j<nsd; j++)
      laplawf += derxy(j,vi)*diffflux(j);
    double v = - laplawf*timefacfac;

    for (int ui=0; ui<nen; ++ui)
    {
      for (int idim=0; idim<nsd; ++idim)
      {
        const int fui = ui*nsd+idim;
        mymat(fvi,fui) += v*derxy(idim,ui);
      }
    }

  }

  //----------------------------------------------------------------
  // standard Galerkin terms  -- "shapederivatives" diffusive term
  //----------------------------------------------------------------
  // see scatra_ele_calc_OD.cpp

  // TODO: anisotropic difftensor and
  //       non-constant viscosity (because of pressure gradient, probably not really necessary)
  const double v = difftensor(0,0)*timefacfac/det;

  //gradient of pressure w.r.t. reference coordinates
  LINALG::Matrix<nsd,1> refgradpres(true);
  refgradpres.Clear();

    //gradient of phi w.r.t. reference coordinates
  std::vector<LINALG::Matrix<nsd,1> > refgradphi(numdofpernode,LINALG::Matrix<nsd,1>(true)); //static LINALG::Matrix<nsd,1> refgradphi;
  for (int idof=0; idof<numdofpernode; ++idof)
    refgradphi[idof].Multiply(xjm,gradphi[idof]);

    // compute the pressure gradient from the phi gradients
    for (int idof=0; idof<numdofpernode; ++idof)
      refgradpres.Update(phasemanager.PressureDeriv(curphase,idof),refgradphi[idof],1.0);

  if(nsd == 3)
  {
    const double xjm_0_0   = xjm(0, 0);
    const double xjm_0_1   = xjm(0, 1);
    const double xjm_0_2   = xjm(0, 2);
    const double xjm_1_0   = xjm(1, 0);
    const double xjm_1_1   = xjm(1, 1);
    const double xjm_1_2   = xjm(1, 2);
    const double xjm_2_0   = xjm(2, 0);
    const double xjm_2_1   = xjm(2, 1);
    const double xjm_2_2   = xjm(2, 2);

    const double gradpres_0   = gradpres(0);
    const double gradpres_1   = gradpres(1);
    const double gradpres_2   = gradpres(2);

    for (unsigned vi = 0; vi < nen; ++vi)
    {
      const double deriv_vi_0   = deriv(0,vi);
      const double deriv_vi_1   = deriv(1,vi);
      const double deriv_vi_2   = deriv(2,vi);

      const int fvi = vi*numdofpernode+phasetoadd;

      for (unsigned ui = 0; ui < nen; ++ui)
      {
        const double v00 = + gradpres_1 * (
                                              deriv_vi_0 * (deriv(2, ui)*xjm_1_2 - deriv(1, ui)*xjm_2_2)
                                            + deriv_vi_1 * (deriv(0, ui)*xjm_2_2 - deriv(2, ui)*xjm_0_2)
                                            + deriv_vi_2 * (deriv(1, ui)*xjm_0_2 - deriv(0, ui)*xjm_1_2)
                                         )
                           + gradpres_2 * (
                                              deriv_vi_0 * (deriv(1, ui)*xjm_2_1 - deriv(2, ui)*xjm_1_1)
                                            + deriv_vi_1 * (deriv(2, ui)*xjm_0_1 - deriv(0, ui)*xjm_2_1)
                                            + deriv_vi_2 * (deriv(0, ui)*xjm_1_1 - deriv(1, ui)*xjm_0_1)
                                         );
        const double v01 = + gradpres_0 * (
                                              deriv_vi_0 * (deriv(1, ui)*xjm_2_2 - deriv(2, ui)*xjm_1_2)
                                            + deriv_vi_1 * (deriv(2, ui)*xjm_0_2 - deriv(0, ui)*xjm_2_2)
                                            + deriv_vi_2 * (deriv(0, ui)*xjm_1_2 - deriv(1, ui)*xjm_0_2)
                                         )
                           + gradpres_2 * (
                                              deriv_vi_0 * (deriv(2, ui)*xjm_1_0 - deriv(1, ui)*xjm_2_0)
                                            + deriv_vi_1 * (deriv(0, ui)*xjm_2_0 - deriv(2, ui)*xjm_0_0)
                                            + deriv_vi_2 * (deriv(1, ui)*xjm_0_0 - deriv(0, ui)*xjm_1_0)
        );
        const double v02 = + gradpres_0 * (
                                              deriv_vi_0 * (deriv(2, ui)*xjm_1_1 - deriv(1, ui)*xjm_2_1)
                                            + deriv_vi_1 * (deriv(0, ui)*xjm_2_1 - deriv(2, ui)*xjm_0_1)
                                            + deriv_vi_2 * (deriv(1, ui)*xjm_0_1 - deriv(0, ui)*xjm_1_1)
                                         )
                           + gradpres_1 * (
                                              deriv_vi_0 * (deriv(1, ui)*xjm_2_0 - deriv(2, ui)*xjm_1_0)
                                            + deriv_vi_1 * (deriv(2, ui)*xjm_0_0 - deriv(0, ui)*xjm_2_0)
                                            + deriv_vi_2 * (deriv(0, ui)*xjm_1_0 - deriv(1, ui)*xjm_0_0)
                                         );

        mymat(fvi, ui * nsd + 0) += v * v00;
        mymat(fvi, ui * nsd + 1) += v * v01;
        mymat(fvi, ui * nsd + 2) += v * v02;
      }
    }

    const double refgradpres_0   = refgradpres(0);
    const double refgradpres_1   = refgradpres(1);
    const double refgradpres_2   = refgradpres(2);

    for (unsigned vi = 0; vi < nen; ++vi)
    {
      const double derxy_vi_0   = derxy(0,vi);
      const double derxy_vi_1   = derxy(1,vi);
      const double derxy_vi_2   = derxy(2,vi);

      const int fvi = vi*numdofpernode+phasetoadd;

      for (unsigned ui = 0; ui < nen; ++ui)
      {
        const double v00 = + derxy_vi_1  * (
                                                refgradpres_0 * (deriv(2, ui)*xjm_1_2 - deriv(1, ui)*xjm_2_2)
                                              + refgradpres_1 * (deriv(0, ui)*xjm_2_2 - deriv(2, ui)*xjm_0_2)
                                              + refgradpres_2 * (deriv(1, ui)*xjm_0_2 - deriv(0, ui)*xjm_1_2)
                                           )
                           + derxy_vi_2  * (
                                                refgradpres_0 * (deriv(1, ui)*xjm_2_1 - deriv(2, ui)*xjm_1_1)
                                              + refgradpres_1 * (deriv(2, ui)*xjm_0_1 - deriv(0, ui)*xjm_2_1)
                                              + refgradpres_2 * (deriv(0, ui)*xjm_1_1 - deriv(1, ui)*xjm_0_1)
                                           );
        const double v01 = + derxy_vi_0  * (
                                                refgradpres_0 * (deriv(1, ui)*xjm_2_2 - deriv(2, ui)*xjm_1_2)
                                              + refgradpres_1 * (deriv(2, ui)*xjm_0_2 - deriv(0, ui)*xjm_2_2)
                                              + refgradpres_2 * (deriv(0, ui)*xjm_1_2 - deriv(1, ui)*xjm_0_2)
                                           )
                           + derxy_vi_2  * (
                                                refgradpres_0 * (deriv(2, ui)*xjm_1_0 - deriv(1, ui)*xjm_2_0)
                                              + refgradpres_1 * (deriv(0, ui)*xjm_2_0 - deriv(2, ui)*xjm_0_0)
                                              + refgradpres_2 * (deriv(1, ui)*xjm_0_0 - deriv(0, ui)*xjm_1_0)
                                           );
        const double v02 = + derxy_vi_0  * (
                                                refgradpres_0 * (deriv(2, ui)*xjm_1_1 - deriv(1, ui)*xjm_2_1)
                                              + refgradpres_1 * (deriv(0, ui)*xjm_2_1 - deriv(2, ui)*xjm_0_1)
                                              + refgradpres_2 * (deriv(1, ui)*xjm_0_1 - deriv(0, ui)*xjm_1_1)
                                           )
                           + derxy_vi_1  * (
                                                refgradpres_0 * (deriv(1, ui)*xjm_2_0 - deriv(2, ui)*xjm_1_0)
                                              + refgradpres_1 * (deriv(2, ui)*xjm_0_0 - deriv(0, ui)*xjm_2_0)
                                              + refgradpres_2 * (deriv(0, ui)*xjm_1_0 - deriv(1, ui)*xjm_0_0)
                                           );

        mymat(fvi, ui * nsd + 0) += v * v00;
        mymat(fvi, ui * nsd + 1) += v * v01;
        mymat(fvi, ui * nsd + 2) += v * v02;
      }
    }
  }
  else if(nsd == 2)
  {
    {

      const double gradpres_0   = gradpres(0);
      const double gradpres_1   = gradpres(1);

      for (int vi = 0; vi < nen; ++vi)
      {
        const double deriv_vi_0   = deriv(0,vi);
        const double deriv_vi_1   = deriv(1,vi);

        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui = 0; ui < nen; ++ui)
        {
          const double v00 = + gradpres_1 * (
                                              - deriv_vi_0 * deriv(1, ui)
                                              + deriv_vi_1 * deriv(0, ui)
                                           );
          const double v01 = + gradpres_0 * (
                                                deriv_vi_0 * deriv(1, ui)
                                              - deriv_vi_1 * deriv(0, ui)
                                           );

          mymat(fvi, ui * nsd + 0) += v * v00;
          mymat(fvi, ui * nsd + 1) += v * v01;
        }
      }
    }

      const double refgradpres_0   = refgradpres(0);
      const double refgradpres_1   = refgradpres(1);

      for (unsigned vi = 0; vi < nen; ++vi)
      {
        const double derxy_vi_0   = derxy(0,vi);
        const double derxy_vi_1   = derxy(1,vi);

        const int fvi = vi*numdofpernode+phasetoadd;

        for (unsigned ui = 0; ui < nen; ++ui)
        {
          const double v00 = + derxy_vi_1  * (
                                                - refgradpres_0 * deriv(1, ui)
                                                + refgradpres_1 * deriv(0, ui)
                                             );
          const double v01 = + derxy_vi_0  * (
                                                  refgradpres_0 * deriv(1, ui)
                                                - refgradpres_1 * deriv(0, ui)
                                             );

          mymat(fvi, ui * nsd + 0) += v * v00;
          mymat(fvi, ui * nsd + 1) += v * v01;
        }
      }
    }
    else
      dserror("shapederivatives not implemented for 1D!");

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorDiff<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // we do not need the matrix if we calculate the initial time derivative
  if(!inittimederiv)
  {
    // get matrix to fill
    Epetra_SerialDenseMatrix& mymat = *elemat[0];

    // TODO a constant density is assumed here
    double scaledtimefacfac =timefacfac/phasemanager.Density(curphase);

    //----------------------------------------------------------------
    // reaction terms
    //----------------------------------------------------------------
    for (int vi=0; vi<nen; ++vi)
    {
      const double v = scaledtimefacfac*funct(vi);
      const int fvi = vi*numdofpernode+phasetoadd;

      for (int ui=0; ui<nen; ++ui)
      {
        const double vfunct = v*funct(ui);
        for (int idof=0; idof<numdofpernode; ++idof)
        {
          const int fui = ui*numdofpernode+idof;

          mymat(fvi,fui) += vfunct*phasemanager.ReacDeriv(curphase,idof);
        }
      }
    }
  } // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  if(not phasemanager.IsReactive(curphase))
    return;

  // TODO a constant density is assumed here
  double scale = 1.0/phasemanager.Density(curphase);

  double vrhs = scale*rhsfac*phasemanager.ReacTerm(curphase);

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    myvec[fvi] -= vrhs*funct(vi);
  }
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  if(not phasemanager.IsReactive(curphase))
    return;

  double scale = 1.0/phasemanager.Density(curphase);
  double vrhs = 0.0;

  // linearization of porosity (may appear in reaction term)
  //-------------dreac/dd = dreac/dporosity * dporosity/dd = dreac/dporosity * dporosity/dJ * dJ/dd = dreac/dporosity * dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1

  if(phasemanager.PorosityDependsOnStruct())
  {
    vrhs = timefacfac*scale*phasemanager.ReacDerivPorosity(curphase)*phasemanager.JacobianDefGrad()
        *phasemanager.PorosityDerivWrtJacobianDefGrad();

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;
      const double v = funct(vi)*vrhs;

      for (int ui=0; ui<nen; ++ui)
      {
        for (int idim=0; idim<nsd; ++idim)
        {
          const int fui = ui*nsd+idim;

          mymat(fvi,fui) += v*derxy(idim,ui);
        }
      }
    }
  }

  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space, i.e. det(dx/ds)
  // in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac * N_x
  // TODO a constant density is assumed here

  vrhs = scale*timefacfac*phasemanager.ReacTerm(curphase);

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    const double v = vrhs*funct(vi);

    for (int ui=0; ui<nen; ++ui)
    {
      for (int idim=0; idim<nsd; ++idim)
      {
        const int fui = ui*nsd+idim;
        mymat(fvi,fui) += v*derxy(idim,ui);
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorReac<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  if(not phasemanager.IsReactive(curphase))
    return;

  const int numscal = phasemanager.NumScal();

  double vrhs = 1.0/phasemanager.Density(curphase)*timefacfac;

  // linearization of reaction term w.r.t scalars
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    const double v = vrhs*funct(vi);

    for (int ui=0; ui<nen; ++ui)
    {
      const double vfunct = v*funct(ui);
      for (int iscal=0; iscal<numscal; ++iscal)
      {
        const int fui = ui*numscal+iscal;
        mymat(fvi,fui) += vfunct*phasemanager.ReacDerivScalar(curphase,iscal);
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
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // in case of an incompressible fluid phase 'curphase' (1/K = 0) the term cancels out
  if(phasemanager.IncompressibleFluidPhase(curphase))
    return;

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  // saturation
  const double saturation = phasemanager.Saturation(curphase);

  // inverse bulk modulus of phase (compressibility)
  const double invbulkmodulus = phasemanager.InvBulkmodulus(curphase);

  // pre factor
  const double facfacmass = fac*phasemanager.Porosity()*saturation*invbulkmodulus;
  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  for (int vi=0; vi<nen; ++vi)
  {
    const double v = facfacmass*funct(vi);
    const int fvi = vi*numdofpernode+phasetoadd;

    for (int ui=0; ui<nen; ++ui)
    {
      const double vfunct = v*funct(ui);
      for (int idof=0; idof<numdofpernode; ++idof)
      {
        const int fui = ui*numdofpernode+idof;

        mymat(fvi,fui) += vfunct*phasemanager.PressureDeriv(curphase,idof);
      }
    }
  }

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if(!inittimederiv)
  {
    //----------------------------------------------------------------
    // linearization of saturation w.r.t. dof
    //----------------------------------------------------------------
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      double hist=0.0;
      //if(curphase==phasetoadd) // bug fix??
        hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass2 = fac*phasemanager.PressureDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

      for (int idof=0; idof<numdofpernode; ++idof)
      {
        if(idof!=phasetoadd)
          facfacmass2 += timefacfac*phasemanager.PressureDeriv(curphase,idof)*(*variablemanager.Phidtnp())[idof];
      }

      facfacmass2 *= phasemanager.Porosity()*invbulkmodulus;

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = facfacmass2*funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;

            mymat(fvi,fui) +=
                vfunct*phasemanager.SaturationDeriv(curphase,idof);
          }
        }
      }
    }

    //----------------------------------------------------------------
    // linearization of porosity w.r.t. dof
    //----------------------------------------------------------------
    if(phasemanager.PorosityDependsOnFluid())
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      double hist=0.0;
      //if(curphase==phasetoadd)  //bugfix??
        hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass = fac*phasemanager.PressureDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

      for (int idof=0; idof<numdofpernode; ++idof)
      {
        if(idof!=phasetoadd)
          facfacmass += timefacfac*phasemanager.PressureDeriv(curphase,idof)*(*variablemanager.Phidtnp())[idof];
      }

      facfacmass *= saturation*invbulkmodulus;

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = facfacmass*funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;

            mymat(fvi,fui) +=
                vfunct*phasemanager.PorosityDeriv(idof);
          }
        }
      }
    }
  } // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{

  // in case of an incompressible fluid phase 'curphase' (1/K = 0) the term cancels out
  if(phasemanager.IncompressibleFluidPhase(curphase))
    return;

  // for the initial time derivative calculation no transient terms enter the rhs
  if(!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    double vtrans = GetRhsTrans(
        curphase,
        phasetoadd,
        numdofpernode,
        phasemanager,
        variablemanager,
        rhsfac,
        fac
      );

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;

      myvec[fvi] -= vtrans*funct(vi);
    }
  } // !inittimederiv
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{

  // in case of an incompressible fluid phase 'curphase' (1/K = 0) the term cancels out
  if(phasemanager.IncompressibleFluidPhase(curphase))
    return;

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  double vtrans = GetRhsTrans(
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      timefacfac,
      fac
    );

  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space, i.e. det(dx/ds)
  // in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac * N_x
  //              fac        = J              --> d(fac)/dd        = fac * N_x
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    const double v = funct(vi)*vtrans;

    for (int ui=0; ui<nen; ++ui)
    {
      for (int idim=0; idim<nsd; ++idim)
      {
        const int fui = ui*nsd+idim;

        mymat(fvi,fui) += v*derxy(idim,ui);
      }
    }
  }

  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  // in our case: vtrans is scaled with porosity in GetRhsTrans --> scale it with 1.0/porosity here

  if(phasemanager.PorosityDependsOnStruct())
  {
    vtrans *= 1.0/phasemanager.Porosity()*phasemanager.JacobianDefGrad()*phasemanager.PorosityDerivWrtJacobianDefGrad();

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;
      const double v = funct(vi)*vtrans;

      for (int ui=0; ui<nen; ++ui)
      {
        for (int idim=0; idim<nsd; ++idim)
        {
          const int fui = ui*nsd+idim;

          mymat(fvi,fui) += v*derxy(idim,ui);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate transient term at GP                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassPressure<nsd,nen>::GetRhsTrans(
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac
    )
{
  // read data from managers
  double hist = 0.0;
  //if(curphase==phasetoadd) //bugfix??
    hist = (*variablemanager.Hist())[phasetoadd];
  //std::cout << "hist = " << hist << std::endl;
  const double  porosity = phasemanager.Porosity();
  const std::vector<double>&  phinp = *variablemanager.Phinp();
  const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

  // saturation
  const double saturation = phasemanager.Saturation(curphase);

  // inverse bulk modulus of phase (compressibility)
  const double invbulkmodulus = phasemanager.InvBulkmodulus(curphase);

  double vtrans = 0.0;

  //TODO check for Genalpha
  // compute scalar at integration point
  vtrans = fac*phasemanager.PressureDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);
  for (int idof=0; idof<numdofpernode; ++idof)
    if(idof!=phasetoadd)
      vtrans += rhsfac*phasemanager.PressureDeriv(curphase,idof)*phidtnp[idof];

  vtrans *= porosity*saturation*invbulkmodulus;

  return vtrans;

}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // in case of an incompressible solid (1/K_s = 0) the term cancels out
  if(phasemanager.IncompressibleSolid())
    return;

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  //  get inverse bulkmodulus (=compressiblity)
  // TODO linearization of bulkmodulus
  const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

  //----------------------------------------------------------------
  // standard Galerkin transient term
  //----------------------------------------------------------------
  {
    const double facfacmass = fac*(1.0-phasemanager.Porosity())*invsolidbulkmodulus;
    for (int vi=0; vi<nen; ++vi)
    {
      const double v = facfacmass*funct(vi);
      const int fvi = vi*numdofpernode+phasetoadd;

      for (int ui=0; ui<nen; ++ui)
      {
        const double vfunct = v*funct(ui);
        for (int idof=0; idof<numdofpernode; ++idof)
        {
          const int fui = ui*numdofpernode+idof;

          mymat(fvi,fui) += vfunct*phasemanager.SolidPressureDeriv(idof);
        }
      }
    }
  }

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if(!inittimederiv)
  {
    //----------------------------------------------------------------
    // linearization of solid pressure derivative w.r.t. dof
    //----------------------------------------------------------------
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
      double hist=0.0;
      //if(curphase==phasetoadd) //bugfix??
      hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass3 = (1.0-phasemanager.Porosity())*invsolidbulkmodulus;

      std::vector<double> val(numdofpernode,0.0);

      for (int idof=0; idof<numdofpernode; ++idof)
      {
        if(idof==phasetoadd)
          for (int jdof=0; jdof<numdofpernode; ++jdof)
            val[jdof]+= fac*phasemanager.SolidPressureDerivDeriv(idof,jdof)*(phinp[idof]-hist);
        else
          for (int jdof=0; jdof<numdofpernode; ++jdof)
            val[jdof]+= timefacfac*phasemanager.SolidPressureDerivDeriv(idof,jdof)*(phidtnp[idof]);
      }

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = facfacmass3*funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;

            mymat(fvi,fui) +=vfunct*val[idof];
          }
        }
      }
    }

    //----------------------------------------------------------------
    // linearization of porosity w.r.t. dof
    //----------------------------------------------------------------
    if(phasemanager.PorosityDependsOnFluid())
    {
      const std::vector<double>& phinp = *variablemanager.Phinp();
      const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
      double hist=0.0;
      //if(curphase==phasetoadd) //bugfix??
        hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass3 = -1.0*invsolidbulkmodulus;

      std::vector<double> val(numdofpernode,0.0);

      for (int idof=0; idof<numdofpernode; ++idof)
      {
        const double solidpressurederiv = phasemanager.SolidPressureDeriv(idof);
        if(idof==phasetoadd)
          for (int jdof=0; jdof<numdofpernode; ++jdof)
            val[jdof]+= fac*solidpressurederiv*phasemanager.PorosityDeriv(jdof)*(phinp[idof]-hist);
        else
          for (int jdof=0; jdof<numdofpernode; ++jdof)
            val[jdof]+= timefacfac*solidpressurederiv*phasemanager.PorosityDeriv(jdof)*(phidtnp[idof]);
      }

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = facfacmass3*funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;

            mymat(fvi,fui) +=vfunct*val[idof];
          }
        }
      }
    }
  } // !inittimederiv
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // in case of an incompressible solid (1/K_s = 0) the term cancels out
  if(phasemanager.IncompressibleSolid())
    return;

  // for the initial time derivative calculation no transient terms enter the rhs
  if(!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    double vtrans = GetRhsTrans(
        curphase,
        phasetoadd,
        numdofpernode,
        phasemanager,
        variablemanager,
        rhsfac,
        fac
      );

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;

      myvec[fvi] -= vtrans*funct(vi);
    }
  }// !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{

  // in case of an incompressible solid (1/K_s = 0) the term cancels out
  if(phasemanager.IncompressibleSolid())
    return;

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  double vtrans = GetRhsTrans(
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      timefacfac,
      fac
    );

  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space, i.e. det(dx/ds)
  // in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac * N_x
  //              fac        = J              --> d(fac)/dd        = fac * N_x
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    const double v = funct(vi)*vtrans;

    for (int ui=0; ui<nen; ++ui)
    {
      for (int idim=0; idim<nsd; ++idim)
      {
        const int fui = ui*nsd+idim;

        mymat(fvi,fui) += v*derxy(idim,ui);
      }
    }
  }

  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  // in our case: vtrans is scaled with (1.0-porosity) in GetRhsTrans --> scale it with 1.0/(1.0-porosity) here

  if(phasemanager.PorosityDependsOnStruct())
  {
    vtrans *= -1.0/(1.0-phasemanager.Porosity())*phasemanager.JacobianDefGrad()*phasemanager.PorosityDerivWrtJacobianDefGrad();

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;
      const double v = funct(vi)*vtrans;

      for (int ui=0; ui<nen; ++ui)
      {
        for (int idim=0; idim<nsd; ++idim)
        {
          const int fui = ui*nsd+idim;

          mymat(fvi,fui) += v*derxy(idim,ui);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate transient term at GP                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressure<nsd,nen>::GetRhsTrans(
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac
    )
{
  // read data from managers
  double hist = 0.0;
  //if(curphase==phasetoadd) //bugfix??
    hist = (*variablemanager.Hist())[phasetoadd];
  const double  porosity = phasemanager.Porosity();
  const std::vector<double>&  phinp = *variablemanager.Phinp();
  const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

  //  get inverse bulkmodulus (=compressiblity)
  const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

  //TODO check genalpha
  // compute scalar at integration point
  double vtrans = fac*phasemanager.SolidPressureDeriv(phasetoadd)*(phinp[phasetoadd]-hist);

  for (int idof=0; idof<numdofpernode; ++idof)
  {
    if(idof!=phasetoadd)
    {
      vtrans += rhsfac*phasemanager.SolidPressureDeriv(idof)*phidtnp[idof];
    }
  }

  vtrans *= (1.0-porosity)*invsolidbulkmodulus;

  return vtrans;

}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // call base class with scaled factors
  EvaluatorMassSolidPressure<nsd,nen>::EvaluateMatrixAndAssemble(
      elemat,
      funct,
      derxy,
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      phasemanager.Saturation(curphase)*timefacfac,
      phasemanager.Saturation(curphase)*fac,
      inittimederiv
    );

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if(!inittimederiv)
  {
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
      double hist=0.0;
      //if(curphase==phasetoadd) //bugfix??
        hist = (*variablemanager.Hist())[phasetoadd];

      double facfacmass = fac*phasemanager.SolidPressureDeriv(phasetoadd)*(phinp[phasetoadd]-hist);

      for (int idof=0; idof<numdofpernode; ++idof)
      {
        if(idof!=phasetoadd)
          facfacmass += timefacfac*phasemanager.SolidPressureDeriv(idof)*phidtnp[idof];
      }

      facfacmass *= (1.0-phasemanager.Porosity())*invsolidbulkmodulus;

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = facfacmass*funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;

            mymat(fvi,fui) +=
                vfunct*phasemanager.SaturationDeriv(curphase,idof);
          }
        }
      }
    }
  }// !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // call base class with scaled factors
  EvaluatorMassSolidPressure<nsd,nen>::EvaluateVectorAndAssemble(
      elevec,
      funct,
      derxy,
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      phasemanager.Saturation(curphase)*rhsfac,
      phasemanager.Saturation(curphase)*fac,
      inittimederiv);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  // call base class with scaled factors
  EvaluatorMassSolidPressure<nsd,nen>::EvaluateMatrixODStructAndAssemble(
      elemat,
      funct,
      deriv,
      derxy,
      xjm,
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      phasemanager.Saturation(curphase)*timefacfac,
      phasemanager.Saturation(curphase)*fac,
      det);

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSolidPressureSat<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  //----------------------------------------------------------------
  // linearization of saturation w.r.t. dof
  //----------------------------------------------------------------
  {
    const double facfacmass = fac*phasemanager.Porosity();
    for (int vi=0; vi<nen; ++vi)
    {
      const double v = facfacmass*funct(vi);
      const int fvi = vi*numdofpernode+phasetoadd;

      for (int ui=0; ui<nen; ++ui)
      {
        const double vfunct = v*funct(ui);
        for (int idof=0; idof<numdofpernode; ++idof)
        {
          const int fui = ui*numdofpernode+idof;

          mymat(fvi,fui) += vfunct*phasemanager.SaturationDeriv(curphase,idof);
        }
      }
    }
  }

  // for the initial time derivative calculation we only need the standard Galerkin transient term
  // and no additional linearizations
  if(!inittimederiv)
  {
    //----------------------------------------------------------------
    // linearization of porosity w.r.t. dof
    //----------------------------------------------------------------
    if(phasemanager.PorosityDependsOnFluid())
    {
      // read data from manager
      double hist = 0.0;
      //if(curphase==phasetoadd) //bugfix??
        hist = (*variablemanager.Hist())[phasetoadd];
      const std::vector<double>&  phinp = *variablemanager.Phinp();
      const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

      //TODO genalpha
      // compute scalar at integration point
      double facfacmass = fac*phasemanager.SaturationDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

      for (int idof=0; idof<numdofpernode; ++idof)
        if(phasetoadd!=idof)
          facfacmass += timefacfac*phasemanager.SaturationDeriv(curphase,idof)*phidtnp[idof];

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = facfacmass*funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui);
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            const int fui = ui*numdofpernode+idof;

            mymat(fvi,fui) += vfunct*phasemanager.PorosityDeriv(idof);
          }
        }
      }
    }

    //----------------------------------------------------------------
    // linearization of derivative of saturation w.r.t. dof
    //----------------------------------------------------------------
    {

      // read data from manager
      double hist = 0.0;
      hist = (*variablemanager.Hist())[phasetoadd];
      const std::vector<double>&  phinp = *variablemanager.Phinp();
      const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

      /*for (int iphase=0; iphase < numdofpernode; iphase++)
      {
        std::cout << iphase << " =====================================" << std::endl;
        for (int jphase = 0; jphase < numdofpernode; jphase++)
        {
          for (int kphase = 0; kphase < numdofpernode; kphase++)
          {
            std::cout << std::setprecision(8) << phasemanager.SaturationDerivDeriv(iphase,jphase,kphase) << "  ";
          }
          std::cout << "\n";
        }
      }*/

      for (int vi=0; vi<nen; ++vi)
      {
        const double v = funct(vi);
        const int fvi = vi*numdofpernode+phasetoadd;

        for (int ui=0; ui<nen; ++ui)
        {
          const double vfunct = v*funct(ui)*phasemanager.Porosity();
          for (int idof=0; idof<numdofpernode; ++idof)
          {
            if(idof == phasetoadd)
            {
              for (int jdof=0; jdof<numdofpernode; ++jdof)
              {
                const int fui = ui*numdofpernode+jdof;
                mymat(fvi,fui) += fac*vfunct*phasemanager.SaturationDerivDeriv(curphase,phasetoadd,jdof)*(phinp[phasetoadd]-hist);
              }
            }
            else
            {
              for (int jdof=0; jdof<numdofpernode; ++jdof)
              {
                const int fui = ui*numdofpernode+jdof;
                mymat(fvi,fui) += timefacfac*vfunct*phasemanager.SaturationDerivDeriv(curphase,idof,jdof)*(phidtnp[idof]);
              }
            }
          }
        }
      }
    }
  } // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // for the initial time derivative calculation no transient terms enter the rhs
  if(!inittimederiv)
  {
    // get vector to fill
    Epetra_SerialDenseVector& myvec = *elevec[0];

    double vtrans = GetRhsTrans(
        curphase,
        phasetoadd,
        numdofpernode,
        phasemanager,
        variablemanager,
        rhsfac,
        fac
      );

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;

      myvec[fvi] -= vtrans*funct(vi);
    }
  } // !inittimederiv

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{

  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  double vtrans = GetRhsTrans(
      curphase,
      phasetoadd,
      numdofpernode,
      phasemanager,
      variablemanager,
      timefacfac,
      fac
    );

  // linearization of mesh motion
  //------------------------------------------------dJ/dd = dJ/dF : dF/dd = J * F^-T . N_{\psi} = J * N_x
  // J denotes the determinant of the Jacobian of the mapping between current and parameter space, i.e. det(dx/ds)
  // in our case: timefacfac = J * dt * theta --> d(timefacfac)/dd = timefacfac * N_x
  //              fac        = J              --> d(fac)/dd        = fac * N_x
  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;
    const double v = funct(vi)*vtrans;

    for (int ui=0; ui<nen; ++ui)
    {
      for (int idim=0; idim<nsd; ++idim)
      {
        const int fui = ui*nsd+idim;

        mymat(fvi,fui) += v*derxy(idim,ui);
      }
    }
  }

  // linearization of porosity
  //------------------------------------------------dporosity/dd = dporosity/dJ * dJ/dd = dporosity/dJ * J * N_x
  // J denotes the determinant of the deformation gradient, i.e. det F = det ( d x / d X ) = det (dx/ds) * ( det(dX/ds) )^-1
  // in our case: vtrans is scaled with porosity in GetRhsTrans --> scale it with 1.0/porosity here

  if(phasemanager.PorosityDependsOnStruct())
  {
    vtrans *= 1.0/phasemanager.Porosity()*phasemanager.JacobianDefGrad()*phasemanager.PorosityDerivWrtJacobianDefGrad();

    for (int vi=0; vi<nen; ++vi)
    {
      const int fvi = vi*numdofpernode+phasetoadd;
      const double v = funct(vi)*vtrans;

      for (int ui=0; ui<nen; ++ui)
      {
        for (int idim=0; idim<nsd; ++idim)
        {
          const int fui = ui*nsd+idim;

          mymat(fvi,fui) += v*derxy(idim,ui);
        }
      }
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate transient term at GP                       kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
double DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorMassSaturation<nsd,nen>::GetRhsTrans(
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac
    )
{
  // read data from manager
  double hist = 0.0;
  //if(curphase==phasetoadd) //bugfix??
    hist = (*variablemanager.Hist())[phasetoadd];

  const double  porosity = phasemanager.Porosity();
  const std::vector<double>&  phinp = *variablemanager.Phinp();
  const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

  //TODO genalpha
  // compute scalar at integration point
  double vtrans = fac*phasemanager.SaturationDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

  for (int idof=0; idof<numdofpernode; ++idof)
    if(phasetoadd!=idof)
      vtrans += rhsfac*phasemanager.SaturationDeriv(curphase,idof)*phidtnp[idof];
  //double vtrans = fac*(phasemanager.Saturation(phasetoadd)-phasemanager.SaturationHist(phasetoadd,histvec));

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
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get vectors to be filled
  // pressure
  Epetra_SerialDenseVector&     pressure = *elevec[0];
  // saturation
  Epetra_SerialDenseVector&     saturation = *elevec[1];
  //counter
  Epetra_SerialDenseVector&     counter = *elevec[2];

  for(int inode=0;inode<nen;inode++)
  {
    // save the pressure value
    pressure[inode*numdofpernode+curphase] += fac*funct(inode)*phasemanager.Pressure(curphase);
    // save the saturation value
    saturation[inode*numdofpernode+curphase] += fac*funct(inode)*phasemanager.Saturation(curphase);
    // mark the evaluated node
    counter[inode*numdofpernode+curphase] += fac*funct(inode);
  }
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPressureAndSaturation<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get vectors to be filled
  Epetra_SerialDenseVector&     solidpressure = *elevec[0];
  Epetra_SerialDenseVector&     counter = *elevec[1];

  for(int inode=0;inode<nen;inode++)
  {
    // save the pressure value
    solidpressure[inode] += fac*funct(inode)*phasemanager.SolidPressure();
    // mark the evaluated node
    counter[inode] += fac*funct(inode);
  }
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorSolidPressure<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                             kremheller 04/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // do nothing
  return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                 kremheller 04/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  // get vectors to be filled
  Epetra_SerialDenseVector&     porosity = *elevec[0];
  Epetra_SerialDenseVector&     counter = *elevec[1];

  for(int inode=0;inode<nen;inode++)
  {
    // save the porosity value
    porosity[inode] += fac*funct(inode)*phasemanager.Porosity();
    // mark the evaluated node
    counter[inode] += fac*funct(inode);
  }
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 04/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorPorosity<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // get matrixes to fill
  Epetra_SerialDenseMatrix&     linearization = *elemat[0];
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
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  //nothing to do
   return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxLinearization<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 * **********************************************************************
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate element matrix                                   vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd,nen>::EvaluateMatrixAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    bool                                                        inittimederiv
  )
{
  // get matrixes to fill
   Epetra_SerialDenseMatrix&     rhs = *elemat[1];

   const std::vector<LINALG::Matrix<nsd,1> >& gradphi = *variablemanager.GradPhinp();

   // current pressure gradient
   LINALG::Matrix<nsd,1> gradpres(true);
   gradpres.Clear();

   // compute the pressure gradient from the phi gradients
   for (int idof=0; idof<numdofpernode; ++idof)
     gradpres.Update(phasemanager.PressureDeriv(curphase,idof),gradphi[idof],1.0);

   double abspressgrad=0.0;
   for (int i=0; i<nsd; i++)
     abspressgrad+=gradpres(i)*gradpres(i);
   abspressgrad=sqrt(abspressgrad);

   // diffusion tensor
   LINALG::Matrix<nsd,nsd> difftensor(true);
   phasemanager.PermeabilityTensor(curphase,difftensor);
   difftensor.Scale(phasemanager.RelPermeability(curphase)/phasemanager.DynViscosity(curphase, abspressgrad));

   // diffusive flux
   static LINALG::Matrix<nsd,1> diffflux(true);
   diffflux.Multiply(-1.0,difftensor,gradpres);

   // Compute element vectors. For L2-Projection
   for (int node_i=0;node_i<nen;node_i++)
   {
     for (int j=0;j<nsd;j++)
     {
       rhs(node_i,nsd*curphase+j) += funct(node_i) * fac * diffflux(j);
     }
   }
   return;
}

/*----------------------------------------------------------------------*
 | evaluate RHS vector                                      vuong 09/16 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd,nen>::EvaluateVectorAndAssemble(
    std::vector<Epetra_SerialDenseVector*>&                     elevec,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      rhsfac,
    double                                                      fac,
    bool                                                        inittimederiv
    )
{
  //nothing to do
   return;
};

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with structure kremheller 03/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd,nen>::EvaluateMatrixODStructAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              deriv,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    const LINALG::Matrix<nsd,nsd>&                              xjm,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac,
    double                                                      det
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 | evaluate off-diagonal coupling matrix with scatra   kremheller 06/17 |
 *----------------------------------------------------------------------*/
template <int nsd, int nen>
void DRT::ELEMENTS::POROFLUIDEVALUATOR::ReconstructFluxRHS<nsd,nen>::EvaluateMatrixODScatraAndAssemble(
    std::vector<Epetra_SerialDenseMatrix*>&                     elemat,
    const LINALG::Matrix<nen,1>&                                funct,
    const LINALG::Matrix<nsd,nen>&                              derxy,
    int                                                         curphase,
    int                                                         phasetoadd,
    int                                                         numdofpernode,
    const POROFLUIDMANAGER::PhaseManagerInterface&              phasemanager,
    const POROFLUIDMANAGER::VariableManagerInterface<nsd,nen>&  variablemanager,
    double                                                      timefacfac,
    double                                                      fac
  )
{
  //nothing to do
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes

// 1D elements
// line 2
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<1,2>;

// line 3
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<1,3>;

// 2D elements
// tri3
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<2,3>;
// quad4
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<2,4>;

// quad9 and nurbs9
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<2,9>;

// 3D elements
// hex8
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3,8>;

// hex27
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3,27>;
// tet4
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3,4>;
// tet10
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3,10>;
// pyramid5
template class DRT::ELEMENTS::POROFLUIDEVALUATOR::EvaluatorInterface<3,5>;
