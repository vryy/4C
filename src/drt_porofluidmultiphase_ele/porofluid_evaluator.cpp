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

  // determine action
  switch (action)
  {
  // calculate true pressures and saturation
  case POROFLUIDMULTIPHASE::calc_mat_and_rhs:
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
//      assembler = Teuchos::rcp(new AssembleStandard(iphase));
//      tmpevaluator = Teuchos::rcp(new EvaluatorConv<nsd, nen>(assembler,iphase));
//      evaluator_phase->AddEvaluator(tmpevaluator);

      // add evaluator for the convective conservative term (w, S \nabla \cdot v )
      if (para.IsAle())
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
        tmpevaluator = Teuchos::rcp(new EvaluatorSatDivVel<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);
      }

      // add evaluator for the diffusive term (\nabla w, K \nabla p)
      // the diffusive term is also assembled into the last phase
      assembler = Teuchos::rcp(new AssembleAlsoIntoOtherPhase(curphase, numdofpernode-1));
      tmpevaluator = Teuchos::rcp(new EvaluatorDiff<nsd, nen>(assembler,curphase));
      evaluator_phase->AddEvaluator(tmpevaluator);

      // add evaluator for the reactive term
      if (phasemanager.IsReactive(curphase))
      {
        // the reactive term is also assembled into the last phase
        assembler = Teuchos::rcp(
            new AssembleAlsoIntoOtherPhase(curphase, numdofpernode-1));
        tmpevaluator = Teuchos::rcp(new EvaluatorReac<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);
      }

      // add evaluators for the instationary terms
      if (not para.IsStationary())
      {
        // add evaluator for the instationary pressure term
        // the term is also assembled into the last phase
        assembler = Teuchos::rcp(
            new AssembleAlsoIntoOtherPhase(curphase, numdofpernode-1));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassPressure<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);

        // add evaluator for the instationary solid pressure term
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassSolidPressureSat<nsd, nen>(assembler,curphase));
        evaluator_phase->AddEvaluator(tmpevaluator);

        // add evaluator for the instationary saturation term
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
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
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
        tmpevaluator = Teuchos::rcp(new EvaluatorDivVel<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);
      }

      // add evaluator for the diffusive term (\nabla w, K \nabla p)
      assembler = Teuchos::rcp(new AssembleStandard(curphase));
      tmpevaluator = Teuchos::rcp(new EvaluatorDiff<nsd, nen>(assembler,curphase));
      evaluator_lastphase->AddEvaluator(tmpevaluator);

      // add evaluator for the reactive term
      if (phasemanager.IsReactive(curphase))
      {
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
        tmpevaluator = Teuchos::rcp(new EvaluatorReac<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);
      }

      // add evaluators for the instationary terms
      if (not para.IsStationary())
      {
        // add evaluator for the instationary pressure term
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
        tmpevaluator = Teuchos::rcp(new EvaluatorMassPressure<nsd, nen>(assembler,curphase));
        evaluator_lastphase->AddEvaluator(tmpevaluator);

        // add evaluator for the instationary solid pressure term
        assembler = Teuchos::rcp(new AssembleStandard(curphase));
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
      assembler = Teuchos::rcp(new AssembleStandard(iphase));
      tmpevaluator = Teuchos::rcp(new EvaluatorPressureAndSaturation<nsd, nen>(assembler,iphase));
      evaluator_multiphase->AddEvaluator(tmpevaluator);
    }
    evaluator = evaluator_multiphase;

    break;
  }
  case POROFLUIDMULTIPHASE::calc_solidpressure:
  {
    Teuchos::RCP<AssembleInterface> assembler = Teuchos::rcp(new AssembleStandard(-1));
    evaluator = Teuchos::rcp(new EvaluatorSolidPressure<nsd, nen>(assembler,-1));

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

    assembler = Teuchos::rcp(new AssembleStandard(-1));
    tmpevaluator = Teuchos::rcp(new ReconstructFluxLinearization<nsd, nen>(assembler,-1));
    evaluator_multiphase->AddEvaluator(tmpevaluator);

    // build evaluators for all but last phase
    for (int iphase = 0; iphase < numdofpernode; iphase++)
    {
      assembler = Teuchos::rcp(new AssembleStandard(iphase));
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
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
      fac
    );

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
    double                                                      fac
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
      phasemanager.Saturation(curphase)*fac);

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
    double                                                      fac
  )
{
  // get matrix to fill
  Epetra_SerialDenseMatrix& mymat = *elemat[0];

  // diffusion tensor
  LINALG::Matrix<nsd,nsd> difftensor(true);
  phasemanager.DiffTensor(curphase,difftensor);

  static LINALG::Matrix<nsd,nen> diffflux(true);
  diffflux.Multiply(difftensor,derxy);
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
      }
    }
  }
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
    double                                                      fac
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  const std::vector<LINALG::Matrix<nsd,1> >& gradphi = *variablemanager.GradPhinp();

  // diffusion tensor
  LINALG::Matrix<nsd,nsd> difftensor(true);
  phasemanager.DiffTensor(curphase,difftensor);

  // current pressure gradient
  LINALG::Matrix<nsd,1> gradpres(true);
  gradpres.Clear();

  // compute the pressure gradient from the phi gradients
  for (int idof=0; idof<numdofpernode; ++idof)
    gradpres.Update(phasemanager.PressureDeriv(curphase,idof),gradphi[idof],1.0);

  // diffusive flux
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
    double                                                      fac
  )
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
    double                                                      fac
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
    double                                                      fac
  )
{
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

  //----------------------------------------------------------------
  // linearization of saturation w.r.t. dof
  //----------------------------------------------------------------
  {
    const std::vector<double>& phinp = *variablemanager.Phinp();
    double hist=0.0;
    if(curphase==phasetoadd)
      hist = (*variablemanager.Hist())[curphase];

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
  {
    const std::vector<double>& phinp = *variablemanager.Phinp();
    double hist=0.0;
    if(curphase==phasetoadd)
      hist = (*variablemanager.Hist())[curphase];

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
    double                                                      fac
    )
{

  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  // read data from managers
  double hist = 0.0;
  if(curphase==phasetoadd)
    hist = (*variablemanager.Hist())[curphase];
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

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    myvec[fvi] -= vtrans*funct(vi);
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
    double                                                      fac
  )
{
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

  //----------------------------------------------------------------
  // linearization of solid pressure derivative w.r.t. dof
  //----------------------------------------------------------------
  {
    const std::vector<double>& phinp = *variablemanager.Phinp();
    const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
    double hist=0.0;
    if(curphase==phasetoadd)
      hist = (*variablemanager.Hist())[curphase];

    double facfacmass3 = (1.0-phasemanager.Porosity())*invsolidbulkmodulus;

    std::vector<double> val(numdofpernode,0.0);

    for (int idof=0; idof<numdofpernode; ++idof)
    {
      if(idof==curphase)
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
  {
    const std::vector<double>& phinp = *variablemanager.Phinp();
    const std::vector<double>& phidtnp = *variablemanager.Phidtnp();
    double hist=0.0;
    if(curphase==phasetoadd)
      hist = (*variablemanager.Hist())[curphase];

    double facfacmass3 = -1.0*invsolidbulkmodulus;

    std::vector<double> val(numdofpernode,0.0);

    for (int idof=0; idof<numdofpernode; ++idof)
    {
      const double solidpressurederiv = phasemanager.SolidPressureDeriv(idof);
      if(idof==curphase)
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
    double                                                      fac
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  // read data from managers
  double hist = 0.0;
  if(curphase==phasetoadd)
    hist = (*variablemanager.Hist())[curphase];
  const double  porosity = phasemanager.Porosity();
  const std::vector<double>&  phinp = *variablemanager.Phinp();
  const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

  //  get inverse bulkmodulus (=compressiblity)
  const double invsolidbulkmodulus = phasemanager.InvBulkmodulusSolid();

  //TODO check genalpha
  // compute scalar at integration point
  double vtrans = fac*phasemanager.SolidPressureDeriv(curphase)*(phinp[curphase]-hist);

  for (int idof=0; idof<numdofpernode; ++idof)
  {
    if(idof!=curphase)
    {
      vtrans += rhsfac*phasemanager.SolidPressureDeriv(idof)*phidtnp[idof];
    }
  }

  vtrans *= (1.0-porosity)*invsolidbulkmodulus;

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    myvec[fvi] -= vtrans*funct(vi);
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
    double                                                      fac
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
      phasemanager.Saturation(curphase)*fac
    );

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
    if(curphase==phasetoadd)
      hist = (*variablemanager.Hist())[curphase];

    double facfacmass = fac*phasemanager.SolidPressureDeriv(curphase)*(phinp[curphase]-hist);

    for (int idof=0; idof<numdofpernode; ++idof)
    {
      if(idof!=curphase)
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
    double                                                      fac
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
      phasemanager.Saturation(curphase)*fac);

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
    double                                                      fac
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

  //----------------------------------------------------------------
  // linearization of porosity w.r.t. dof
  //----------------------------------------------------------------
  {
    // read data from manager
    double hist = 0.0;
    if(curphase==phasetoadd)
      hist = (*variablemanager.Hist())[curphase];
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
    double                                                      fac
    )
{
  // get matrix to fill
  Epetra_SerialDenseVector& myvec = *elevec[0];

  // read data from manager
  double hist = 0.0;
  if(curphase==phasetoadd)
    hist = (*variablemanager.Hist())[curphase];
  const double  porosity = phasemanager.Porosity();
  const std::vector<double>&  phinp = *variablemanager.Phinp();
  const std::vector<double>&  phidtnp = *variablemanager.Phidtnp();

  //TODO genalpha
  // compute scalar at integration point
  double vtrans = fac*phasemanager.SaturationDeriv(curphase,phasetoadd)*(phinp[phasetoadd]-hist);

  for (int idof=0; idof<numdofpernode; ++idof)
    if(phasetoadd!=idof)
      vtrans += rhsfac*phasemanager.SaturationDeriv(curphase,idof)*phidtnp[idof];

  vtrans *= porosity;

  for (int vi=0; vi<nen; ++vi)
  {
    const int fvi = vi*numdofpernode+phasetoadd;

    myvec[fvi] -= vtrans*funct(vi);
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
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
    double                                                      fac
    )
{
  //nothing to do
   return;
};

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
    double                                                      fac
  )
{
  // get matrixes to fill
   Epetra_SerialDenseMatrix&     rhs = *elemat[1];

   // diffusion tensor
   LINALG::Matrix<nsd,nsd> difftensor(true);
   phasemanager.DiffTensor(curphase,difftensor);

   // current pressure gradient
   LINALG::Matrix<nsd,1> gradpres(true);
   gradpres.Clear();

   // compute the pressure gradient from the phi gradients
   for (int idof=0; idof<numdofpernode; ++idof)
   {
     gradpres.Update(phasemanager.PressureDeriv(curphase,idof),(*variablemanager.GradPhinp())[idof],1.0);
   }

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
    double                                                      fac
    )
{
  //nothing to do
   return;
};


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
