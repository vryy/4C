/*----------------------------------------------------------------------*/
/*!
\file scatra_timint_convcheck_strategies.cpp

\brief strategies for Newton-Raphson convergence check for scalar transport problems

To keep the scalar transport time integrator class and derived classes as plain as possible,
the convergence check for the Newton-Raphson iteration has been encapsulated within separate
strategy classes. Every specific convergence check strategy (e.g., for standard scalar transport
or electrochemistry problems with or without scatra-scatra interface coupling involving Lagrange
multipliers) computes, checks, and outputs different relevant vector norms and is implemented
in a subclass derived from an abstract, purely virtual interface class.

\level 1

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
 */
/*----------------------------------------------------------------------*/
#include "scatra_timint_convcheck_strategies.H"
#include "scatra_timint_elch.H"
#include "scatra_timint_implicit.H"
#include "scatra_timint_meshtying_strategy_s2i.H"

#include "../drt_lib/drt_discret.H"

#include "../linalg/linalg_mapextractor.H"

#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*
 | constructor                                               fang 02/16 |
 *----------------------------------------------------------------------*/
SCATRA::ConvCheckStrategyBase::ConvCheckStrategyBase(
    const Teuchos::ParameterList&   parameters   //!< parameter list for Newton-Raphson iteration
    ) :
itmax_(parameters.get<int>("ITEMAX")),
ittol_(parameters.get<double>("CONVTOL")),
abstolres_(parameters.get<double>("ABSTOLRES")),
itmax_outer_(parameters.get<int>("ITEMAX_OUTER")),
ittol_outer_(parameters.get<double>("CONVTOL_OUTER"))
{
  return;
}


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 02/16 |
 *----------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyStd::AbortNonlinIter(
    const ScaTraTimIntImpl&   scatratimint,   //!< scalar transport time integrator
    double&                   actresidual     //!< return maximum current residual value
    ) const
{
  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current Newton-Raphson iteration step
  const int itnum = scatratimint.IterNum();

  // compute L2 norm of concentration state vector
  double conc_state_L2(0.0);
  scatratimint.Phinp()->Norm2(&conc_state_L2);

  // compute L2 norm of concentration residual vector
  double conc_res_L2(0.);
  scatratimint.Residual()->Norm2(&conc_res_L2);

  // compute infinity norm of concentration residual vector
  double conc_res_inf(0.);
  scatratimint.Residual()->NormInf(&conc_res_inf);

  // compute L2 norm of concentration increment vector
  double conc_inc_L2(0.);
  scatratimint.Increment()->Norm2(&conc_inc_L2);

  // safety checks
  if(std::isnan(conc_state_L2) or std::isnan(conc_res_L2) or std::isnan(conc_inc_L2))
    dserror("Calculated vector norm for concentration is not a number!");
  if(std::isinf(conc_state_L2) or std::isinf(conc_res_L2) or std::isinf(conc_inc_L2))
    dserror("Calculated vector norm for concentration is infinity!");

  // care for the case that nothing really happens in the concentration field
  if(conc_state_L2 < 1.e-5)
    conc_state_L2 = 1.;

  // special case: very first iteration step --> solution increment is not yet available
  if(itnum == 1)
  {
    if(mypid == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- con-res ---|-- con-inc ---|-- con-res-inf ---|" << std::endl;

      // print first line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (      --     ,te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;
    }
  }

  // ordinary case: later iteration steps --> solution increment can be printed and convergence check should be done
  else
  {
    if(mypid == 0)
      // print current line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_inc_L2/conc_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (ts="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtSolve() << ",te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;

    // convergence check
    if(conc_res_L2 <= ittol_ and conc_inc_L2/conc_state_L2 <= ittol_)
    {
      if(mypid == 0)
        // print finish line of convergence table to screen
        std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;

      return true;
    }
  }

  // abort iteration when there is nothing more to do --> better robustness
  // absolute tolerance determines whether residual is already zero
  // prevents additional solver calls that will not improve the solution anymore
  if(conc_res_L2 < abstolres_)
  {
    if(mypid == 0)
      // print finish line of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;

    return true;
  }

  // output warning in case maximum number of iteration steps is reached without convergence, and proceed to next time step
  if(itnum == itmax_)
  {
    if (mypid == 0)
    {
      std::cout << "+---------------------------------------------------------------+" << std::endl;
      std::cout << "|       >>>>>> Newton-Raphson iteration did not converge!       |" << std::endl;
      std::cout << "+---------------------------------------------------------------+" << std::endl;
    }

    return true;
  }

  // return maximum residual value for adaptivity of linear solver tolerance
  actresidual = std::max(conc_res_L2,conc_inc_L2/conc_state_L2);

  // proceed with next iteration step
  return false;
} // SCATRA::ConvCheckStrategyStd::AbortNonlinIter()


/*--------------------------------------------------------------------------------------------*
 | perform convergence check for outer iteration in partitioned coupling schemes   fang 08/17 |
 *--------------------------------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyStd::AbortOuterIter(
    const ScaTraTimIntImpl&   scatratimint   //!< scalar transport time integrator
    ) const
{
  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current outer iteration step
  const int itnum = scatratimint.IterNumOuter();

  // compute vector norms
  double L2_phinp(0.);
  double L2_phinp_inc(0.);
  scatratimint.Phinp()->Norm2(&L2_phinp);
  scatratimint.PhinpInc()->Norm2(&L2_phinp_inc);
  if(L2_phinp < 1.e-10)
    L2_phinp = 1.;

  // print convergence status
  if(mypid == 0)
  {
    std::cout << "|                                  OUTER ITERATION                                |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;
    std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_outer_ << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << ittol_outer_ << "[L_2 ]  |      --      | "
              << std::setw(10) << std::setprecision(3) << std::scientific << L2_phinp_inc/L2_phinp << "   |        --        |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;
  }

  // convergence check
  if(L2_phinp_inc/L2_phinp <= ittol_outer_)
  {
    if(mypid == 0)
    {
      std::cout << "|            OUTER ITERATION LOOP CONVERGED AFTER ITERATION"
                << std::setw(4) << itnum << "/"
                << std::setw(4) << itmax_outer_ << " !            |" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;
    }

    return true;
  }

  // throw error in case maximum number of iteration steps is reached without convergence
  else if(itnum == itmax_outer_)
  {
    if(mypid == 0)
    {
      std::cout << "|         >>>>>> not converged within maximum number of iteration steps!          |" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+------------------+" << std::endl;
    }

    dserror("Outer iteration did not converge within maximum number of iteration steps!");

    return true;
  }

  // proceed with next outer iteration step
  return false;
} // SCATRA::ConvCheckStrategyStd::AbortOuterIter


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 04/16 |
 *----------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyStdMicroScale::AbortNonlinIter(
    const ScaTraTimIntImpl&   scatratimint,   //!< scalar transport time integrator
    double&                   actresidual     //!< return maximum current residual value
    ) const
{
  // extract current Newton-Raphson iteration step
  const int itnum = scatratimint.IterNum();

  // compute L2 norm of concentration state vector
  double conc_state_L2(0.0);
  scatratimint.Phinp()->Norm2(&conc_state_L2);

  // compute L2 norm of concentration residual vector
  double conc_res_L2(0.);
  scatratimint.Residual()->Norm2(&conc_res_L2);

  // compute infinity norm of concentration residual vector
  double conc_res_inf(0.);
  scatratimint.Residual()->NormInf(&conc_res_inf);

  // compute L2 norm of concentration increment vector
  double conc_inc_L2(0.);
  scatratimint.Increment()->Norm2(&conc_inc_L2);

  // safety checks
  if(std::isnan(conc_state_L2) or std::isnan(conc_res_L2) or std::isnan(conc_inc_L2))
    dserror("Calculated vector norm for concentration is not a number!");
  if(std::isinf(conc_state_L2) or std::isinf(conc_res_L2) or std::isinf(conc_inc_L2))
    dserror("Calculated vector norm for concentration is infinity!");

  // care for the case that nothing really happens in the concentration field
  if(conc_state_L2 < 1.e-5)
    conc_state_L2 = 1.;

  // convergence check
  if((itnum > 1 and conc_res_L2 <= ittol_ and conc_inc_L2/conc_state_L2 <= ittol_) or conc_res_L2 < abstolres_)
    return true; // converged

  // throw error in case maximum number of iteration steps is reached without convergence
  if(itnum == itmax_)
    dserror("Newton-Raphson algorithm on micro scale did not converge!");

  // return maximum residual value for adaptivity of linear solver tolerance
  actresidual = std::max(conc_res_L2,conc_inc_L2/conc_state_L2);

  // proceed with next iteration step
  return false;
} // SCATRA::ConvCheckStrategyStdMicroScale::AbortNonlinIter()


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 02/16 |
 *----------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyStdElch::AbortNonlinIter(
    const ScaTraTimIntImpl&   scatratimint,   //!< scalar transport time integrator
    double&                   actresidual     //!< return maximum current residual value
    ) const
{
  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current Newton-Raphson iteration step
  const int itnum = scatratimint.IterNum();

  // compute L2 norm of concentration state vector
  Teuchos::RCP<Epetra_Vector> conc_vector = scatratimint.Splitter()->ExtractOtherVector(scatratimint.Phinp());
  double conc_state_L2(0.0);
  conc_vector->Norm2(&conc_state_L2);

  // compute L2 norm of concentration residual vector
  scatratimint.Splitter()->ExtractOtherVector(scatratimint.Residual(),conc_vector);
  double conc_res_L2(0.);
  conc_vector->Norm2(&conc_res_L2);

  // compute infinity norm of concentration residual vector
  double conc_res_inf(0.);
  conc_vector->NormInf(&conc_res_inf);

  // compute L2 norm of concentration increment vector
  scatratimint.Splitter()->ExtractOtherVector(scatratimint.Increment(),conc_vector);
  double conc_inc_L2(0.);
  conc_vector->Norm2(&conc_inc_L2);

  // compute L2 norm of electric potential state vector
  Teuchos::RCP<Epetra_Vector> pot_vector = scatratimint.Splitter()->ExtractCondVector(scatratimint.Phinp());
  double pot_state_L2(0.0);
  pot_vector->Norm2(&pot_state_L2);

  // compute L2 norm of electric potential residual vector
  scatratimint.Splitter()->ExtractCondVector(scatratimint.Residual(),pot_vector);
  double pot_res_L2(0.);
  pot_vector->Norm2(&pot_res_L2);

  // compute L2 norm of electric potential increment vector
  scatratimint.Splitter()->ExtractCondVector(scatratimint.Increment(),pot_vector);
  double pot_inc_L2(0.);
  pot_vector->Norm2(&pot_inc_L2);

  // safety checks
  if(std::isnan(conc_state_L2) or std::isnan(conc_res_L2) or std::isnan(conc_inc_L2))
    dserror("Calculated vector norm for concentration is not a number!");
  if(std::isinf(conc_state_L2) or std::isinf(conc_res_L2) or std::isinf(conc_inc_L2))
    dserror("Calculated vector norm for concentration is infinity!");
  if(std::isnan(pot_state_L2) or std::isnan(pot_res_L2) or std::isnan(pot_inc_L2))
    dserror("Calculated vector norm for electric potential is not a number!");
  if(std::isinf(pot_state_L2) or std::isinf(pot_res_L2) or std::isinf(pot_inc_L2))
    dserror("Calculated vector norm for electric potential is infinity!");

  // care for the case that nothing really happens in the concentration or electric potential fields
  if(conc_state_L2 < 1.e-5)
    conc_state_L2 = 1.;
  if(pot_state_L2 < 1.e-5)
    pot_state_L2 = 1.;

  // special case: very first iteration step --> solution increment is not yet available
  if(itnum == 1)
  {
    if(mypid == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- con-res ---|-- pot-res ---|-- con-inc ---|-- pot-inc ---|-- con-res-inf ---|" << std::endl;

      // print first line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << pot_res_L2 << "   |      --      |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (      --     ,te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;
    }
  }

  // ordinary case: later iteration steps --> solution increment can be printed and convergence check should be done
  else
  {
    if(mypid == 0)
      // print current line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << pot_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_inc_L2/conc_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << pot_inc_L2/pot_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (ts="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtSolve() << ",te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;

    // convergence check
    if(conc_res_L2 <= ittol_ and conc_inc_L2/conc_state_L2 <= ittol_ and pot_res_L2 <= ittol_ and pot_inc_L2/pot_state_L2 <= ittol_)
    {
      if(mypid == 0)
        // print finish line of convergence table to screen
        std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

      return true;
    }
  }

  // abort iteration when there is nothing more to do --> better robustness
  // absolute tolerance determines whether residual is already zero
  // prevents additional solver calls that will not improve the solution anymore
  if(conc_res_L2 < abstolres_ and pot_res_L2 < abstolres_)
  {
    if(mypid == 0)
      // print finish line of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

    return true;
  }

  // output warning in case maximum number of iteration steps is reached without convergence, and proceed to next time step
  if(itnum == itmax_)
  {
    if (mypid == 0)
    {
      std::cout << "+---------------------------------------------------------------+" << std::endl;
      std::cout << "|       >>>>>> Newton-Raphson iteration did not converge!       |" << std::endl;
      std::cout << "+---------------------------------------------------------------+" << std::endl;
    }

    return true;
  }

  // return maximum residual value for adaptivity of linear solver tolerance
  actresidual = std::max(conc_res_L2,conc_inc_L2/conc_state_L2);
  actresidual = std::max(actresidual,pot_res_L2);
  actresidual = std::max(actresidual,pot_inc_L2/pot_state_L2);

  // proceed with next iteration step
  return false;
} // SCATRA::ConvCheckStrategyStdElch::AbortNonlinIter()


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 02/16 |
 *----------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyS2ILM::AbortNonlinIter(
    const ScaTraTimIntImpl&   scatratimint,   //!< scalar transport time integrator
    double&                   actresidual     //!< return maximum current residual value
    ) const
{
  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current Newton-Raphson iteration step
  const int itnum = scatratimint.IterNum();

  // compute L2 norm of concentration state vector
  double conc_state_L2(0.0);
  scatratimint.Phinp()->Norm2(&conc_state_L2);

  // compute L2 norm of concentration residual vector
  double conc_res_L2(0.);
  scatratimint.Residual()->Norm2(&conc_res_L2);

  // compute infinity norm of concentration residual vector
  double conc_res_inf(0.);
  scatratimint.Residual()->NormInf(&conc_res_inf);

  // compute L2 norm of concentration increment vector
  double conc_inc_L2(0.);
  scatratimint.Increment()->Norm2(&conc_inc_L2);

  // extract meshtying strategy from scalar transport time integrator
  const Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategys2i = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(scatratimint.Strategy());
  if(meshtyingstrategys2i == Teuchos::null)
    dserror("Invalid scalar transport meshtying strategy!");

  // compute L2 norm of Lagrange multiplier state vector
  double lm_state_L2(0.0);
  meshtyingstrategys2i->LM()->Norm2(&lm_state_L2);

  // compute L2 norm of Lagrange multiplier residual vector
  double lm_res_L2(0.);
  meshtyingstrategys2i->LMResidual()->Norm2(&lm_res_L2);

  // compute L2 norm of Lagrange multiplier increment vector
  double lm_inc_L2(0.);
  meshtyingstrategys2i->LMIncrement()->Norm2(&lm_inc_L2);

  // safety checks
  if(std::isnan(conc_state_L2) or std::isnan(conc_res_L2) or std::isnan(conc_inc_L2))
    dserror("Calculated vector norm for concentration is not a number!");
  if(std::isinf(conc_state_L2) or std::isinf(conc_res_L2) or std::isinf(conc_inc_L2))
    dserror("Calculated vector norm for concentration is infinity!");
  if(std::isnan(lm_state_L2) or std::isnan(lm_res_L2) or std::isnan(lm_inc_L2))
    dserror("Calculated vector norm for Lagrange multipliers is not a number!");
  if(std::isinf(lm_state_L2) or std::isinf(lm_res_L2) or std::isinf(lm_inc_L2))
    dserror("Calculated vector norm for Lagrange multipliers is infinity!");

  // care for the case that nothing really happens in the concentration or Lagrange multiplier fields
  if(conc_state_L2 < 1.e-5)
    conc_state_L2 = 1.;
  if(lm_state_L2 < 1.e-5)
    lm_state_L2 = 1.;

  // special case: very first iteration step --> solution increment is not yet available
  if(itnum == 1)
  {
    if(mypid == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- con-res ---|--- lm-res ---|-- con-inc ---|--- lm-inc ---|-- con-res-inf ---|" << std::endl;

      // print first line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << lm_res_L2 << "   |      --      |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (      --     ,te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;
    }
  }

  // ordinary case: later iteration steps --> solution increment can be printed and convergence check should be done
  else
  {
    if(mypid == 0)
      // print current line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << lm_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_inc_L2/conc_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << lm_inc_L2/lm_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (ts="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtSolve() << ",te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;

    // convergence check
    if(conc_res_L2 <= ittol_ and conc_inc_L2/conc_state_L2 <= ittol_ and lm_res_L2 <= ittol_ and lm_inc_L2/lm_state_L2 <= ittol_)
    {
      if(mypid == 0)
        // print finish line of convergence table to screen
        std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

      return true;
    }
  }

  // abort iteration when there is nothing more to do --> better robustness
  // absolute tolerance determines whether residual is already zero
  // prevents additional solver calls that will not improve the solution anymore
  if(conc_res_L2 < abstolres_ and lm_res_L2 < abstolres_)
  {
    if(mypid == 0)
      // print finish line of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

    return true;
  }

  // output warning in case maximum number of iteration steps is reached without convergence, and proceed to next time step
  if(itnum == itmax_)
  {
    if (mypid == 0)
    {
      std::cout << "+---------------------------------------------------------------+" << std::endl;
      std::cout << "|       >>>>>> Newton-Raphson iteration did not converge!       |" << std::endl;
      std::cout << "+---------------------------------------------------------------+" << std::endl;
    }

    return true;
  }

  // return maximum residual value for adaptivity of linear solver tolerance
  actresidual = std::max(conc_res_L2,conc_inc_L2/conc_state_L2);
  actresidual = std::max(actresidual,lm_res_L2);
  actresidual = std::max(actresidual,lm_inc_L2/lm_state_L2);

  // proceed with next iteration step
  return false;
} // SCATRA::ConvCheckStrategyS2ILM::AbortNonlinIter()


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 02/16 |
 *----------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyS2ILMElch::AbortNonlinIter(
    const ScaTraTimIntImpl&   scatratimint,   //!< scalar transport time integrator
    double&                   actresidual     //!< return maximum current residual value
    ) const
{
  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current Newton-Raphson iteration step
  const int itnum = scatratimint.IterNum();

  // compute L2 norm of concentration state vector
  Teuchos::RCP<Epetra_Vector> conc_vector = scatratimint.Splitter()->ExtractOtherVector(scatratimint.Phinp());
  double conc_state_L2(0.0);
  conc_vector->Norm2(&conc_state_L2);

  // compute L2 norm of concentration residual vector
  scatratimint.Splitter()->ExtractOtherVector(scatratimint.Residual(),conc_vector);
  double conc_res_L2(0.);
  conc_vector->Norm2(&conc_res_L2);

  // compute infinity norm of concentration residual vector
  double conc_res_inf(0.);
  conc_vector->NormInf(&conc_res_inf);

  // compute L2 norm of concentration increment vector
  scatratimint.Splitter()->ExtractOtherVector(scatratimint.Increment(),conc_vector);
  double conc_inc_L2(0.);
  conc_vector->Norm2(&conc_inc_L2);

  // compute L2 norm of electric potential state vector
  Teuchos::RCP<Epetra_Vector> pot_vector = scatratimint.Splitter()->ExtractCondVector(scatratimint.Phinp());
  double pot_state_L2(0.0);
  pot_vector->Norm2(&pot_state_L2);

  // compute L2 norm of electric potential residual vector
  scatratimint.Splitter()->ExtractCondVector(scatratimint.Residual(),pot_vector);
  double pot_res_L2(0.);
  pot_vector->Norm2(&pot_res_L2);

  // compute L2 norm of electric potential increment vector
  scatratimint.Splitter()->ExtractCondVector(scatratimint.Increment(),pot_vector);
  double pot_inc_L2(0.);
  pot_vector->Norm2(&pot_inc_L2);

  // extract meshtying strategy from scalar transport time integrator
  const Teuchos::RCP<const SCATRA::MeshtyingStrategyS2I> meshtyingstrategys2i = Teuchos::rcp_dynamic_cast<const SCATRA::MeshtyingStrategyS2I>(scatratimint.Strategy());
  if(meshtyingstrategys2i == Teuchos::null)
    dserror("Invalid scalar transport meshtying strategy!");

  // compute L2 norm of Lagrange multiplier state vector
  double lm_state_L2(0.0);
  meshtyingstrategys2i->LM()->Norm2(&lm_state_L2);

  // compute L2 norm of Lagrange multiplier residual vector
  double lm_res_L2(0.);
  meshtyingstrategys2i->LMResidual()->Norm2(&lm_res_L2);

  // compute L2 norm of Lagrange multiplier increment vector
  double lm_inc_L2(0.);
  meshtyingstrategys2i->LMIncrement()->Norm2(&lm_inc_L2);

  // safety checks
  if(std::isnan(conc_state_L2) or std::isnan(conc_res_L2) or std::isnan(conc_inc_L2))
    dserror("Calculated vector norm for concentration is not a number!");
  if(std::isinf(conc_state_L2) or std::isinf(conc_res_L2) or std::isinf(conc_inc_L2))
    dserror("Calculated vector norm for concentration is infinity!");
  if(std::isnan(pot_state_L2) or std::isnan(pot_res_L2) or std::isnan(pot_inc_L2))
    dserror("Calculated vector norm for electric potential is not a number!");
  if(std::isinf(pot_state_L2) or std::isinf(pot_res_L2) or std::isinf(pot_inc_L2))
    dserror("Calculated vector norm for electric potential is infinity!");
  if(std::isnan(lm_state_L2) or std::isnan(lm_res_L2) or std::isnan(lm_inc_L2))
    dserror("Calculated vector norm for Lagrange multipliers is not a number!");
  if(std::isinf(lm_state_L2) or std::isinf(lm_res_L2) or std::isinf(lm_inc_L2))
    dserror("Calculated vector norm for Lagrange multipliers is infinity!");

  // care for the case that nothing really happens in the concentration, electric potential, or Lagrange multiplier fields
  if(conc_state_L2 < 1.e-5)
    conc_state_L2 = 1.;
  if(pot_state_L2 < 1.e-5)
    pot_state_L2 = 1.;
  if(lm_state_L2 < 1.e-5)
    lm_state_L2 = 1.;

  // special case: very first iteration step --> solution increment is not yet available
  if(itnum == 1)
  {
    if(mypid == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- con-res ---|-- pot-res ---|--- lm-res ---|-- con-inc ---|-- pot-inc ---|--- lm-inc ---|-- con-res-inf ---|" << std::endl;

      // print first line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << pot_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << lm_res_L2 << "   |      --      |      --      |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (      --     ,te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;
    }
  }

  // ordinary case: later iteration steps --> solution increment can be printed and convergence check should be done
  else
  {
    if(mypid == 0)
      // print current line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << pot_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << lm_res_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_inc_L2/conc_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << pot_inc_L2/pot_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << lm_inc_L2/lm_state_L2 << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << conc_res_inf << "       | (ts="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtSolve() << ",te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;

    // convergence check
    if(conc_res_L2 <= ittol_ and conc_inc_L2/conc_state_L2 <= ittol_ and pot_res_L2 <= ittol_ and pot_inc_L2/pot_state_L2 <= ittol_ and lm_res_L2 <= ittol_ and lm_inc_L2/lm_state_L2 <= ittol_)
    {
      if(mypid == 0)
        // print finish line of convergence table to screen
        std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

      return true;
    }
  }

  // abort iteration when there is nothing more to do --> better robustness
  // absolute tolerance determines whether residual is already zero
  // prevents additional solver calls that will not improve the solution anymore
  if(conc_res_L2 < abstolres_ and pot_res_L2 < abstolres_ and lm_res_L2 < abstolres_)
  {
    if(mypid == 0)
      // print finish line of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

    return true;
  }

  // output warning in case maximum number of iteration steps is reached without convergence, and proceed to next time step
  if(itnum == itmax_)
  {
    if (mypid == 0)
    {
      std::cout << "+---------------------------------------------------------------+" << std::endl;
      std::cout << "|       >>>>>> Newton-Raphson iteration did not converge!       |" << std::endl;
      std::cout << "+---------------------------------------------------------------+" << std::endl;
    }

    return true;
  }

  // return maximum residual value for adaptivity of linear solver tolerance
  actresidual = std::max(conc_res_L2,conc_inc_L2/conc_state_L2);
  actresidual = std::max(actresidual,pot_res_L2);
  actresidual = std::max(actresidual,pot_inc_L2/pot_state_L2);
  actresidual = std::max(actresidual,lm_res_L2);
  actresidual = std::max(actresidual,lm_inc_L2/lm_state_L2);

  // proceed with next iteration step
  return false;
} // SCATRA::ConvCheckStrategyS2ILMElch::AbortNonlinIter()


/*----------------------------------------------------------------------*
 | perform convergence check for Newton-Raphson iteration    fang 08/17 |
 *----------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyStdMacroScaleElch::AbortNonlinIter(
    const ScaTraTimIntImpl&   scatratimint,   //!< scalar transport time integrator
    double&                   actresidual     //!< return maximum current residual value
    ) const
{
  // cast scalar transport time integrator
  const ScaTraTimIntElch* const elchtimint = dynamic_cast<const ScaTraTimIntElch*>(&scatratimint);
  if(elchtimint == NULL)
    dserror("Cast of scalar transport time integrator failed!");

  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current Newton-Raphson iteration step
  const int itnum = scatratimint.IterNum();

  // compute L2 norm of state vector associated with electrolyte concentration
  const Teuchos::RCP<Epetra_Vector> vector_conc_el = elchtimint->SplitterMacro()->ExtractVector(scatratimint.Phinp(),0);
  double L2_state_conc_el(0.);
  vector_conc_el->Norm2(&L2_state_conc_el);

  // compute L2 norm of residual vector associated with electrolyte concentration
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.Residual(),0,vector_conc_el);
  double L2_res_conc_el(0.);
  vector_conc_el->Norm2(&L2_res_conc_el);

  // compute infinity norm of residual vector associated with electrolyte concentration
  double inf_res_conc_el(0.);
  vector_conc_el->NormInf(&inf_res_conc_el);

  // compute L2 norm of increment vector associated with electrolyte concentration
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.Increment(),0,vector_conc_el);
  double L2_inc_conc_el(0.);
  vector_conc_el->Norm2(&L2_inc_conc_el);

  // compute L2 norm of state vector associated with electrolyte potential
  const Teuchos::RCP<Epetra_Vector> vector_pot_el = elchtimint->SplitterMacro()->ExtractVector(scatratimint.Phinp(),1);
  double L2_state_pot_el(0.);
  vector_pot_el->Norm2(&L2_state_pot_el);

  // compute L2 norm of residual vector associated with electrolyte potential
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.Residual(),1,vector_pot_el);
  double L2_res_pot_el(0.);
  vector_pot_el->Norm2(&L2_res_pot_el);

  // compute L2 norm of increment vector associated with electrolyte potential
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.Increment(),1,vector_pot_el);
  double L2_inc_pot_el(0.);
  vector_pot_el->Norm2(&L2_inc_pot_el);

  // compute L2 norm of state vector associated with electrode potential
  const Teuchos::RCP<Epetra_Vector> vector_pot_ed = elchtimint->SplitterMacro()->ExtractVector(scatratimint.Phinp(),2);
  double L2_state_pot_ed(0.);
  vector_pot_ed->Norm2(&L2_state_pot_ed);

  // compute L2 norm of residual vector associated with electrode potential
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.Residual(),2,vector_pot_ed);
  double L2_res_pot_ed(0.);
  vector_pot_ed->Norm2(&L2_res_pot_ed);

  // compute L2 norm of increment vector associated with electrode potential
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.Increment(),2,vector_pot_ed);
  double L2_inc_pot_ed(0.);
  vector_pot_ed->Norm2(&L2_inc_pot_ed);

  // safety checks
  if(std::isnan(L2_state_conc_el) or std::isnan(L2_res_conc_el) or std::isnan(L2_inc_conc_el))
    dserror("Calculated vector norm for electrolyte concentration is not a number!");
  if(std::isinf(L2_state_conc_el) or std::isinf(L2_res_conc_el) or std::isinf(L2_inc_conc_el))
    dserror("Calculated vector norm for electrolyte concentration is infinity!");
  if(std::isnan(L2_state_pot_el) or std::isnan(L2_res_pot_el) or std::isnan(L2_inc_pot_el))
    dserror("Calculated vector norm for electrolyte potential is not a number!");
  if(std::isinf(L2_state_pot_el) or std::isinf(L2_res_pot_el) or std::isinf(L2_inc_pot_el))
    dserror("Calculated vector norm for electrolyte potential is infinity!");
  if(std::isnan(L2_state_pot_ed) or std::isnan(L2_res_pot_ed) or std::isnan(L2_inc_pot_ed))
    dserror("Calculated vector norm for electrode potential is not a number!");
  if(std::isinf(L2_state_pot_ed) or std::isinf(L2_res_pot_ed) or std::isinf(L2_inc_pot_ed))
    dserror("Calculated vector norm for electrode potential is infinity!");

  // care for the case that nothing really happens in the electrolyte concentration, electrolyte potential, or electrode potential fields
  if(L2_state_conc_el < 1.e-5)
    L2_state_conc_el = 1.;
  if(L2_state_pot_el < 1.e-5)
    L2_state_pot_el = 1.;
  if(L2_state_pot_ed < 1.e-5)
    L2_state_pot_ed = 1.;

  // special case: very first iteration step --> solution increment is not yet available
  if(itnum == 1)
  {
    if(mypid == 0)
    {
      // print header of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
      std::cout << "|- step/max -|- tol      [norm] -|-- con-res ---|- pot-el-res -|- pot-ed-res -|-- con-inc ---|- pot-el-inc -|- pot-ed-inc -|-- con-res-inf ---|" << std::endl;

      // print first line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_res_conc_el << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_res_pot_el << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_res_pot_ed << "   |      --      |      --      |      --      | "
                << std::setw(10) << std::setprecision(3) << std::scientific << inf_res_conc_el << "       | (      --     ,te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;
    }
  }

  // ordinary case: later iteration steps --> solution increment can be printed and convergence check should be done
  else
  {
    if(mypid == 0)
      // print current line of convergence table to screen
      std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_ << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << ittol_ << "[L_2 ]  | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_res_conc_el << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_res_pot_el << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_res_pot_ed << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_inc_conc_el/L2_state_conc_el << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_inc_pot_el/L2_state_pot_el << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << L2_inc_pot_ed/L2_state_pot_ed << "   | "
                << std::setw(10) << std::setprecision(3) << std::scientific << inf_res_conc_el << "       | (ts="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtSolve() << ",te="
                << std::setw(10) << std::setprecision(3) << std::scientific << scatratimint.DtEle() << ")" << std::endl;

    // convergence check
    if(L2_res_conc_el <= ittol_ and L2_inc_conc_el/L2_state_conc_el <= ittol_ and L2_res_pot_el <= ittol_ and L2_inc_pot_el/L2_state_pot_el <= ittol_ and L2_res_pot_ed <= ittol_ and L2_inc_pot_ed/L2_state_pot_ed <= ittol_)
    {
      if(mypid == 0)
        // print finish line of convergence table to screen
        std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

      return true;
    }
  }

  // abort iteration when there is nothing more to do --> better robustness
  // absolute tolerance determines whether residual is already zero
  // prevents additional solver calls that will not improve the solution anymore
  if(L2_res_conc_el < abstolres_ and L2_res_pot_el < abstolres_ and L2_res_pot_ed < abstolres_)
  {
    if(mypid == 0)
      // print finish line of convergence table to screen
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;

    return true;
  }

  // output warning in case maximum number of iteration steps is reached without convergence, and proceed to next time step
  if(itnum == itmax_)
  {
    if (mypid == 0)
    {
      std::cout << "+---------------------------------------------------------------+" << std::endl;
      std::cout << "|       >>>>>> Newton-Raphson iteration did not converge!       |" << std::endl;
      std::cout << "+---------------------------------------------------------------+" << std::endl;
    }

    return true;
  }

  // return maximum residual value for adaptivity of linear solver tolerance
  actresidual = std::max(L2_res_conc_el,L2_inc_conc_el/L2_state_conc_el);
  actresidual = std::max(actresidual,L2_res_pot_el);
  actresidual = std::max(actresidual,L2_inc_pot_el/L2_state_pot_el);
  actresidual = std::max(actresidual,L2_res_pot_ed);
  actresidual = std::max(actresidual,L2_inc_pot_ed/L2_state_pot_ed);

  // proceed with next iteration step
  return false;
} // SCATRA::ConvCheckStrategyStdMacroScaleElch::AbortNonlinIter


/*--------------------------------------------------------------------------------------------*
 | perform convergence check for outer iteration in partitioned coupling schemes   fang 08/17 |
 *--------------------------------------------------------------------------------------------*/
bool SCATRA::ConvCheckStrategyStdMacroScaleElch::AbortOuterIter(
    const ScaTraTimIntImpl&   scatratimint   //!< scalar transport time integrator
    ) const
{
  // cast scalar transport time integrator
  const ScaTraTimIntElch* const elchtimint = dynamic_cast<const ScaTraTimIntElch*>(&scatratimint);
  if(elchtimint == NULL)
    dserror("Cast of scalar transport time integrator failed!");

  // extract processor ID
  const int mypid = scatratimint.Discretization()->Comm().MyPID();

  // extract current outer iteration step
  const int itnum = scatratimint.IterNumOuter();

  // compute vector norms
  const Teuchos::RCP<Epetra_Vector> vector_conc_el = elchtimint->SplitterMacro()->ExtractVector(scatratimint.Phinp(),0);
  double L2_state_conc_el(0.);
  vector_conc_el->Norm2(&L2_state_conc_el);
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.PhinpInc(),0,vector_conc_el);
  double L2_inc_conc_el(0.);
  vector_conc_el->Norm2(&L2_inc_conc_el);
  const Teuchos::RCP<Epetra_Vector> vector_pot_el = elchtimint->SplitterMacro()->ExtractVector(scatratimint.Phinp(),1);
  double L2_state_pot_el(0.);
  vector_pot_el->Norm2(&L2_state_pot_el);
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.PhinpInc(),1,vector_pot_el);
  double L2_inc_pot_el(0.);
  vector_pot_el->Norm2(&L2_inc_pot_el);
  const Teuchos::RCP<Epetra_Vector> vector_pot_ed = elchtimint->SplitterMacro()->ExtractVector(scatratimint.Phinp(),2);
  double L2_state_pot_ed(0.);
  vector_pot_ed->Norm2(&L2_state_pot_ed);
  elchtimint->SplitterMacro()->ExtractVector(scatratimint.PhinpInc(),2,vector_pot_ed);
  double L2_inc_pot_ed(0.);
  vector_pot_ed->Norm2(&L2_inc_pot_ed);
  if(L2_state_conc_el < 1.e-10)
    L2_state_conc_el = 1.;
  if(L2_state_pot_el < 1.e-10)
    L2_state_pot_el = 1.;
  if(L2_state_pot_ed < 1.e-10)
    L2_state_pot_ed = 1.;

  // print convergence status
  if(mypid == 0)
  {
    std::cout << "|                                                               OUTER ITERATION                                                               |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
    std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itmax_outer_ << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << ittol_outer_ << "[L_2 ]  |      --      |      --      |      --      | "
              << std::setw(10) << std::setprecision(3) << std::scientific << L2_inc_conc_el/L2_state_conc_el << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << L2_inc_pot_el/L2_state_pot_el << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << L2_inc_pot_ed/L2_state_pot_ed << "   |        --        |" << std::endl;
    std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
  }

  // convergence check
  if(L2_inc_conc_el/L2_state_conc_el <= ittol_outer_ and L2_inc_pot_el/L2_state_pot_el <= ittol_outer_ and L2_inc_pot_ed/L2_state_pot_ed <= ittol_outer_)
  {
    if(mypid == 0)
    {
      std::cout << "|                                          OUTER ITERATION LOOP CONVERGED AFTER ITERATION"
                << std::setw(4) << itnum << "/"
                << std::setw(4) << itmax_outer_ << " !                                          |" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
    }

    return true;
  }

  // throw error in case maximum number of iteration steps is reached without convergence
  else if(itnum == itmax_outer_)
  {
    if(mypid == 0)
    {
      std::cout << "|                                       >>>>>> not converged within maximum number of iteration steps!                                        |" << std::endl;
      std::cout << "+------------+-------------------+--------------+--------------+--------------+--------------+--------------+--------------+------------------+" << std::endl;
    }

    dserror("Outer iteration did not converge within maximum number of iteration steps!");

    return true;
  }

  // proceed with next outer iteration step
  return false;
} // SCATRA::ConvCheckStrategyStdMacroScaleElch::AbortOuterIter
