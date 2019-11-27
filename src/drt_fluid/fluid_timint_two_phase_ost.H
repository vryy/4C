/*-----------------------------------------------------------*/
/*! \file

\brief One-step theta time integration for two-phase flow

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_TWO_PHASE_OST_H
#define FLD_TIMINT_TWO_PHASE_OST_H


#include "fluid_timint_two_phase.H"
#include "fluid_timint_ost.H"

#include "../linalg/linalg_utils_sparse_algebra_math.H"

namespace FLD
{
  class TimIntTwoPhaseOst : public TimIntOneStepTheta, public TimIntTwoPhase
  {
   public:
    /// Standard Constructor
    TimIntTwoPhaseOst(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntTwoPhaseOst();

    /*!
    \brief initialization

    */
    virtual void Init();


    /*!
    \brief set scalar fields at previous time step within outer iteration loop (needed for new OST)

    */
    void SetIterScalarFieldsn(Teuchos::RCP<const Epetra_Vector> curvaturen,
        Teuchos::RCP<const Epetra_MultiVector> gradphin,
        Teuchos::RCP<DRT::Discretization> scatradis);

    /*!
      \brief Add problem dependent vector.

      */
    virtual void AddProblemDependentVectors();

   protected:
    Teuchos::RCP<Epetra_Vector> tpf_gradphi_curvn_;

   private:
  };  // class TimIntTwoPhaseOst

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_TWO_PHASE_OST_H
