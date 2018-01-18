/*!----------------------------------------------------------------------
\file fluid_timint_two_phase_ost.H
\brief

\level 2

<pre>
\maintainer Christoph Ager & Martin Kronbichler
             {ager,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_TWO_PHASE_OST_H
#define FLD_TIMINT_TWO_PHASE_OST_H


#include "fluid_timint_two_phase.H"
#include "fluid_timint_ost.H"

#include "../linalg/linalg_utils.H"

namespace FLD
{

class TimIntTwoPhaseOst : public TimIntOneStepTheta, public TimIntTwoPhase
{
public:
  /// Standard Constructor
  TimIntTwoPhaseOst(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

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

  Teuchos::RCP<Epetra_Vector>              tpf_gradphi_curvn_;

private:


}; // class TimIntTwoPhaseOst

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_TWO_PHASE_OST_H
