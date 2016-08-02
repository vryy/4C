/*!----------------------------------------------------------------------
\file fluid_timint_two_phase.H

\brief TimIntTwoPhase

\level 2
<pre>
\maintainer Magnus Winter & Martin Kronbichler
             {winter,kronbichler}@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15236/-235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_TWO_PHASE_H
#define FLD_TIMINT_TWO_PHASE_H


#include "fluidimplicitintegration.H"


namespace FLD
{

class TimIntTwoPhase : public virtual FluidImplicitTimeInt
{
//  friend class TimIntGenAlpha; //Exists in LOMA, not needed for TimIntTwoPhase though.
public:
  /// Standard Constructor
  TimIntTwoPhase(const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&                      solver,
      const Teuchos::RCP<Teuchos::ParameterList>&              params,
      const Teuchos::RCP<IO::DiscretizationWriter>&            output,
      bool                                                     alefluid = false);

  /// Destructor
  virtual ~TimIntTwoPhase();

  /*!
  \brief initialization

  */
  virtual void Init();

  /*!
  \brief set scalar fields within outer iteration loop

  */
  void SetIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
                         Teuchos::RCP<const Epetra_Vector> scalaram,
                         Teuchos::RCP<const Epetra_Vector> scalardtam,
                         Teuchos::RCP<const Epetra_Vector> fsscalaraf,
                         Teuchos::RCP<const Epetra_Vector> curvatureaf,
                         Teuchos::RCP<const Epetra_MultiVector> gradphiaf,
                         Teuchos::RCP<DRT::Discretization> scatradis);

  /*!
  \brief set scalar fields

  */
  void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp,
                      Teuchos::RCP<const Epetra_Vector> curvatureanp,
                      Teuchos::RCP<const Epetra_MultiVector> gradphinp,
                      Teuchos::RCP<const Epetra_Vector> scatraresidual,
                      Teuchos::RCP<DRT::Discretization> scatradis,
                      const int                whichscalar = -1);

  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element*/
  void SetElementCustomParameter();


  /*!
    \brief Add problem dependent vector.

    */
  virtual void AddProblemDependentVectors();


protected:

  ///Curvature of phi and gradient of phi at node points.
  Teuchos::RCP<Epetra_Vector>              tpf_gradphi_curvaf_;

private:


}; // class TimIntTwoPhase

} // namespace FLD


#endif  // #ifndef FLD_TIMINT_TWO_PHASE_H
