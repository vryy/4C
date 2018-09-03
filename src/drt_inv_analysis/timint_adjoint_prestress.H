/*----------------------------------------------------------------------*/
/*!
 \file timint_adjoint_prestress.H
 \brief Adjoint time integration for a quasi static hyperelasticity including prestressing

 <pre>
 Maintainer: Sebastian Kehl
 kehl@mhpc.mw.tum.de
 089 - 289-10361
 </pre>
 */

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef TIMINT_ADJOINT_PRESTRESS_H
#define TIMINT_ADJOINT_PRESTRESS_H

/*----------------------------------------------------------------------*/
/* headers */
#include "timint_adjoint.H"

#include "../drt_inpar/inpar_structure.H"

/*----------------------------------------------------------------------*/
namespace STR
{
  class TimIntAdjointPrestress : public TimIntAdjoint
  {
   public:
    //! constructor
    TimIntAdjointPrestress(Teuchos::RCP<DRT::Discretization> dis);

    //! destructor
    virtual ~TimIntAdjointPrestress() {}

    //! interface to the calling optimization routine to set up the rhs and the primal variables
    /*! \param In
        rhs - rigth hand side for the adjoint equation.
        \param In
        mtime - vector with time for each entry in rhs.
        \param In
        dis - primal solution to the adjoint problem.
        \param In
        mtime - vector with time for each entry in the primal solution.
        \param In
    */
    virtual void SetupAdjoint(Teuchos::RCP<Epetra_MultiVector> rhs, std::vector<double> mtime,
        Teuchos::RCP<Epetra_MultiVector> dis, std::vector<double> time);

    //! dual variables corresponding to the "prestressed adjoint equation"
    virtual Teuchos::RCP<Epetra_MultiVector> ExtractPrestressSolution() { return disdualp_; }

   protected:
    //! solve the linear system
    virtual void Solve();

    //! update state
    virtual void UpdateStepState();

    //! dual variables (lagrange multiplier)
    Teuchos::RCP<Epetra_MultiVector> disdualp_;

    //! disdualn_ at step n
    Teuchos::RCP<Epetra_Vector> disdualnp_;

    //! rhs with additional terms for the prestressing
    Teuchos::RCP<Epetra_Vector> rhsnp_;

    //! stiffness matrix with neumann terms
    Teuchos::RCP<LINALG::SparseOperator> stiffp_;

    //! prestress type
    INPAR::STR::PreStress pstype_;

    //! prestress time
    double pstime_;
  };

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif
