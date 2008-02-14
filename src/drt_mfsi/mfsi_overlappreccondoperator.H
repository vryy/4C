
#ifdef CCADISCRET

#ifndef MFSI_OVERLAPPRECCONDOPERATOR_H
#define MFSI_OVERLAPPRECCONDOPERATOR_H

#include <Thyra_LinearOpBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_DefaultBlockedLinearOp.hpp>

namespace MFSI
{

  /// internal predictor operator for overlapping monolithic FSI
  /*!

    This is a block Gauss-Seidel preconditioner operator that works on the
    lower triangle of a given 3x3 FSI block matrix. The magic happens in
    apply().

    To the outside world this is a simple linear operator. Internally the
    linear structure, fluid and ale fields are solved, with structure rhs on
    both fluid and ale fields. Thus the coupling from structure to fluid is
    taken into account, whereas the fluid to structure coupling is neglected.

    \author u.kue
    \date 01/08
   */
  class OverlappingPCOperator : public Thyra::LinearOpBase<double,double>
  {
  public:

    /// constructor with field solver objects and 3x3 block matrix
    OverlappingPCOperator(
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > structure_lows,
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > fluid_lows,
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > ale_lows,
      Teuchos::RCP<const Thyra::DefaultBlockedLinearOp<double> > blockFsiOp
      );

    /// Return a smart pointer for the range space for this operator.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > range() const;

    /// Return a smart pointer for the domain space for this operator.
    Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domain() const;

    /// Apply the forward non-conjugate or conjugate linear operator to a multi-vector : Y = alpha*M*X + beta*Y.
    void apply(Thyra::EConj conj,
               const Thyra::MultiVectorBase<double>& X,
               Thyra::MultiVectorBase<double>* Y,
               double alpha=1,
               double beta=0) const;

  private:

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > structure_lows_;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > fluid_lows_;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > ale_lows_;
    Teuchos::RCP<const Thyra::DefaultBlockedLinearOp<double> > blockFsiOp_;
  };
}

#endif
#endif
