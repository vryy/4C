
#ifdef CCADISCRET

#ifndef MFSI_COUPLINGOPERATOR_H
#define MFSI_COUPLINGOPERATOR_H

#include <Thyra_SingleScalarEuclideanLinearOpBase.hpp>
#include <Thyra_EpetraLinearOp.hpp>

#include "mfsi_coupling.H"

namespace MFSI
{

  /// Thyra based wrapper to a Epetra_Operator that does the
  /// transformation of applied vectors on the fly.
  class CouplingOperator : virtual public Thyra::SingleScalarEuclideanLinearOpBase<double>
  {
  public:

    /** @name Public types */
    //@{

    /** \brief . */
    typedef double Scalar;

    //@}

    CouplingOperator();

    CouplingOperator(const Teuchos::RCP<Epetra_Operator> &op);

    /** \brief Initialize
     *
     * Preconditions:<ul>
     * <li> <tt>op.get() != NULL</tt> (throw <tt>std::invalid_argument</tt>)
     * </ul>
     *
     * Postconditions:<ul>
     * <li> <tt>this->epetra_op().get() == op.get()</tt>
     * <li> [<tt>spmdRange.get() != NULL</tt>] <tt>this->spmdRange().get() == spmdRange.get()</tt>
     * <li> [<tt>spmdDomain.get() != NULL</tt>] <tt>this->spmdDomain().get() == spmdDomain.get()</tt>
     * <li> [<tt>spmdRange.get() == NULL</tt>] <tt>this->spmdRange().get() != NULL</tt>
     * <li> [<tt>spmdDomain.get() == NULL</tt>] <tt>this->spmdDomain().get() != NULL</tt>
     * <li> <tt>this->opSupported(NOTRANS) == true</tt>
     * <li> <tt>this->opSupported(TRNAS) == adjointSupport==EPETRA_OP_ADJOINT_SUPPORTED</tt>
     * </ul>
     */
    void initialize(const Teuchos::RCP<Epetra_Operator> &op);

    /** \brief Set to uninitialized and optionally return the current state.
     *
     * Postconditions:<ul>
     * <li> <tt>this->domain().get() == NULL</tt>
     * <li> <tt>this->range().get() == NULL</tt>
     * </ul>
     */
    void uninitialize(Teuchos::RCP<Epetra_Operator> *op);

    /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the range.
     *
     * Postconditions:<ul>
     * <li> [<tt>this->range().get() != NULL</tt>] <tt>return.get() != NULL</tt>
     * <li> [<tt>this->range().get() == NULL</tt>] <tt>return.get() == NULL</tt>
     * </ul>
     */
    Teuchos::RCP< const Thyra::SpmdVectorSpaceBase<Scalar> > spmdRange() const;

    /** \brief Return a smart pointer to the SpmdVectorSpaceBase object for the domain.
     *
     * Postconditions:<ul>
     * <li> [<tt>this->domain().get() != NULL</tt>] <tt>return.get() != NULL</tt>
     * <li> [<tt>this->domain().get() == NULL</tt>] <tt>return.get() == NULL</tt>
     * </ul>
     */
    Teuchos::RCP< const Thyra::SpmdVectorSpaceBase<Scalar> > spmdDomain() const;

    /** @name Overridden from SingleScalarLinearOpBase */
    //@{

    /** \brief . */
    bool opSupported(Thyra::ETransp M_trans) const;

    //@}

    /** @name Overridden from EuclideanLinearOpBase */
    //@{

    /// Returns <tt>this->spmdRange()</tt>
    Teuchos::RCP< const Thyra::ScalarProdVectorSpaceBase<Scalar> > rangeScalarProdVecSpc() const;

    /// Returns <tt>this->spmdDomain()</tt>
    Teuchos::RCP< const Thyra::ScalarProdVectorSpaceBase<Scalar> > domainScalarProdVecSpc() const;

    /** \brief . */
    void euclideanApply(
      const Thyra::ETransp M_trans,
      const Thyra::MultiVectorBase<Scalar> &X,
      Thyra::MultiVectorBase<Scalar> *Y,
      const Scalar alpha,
      const Scalar beta
      ) const;

    //@}

    /** @name Overridden from LinearOpBase */
    //@{

    /** \brief . */
    Teuchos::RCP<const Thyra::LinearOpBase<Scalar> > clone() const;

    //@}

    /** \name Overridden from Teuchos::Describable */
    //@{

    /** \brief . */
    std::string description() const;

    /** \brief . */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel) const;

    //@}

  protected:

  private:

    ///
    Teuchos::RCP<Thyra::EpetraLinearOp> op_;
  };


  /** \brief Dynamically allocate a CouplingOperator to wrap a
   * Epetra_Operator object.
   *
   * \relates CouplingOperator
   */
  inline
  Teuchos::RCP<CouplingOperator>
  nonconstCouplingOp(
    const Teuchos::RCP<Epetra_Operator> &op
    )
  {
    return Teuchos::rcp(new CouplingOperator(op));
  }

}

#endif
#endif
