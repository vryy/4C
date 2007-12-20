#ifdef CCADISCRET

#ifndef MFSI_OVERLAPPRECCOND_H
#define MFSI_OVERLAPPRECCOND_H

#include <Thyra_PreconditionerFactoryBase.hpp>
#include <Thyra_PreconditionerBase.hpp>

#include <Thyra_LinearOpWithSolveFactoryBase.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOperatorDecl.hpp>


namespace MFSI
{

  class OverlappingPCFactory : public Thyra::PreconditionerFactoryBase<double>
  {
  public:

    /** @name Constructors/initializers/accessors */
    //@{

    OverlappingPCFactory(
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  structure,
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  interface,
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  fluid,
      Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > const&  ale
      );

    //@}

    /** @name Overridden from PreconditionerFactoryBase */
    //@{

    /** \brief Check that a <tt>LinearOperator</tt> object is compatible with
     * <tt>*this</tt> factory object.
     */
    bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc ) const;

    /** \brief Create an (uninitialized) <tt>LinearOperator</tt>
     * object to be initalized as the preconditioner later in
     * <tt>this->initializePrecOp()</tt>.
     *
     * Note that on output <tt>return->domain().get()==NULL</tt> may
     * be true which means that the operator is not fully
     * initialized.  In fact, the output operator object is not
     * guaranteed to be fully initialized until after it is passed
     * through <tt>this->initializePrecOp()</tt>.
     */
    Teuchos::RCP<Thyra::PreconditionerBase<double> > createPrec() const;

    /** \brief Initialize the LSCPreconditioner object */
    void initializePrec(const Teuchos::RCP<const Thyra::LinearOpSourceBase<double> >
                        &fwdOpSrc,
                        Thyra::PreconditionerBase<double> *precOp,
                        const Thyra::ESupportSolveUse
                        supportSolveUse = Thyra::SUPPORT_SOLVE_UNSPECIFIED) const;

    /** \brief Uninitialize the LSCPreconditioner object */
    void uninitializePrec(Thyra::PreconditionerBase<double> *prec,
                          Teuchos::RCP<const Thyra::LinearOpSourceBase<double> >  *fwdOpSrc,
                          Thyra::ESupportSolveUse *supportSolveUse = NULL) const;
    //@}

    virtual void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>&);

    virtual Teuchos::RCP<Teuchos::ParameterList> getParameterList();

    virtual Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  private:

    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > structure_;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > interface_;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > fluid_;
    Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > ale_;
  };
}

#endif
#endif
