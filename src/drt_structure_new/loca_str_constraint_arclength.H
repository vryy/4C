/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX for arc-length
continuation.


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_DRT_STRUCTURE_NEW_LOCA_STR_CONSTRAINT_ARCLENGTH_H_
#define SRC_DRT_STRUCTURE_NEW_LOCA_STR_CONSTRAINT_ARCLENGTH_H_

#include <LOCA_MultiContinuation_ConstraintInterfaceMVDX.H>

// forward declaration
namespace STR
{
  namespace TIMINT
  {
    class BaseDataGlobalState;
  }  // namespace TIMINT
}  // namespace STR
namespace LOCA
{
  class ParameterVector;
  namespace STR
  {
    namespace MultiContinuation
    {
      /*! \brief Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX
       *  for arc-length continuation.
       *
       *  This class implements the arc-length constraint equation for
       *  standard-arc-length continuation:
       *  \f[
       *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
       *         \|x-x^\ast\|^{2} + (p-p^\ast)^2 \|f_{\mathrm{ext}}\|^{2} - (\Delta s)^2.
       *  \f]
       *  In contrast to the pseudo arc-length continuation method (see
       *  LOCA::MultiContinuation::ArcLengthConstraint) we do not keep the tangential
       *  direction (see \f$ v_x\f$ and \f$ v_p \f$) constant after the predictor step,
       *  instead we update them in each Newton iteration.
       *
       *  The derivative of g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s) w.r.t. x is
       *  \f[
       *    \nabla_{x} g = 2 (x - x^\ast)
       *  \f]
       *  and w.r.t. p
       *  \f[
       *    \nabla_{p} g = 2 \|f_{\mathrm{ext}}\|^{2} (p - p^\ast).
       *  \f] */
      class ArcLengthConstraint : public LOCA::MultiContinuation::ConstraintInterfaceMVDX
      {
       public:
        //! constructor
        ArcLengthConstraint();

        //! copy constructor
        ArcLengthConstraint(const ArcLengthConstraint& source, NOX::CopyType type);

        //! destructor
        virtual ~ArcLengthConstraint(){};

        //! LOCA::MultiContinuation::ConstraintInterface methods
        //!@{
        //! Copy
        virtual void copy(const ConstraintInterface& source);

        //! Cloning function
        virtual Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface> clone(
            NOX::CopyType type = NOX::DeepCopy) const;

        //! Return number of constraints
        virtual int numConstraints() const;

        //! Set the solution vector to x.
        virtual void setX(const NOX::Abstract::Vector& x);

        //! Sets parameter indexed by paramID
        virtual void setParam(int paramID, double val);

        //! Sets parameters indexed by paramIDs
        virtual void setParams(
            const std::vector<int>& paramIDs, const NOX::Abstract::MultiVector::DenseMatrix& vals);

        //! Compute constraint residuals
        virtual NOX::Abstract::Group::ReturnType computeConstraints();

        //! Compute derivative of constraints w.r.t. solution vector x
        virtual NOX::Abstract::Group::ReturnType computeDX();

        /*! \brief Compute derivative of constraints w.r.t. supplied parameters.
         *
         *  The first column of \c dgdp should be filled with the constraint
         *  residuals \f$g\f$ if \c isValidG is \c false.  If \c isValidG is
         *  \c true, then the \c dgdp contains \f$g\f$ on input. */
        virtual NOX::Abstract::Group::ReturnType computeDP(const std::vector<int>& paramIDs,
            NOX::Abstract::MultiVector::DenseMatrix& dgdp, bool isValidG);

        //! Return \c true if constraint residuals are valid
        virtual bool isConstraints() const;

        //! Return \c true if derivative of constraint w.r.t. x is valid
        virtual bool isDX() const;

        //! Return constraint residuals
        virtual const NOX::Abstract::MultiVector::DenseMatrix& getConstraints() const;

        /*! \brief Return \c true if solution component of constraint
         * derivatives is zero */
        virtual bool isDXZero() const;

        /*! \brief Perform any preprocessing before a continuation step starts.
         *
         * The \c stepStatus argument indicates whether the previous step was
         * successful.  The default implementation is empty. */
        virtual void preProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus);

        /*! \brief Perform any postprocessing after a continuation step finishes.
         *
         * The \c stepStatus argument indicates whether the step was
         * successful.  The default implementation is empty. */
        virtual void postProcessContinuationStep(LOCA::Abstract::Iterator::StepStatus stepStatus);
        //!@}

        //! MultiVector Constraint interface methods
        //!@{
        /*! \brief Return solution component of constraint derivatives
         *
         *  May return NULL if constraint derivative is zero */
        virtual const NOX::Abstract::MultiVector* getDX() const;
        //!@}

        //! Initialize
        void Init(const Teuchos::RCP<const ::STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
            const Teuchos::RCP<LOCA::ParameterVector>& loca_param_vec_ptr,
            const Teuchos::RCP<NOX::Abstract::Vector>& xvector_ptr);

        //! Setup
        void Setup();

       protected:
        //! check the isinit_ flag
        void CheckInit() const;

        //! check the isinit_ and issetup_
        void CheckInitSetup() const;

        //! reset the isvalid flags
        void resetIsValid();

        //! returns the initialization indicator
        const bool& IsInit() const { return isinit_; };

        //! returns the setup indicator
        const bool& IsSetup() const { return issetup_; };

       protected:
        //! flagg which indicates if init has been called
        bool isinit_;

        //! flagg which indicates if setup has been called
        bool issetup_;

        //! flagg which indicates if the current constraints are valid
        bool isvalid_constraints_;

        //! flagg which indicates if the current derivatives w.r.t. x are valid
        bool isvalid_dx_;

        //! number of constraints
        //! ToDo Currently hard coded.
        int numconstraints_;

        //! Constraint values
        NOX::Abstract::MultiVector::DenseMatrix constraints_;

        //! loca parameter vector
        LOCA::ParameterVector loca_param_vec_;

        //! global state data container (read-only access)
        Teuchos::RCP<const ::STR::TIMINT::BaseDataGlobalState> gstate_ptr_;

        //! solution multi vector
        Teuchos::RCP<NOX::Abstract::MultiVector> xvector_ptr_;
      };
    }  // namespace MultiContinuation
  }    // namespace STR
}  // namespace LOCA


#endif /* SRC_DRT_STRUCTURE_NEW_LOCA_STR_CONSTRAINT_ARCLENGTH_H_ */
