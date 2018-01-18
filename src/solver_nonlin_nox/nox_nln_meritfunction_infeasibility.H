/*-----------------------------------------------------------*/
/*!
\file nox_nln_meritfunction_infeasibility.H

\brief Implementation of the infeasibility merit function for
       constrained problems. Especially useful for the filter method.

\maintainer Michael Hiermeier

\date Apr 19, 2017

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_SOLVER_NONLIN_NOX_NOX_NLN_MERITFUNCTION_INFEASIBILITY_H_
#define SRC_SOLVER_NONLIN_NOX_NOX_NLN_MERITFUNCTION_INFEASIBILITY_H_

#include <NOX_MeritFunction_Generic.H>

namespace Teuchos {
  class ParameterList;
} // namespace Teuchos
namespace NOX {
  class Utils;
namespace NLN {
namespace MeritFunction {
  enum MeritFctName : int;

  class Infeasibility : public NOX::MeritFunction::Generic
  {
    enum Type
    {
      type_vague,     //!< undefined type
      type_two_norm   //!< use a L2-norm of the infeasibility vector
    };

  public:
    /// constructor
    Infeasibility(
        const Teuchos::ParameterList& params,
        const NOX::Utils& u );

    virtual ~Infeasibility() {};

    //! Computes the merit function, \f$ f(x) \f$.
    virtual double computef(const NOX::Abstract::Group& grp) const;

    /*! Computes the gradient of the merit function, \f$ \nabla f \f$, and
     *  returns the result in the \c result vector. */
    virtual void computeGradient(
        const NOX::Abstract::Group& group,
        NOX::Abstract::Vector& result) const;

    /*! Computes the inner product of the given direction and the gradient
     *  associated with the merit function. Returns the steepest descent
     *  direction in the \c result vector. */
    virtual double computeSlope(
        const NOX::Abstract::Vector& dir,
        const NOX::Abstract::Group& grp) const;

    //! Compute the quadratic model,\f$ m(d) \f$, for the given merit function.
    virtual double computeQuadraticModel(
        const NOX::Abstract::Vector& dir,
        const NOX::Abstract::Group& grp) const;

    /*! Computes the vector in the steepest descent direction that minimizes
     *  the quadratic model. */
    virtual void computeQuadraticMinimizer(
        const NOX::Abstract::Group& grp,
        NOX::Abstract::Vector& result) const;

    //! Returns the name of the merit function.
    virtual const std::string& name() const;

    //! return the name of the merit function as enumerator
    enum MeritFctName Type() const;

  private:
    /// \brief Get a list of currently supported infeasibility merit function types
    /** This list is a sub-list of the merit function enumerator list.
     *
     *  \author hiermeier \date 12/17 */
    std::map<std::string,MeritFctName> GetSupportedTypeList() const;

    /// Set the infeasibility merit function type
    void SetType( const std::string& type_name );

  private:
//    const NOX::Utils& utils_;

    enum MeritFctName infeasibility_type_;

    std::string meritFunctionName_;
  };
}
}
}


#endif /* SRC_SOLVER_NONLIN_NOX_NOX_NLN_MERITFUNCTION_INFEASIBILITY_H_ */
