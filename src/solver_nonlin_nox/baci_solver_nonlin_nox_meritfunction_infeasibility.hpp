/*-----------------------------------------------------------*/
/*! \file

\brief Implementation of the infeasibility merit function for
       constrained problems. Especially useful for the filter method.



\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_SOLVER_NONLIN_NOX_MERITFUNCTION_INFEASIBILITY_HPP
#define BACI_SOLVER_NONLIN_NOX_MERITFUNCTION_INFEASIBILITY_HPP

#include "baci_config.hpp"

#include "baci_solver_nonlin_nox_forward_decl.hpp"

#include <NOX_MeritFunction_Generic.H>

#include <map>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace NLN
  {
    namespace MeritFunction
    {
      enum MeritFctName : int;

      class Infeasibility : public ::NOX::MeritFunction::Generic
      {
        enum Type
        {
          type_vague,    //!< undefined type
          type_two_norm  //!< use a L2-norm of the infeasibility vector
        };

       public:
        /// constructor
        Infeasibility(const Teuchos::ParameterList& params, const ::NOX::Utils& u);

        //! Computes the merit function, \f$ f(x) \f$.
        double computef(const ::NOX::Abstract::Group& grp) const override;

        /*! Computes the gradient of the merit function, \f$ \nabla f \f$, and
         *  returns the result in the \c result vector. */
        void computeGradient(
            const ::NOX::Abstract::Group& group, ::NOX::Abstract::Vector& result) const override;

        /*! Computes the inner product of the given direction and the gradient
         *  associated with the merit function. Returns the steepest descent
         *  direction in the \c result vector. */
        double computeSlope(
            const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const override;

        //! Compute the quadratic model,\f$ m(d) \f$, for the given merit function.
        double computeQuadraticModel(
            const ::NOX::Abstract::Vector& dir, const ::NOX::Abstract::Group& grp) const override;

        /*! Computes the vector in the steepest descent direction that minimizes
         *  the quadratic model. */
        void computeQuadraticMinimizer(
            const ::NOX::Abstract::Group& grp, ::NOX::Abstract::Vector& result) const override;

        //! Returns the name of the merit function.
        const std::string& name() const override;

        //! return the name of the merit function as enumerator
        enum MeritFctName Type() const;

       private:
        /// \brief Get a list of currently supported infeasibility merit function types
        /** This list is a sub-list of the merit function enumerator list.
         *
         *  \author hiermeier \date 12/17 */
        std::map<std::string, MeritFctName> GetSupportedTypeList() const;

        /// Set the infeasibility merit function type
        void SetType(const std::string& type_name);

       private:
        //    const ::NOX::Utils& utils_;

        enum MeritFctName infeasibility_type_;

        std::string meritFunctionName_;
      };
    }  // namespace MeritFunction
  }    // namespace NLN
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif  // SOLVER_NONLIN_NOX_MERITFUNCTION_INFEASIBILITY_H
