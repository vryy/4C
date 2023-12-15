/*----------------------------------------------------------------------*/
/*! \file

\brief Calculates vector extrapolation via epsilon algorithm.


\level 1
*/
/*---------------------------------------------------------------------*/

#ifndef BACI_FSI_NOX_EPSILON_H
#define BACI_FSI_NOX_EPSILON_H

#include <NOX_Direction_Generic.H>  // base class
#include <NOX_Utils.H>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace FSI
  {
    //! %Calculates vector extrapolation via epsilon algorithm.
    class EpsilonExtrapolation : public ::NOX::Direction::Generic
    {
     public:
      //! Constructor
      EpsilonExtrapolation(const Teuchos::RCP<::NOX::Utils>& utils, Teuchos::ParameterList& params);


      // derived
      bool reset(
          const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) override;

      // derived
      bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
          const ::NOX::Solver::Generic& solver) override;

      // derived
      bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& group,
          const ::NOX::Solver::LineSearchBased& solver) override;

     private:
      //! Print error message and throw error
      void throwError(const std::string& functionName, const std::string& errorMsg);

     private:
      //! Printing Utils
      Teuchos::RCP<::NOX::Utils> utils_;

      int kmax_;
      double omega_;
      // double eps_;
      int maxcol_;
    };


  }  // namespace FSI
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
