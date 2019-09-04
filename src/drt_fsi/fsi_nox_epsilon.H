/*----------------------------------------------------------------------*/
/*! \file

\brief Calculates vector extrapolation via epsilon algorithm.

\maintainer Matthias Mayr

\level 1
*/
/*---------------------------------------------------------------------*/

#ifndef FSI_NOX_EPSILON_H
#define FSI_NOX_EPSILON_H

#include <NOX_Direction_Generic.H>  // base class
#include <NOX_Utils.H>
#include <Teuchos_RCP.hpp>

namespace NOX
{
  namespace FSI
  {
    //! %Calculates vector extrapolation via epsilon algorithm.
    class EpsilonExtrapolation : public NOX::Direction::Generic
    {
     public:
      //! Constructor
      EpsilonExtrapolation(const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params);

      //! Destructor
      ~EpsilonExtrapolation();

      // derived
      bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

      // derived
      bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
          const NOX::Solver::Generic& solver);

      // derived
      virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
          const NOX::Solver::LineSearchBased& solver);

     private:
      //! Print error message and throw error
      void throwError(const std::string& functionName, const std::string& errorMsg);

     private:
      //! Printing Utils
      Teuchos::RCP<NOX::Utils> utils_;

      int kmax_;
      double omega_;
      // double eps_;
      int maxcol_;
    };


  }  // namespace FSI
}  // namespace NOX

#endif
