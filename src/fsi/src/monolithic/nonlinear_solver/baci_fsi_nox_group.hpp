/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of NOX::Group for FSI

\level 1

*/

/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_NOX_GROUP_HPP
#define FOUR_C_FSI_NOX_GROUP_HPP

#include "baci_config.hpp"

#include <NOX_Epetra_Group.H>

BACI_NAMESPACE_OPEN

// forward declarations
namespace FSI
{
  class MonolithicInterface;
}

namespace NOX
{
  namespace FSI
  {
    /// Special NOX group that always sets Jacobian and RHS at the same time.
    class Group : public ::NOX::Epetra::Group
    {
     public:
      Group(BACI::FSI::MonolithicInterface& mfsi,                     ///< monolithic FSI interface
          Teuchos::ParameterList& printParams,                        ///< printing parameters
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i,  ///< NOX interface
          const ::NOX::Epetra::Vector& x,                             ///< initial guess
          const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys     ///< linear system
      );

      /// fetch the known Jacobian and RHS from the field solvers
      void CaptureSystemState();

      /// calculate the RHS vector
      ::NOX::Abstract::Group::ReturnType computeF() override;

      /// calculate the Jacobian matrix
      ::NOX::Abstract::Group::ReturnType computeJacobian() override;

      ::NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p) override;

     private:
      BACI::FSI::MonolithicInterface& mfsi_;
    };
  }  // namespace FSI
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
