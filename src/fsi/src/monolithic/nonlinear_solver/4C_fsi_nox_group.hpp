// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FSI_NOX_GROUP_HPP
#define FOUR_C_FSI_NOX_GROUP_HPP

#include "4C_config.hpp"

#include <NOX_Epetra_Group.H>

FOUR_C_NAMESPACE_OPEN

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
      Group(FourC::FSI::MonolithicInterface& mfsi,                    ///< monolithic FSI interface
          Teuchos::ParameterList& printParams,                        ///< printing parameters
          const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i,  ///< NOX interface
          const ::NOX::Epetra::Vector& x,                             ///< initial guess
          const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys     ///< linear system
      );

      /// fetch the known Jacobian and RHS from the field solvers
      void capture_system_state();

      /// calculate the RHS vector
      ::NOX::Abstract::Group::ReturnType computeF() override;

      /// calculate the Jacobian matrix
      ::NOX::Abstract::Group::ReturnType computeJacobian() override;

      ::NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p) override;

     private:
      FourC::FSI::MonolithicInterface& mfsi_;
    };
  }  // namespace FSI
}  // namespace NOX

FOUR_C_NAMESPACE_CLOSE

#endif
