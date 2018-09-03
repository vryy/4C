/*----------------------------------------------------------------------*/
/*!
\file fsi_nox_group.H

\brief Implementation of NOX::Group for FSI

<pre>
Maintainer: Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

#ifndef FSI_NOX_GROUP_H
#define FSI_NOX_GROUP_H

#include <NOX_Epetra_Group.H>

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
    class Group : public NOX::Epetra::Group
    {
     public:
      Group(::FSI::MonolithicInterface& mfsi,                       ///< monolithic FSI interface
          Teuchos::ParameterList& printParams,                      ///< printing parameters
          const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  ///< NOX interface
          const NOX::Epetra::Vector& x,                             ///< initial guess
          const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys     ///< linear system
      );

      /// fetch the known Jacobian and RHS from the field solvers
      void CaptureSystemState();

      /// calculate the RHS vector
      virtual NOX::Abstract::Group::ReturnType computeF();

      /// calculate the Jacobian matrix
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p);

     private:
      ::FSI::MonolithicInterface& mfsi_;
    };
  }  // namespace FSI
}  // namespace NOX

#endif
