/*-----------------------------------------------------------*/
/*!
\file loca_nln_group.H

\brief LOCA specific group

\maintainer Michael Hiermeier

\date Nov 16, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef SRC_LOCA_CONTINUATION_LOCA_NLN_GROUP_H_
#define SRC_LOCA_CONTINUATION_LOCA_NLN_GROUP_H_

#include <LOCA_Epetra_Group.H>
#include "../solver_nonlin_nox/nox_nln_group.H"

// forward declaration
namespace Teuchos
{
  class ParameterList;
}  // namespace Teuchos
namespace LOCA
{
  namespace NLN
  {
    class Group : public virtual NOX::NLN::Group, public virtual LOCA::Epetra::Group
    {
     public:
      //! standard constructor enabling most %LOCA support
      Group(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata, Teuchos::ParameterList& grp_options,
          Teuchos::ParameterList& print_params,
          const Teuchos::RCP<LOCA::Epetra::Interface::Required>& ireq,
          NOX::Epetra::Vector& initial_guess, const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys,
          const LOCA::ParameterVector& loca_params);

      //! Constructor with time-dependent interface and shifted linear system
      Group(const Teuchos::RCP<LOCA::GlobalData>& loca_gdata, Teuchos::ParameterList& grp_options,
          Teuchos::ParameterList& print_params,
          const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent>& itime,
          NOX::Epetra::Vector& initial_guess, const Teuchos::RCP<NOX::Epetra::LinearSystem>& linsys,
          const Teuchos::RCP<NOX::Epetra::LinearSystem>& shifted_linsys,
          const LOCA::ParameterVector& loca_params);

      /*!\brief Copy constructor.
       *
       * If type is DeepCopy, takes ownership of valid shared Jacobian and shared
       * preconditioning matrix. */
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! derived
      virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source);
      virtual NOX::Abstract::Group& operator=(const NOX::Epetra::Group& source);

      //! destructor
      virtual ~Group(){};

      //! derived
      Teuchos::RCP<NOX::Abstract::Group> clone(NOX::CopyType type) const;

      //! derived
      NOX::Abstract::Group::ReturnType computeF();

      //! derived
      NOX::Abstract::Group::ReturnType computeFandJacobian();

      //! derived
      void resetIsValid() override;

    };  // class Group
  }     // namespace NLN
}  // namespace LOCA


#endif /* SRC_LOCA_CONTINUATION_LOCA_NLN_GROUP_H_ */
