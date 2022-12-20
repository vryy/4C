/*----------------------------------------------------------------------------*/
/*! \file
\brief Definitions of NOX group for non-linear solution techniques
       used within implicit structural time integration
\level 1
*/

/*----------------------------------------------------------------------------*/
#ifndef STRTIMINT_NOXGROUP_H
#define STRTIMINT_NOXGROUP_H

/*----------------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Group.H>
#include "strtimint_impl.H"

/*----------------------------------------------------------------------------*/
/* forward declarations */
namespace STR
{
  class TimIntImpl;
}

/*----------------------------------------------------------------------------*/
/* belongs to NOX name space of Trilinos */
namespace NOX
{
  /*--------------------------------------------------------------------------*/
  namespace STR
  {
    /*==================================================================*/
    /*!
     * \brief Special NOX::Epetra::Group that always //ToDo (mayr) That's not correct anymore.
     *        sets Jacobian and RHS at the same time.
     *
     * \author bborn
     * \date 11/08
     */
    class Group : public NOX::Epetra::Group
    {
     public:
      //! Constructor
      Group(::STR::TimIntImpl& sti,             //!< time integrator
          Teuchos::ParameterList& printParams,  //!< printing parameters
          const Teuchos::RCP<NOX::Epetra::Interface::Required>&
              i,                         //!< basically the NOXified time integrator
          const NOX::Epetra::Vector& x,  //!< current solution vector
          const Teuchos::RCP<NOX::Epetra::LinearSystem>&
              linSys  //!< linear system, matrix and RHS etc.
      );

      //! Compute and store \f$F(x)\f$
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! Compute and store Jacobian \f$\frac{\partial F(x)}{\partial x}\f$
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

     private:
      //! structural time integrator
      //! HINT: Currently this variable is unused within the group
      //::STR::TimIntImpl& sti_;
    };

  }  // namespace STR

}  // namespace NOX

/*----------------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_NOXGROUP_H
