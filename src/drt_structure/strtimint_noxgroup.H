/*----------------------------------------------------------------------*/
/*!
\file strtimint_noxgroup.H
\brief Definitions of NOX group for non-linear solution techniques
       used within implicit structural time integration

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_NOXGROUP_H
#define STRTIMINT_NOXGROUP_H

/*----------------------------------------------------------------------*/
/* headers */
#include <NOX_Epetra_Group.H>
#include "strtimint_impl.H"

/*----------------------------------------------------------------------*/
/* foward declarations */
namespace STR
{
  class TimIntImpl;
}

/*----------------------------------------------------------------------*/
/* belongs to NOX name space of Trilinos */
namespace NOX
{
  /*--------------------------------------------------------------------*/
  /* belongs to STR name space of Baci/DRT */
  namespace STR
  {

    /*==================================================================*/
    /*!
     * \brief Special NOX::Epetra::NOX group that always 
     *        sets Jacobian and RHS at the same time.
     *
     * \author bborn
     * \date 11/08
     */
    class Group : public NOX::Epetra::Group
    {
    public:
      
      //! Constructor
      Group
      (
        ::STR::TimIntImpl& sti,  //!< time integrator
        Teuchos::ParameterList& printParams,  //!< printing parameters
        const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  //!< basically time integrator 
        const NOX::Epetra::Vector& x,  //!< current solution vector
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys  //!< linear system, i.e. matrix and RHS
      );
      
      virtual NOX::Abstract::Group::ReturnType computeF();
      
      virtual NOX::Abstract::Group::ReturnType computeJacobian();
      
      NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& p);
      
    private:
      
      //! structural time integrator
      //! HINT: Currently this not unused within the group
      ::STR::TimIntImpl& sti_;
    };

  }  // namespace STR

}  // namespace NOX

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_NOXGROUP_H
#endif  // #ifdef CCADISCRET
