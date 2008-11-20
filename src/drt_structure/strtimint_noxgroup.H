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
  //! NOX::STR holds local extensions to NOX for structural time integration
  namespace STR
  {

    /*==================================================================*/
    /*!
     * \brief Special NOX::Epetra::Group that always 
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
        const Teuchos::RCP<NOX::Epetra::Interface::Required>& i,  //!< basically the NOXified 
                                                                  //!< time integrator 
        const NOX::Epetra::Vector& x,  //!< current solution vector
        const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys  //!< linear system, 
                                                               //!< i.e. matrix and RHS etc.
      );

      //! Compute and store \f$F(x)\f$
      virtual NOX::Abstract::Group::ReturnType computeF();
      
      //! Compute and store Jacobian \f$\frac{\partial F(x)}{\partial x}\f$
      virtual NOX::Abstract::Group::ReturnType computeJacobian();
      
      //! Compute the Newton direction, using parameters for the linear solvexs
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
