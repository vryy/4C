/*----------------------------------------------------------------------*/
/*!
\file strtimint_statics.H
\brief Static Prestress  analysis

<pre>
Maintainer: Sebastian Kehl
            kehl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_PRESTRESS_H
#define STRTIMINT_PRESTRESS_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_statics.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */

namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Static Prestress analysis
   *
   * This is the prestress version of the static analysis inside the structural dynamics section.
   *
   * Regarding this matter, please direct any complaints to Michael Gee.
   *
   * \authors amaier/kehl
   * \date 03/12
   */

  class TimIntPrestress : public TimIntStatics
  {

  public:

    //! @name Construction
    //@{
    //! Constructor
    TimIntPrestress
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,  //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Update element
    void UpdateStepElement();

    //! Write external and internal forces for dynamic restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);
    //@}

  };  // class TimIntPrestress
}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_STATICS_H
