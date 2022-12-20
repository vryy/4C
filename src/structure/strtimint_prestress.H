/*----------------------------------------------------------------------*/
/*! \file
\brief Static Prestress  analysis
\level 2

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
    TimIntPrestress(const Teuchos::ParameterList& timeparams,  //!< ioflags
        const Teuchos::ParameterList& ioparams,                //!< ioflags
        const Teuchos::ParameterList& sdynparams,              //!< input parameters
        const Teuchos::ParameterList& xparams,                 //!< extra flags
        const Teuchos::RCP<DRT::Discretization>& actdis,       //!< current discretisation
        const Teuchos::RCP<LINALG::Solver>& solver,            //!< the solver
        const Teuchos::RCP<LINALG::Solver>& contactsolver,     //!< the solver for contact meshtying
        const Teuchos::RCP<IO::DiscretizationWriter>& output   //!< the output
    );

    void Setup() override;

    //! Update element
    void UpdateStepElement() override;
    //@}

  };  // class TimIntPrestress
}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_STATICS_H
