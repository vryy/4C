/*----------------------------------------------------------------------*/
/*! \file
\brief Static Prestress  analysis
\level 2

*/

/*----------------------------------------------------------------------*/
#ifndef FOUR_C_STRUCTURE_TIMINT_PRESTRESS_HPP
#define FOUR_C_STRUCTURE_TIMINT_PRESTRESS_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_structure_timint_statics.hpp"

FOUR_C_NAMESPACE_OPEN

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
        const Teuchos::RCP<Core::FE::Discretization>& actdis,  //!< current discretisation
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,      //!< the solver
        const Teuchos::RCP<Core::LinAlg::Solver>&
            contactsolver,  //!< the solver for contact meshtying
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output  //!< the output
    );

    void Setup() override;

    //! Update element
    void UpdateStepElement() override;
    //@}

  };  // class TimIntPrestress
}  // namespace STR

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
