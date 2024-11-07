// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_TIMINT_PRESTRESS_HPP
#define FOUR_C_STRUCTURE_TIMINT_PRESTRESS_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_structure_timint_statics.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */

namespace Solid
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
    TimIntPrestress(const Teuchos::ParameterList& timeparams,     //!< ioflags
        const Teuchos::ParameterList& ioparams,                   //!< ioflags
        const Teuchos::ParameterList& sdynparams,                 //!< input parameters
        const Teuchos::ParameterList& xparams,                    //!< extra flags
        const std::shared_ptr<Core::FE::Discretization>& actdis,  //!< current discretisation
        const std::shared_ptr<Core::LinAlg::Solver>& solver,      //!< the solver
        const std::shared_ptr<Core::LinAlg::Solver>&
            contactsolver,  //!< the solver for contact meshtying
        const std::shared_ptr<Core::IO::DiscretizationWriter>& output  //!< the output
    );

    void setup() override;

    //! Update element
    void update_step_element() override;
    //@}

  };  // class TimIntPrestress
}  // namespace Solid

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
