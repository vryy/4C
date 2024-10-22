// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_THERMO_TIMINT_EXPLEULER_HPP
#define FOUR_C_THERMO_TIMINT_EXPLEULER_HPP


/*----------------------------------------------------------------------*
 | headers                                                   dano 01/12 |
 *----------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_thermo_timint_expl.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                    bborn 01/12 |
 *----------------------------------------------------------------------*/
namespace Thermo
{
  /*====================================================================*/
  /*!
   * \brief Forward Euler
   *        explicit time integrator
   * \author dano
   * \date 01/12
   */
  class TimIntExplEuler : public TimIntExpl
  {
   public:
    //! @name Life
    //@{

    //! Constructor
    TimIntExplEuler(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& tdynparams,            //!< input parameters
        const Teuchos::ParameterList& xparams,               //!< extra flags
        Teuchos::RCP<Core::FE::Discretization> actdis,       //!< current discretisation
        Teuchos::RCP<Core::LinAlg::Solver> solver,           //!< the solver
        Teuchos::RCP<Core::IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

    //! Empty constructor
    TimIntExplEuler() : TimIntExpl() { ; }

    //! Copy constructor
    TimIntExplEuler(const TimIntExplEuler& old) : TimIntExpl(old) { ; }

    //! Resize #TimIntMStep<T> multi-step quantities
    void resize_m_step() override { FOUR_C_THROW("not a multistep method"); }

    //@}

    //! @name Actions
    //@{

    //! Do time integration of single step
    void integrate_step() override;

    //! Update configuration after time step
    //!
    //! Thus the 'last' converged is lost and a reset of the time step
    //! becomes impossible. We are ready and keen awating the next time step.
    void update_step_state() override;

    //! Update Element
    void update_step_element() override;

    //@}

    //! @name Attribute access functions
    //@{

    //! Return time integrator name
    enum Inpar::Thermo::DynamicType method_name() const override
    {
      return Inpar::Thermo::dyna_expleuler;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int method_steps() override { return 1; }

    //! Give local order of accuracy of temperature part
    int method_order_of_accuracy() override { return 1; }

    //! Return linear error coefficient of temperatures
    double method_lin_err_coeff() override
    {
      FOUR_C_THROW("no time adaptivity possible.");
      return 0.0;
    }

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> fext() override { return fextn_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> fext_new()
    {
      FOUR_C_THROW("FextNew() not available in ExplEuler");
      return Teuchos::null;
    }

    //! Read and set restart for forces
    void read_restart_force() override;

    //! Write internal and external forces for restart
    void write_restart_force(Teuchos::RCP<Core::IO::DiscretizationWriter> output) override;

    //@}

   protected:
    //! @name Global forces at \f$t_{n+1}\f$
    //@{
    Teuchos::RCP<Core::LinAlg::Vector<double>> fextn_;  //!< external force
                                                        //!< \f$F_{int;n+1}\f$
    Teuchos::RCP<Core::LinAlg::Vector<double>> fintn_;  //!< internal force
                                                        //!< \f$F_{int;n+1}\f$
    //@}

  };  // class TimIntExplEuler

}  // namespace Thermo

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
