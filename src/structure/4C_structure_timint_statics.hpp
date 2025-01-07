// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_STRUCTURE_TIMINT_STATICS_HPP
#define FOUR_C_STRUCTURE_TIMINT_STATICS_HPP

/*----------------------------------------------------------------------*/
/* headers */
#include "4C_config.hpp"

#include "4C_structure_timint_impl.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace Solid
{
  /*====================================================================*/
  /*!
   * \brief Static analysis
   *
   * This static analysis inside the structural dynamics methods appears
   * slightly displaced, however, it comes in handy in case of
   * fluid-structure-interaction, which is built upon structural
   * dynamics.
   *
   * Regarding this matter, please direct any complaints to Axel Gerstenberger.
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntStatics : public TimIntImpl
  {
   public:
    //! @name Construction
    //@{

    //! Constructor
    TimIntStatics(const Teuchos::ParameterList& timeparams,   //!< ioflags
        const Teuchos::ParameterList& ioparams,               //!< ioflags
        const Teuchos::ParameterList& sdynparams,             //!< input parameters
        const Teuchos::ParameterList& xparams,                //!< extra flags
        std::shared_ptr<Core::FE::Discretization> actdis,     //!< current discretisation
        std::shared_ptr<Core::LinAlg::Solver> solver,         //!< the solver
        std::shared_ptr<Core::LinAlg::Solver> contactsolver,  //!< the solver for contact meshtying
        std::shared_ptr<Core::IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref setup().

    \warning none
    \return bool
    \date 08/16
    \author rauch  */
    void init(const Teuchos::ParameterList& timeparams, const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams, std::shared_ptr<Core::FE::Discretization> actdis,
        std::shared_ptr<Core::LinAlg::Solver> solver) override;

    /*! \brief Setup all class internal objects and members

     setup() is not supposed to have any input arguments !

     Must only be called after init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    void setup() override;

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here except when doing optimization
    void resize_m_step() override { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum Inpar::Solid::DynamicType method_name() const override
    {
      return Inpar::Solid::dyna_statics;
    }

    //! Provide number of steps, a single-step method returns 1
    int method_steps() const override { return 1; }

    //! Give local order of accuracy of displacement part
    int method_order_of_accuracy_dis() const override
    {
      FOUR_C_THROW("Sensible to ask?");
      return 0;
    }

    //! Give local order of accuracy of velocity part
    int method_order_of_accuracy_vel() const override
    {
      FOUR_C_THROW("Sensible to ask?");
      return 0;
    }

    //! Return linear error coefficient of displacements
    double method_lin_err_coeff_dis() const override
    {
      FOUR_C_THROW("Sensible to ask?");
      return 0.0;
    }

    //! Return linear error coefficient of velocities
    double method_lin_err_coeff_vel() const override
    {
      FOUR_C_THROW("Sensible to ask?");
      return 0.0;
    }

    //! return time integration factor
    double tim_int_param() const override { return 0.0; }

    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void predict_const_dis_consist_vel_acc() override;

    //! Consistent predictor with constant velocities,
    //! extrapolated displacements and consistent accelerations
    //! For quasi-static problems this is equivalent to
    //! a linear extrapolation of the displacement field.
    //! In the first step we do TangDis
    void predict_const_vel_consist_acc() override;

    //! Consistent predictor with constant accelerations
    //! and extrapolated velocities and displacements
    //! For quasi-static problems this is equivalent to
    //! a quadratic extrapolation of the displacement field.
    //! In the first step we do TangDis, in the second ConstVel
    void predict_const_acc() override;

    //! Create force residual #fres_ and its stiffness #stiff_
    void evaluate_force_stiff_residual(Teuchos::ParameterList& params) final;

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with solve_relaxation_linear
    void evaluate_force_stiff_residual_relax(Teuchos::ParameterList& params) override;

    //! Evaluate residual #fres_
    void evaluate_force_residual() override;

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double calc_ref_norm_force() override;

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #acc_ from scratch by using the
    //! newly updated #disn_.
    //! This method is necessary for certain predictors (like #predict_const_dis_consist_vel_acc)
    void update_iter_incrementally() override;

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    void update_iter_iteratively() override;

    //! Update step
    void update_step_state() override;

    //! Update element
    void update_step_element() override;

    //! Read and set restart for forces
    void read_restart_force() override;

    //! Write internal and external forces for restart
    void write_restart_force(std::shared_ptr<Core::IO::DiscretizationWriter> output) override;
    //@}

    void apply_dirichlet_bc(const double time,              //!< at time
        std::shared_ptr<Core::LinAlg::Vector<double>> dis,  //!< displacements
                                                            //!< (may be nullptr)
        std::shared_ptr<Core::LinAlg::Vector<double>> vel,  //!< velocities
                                                            //!< (may be nullptr)
        std::shared_ptr<Core::LinAlg::Vector<double>> acc,  //!< accelerations
                                                            //!< (may be nullptr)
        bool recreatemap                                    //!< recreate mapextractor/toggle-vector
                                                            //!< which stores the DOF IDs subjected
                                                            //!< to Dirichlet BCs
        //!< This needs to be true if the bounded DOFs
        //!< have been changed.
        ) override;

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> fext() override
    {
      FOUR_C_THROW("Statics: no external forces at t_n available");
      return nullptr;
    }

    //! Return external force \f$F_{ext,n+1}\f$
    std::shared_ptr<Core::LinAlg::Vector<double>> fext_new() override { return fextn_; }

    //@}

   protected:
    //! equal operator is hidden
    TimIntStatics operator=(const TimIntStatics& old);

    //! copy constructor is hidden
    TimIntStatics(const TimIntStatics& old);

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    std::shared_ptr<Core::LinAlg::Vector<double>> fint_;  //!< internal force at \f$t_n\f$

    std::shared_ptr<Core::LinAlg::Vector<double>> fintn_;  //!< internal force at \f$t_{n+1}\f$

    std::shared_ptr<Core::LinAlg::Vector<double>> fext_;  //!< internal force at \f$t_n\f$

    std::shared_ptr<Core::LinAlg::Vector<double>> fextn_;  //!< external force at \f$t_{n+1}\f$
    //@}

  };  // class TimIntStatics

}  // namespace Solid

/*----------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
