/*-----------------------------------------------------------*/
/*! \file

\brief This file contains the declaration of a time integrator for prestressing

\level 3


*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_STRUCTURE_NEW_IMPL_PRESTRESS_HPP
#define FOUR_C_STRUCTURE_NEW_IMPL_PRESTRESS_HPP

#include "4C_config.hpp"

#include "4C_structure_new_impl_statics.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Solid
{
  namespace IMPLICIT
  {
    /*!
     * \brief The time integrator that also writes dynamic forces to allow a dynamic restart
     */
    class PreStress : public Statics
    {
     public:
      //! constructor
      PreStress();

      void write_restart(
          Core::IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      void update_step_state() override;

      void update_step_element() override;

      /*!
       * \brief During MULF and material iterative prestressing, the displacements resetted after
       * each timestep, which is done here.
       */
      void post_update() override;

      /*!
       * \brief Stop the simulation in case of a converged prestress with the material iterative
       * prestressing method
       *
       * \return true If a converged prestress state is obtained with material iterative
       * prestressing method
       * \return false otherwise
       */
      bool early_stopping() const override;

      void post_time_loop() override;

     private:
      /*!
       * \brief Returns whether the material iterative prestress algorithm reached a converged state
       */
      bool is_material_iterative_prestress_converged() const;

      /// Absolute displacement norm used for the convergence check of the prestress state
      double absolute_displacement_norm_;
    };
  }  // namespace IMPLICIT
}  // namespace Solid

FOUR_C_NAMESPACE_CLOSE

#endif
