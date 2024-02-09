/*-----------------------------------------------------------*/
/*! \file

\brief This file contains the declaration of a time integrator for prestressing

\level 3


*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_IMPL_PRESTRESS_HPP
#define BACI_STRUCTURE_NEW_IMPL_PRESTRESS_HPP

#include "baci_config.hpp"

#include "baci_structure_new_impl_statics.hpp"

BACI_NAMESPACE_OPEN

namespace STR
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

      void WriteRestart(
          IO::DiscretizationWriter& iowriter, const bool& forced_writerestart) const override;

      void UpdateStepState() override;

      void UpdateStepElement() override;

      /*!
       * \brief During MULF and material iterative prestressing, the displacements resetted after
       * each timestep, which is done here.
       */
      void PostUpdate() override;

      /*!
       * \brief Stop the simulation in case of a converged prestress with the material iterative
       * prestressing method
       *
       * \return true If a converged prestress state is obtained with material iterative
       * prestressing method
       * \return false otherwise
       */
      bool EarlyStopping() const override;

      void PostTimeLoop() override;

     private:
      /*!
       * \brief Returns whether the material iterative prestress algorithm reached a converged state
       */
      bool IsMaterialIterativePrestressConverged() const;

      /// Absolute displacement norm used for the convergence check of the prestress state
      double absoluteDisplacementNorm_;
    };
  }  // namespace IMPLICIT
}  // namespace STR


BACI_NAMESPACE_CLOSE

#endif  // STRUCTURE_NEW_IMPL_PRESTRESS_H
