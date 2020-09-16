/*-----------------------------------------------------------*/
/*! \file

\brief This file contains the declaration of a time integrator for prestressing

\level 3


*/
/*-----------------------------------------------------------*/

#ifndef STR_IMPL_PRESTRESS_H_
#define STR_IMPL_PRESTRESS_H_

#include "str_impl_statics.H"

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

     private:
      /// Absolute displacement norm used for the convergence check of the prestress state
      double absoluteDisplacementNorm_;
    };
  }  // namespace IMPLICIT
}  // namespace STR


#endif /* STR_IMPL_PRESTRESS_H_ */
