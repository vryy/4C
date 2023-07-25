/*-----------------------------------------------------------*/
/*! \file

\brief Generic class for all explicit time integrators


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef BACI_STRUCTURE_NEW_EXPL_GENERIC_H
#define BACI_STRUCTURE_NEW_EXPL_GENERIC_H

#include "baci_structure_new_integrator.H"

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief A generic fully explicit time integrator for solid dynamics
     *
     * This serves as a base class for all fully explicit time integration schemes.
     */
    class Generic : public STR::Integrator
    {
     public:
      //! constructor
      Generic(){};

      //! destructor
      ~Generic() override = default;

      //! @name Monolithic update routines
      //! @{
      //! things that should be done before updating (derived)
      void PreUpdate() override{/* do nothing for now */};

      virtual void OutputStepstate() { ; };

      //! things that should be done after updating (derived)
      void PostUpdate() override{/* do nothing for now */};
      //! @}
    };
  }  // namespace EXPLICIT
}  // namespace STR


#endif  // STRUCTURE_NEW_EXPL_GENERIC_H
