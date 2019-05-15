/*-----------------------------------------------------------*/
/*!

\brief Generic class for all explicit time integrators

\maintainer Anh-Tu Vuong

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_EXPL_GENERIC_H_
#define STR_EXPL_GENERIC_H_

#include "str_integrator.H"

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
      virtual ~Generic(){};

      //! Setup (has to be implemented by the derived classes)
      virtual void Setup() = 0;

      //! Set state variables (derived)
      virtual void SetState(const Epetra_Vector& x) = 0;

      //! return integration factor (derived)
      virtual bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) = 0;

      //! (derived)
      double GetIntParam() const = 0;

      //! @name Monolithic update routines
      //! @{
      //! things that should be done before updating (derived)
      virtual void PreUpdate(){/* do nothing for now */};

      //! Update configuration after time step (derived)
      virtual void UpdateStepState() = 0;

      //! Update everything on element level after time step and after output (derived)
      virtual void UpdateStepElement() = 0;

      virtual void OutputStepstate() { ; };

      //! things that should be done after updating (derived)
      virtual void PostUpdate(){/* do nothing for now */};
      //! @}
    };
  }  // namespace EXPLICIT
}  // namespace STR


#endif /* STR_EXPL_GENERIC_H_ */
