/*-----------------------------------------------------------*/
/*! \file

\brief Central differences time integration for solid dynamics


\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STRUCTURE_NEW_EXPL_CENTRDIFF_H
#define STRUCTURE_NEW_EXPL_CENTRDIFF_H

#include "structure_new_expl_generic.H"

namespace STR
{
  namespace EXPLICIT
  {
    /*! \brief Central differences time integration for solid dynamics
     *
     */
    class CentrDiff : public Generic
    {
     public:
      //! constructor
      CentrDiff(){};

      //! destructor
      ~CentrDiff() override = default;

      //! Setup class variables (derived)
      void Setup() override{};

      //! Set state variables (derived)
      void SetState(const Epetra_Vector& x) override{};

      //! Apply the right hand side only (derived)
      bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) override { return false; };

      //! return integration factor (derived)
      double GetIntParam() const override { return -1.0; }

      //! Update configuration after time step (derived)
      void UpdateStepState() override{};

      //! Update everything on element level after time step and after output (derived)
      void UpdateStepElement() override{};
    };
  }  // namespace EXPLICIT
}  // namespace STR


#endif  // STRUCTURE_NEW_EXPL_CENTRDIFF_H
