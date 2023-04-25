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
      virtual ~CentrDiff(){};

      //! Setup class variables (derived)
      virtual void Setup(){};

      //! Set state variables (derived)
      virtual void SetState(const Epetra_Vector& x){};

      //! Apply the right hand side only (derived)
      virtual bool ApplyForce(const Epetra_Vector& x, Epetra_Vector& f) { return false; };

      //! return integration factor (derived)
      virtual double GetIntParam() const { return -1.0; }

      //! Update configuration after time step (derived)
      virtual void UpdateStepState(){};

      //! Update everything on element level after time step and after output (derived)
      virtual void UpdateStepElement(){};
    };
  }  // namespace EXPLICIT
}  // namespace STR


#endif  // STRUCTURE_NEW_EXPL_CENTRDIFF_H
