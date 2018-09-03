/*-----------------------------------------------------------*/
/*!
\file str_expl_ab2.H

\maintainer Michael Hiermeier

\date Aug 13, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef STR_EXPL_AB2_H_
#define STR_EXPL_AB2_H_

#include "str_expl_generic.H"

namespace STR
{
  namespace EXPLICIT
  {
    class AdamsBashforth2 : public Generic
    {
     public:
      //! constructor
      AdamsBashforth2(){};

      //! destructor
      virtual ~AdamsBashforth2(){};

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



#endif /* STR_EXPL_AB2_H_ */
