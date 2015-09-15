/*
 * str_timint_explicit.H
 *
 *  Created on: Aug 12, 2015
 *      Author: farah
 */

#ifndef STR_TIMINT_EXPLICIT_H_
#define STR_TIMINT_EXPLICIT_H_


#include "str_timint_base.H"

namespace STR {
namespace TIMINT {
  class Explicit : public Base
  {
  public:
    //! constructor
    Explicit();

    //! destructor
    virtual ~Explicit() {};

    virtual int IntegrateStep(){return 0;};

    //! RHS() is not available for explicit time integrators
    virtual Teuchos::RCP<const Epetra_Vector> RHS()
    {
      dserror("RHS() is not available for explicit time integrators");
      return Teuchos::null;
    }

    //! Prepare time step
    virtual void PrepareTimeStep()
    {
      // update end time \f$t_{n+1}\f$ of this time step to cope with time step size adaptivity
      DataGlobalStatePtr()->SetTimenp(DataGlobalStatePtr()->GetUpdatedTime());

      // things that need to be done before Predict
      PrePredict();

      // TODO prepare contact for new time step
      // PrepareStepContact();

      // things that need to be done after Predict
      PostPredict();

      return;
    }

    //!  Evaluate routine for coupled problems with monolithic approach
    virtual void Evaluate(
      Teuchos::RCP<const Epetra_Vector> disiterinc ///< iterative solution increment
    )
    {
      dserror("All monolithically coupled problems work with implicit time "
          "integration schemes. Thus, calling Evaluate() in an explicit scheme "
          "is not possible.");
    }

    //! Update routine for coupled problems with monolithic approach
    // TODO: do we need this?
    virtual void Update()
    {
      PreUpdate();
      //UpdateStepState();
      //UpdateStepTime();
      //UpdateStepElement();
      PostUpdate();
      return;
    }

    //! Update routine for coupled problems with monolithic approach with time adaptivity
    virtual void Update(const double endtime)
    {
      dserror("Not implemented. No time adaptivity available for explicit time integration.");
    }

    /// has to be renamed either here or PrintStep()
    virtual void Output(bool forced_writerestart)
    {
      PreOutput();
      //OutputStep(forced_writerestart);
      // write Gmsh output
      //writeGmshStrucOutputStep();
      return;
    }

  protected:
    virtual void Setup();
  };
} // namespace TimInt
} // namespace STR

#endif /* STR_TIMINT_EXPLICIT_H_ */
