/*----------------------------------------------------------------------*/
/*! \file
\brief Objective function interface

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/

#ifndef OBJECTIVE_FUNCT_H
#define OBJECTIVE_FUNCT_H

#include "Teuchos_RCP.hpp"
#include "Epetra_MultiVector.h"

namespace INVANA
{
  class MatParManager;


  /*!
  \brief objective function interface

  This is a general interface type class for objective functions to be derived from. Since
  an objective function can actually be anything you like this is just an interface. It only
  makes sure you have to implement the "Evaluate" and "EvaluateGradient" routine which should
  have a meaning to any type of objective function.

  */
  class ObjectiveFunct
  {
   public:
    //! standard destructor
    virtual ~ObjectiveFunct(){};

    //! evaluate value of the objective function at time "time"
    virtual void Evaluate(Teuchos::RCP<Epetra_Vector> disp, double time, double& val) = 0;

    //! evaluate gradient of objective function w.r.t displacements/state at time "time"
    virtual void EvaluateGradient(
        Teuchos::RCP<Epetra_Vector> disp, double time, Teuchos::RCP<Epetra_Vector> gradient) = 0;

    //! measured time steps
    virtual std::vector<double> MeasuredTime() = 0;

    //! query for step of measurements of this time step
    /*!
      \param In
      time - timestep to be found with accuracy up to 1.0e-10
    */
    //! returns step of corresponding measurements; otherwise returns -1
    virtual int FindStep(double time) = 0;

    /*! Get scaling factor
     *
     */
    virtual double GetScaleFac() = 0;

  };  // class ObjectiveFunct
}  // namespace INVANA

#endif  // OBJECTIVE_FUNCT_H_
