/*======================================================================*/
/*!
\file timeadaptivity.H

\class TimeAdaptivity

\brief Time step size adaptivity

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* headers */
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef TIMEADAPTIVITY_H
#define TIMEADAPTIVITY_H

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../io/io_drt.H"
#include "../drt_contact/drt_contact_manager.H"
#include "../drt_lib/linalg_mfoperator.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include <NOX.H>
#include <NOX_Epetra.H>
#include "ml_common.h"
#include "ml_include.h"
#include "nlnml_preconditioner.H"
#include "nlnml_linearsystem.H"
#include "nlnml_finelevelnoxinterface.H"

/*----------------------------------------------------------------------*/
/* namespaces */
using namespace std;
using namespace Teuchos;

/*======================================================================*/
/*!
\class TimeAdaptivity

\brief Time step size adaptivity

\author bborn
\date 10/07
*/
class TimeAdaptivity 
{

public:

   /*!
     \brief Norm for error vector
   */
   enum TAErrNorm
   { 
      norm_vague = 0,  // undetermined norm
      norm_l1,  // L1/linear/average norm
      norm_l2,   // L2/Euclidean norm
      norm_rms,  // root mean square (RMS) norm
      norm_inf   // Maximum/infinity norm
   };

   /*!
     \brief Ordinary Constructor with parameters
   */
   explicit TimeAdaptivity
   (
      double timeinitial,
      double timefinal,
      int timestepinitial,
      int timestepfinal,
      double stepsizeinitial,
      //
      double stepsizemax,
      double stepsizemin,
      double sizeratiomax,
      double sizeratiomin,
      double sizeratioscale,
      TAErrNorm errnorm,
      double errtol,
      int errorder,
      int adaptstepmax,
      //
      DRT::Discretization& dis,
      LINALG::Solver& solver,
      IO::DiscretizationWriter& output
   );

   /*!
     \brief Destructor
   */
   virtual ~TimeAdaptivity();

   /*!
      \brief Set parameter list
    */
   //void SetParaList(ParameterList& params);

   /*!
     \brief Set-up integration in time

     see Integrate() underneath
   */
   //virtual void SetUpIntegration
   //(
   //   DRT::Discretization& discret,
   //   LINALG::Solver& solver,
   //   IO::DiscretizationWriter& output
   //) = 0;

   /*!
     \brief Integrate in time
            
     This is the key method here, i.e. the time integration
     algorithm. It has to be implemented by derived classes.
   */
   // Here, a pure virtual method makes only sense, when we have a
   // base class for all time integration schemes: TimeIntegrator.
   // In this case StrGenAlpha is derived from this class.
   // Now, we can demand of any object derived form TimeAdaptivity
   // to offer an Integrate method for any TimeIntegrator.
   //virtual void Integrate(TimeIntegrator::TimeIntegrator& timint) = 0;

   /*
     \brief Indicate error and determine new step size
   */
   void Indicate
   (
      bool& accepted,
      double& stpsiznew
   );

   /*!
     \brief Print error norm string
   */
   string PrintErrNorm() const;

   /*!
     \brief Print time adapting constants
   */
   void PrintConstants(std::ostream& str) const;

   /*!
     \brief Print time adapting variables
   */
   void PrintVariables(std::ostream& str) const;

   /*!
     \brief Print time adapting parameters:TimeIntegrator
   */
   void Print(std::ostream& str) const;


protected:

   // not wanted: = operator
   //TimeAdaptivity operator = (const TimeAdaptivity& old);
   // not wanted: copy constructor
   TimeAdaptivity(const TimeAdaptivity& old);

   // plain time integration constants
   double timeinitial_;  // initial time: t_0
   double timefinal_;  // final time
   int timestepinitial_;  // initial time step index: 0 (often)
   int timestepfinal_;  // maximum time step: n_max
   double stepsizeinitial_;  // initial step size: dt_n

   // adaptive time integration constants
   double stepsizemax_;
   double stepsizemin_;
   double sizeratiomax_;
   double sizeratiomin_;
   double sizeratioscale_;
   enum TAErrNorm errnorm_;
   double errtol_;
   int errorder_;
   int adaptstepmax_;

   // link to discretisation
   DRT::Discretization& discret_;
   LINALG::Solver& solver_;
   IO::DiscretizationWriter& output_;
 
   // more parameters
   //ParameterList* params_;

   // plain time integration variables
   double time_;
   int timestep_;

   // adaptive time integration variables
   double stepsize_;
   RefCountPtr<Epetra_Vector> locdiserrn_;
   int adaptstep_;


}; // class TimeAdaptivity


/*======================================================================*/
/*!
\brief Out stream operator
\author bborn
\date 10/07
*/
std::ostream& operator<<
(
   ostream& str, 
   const TimeAdaptivity& ta
);



#endif  // #ifndef TIMEADAPTIVITY_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
