/*======================================================================*/
/*!
\file strutimada.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/bornemann
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifdef CCADISCRET
#ifndef STRUTIMADA_H
#define STRUTIMADA_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_io/io.H"

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

#include "strutimint.H"


/*======================================================================*/
/*!
\brief Front-end for time step size adaptivity in structural dynamics

\author bborn
\date 07/08
*/
class StruTimAda
{

public:

  //! List of all available auxiliar schemes
  enum NameEnum
  {
    name_zienxie,  //!< Zienkiewicz-Xie method
    name_ab2  //!< Adams-Bashforth 2nd order
  };

  //! Constructor
  StruTimAda
  (
    const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
    const Teuchos::ParameterList& tap,  //!< adaptive input flags
    Teuchos::RCP<StruTimInt> tis  //!< marching time integrator
  );

  //! Destructor
  virtual ~StruTimAda() { ; }

  //! @name Actions
  //@{

  //! Integrate in time
  //! This is the key method here, i.e. the time integration
  //!  algorithm.
  void Integrate();

  //! Indicate error and determine new step size
  void Indicate
  (
    bool& accepted,  //!< true=accpeted, false=not accepted
    double& stpsiznew  //!< step size prediction for next or
                       //!< step repetition
  );

  //! Provide local discretisation error specific to
  //! marching/auxiliar time integrator pair
  virtual void EvaluateLocalErrorVector() = 0;

  //@}

  //! @name Output
  //@{

  //! Print error norm string
  std::string PrintErrNorm() const;

  //! Print time adapting constants
  void PrintConstants
  (
    std::ostream& str  //!< output stream
  ) const;

  //! Print time adapting variables
  void PrintVariables
  (
    std::ostream& str  //!< output stream
  ) const;

  //! Print time adapting parameters:TimeIntegrator
  void Print
  (
    std::ostream& str  //!< output stream
  ) const;

  //@}

  //! @name Attributes
  //@{

  //! Provide the name
  virtual enum NameEnum MethodName() = 0;

  //! Provide the name as string
  virtual std::string MethodTitle() = 0;

  //! Provide local order of accuracy
  virtual int MethodOrderOfAccuracy() = 0;

  //@}

protected:

  //! not wanted: = operator
  //TimeAdaptivity operator = (const TimeAdaptivity& old);

  //! not wanted: copy constructor
  StruTimAda(const StruTimAda& old);

  //! @name General purpose algorithm members
  //@{
  Teuchos::RCP<StruTimInt> tis_;  //!< marching time integrator
  Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
  int mypid_;  //!< processor ID
  Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver
  Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
  //@}

  //! @name Plain time integration constants
  //@{
  double timeinitial_;  //!< initial time: t_0
  double timefinal_;  //!< final time
  int timestepinitial_;  //!< initial time step index: 0 (often)
  int timestepfinal_;  //!< maximum time step: n_max
  double stepsizeinitial_;  //!< initial step size: dt_n
  //@}

  //! @name Adaptive time integration constants
  //@{
  double stepsizemax_;  //!< maximum time step size (upper limit)
  double stepsizemin_;  //!< minimum time step size (lower limit)
  double sizeratiomax_;
  double sizeratiomin_;
  double sizeratioscale_;
  enum StruTimIntVector::NormEnum errnorm_;
  double errtol_;
  int errorder_;
  int adaptstepmax_;
  //@}

  //! @name plain time integration variables
  //@{
  double time_;  //!< current time \f$t_n\f$
  int timestep_;  //!< current time step \f$n\f$
  //@}

  //! @name Adaptive time integration variables
  //@{
  double stepsizepre_;  //!< previous time step size \f$\Delta t_{n-1}\f$
  double stepsize_;  //!< current time step size \f$\Delta t_n\f$
  Teuchos::RCP<Epetra_Vector> locerrn_;  //!< current local error 
                                         //!< estimation \f$l_n\f$
  int adaptstep_;
  //@}
};


/*======================================================================*/
/*!
\brief Out stream operator
\author bborn
\date 10/07
*/
std::ostream& operator<<
(
  std::ostream& str,
  const StruTimAda::StruTimAda& ta
);


/*----------------------------------------------------------------------*/
#endif  // #ifndef STRUTIMADA_H
#endif  // #ifdef CCADISCRET
