/*======================================================================*/
/*!
\file strtimada.H

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
#ifndef STRTIMADA_H
#define STRTIMADA_H

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

#include "strtimint.H"

/*----------------------------------------------------------------------*/
/* belongs to structure dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Front-end for time step size adaptivity in structural dynamics
   *
   * \author bborn
   * \date 07/08
   */
  class TimAda
  {

  public:

    //! List of all available auxiliar schemes
    enum NameEnum
    {
      name_zienxie,  //!< Zienkiewicz-Xie method
      name_ab2  //!< Adams-Bashforth 2nd order
    };

    //! Provide the name as string
    static std::string MapNameEnumToString
    (
      const enum NameEnum term  //!< the enum
    )
    {
      switch (term)
      {
      case name_zienxie :
        return "ZienkiewiczXie";
        break;
      case name_ab2 :
        return "AdamsBashforth2";
        break;
      default :
        dserror("Cannot cope with name enum %d", term);
        return "";
        break;
      }
    }

    //! List type of local error control
    enum CtrlEnum
    {
      ctrl_dis,  //!< check only displacements
      ctrl_vel,  //!< check only velocities
      ctrl_dis_and_vel  //!< check displacements and velocities
    };

    //! Constructor
    TimAda
    (
      const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
      const Teuchos::ParameterList& tap,  //!< adaptive input flags
      Teuchos::RCP<TimInt> tis  //!< marching time integrator
    );

    //! Destructor
    virtual ~TimAda() { ; }

    //! @name Actions
    //@{

    //! Integrate in time
    //! This is the key method here, i.e. the time integration
    //!  algorithm.
    void Integrate();

    //! Make one step with auxiliar scheme
    //!
    //! Afterwards, we are going to hold 
    //! \f$D_{n+1}^{AUX}\f$ on #locdiserrn_
    //! and \f$V_{n+1}^{AUX}\f$ on #locvelerrn_
    virtual void IntegrateStepAuxiliar() = 0;

    //! Indicate error and determine new step size
    void Indicate
    (
      bool& accepted,  //!< true=accpeted, false=not accepted
      double& stpsiznew  //!< step size prediction for next or
                         //!< step repetition
    );

    //! Provide local discretisation error vector
    //! of displacements/velocities specific to
    //! marching/auxiliar time integrator pair
    void EvaluateLocalErrorDis();

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

    //! Modify step size to hit precisely output period
    void SizeForOutput();

    //! Output to file(s)
    void OutputPeriod();

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    virtual enum NameEnum MethodName() = 0;

    //! Provide the name as string
    std::string MethodTitle()
    {
      return MapNameEnumToString(MethodName());
    }

    //! Provide local order of accuracy based upon linear test equation
    virtual int MethodOrderOfAccuracy() = 0;

    //@}

  protected:

    //! not wanted: copy constructor
    TimAda(const TimAda& old);

    //! A revolutionary routine to get --well-- the sign of a number.
    static int Sign
    (
      const double number
    )
    {
      return (number==0.0) ? 0 : (number>0.0) ? +1 : -1;
    }

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<TimInt> sti_;  //!< marching time integrator
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
    int mypid_;  //!< processor ID
    Teuchos::RCP<LINALG::Solver> solver_;  //!< linear algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    //@}



    //! @name Plain time integration constants
    //@{
    double timeinitial_;  //!< initial time: t_0
    double timefinal_;  //!< final time
    int timedirect_;  //!< +1: in positive, -1: in negative time direction
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
    enum CtrlEnum errctrl_;  //!< type of check
    enum TimIntVector::NormEnum errnorm_;
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
    Teuchos::RCP<Epetra_Vector> locerrdisn_;  //!< current local disp. error 
                                              //!< estimation \f$l_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> locerrveln_;  //!< current local vel. error 
                                              //!< estimation \f$\dot{l}_{n+1}\f$
    int adaptstep_;
    //@}

    //! @name Output settings
    //@{
    bool outsys_;  //!< do it this step: write system to file
    bool outstr_;  //!< do it this step: write stress/strain to file
    bool outrest_;  //!< do it this step: write restart data to file
    double outsysperiod_;  //!< print system (dis,vel,acc,...) 
                           //!< every given period of time
    double outstrperiod_;  //!< print stress/strain every given 
                           //!< period of time
    double outrestperiod_;  //!< print restart every given 
                            //!< period of time
    double outsystime_;  //!< next output time point for system
    double outstrtime_;  //!< next output time point for stress/strain
    double outresttime_;  //!< next output time point for restart
    //@}

  }; // class TimAda

}  // namespace STR

/*======================================================================*/
/*!
\brief Out stream operator
\author bborn
\date 10/07
*/
std::ostream& operator<<
(
  std::ostream& str,
  const STR::TimAda::TimAda& ta
);


/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMADA_H
#endif  // #ifdef CCADISCRET
