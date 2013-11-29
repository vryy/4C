/*======================================================================*/
/*!
\file strtimada.H
\brief Time step adaptivity front-end for structural dynamics

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>
*/

/*----------------------------------------------------------------------*/
/* definitions */
#ifndef STRTIMADA_H
#define STRTIMADA_H

/*----------------------------------------------------------------------*/
/* headers */
#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include <Teuchos_StandardParameterEntryValidators.hpp>
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

#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

namespace LINALG
{
  class Solver;
}

namespace IO
{
  class DiscretizationWriter;
}


/*----------------------------------------------------------------------*/
/* belongs to structure dynamics namespace */
namespace STR
{
  // forward declarations
  class TimInt;

  /*====================================================================*/
  /*!
   * \brief Front-end for time step size adaptivity in structural dynamics
   *
   * Time step size adaptivity is based on <e>a posteriori</e> error indicators/
   * estimators. Therefore, the solution obtained with the marching time integrator
   * is compared to one obtained with an auxiliary time integrator. Based on this
   * estimation of the local truncation error, a new time step size is suggested.
   *
   * \author bborn
   * \date 07/08
   */
  class TimAda
  {

  public:

    //! Provide the name as std::string
    static std::string MapKindEnumToString
    (
      const enum INPAR::STR::TimAdaKind term  //!< the enum
    )
    {
      switch (term)
      {
      case INPAR::STR::timada_kind_zienxie :
        return "ZienkiewiczXie";
        break;
      case INPAR::STR::timada_kind_ab2 :
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

    //! Type of adaptivity algorithm
    enum AdaEnum
    {
      ada_vague,  //!< algorithm is unknown
      ada_upward,  //!< of upward type, i.e. auxiliar scheme has \b higher order
                   //!< of accuracy than marching scheme
      ada_downward,  //!< of downward type, i.e. auxiliar scheme has \b lower order
                     //!< of accuracy than marching scheme
      ada_orderequal,  //!< of equal order type, i.e. auxiliar scheme has the \b some
                       //!< order of accuracy like the marching method
      ada_ident  //!< auxiliar scheme is \b identical to marching scheme
    };

    //! Constructor
    TimAda
    (
      const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
      const Teuchos::ParameterList& tap,  //!< adaptive input flags
      Teuchos::RCP<TimInt> tis  //!< marching time integrator
    );

    //! Destructor
    virtual ~TimAda()
    {
      DetachFileStepSize();
    }

    //! @name Actions
    //@{

    //! Integrate in time
    //! This is the key method here, i.e. the time integration
    //!  algorithm.
    int Integrate();

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

    //! Print error norm std::string
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

    //! Prepare output to file(s)
    void PrepareOutputPeriod();

    //! Output to file(s)
    void OutputPeriod();

    //!  Update output periods
    void UpdatePeriod();

    //! Check whether step size output file is attached
    bool AttachedFileStepSize()
    {
      if (outsizefile_) return true;
      else return false;
    }

    //! Attach file handle for step size file #outsizefile_
    void AttachFileStepSize();

    //! Detach file handle for step size file #outsizefile_
    void DetachFileStepSize()
    {
      if (outsizefile_) delete outsizefile_;
    }

    //! Write step size
    void OutputStepSize();

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    virtual enum INPAR::STR::TimAdaKind MethodName() const = 0;

    //! Provide the name as std::string
    std::string MethodTitle() const
    {
      return MapKindEnumToString(MethodName());
    }

    //! Provide local order of accuracy based upon linear test equation
    //! for displacements
    virtual int MethodOrderOfAccuracyDis() const = 0;

    //! Provide local order of accuracy based upon linear test equation
    //! for velocities
    virtual int MethodOrderOfAccuracyVel() const = 0;

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const = 0;

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const = 0;

    //! Provide type of algorithm
    virtual enum AdaEnum MethodAdaptDis() const = 0;

    //@}

    //! use contact solver or not
    bool UseContactSolver() { return false; };

  protected:

    //! not wanted: copy constructor
    TimAda(const TimAda& old);

    //! A revolutionary routine to get --well-- the sign of a number.
    static int Sign
    (
      const double number  //!< a real number
    )
    {
      return (number==0.0) ? 0 : (number>0.0) ? +1 : -1;
    }

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<TimInt> sti_;  //!< marching time integrator
    Teuchos::RCP<DRT::Discretization> discret_;  //!< attached discretisation
    int myrank_;  //!< processor ID
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
    double sizeratiomax_;  //!< maximally permitted increase of current step size
                           //!< relative to last converged one
    double sizeratiomin_;  //!< minimally permitted increase
                           //!< (or maximally permitted decrease)
                           //!< of current step size relative to last converged one
    double sizeratioscale_;  //!< safety factor, should be lower than 1.0
    enum CtrlEnum errctrl_;  //!< type of control, see #CtrlEnum
    enum INPAR::STR::VectorNorm errnorm_;  //!< norm for local error vector
    double errtol_;  //!< target local error tolerance
    int errorder_;  //!< order of local error indication
    int adaptstepmax_;  //!< maximally permitted trials to find tolerable step size
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
    int adaptstep_;  //!< trial counter, cf. #adaptstepmax_
    //@}

    //! @name Output settings
    //@{
    bool outsys_;  //!< do it this step: write system to file
    bool outstr_;  //!< do it this step: write stress/strain to file
    bool outene_;  //!< do it this step: write energy to file
    bool outrest_;  //!< do it this step: write restart data to file
    double outsysperiod_;  //!< print system (dis,vel,acc,...)
                           //!< every given period of time
    double outstrperiod_;  //!< print stress/strain every given
                           //!< period of time
    double outeneperiod_;  //!< print energies every given
                           //!< period of time
    double outrestperiod_;  //!< print restart every given
                            //!< period of time
    int outsizeevery_;  //!< print step size every given step
    double outsystime_;  //!< next output time point for system
    double outstrtime_;  //!< next output time point for stress/strain
    double outenetime_;  //!< next output time point for energy
    double outresttime_;  //!< next output time point for restart
    std::ofstream* outsizefile_;  //!< outputfile for step sizes
    //@}

  }; // class TimAda

}  // namespace STR

/*======================================================================*/
/*!
 * \brief Out stream inserter for STR::TimAda
 *
 * \author bborn \date 10/07
 */
std::ostream& operator<<
(
  std::ostream& str,
  const STR::TimAda& ta
);


/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMADA_H
