/*======================================================================*/
/*!
\brief Time step adaptivity front-end for structural dynamics
\level 1

\maintainer Matthias Mayr
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

#include "../drt_inpar/inpar_structure.H"

// forward declarations
namespace DRT
{
  class Discretization;
}

class Epetra_Vector;

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
    static std::string MapKindEnumToString(const enum INPAR::STR::TimAdaKind term  //!< the enum
    )
    {
      switch (term)
      {
        case INPAR::STR::timada_kind_zienxie:
          return "ZienkiewiczXie";
          break;
        case INPAR::STR::timada_kind_ab2:
          return "AdamsBashforth2";
          break;
        case INPAR::STR::timada_kind_expleuler:
          return "ExplicitEuler";
          break;
        case INPAR::STR::timada_kind_centraldiff:
          return "CentralDifference";
          break;
        default:
          dserror("Cannot cope with name enum %d", term);
          return "";
          break;
      }
      return "";  // make compiler happy
    }

    //! List type of local error control
    enum CtrlEnum
    {
      ctrl_dis,         //!< check only displacements
      ctrl_vel,         //!< check only velocities
      ctrl_dis_and_vel  //!< check displacements and velocities
    };

    //! Type of adaptivity algorithm
    enum AdaEnum
    {
      ada_vague,     //!< algorithm is unknown
      ada_upward,    //!< of upward type, i.e. auxiliary scheme has \b higher order of accuracy than
                     //!< marching scheme
      ada_downward,  //!< of downward type, i.e. auxiliary scheme has \b lower order of accuracy
                     //!< than marching scheme
      ada_orderequal,  //!< of equal order type, i.e. auxiliary scheme has the \b same order of
                       //!< accuracy like the marching method
      ada_ident        //!< auxiliary scheme is \b identical to marching scheme
    };

    //! Constructor
    TimAda(const Teuchos::ParameterList& sdyn,  //!< TIS input parameters
        const Teuchos::ParameterList& tap,      //!< adaptive input flags
        Teuchos::RCP<TimInt> tis                //!< marching time integrator
    );

    //! Destructor
    virtual ~TimAda() { ; }

    //! @name Actions
    //@{

    /*! \brief Integrate in time
     * This is the key method here, i.e. the time integration algorithm.
     */
    int Integrate();

    /*! Finalize the class initialization
     * Merge() and ResizeMStep() need to be called after(!) both Init()
     * and Setup() have been called on both the marching time integrator
     * and the auxiliary time integrator if existing (popp 01/2017).
     */
    virtual void Init(Teuchos::RCP<TimInt>& sti) = 0;

    /*! \brief Make one step with auxiliary scheme
     *
     *  Afterwards, the auxiliary solutions are stored in the local error
     *  vectors, ie:
     *  - \f$D_{n+1}^{AUX}\f$ in #locdiserrn_
     *  - \f$V_{n+1}^{AUX}\f$ in #locvelerrn_
     */
    virtual void IntegrateStepAuxiliar() = 0;

    //! Indicate error and determine new step size
    void Indicate(bool& accepted,  //!< true=accepted, false=not accepted
        double& stpsiznew          //!< step size prediction for next step or step repetition
    );

    /*! \brief Compute local discretisation error
     *
     *  Compute the local discretisation error vector of displacements/velocities
     *  specific to marching/auxiliary time integrator pair.
     *
     *  \note Solution of auxiliary step is already stored in ##locdiserrn_/#locvelerrn_,
     *  so we just need to subtract the solution of the marching TIS.
     */
    void EvaluateLocalErrorDis();

    /*! \brief Calculate time step size
     *
     *  Using the ratio of the desired tolerance \f$tol\f$ (#errtol_) to the
     *  estimated local discretization error, an optimal scaling
     *  factor \f$\kappa_{opt}\f$ is computed, such that the user given error
     *  tolerance is met 'exactly'.
     *  \f[
     *    \kappa_{opt} = \left(\frac{tol}{\vert error\vert}\right)^\frac{1}{p+1}
     *  \f]
     *  To reduce the number of time step repetitions, the scaling factor is
     *  reduced by a safety factor \f$\kappa_{safe} \in [0, 1]\f$ (#sizeratioscale_)
     *  to hopefully keep the achieved local discretization error a little bit
     *  below the tolerance.
     *
     *  Starting with the current time step size \f$\Delta t_{curr}\f$ (#stepsize_),
     *  the new time step size is computed as
     *  \f[
     *    \Delta t_{new} = \kappa_{opt} \cdot \kappa_{safe} \cdot \Delta t_{curr}
     *  \f]
     *
     *  Now, we update the actual scaling factor
     *  \f$\kappa_{eff} = \Delta t_{new} / \Delta t^{n-1}\f$,
     *  limit it by upper and lower bounds (#sizeratiomax_ and #sizeratiomin_)
     *  and recompute the new time step size, if necessary. Finally, we make sure
     *  that the new time step size also satisfies upper and lower bounds (#stepsizemax_
     *  and #stepsizemin_).
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual double CalculateDt(const double norm  ///< current norm of local discretization error
    );

    /*! \brief Prepare repetition of current time step
     *
     *  Print to screen and reset certain quantities in case that the current time
     *  step has to be repeated.
     *
     *  \author mayr.mt \date 12/2013
     */
    virtual void ResetStep();

    //@}

    //! @name Access routines
    //@{

    //! get the vector of the local discretization error
    Teuchos::RCP<Epetra_Vector>& LocErrDis() { return locerrdisn_; }

    //@}

    //! @name Output
    //@{

    //! Print error norm string
    std::string PrintErrNorm() const;

    //! Print time adapting constants
    void PrintConstants(std::ostream& str  //!< output stream
        ) const;

    //! Print time adapting variables
    void PrintVariables(std::ostream& str  //!< output stream
        ) const;

    //! Print time adapting parameters:TimeIntegrator
    void Print(std::ostream& str  //!< output stream
        ) const;

    //! Modify step size to hit precisely output period
    void SizeForOutput();

    //! Prepare output to file(s)
    void PrepareOutputPeriod();

    //! Output to file(s)
    void OutputPeriod();

    //!  Update output periods
    void UpdatePeriod();

    //! Set new step size
    virtual void SetDt(const double dtnew);

    //! Update step size
    virtual void UpdateStepSize();

    //! Update step size
    virtual void UpdateStepSize(const double stpsiznew);

    //! Access to current time step size
    virtual double Dt() const { return stepsize_; }

    //! Access to target time \f$t_{n+1}\f$ of current time step
    virtual double Time() const { return time_ + stepsize_; }

    //! Check whether step size output file is attached
    bool AttachedFileStepSize()
    {
      if (not outsizefile_.is_null())
        return true;
      else
        return false;
    }

    //! Attach file handle for step size file #outsizefile_
    void AttachFileStepSize();


    //! Write step size
    void OutputStepSize();

    //@}

    //! @name Attributes
    //@{

    //! Provide the name
    virtual enum INPAR::STR::TimAdaKind MethodName() const = 0;

    //! Provide the name as std::string
    std::string MethodTitle() const { return MapKindEnumToString(MethodName()); }

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
    static int Sign(const double number  //!< a real number
    )
    {
      return (number == 0.0) ? 0 : (number > 0.0) ? +1 : -1;
    }

    //! @name General purpose algorithm members
    //@{
    Teuchos::RCP<TimInt> sti_;                       //!< marching time integrator
    Teuchos::RCP<DRT::Discretization> discret_;      //!< attached discretisation
    int myrank_;                                     //!< processor ID
    Teuchos::RCP<LINALG::Solver> solver_;            //!< linear algebraic solver
    Teuchos::RCP<IO::DiscretizationWriter> output_;  //!< binary output
    //@}

    //! @name Plain time integration constants
    //@{
    double timeinitial_;      //!< initial time: t_0
    double timefinal_;        //!< final time
    int timedirect_;          //!< +1: in positive, -1: in negative time direction
    int timestepinitial_;     //!< initial time step index: 0 (often)
    int timestepfinal_;       //!< maximum time step: n_max
    double stepsizeinitial_;  //!< initial step size: dt_n
    //@}

    //! @name Adaptive time integration constants
    //@{
    double stepsizemax_;                   //!< maximum time step size (upper limit)
    double stepsizemin_;                   //!< minimum time step size (lower limit)
    double sizeratiomax_;                  //!< maximally permitted increase of current step size
                                           //!< relative to last converged one
    double sizeratiomin_;                  //!< minimally permitted increase
                                           //!< (or maximally permitted decrease)
                                           //!< of current step size relative to last converged one
    double sizeratioscale_;                //!< safety factor, should be lower than 1.0
    enum CtrlEnum errctrl_;                //!< type of control, see #CtrlEnum
    enum INPAR::STR::VectorNorm errnorm_;  //!< norm for local error vector
    double errtol_;                        //!< target local error tolerance
    int errorder_;                         //!< order of local error indication
    int adaptstepmax_;  //!< maximally permitted trials to find tolerable step size
    //@}

    //! @name plain time integration variables
    //@{
    double time_;   //!< current time \f$t_n\f$
    int timestep_;  //!< current time step \f$n\f$
    //@}

    //! @name Adaptive time integration variables
    //@{
    double stepsizepre_;                      //!< previous time step size \f$\Delta t_{n-1}\f$
    double stepsize_;                         //!< current time step size \f$\Delta t_n\f$
    Teuchos::RCP<Epetra_Vector> locerrdisn_;  //!< current local disp. error
                                              //!< estimation \f$l_{n+1}\f$
    Teuchos::RCP<Epetra_Vector> locerrveln_;  //!< current local vel. error
                                              //!< estimation \f$\dot{l}_{n+1}\f$
    int adaptstep_;                           //!< trial counter, cf. #adaptstepmax_
    //@}

    //! @name Output settings
    //@{
    bool outsys_;                              //!< do it this step: write system to file
    bool outstr_;                              //!< do it this step: write stress/strain to file
    bool outene_;                              //!< do it this step: write energy to file
    bool outrest_;                             //!< do it this step: write restart data to file
    double outsysperiod_;                      //!< print system (dis,vel,acc,...)
                                               //!< every given period of time
    double outstrperiod_;                      //!< print stress/strain every given
                                               //!< period of time
    double outeneperiod_;                      //!< print energies every given
                                               //!< period of time
    double outrestperiod_;                     //!< print restart every given
                                               //!< period of time
    int outsizeevery_;                         //!< print step size every given step
    double outsystime_;                        //!< next output time point for system
    double outstrtime_;                        //!< next output time point for stress/strain
    double outenetime_;                        //!< next output time point for energy
    double outresttime_;                       //!< next output time point for restart
    Teuchos::RCP<std::ofstream> outsizefile_;  //!< outputfile for step sizes
    //@}

  };  // class TimAda

}  // namespace STR

/*======================================================================*/
/*!
 * \brief Out stream inserter for STR::TimAda
 *
 * \author bborn \date 10/07
 */
std::ostream& operator<<(std::ostream& str, const STR::TimAda& ta);


/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMADA_H
