/*----------------------------------------------------------------------*/
/*!
\file acou_timeint.H

\brief Acoustics Base Algorithm for implicit and explicit time integration

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*----------------------------------------------------------------------*/

#ifndef ACOU_TIMINT_H
#define ACOU_TIMINT_H

/*----------------------------------------------------------------------*
 | headers                                               schoeder 03/15 |
 *----------------------------------------------------------------------*/

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Epetra_Operator.h>
#include <Epetra_CrsGraph.h>

#include "../drt_inpar/inpar_acou.H"

/*----------------------------------------------------------------------*
 | forward declarations                                  schoeder 03/15 |
 *----------------------------------------------------------------------*/
namespace LINALG
{
  class SparseOperator;
  class MapExtractor;
  class Solver;
}  // namespace LINALG
namespace DRT
{
  class Discretization;
  class DiscretizationHDG;
  class Node;
  class ResultTest;
}  // namespace DRT
namespace IO
{
  class DiscretizationWriter;
}

namespace SCATRA
{
  class TimIntStationary;
}
/*----------------------------------------------------------------------*
 | general time intergrator for acoustics                schoeder 03/15 |
 *----------------------------------------------------------------------*/
namespace ACOU
{
  class PATMonitorManager;

  class AcouTimeInt
  {
   public:
    // constructor
    AcouTimeInt(const Teuchos::RCP<DRT::DiscretizationHDG>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output);

    /// virtual destructor
    virtual ~AcouTimeInt();

    /// write restart
    virtual void WriteRestart();

    /// read restart information for given time step
    virtual void ReadRestart(int step);

    /// set initial field to zero (needed for adjoint problem)
    virtual void SetInitialZeroField();

    /// set initial field by given function
    virtual void SetInitialField(int startfuncno) = 0;

    /// set initial field by given light distribution
    // virtual void SetInitialPhotoAcousticField(Teuchos::RCP<SCATRA::TimIntStationary> scatraalgo);
    virtual void SetInitialPhotoAcousticField(Teuchos::RCP<Epetra_Vector> light,
        Teuchos::RCP<DRT::Discretization> scatradis, bool meshconform);

    /// Output nodal internal pressure field at new time step
    virtual void NodalPressureField(Teuchos::RCP<Epetra_Vector> outvec) = 0;

    /// output routine: not inherited since it does something completely different than the one for
    /// the implicit integration
    virtual void Output();

    // output of material parameters for inverse analysis
    void OutputDensityAndSpeedOfSound();

    /// init a monitor file
    void InitMonitorFile();

    /// fill monitor file
    void FillMonitorFile(Teuchos::RCP<Epetra_Vector> ip);

    /// print out user information
    void PrintInformationToScreen();

    /// run a complete simulation
    virtual void Integrate() = 0;

    /// increment time and step value
    void IncrementTimeAndStep()
    {
      step_ += 1;
      time_ += dtp_;
    }

    double Time() { return time_; }
    int Step() { return step_; }
    int UpRes() { return upres_; }
    double TimeStep() { return dtp_; }

    /// create result test for encapsulated structure algorithm
    virtual Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    /// returns pointer to the discretization
    Teuchos::RCP<DRT::Discretization> Discretization();

    /// Output nodal internal psi field at new time step
    virtual void NodalPsiField(Teuchos::RCP<Epetra_Vector> outvec);

    /// return the name of the time integration scheme
    virtual std::string Name() = 0;

    /// return the adjoint source vector
    const Teuchos::RCP<PATMonitorManager> MonitorManager() { return monitor_manager_; }

    bool iswithPMLS() { return withpmls_; }

    bool isWriteMonitor() { return writemonitor_; }

    // get values for speed of sound gradient
    virtual double GetSoSGradient(int colid) = 0;

    // get values for density gradient
    virtual double GetDensityGradient(int colid) = 0;

   protected:
    /// discretization, solver, parameter list and output
    Teuchos::RCP<DRT::DiscretizationHDG> discret_;
    Teuchos::RCP<LINALG::Solver> solver_;
    Teuchos::RCP<Teuchos::ParameterList> params_;
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    INPAR::ACOU::DynamicType dyna_;   /// time integration scheme
    INPAR::ACOU::PhysicalType phys_;  /// physical type

    bool invana_;  /// flag indicating whether inverse analysis is performed
    bool reduction_;
    bool padaptivity_;  /// p-adaptivity flag
    bool adjoint_;      /// flag indicating whether this is standard or an adjoint solve
    bool withpmls_;
    int myrank_;         /// processor id
    bool writemonitor_;  /// flag indicating if a monitor file is written
    bool writestress_;   /// flag indicating if a stresses or -if no- the displacment gradient is
                         /// written to the output

    double time_;  /// physical time
    int step_;     /// time step
    int restart_;  /// restart step

    double maxtime_;  /// maximum time
    int stepmax_;     /// maximum step

    int uprestart_;  /// write restart data every uprestart_ steps
    int upres_;      /// write output every upres_ steps

    int numdim_;  /// number of spatial dimensions
    double dtp_;  /// time step size

    //! @name trace at time n+1
    Teuchos::RCP<Epetra_Vector> velnp_;

    //! monitor manager for inverse analysis
    Teuchos::RCP<PATMonitorManager> monitor_manager_;
  };

}  // namespace ACOU

#endif /* ACOU_TIMINT_H */
