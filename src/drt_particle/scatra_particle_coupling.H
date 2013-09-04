/*----------------------------------------------------------------------*/
/*!
\file scatra_particle_coupling.H

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 11/12 |
 *----------------------------------------------------------------------*/
#ifndef SCATRA_PARTICLE_COUPLING_H
#define SCATRA_PARTICLE_COUPLING_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "cavitation_algorithm.H"
#include "../drt_inpar/inpar_cavitation.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
namespace SCATRA
{
  class ScaTraTimIntImpl;
}
namespace DRT
{
  class Discretization;
  class Element;
}

/*----------------------------------------------------------------------*
 | cavitation algorithm                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
/// CAVITATION: Cavitation problem
namespace PARTICLE
{
  // forward declaration
//  class BubbleSource;

  class ScatraParticleCoupling : public PARTICLE::Algorithm
  {
  public:

    ScatraParticleCoupling(
      Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra,   ///< underlying scatra field
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      const Teuchos::ParameterList& params   ///< problem parameters
      );

    /// do initialization of problem
    virtual void Init(bool restarted);

    /// read restart
    virtual void ReadRestart(
      int restart  ///< time step from which to restart
      );

    /// setup of the system
    virtual void SetupSystem();

    /// outer time loop for cavitation problem
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(
      const Epetra_Comm& comm  ///< local epetra communicator
      );


  protected:

    /// prepare time step
    virtual void PrepareTimeStep();

    /// solve the current particle time step
    virtual void Integrate();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles();

    /// insert particles into domain
//    void ParticleInflow();

    /// update the current time step
    virtual void Update();

    /// output particle time step
    virtual void Output();

    /// find XAABB and divide into bins using the cutoff radius
    virtual void CreateBins();

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcs(
      std::map<int, std::set<int> >& fluideles  ///< map from bin id to fluid elements in this bin
      );

    /// setup of bubble sources in each bin
//    virtual void BuildBubbleInflowCondition();

    /// setup ghosting of bins, particles & underlying fluid
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap,  ///< rowmap of bins
      std::map<int, std::set<int> >& fluideles   ///< map from bin id to fluid elements in this bin
      );

     //@}

  private:

    //! @name scatra related stuff
    //@{

    Teuchos::RCP<SCATRA::ScaTraTimIntImpl> scatra_;     ///< underlying fluid time integrator

    Teuchos::RCP<DRT::Discretization> scatradis_;     ///< underlying fluid discretization

//    std::map<int, std::list<Teuchos::RCP<BubbleSource> > > bubble_source_;   ///< store particle sources for each bin

    //@}

  }; // Algorithm


//  class BubbleSource
//  {
//  public:
//
//    BubbleSource(
//      int bubbleinflowid,   ///< unique id over all inflow conditions
//      std::vector<double> inflow_position,  ///< position of particle source
//      std::vector<double> inflow_vel,   ///< inflow velocity
//      int inflow_vel_curve,   ///< inflow velocity can be superposed with a time curve
//      double inflow_radius,   ///< initial radius of particle
//      double inflow_freq   ///< frequency of entering particles
//      );
//
//    // internal variables
//    int inflowid_;
//    std::vector<double> inflow_position_;
//    std::vector<double> inflow_vel_;
//    int inflow_vel_curve_;
//    double inflow_radius_;
//    double inflow_freq_;
//
//  }; // BubbleSource

} // namespace CAVITATION


/*----------------------------------------------------------------------*/
#endif  // CAVITATION_ALGORITHM_H
