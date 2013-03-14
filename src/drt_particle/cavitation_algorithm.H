/*----------------------------------------------------------------------*/
/*!
\file cavitation_algorithm.H

\brief Algorithm to control cavitation simulations

<pre>
Maintainer: Georg Hammerl
            hammerl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-152537
</pre>
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 11/12 |
 *----------------------------------------------------------------------*/
#ifndef CAVITATION_ALGORITHM_H
#define CAVITATION_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 11/12 |
 *----------------------------------------------------------------------*/
#include "particle_algorithm.H"

#include "Epetra_MpiComm.h"
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
namespace ADAPTER
{
  class Fluid;
}

/*----------------------------------------------------------------------*
 | cavitation algorithm                                     ghamm 11/12 |
 *----------------------------------------------------------------------*/
/// CAVITATION: Cavitation problem
namespace CAVITATION
{
  // forward declaration
  class BubbleSource;

  class Algorithm : public PARTICLE::Algorithm
  {
  public:

    Algorithm(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      const Teuchos::ParameterList& params   ///< problem parameters
      );

    virtual ~Algorithm(){};

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
    void ParticleInflow();

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
    virtual void BuildBubbleInflowCondition();

    /// setup ghosting of bins, particles & underlying fluid
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap,  ///< rowmap of bins
      std::map<int, std::set<int> >& fluideles   ///< map from bin id to fluid elements in this bin
      );

    /// get adjacent bins to corner, where ijk is in 1st octant
    std::vector<int> AdjacentBinstoCorner(
      int* ijk   ///< ijk is 1st octant
      );


  private:
    //! @name fluid related stuff
    //@{

    Teuchos::RCP<DRT::Discretization> fluiddis_;     ///< underlying fluid discretization

    Teuchos::RCP<ADAPTER::Fluid> fluid_;     ///< underlying fluid time integrator

    //! Info: Currently it is assumed that bins and fluid elements are fixed in space (initial content is stored here)
    //! Future: Dynamical registering of fluid elements in bins when they are moving
    std::map<int, std::set<int> > extendedfluidghosting_;   ///< map that contains all fluid element ids of each bin

    std::map<int, std::list<Teuchos::RCP<BubbleSource> > > bubble_source_;   ///< store particle sources for each bin

    //@}

  }; // Algorithm


  class BubbleSource
  {
  public:

    BubbleSource(
      int bubbleinflowid,   ///< unique id over all inflow conditions
      std::vector<double> inflow_position,  ///< position of particle source
      std::vector<double> inflow_vel,   ///< inflow velocity
      double inflow_radius,   ///< initial radius of particle
      double inflow_freq   ///< frequency of entering particles
      );

    // internal variables
    int inflowid_;
    std::vector<double> inflow_position_;
    std::vector<double> inflow_vel_;
    double inflow_radius_;
    double inflow_freq_;

  }; // BubbleSource

} // namespace CAVITATION


/*----------------------------------------------------------------------*/
#endif  // CAVITATION_ALGORITHM_H
