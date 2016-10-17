/*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.H

\brief Algorithm to control particle simulations

\level 2

\maintainer Georg Hammerl
*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | definitions                                              ghamm 09/12 |
 *----------------------------------------------------------------------*/
#ifndef PARTICLE_ALGORITHM_H
#define PARTICLE_ALGORITHM_H

/*----------------------------------------------------------------------*
 | headers                                                  ghamm 09/12 |
 *----------------------------------------------------------------------*/
#include "../drt_adapter/adapter_algorithmbase.H"
#include "binning_strategy.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_element.H"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <vector>


/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 09/12 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Discretization;
  class Node;
}
namespace ADAPTER
{
  class Structure;
  class Particle;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class MapExtractor;
}

namespace MAT
{
  namespace PAR
  {
    class ParticleMat;
    class ExtParticleMat;
  }
}

/*----------------------------------------------------------------------*
 | particle algorithm                                       ghamm 09/12 |
 *----------------------------------------------------------------------*/
/// PARTICLE: Particle simulation
namespace PARTICLE
{
  // forward declaration
  class HeatSource;

  class Algorithm : public ADAPTER::AlgorithmBase, public BINSTRATEGY::BinningStrategy
  {
  public:

    Algorithm(
      const Epetra_Comm& lcomm,  ///< local epetra communicator
      const Teuchos::ParameterList& params   ///< problem parameters
      );

    virtual ~Algorithm(){};

    /// do initialization of problem, reduced version for restart
    virtual void Init(bool restarted);

    /// set up pointers to material bundles
    virtual void InitMaterials();

    /// read restart
    virtual void ReadRestart(
      int restart  ///< time step from which to restart
      );

    /// setup of the system
    virtual void SetupSystem();

    /// outer time loop for particle problem
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(
      const Epetra_Comm& comm  ///< local epetra communicator
      );

    /// get wall discretization
    virtual Teuchos::RCP<DRT::Discretization> WallDiscret() { return particlewalldis_; }

    /// get underlyling structural problem
    virtual Teuchos::RCP<ADAPTER::Structure> Structure() { return structure_; }

    /// return gravity acceleration
    virtual LINALG::Matrix<3,1> GetGravityAcc() {return gravity_acc_ ; }

    /// particles are checked whether they have moved out of their current bin
    /// and transferred if necessary
    virtual void TransferParticles(const bool updatestates, const bool ghosting = true);

    /// get bins in range and call FindNeighbors
    virtual void GetNeighbouringParticlesAndWalls(
                    DRT::Node* CurrentParticle,
                    std::list<DRT::Node*>& neighboring_particles,
                    std::set<DRT::Element*>& neighboring_walls);

    virtual const Teuchos::RCP<ADAPTER::Particle> AdapterParticle() { return particles_; }

    /// map that links bins to heat sources
    virtual std::map<int,std::list<Teuchos::RCP<HeatSource> > > & Bins2HeatSources() {return bins2heatSources_; }

    /// list of heat sources
    virtual const std::map<int,Teuchos::RCP<HeatSource> > HeatSources() {return heatSources_; }

    /// contact type: meshfree, particle-thermo, etc.
    INPAR::PARTICLE::ParticleInteractions ParticleInteractionType() {return particleInteractionType_; }

    /// Return the pointer to the particle material
    const MAT::PAR::ParticleMat* ParticleMat() { return particleMat_; }

    /// Return the pointer to the extended particle material
    const MAT::PAR::ExtParticleMat* ExtParticleMat() { return extParticleMat_; }

  protected:

    /// prepare time step
    virtual void PrepareTimeStep();

    /// solve the current particle time step
    virtual void Integrate();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles(bool init = false);

    /// update the current time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// output particle time step
    virtual void Output(bool forced_writerestart = false);

    /// update state vectors which do not live in the time integrator to new layout
    virtual void UpdateStates(){};

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcs();

    // bins are redistributed according to the current fill level with particles
    virtual void DynamicLoadBalancing();

    /// create graph for bin distribution
    Teuchos::RCP<const Epetra_CrsGraph> CreateGraph();

    /// rough safety check for proper bin size
    virtual void BinSizeSafetyCheck(const double dt);

    /// round robin loop to fill particles into its correct bin on according proc
    virtual void FillParticlesIntoBinsRoundRobin(
      std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>& homelessparticles  ///< set of homeless particles
      );

    /// fill particles into their correct bin on according proc
    virtual void FillParticlesIntoBinsRemoteIdList(
      std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>& homelessparticles   ///< set of homeless particles
      );

    /// node is placed into the correct row bin or put into the set of homeless particles
    virtual bool PlaceNodeCorrectly(
      Teuchos::RCP<DRT::Node> node,  ///< node to be placed
      const double* currpos,  ///< current position of this node
      std::set<Teuchos::RCP<DRT::Node>, BINSTRATEGY::Less>& homelessparticles  ///< set of homeless particles
      );

    /// setup ghosting of bins and particles
    virtual void SetupGhosting(
      Teuchos::RCP<Epetra_Map> binrowmap  ///< rowmap of bins
      );

    /// particle walls are set up and assigned to bins
    virtual void SetupParticleWalls(
      Teuchos::RCP<DRT::Discretization> basediscret
      );

    /// build connectivity from particle wall elements to bins
    virtual void BuildElementToBinPointers(
      bool wallpointer   ///< bool whether wall element pointer need rebuild
      );

    /// bins are distributed to the processors
    virtual Teuchos::RCP<Epetra_Map> DistributeBinsToProcsBasedOnUnderlyingDiscret(
      Teuchos::RCP<DRT::Discretization> underlyingdis,  ///< underlying discretization on which distribution depends on
      std::map<int, std::set<int> >& fluideles  ///< map from bin id to row fluid elements in this bin
      );

    /// assign wall elements to bins
    virtual void AssignWallElesToBins();

    /// set up heat sources
    virtual void SetUpHeatSources();

    /// update of the map bins->heat Sources
    virtual void UpdateHeatSourcesConnectivity(bool forceRestart);

    /// update connectivity: mask to sort various connectivity updates
    virtual void UpdateConnectivity();

    /// dismember particle using the dismembering radius (it shall be smaller than the oldParticleRadius)
    virtual void ParticleDismemberer();

    /// mass and density updater for particle dismemberer
    void MassDensityUpdaterForParticleDismemberer(
        Teuchos::RCP<Epetra_Vector> &mass,
        Teuchos::RCP<Epetra_Vector> &density,
        Teuchos::RCP<Epetra_Vector> &radius,
        const int &lidNode_new,
        const int &lidNode_old,
        const int &nlist);

    /// small function for the dismember
    virtual int ComputeSemiLengthInParticlesForParticleDismemberer(const double &oldRadius,const double &semiStep);

    /// get particles and wall elements in given bins
    void GetBinContent(
      std::list<DRT::Node*>& particles,
      std::set<DRT::Element*>& walls,
      std::vector<int>& binIds
      );

    //! @name particle related stuff
    //@{

    Teuchos::RCP<ADAPTER::Particle> particles_;     ///< particle time integration

    Teuchos::RCP<Epetra_Map> bincolmap_;     ///< colmap of bins

    LINALG::Matrix<3,1> gravity_acc_;   ///< acceleration due to gravity

    int writeresultsevery_;   ///< output interval

    Teuchos::RCP<ADAPTER::Structure> structure_;     ///< structural (wall) time integration

    Teuchos::RCP<DRT::Discretization> particlewalldis_;     ///< fully redundant particle wall discretization

    Teuchos::RCP<LINALG::MapExtractor> wallextractor_;   ///< map extractor for particle wall

    const bool moving_walls_;   ///< flag whether walls are moving


    std::map<int,Teuchos::RCP<HeatSource> > heatSources_; ///< heat source storing

    std::map<int,std::list<Teuchos::RCP<HeatSource> > > bins2heatSources_;   ///< heat sources 2 bins map

    const INPAR::PARTICLE::ParticleInteractions particleInteractionType_;

    const MAT::PAR::ParticleMat* particleMat_;  ///< pointer to the extParticleMat material

    const MAT::PAR::ExtParticleMat* extParticleMat_;  ///< pointer to the extParticleMat material

    //@}

  }; // Algorithm

  struct HeatSource
  {
  public:

  HeatSource(
    bool active,                            ///< flag that simplifies the connectivity algorithm
    const int id,                           ///< unique id over all heat sources conditions
    const std::vector<double> minVerZone,   ///< heat zone, vertex 0
    const std::vector<double> maxVerZone,   ///< heat zone, vertex 1
    const double QDot,                      ///< intake q/s
    const double Tstart,                    ///< t start
    const double Tend) :                    ///< t end
       active_(active),
       id_(id),
       minVerZone_(minVerZone),
       maxVerZone_(maxVerZone),
       QDot_(QDot),
       Tstart_(Tstart),
       Tend_(Tend) {}

  // internal variables
  bool active_;                             ///< flag that simplifies the connectivity algorithm
  const int id_;                            ///< unique id over all heat sources conditions
  const std::vector<double> minVerZone_;    ///< heat zone, vertex 0
  const std::vector<double> maxVerZone_;    ///< heat zone, vertex 1
  const double QDot_;                       ///< intake heat per second per unit volume ( q/(s * V) )
  const double Tstart_;                     ///< t start
  const double Tend_;                       ///< t end
  }; // HeatSource

  struct homelessParticleTemp
  {
  public:
  int lidNode_old;
  int lidDof_old;
  std::vector<double> pos;
  };

} // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // PARTICLE_ALGORITHM_H
