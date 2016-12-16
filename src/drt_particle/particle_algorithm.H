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
#include "particle_handler.H"
#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_inpar/inpar_binstrategy.H"

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
  class Element;
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
namespace BINSTRATEGY
{
  class BinningStrategy;
}

/*----------------------------------------------------------------------*
 | particle algorithm                                       ghamm 09/12 |
 *----------------------------------------------------------------------*/
/// PARTICLE: Particle simulation
namespace PARTICLE
{
  // forward declaration
  class HeatSource;
  class Rendering;

  class Algorithm : public ADAPTER::AlgorithmBase, public PARTICLE::ParticleHandler
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
    /// and transferred if necessary, update states and do ghosting if needed
    virtual Teuchos::RCP<std::list<int> > TransferParticles(const bool updatestates,
                                   const bool ghosting = true);

    /// get bins in range and call FindNeighbors
    virtual void GetNeighbouringItems(
                    DRT::Node* particle,
                    std::list<DRT::Node*>& neighboring_particles,
                    std::set<DRT::Element*>& neighboring_walls,
                    const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources = Teuchos::null);

    /// get bins in range and call FindNeighbors (bin version)
    virtual void GetNeighbouringItems(
                    const int binId,
                    std::list<DRT::Node*>& neighboring_particles,
                    std::set<DRT::Element*>& neighboring_walls,
                    const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources = Teuchos::null);

    /// link to timInt
    virtual const Teuchos::RCP<ADAPTER::Particle> AdapterParticle() { return particles_; }

    /// map that links bins to heat sources
    virtual std::map<int,std::list<Teuchos::RCP<HeatSource> > > & Bins2HeatSources() {return bins2heatSources_; }

    /// list of heat sources
    virtual const std::list< Teuchos::RCP<HeatSource> > HeatSources() {return heatSources_; }

    /// contact type: meshfree, particle-thermo, etc.
    INPAR::PARTICLE::ParticleInteractions ParticleInteractionType() {return particleInteractionType_; }

    /// Return the pointer to the particle material
    const MAT::PAR::ParticleMat* ParticleMat() { return particleMat_; }

    /// Return the pointer to the extended particle material
    const MAT::PAR::ExtParticleMat* ExtParticleMat() { return extParticleMat_; }

    virtual Teuchos::RCP<Rendering> GetRendering() { return rendering_; }

  protected:

    /// prepare time step
    virtual void PrepareTimeStep();

    /// solve the current particle time step
    virtual void Integrate();

    /// force on particles is calculated and applied
    virtual void CalculateAndApplyForcesToParticles(bool init = false);

    /// accelerations on particleMeshFree are calculated and applied
    virtual void CalculateAndApplyAccelerationsToParticles(bool init = false);

    /// update the current time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// output particle time step
    virtual void Output(bool forced_writerestart = false);

    /// update state vectors which do not live in the time integrator to new layout
    virtual void UpdateStates(){};

    // bins are redistributed according to the current fill level with particles
    virtual void DynamicLoadBalancing();

    /// build connectivity from particle wall elements to bins
    virtual void BuildElementToBinPointers(
      bool wallpointer   ///< bool whether wall element pointer need rebuild
      );

    /// rough safety check for proper bin size
    virtual void BinSizeSafetyCheck(const double dt);

    /// particle walls are set up and assigned to bins
    virtual void SetupParticleWalls(
      Teuchos::RCP<DRT::Discretization> basediscret
      );

    /// assign wall elements to bins
    virtual void AssignWallElesToBins();

    /// set up heat sources
    virtual void SetUpHeatSources();

    /// update of the map bins->heat Sources
    virtual void UpdateHeatSourcesConnectivity(bool forceRestart);

    /// update connectivity: mask to sort various connectivity updates
    virtual void UpdateConnectivity();

    /// get particles and wall elements in given bins
    void GetBinContent(
      std::list<DRT::Node*>& particles,
      std::set<DRT::Element*>& walls,
      const Teuchos::RCP<std::set<Teuchos::RCP<HeatSource>, BINSTRATEGY::Less> > neighboring_heatSources,
      std::vector<int>& binIds);

    /// dismember particle using the dismembering radius (it shall be smaller than the oldParticleRadius)
    virtual void ParticleDismemberer();

    /// mass and density updater for particle dismemberer
    virtual void MassDensityUpdaterForParticleDismemberer(
        Teuchos::RCP<Epetra_Vector> &mass,
        Teuchos::RCP<Epetra_Vector> &density,
        Teuchos::RCP<Epetra_Vector> &radius,
        const int &lidNode_new,
        const int &lidNode_old,
        const int &nlist);

    /// small function for the dismember
    virtual int ComputeSemiLengthInParticlesForParticleDismemberer(const double &oldRadius,const double &semiStep)
      {return (oldRadius-semiStep)/(2 * semiStep) + 2; }

    /// compute thermodynamic expansion - new densities and radii
    virtual void ThermalExpansion();

    //! @name particle related stuff
    //@{

    Teuchos::RCP<ADAPTER::Particle> particles_;     ///< particle time integration

    LINALG::Matrix<3,1> gravity_acc_;   ///< acceleration due to gravity

    int writeresultsevery_;   ///< output interval

    Teuchos::RCP<ADAPTER::Structure> structure_;     ///< structural (wall) time integration

    Teuchos::RCP<DRT::Discretization> particlewalldis_;     ///< fully redundant particle wall discretization

    Teuchos::RCP<Epetra_Map> particlewallelecolmap_standardghosting_;     ///< wall ele col map containing standard ghosting

    Teuchos::RCP<LINALG::MapExtractor> wallextractor_;   ///< map extractor for particle wall

    const bool moving_walls_;   ///< flag whether walls are moving

    const int transfer_every_;   ///< interval to call particle transfer

    std::list<Teuchos::RCP<HeatSource> > heatSources_; ///< heat source storing

    std::map<int,std::list<Teuchos::RCP<HeatSource> > > bins2heatSources_;   ///< bins 2 heat sources map

    const INPAR::PARTICLE::ParticleInteractions particleInteractionType_;

    const MAT::PAR::ParticleMat* particleMat_;  ///< pointer to the extParticleMat material

    const MAT::PAR::ExtParticleMat* extParticleMat_;  ///< pointer to the extParticleMat material

    INPAR::BINSTRATEGY::BinContent bin_surfcontent_;

    Teuchos::RCP<Rendering> rendering_; /// handler of the rendering discretization (used for meshFree dynamics)

    //@}

  }; // Algorithm

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
