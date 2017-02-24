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
#include "../drt_particle/binning_strategy_utils.H"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <vector>
#include <boost/unordered_map.hpp>

//TODO: Replace these defines by input parameters
#define PARTICLE_DIM 3                        //space dimension
//#define PARTICLE_CONSISTENT_VOLUME 0.2      //Value of problem/domain volume in order to apply SPH-consistent calculation of initial particle mass
#define PARTICLE_DELTADENSFAC 1.0            //Scale factor for reference density in equation of state

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

    /// set structural states
    virtual void SetStructStates(
        Teuchos::RCP<const Epetra_Vector> structdispn,   ///< structural displacements at t_n
        Teuchos::RCP<const Epetra_Vector> structdispnp,  ///< structural displacements at t_n+1
        Teuchos::RCP<const Epetra_Vector> structvelnp    ///< structural velocities at t_n
        );

    /// return gravity acceleration
    virtual LINALG::Matrix<3,1> GetGravityAcc() {return gravity_acc_ ; }

    /// particles are checked whether they have moved out of their current bin
    /// and transferred if necessary, update states and do ghosting if needed
    virtual Teuchos::RCP<std::list<int> > TransferParticles(const bool updatestates,
                                   const bool ghosting = true);

    /// get bins in range and call FindNeighbors
    virtual void GetNeighbouringItems(
                    DRT::Node* particle,
                    std::list<DRT::Node*>& neighboursLinf_p,
                    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
                    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs = Teuchos::null);

    /// get bins in range and call FindNeighbors (bin version)
    virtual void GetNeighbouringItems(
                    const int binId,
                    std::list<DRT::Node*>& neighboursLinf_p,
                    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
                    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs = Teuchos::null);

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

    /// set up walls
    virtual void SetUpWallDiscret();

    inline virtual void CalculateAndApplyForcesToParticles(bool init) { ; };

    /// bins are redistributed according to the current fill level with particles
    virtual void DynamicLoadBalancing();

    /// update connectivity: mask to sort various connectivity updates
    virtual void UpdateConnectivity();

    /// prepare time step
    virtual void PrepareTimeStep(bool print_header = true);

    /// output particle time step
    virtual void Output(bool forced_writerestart = false);

    /// adaptions for Normal_DEM_thermo
    virtual void NormDemThermoAdapt();

  protected:

    /// solve the current particle time step
    virtual void Integrate();

    /// update the current time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// update state vectors which do not live in the time integrator to new layout
    virtual void UpdateStates(){};

    /// build connectivity from particle wall elements to bins
    virtual void BuildElementToBinPointers(
      bool wallpointer   ///< bool whether wall element pointer need rebuild
      );

    /// rough safety check for proper bin size
    virtual void BinSizeSafetyCheck(const double dt);

    /// access structure and setup particle wall
    virtual void AccessStructure();

    /// particle walls are set up and assigned to bins
    virtual void SetupParticleWalls(
      Teuchos::RCP<DRT::Discretization> basediscret,
      const std::string elename
      );

    /// assign wall elements to bins
    virtual void AssignWallElesToBins();

    /// set up heat sources
    virtual void SetUpHeatSources();

    /// update of the map bins->heat Sources
    virtual void UpdateHeatSourcesConnectivity(bool forceRestart);

    /// get particles and wall elements in given bins
    void GetBinContent(
      std::list<DRT::Node*>& bin_p,
      boost::unordered_map<int, DRT::Element*>& bin_w,
      const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > bin_hs,
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

    Teuchos::RCP<DRT::Discretization> particlewalldis_;     ///< fully redundant particle wall discretization

    Teuchos::RCP<Epetra_Map> particlewallelecolmap_standardghosting_;     ///< wall ele col map containing standard ghosting

    Teuchos::RCP<LINALG::MapExtractor> wallextractor_;   ///< map extractor for particle wall

    Teuchos::RCP<const Epetra_Vector> structdispn_; ///< structural displacements at t_n

    Teuchos::RCP<const Epetra_Vector> structdispnp_; ///< structural displacements at t_n+1

    Teuchos::RCP<const Epetra_Vector> structvelnp_; ///< structural velocities at t_n

    const bool moving_walls_;   ///< flag whether walls are moving

    const int transfer_every_;   ///< interval to call particle transfer

    std::list<Teuchos::RCP<HeatSource> > heatSources_; ///< heat source storing

    std::map<int,std::list<Teuchos::RCP<HeatSource> > > bins2heatSources_;   ///< bins 2 heat sources map

    const INPAR::PARTICLE::ParticleInteractions particleInteractionType_;

    const MAT::PAR::ParticleMat* particleMat_;  ///< pointer to the extParticleMat material

    const MAT::PAR::ExtParticleMat* extParticleMat_;  ///< pointer to the extParticleMat material

    BINSTRATEGY::UTILS::BinContentType bin_wallcontent_;

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
