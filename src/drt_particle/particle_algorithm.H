 /*----------------------------------------------------------------------*/
/*!
\file particle_algorithm.H

\brief Algorithm to control particle simulations

\level 2

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/

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
#include "../drt_binstrategy/binning_strategy_utils.H"

#include <Epetra_MpiComm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_CrsGraph.h>
#include <vector>
#include <boost/unordered_map.hpp>

//Usefull defines (tolerances and debugging tools)
#define PARTICLE_TENSILESAFETYFAC 0.0000001         //Factor prescribing the minimal particle distance as compared to the average initial particle distance (tensile instabilites)
//#define PARTICLE_OVERLAPPINGNEIGHBORS             //Determine also fully overlapping variant of neighbors vectors
//#define PARTICLE_ONLYLAPLACETERM                  //Apply only Laplace term in viscous force calculation according to variant 1 (this is sufficient for incompressible flows)
//#define PARTICLE_MINMAXOUTPUT_EVRY 100            //Ouput of minimal/maximal pressures, particle distances etc. every PARTICLE_MINMAXOUTPUT_EVRY time steps
//#define PARTICLE_VELDIFFBOUNDARY                  //Add boundary particle contributions to the velocity difference term (tensor A)
#define PARTICLE_ARTVISCBOUNDARY                    //Apply artificial viscosity also for boundary particle interaction (in contrast to Adami et al. 2012)
#define PARTICLE_COLORLIMIT 0.9                     //Color field limit value to identify free-surface particle
//#define PARTICLE_NORESTARTACC                       //Calculate initial acceleration new in case of a restart (don't use 'restart' value)


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
    virtual Teuchos::RCP<DRT::Discretization> WallDiscret() const { return particlewalldis_; }

    /// get wall extractor
    virtual Teuchos::RCP<LINALG::MapExtractor> GetWallExtractor() { return wallextractor_; }

    /// set wall states
    virtual void SetWallStates(
        Teuchos::RCP<const Epetra_Vector> walldispn,   ///< wall displacements at \f$t_{n}\f$
        Teuchos::RCP<const Epetra_Vector> walldispnp,  ///< wall displacements at \f$t_{n+1}\f$
        Teuchos::RCP<const Epetra_Vector> wallvelnp    ///< wall velocities at \f$t_{n+1}\f$
        );

    /// set up walls
    virtual void SetUpWallDiscret();

    /// return gravity acceleration
    virtual LINALG::Matrix<3,1> GetGravityAcc(const double time=-1.0);

    /// particles are checked whether they have moved out of their current bin
    /// and transferred if necessary, update states and do ghosting if needed
    virtual Teuchos::RCP<std::list<int> > TransferParticles(const bool updatestates,
                                   const bool ghosting = true);

    /// update of the map bins->heat Sources
    virtual void UpdateHeatSourcesConnectivity(bool forceRestart);

    /// get bins in range and call FindNeighbors
    virtual void GetNeighbouringItems(
                    DRT::Node* particle,
                    std::list<DRT::Node*>& neighboursLinf_p,
                    boost::unordered_map<int, DRT::Element*>& neighboursLinf_w,
                    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs = Teuchos::null) const;

    /// get bins in range and call FindNeighbors (bin version)
    virtual void GetNeighbouringItems(
                    const int binId,
                    std::list<DRT::Node*>& neighboursLinf_p,
                    boost::unordered_map<int, DRT::Element*>* neighboursLinf_w = NULL,
                    const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > neighboursLinf_hs = Teuchos::null) const;

    /// link to timInt
    virtual const Teuchos::RCP<ADAPTER::Particle> AdapterParticle() { return particles_; }

    /// map that links bins to heat sources
    virtual std::map<int,std::list<Teuchos::RCP<HeatSource> > > & Bins2HeatSources() {return bins2heatSources_; }

    /// list of heat sources
    virtual const std::list< Teuchos::RCP<HeatSource> > HeatSources() {return heatSources_; }

    /// contact type: meshfree, particle-thermo, etc.
    INPAR::PARTICLE::ParticleInteractions ParticleInteractionType() const {return particleInteractionType_; }

    /// Return the pointer to the particle material
    const MAT::PAR::ParticleMat* ParticleMat() const { return particleMat_; }

    /// Return the pointer to the extended particle material
    const MAT::PAR::ExtParticleMat* ExtParticleMat() const { return extParticleMat_; }

    /// get meshfree rendering handler
    virtual Teuchos::RCP<Rendering> GetRendering() { return rendering_; }

    /// set meshfree rendering handler
    virtual void SetRendering(Teuchos::RCP<Rendering> rendering) { rendering_ = rendering; }

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

    /// build extended bin column map
    virtual Teuchos::RCP<Epetra_Map> ExtendedBinColMap();

    /// ghosting of an additional second bin layer
    virtual void AddLayerGhosting(
        std::set<int>& colbins
        );

    /// ghosting of bins in proximity of ghosted boundary particles
    virtual void BdryParticleGhosting(
        std::set<int>& colbins
        );

    /// ghosting of bins in proximity of ghosted wall elements
    virtual void WallElementGhosting(
        std::set<int>& colbins
        );

    /// build connectivity from particle wall elements to bins
    virtual void BuildElementToBinPointers(
      bool wallpointer   ///< bool whether wall element pointer need rebuild
      );

    /// assign wall elements and gids to bins
    virtual void AssignWallElesAndGidsToBins();

    /// rough safety check for proper bin size
    virtual void BinSizeSafetyCheck(const double dt);

    /// access structure and setup particle wall
    virtual void AccessStructure();

    /// particle walls are set up and assigned to bins
    virtual void SetupParticleWalls(
      Teuchos::RCP<DRT::Discretization> basediscret,
      const std::string elename
      );

    /// relate wall gids to bin ids
    virtual void RelateWallGidsToBinIds();

    /// set up heat sources
    virtual void SetUpHeatSources();

    /// get particles and wall elements in given bins
    void GetBinContent(
      std::list<DRT::Node*>& bin_p,
      boost::unordered_map<int, DRT::Element*>* bin_w,
      const Teuchos::RCP<boost::unordered_map<int , Teuchos::RCP<HeatSource> > > bin_hs,
      std::vector<int>& binIds) const;

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

    Teuchos::RCP<const Epetra_Vector> walldispn_; ///< wall displacements at \f$t_{n}\f$

    Teuchos::RCP<const Epetra_Vector> walldispnp_; ///< wall displacements at \f$t_{n+1}\f$

    Teuchos::RCP<const Epetra_Vector> wallvelnp_; ///< wall velocities at \f$t_{n+1}\f$

    const bool moving_walls_;   ///< flag whether walls are moving

    const int transfer_every_;   ///< interval to call particle transfer

    std::list<Teuchos::RCP<HeatSource> > heatSources_; ///< heat source storing

    std::map<int,std::list<Teuchos::RCP<HeatSource> > > bins2heatSources_;   ///< bins 2 heat sources map

    const INPAR::PARTICLE::ParticleInteractions particleInteractionType_;

    const INPAR::PARTICLE::ExtendedGhosting extendedGhosting_;

    const MAT::PAR::ParticleMat* particleMat_;  ///< pointer to the extParticleMat material

    const MAT::PAR::ExtParticleMat* extParticleMat_;  ///< pointer to the extParticleMat material

    BINSTRATEGY::UTILS::BinContentType bin_wallcontent_;

    Teuchos::RCP<Rendering> rendering_; ///< handler of the rendering discretization (used for meshFree dynamics)

    std::map<int, std::set<int> > relwallgidtobinids_; ///< map relating wall gids to bin ids

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
