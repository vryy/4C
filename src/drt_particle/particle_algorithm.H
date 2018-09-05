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
#include <boost/unordered_map.hpp>

//*********************************SPH and
// DEM************************************************************************************************************************************
//#define PARTICLE_MINMAXOUTPUT_EVRY 25            //Ouput of minimal/maximal pressures, particle
// distances etc. every PARTICLE_MINMAXOUTPUT_EVRY time steps

//*********************************SPH********************************************************************************************************************************************
// Usefull defines (tolerances and constants)
#define PARTICLE_TENSILESAFETYFAC \
  0.0000001  // Factor prescribing the minimal particle distance as compared to the average initial
             // particle distance (tensile instabilites)
#define PARTICLE_COLORLIMIT 0.9  // Color field limit value to identify free-surface particle
#define PARTICLE_COLORGRADLIMIT \
  0.1  // Color gradient limit for the evaluation of surface tension forces
#define PARTICLE_SFEXP \
  5  // Additional weighting of surface tension forces (continuum model according to Adami et al.
     // 2010) for free-surface flow
#define PARTICLE_P0ZHANG \
  3.0  // Factor between h and \tilde{h} of smoothing kernel applied in background pressure term for
       // free-surface flow according to Zhang et al. 2017

#define PARTICLE_REINITSHIFT \
  2  // Define shifting criterion between free-surface and interior particles in case of density
     // reinitialization (default: 2) 1: Apply shifting criterion for free-surface particles
     // according to free-surface particle type (fspType_) 2: Apply shifting criterion for
     // free-surface particles according to colorfield<1 (see Zhang et al. 2017) 3: Apply
     // reinitialization everywhere

#define PARTICLE_ST_TRIPPLEPOINTSMOOTHFAC \
  2.0  // Smoothing factor for prescription of normal vector via static contact angle
#define PARTICLE_ST_BOUNDARYSCALEFAC \
  1.0  // Increase the surface tension forces in the triple zone by scaling
       // factor 1.0+PARTICLE_ST_BOUNDARYSCALEFAC
#define PARTICLE_ST_CFGEXTRAPLATION  // Extrapolate direction of CFG vector
//#define PARTICLE_ONLYPHASE1ACC                      //Apply surface tension forces only to the
// phase with higher density

// Usefull defines (additional functionalities and debugging tools)
//#define PARTICLE_ONLYLAPLACETERM                  //Apply only Laplace term in viscous force
// calculation according to variant 1 (this is sufficient for incompressible flows) #define
// PARTICLE_BOUNDARYDENSITY                  //Determine the density of boundary particles required
// for colorfield and density reinitialization on the basis of pressure extrapolation (similar to
// Adami et al. 2012) #define PARTICLE_WRITECOLORFIELD                  //Output of colorField_ and
// fspType_ (free-surface particle type) for paraview

//*********************************DEM********************************************************************************************************************************************
#define PARTICLE_BINNINGCONST 5  // Constant required for binning strategy
//#define PARTICLE_NOTIMESTEPDSERROR                  //Only throw warning (no dserror), if yield
// strength criterion is violated!

/*----------------------------------------------------------------------*
 | forward declarations                                     ghamm 09/12 |
 *----------------------------------------------------------------------*/
namespace DRT
{
  class Element;
  class Discretization;
  class Node;
}  // namespace DRT
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
  }  // namespace PAR
}  // namespace MAT
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
  class Rendering;

  class Algorithm : public ADAPTER::AlgorithmBase, public PARTICLE::ParticleHandler
  {
   public:
    Algorithm(const Epetra_Comm& lcomm,       ///< local epetra communicator
        const Teuchos::ParameterList& params  ///< problem parameters
    );

    virtual ~Algorithm(){};

    /// do initialization of problem, reduced version for restart
    virtual void Init(bool restarted);

    /// set up pointers to material bundles
    virtual void InitMaterials();

    /// read restart
    virtual void ReadRestart(int restart  ///< time step from which to restart
    );

    /// setup of the system
    virtual void SetupSystem();

    /// outer time loop for particle problem
    virtual void Timeloop();

    /// test results (if necessary)
    virtual void TestResults(const Epetra_Comm& comm  ///< local epetra communicator
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
    virtual LINALG::Matrix<3, 1> GetGravityAcc(const double time = -1.0);

    /// bins are redistributed according to the current fill level with particles
    virtual void DynamicLoadBalancing();

    /// update connectivity: mask to sort various connectivity updates
    virtual void UpdateConnectivity();

    /// set position of particle nodes
    virtual void SetParticleNodePos();

    /// movement of particles since last redistribution is checked
    /// in order to determine if a new redistribution of particles is needed
    virtual bool CheckAdaptiveRepartition();

    /// particles are checked whether they have moved out of their current bin
    /// and transferred if necessary, update states and do ghosting if needed
    virtual Teuchos::RCP<std::list<int>> TransferParticles(
        const bool updatestates, const bool ghosting = true);

    /// return importer for dof based vectors
    virtual const Teuchos::RCP<Epetra_Import> DofImporter() const { return dofimporter_; }

    /// return importer for node based vectors
    virtual const Teuchos::RCP<Epetra_Import> NodeImporter() const { return nodeimporter_; }

    /// get bins in range and call FindNeighbors
    virtual void GetNeighbouringItems(DRT::Node* particle, std::list<DRT::Node*>& neighboursLinf_p,
        boost::unordered_map<int, DRT::Element*>& neighboursLinf_w) const;

    /// get bins in range and call FindNeighbors (bin version)
    virtual void GetNeighbouringItems(const int binId, std::list<DRT::Node*>& neighboursLinf_p,
        boost::unordered_map<int, DRT::Element*>* neighboursLinf_w = NULL) const;

    /// link to time integration
    virtual const Teuchos::RCP<ADAPTER::Particle> AdapterParticle() const { return particles_; }

    /// get particle interaction type
    INPAR::PARTICLE::ParticleInteractions ParticleInteractionType() const
    {
      return particleInteractionType_;
    }

    /// return vector of pointers to particle material parameters
    std::vector<const MAT::PAR::ParticleMat*> ParticleMat() const { return particleMat_; }

    /// get SPH rendering handler
    virtual Teuchos::RCP<Rendering> GetRendering() { return rendering_; }

    /// set SPH rendering handler
    virtual void SetRendering(Teuchos::RCP<Rendering> rendering) { rendering_ = rendering; }

    inline virtual void CalculateAndApplyForcesToParticles(bool init) { ; };

    /// prepare time step
    virtual void PrepareTimeStep(bool print_header = true);

    /// output particle time step
    virtual void Output(bool forced_writerestart = false);

   protected:
    /// solve the current particle time step
    virtual void Integrate();

    /// update the current time step
    virtual void Update();

    /// calculate stresses, strains, energies
    virtual void PrepareOutput();

    /// update state vectors which do not live in the time integrator to new layout
    virtual void UpdateStates(){};

    /// setup importer for dof and node based vectors
    virtual void SetupImporter();

    /// build extended bin column map
    virtual Teuchos::RCP<Epetra_Map> ExtendedBinColMap();

    /// ghosting of an additional second bin layer
    virtual void AddLayerGhosting(std::set<int>& colbins);

    /// ghosting of bins in proximity of ghosted boundary particles
    virtual void BdryParticleGhosting(std::set<int>& colbins);

    /// ghosting of bins in proximity of ghosted wall elements
    virtual void WallElementGhosting(std::set<int>& colbins);

    /// build connectivity from particle wall elements to bins
    virtual void BuildElementToBinPointers(
        bool wallpointer  ///< bool whether wall element pointer need rebuild
    );

    /// assign wall elements and gids to bins
    virtual void AssignWallElesAndGidsToBins();

    /// safety check for proper bin size
    virtual void BinSizeSafetyCheck();

    /// access structure and setup particle wall
    virtual void AccessStructure();

    /// particle walls are set up and assigned to bins
    virtual void SetupParticleWalls(
        Teuchos::RCP<DRT::Discretization> basediscret, const std::string elename);

    /// relate wall gids to bin ids
    virtual void RelateWallGidsToBinIds();

    /// get particles and wall elements in given bins
    void GetBinContent(std::list<DRT::Node*>& bin_p,
        boost::unordered_map<int, DRT::Element*>* bin_w, std::vector<int>& binIds) const;

    /// determine particle interaction distance
    virtual double ParticleInteractionDistance();

    //! @name particle related stuff
    //@{

    /// particle time integration
    Teuchos::RCP<ADAPTER::Particle> particles_;

    ///< importer for dof based vectors
    Teuchos::RCP<Epetra_Import> dofimporter_;

    ///< importer for node based vectors
    Teuchos::RCP<Epetra_Import> nodeimporter_;

    /// acceleration due to gravity
    LINALG::Matrix<3, 1> gravity_acc_;

    /// output interval
    int writeresultsevery_;

    /// fully redundant particle wall discretization
    Teuchos::RCP<DRT::Discretization> particlewalldis_;

    /// wall ele col map containing standard ghosting
    Teuchos::RCP<Epetra_Map> particlewallelecolmap_standardghosting_;

    /// map extractor for particle wall
    Teuchos::RCP<LINALG::MapExtractor> wallextractor_;

    /// wall displacements at \f$t_{n}\f$
    Teuchos::RCP<const Epetra_Vector> walldispn_;

    /// wall displacements at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> walldispnp_;

    /// wall velocities at \f$t_{n+1}\f$
    Teuchos::RCP<const Epetra_Vector> wallvelnp_;

    /// flag whether walls are moving
    const bool moving_walls_;

    /// type of repartitioning strategy
    const INPAR::PARTICLE::RepartitionStrategy rep_strategy_;

    /// displacement state of particles since last redistribution
    Teuchos::RCP<Epetra_Vector> dis_at_last_redistr_;

    /// type of particle interaction
    const INPAR::PARTICLE::ParticleInteractions particleInteractionType_;

    /// extended ghosting strategy
    const INPAR::PARTICLE::ExtendedGhosting extendedGhosting_;

    /// vector of pointers to particle material parameters
    std::vector<const MAT::PAR::ParticleMat*> particleMat_;

    BINSTRATEGY::UTILS::BinContentType bin_wallcontent_;

    /// handler of the rendering discretization (used for SPH dynamics)
    Teuchos::RCP<Rendering> rendering_;

    /// map relating wall gids to bin ids
    std::map<int, std::set<int>> relwallgidtobinids_;

    //@}

  };  // Algorithm

}  // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // PARTICLE_ALGORITHM_H
