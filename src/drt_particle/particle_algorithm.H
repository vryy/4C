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

// Ouput of minimal/maximal pressures, particle distances etc. every PARTICLE_MINMAXOUTPUT_EVRY time
// steps
// #define PARTICLE_MINMAXOUTPUT_EVRY 25

// Only throw warning (no dserror), if yield strength criterion is violated!
// #define PARTICLE_NOTIMESTEPDSERROR

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
  }
}  // namespace MAT

namespace BINSTRATEGY
{
  class BinningStrategy;
}

/*----------------------------------------------------------------------*
 | particle algorithm                                       ghamm 09/12 |
 *----------------------------------------------------------------------*/
namespace PARTICLE
{
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
    INPAR::PARTICLEOLD::ParticleInteractions ParticleInteractionType() const
    {
      return particleInteractionType_;
    }

    /// return vector of pointers to particle material parameters
    std::vector<const MAT::PAR::ParticleMat*> ParticleMat() const { return particleMat_; }

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
    const INPAR::PARTICLEOLD::RepartitionStrategy rep_strategy_;

    /// displacement state of particles since last redistribution
    Teuchos::RCP<Epetra_Vector> dis_at_last_redistr_;

    /// type of particle interaction
    const INPAR::PARTICLEOLD::ParticleInteractions particleInteractionType_;

    /// extended ghosting strategy
    const INPAR::PARTICLEOLD::ExtendedGhosting extendedGhosting_;

    /// vector of pointers to particle material parameters
    std::vector<const MAT::PAR::ParticleMat*> particleMat_;

    BINSTRATEGY::UTILS::BinContentType bin_wallcontent_;

    /// map relating wall gids to bin ids
    std::map<int, std::set<int>> relwallgidtobinids_;

    //@}

  };  // Algorithm

}  // namespace PARTICLE


/*----------------------------------------------------------------------*/
#endif  // PARTICLE_ALGORITHM_H
