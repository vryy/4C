/*---------------------------------------------------------------------------*/
/*!
\file particle_engine.H

\brief particle engine to control particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ENGINE_H
#define PARTICLE_ENGINE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_engine_interface.H"

#include "particle_enums.H"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainer;
  class ParticleContainerBundle;
  class ParticleObject;
  class ParticleRuntimeVtpWriter;
}  // namespace PARTICLEENGINE

namespace BINSTRATEGY
{
  class BinningStrategy;
}

namespace IO
{
  class DiscretizationReader;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngine : public ParticleEngineInterface
  {
   public:
    //! constructor
    explicit ParticleEngine(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~ParticleEngine();

    //! init particle engine
    void Init();

    //! setup particle engine
    void Setup(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    //! write restart of particle engine
    void WriteRestart(const int step, const double time) const;

    //! read restart of particle engine
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader,
        std::vector<ParticleObjShrdPtr>& particlestoread) const;

    //! write particle runtime vtp output
    void WriteParticleRuntimeVtpOutput(const int step, const double time) const;

    //! erase particles outside bounding box
    void EraseParticlesOutsideBoundingBox(std::vector<ParticleObjShrdPtr>& particlestocheck);

    //! distribute particles to owning processor
    void DistributeParticles(std::vector<ParticleObjShrdPtr>& particlestodistribute);

    //! transfer particles to new bins and processors
    void TransferParticles();

    //! ghost particles on other processors
    void GhostParticles();

    //! refresh particles being ghosted on other processors
    void RefreshParticles() const;

    //! refresh particles of specific states and types
    void RefreshParticlesOfSpecificStatesAndTypes(
        const StatesOfTypesToRefresh& particlestatestotypes) const override;

    //! dynamic load balancing
    void DynamicLoadBalancing();

    //! change type of particles
    void TypeChangeParticles(std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) override;

    //! build particle to particle neighbors
    void BuildParticleToParticleNeighbors();

    //! build global id to local index map
    void BuildGlobalIDToLocalIndexMap();

    //! check for valid particle connectivity
    bool HaveValidParticleConnectivity() const;

    //! get binning strategy
    const std::shared_ptr<BINSTRATEGY::BinningStrategy> GetBinningStrategy() const
    {
      return binstrategy_;
    };

    //! get bin row map
    const Teuchos::RCP<Epetra_Map> GetBinRowMap() const { return binrowmap_; };

    //! get bin col map
    const Teuchos::RCP<Epetra_Map> GetBinColMap() const { return bincolmap_; };

    //! get particle container bundle
    const ParticleContainerBundleShrdPtr GetParticleContainerBundle() const override
    {
      return particlecontainerbundle_;
    };

    //! get reference to (owned and ghosted) particles to bins
    const ParticlesToBins& GetParticlesToBins() const;

    //! get reference to potential particle neighbors
    const PotentialParticleNeighbors& GetPotentialParticleNeighbors() const override;

    //! get local index in specific particle container
    const LocalIndexTupleShrdPtr GetLocalIndexInSpecificContainer(int globalid) const override;

    //! return bin size
    const double* BinSize() const override;

    //! return minimum relevant bin size
    double MinBinSize() const override { return minbinsize_; };

    //! return flag whether pbc are applied
    bool HavePBC(const int dim = -1) const override;

    //! return delta for pbc in x, y, or z direction
    double PBCDelta(const int dim) const override;

    //! get bounding box dimensions
    LINALG::Matrix<3, 2>& XAABB() const override;

    //! get distance between particles considering periodic boundaries
    void DistanceBetweenParticles(
        const double* pos_i, const double* pos_j, double* r_ji) const override;

    //! create binning discretization reader
    const std::shared_ptr<IO::DiscretizationReader> BinDisReader(int restartstep) const;

    //! get number of particles on this processors
    int GetNumberOfParticles() const override;

    //! get number of particles on this processor of specific type
    int GetNumberOfParticlesOfSpecificType(const TypeEnum typeEnum) const override;

    //! write binning discretization output (debug feature)
    void WriteBinDisOutput(const int step, const double time) const;

   private:
    //! init binning strategy
    void InitBinningStrategy();

    //! setup binning strategy
    void SetupBinningStrategy();

    //! setup ghosting of bins
    void SetupBinGhosting();

    //! init particle container bundle
    void InitParticleContainerBundle();

    //! setup particle container bundle
    void SetupParticleContainerBundle(
        const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes) const;

    //! setup data storage
    void SetupDataStorage(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    //! init particle runtime vtp writer
    void InitParticleVtpWriter();

    //! setup particle runtime vtp writer
    void SetupParticleVtpWriter() const;

    //! setup particle type weights for dynamic load balancing
    void SetupTypeWeights();

    //! determine bin distribution dependent maps/sets
    void DetermineBinDisDependentMapsAndSets();

    //! determine ghosting dependent maps/sets for communication
    void DetermineGhostingDependentMapsAndSets();

    //! relate half neighboring bins to owned bins
    void RelateHalfNeighboringBinsToOwnedBins();

    //! check particles for periodic boundaries/leaving domain
    void CheckParticlesAtBoundaries(std::vector<std::set<int>>& particlestoremove) const;

    //! determine particles that need to be distributed
    void DetermineParticlesToBeDistributed(std::vector<ParticleObjShrdPtr>& particlestodistribute,
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestokeep) const;

    //! determine particles that need to be transfered
    void DetermineParticlesToBeTransfered(std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    //! determine particles that need to be ghosted
    void DetermineParticlesToBeGhosted(
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    //! determine particles that need to be refreshed
    void DetermineParticlesToBeRefreshed(
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    //! determine particles that need to be refreshed
    void DetermineSpecificStatesOfParticlesOfSpecificTypesToBeRefreshed(
        const StatesOfTypesToRefresh& particlestatestotypes,
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    //! communicate particles
    void CommunicateParticles(std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoreceive) const;

    //! communicate and build map for direct ghosting
    void CommunicateDirectGhostingMap(
        std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting);

    //! insert owned particles received from other processors
    void InsertOwnedParticles(
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert);

    //! insert ghosted particles received from other processors
    void InsertGhostedParticles(
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert,
        std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting);

    //! insert refreshed particles received from other processors
    void InsertRefreshedParticles(
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) const;

    //! remove particles from containers
    void RemoveParticlesFromContainers(std::vector<std::set<int>>& particlestoremove);

    //! store particle positions after transfer of particles
    void StorePositionsAfterParticleTransfer();

    //! relate owned particles to bins
    void RelateOwnedParticlesToBins();

    //! determine minimum relevant bin size
    void DetermineMinRelevantBinSize();

    //! determine bin weights needed for repartitioning
    void DetermineBinWeights();

    //! invalidate particle safety flags
    void InvalidateParticleSafetyFlags();

    //! communication
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! binning strategy
    std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy_;

    //! distribution of row bins
    Teuchos::RCP<Epetra_Map> binrowmap_;

    //! distribution of col bins
    Teuchos::RCP<Epetra_Map> bincolmap_;

    //! minimum relevant bin size
    double minbinsize_;

    //! size of vectors indexed by particle types
    int typevectorsize_;

    //! vector of bin center coordinates
    Teuchos::RCP<Epetra_MultiVector> bincenters_;

    //! vector of bin weights
    Teuchos::RCP<Epetra_MultiVector> binweights_;

    //! vector of particle type weights for dynamic load balancing
    std::vector<double> typeweights_;

    //! particle container bundle
    ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! particle runtime vtp writer
    std::unique_ptr<PARTICLEENGINE::ParticleRuntimeVtpWriter> particlevtpwriter_;

    //! relate (owned and ghosted) particles to bins
    ParticlesToBins particlestobins_;

    //! relate potential particle neighbors of all types and statuses
    PotentialParticleNeighbors potentialparticleneighbors_;

    //! maps particle global ids to local index in specific particle container
    std::unordered_map<int, LocalIndexTupleShrdPtr> globalidtolocalindex_;

    //! relate local index of owned particles to other processors local index of ghosted particles
    std::vector<std::map<int, std::vector<std::pair<int, int>>>> directghostingtargets_;

    //! relate half surrounding neighboring bins (including owned bin itself) to owned bins
    std::vector<std::set<int>> halfneighboringbinstobins_;

    //! flag denoting valid relation of owned particles to bins
    bool validownedparticles_;

    //! flag denoting valid relation of ghosted particles to bins
    bool validghostedparticles_;

    //! flag denoting valid relation of particle neighbors
    bool validparticleneighbors_;

    //! flag denoting validity of map relating particle global ids to local index
    bool validglobalidtolocalindex_;

    //! flag denoting validity of direct ghosting
    bool validdirectghosting_;

    //! flag denoting valid relation of half surrounding neighboring bins to owned bins
    bool validhalfneighboringbins_;

    //! owned bins at an open boundary or a periodic boundary
    std::set<int> boundarybins_;

    //! owned bins touched by other processors
    std::set<int> touchedbins_;

    //! maps bins of surrounding first layer to owning processors
    std::map<int, int> firstlayerbinsownedby_;

    //! bins being ghosted on this processor
    std::set<int> ghostedbins_;

    //! maps bins on this processor to processors ghosting that bins
    std::map<int, std::set<int>> thisbinsghostedby_;
  };

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
