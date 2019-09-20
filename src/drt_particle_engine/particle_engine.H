/*---------------------------------------------------------------------------*/
/*! \file
\brief particle engine to control particle problem

\level 1

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ENGINE_H
#define PARTICLE_ENGINE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_engine_interface.H"

#include <Epetra_Comm.h>
#include <Epetra_Vector.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
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
  class DiscretizationWriter;
  class DiscretizationReader;
}  // namespace IO

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  /*!
   * \brief particle engine to control particle problem
   *
   * The particle engine is responsible for the parallel distribution of particles to all
   * processors and keeps track of the dynamic change in position of all particles, including
   * communication of particles to other processors if required. In addition potential particle
   * neighbor pair relations are build.
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  class ParticleEngine final : public ParticleEngineInterface
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit ParticleEngine(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~ParticleEngine() override;

    /*!
     * \brief init particle engine
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void Init();

    /*!
     * \brief setup particle engine
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] particlestatestotypes particle types and corresponding particle states
     */
    void Setup(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    /*!
     * \brief write restart of particle engine
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    void WriteRestart(const int step, const double time) const;

    /*!
     * \brief read restart of particle engine
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in]  reader          discretization reader
     * \param[out] particlestoread particle objects read in from restart
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader,
        std::vector<ParticleObjShrdPtr>& particlestoread) const;

    /*!
     * \brief write particle runtime output
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] step output step
     * \param[in] time output time
     */
    void WriteParticleRuntimeOutput(const int step, const double time) const;

    /*!
     * \brief erase particles outside bounding box
     *
     * All particles handed in are checked if they are located inside the bounding box. Particles
     * that are outside the bounding box are erased.
     *
     * \author Sebastian Fuchs \date 09/2018
     *
     * \param[in] particlestocheck particles to be checked if located in bounding box
     */
    void EraseParticlesOutsideBoundingBox(std::vector<ParticleObjShrdPtr>& particlestocheck);

    /*!
     * \brief distribute particles to owning processor
     *
     * All particle objects handed in are checked for their owning processor and send to that owning
     * processor. This method is primarily used for initialization or after redistribution of the
     * underlying bins.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] particlestodistribute particles to be distributed to their owning processor
     */
    void DistributeParticles(std::vector<ParticleObjShrdPtr>& particlestodistribute);

    /*!
     * \brief transfer particles to new bins and processors
     *
     * Check all owned particle for their current position. If a particle leaves the spatial domain
     * of its currently owning processor it either left the computational domain and is removed or
     * it entered the spatial domain of a neighboring processor and is transfered to this new owning
     * processor.
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void TransferParticles();

    /*!
     * \brief ghost particles on other processors
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void GhostParticles();

    /*!
     * \brief refresh particles being ghosted on other processors
     *
     * Make use of the direct ghosting targets map and send particles that are ghosted on other
     * processors directly to those processors together with the information of the local index in
     * the respective ghosted particle container.
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void RefreshParticles() const;

    /*!
     * \brief refresh specific states of particles of specific types
     *
     * Communicate the information of specific states of particles of specific types from owned
     * particles to processors ghosting that respective particles.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] particlestatestotypes particle types and corresponding particle states to be
     *                                  refreshed
     */
    void RefreshParticlesOfSpecificStatesAndTypes(
        const StatesOfTypesToRefresh& particlestatestotypes) const override;

    /*!
     * \brief dynamic load balancing
     *
     * Dynamically balance the work load over all processors based on the total number of particles
     * on each processor. This results in a new distribution of bins and requires to rebuild all
     * dependent maps and sets.
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void DynamicLoadBalancing();

    /*!
     * \brief change type of particles
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] particlestoremove particles to be removed from containers on this processor
     * \param[in] particlestoinsert particles to be inserted into containers on this processor
     */
    void TypeChangeParticles(std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) override;

    /*!
     * \brief build particle to particle neighbors
     *
     * Build potential particle to particle neighbor pairs storing each pair only once.
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void BuildParticleToParticleNeighbors();

    /*!
     * \brief build global id to local index map
     *
     * Relate the global id of all particles to the local index as stored in the particle container.
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void BuildGlobalIDToLocalIndexMap();

    /*!
     * \brief check for valid particle connectivity
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    bool HaveValidParticleConnectivity() const;

    /*!
     * \brief check for valid particle neighbors
     *
     * \author Sebastian Fuchs \date 04/2019
     */
    bool HaveValidParticleNeighbors() const;

    /*!
     * \brief get binning strategy
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \return binning strategy
     */
    const std::shared_ptr<BINSTRATEGY::BinningStrategy> GetBinningStrategy() const
    {
      return binstrategy_;
    };

    /*!
     * \brief get bin row map
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \return bin row map
     */
    const Teuchos::RCP<Epetra_Map> GetBinRowMap() const { return binrowmap_; };

    /*!
     * \brief get bin column map
     *
     * \author Sebastian Fuchs \date 02/2019
     *
     * \return bin column map
     */
    const Teuchos::RCP<Epetra_Map> GetBinColMap() const { return bincolmap_; };

    const ParticleContainerBundleShrdPtr GetParticleContainerBundle() const override
    {
      return particlecontainerbundle_;
    };

    /*!
     * \brief get reference to relation of (owned and ghosted) particles to bins
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \return relation of (owned and ghosted) particles to bins
     */
    const ParticlesToBins& GetParticlesToBins() const;

    const PotentialParticleNeighbors& GetPotentialParticleNeighbors() const override;

    const std::vector<std::vector<int>>& GetCommunicatedParticleTargets() const override
    {
      return communicatedparticletargets_;
    };

    const LocalIndexTupleShrdPtr GetLocalIndexInSpecificContainer(int globalid) const override;

    const std::shared_ptr<IO::DiscretizationWriter> GetBinDiscretizationWriter() const override;

    /*!
     * \brief relate all particles to all processors
     *
     * Relate the global id of all particles to the global id of the owning processor stored in a
     * vector indexed by the particle global id. In case a particle for a certain global id does not
     * exist the value is set to -1.
     *
     * \note This method is computational expensive and should be used only for initialization.
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[out] particlestoproc relate global id of particles to global id of processor
     */
    void RelateAllParticlesToAllProcs(std::vector<int>& particlestoproc) const override;

    const double* BinSize() const override;

    double MinBinSize() const override { return minbinsize_; };

    bool HavePBC(const int dim = -1) const override;

    double PBCDelta(const int dim) const override;

    LINALG::Matrix<3, 2>& XAABB() const override;

    /*!
     * \brief get distance between particles considering periodic boundaries
     *
     * Compute the distance vector between two given particles. In the case of periodic boundary
     * conditions it is checked if the distance over the periodic boundary is closer than inside the
     * domain.
     *
     * \author Sebastian Fuchs \date 08/2018
     *
     * \param[in]  pos_i pointer to position of particle i
     * \param[in]  pos_j pointer to position of particle j
     * \param[out] r_ji  vector from particle i to j
     */
    void DistanceBetweenParticles(
        const double* pos_i, const double* pos_j, double* r_ji) const override;

    /*!
     * \brief create binning discretization reader
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in] restartstep restart step
     *
     * \return discretization reader
     */
    const std::shared_ptr<IO::DiscretizationReader> BinDisReader(int restartstep) const;

    int GetNumberOfParticles() const override;

    int GetNumberOfParticlesOfSpecificType(const TypeEnum typeEnum) const override;

    /*!
     * \brief write binning discretization output
     *
     * \note The binning discretization output is solely meant as a debug feature as this method is
     * computational expensive!
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] step output step
     * \param[in] time output time
     */
    void WriteBinDisOutput(const int step, const double time) const;

   private:
    //! \name init and setup methods
    //! @{

    /*!
     * \brief init binning strategy
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void InitBinningStrategy();

    /*!
     * \brief setup binning strategy
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void SetupBinningStrategy();

    /*!
     * \brief setup ghosting of bins
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void SetupBinGhosting();

    /*!
     * \brief init particle container bundle
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void InitParticleContainerBundle();

    /*!
     * \brief setup particle container bundle
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] particlestatestotypes particle types and corresponding particle states
     */
    void SetupParticleContainerBundle(
        const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes) const;

    /*!
     * \brief setup data storage
     *
     * \author Sebastian Fuchs \date 12/2018
     *
     * \param[in] particlestatestotypes particle types and corresponding particle states
     */
    void SetupDataStorage(const std::map<TypeEnum, std::set<StateEnum>>& particlestatestotypes);

    /*!
     * \brief init particle runtime vtp writer
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void InitParticleVtpWriter();

    /*!
     * \brief setup particle runtime vtp writer
     *
     * \author Sebastian Fuchs \date 04/2018
     */
    void SetupParticleVtpWriter() const;

    /*!
     * \brief setup particle type weights for dynamic load balancing
     *
     * \author Sebastian Fuchs \date 03/2019
     */
    void SetupTypeWeights();

    //! @}

    //! \name bin distribution and ghosting
    //! @{

    /*!
     * \brief determine bin distribution dependent maps/sets
     *
     * Determine the following maps respectively sets that depend on the parallel distribution of
     * the bins to the processors:
     * - set of owned bins at an open boundary or a periodic boundary
     * - set of owned bins touched by other processors
     * - map of non-owned neighboring bins related to global id of neighboring owning processor
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void DetermineBinDisDependentMapsAndSets();

    /*!
     * \brief determine ghosting dependent maps/sets for communication
     *
     * Determine the following maps respectively sets that depend on the ghosting of the bins:
     * - set of bins being ghosted on this processor
     * - map of owned bins related to global ids of processors ghosting that bins
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void DetermineGhostingDependentMapsAndSets();

    /*!
     * \brief relate half neighboring bins to owned bins
     *
     * In order to evaluate potential particle neighbor pairs all neighboring bin pairs need to be
     * processed. It is sufficient to find each particle neighbor pair once and evaluate the
     * contribution of both particles onto each other, thus it also is sufficient to relate only
     * half of the surrounding neighboring bins to an owned bin following a specific stencil. In
     * addition, all ghosted bins and also the bin itself are added to that list in order to not
     * miss any potential particle pairs.
     *
     * \author Sebastian Fuchs \date 01/2019
     */
    void RelateHalfNeighboringBinsToOwnedBins();

    //! @}

    /*!
     * \brief check particles for periodic boundaries/leaving domain
     *
     * Check the current position of owned particles in the proximity of the boundary of the
     * computational domain. In case a particle left the computational domain it is removed from the
     * containers. In case a particle travels over a periodic boundary it is shifted by the periodic
     * length to be re-injected at the opposite site.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[out] particlestoremove particles to be removed from containers
     */
    void CheckParticlesAtBoundaries(std::vector<std::set<int>>& particlestoremove) const;

    /*!
     * \brief determine particles that need to be distributed
     *
     * Find the owning processors of all particle objects handed in.
     *
     * \author Sebastian Fuchs \date 04/2018
     *
     * \param[in]  particlestodistribute particles to be distributed to their owning processor
     * \param[out] particlestosend       particles to be send to other processors
     * \param[out] particlestokeep       particles to be kept at this processor
     */
    void DetermineParticlesToBeDistributed(std::vector<ParticleObjShrdPtr>& particlestodistribute,
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestokeep);

    /*!
     * \brief determine particles that need to be transfered
     *
     * Check the current position of owned particles close to the spatial domain of neighboring
     * processors. In case a particle left the spatial domain of its current owning processor it
     * needs to be transfered to its new owning processor. This requires that particles travel not
     * more than one bin layer as each processor only knows the owning processor of the first layer
     * of bins surrounding its own spatial domain.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[out] particlestoremove particles to be removed from containers on this processor
     * \param[out] particlestosend   particles to be send to other processors
     */
    void DetermineParticlesToBeTransfered(std::vector<std::set<int>>& particlestoremove,
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend);

    /*!
     * \brief determine particles that need to be ghosted
     *
     * Determine all particles that need to be ghosted on other processors based on the information
     * of owned bins, that are ghosted by other processors.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[out] particlestosend particles to be send to other processors
     */
    void DetermineParticlesToBeGhosted(
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    /*!
     * \brief determine particles that need to be refreshed
     *
     * Determine particles that need to be refreshed on other processors based on the direct
     * ghosting targets map.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[out] particlestosend particles to be send to other processors
     */
    void DetermineParticlesToBeRefreshed(
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    /*!
     * \brief determine specific states of particles of specific type that need to be refreshed
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in]  particlestatestotypes particle types and corresponding particle states to be
     *                                   refreshed
     * \param[out] particlestosend       particles to be send to other processors
     */
    void DetermineSpecificStatesOfParticlesOfSpecificTypesToBeRefreshed(
        const StatesOfTypesToRefresh& particlestatestotypes,
        std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend) const;

    /*!
     * \brief communicate particles
     *
     * This method communicates particles to be send to other processors and receives particles from
     * other processors.
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in]  particlestosend    particles to be send to other processors
     * \param[out] particlestoreceive particles to be received on this processor
     */
    void CommunicateParticles(std::vector<std::vector<ParticleObjShrdPtr>>& particlestosend,
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoreceive) const;

    /*!
     * \brief communicate and build map for direct ghosting
     *
     * Communicate the information at which local index in the particle container a particle at the
     * processor ghosting that particle is inserted to the owning processors.
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] directghosting direct ghosting information
     */
    void CommunicateDirectGhostingMap(
        std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting);

    /*!
     * \brief insert owned particles received from other processors
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] particlestoinsert particles to be inserted into containers on this processor
     */
    void InsertOwnedParticles(
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert);

    /*!
     * \brief insert ghosted particles received from other processors
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in]  particlestoinsert particles to be inserted into containers on this processor
     * \param[out] directghosting    direct ghosting information
     */
    void InsertGhostedParticles(
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert,
        std::map<int, std::map<TypeEnum, std::map<int, std::pair<int, int>>>>& directghosting);

    /*!
     * \brief insert refreshed particles received from other processors
     *
     * \author Sebastian Fuchs \date 05/2018
     *
     * \param[in] particlestoinsert particles to be inserted into containers on this processor
     */
    void InsertRefreshedParticles(
        std::vector<std::vector<std::pair<int, ParticleObjShrdPtr>>>& particlestoinsert) const;

    /*!
     * \brief remove particles from containers
     *
     * \author Sebastian Fuchs \date 03/2018
     *
     * \param[in] particlestoremove particles to be removed from containers on this processor
     */
    void RemoveParticlesFromContainers(std::vector<std::set<int>>& particlestoremove);

    /*!
     * \brief store particle positions after transfer of particles
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void StorePositionsAfterParticleTransfer();

    /*!
     * \brief relate owned particles to bins
     *
     * \author Sebastian Fuchs \date 03/2018
     */
    void RelateOwnedParticlesToBins();

    /*!
     * \brief determine minimum relevant bin size
     *
     * \author Sebastian Fuchs \date 08/2018
     */
    void DetermineMinRelevantBinSize();

    /*!
     * \brief determine bin weights needed for repartitioning
     *
     * \author Sebastian Fuchs \date 05/2018
     */
    void DetermineBinWeights();

    /*!
     * \brief invalidate particle safety flags
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void InvalidateParticleSafetyFlags();

    //! communicator
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! binning strategy
    std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy_;

    //! distribution of row bins
    Teuchos::RCP<Epetra_Map> binrowmap_;

    //! distribution of column bins
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

    //! owned particles being communicated (transfered/distributed) to target processors
    std::vector<std::vector<int>> communicatedparticletargets_;

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
