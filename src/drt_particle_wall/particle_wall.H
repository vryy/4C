/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_H
#define PARTICLE_WALL_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_wall_interface.H"

#include "../drt_particle_engine/particle_typedefs.H"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEWALL
{
  class WallDataState;
  class WallDiscretizationRuntimeVtuWriter;
}  // namespace PARTICLEWALL

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

namespace BINSTRATEGY
{
  class BinningStrategy;
}

namespace DRT
{
  class Discretization;
  class Element;
}  // namespace DRT

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  class WallHandlerBase : public WallHandlerInterface
  {
   public:
    //! constructor
    explicit WallHandlerBase(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~WallHandlerBase() override;

    //! init wall handler
    virtual void Init();

    //! setup wall handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy);

    //! write restart of wall handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of wall handler
    virtual void ReadRestart(const int restartstep);

    //! insert wall handler dependent states of all particle types
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const final;

    //! write wall runtime output
    virtual void WriteWallRuntimeOutput(const int step, const double time) const final;

    //! update bin row and column map
    virtual void UpdateBinRowAndColMap(
        const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap) final;

    //! check that wall nodes are located in bounding box
    virtual void CheckWallNodesLocatedInBoundingBox() const final;

    //! get max wall position increment since last transfer
    virtual void GetMaxWallPositionIncrement(double& allprocmaxpositionincrement) const final;

    //! distribute wall elements and nodes
    virtual void DistributeWallElementsAndNodes() = 0;

    //! transfer wall elements and nodes
    virtual void TransferWallElementsAndNodes() = 0;

    //! relate bins to column wall elements
    virtual void RelateBinsToColWallEles() final;

    //! build particle to wall neighbors
    virtual void BuildParticleToWallNeighbors(
        const PARTICLEENGINE::ParticlesToBins& particlestobins) final;

    //! check for valid wall neighbors
    virtual bool HaveValidWallNeighbors() const final { return validwallneighbors_; };

    //! get wall discretization
    Teuchos::RCP<const DRT::Discretization> GetWallDiscretization() const final
    {
      return walldiscretization_;
    };

    //! get wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> GetWallDataState() const final
    {
      return walldatastate_;
    }

    //! get reference to potential wall neighbors
    const PARTICLEENGINE::PotentialWallNeighbors& GetPotentialWallNeighbors() const final;

    //! determine nodal positions of column wall element
    void DetermineColWallEleNodalPos(
        DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& colelenodalpos) const final;

   private:
    //! init wall discretization
    virtual void InitWallDiscretization() final;

    //! setup wall discretization
    virtual void SetupWallDiscretization() const = 0;

    //! init wall data state container
    virtual void InitWallDataState() final;

    //! init wall discretization runtime vtu writer
    virtual void InitWallDiscretizationRuntimeVtuWriter() final;

    //! relate bins to column wall elements
    PARTICLEENGINE::BinsToColWallEles binstocolwalleles_;

    //! relate potential neighboring column wall elements to particles of all types
    PARTICLEENGINE::PotentialWallNeighbors potentialwallneighbors_;

    //! wall discretization runtime vtu writer
    std::unique_ptr<PARTICLEWALL::WallDiscretizationRuntimeVtuWriter>
        walldiscretizationruntimevtuwriter_;

   protected:
    //! communication
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! binning strategy
    std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy_;

    //! distribution of row bins
    Teuchos::RCP<Epetra_Map> binrowmap_;

    //! distribution of col bins
    Teuchos::RCP<Epetra_Map> bincolmap_;

    //! wall discretization
    Teuchos::RCP<DRT::Discretization> walldiscretization_;

    //! wall data state container
    std::shared_ptr<PARTICLEWALL::WallDataState> walldatastate_;

    //! flag denoting valid relation of bins to column wall elements
    bool validwallelements_;

    //! flag denoting valid relation of wall neighbors
    bool validwallneighbors_;
  };

  class WallHandlerDiscretCondition final : public WallHandlerBase
  {
   public:
    //! constructor
    explicit WallHandlerDiscretCondition(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~WallHandlerDiscretCondition() override{};

    //! distribute wall elements and nodes
    void DistributeWallElementsAndNodes() override;

    //! transfer wall elements and nodes
    void TransferWallElementsAndNodes() override;

   private:
    //! extend wall element ghosting
    void ExtendWallElementGhosting(std::map<int, std::set<int>>& bintorowelemap);

    //! setup wall discretization
    void SetupWallDiscretization() const override;
  };

  class WallHandlerBoundingBox final : public WallHandlerBase
  {
   public:
    //! constructor
    explicit WallHandlerBoundingBox(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~WallHandlerBoundingBox() override{};

    //! distribute wall elements and nodes
    void DistributeWallElementsAndNodes() override;

    //! transfer wall elements and nodes
    void TransferWallElementsAndNodes() override;

   private:
    //! setup wall discretization
    void SetupWallDiscretization() const override;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
