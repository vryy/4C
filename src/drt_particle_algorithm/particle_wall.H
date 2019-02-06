/*---------------------------------------------------------------------------*/
/*!
\file particle_wall.H

\brief particle wall handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

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

#include "../linalg/linalg_fixedsizematrix.H"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>

#include "../linalg/linalg_utils.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 10/2018 |
 *---------------------------------------------------------------------------*/
class DiscretizationRuntimeVtuWriter;

namespace IO
{
  class DiscretizationReader;
}

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
namespace PARTICLEALGORITHM
{
  class WallHandlerBase : public WallHandlerInterface
  {
   public:
    //! constructor
    explicit WallHandlerBase(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~WallHandlerBase();

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

    //! write wall runtime vtu output
    virtual void WriteWallRuntimeVtuOutput(const int step, const double time) const;

    //! update bin row and column map
    virtual void UpdateBinRowAndColMap(
        const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap) final;

    //! distribute wall elements and nodes
    virtual void DistributeWallElementsAndNodes() = 0;

    //! relate bins to column wall elements
    virtual void RelateBinsToColWallEles() final;

    //! build particle to wall neighbors
    virtual void BuildParticleToWallNeighbors(
        const PARTICLEENGINE::ParticlesToBins& particlestobins) final;

    //! get reference to potential wall neighbors
    virtual const PARTICLEENGINE::PotentialWallNeighbors& GetPotentialWallNeighbors()
        const override final;

   private:
    //! init wall discretization
    virtual void InitWallDiscretization() final;

    //! setup wall discretization
    virtual void SetupWallDiscretization() const = 0;

    //! init wall discretization runtime vtu writer
    virtual void InitWallVtuWriter() final;

    //! setup wall discretization runtime vtu writer
    virtual void SetupWallVtuWriter() final;

    //! determine nodal positions of column wall element
    virtual void DetermineColWallEleNodalPos(
        DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& colelenodalpos) final;

    //! setup time of runtime vtu writer
    double setuptime_;

    //! relate bins to column wall elements
    PARTICLEENGINE::BinsToColWallEles binstocolwalleles_;

    //! relate potential neighboring column wall elements to particles of all types
    PARTICLEENGINE::PotentialWallNeighbors potentialwallneighbors_;

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

    //! wall discretization runtime vtu writer
    std::unique_ptr<DiscretizationRuntimeVtuWriter> wallvtuwriter_;

    //! flag denoting valid relation of bins to column wall elements
    bool validwallelements_;

    //! flag denoting valid relation of wall neighbors
    bool validwallneighbors_;
  };

  class WallHandlerDiscretCondition : public WallHandlerBase
  {
   public:
    //! constructor
    explicit WallHandlerDiscretCondition(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~WallHandlerDiscretCondition(){};

    //! distribute wall elements and nodes
    virtual void DistributeWallElementsAndNodes() override final;

   private:
    //! setup wall discretization
    virtual void SetupWallDiscretization() const override;
  };

  class WallHandlerBoundingBox : public WallHandlerBase
  {
   public:
    //! constructor
    explicit WallHandlerBoundingBox(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~WallHandlerBoundingBox(){};

    //! distribute wall elements and nodes
    virtual void DistributeWallElementsAndNodes() override final;

   private:
    //! setup wall discretization
    virtual void SetupWallDiscretization() const override;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
