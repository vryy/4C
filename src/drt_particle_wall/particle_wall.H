/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall handler for particle simulations

\level 2

\maintainer Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_H
#define PARTICLE_WALL_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_wall_interface.H"

#include "../drt_particle_engine/particle_typedefs.H"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
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
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEWALL
{
  /*!
   * \brief particle wall handler for particle problem
   *
   * The particle wall handler is responsible for the parallel distribution of all wall elements to
   * all processors following the bin distribution as handed in from the particle engine. In
   * addition potential particle wall neighbor pair relations are build.
   *
   * \author Sebastian Fuchs \date 10/2018
   */
  class WallHandlerBase : public WallHandlerInterface
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit WallHandlerBase(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~WallHandlerBase() override;

    /*!
     * \brief init wall handler
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \param[in] binstrategy binning strategy
     */
    virtual void Init(const std::shared_ptr<BINSTRATEGY::BinningStrategy> binstrategy);

    /*!
     * \brief setup wall handler
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \param[in] particlestatestotypes particle types and corresponding particle states
     */
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief write restart of wall handler
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    virtual void WriteRestart(const int step, const double time) const;

    /*!
     * \brief read restart of wall handler
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] restartstep restart step
     */
    virtual void ReadRestart(const int restartstep);

    /*!
     * \brief insert wall handler dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 05/2019
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const final;

    /*!
     * \brief write wall runtime output
     *
     * \author Sebastian Fuchs \date 08/2019
     *
     * \param[in] step output step
     * \param[in] time output time
     */
    virtual void WriteWallRuntimeOutput(const int step, const double time) const final;

    /*!
     * \brief update bin row and column map
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] binrowmap bin row map
     * \param[in] bincolmap bin column map
     */
    virtual void UpdateBinRowAndColMap(
        const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap) final;

    /*!
     * \brief check that wall nodes are located in bounding box
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    virtual void CheckWallNodesLocatedInBoundingBox() const final;

    /*!
     * \brief get maximum wall position increment since last transfer
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[out] allprocmaxpositionincrement maximum wall position increment
     */
    virtual void GetMaxWallPositionIncrement(double& allprocmaxpositionincrement) const final;

    //! distribute wall elements and nodes
    virtual void DistributeWallElementsAndNodes() = 0;

    //! transfer wall elements and nodes
    virtual void TransferWallElementsAndNodes() = 0;

    /*!
     * \brief relate bins to column wall elements
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    virtual void RelateBinsToColWallEles() final;

    /*!
     * \brief build particle to wall neighbors
     *
     * Build potential particle to wall neighbor pairs.
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] particlestobins relation of (owned and ghosted) particles to bins
     */
    virtual void BuildParticleToWallNeighbors(
        const PARTICLEENGINE::ParticlesToBins& particlestobins) final;

    /*!
     * \brief check for valid wall neighbors
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    virtual bool HaveValidWallNeighbors() const final { return validwallneighbors_; };

    Teuchos::RCP<const DRT::Discretization> GetWallDiscretization() const final
    {
      return walldiscretization_;
    };

    std::shared_ptr<PARTICLEWALL::WallDataState> GetWallDataState() const final
    {
      return walldatastate_;
    }

    const PARTICLEENGINE::PotentialWallNeighbors& GetPotentialWallNeighbors() const final;

    /*!
     * \brief determine nodal positions of column wall element
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param ele[in]             column wall element
     * \param colelenodalpos[out] current nodal position
     */
    void DetermineColWallEleNodalPos(
        DRT::Element* ele, std::map<int, LINALG::Matrix<3, 1>>& colelenodalpos) const final;

   private:
    //! \name init and setup methods
    //! @{

    //! init wall discretization
    virtual void InitWallDiscretization() = 0;

    //! setup wall discretization
    virtual void SetupWallDiscretization() const = 0;

    /*!
     * \brief init wall data state container
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    virtual void InitWallDataState() final;

    /*!
     * \brief init wall discretization runtime vtu writer
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    virtual void InitWallDiscretizationRuntimeVtuWriter() final;

    //! @}

    //! relate bins to column wall elements
    PARTICLEENGINE::BinsToColWallEles binstocolwalleles_;

    //! relate potential neighboring column wall elements to particles of all types
    PARTICLEENGINE::PotentialWallNeighbors potentialwallneighbors_;

    //! wall discretization runtime vtu writer
    std::unique_ptr<PARTICLEWALL::WallDiscretizationRuntimeVtuWriter>
        walldiscretizationruntimevtuwriter_;

   protected:
    //! create wall discretization
    virtual void CreateWallDiscretization() final;

    //! communicator
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

  /*!
   * \brief particle wall handler with wall from condition on discretization
   *
   * Particle wall handler with wall discretization generated from a surface condition on a
   * structure discretization.
   *
   * \author Sebastian Fuchs \date 10/2018
   */
  class WallHandlerDiscretCondition final : public WallHandlerBase
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit WallHandlerDiscretCondition(
        const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~WallHandlerDiscretCondition() override{};

    /*!
     * \brief distribute wall elements and nodes
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void DistributeWallElementsAndNodes() override;

    /*!
     * \brief transfer wall elements and nodes
     *
     * \author Sebastian Fuchs \date 03/2019
     */
    void TransferWallElementsAndNodes() override;

   private:
    /*!
     * \brief extend wall element ghosting
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[in] bintorowelemap bin to row wall element distribution
     */
    void ExtendWallElementGhosting(std::map<int, std::set<int>>& bintorowelemap);

    /*!
     * \brief init wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void InitWallDiscretization() override;

    /*!
     * \brief setup wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void SetupWallDiscretization() const override;
  };

  /*!
   * \brief particle wall handler with wall from bounding box
   *
   * Particle wall handler with wall discretization generated from bounding box of binning strategy.
   *
   * \author Sebastian Fuchs \date 10/2018
   */
  class WallHandlerBoundingBox final : public WallHandlerBase
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \param[in] comm   communicator
     * \param[in] params particle simulation parameter list
     */
    explicit WallHandlerBoundingBox(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! destructor
    ~WallHandlerBoundingBox() override{};

    /*!
     * \brief distribute wall elements and nodes
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void DistributeWallElementsAndNodes() override;

    /*!
     * \brief transfer wall elements and nodes
     *
     * \author Sebastian Fuchs \date 03/2019
     */
    void TransferWallElementsAndNodes() override;

   private:
    /*!
     * \brief init wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void InitWallDiscretization() override;

    /*!
     * \brief setup wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void SetupWallDiscretization() const override;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
#endif
