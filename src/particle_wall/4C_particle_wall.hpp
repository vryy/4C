/*---------------------------------------------------------------------------*/
/*! \file
\brief particle wall handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_WALL_HPP
#define FOUR_C_PARTICLE_WALL_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"
#include "4C_particle_wall_interface.hpp"

#include <Epetra_Comm.h>
#include <Epetra_Map.h>
#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
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

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::Elements
{
  class Element;
}

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

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 10/2018
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
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
     * \param[in] restart_time          restart time of the simulation
     */
    virtual void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        double restart_time);

    /*!
     * \brief write restart of wall handler
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] step restart step
     * \param[in] time restart time
     */
    virtual void write_restart(const int step, const double time) const;

    /*!
     * \brief read restart of wall handler
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] restartstep restart step
     */
    virtual void read_restart(const int restartstep);

    /*!
     * \brief insert wall handler dependent states of all particle types
     *
     * \author Sebastian Fuchs \date 05/2019
     *
     * \param[out] particlestatestotypes map of particle types and corresponding states
     */
    virtual void insert_particle_states_of_particle_types(
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
    virtual void write_wall_runtime_output(const int step, const double time) const final;

    /*!
     * \brief update bin row and column map
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] binrowmap bin row map
     * \param[in] bincolmap bin column map
     */
    virtual void update_bin_row_and_col_map(
        const Teuchos::RCP<Epetra_Map> binrowmap, const Teuchos::RCP<Epetra_Map> bincolmap) final;

    /*!
     * \brief check that wall nodes are located in bounding box
     *
     * \author Sebastian Fuchs \date 08/2019
     */
    virtual void check_wall_nodes_located_in_bounding_box() const final;

    /*!
     * \brief get maximum wall position increment since last transfer
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[out] allprocmaxpositionincrement maximum wall position increment
     */
    virtual void get_max_wall_position_increment(double& allprocmaxpositionincrement) const final;

    //! distribute wall elements and nodes
    virtual void distribute_wall_elements_and_nodes() = 0;

    //! transfer wall elements and nodes
    virtual void transfer_wall_elements_and_nodes() = 0;

    /*!
     * \brief relate bins to column wall elements
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    virtual void relate_bins_to_col_wall_eles() final;

    /*!
     * \brief build particle to wall neighbors
     *
     * Build potential particle to wall neighbor pairs.
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param[in] particlestobins relation of (owned and ghosted) particles to bins
     */
    virtual void build_particle_to_wall_neighbors(
        const PARTICLEENGINE::ParticlesToBins& particlestobins) final;

    /*!
     * \brief check for valid wall neighbors
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    virtual bool have_valid_wall_neighbors() const final { return validwallneighbors_; };

    Teuchos::RCP<const Core::FE::Discretization> get_wall_discretization() const final
    {
      return walldiscretization_;
    };

    std::shared_ptr<PARTICLEWALL::WallDataState> GetWallDataState() const final
    {
      return walldatastate_;
    }

    const PARTICLEENGINE::PotentialWallNeighbors& get_potential_wall_neighbors() const final;

    /*!
     * \brief determine nodal positions of column wall element
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \param ele[in]             column wall element
     * \param colelenodalpos[out] current nodal position
     */
    void determine_col_wall_ele_nodal_pos(Core::Elements::Element* ele,
        std::map<int, Core::LinAlg::Matrix<3, 1>>& colelenodalpos) const final;

   private:
    //! \name init and setup methods
    //! @{

    //! init wall discretization
    virtual void init_wall_discretization() = 0;

    //! setup wall discretization
    virtual void setup_wall_discretization() const = 0;

    /*!
     * \brief init wall data state container
     *
     * \author Sebastian Fuchs \date 05/2019
     */
    virtual void init_wall_data_state() final;

    /*!
     * \brief create wall discretization runtime vtu writer
     *
     * \author Sebastian Fuchs \date 11/2018
     *
     * \note this has to be called during setup() as it requires the restart to be read already
     *
     * \param[in] restart_time   restart time of the simulation
     */
    virtual void create_wall_discretization_runtime_vtu_writer(double restart_time) final;

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
    virtual void create_wall_discretization() final;

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
    Teuchos::RCP<Core::FE::Discretization> walldiscretization_;

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

    /*!
     * \brief distribute wall elements and nodes
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void distribute_wall_elements_and_nodes() override;

    /*!
     * \brief transfer wall elements and nodes
     *
     * \author Sebastian Fuchs \date 03/2019
     */
    void transfer_wall_elements_and_nodes() override;

   private:
    /*!
     * \brief extend wall element ghosting
     *
     * \author Sebastian Fuchs \date 03/2019
     *
     * \param[in] bintorowelemap bin to row wall element distribution
     */
    void extend_wall_element_ghosting(std::map<int, std::set<int>>& bintorowelemap);

    /*!
     * \brief init wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void init_wall_discretization() override;

    /*!
     * \brief setup wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void setup_wall_discretization() const override;
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

    /*!
     * \brief distribute wall elements and nodes
     *
     * \author Sebastian Fuchs \date 11/2018
     */
    void distribute_wall_elements_and_nodes() override;

    /*!
     * \brief transfer wall elements and nodes
     *
     * \author Sebastian Fuchs \date 03/2019
     */
    void transfer_wall_elements_and_nodes() override;

   private:
    /*!
     * \brief init wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void init_wall_discretization() override;

    /*!
     * \brief setup wall discretization
     *
     * \author Sebastian Fuchs \date 10/2018
     */
    void setup_wall_discretization() const override;
  };

}  // namespace PARTICLEWALL

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
