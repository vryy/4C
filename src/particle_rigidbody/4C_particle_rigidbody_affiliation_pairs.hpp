/*---------------------------------------------------------------------------*/
/*! \file
\brief affiliation pair handler for rigid bodies
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_RIGIDBODY_AFFILIATION_PAIRS_HPP
#define FOUR_C_PARTICLE_RIGIDBODY_AFFILIATION_PAIRS_HPP

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_config.hpp"

#include <Epetra_Comm.h>

#include <memory>
#include <unordered_map>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace Core::IO
{
  class DiscretizationReader;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleRigidBody
{
  /*!
   * \brief affiliation pair handler for rigid bodies
   *
   * The affiliation pair handler relates the global ids of rigid particles to the corresponding
   * global ids of rigid bodies.
   *
   * \author Sebastian Fuchs \date 08/2020
   */
  class RigidBodyAffiliationPairs final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] comm   communicator
     */
    explicit RigidBodyAffiliationPairs(const Epetra_Comm& comm);

    /*!
     * \brief init affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void init();

    /*!
     * \brief setup affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief write restart of affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void write_restart() const;

    /*!
     * \brief read restart of affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] reader discretization reader
     */
    void read_restart(const std::shared_ptr<Core::IO::DiscretizationReader> reader);

    /*!
     * \brief get reference to affiliation pair data
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \return reference to affiliation pair data
     */
    inline std::unordered_map<int, int>& get_ref_to_affiliation_pair_data()
    {
      return affiliationdata_;
    };

    /*!
     * \brief distribute affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void distribute_affiliation_pairs();

    /*!
     * \brief communicate affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void communicate_affiliation_pairs();

   private:
    /*!
     * \brief communicate specific affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void communicate_specific_affiliation_pairs(
        const std::vector<std::vector<int>>& particletargets);

    /*!
     * \brief pack all affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] buffer buffer containing affiliation data
     */
    void pack_all_affiliation_pairs(std::vector<char>& buffer) const;

    /*!
     * \brief unpack affiliation pairs
     *
     * Unpack affiliation pairs relating rigid particles to rigid bodies.
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] buffer buffer containing affiliation data
     */
    void unpack_affiliation_pairs(const std::vector<char>& buffer);

    /*!
     * \brief add affiliation pair to buffer
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in,out] buffer    buffer containing affiliation data
     * \param[in]     globalid  global id of rigid particle
     * \param[in]     rigidbody rigid body
     */
    void add_affiliation_pair_to_buffer(
        std::vector<char>& buffer, int globalid, int rigidbody) const;

    //! communicator
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! affiliation pair data relating rigid particles to rigid bodies
    std::unordered_map<int, int> affiliationdata_;
  };
}  // namespace ParticleRigidBody

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
