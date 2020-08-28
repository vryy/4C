/*---------------------------------------------------------------------------*/
/*! \file
\brief affiliation pair handler for rigid bodies
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_RIGIDBODY_AFFILIATION_PAIRS_H
#define PARTICLE_RIGIDBODY_AFFILIATION_PAIRS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include <Epetra_Comm.h>
#include <vector>
#include <unordered_map>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
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
namespace PARTICLERIGIDBODY
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

    //! destructor
    ~RigidBodyAffiliationPairs() = default;

    /*!
     * \brief init affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void Init();

    /*!
     * \brief setup affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    /*!
     * \brief write restart of affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void WriteRestart() const;

    /*!
     * \brief read restart of affiliation pair handler
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    /*!
     * \brief get reference to affiliation pair data
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \return reference to affiliation pair data
     */
    inline std::unordered_map<int, int>& GetRefToAffiliationPairData() { return affiliationdata_; };

    /*!
     * \brief distribute affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void DistributeAffiliationPairs();

    /*!
     * \brief communicate affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void CommunicateAffiliationPairs();

   private:
    /*!
     * \brief communicate specific affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     */
    void CommunicateSpecificAffiliationPairs(const std::vector<std::vector<int>>& particletargets);

    /*!
     * \brief pack all affiliation pairs
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] buffer buffer containing affiliation data
     */
    void PackAllAffiliationPairs(std::vector<char>& buffer) const;

    /*!
     * \brief unpack affiliation pairs
     *
     * Unpack affiliation pairs relating rigid particles to rigid bodies.
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in] buffer buffer containing affiliation data
     */
    void UnpackAffiliationPairs(const std::vector<char>& buffer);

    /*!
     * \brief add affiliation pair to buffer
     *
     * \author Sebastian Fuchs \date 08/2020
     *
     * \param[in,out] buffer    buffer containing affiliation data
     * \param[in]     globalid  global id of rigid particle
     * \param[in]     rigidbody rigid body
     */
    void AddAffiliationPairToBuffer(std::vector<char>& buffer, int globalid, int rigidbody) const;

    //! communicator
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! affiliation pair data relating rigid particles to rigid bodies
    std::unordered_map<int, int> affiliationdata_;
  };
}  // namespace PARTICLERIGIDBODY

/*---------------------------------------------------------------------------*/
#endif
