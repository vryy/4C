/*---------------------------------------------------------------------------*/
/*! \file
\brief unique global identifier handler for particle simulations
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_UNIQUE_GLOBAL_ID_HANDLER_H
#define PARTICLE_UNIQUE_GLOBAL_ID_HANDLER_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "Epetra_Comm.h"

#include <map>
#include <vector>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
class RuntimeVtpWriter;

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
   * \brief unique global identifier handler
   *
   * A class that handles unique global identifiers (or also denoted global ids) for objects (e.g.,
   * particles) in particle simulations, meaning assignment of global ids to newly generated
   * objects, collecting of global ids of removed objects, and re-usage of free global ids.
   *
   * \note This requires communication among the processors in order to ensure that all global ids
   * are unique.
   *
   * \author Sebastian Fuchs \date 11/2019
   */
  class UniqueGlobalIdHandler final
  {
   public:
    /*!
     * \brief constructor
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] comm communicator
     */
    explicit UniqueGlobalIdHandler(const Epetra_Comm& comm, const std::string& objectname);

    //! destructor
    ~UniqueGlobalIdHandler() = default;

    /*!
     * \brief init unique global identifier handler
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void Init();

    /*!
     * \brief setup unique global identifier handler
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void Setup();

    /*!
     * \brief write restart of unique global identifier handler
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] writer discretization writer
     */
    void WriteRestart(std::shared_ptr<IO::DiscretizationWriter> writer) const;

    /*!
     * \brief read restart of unique global identifier handler
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] reader discretization reader
     */
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    /*!
     * \brief return maximum global id
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \return maximum global id
     */
    inline int GetMaxGlobalId() { return maxglobalid_; }

    /*!
     * \brief return number of reusable global ids on this processor
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \return number of reusable global ids on this processor
     */
    inline int GetNumberOfReusableGlobalIds() { return reusableglobalids_.size(); }

    /*!
     * \brief insert freed global id
     *
     * The global id of an object that is removed from the simulation is freed and stored herein to
     * be reassigned in case of new object generation.
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in] freedglobalid freed global id
     */
    inline void InsertFreedGlobalId(int freedglobalid)
    {
      reusableglobalids_.push_back(freedglobalid);
    }

    /*!
     * \brief draw requested number of global ids
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in,out] requesteduniqueglobalids requested unique global ids
     */
    void DrawRequestedNumberOfGlobalIds(std::vector<int>& requesteduniqueglobalids);

   private:
    /*!
     * \brief gather reusable global ids from all processors on master processor
     *
     * Gather all reusable global ids from all processors on the master processor. The gathered
     * reusable global ids are sorted on the master processor. On all non-master processors the
     * storage of reusable global ids is cleared.
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    void GatherReusableGlobalIdsFromAllProcsOnMasterProc();

    /*!
     * \brief prepare requested global ids for all processors
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in]     numberofrequestedgids number of requested global ids on this processor
     * \param[in,out] preparedglobalids     prepared global ids
     */
    void PrepareRequestedGlobalIdsForAllProcs(
        int numberofrequestedgids, std::map<int, std::vector<int>>& preparedglobalids);

    /*!
     * \brief extract requested global ids on master processor
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in]     preparedglobalids        prepared global ids
     * \param[in,out] requesteduniqueglobalids requested unique global ids
     */
    void ExtractRequestedGlobalIdsOnMasterProc(std::map<int, std::vector<int>>& preparedglobalids,
        std::vector<int>& requesteduniqueglobalids) const;

    /*!
     * \brief distribute requested global ids from master processor to all processors
     *
     * \author Sebastian Fuchs \date 11/2019
     *
     * \param[in]     tobesendglobalids        global ids to be send to each processor
     * \param[in,out] requesteduniqueglobalids requested unique global ids
     */
    void DistributeRequestedGlobalIdsFromMasterProcToAllProcs(
        std::map<int, std::vector<int>>& tobesendglobalids,
        std::vector<int>& requesteduniqueglobalids) const;

    //! communicator
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! master processor id
    const int masterrank_;

    //! object name
    const std::string objectname_;

    /*!
     * \brief maximum global id
     *
     * This is the maximum global id either actively used in the simulation by an object or as
     * present in the reusable global ids.
     *
     * \note This maximum global id is solely incremented on the master processor but broadcasted to
     * all other processors.
     *
     * \note The global ids are assigned starting from zero.
     *
     * \author Sebastian Fuchs \date 11/2019
     */
    int maxglobalid_;

    //! reusable global ids
    std::vector<int> reusableglobalids_;
  };
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
