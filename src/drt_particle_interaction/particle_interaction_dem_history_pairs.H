/*---------------------------------------------------------------------------*/
/*! \file
\brief history pair handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_HISTORY_PAIRS_H
#define PARTICLE_INTERACTION_DEM_HISTORY_PAIRS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "particle_interaction_dem_history_pair_struct.H"

#include <Epetra_Comm.h>
#include <unordered_map>

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
 | type definitions                                           sfuchs 03/2019 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  typedef std::pair<bool, PARTICLEINTERACTION::DEMHistoryPairTangential>
      TouchedDEMHistoryPairTangential;
  typedef std::unordered_map<int, std::unordered_map<int, TouchedDEMHistoryPairTangential>>
      DEMHistoryPairTangentialData;

  typedef std::pair<bool, PARTICLEINTERACTION::DEMHistoryPairRolling> TouchedDEMHistoryPairRolling;
  typedef std::unordered_map<int, std::unordered_map<int, TouchedDEMHistoryPairRolling>>
      DEMHistoryPairRollingData;

  typedef std::pair<bool, PARTICLEINTERACTION::DEMHistoryPairAdhesion>
      TouchedDEMHistoryPairAdhesion;
  typedef std::unordered_map<int, std::unordered_map<int, TouchedDEMHistoryPairAdhesion>>
      DEMHistoryPairAdhesionData;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMHistoryPairs final
  {
   public:
    //! constructor
    explicit DEMHistoryPairs(const Epetra_Comm& comm);

    //! destructor
    ~DEMHistoryPairs(){};

    //! init history pair handler
    void Init();

    //! setup history pair handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface);

    //! write restart of history pair handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of history pair handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to particle tangential history pair data
    inline DEMHistoryPairTangentialData& GetRefToParticleTangentialHistoryData()
    {
      return particletangentialhistorydata_;
    };

    //! get reference to particle-wall tangential history pair data
    inline DEMHistoryPairTangentialData& GetRefToParticleWallTangentialHistoryData()
    {
      return particlewalltangentialhistorydata_;
    };

    //! get reference to particle rolling history pair data
    inline DEMHistoryPairRollingData& GetRefToParticleRollingHistoryData()
    {
      return particlerollinghistorydata_;
    };

    //! get reference to particle-wall rolling history pair data
    inline DEMHistoryPairRollingData& GetRefToParticleWallRollingHistoryData()
    {
      return particlewallrollinghistorydata_;
    };

    //! get reference to particle adhesion history pair data
    inline DEMHistoryPairAdhesionData& GetRefToParticleAdhesionHistoryData()
    {
      return particleadhesionhistorydata_;
    };

    //! get reference to particle-wall adhesion history pair data
    inline DEMHistoryPairAdhesionData& GetRefToParticleWallAdhesionHistoryData()
    {
      return particlewalladhesionhistorydata_;
    };

    //! distribute history pairs
    void DistributeHistoryPairs();

    //! communicate history pairs
    void CommunicateHistoryPairs();

    //! update history pairs
    void UpdateHistoryPairs();

   private:
    //! communicate specific history pairs
    template <typename historypairtype>
    void CommunicateSpecificHistoryPairs(const std::vector<std::vector<int>>& particletargets,
        std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>&
            historydata);

    //! erase untouched history pairs
    template <typename historypairtype>
    void EraseUntouchedHistoryPairs(
        std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>&
            historydata);

    //! pack all history pairs
    template <typename historypairtype>
    void PackAllHistoryPairs(std::vector<char>& buffer,
        const std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>&
            historydata) const;

    //! unpack history pairs
    template <typename historypairtype>
    void UnpackHistoryPairs(const std::vector<char>& buffer,
        std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>&
            historydata);

    //! add history pair to buffer
    template <typename historypairtype>
    void AddHistoryPairToBuffer(std::vector<char>& buffer, int globalid_i, int globalid_j,
        const historypairtype& historypair) const;

    //! communication
    const Epetra_Comm& comm_;

    //! particle tangential history pair data
    DEMHistoryPairTangentialData particletangentialhistorydata_;

    //! particle-wall tangential history pair data
    DEMHistoryPairTangentialData particlewalltangentialhistorydata_;

    //! particle rolling history pair data
    DEMHistoryPairRollingData particlerollinghistorydata_;

    //! particle-wall rolling history pair data
    DEMHistoryPairRollingData particlewallrollinghistorydata_;

    //! particle adhesion history pair data
    DEMHistoryPairAdhesionData particleadhesionhistorydata_;

    //! particle-wall adhesion history pair data
    DEMHistoryPairAdhesionData particlewalladhesionhistorydata_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
