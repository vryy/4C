/*---------------------------------------------------------------------------*/
/*! \file
\brief history pair handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_dem_history_pairs.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_io.hpp"
#include "4C_particle_engine_communication_utils.hpp"
#include "4C_particle_engine_interface.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::DEMHistoryPairs::DEMHistoryPairs(const Epetra_Comm& comm) : comm_(comm)
{
  // empty constructor
}

void PARTICLEINTERACTION::DEMHistoryPairs::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::DEMHistoryPairs::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEINTERACTION::DEMHistoryPairs::WriteRestart() const
{
  // get bin discretization writer
  std::shared_ptr<IO::DiscretizationWriter> binwriter =
      particleengineinterface_->get_bin_discretization_writer();

  // prepare buffer
  Teuchos::RCP<std::vector<char>> buffer;

  // particle tangential history data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    if (not particletangentialhistorydata_.empty())
      PackAllHistoryPairs(*buffer, particletangentialhistorydata_);

    binwriter->WriteCharVector("ParticleTangentialHistoryData", buffer);
  }

  // particle-wall tangential history pair data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    if (not particlewalltangentialhistorydata_.empty())
      PackAllHistoryPairs(*buffer, particlewalltangentialhistorydata_);

    binwriter->WriteCharVector("ParticleWallTangentialHistoryData", buffer);
  }

  // particle rolling history data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    if (not particlerollinghistorydata_.empty())
      PackAllHistoryPairs(*buffer, particlerollinghistorydata_);

    binwriter->WriteCharVector("ParticleRollingHistoryData", buffer);
  }

  // particle-wall rolling history pair data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    if (not particlewallrollinghistorydata_.empty())
      PackAllHistoryPairs(*buffer, particlewallrollinghistorydata_);

    binwriter->WriteCharVector("ParticleWallRollingHistoryData", buffer);
  }

  // particle adhesion history data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    if (not particleadhesionhistorydata_.empty())
      PackAllHistoryPairs(*buffer, particleadhesionhistorydata_);

    binwriter->WriteCharVector("ParticleAdhesionHistoryData", buffer);
  }

  // particle-wall adhesion history pair data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    if (not particlewalladhesionhistorydata_.empty())
      PackAllHistoryPairs(*buffer, particlewalladhesionhistorydata_);

    binwriter->WriteCharVector("ParticleWallAdhesionHistoryData", buffer);
  }
}

void PARTICLEINTERACTION::DEMHistoryPairs::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // prepare buffer
  Teuchos::RCP<std::vector<char>> buffer;

  // particle tangential history data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    reader->ReadCharVector(buffer, "ParticleTangentialHistoryData");

    if (buffer->size() > 0) UnpackHistoryPairs(*buffer, particletangentialhistorydata_);
  }

  // particle-wall tangential history pair data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    reader->ReadCharVector(buffer, "ParticleWallTangentialHistoryData");

    if (buffer->size() > 0) UnpackHistoryPairs(*buffer, particlewalltangentialhistorydata_);
  }

  // particle rolling history data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    reader->ReadCharVector(buffer, "ParticleRollingHistoryData");

    if (buffer->size() > 0) UnpackHistoryPairs(*buffer, particlerollinghistorydata_);
  }

  // particle-wall rolling history pair data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    reader->ReadCharVector(buffer, "ParticleWallRollingHistoryData");

    if (buffer->size() > 0) UnpackHistoryPairs(*buffer, particlewallrollinghistorydata_);
  }

  // particle adhesion history data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    reader->ReadCharVector(buffer, "ParticleAdhesionHistoryData");

    if (buffer->size() > 0) UnpackHistoryPairs(*buffer, particleadhesionhistorydata_);
  }

  // particle-wall adhesion history pair data
  {
    buffer = Teuchos::rcp(new std::vector<char>);

    reader->ReadCharVector(buffer, "ParticleWallAdhesionHistoryData");

    if (buffer->size() > 0) UnpackHistoryPairs(*buffer, particlewalladhesionhistorydata_);
  }
}

void PARTICLEINTERACTION::DEMHistoryPairs::distribute_history_pairs()
{
  // relate all particles to all processors
  std::vector<int> particlestoproc(0);
  particleengineinterface_->relate_all_particles_to_all_procs(particlestoproc);

  // allocate memory
  std::vector<std::vector<int>> particletargets(comm_.NumProc());

  // iterate over all particle global ids
  for (int gid = 0; gid < static_cast<int>(particlestoproc.size()); ++gid)
  {
    // processor id of current particle
    const int currproc = particlestoproc[gid];

    // no need to send history pairs
    if (currproc == comm_.MyPID()) continue;

    // no particle with current global id in simulation
    if (currproc < 0) continue;

    // append current particle global id
    particletargets[currproc].push_back(gid);
  }

  // communicate specific history pairs
  communicate_specific_history_pairs(particletargets, particletangentialhistorydata_);
  communicate_specific_history_pairs(particletargets, particlewalltangentialhistorydata_);

  communicate_specific_history_pairs(particletargets, particlerollinghistorydata_);
  communicate_specific_history_pairs(particletargets, particlewallrollinghistorydata_);

  communicate_specific_history_pairs(particletargets, particleadhesionhistorydata_);
  communicate_specific_history_pairs(particletargets, particlewalladhesionhistorydata_);
}

void PARTICLEINTERACTION::DEMHistoryPairs::communicate_history_pairs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMHistoryPairs::communicate_history_pairs");

  // get reference to particles being communicated to target processors
  const std::vector<std::vector<int>>& particletargets =
      particleengineinterface_->get_communicated_particle_targets();

  // communicate specific history pairs
  communicate_specific_history_pairs(particletargets, particletangentialhistorydata_);
  communicate_specific_history_pairs(particletargets, particlewalltangentialhistorydata_);

  communicate_specific_history_pairs(particletargets, particlerollinghistorydata_);
  communicate_specific_history_pairs(particletargets, particlewallrollinghistorydata_);

  communicate_specific_history_pairs(particletargets, particleadhesionhistorydata_);
  communicate_specific_history_pairs(particletargets, particlewalladhesionhistorydata_);
}

void PARTICLEINTERACTION::DEMHistoryPairs::UpdateHistoryPairs()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::DEMHistoryPairs::UpdateHistoryPairs");

  // erase untouched history pairs
  if (not particletangentialhistorydata_.empty())
    erase_untouched_history_pairs(particletangentialhistorydata_);

  if (not particlewalltangentialhistorydata_.empty())
    erase_untouched_history_pairs(particlewalltangentialhistorydata_);

  if (not particlerollinghistorydata_.empty())
    erase_untouched_history_pairs(particlerollinghistorydata_);

  if (not particlewallrollinghistorydata_.empty())
    erase_untouched_history_pairs(particlewallrollinghistorydata_);

  if (not particleadhesionhistorydata_.empty())
    erase_untouched_history_pairs(particleadhesionhistorydata_);

  if (not particlewalladhesionhistorydata_.empty())
    erase_untouched_history_pairs(particlewalladhesionhistorydata_);
}

template <typename historypairtype>
void PARTICLEINTERACTION::DEMHistoryPairs::communicate_specific_history_pairs(
    const std::vector<std::vector<int>>& particletargets,
    std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>& historydata)
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // pack history pairs
  for (int torank = 0; torank < comm_.NumProc(); ++torank)
  {
    if (particletargets[torank].empty()) continue;

    for (int globalid : particletargets[torank])
    {
      // no history pairs for current global id
      if (not historydata.count(globalid)) continue;

      for (auto& it_j : historydata[globalid])
      {
        // get reference to history pair
        const historypairtype& historypair = (it_j.second).second;

        // add history pair to buffer
        add_history_pair_to_buffer(sdata[torank], globalid, it_j.first, historypair);
      }
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack history pairs
  for (auto& p : rdata) UnpackHistoryPairs(p.second, historydata);
}

template <typename historypairtype>
void PARTICLEINTERACTION::DEMHistoryPairs::erase_untouched_history_pairs(
    std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>& historydata)
{
  // iterate over nested unordered maps of stored history pairs
  for (auto it_i = historydata.begin(); it_i != historydata.end();)
  {
    for (auto it_j = (it_i->second).begin(); it_j != (it_i->second).end();)
    {
      // remove untouched history pair
      if ((it_j->second).first == false) (it_i->second).erase(it_j++);
      // invalidate touched history pair
      else
        ((it_j++)->second).first = false;
    }

    // no history pairs left: erase entry
    if ((it_i->second).empty()) historydata.erase(it_i++);
    // increment iterator
    else
      ++it_i;
  }
}

template <typename historypairtype>
void PARTICLEINTERACTION::DEMHistoryPairs::PackAllHistoryPairs(std::vector<char>& buffer,
    const std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>&
        historydata) const
{
  // iterate over nested unordered maps of stored history pairs
  for (auto& it_i : historydata)
  {
    for (auto& it_j : it_i.second)
    {
      // get reference to history pair
      const historypairtype& historypair = (it_j.second).second;

      // add history pair to buffer
      add_history_pair_to_buffer(buffer, it_i.first, it_j.first, historypair);
    }
  }
}

template <typename historypairtype>
void PARTICLEINTERACTION::DEMHistoryPairs::UnpackHistoryPairs(const std::vector<char>& buffer,
    std::unordered_map<int, std::unordered_map<int, std::pair<bool, historypairtype>>>& historydata)
{
  std::vector<char>::size_type position = 0;
  while (position < buffer.size())
  {
    // get global ids
    int globalid_i = CORE::COMM::ParObject::ExtractInt(position, buffer);
    int globalid_j = CORE::COMM::ParObject::ExtractInt(position, buffer);

    // unpack history pair data
    historypairtype historypair = historypairtype();
    historypair.Unpack(position, buffer);

    // add history pair data
    historydata[globalid_i][globalid_j] = std::make_pair(true, historypair);
  }
  if (position != buffer.size())
    FOUR_C_THROW("mismatch in size of data %d <-> %d", static_cast<int>(buffer.size()), position);
}

template <typename historypairtype>
void PARTICLEINTERACTION::DEMHistoryPairs::add_history_pair_to_buffer(std::vector<char>& buffer,
    int globalid_i, int globalid_j, const historypairtype& historypair) const
{
  CORE::COMM::PackBuffer data;
  data.StartPacking();

  // add global ids
  data.AddtoPack(globalid_i);
  data.AddtoPack(globalid_j);

  // pack history pair data
  historypair.Pack(data);

  // append packed history pair to buffer
  buffer.insert(buffer.end(), data().begin(), data().end());
}

/*---------------------------------------------------------------------------*
 | template instantiations                                                   |
 *---------------------------------------------------------------------------*/
template void PARTICLEINTERACTION::DEMHistoryPairs::communicate_specific_history_pairs<
    PARTICLEINTERACTION::DEMHistoryPairTangential>(
    const std::vector<std::vector<int>>&, DEMHistoryPairTangentialData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::communicate_specific_history_pairs<
    PARTICLEINTERACTION::DEMHistoryPairRolling>(
    const std::vector<std::vector<int>>&, DEMHistoryPairRollingData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::communicate_specific_history_pairs<
    PARTICLEINTERACTION::DEMHistoryPairAdhesion>(
    const std::vector<std::vector<int>>&, DEMHistoryPairAdhesionData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::erase_untouched_history_pairs<
    PARTICLEINTERACTION::DEMHistoryPairTangential>(DEMHistoryPairTangentialData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::erase_untouched_history_pairs<
    PARTICLEINTERACTION::DEMHistoryPairRolling>(DEMHistoryPairRollingData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::erase_untouched_history_pairs<
    PARTICLEINTERACTION::DEMHistoryPairAdhesion>(DEMHistoryPairAdhesionData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::PackAllHistoryPairs<
    PARTICLEINTERACTION::DEMHistoryPairTangential>(
    std::vector<char>&, const DEMHistoryPairTangentialData&) const;

template void PARTICLEINTERACTION::DEMHistoryPairs::PackAllHistoryPairs<
    PARTICLEINTERACTION::DEMHistoryPairRolling>(
    std::vector<char>&, const DEMHistoryPairRollingData&) const;

template void PARTICLEINTERACTION::DEMHistoryPairs::PackAllHistoryPairs<
    PARTICLEINTERACTION::DEMHistoryPairAdhesion>(
    std::vector<char>&, const DEMHistoryPairAdhesionData&) const;

template void PARTICLEINTERACTION::DEMHistoryPairs::UnpackHistoryPairs<
    PARTICLEINTERACTION::DEMHistoryPairTangential>(
    const std::vector<char>&, DEMHistoryPairTangentialData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::UnpackHistoryPairs<
    PARTICLEINTERACTION::DEMHistoryPairRolling>(
    const std::vector<char>&, DEMHistoryPairRollingData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::UnpackHistoryPairs<
    PARTICLEINTERACTION::DEMHistoryPairAdhesion>(
    const std::vector<char>&, DEMHistoryPairAdhesionData&);

template void PARTICLEINTERACTION::DEMHistoryPairs::add_history_pair_to_buffer<
    PARTICLEINTERACTION::DEMHistoryPairTangential>(
    std::vector<char>&, int, int, const PARTICLEINTERACTION::DEMHistoryPairTangential&) const;

template void PARTICLEINTERACTION::DEMHistoryPairs::add_history_pair_to_buffer<
    PARTICLEINTERACTION::DEMHistoryPairRolling>(
    std::vector<char>&, int, int, const PARTICLEINTERACTION::DEMHistoryPairRolling&) const;

template void PARTICLEINTERACTION::DEMHistoryPairs::add_history_pair_to_buffer<
    PARTICLEINTERACTION::DEMHistoryPairAdhesion>(
    std::vector<char>&, int, int, const PARTICLEINTERACTION::DEMHistoryPairAdhesion&) const;

FOUR_C_NAMESPACE_CLOSE
