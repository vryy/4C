/*---------------------------------------------------------------------------*/
/*! \file
\brief affiliation pair handler for rigid bodies
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_rigidbody_affiliation_pairs.hpp"

#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_io.hpp"
#include "4C_particle_engine_communication_utils.hpp"
#include "4C_particle_engine_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLERIGIDBODY::RigidBodyAffiliationPairs::RigidBodyAffiliationPairs(const Epetra_Comm& comm)
    : comm_(comm), myrank_(comm.MyPID())
{
  // empty constructor
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::Init()
{
  // nothing to do
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::write_restart() const
{
  // get bin discretization writer
  std::shared_ptr<CORE::IO::DiscretizationWriter> binwriter =
      particleengineinterface_->get_bin_discretization_writer();

  // prepare buffer
  Teuchos::RCP<std::vector<char>> buffer = Teuchos::rcp(new std::vector<char>);

  // pack all affiliation pairs
  if (not affiliationdata_.empty()) pack_all_affiliation_pairs(*buffer);

  binwriter->WriteCharVector("RigidBodyAffiliationData", buffer);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::read_restart(
    const std::shared_ptr<CORE::IO::DiscretizationReader> reader)
{
  // prepare buffer
  Teuchos::RCP<std::vector<char>> buffer = Teuchos::rcp(new std::vector<char>);

  reader->ReadCharVector(buffer, "RigidBodyAffiliationData");

  // unpack affiliation pairs
  if (buffer->size() > 0) unpack_affiliation_pairs(*buffer);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::distribute_affiliation_pairs()
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

  // communicate specific affiliation pairs
  communicate_specific_affiliation_pairs(particletargets);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::communicate_affiliation_pairs()
{
  // get reference to particles being communicated to target processors
  const std::vector<std::vector<int>>& particletargets =
      particleengineinterface_->get_communicated_particle_targets();

  // communicate specific affiliation pairs
  communicate_specific_affiliation_pairs(particletargets);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::communicate_specific_affiliation_pairs(
    const std::vector<std::vector<int>>& particletargets)
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // pack affiliation pairs
  for (int torank = 0; torank < comm_.NumProc(); ++torank)
  {
    if (particletargets[torank].empty()) continue;

    for (int globalid : particletargets[torank])
    {
      auto it = affiliationdata_.find(globalid);

      // no affiliation pair for current global id
      if (it == affiliationdata_.end()) continue;

      // add affiliation pair to buffer
      add_affiliation_pair_to_buffer(sdata[torank], globalid, it->second);

      // erase affiliation pair
      affiliationdata_.erase(it);
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack affiliation pairs
  for (auto& p : rdata) unpack_affiliation_pairs(p.second);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::pack_all_affiliation_pairs(
    std::vector<char>& buffer) const
{
  // iterate over all affiliation pairs
  for (auto& it : affiliationdata_)
  {
    // add affiliation pair to buffer
    add_affiliation_pair_to_buffer(buffer, it.first, it.second);
  }
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::unpack_affiliation_pairs(
    const std::vector<char>& buffer)
{
  std::vector<char>::size_type position = 0;
  while (position < buffer.size())
  {
    // get affiliation pair
    const int globalid = CORE::COMM::ParObject::ExtractInt(position, buffer);
    const int rigidbody = CORE::COMM::ParObject::ExtractInt(position, buffer);

    // add affiliation pair
    affiliationdata_[globalid] = rigidbody;
  }
  if (position != buffer.size())
    FOUR_C_THROW("mismatch in size of data %d <-> %d", static_cast<int>(buffer.size()), position);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::add_affiliation_pair_to_buffer(
    std::vector<char>& buffer, int globalid, int rigidbody) const
{
  CORE::COMM::PackBuffer data;
  data.StartPacking();

  // add affiliation pair
  data.AddtoPack(globalid);
  data.AddtoPack(rigidbody);

  // append packed affiliation pair to buffer
  buffer.insert(buffer.end(), data().begin(), data().end());
}

FOUR_C_NAMESPACE_CLOSE
