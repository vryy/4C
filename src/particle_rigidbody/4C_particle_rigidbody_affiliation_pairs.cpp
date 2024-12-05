// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_particle_rigidbody_affiliation_pairs.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_comm_pack_buffer.hpp"
#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_io.hpp"
#include "4C_particle_engine_communication_utils.hpp"
#include "4C_particle_engine_interface.hpp"
FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleRigidBody::RigidBodyAffiliationPairs::RigidBodyAffiliationPairs(MPI_Comm comm)
    : comm_(comm), myrank_(Core::Communication::my_mpi_rank(comm))
{
  // empty constructor
}

void ParticleRigidBody::RigidBodyAffiliationPairs::init()
{
  // nothing to do
}

void ParticleRigidBody::RigidBodyAffiliationPairs::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void ParticleRigidBody::RigidBodyAffiliationPairs::write_restart() const
{
  // get bin discretization writer
  std::shared_ptr<Core::IO::DiscretizationWriter> binwriter =
      particleengineinterface_->get_bin_discretization_writer();

  // prepare buffer
  std::shared_ptr<std::vector<char>> buffer = std::make_shared<std::vector<char>>();

  // pack all affiliation pairs
  if (not affiliationdata_.empty()) pack_all_affiliation_pairs(*buffer);

  binwriter->write_char_data("RigidBodyAffiliationData", *buffer);
}

void ParticleRigidBody::RigidBodyAffiliationPairs::read_restart(
    const std::shared_ptr<Core::IO::DiscretizationReader> reader)
{
  // prepare buffer
  std::shared_ptr<std::vector<char>> buffer = std::make_shared<std::vector<char>>();

  reader->read_char_vector(buffer, "RigidBodyAffiliationData");

  // unpack affiliation pairs
  if (buffer->size() > 0) unpack_affiliation_pairs(*buffer);
}

void ParticleRigidBody::RigidBodyAffiliationPairs::distribute_affiliation_pairs()
{
  // relate all particles to all processors
  std::vector<int> particlestoproc(0);
  particleengineinterface_->relate_all_particles_to_all_procs(particlestoproc);

  // allocate memory
  std::vector<std::vector<int>> particletargets(Core::Communication::num_mpi_ranks(comm_));

  // iterate over all particle global ids
  for (int gid = 0; gid < static_cast<int>(particlestoproc.size()); ++gid)
  {
    // processor id of current particle
    const int currproc = particlestoproc[gid];

    // no need to send history pairs
    if (currproc == Core::Communication::my_mpi_rank(comm_)) continue;

    // no particle with current global id in simulation
    if (currproc < 0) continue;

    // append current particle global id
    particletargets[currproc].push_back(gid);
  }

  // communicate specific affiliation pairs
  communicate_specific_affiliation_pairs(particletargets);
}

void ParticleRigidBody::RigidBodyAffiliationPairs::communicate_affiliation_pairs()
{
  // get reference to particles being communicated to target processors
  const std::vector<std::vector<int>>& particletargets =
      particleengineinterface_->get_communicated_particle_targets();

  // communicate specific affiliation pairs
  communicate_specific_affiliation_pairs(particletargets);
}

void ParticleRigidBody::RigidBodyAffiliationPairs::communicate_specific_affiliation_pairs(
    const std::vector<std::vector<int>>& particletargets)
{
  // prepare buffer for sending and receiving
  std::map<int, std::vector<char>> sdata;
  std::map<int, std::vector<char>> rdata;

  // pack affiliation pairs
  for (int torank = 0; torank < Core::Communication::num_mpi_ranks(comm_); ++torank)
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
  PARTICLEENGINE::COMMUNICATION::immediate_recv_blocking_send(comm_, sdata, rdata);

  // unpack affiliation pairs
  for (auto& p : rdata) unpack_affiliation_pairs(p.second);
}

void ParticleRigidBody::RigidBodyAffiliationPairs::pack_all_affiliation_pairs(
    std::vector<char>& buffer) const
{
  // iterate over all affiliation pairs
  for (auto& it : affiliationdata_)
  {
    // add affiliation pair to buffer
    add_affiliation_pair_to_buffer(buffer, it.first, it.second);
  }
}

void ParticleRigidBody::RigidBodyAffiliationPairs::unpack_affiliation_pairs(
    const std::vector<char>& buffer)
{
  Core::Communication::UnpackBuffer data(buffer);
  while (!data.at_end())
  {
    // get affiliation pair
    int globalid;
    extract_from_pack(data, globalid);
    int rigidbody;
    extract_from_pack(data, rigidbody);

    // add affiliation pair
    affiliationdata_[globalid] = rigidbody;
  }
}

void ParticleRigidBody::RigidBodyAffiliationPairs::add_affiliation_pair_to_buffer(
    std::vector<char>& buffer, int globalid, int rigidbody) const
{
  Core::Communication::PackBuffer data;

  // add affiliation pair
  data.add_to_pack(globalid);
  data.add_to_pack(rigidbody);

  // append packed affiliation pair to buffer
  buffer.insert(buffer.end(), data().begin(), data().end());
}

FOUR_C_NAMESPACE_CLOSE
