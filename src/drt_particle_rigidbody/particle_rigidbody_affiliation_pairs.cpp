/*---------------------------------------------------------------------------*/
/*! \file
\brief affiliation pair handler for rigid bodies
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_rigidbody_affiliation_pairs.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_communication_utils.H"

#include "../drt_io/io.H"
#include "../drt_lib/drt_pack_buffer.H"
#include "../drt_lib/drt_parobject.H"

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

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::WriteRestart() const
{
  // get bin discretization writer
  std::shared_ptr<IO::DiscretizationWriter> binwriter =
      particleengineinterface_->GetBinDiscretizationWriter();

  // prepare buffer
  Teuchos::RCP<std::vector<char>> buffer = Teuchos::rcp(new std::vector<char>);

  // pack all affiliation pairs
  if (not affiliationdata_.empty()) PackAllAffiliationPairs(*buffer);

  binwriter->WriteCharVector("RigidBodyAffiliationData", buffer);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // prepare buffer
  Teuchos::RCP<std::vector<char>> buffer = Teuchos::rcp(new std::vector<char>);

  reader->ReadCharVector(buffer, "RigidBodyAffiliationData");

  // unpack affiliation pairs
  if (buffer->size() > 0) UnpackAffiliationPairs(*buffer);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::DistributeAffiliationPairs()
{
  // relate all particles to all processors
  std::vector<int> particlestoproc(0);
  particleengineinterface_->RelateAllParticlesToAllProcs(particlestoproc);

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
  CommunicateSpecificAffiliationPairs(particletargets);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::CommunicateAffiliationPairs()
{
  // get reference to particles being communicated to target processors
  const std::vector<std::vector<int>>& particletargets =
      particleengineinterface_->GetCommunicatedParticleTargets();

  // communicate specific affiliation pairs
  CommunicateSpecificAffiliationPairs(particletargets);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::CommunicateSpecificAffiliationPairs(
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
      AddAffiliationPairToBuffer(sdata[torank], globalid, it->second);

      // erase affiliation pair
      affiliationdata_.erase(it);
    }
  }

  // communicate data via non-buffered send from proc to proc
  PARTICLEENGINE::COMMUNICATION::ImmediateRecvBlockingSend(comm_, sdata, rdata);

  // unpack affiliation pairs
  for (auto& p : rdata) UnpackAffiliationPairs(p.second);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::PackAllAffiliationPairs(
    std::vector<char>& buffer) const
{
  // iterate over all affiliation pairs
  for (auto& it : affiliationdata_)
  {
    // add affiliation pair to buffer
    AddAffiliationPairToBuffer(buffer, it.first, it.second);
  }
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::UnpackAffiliationPairs(
    const std::vector<char>& buffer)
{
  std::vector<char>::size_type position = 0;
  while (position < buffer.size())
  {
    // get affiliation pair
    const int globalid = DRT::ParObject::ExtractInt(position, buffer);
    const int rigidbody = DRT::ParObject::ExtractInt(position, buffer);

    // add affiliation pair
    affiliationdata_[globalid] = rigidbody;
  }
  if (position != buffer.size())
    dserror("mismatch in size of data %d <-> %d", static_cast<int>(buffer.size()), position);
}

void PARTICLERIGIDBODY::RigidBodyAffiliationPairs::AddAffiliationPairToBuffer(
    std::vector<char>& buffer, int globalid, int rigidbody) const
{
  DRT::PackBuffer data;
  data.StartPacking();

  // add affiliation pair
  data.AddtoPack(globalid);
  data.AddtoPack(rigidbody);

  // append packed affiliation pair to buffer
  buffer.insert(buffer.end(), data().begin(), data().end());
}
