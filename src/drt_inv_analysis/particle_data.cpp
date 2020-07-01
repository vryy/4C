/*----------------------------------------------------------------------*/
/*! \file
\brief Data for particles in the Sequential Monte Carlo algortihm

\level 3

*/
/*----------------------------------------------------------------------*/
#include "particle_data.H"

#include "Epetra_Map.h"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"

/*----------------------------------------------------------------------*/
INVANA::ParticleDataType INVANA::ParticleDataType::instance_;

/*----------------------------------------------------------------------*/
DRT::ParObject* INVANA::ParticleDataType::Create(const std::vector<char>& data)
{
  INVANA::ParticleData* my_particle = new INVANA::ParticleData();
  my_particle->Unpack(data);
  return my_particle;
}

/*----------------------------------------------------------------------*/
INVANA::ParticleData::ParticleData()
    : state_(Teuchos::null),
      posterior_(0.0),
      prior_(0.0),
      statechanged_(true),
      lcomm_(DRT::Problem::Instance()->GetNPGroup()->LocalComm())
{
}

/*----------------------------------------------------------------------*/
INVANA::ParticleData::ParticleData(ParticleData& data)
{
  state_ = Teuchos::rcp(new Epetra_Vector(data.GetState()));
  posterior_ = data.GetPosterior();
  prior_ = data.GetPrior();

  // since the above only works with valid state:
  statechanged_ = false;

  lcomm_ = DRT::Problem::Instance()->GetNPGroup()->LocalComm();
}

/*----------------------------------------------------------------------*/
INVANA::ParticleData& INVANA::ParticleData::operator=(ParticleData& data)
{
  state_->Scale(1.0, data.GetState());
  posterior_ = data.GetPosterior();
  prior_ = data.GetPrior();

  // since the above only works with valid state:
  statechanged_ = false;

  // the communicator doesnt need to be copied?

  return *this;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleData::Init(const Epetra_BlockMap& map)
{
  state_ = Teuchos::rcp(new Epetra_Vector(map, false));

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleData::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // ---- extract stuff from state_
  int llength = state_->MyLength();
  int glength = state_->GlobalLength();
  // extract the values
  std::vector<double> state(llength);
  state_->ExtractCopy(&state[0]);
  // extract the gids
  std::vector<int> gids(llength);
  state_->Map().MyGlobalElements(&gids[0]);
  // ---- end

  // pack stuff
  AddtoPack(data, glength);
  AddtoPack(data, llength);
  AddtoPack(data, state);
  AddtoPack(data, gids);
  AddtoPack(data, posterior_);
  AddtoPack(data, prior_);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleData::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  std::vector<double> state;
  std::vector<int> gids;

  // extract stuff
  int glength = ExtractInt(position, data);
  int llength = ExtractInt(position, data);
  state.reserve(llength);
  gids.resize(llength);
  ExtractfromPack(position, data, state);
  ExtractfromPack(position, data, gids);
  posterior_ = ExtractDouble(position, data);
  prior_ = ExtractDouble(position, data);
  statechanged_ = false;  // assume sending only with correctly evaluated state

  //  std::cout << "this procs gids: ";
  //  for (int i=0; i<llength; i++)
  //    std::cout << gids[i] << " ";
  //  std::cout << std::endl;
  //  std::cout << "GLENGHT: " << glength << std::endl;
  //  std::cout << "PROC: " << DRT::Problem::Instance()->GetNPGroup()->GlobalComm()->MyPID() <<
  //  "arrives here" << std::endl;
  // ---- reconstruct state_
  Epetra_Map amap((int)glength, llength, &gids[0], 0, *lcomm_);
  state_ = Teuchos::rcp(new Epetra_Vector(amap, false));
  state_->ReplaceGlobalValues(llength, &state[0], &gids[0]);
  // ---- end

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleData::SetState(const Epetra_Vector& state)
{
  state_->Scale(1.0, state);
  statechanged_ = true;

  return;
}

/*----------------------------------------------------------------------*/
void INVANA::ParticleData::SetData(double& posterior, double& prior)
{
  posterior_ = posterior;
  prior_ = prior;
  statechanged_ = false;

  return;
}
