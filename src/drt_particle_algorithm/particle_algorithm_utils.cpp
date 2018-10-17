/*---------------------------------------------------------------------------*/
/*!
\file particle_algorithm_utils.cpp

\brief utils for particle algorithm

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_algorithm_utils.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | read parameters relating particle types to IDs             sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToIDs(const Teuchos::ParameterList& params,
    const std::string& name, std::map<PARTICLEENGINE::TypeEnum, int>& typetoidmap)
{
  // read from input file
  std::vector<std::string> typetoid;
  std::string word;
  std::istringstream typetoidstream(Teuchos::getNumericStringParameter(params, name));
  while (typetoidstream >> word) typetoid.push_back(word);

  // safety check
  if ((int)typetoid.size() % 2)
    dserror("input of '%s' (vector size = %d) relating particle type to ID is odd!", name.c_str(),
        typetoid.size());

  int numoftypestoids = typetoid.size() / 2;

  std::string typestring;
  std::string idstring;
  for (int i = 0; i < numoftypestoids; ++i)
  {
    typestring = typetoid[2 * i];
    idstring = typetoid[2 * i + 1];

    // get enum of particle types
    PARTICLEENGINE::TypeEnum particleType = PARTICLEENGINE::EnumFromTypeName(typestring);

    // get function id
    int id = -1;
    try
    {
      id = std::stoi(idstring);
    }
    catch (const std::invalid_argument& e)
    {
      dserror(
          "wrong format provided: expecting an integer of corresponding ID following the name of "
          "the particle type!");
    }

    // safety check
    if (id < 0) dserror("negative function ID = %d not possible!", id);

    // insert into map
    auto iterator = typetoidmap.insert(std::make_pair(particleType, id));

    // safety check
    if (not iterator.second)
      dserror("failed inserting pair ('%s', '%d') into map since key already existing!",
          PARTICLEENGINE::EnumToTypeName(particleType).c_str(), id);
  }
}
