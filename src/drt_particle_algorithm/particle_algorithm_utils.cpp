/*---------------------------------------------------------------------------*/
/*! \file
\brief utils for particle algorithm
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_algorithm_utils.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
template <typename valtype>
void PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(const Teuchos::ParameterList& params,
    const std::string& name, std::map<PARTICLEENGINE::TypeEnum, valtype>& typetovalmap)
{
  // read from input file
  std::vector<std::string> typetoval;
  std::string word;
  std::istringstream typetovalstream(Teuchos::getNumericStringParameter(params, name));
  while (typetovalstream >> word) typetoval.push_back(word);

  // default case
  if (typetoval.size() == 1 and typetoval[0] == "none") return;

  // safety check
  if ((int)typetoval.size() % 2)
    dserror("input of '%s' (size = %d) relating particle type to value is odd!", name.c_str(),
        typetoval.size());

  std::string typestring;
  std::string valstring;
  for (int i = 0; i < (int)(typetoval.size() / 2); ++i)
  {
    typestring = typetoval[2 * i];
    valstring = typetoval[2 * i + 1];

    // get enum of particle types
    PARTICLEENGINE::TypeEnum particleType = PARTICLEENGINE::EnumFromTypeName(typestring);

    // get numeric value
    valtype val;
    try
    {
      // standard conversion (double to valtype) via assignment
      val = std::stod(valstring);
    }
    catch (const std::invalid_argument& e)
    {
      dserror("wrong format of numeric value provided following the name of the particle type!");
    }

    // insert into map
    auto iterator = typetovalmap.insert(std::make_pair(particleType, val));

    // safety check
    if (not iterator.second)
      dserror("failed inserting numeric value into map since key '%s' is already existing!",
          PARTICLEENGINE::EnumToTypeName(particleType).c_str());
  }
}

/*---------------------------------------------------------------------------*
 | template instantiations                                                   |
 *---------------------------------------------------------------------------*/
template void PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues<int>(
    const Teuchos::ParameterList& params, const std::string& name,
    std::map<PARTICLEENGINE::TypeEnum, int>& typetovalmap);

template void PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues<double>(
    const Teuchos::ParameterList& params, const std::string& name,
    std::map<PARTICLEENGINE::TypeEnum, double>& typetovalmap);
