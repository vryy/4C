/*----------------------------------------------------------------------*/
/*! \file
\brief Bundle holds all read-in materials of a #GLOBAL::Problem

\level 1

*/

/*----------------------------------------------------------------------*/
/* macros */


/*----------------------------------------------------------------------*/
/* headers */
#include "4C_mat_par_bundle.hpp"

#include "4C_mat_material_factory.hpp"
#include "4C_matelast_summand.hpp"
#include "4C_material_base.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
void MAT::PAR::Bundle::insert(int matid, CORE::UTILS::LazyPtr<CORE::MAT::PAR::Parameter> mat)
{
  matmap_.emplace(matid, mat);
}

/*----------------------------------------------------------------------*/
bool MAT::PAR::Bundle::id_exists(int id) const { return !(matmap_.find(id) == matmap_.end()); }

/*----------------------------------------------------------------------*/
CORE::MAT::PAR::Parameter* MAT::PAR::Bundle::ParameterById(const int num) const
{
  if (matmap_.size() == 0) FOUR_C_THROW("No materials available, num=%d", num);

  if (auto it = matmap_.find(num); it != matmap_.end())
    return it->second.get();
  else
    FOUR_C_THROW("Material 'MAT %d' could not be found", num);
}

/*----------------------------------------------------------------------*/
int MAT::PAR::Bundle::FirstIdByType(const CORE::Materials::MaterialType type) const
{
  if (auto it = std::find_if(matmap_.begin(), matmap_.end(),
          [&type](const auto& m) { return m.second->Type() == type; });
      it != matmap_.end())
    return it->first;
  else
    return -1;
}

FOUR_C_NAMESPACE_CLOSE
