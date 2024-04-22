/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of ion species in electrolyte solutions

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "baci_mat_elchphase.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElchPhase::ElchPhase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      epsilon_(*matdata->Get<double>("EPSILON")),
      tortuosity_(*matdata->Get<double>("TORTUOSITY")),
      nummat_(*matdata->Get<int>("NUMMAT")),
      matids_(*matdata->Get<std::vector<int>>("MATIDS")),
      local_(*matdata->Get<bool>("LOCAL"))
{
  if (nummat_ != (int)matids_.size())
    FOUR_C_THROW(
        "number of phases %d does not fit to size of phase vector %d", nummat_, matids_.size());

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator n;
    // phase
    for (n = matids_.begin(); n != matids_.end(); ++n)
    {
      const int matid = *n;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
      mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(matid, mat));
    }
  }
}


Teuchos::RCP<MAT::Material> MAT::PAR::ElchPhase::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElchPhase(this));
}

MAT::ElchPhaseType MAT::ElchPhaseType::instance_;


CORE::COMM::ParObject* MAT::ElchPhaseType::Create(const std::vector<char>& data)
{
  MAT::ElchPhase* elchphase = new MAT::ElchPhase();
  elchphase->Unpack(data);
  return elchphase;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchPhase::ElchPhase() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchPhase::ElchPhase(MAT::PAR::ElchPhase* params) : params_(params)
{
  // setup of material map
  if (params_->local_)
  {
    SetupMatMap();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::SetupMatMap()
{
  // safety first
  mat_.clear();
  if (not mat_.empty()) FOUR_C_THROW("What's going wrong here?");

  // make sure the referenced materials in material list have quick access parameters

  // here's the recursive creation of materials
  std::vector<int>::const_iterator n;
  for (n = params_->MatIds().begin(); n != params_->MatIds().end(); ++n)
  {
    const int matid = *n;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(matid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::Clear()
{
  params_ = nullptr;
  mat_.clear();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  if (params_->local_)
  {
    // loop map of associated local materials
    if (params_ != nullptr)
    {
      // std::map<int, Teuchos::RCP<MAT::Material> >::const_iterator m;
      std::vector<int>::const_iterator n;
      for (n = params_->MatIds().begin(); n != params_->MatIds().end(); n++)
      {
        (mat_.find(*n))->second->Pack(data);
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElchPhase*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (params_ != nullptr)  // params_ are not accessible in postprocessing mode
  {
    std::vector<int>::const_iterator n;
    for (n = params_->MatIds().begin(); n != params_->MatIds().end(); n++)
    {
      const int actmatid = *n;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actmatid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(actmatid, mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (n = params_->MatIds().begin(); n != params_->MatIds().end(); n++)
      {
        std::vector<char> pbtest;
        ExtractfromPack(position, data, pbtest);
        (mat_.find(*n))->second->Unpack(pbtest);
      }
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
  }  // if (params_ != nullptr)
}

FOUR_C_NAMESPACE_CLOSE
