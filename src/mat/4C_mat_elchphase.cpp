/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of ion species in electrolyte solutions

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_elchphase.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ElchPhase::ElchPhase(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      epsilon_(matdata->Get<double>("EPSILON")),
      tortuosity_(matdata->Get<double>("TORTUOSITY")),
      nummat_(matdata->Get<int>("NUMMAT")),
      matids_(matdata->Get<std::vector<int>>("MATIDS")),
      local_(matdata->Get<bool>("LOCAL"))
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
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(matid);
      mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(matid, mat));
    }
  }
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::ElchPhase::create_material()
{
  return Teuchos::rcp(new Mat::ElchPhase(this));
}

Mat::ElchPhaseType Mat::ElchPhaseType::instance_;


Core::Communication::ParObject* Mat::ElchPhaseType::Create(const std::vector<char>& data)
{
  Mat::ElchPhase* elchphase = new Mat::ElchPhase();
  elchphase->Unpack(data);
  return elchphase;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ElchPhase::ElchPhase() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ElchPhase::ElchPhase(Mat::PAR::ElchPhase* params) : params_(params)
{
  // setup of material map
  if (params_->local_)
  {
    setup_mat_map();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchPhase::setup_mat_map()
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
    Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(matid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(matid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchPhase::clear()
{
  params_ = nullptr;
  mat_.clear();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchPhase::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  if (params_->local_)
  {
    // loop map of associated local materials
    if (params_ != nullptr)
    {
      // std::map<int, Teuchos::RCP<Core::Mat::Material> >::const_iterator m;
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
void Mat::ElchPhase::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ElchPhase*>(mat);
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
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(actmatid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(actmatid, mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (n = params_->MatIds().begin(); n != params_->MatIds().end(); n++)
      {
        std::vector<char> pbtest;
        extract_from_pack(position, data, pbtest);
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
