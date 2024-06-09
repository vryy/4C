/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of species and phases for electrochemistry applications

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "4C_mat_elchmat.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::ElchMat::ElchMat(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : Parameter(matdata),
      numdof_((matdata->Get<int>("NUMDOF"))),
      numscal_((matdata->Get<int>("NUMSCAL"))),
      numphase_(matdata->Get<int>("NUMPHASE")),
      phaseids_(matdata->Get<std::vector<int>>("PHASEIDS")),
      local_(matdata->Get<bool>("LOCAL"))
{
  if (numphase_ != (int)phaseids_.size())
    FOUR_C_THROW(
        "number of phases %d does not fit to size of phase vector %d", numphase_, phaseids_.size());

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator n;
    // phase
    for (n = phaseids_.begin(); n != phaseids_.end(); ++n)
    {
      const int phaseid = *n;
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(phaseid);
      mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(phaseid, mat));
    }
  }
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::ElchMat::create_material()
{
  return Teuchos::rcp(new Mat::ElchMat(this));
}


Mat::ElchMatType Mat::ElchMatType::instance_;


Core::Communication::ParObject* Mat::ElchMatType::Create(const std::vector<char>& data)
{
  Mat::ElchMat* elchmat = new Mat::ElchMat();
  elchmat->Unpack(data);
  return elchmat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ElchMat::ElchMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::ElchMat::ElchMat(Mat::PAR::ElchMat* params) : params_(params)
{
  // setup of material map
  if (params_->local_)
  {
    setup_mat_map();
  }
  // else: material rcps live inside Mat::PAR::MatList
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchMat::setup_mat_map()
{
  // safety first
  mat_.clear();
  if (not mat_.empty()) FOUR_C_THROW("What's going wrong here?");

  // make sure the referenced materials in material list have quick access parameters

  // here's the recursive creation of materials
  std::vector<int>::const_iterator n;
  for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); ++n)
  {
    const int phaseid = *n;
    Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(phaseid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(phaseid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchMat::clear()
{
  params_ = nullptr;
  mat_.clear();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchMat::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->Id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  if (params_ != nullptr and params_->local_)
  {
    // loop map of associated local materials
    std::vector<int>::const_iterator n;
    for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); n++)
      (mat_.find(*n))->second->Pack(data);
  }

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::ElchMat::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid(-1);
  extract_from_pack(position, data, matid);
  params_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<Mat::PAR::ElchMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (params_ != nullptr)  // params_ are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator n;
    for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); n++)
    {
      const int actphaseid = *n;
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(actphaseid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      mat_.insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(actphaseid, mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); n++)
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
