/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of species and phases for electrochemistry applications

\level 2


*/
/*----------------------------------------------------------------------------*/

#include "baci_mat_elchmat.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElchMat::ElchMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      numdof_((*matdata->Get<int>("NUMDOF"))),
      numscal_((*matdata->Get<int>("NUMSCAL"))),
      numphase_(*matdata->Get<int>("NUMPHASE")),
      phaseids_(*matdata->Get<std::vector<int>>("PHASEIDS")),
      local_(*matdata->Get<bool>("LOCAL"))
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
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(phaseid);
      mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(phaseid, mat));
    }
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::ElchMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElchMat(this));
}


MAT::ElchMatType MAT::ElchMatType::instance_;


CORE::COMM::ParObject* MAT::ElchMatType::Create(const std::vector<char>& data)
{
  MAT::ElchMat* elchmat = new MAT::ElchMat();
  elchmat->Unpack(data);
  return elchmat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchMat::ElchMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchMat::ElchMat(MAT::PAR::ElchMat* params) : params_(params)
{
  // setup of material map
  if (params_->local_)
  {
    SetupMatMap();
  }
  // else: material rcps live inside MAT::PAR::MatList
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchMat::SetupMatMap()
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
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(phaseid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(phaseid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchMat::Clear()
{
  params_ = nullptr;
  mat_.clear();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchMat::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ElchMat::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElchMat*>(mat);
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
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actphaseid);
      if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
      mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(actphaseid, mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); n++)
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
