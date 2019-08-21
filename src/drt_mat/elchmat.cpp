/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of species and phases for electrochemistry applications

\level 2

\maintainer Christoph Schmidt

*/
/*----------------------------------------------------------------------------*/

#include <vector>
#include "elchmat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElchMat::ElchMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      numdof_((matdata->GetInt("NUMDOF"))),
      numscal_((matdata->GetInt("NUMSCAL"))),
      numphase_(matdata->GetInt("NUMPHASE")),
      phaseids_(*matdata->Get<std::vector<int>>("PHASEIDS")),
      local_(matdata->GetInt("LOCAL"))
{
  if (numphase_ != (int)phaseids_.size())
    dserror(
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


DRT::ParObject* MAT::ElchMatType::Create(const std::vector<char>& data)
{
  MAT::ElchMat* elchmat = new MAT::ElchMat();
  elchmat->Unpack(data);
  return elchmat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchMat::ElchMat() : params_(NULL) {}


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
  if (not mat_.empty()) dserror("What's going wrong here?");

  // make sure the referenced materials in material list have quick access parameters

  // here's the recursive creation of materials
  std::vector<int>::const_iterator n;
  for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); ++n)
  {
    const int phaseid = *n;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(phaseid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(phaseid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchMat::Clear()
{
  params_ = NULL;
  mat_.clear();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchMat::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data, matid);

  if (params_ != NULL and params_->local_)
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElchMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (params_ != NULL)  // params_ are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator n;
    for (n = params_->PhaseIds().begin(); n != params_->PhaseIds().end(); n++)
    {
      const int actphaseid = *n;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actphaseid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
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
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }  // if (params_ != NULL)
}
