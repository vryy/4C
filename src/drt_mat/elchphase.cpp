/*----------------------------------------------------------------------------*/
/*! \file
\brief material that stores a list of ion species in electrolyte solutions

\level 2


*/
/*----------------------------------------------------------------------------*/

#include <vector>
#include "elchphase.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElchPhase::ElchPhase(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      epsilon_(matdata->GetDouble("EPSILON")),
      tortuosity_(matdata->GetDouble("TORTUOSITY")),
      nummat_(matdata->GetInt("NUMMAT")),
      matids_(*matdata->Get<std::vector<int>>("MATIDS")),
      local_(matdata->GetInt("LOCAL"))
{
  if (nummat_ != (int)matids_.size())
    dserror("number of phases %d does not fit to size of phase vector %d", nummat_, matids_.size());

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


DRT::ParObject* MAT::ElchPhaseType::Create(const std::vector<char>& data)
{
  MAT::ElchPhase* elchphase = new MAT::ElchPhase();
  elchphase->Unpack(data);
  return elchphase;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchPhase::ElchPhase() : params_(NULL) {}


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
  if (not mat_.empty()) dserror("What's going wrong here?");

  // make sure the referenced materials in material list have quick access parameters

  // here's the recursive creation of materials
  std::vector<int>::const_iterator n;
  for (n = params_->MatIds().begin(); n != params_->MatIds().end(); ++n)
  {
    const int matid = *n;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    mat_.insert(std::pair<int, Teuchos::RCP<MAT::Material>>(matid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::Clear()
{
  params_ = NULL;
  mat_.clear();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ElchPhase::Pack(DRT::PackBuffer& data) const
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

  if (params_->local_)
  {
    // loop map of associated local materials
    if (params_ != NULL)
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
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElchPhase*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (params_ != NULL)  // params_ are not accessible in postprocessing mode
  {
    std::vector<int>::const_iterator n;
    for (n = params_->MatIds().begin(); n != params_->MatIds().end(); n++)
    {
      const int actmatid = *n;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actmatid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
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
      dserror("Mismatch in size of data %d <-> %d", data.size(), position);
  }  // if (params_ != NULL)
}
