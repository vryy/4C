/*!----------------------------------------------------------------------*/
/*!
\file elchmat.cpp

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/


#include <vector>
#include "elchmat.H"
#include "diffcond.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElchMat::ElchMat(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  current_((matdata->GetInt("CURRENT"))),
  numspec_(matdata->GetInt("NUMSPEC")),
  specids_(matdata->Get<std::vector<int> >("SPECIDS")),
  numphase_(matdata->GetInt("NUMPHASE")),
  phaseids_(matdata->Get<std::vector<int> >("PHASEIDS")),
  local_(false)
{
  // check if sizes fit
  if (numspec_ != (int)specids_->size())
    dserror("number of species %d does not fit to size of species vector %d", numspec_, specids_->size());

  if (numphase_ != (int)phaseids_->size())
    dserror("number of phases %d does not fit to size of phase vector %d", numphase_, phaseids_->size());

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    // species
    for (m=specids_->begin(); m!=specids_->end(); ++m)
    {
      const int specid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(specid);
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(specid,mat));
    }
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator n;
    // phase
    for (n=phaseids_->begin(); n!=phaseids_->end(); ++n)
    {
      const int phaseid = *n;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(phaseid);
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(phaseid,mat));
    }
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::ElchMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElchMat(this));
}


MAT::ElchMatType MAT::ElchMatType::instance_;


DRT::ParObject* MAT::ElchMatType::Create( const std::vector<char> & data )
{
  MAT::ElchMat* elchmat = new MAT::ElchMat();
  elchmat->Unpack(data);
  return elchmat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchMat::ElchMat()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ElchMat::ElchMat(MAT::PAR::ElchMat* params)
  : params_(params)
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
  std::vector<int>::const_iterator m;
  for (m=params_->SpecIds()->begin(); m!=params_->SpecIds()->end(); ++m)
  {
    const int specid = *m;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(specid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(specid,mat));
  }

  // here's the recursive creation of materials
  std::vector<int>::const_iterator n;
  for (n=params_->PhaseIds()->begin(); n!=params_->PhaseIds()->end(); ++n)
  {
    const int phaseid = *n;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(phaseid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(phaseid,mat));
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
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);

  if (params_->local_)
  {
    // loop map of associated local materials
    if (params_ != NULL)
    {
      //std::map<int, Teuchos::RCP<MAT::Material> >::const_iterator m;
      std::vector<int>::const_iterator m;
      for (m=params_->SpecIds()->begin(); m!=params_->SpecIds()->end(); m++)
      {
        (mat_.find(*m))->second->Pack(data);
      }

      //std::map<int, Teuchos::RCP<MAT::Material> >::const_iterator m;
      std::vector<int>::const_iterator n;
      for (n=params_->PhaseIds()->begin(); n!=params_->PhaseIds()->end(); n++)
      {
        (mat_.find(*n))->second->Pack(data);
      }
    }
  }
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
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid(-1);
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElchMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (params_ != NULL) // params_ are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m=params_->SpecIds()->begin(); m!=params_->SpecIds()->end(); m++)
    {
      const int actspecid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actspecid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(actspecid,mat));
    }

    std::vector<int>::const_iterator n;
    for (n=params_->PhaseIds()->begin(); n!=params_->PhaseIds()->end(); n++)
    {
      const int actphaseid = *n;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actphaseid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(actphaseid,mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (m=params_->SpecIds()->begin(); m!=params_->SpecIds()->end(); m++)
      {
        std::vector<char> pbtest;
        ExtractfromPack(position,data,pbtest);
        (mat_.find(*m))->second->Unpack(pbtest);
      }

      // loop map of associated local materials
      for (n=params_->PhaseIds()->begin(); n!=params_->PhaseIds()->end(); n++)
      {
        std::vector<char> pbtest;
        ExtractfromPack(position,data,pbtest);
        (mat_.find(*m))->second->Unpack(pbtest);
      }
    }
    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",data.size(),position);
  } // if (params_ != NULL)
}


