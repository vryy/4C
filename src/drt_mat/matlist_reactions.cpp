/*!----------------------------------------------------------------------
\file matlist_reactions.cpp

<pre>
Maintainer: Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364
</pre>
*----------------------------------------------------------------------*/


#include <vector>
#include "matlist_reactions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::MatListReactions::MatListReactions(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  local_((matdata->GetInt("LOCAL"))),
  nummat_(matdata->GetInt("NUMMAT")),
  matids_(matdata->Get<std::vector<int> >("MATIDS")),
  numreac_((matdata->GetInt("NUMREAC"))),
  reacids_((matdata->Get<std::vector<int> >("REACIDS")))
{
  // check if sizes fit
  if (nummat_ != (int)matids_->size())
    dserror("number of materials %d does not fit to size of material vector %d", nummat_, matids_->size());

  if (numreac_ != (int)reacids_->size())
      dserror("number of materials %d does not fit to size of material vector %d", nummat_, reacids_->size());

    if (numreac_< 1)
      dserror("if you don't have reactions, use MAT_matlist instead of MAT_matlist_reactions!");

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m=matids_->begin(); m!=matids_->end(); ++m)
    {
      const int matid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(matid,mat));
    }

    for (m=reacids_->begin(); m!=reacids_->end(); ++m)
    {
      const int reacid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(reacid);
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(reacid,mat));
    }
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::MatListReactions::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MatListReactions(this));
}


MAT::MatListReactionsType MAT::MatListReactionsType::instance_;


DRT::ParObject* MAT::MatListReactionsType::Create( const std::vector<char> & data )
{
  MAT::MatListReactions* MatListReactions = new MAT::MatListReactions();
  MatListReactions->Unpack(data);
  return MatListReactions;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MatListReactions::MatListReactions()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MatListReactions::MatListReactions(MAT::PAR::MatListReactions* params)
  : params_(params)
{
  // setup of material map
  if (params_->local_)
  {
    SetupMatMap();
  }
  // else: material Teuchos::rcps live inside MAT::PAR::MatListReactions
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MatListReactions::SetupMatMap()
{
  // safety first
  mat_.clear();
  if (not mat_.empty()) dserror("What's going wrong here?");

  // make sure the referenced materials in material list have quick access parameters

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m=params_->MatIds()->begin(); m!=params_->MatIds()->end(); ++m)
  {
    const int matid = *m;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(matid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(matid,mat));
  }

  for (m=params_->ReacIds()->begin(); m!=params_->ReacIds()->end(); ++m)
  {
    const int reacid = *m;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(reacid);
    if (mat == Teuchos::null) dserror("Failed to allocate this material");
    mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(reacid,mat));
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MatListReactions::Clear()
{
  params_ = NULL;
  mat_.clear();
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MatListReactions::Pack(DRT::PackBuffer& data) const
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
      std::vector<int>::const_iterator m;
      for (m=params_->MatIds()->begin(); m!=params_->MatIds()->end(); m++)
      {
        (mat_.find(*m))->second->Pack(data);
      }

      for (m=params_->ReacIds()->begin(); m!=params_->ReacIds()->end(); m++)
      {
        (mat_.find(*m))->second->Pack(data);
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MatListReactions::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::MatListReactions*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (params_ != NULL) // params_ are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m=params_->MatIds()->begin(); m!=params_->MatIds()->end(); m++)
    {
      const int actmatid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actmatid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(actmatid,mat));
    }

    for (m=params_->ReacIds()->begin(); m!=params_->ReacIds()->end(); m++)
    {
      const int actmatid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(actmatid);
      if (mat == Teuchos::null) dserror("Failed to allocate this material");
      mat_.insert(std::pair<int,Teuchos::RCP<MAT::Material> >(actmatid,mat));
    }

    if (params_->local_)
    {
      // loop map of associated local materials
      for (m=params_->MatIds()->begin(); m!=params_->MatIds()->end(); m++)
      {
        std::vector<char> pbtest;
        ExtractfromPack(position,data,pbtest);
        (mat_.find(*m))->second->Unpack(pbtest);
      }

      for (m=params_->ReacIds()->begin(); m!=params_->ReacIds()->end(); m++)
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
