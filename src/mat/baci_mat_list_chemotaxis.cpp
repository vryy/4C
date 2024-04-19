/*----------------------------------------------------------------------*/
/*! \file
 \brief

This file contains the material for chemotactic scalars. It derives from MAT_matlist
and adds everything to supervise all the MAT_scatra_chemotaxis materials. The chemotaxation
itself is defined inside the MAT_scatra_chemotaxis materials. So MAT_matlist_chemotaxis
is just a "control instance".


\level 3
*----------------------------------------------------------------------*/


#include "baci_mat_list_chemotaxis.hpp"

#include "baci_global_data.hpp"
#include "baci_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | standard constructor                                      thon 06/15 |
 *----------------------------------------------------------------------*/
MAT::PAR::MatListChemotaxis::MatListChemotaxis(Teuchos::RCP<MAT::PAR::Material> matdata)
    : MatList(matdata),
      numpair_((*matdata->Get<int>("NUMPAIR"))),
      pairids_((matdata->Get<std::vector<int>>("PAIRIDS")))
{
  // check if sizes fit
  if (numpair_ != (int)pairids_->size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", nummat_,
        pairids_->size());

  if (numpair_ < 1)
    FOUR_C_THROW(
        "If you don't have chemotactic pairs, use MAT_matlist instead of MAT_matlist_chemotaxis!");

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = pairids_->begin(); m != pairids_->end(); ++m)
    {
      const int pairid = *m;
      Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(pairid);
      MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<MAT::Material>>(pairid, mat));
    }
  }
}

Teuchos::RCP<MAT::Material> MAT::PAR::MatListChemotaxis::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MatListChemotaxis(this));
}


MAT::MatListChemotaxisType MAT::MatListChemotaxisType::instance_;


CORE::COMM::ParObject* MAT::MatListChemotaxisType::Create(const std::vector<char>& data)
{
  MAT::MatListChemotaxis* MatListChemotaxis = new MAT::MatListChemotaxis();
  MatListChemotaxis->Unpack(data);
  return MatListChemotaxis;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           thon 06/15 |
 *----------------------------------------------------------------------*/
MAT::MatListChemotaxis::MatListChemotaxis() : MatList(), paramschemo_(nullptr) {}


/*----------------------------------------------------------------------*
 | construct the material object given material paramete     thon 06/15 |
 *----------------------------------------------------------------------*/
MAT::MatListChemotaxis::MatListChemotaxis(MAT::PAR::MatListChemotaxis* params)
    : MatList(params), paramschemo_(params)
{
  // setup of material map
  if (paramschemo_->local_)
  {
    SetupMatMap();
  }
}


/*----------------------------------------------------------------------*
 | setup of material map                                     thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemotaxis::SetupMatMap()
{
  // We just have to add the chemotactic materials, since the rest is already done in
  // MAT::MatList::SetupMatMap() called from the MatList constructor

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m = paramschemo_->PairIds()->begin(); m != paramschemo_->PairIds()->end(); ++m)
  {
    const int pairid = *m;
    Teuchos::RCP<MAT::Material> mat = MAT::Material::Factory(pairid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    MaterialMapWrite()->insert(std::pair<int, Teuchos::RCP<MAT::Material>>(pairid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*
 | reset everything                                          thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemotaxis::Clear()
{
  paramschemo_ = nullptr;
  return;
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemotaxis::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (paramschemo_ != nullptr) matid = paramschemo_->Id();  // in case we are in post-process mode

  AddtoPack(data, matid);

  // Pack base class material
  MAT::MatList::Pack(data);
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemotaxis::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  Clear();

  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover paramsreac_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  paramschemo_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramschemo_ = dynamic_cast<MAT::PAR::MatListChemotaxis*>(mat);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position, data, basedata);
  MAT::MatList::Unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | reaction ID by Index                                      thon 06/15 |
 *----------------------------------------------------------------------*/
int MAT::MatListChemotaxis::PairID(const unsigned index) const
{
  if ((int)index < paramschemo_->numpair_)
    return paramschemo_->pairids_->at(index);
  else
  {
    FOUR_C_THROW("Index too large");
    return -1;
  }
}

FOUR_C_NAMESPACE_CLOSE
