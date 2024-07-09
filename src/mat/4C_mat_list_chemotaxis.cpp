/*----------------------------------------------------------------------*/
/*! \file
 \brief

This file contains the material for chemotactic scalars. It derives from MAT_matlist
and adds everything to supervise all the MAT_scatra_chemotaxis materials. The chemotaxation
itself is defined inside the MAT_scatra_chemotaxis materials. So MAT_matlist_chemotaxis
is just a "control instance".


\level 3
*----------------------------------------------------------------------*/


#include "4C_mat_list_chemotaxis.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | standard constructor                                      thon 06/15 |
 *----------------------------------------------------------------------*/
Mat::PAR::MatListChemotaxis::MatListChemotaxis(const Core::Mat::PAR::Parameter::Data& matdata)
    : MatList(matdata),
      numpair_((matdata.parameters.get<int>("NUMPAIR"))),
      pairids_((matdata.parameters.get<std::vector<int>>("PAIRIDS")))
{
  // check if sizes fit
  if (numpair_ != (int)pairids_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", nummat_,
        pairids_.size());

  if (numpair_ < 1)
    FOUR_C_THROW(
        "If you don't have chemotactic pairs, use MAT_matlist instead of MAT_matlist_chemotaxis!");

  if (not local_)
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;
    for (m = pairids_.begin(); m != pairids_.end(); ++m)
    {
      const int pairid = *m;
      Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(pairid);
      material_map_write()->insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(pairid, mat));
    }
  }
}

Teuchos::RCP<Core::Mat::Material> Mat::PAR::MatListChemotaxis::create_material()
{
  return Teuchos::rcp(new Mat::MatListChemotaxis(this));
}


Mat::MatListChemotaxisType Mat::MatListChemotaxisType::instance_;


Core::Communication::ParObject* Mat::MatListChemotaxisType::create(const std::vector<char>& data)
{
  Mat::MatListChemotaxis* MatListChemotaxis = new Mat::MatListChemotaxis();
  MatListChemotaxis->unpack(data);
  return MatListChemotaxis;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           thon 06/15 |
 *----------------------------------------------------------------------*/
Mat::MatListChemotaxis::MatListChemotaxis() : MatList(), paramschemo_(nullptr) {}


/*----------------------------------------------------------------------*
 | construct the material object given material paramete     thon 06/15 |
 *----------------------------------------------------------------------*/
Mat::MatListChemotaxis::MatListChemotaxis(Mat::PAR::MatListChemotaxis* params)
    : MatList(params), paramschemo_(params)
{
  // setup of material map
  if (paramschemo_->local_)
  {
    setup_mat_map();
  }
}


/*----------------------------------------------------------------------*
 | setup of material map                                     thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemotaxis::setup_mat_map()
{
  // We just have to add the chemotactic materials, since the rest is already done in
  // Mat::MatList::setup_mat_map() called from the MatList constructor

  // here's the recursive creation of materials
  std::vector<int>::const_iterator m;
  for (m = paramschemo_->pair_ids()->begin(); m != paramschemo_->pair_ids()->end(); ++m)
  {
    const int pairid = *m;
    Teuchos::RCP<Core::Mat::Material> mat = Mat::Factory(pairid);
    if (mat == Teuchos::null) FOUR_C_THROW("Failed to allocate this material");
    material_map_write()->insert(std::pair<int, Teuchos::RCP<Core::Mat::Material>>(pairid, mat));
  }
  return;
}


/*----------------------------------------------------------------------*
 | reset everything                                          thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemotaxis::clear()
{
  paramschemo_ = nullptr;
  return;
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemotaxis::pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (paramschemo_ != nullptr) matid = paramschemo_->id();  // in case we are in post-process mode

  add_to_pack(data, matid);

  // Pack base class material
  Mat::MatList::pack(data);
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemotaxis::unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, unique_par_object_id());

  // matid and recover paramsreac_
  int matid(-1);
  extract_from_pack(position, data, matid);
  paramschemo_ = nullptr;
  if (Global::Problem::instance()->materials() != Teuchos::null)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramschemo_ = dynamic_cast<Mat::PAR::MatListChemotaxis*>(mat);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  // extract base class material
  std::vector<char> basedata(0);
  Mat::MatList::extract_from_pack(position, data, basedata);
  Mat::MatList::unpack(basedata);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}


/*----------------------------------------------------------------------*
 | reaction ID by Index                                      thon 06/15 |
 *----------------------------------------------------------------------*/
int Mat::MatListChemotaxis::pair_id(const unsigned index) const
{
  if ((int)index < paramschemo_->numpair_)
    return paramschemo_->pairids_.at(index);
  else
  {
    FOUR_C_THROW("Index too large");
    return -1;
  }
}

FOUR_C_NAMESPACE_CLOSE
