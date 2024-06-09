/*----------------------------------------------------------------------*/
/*! \file
 \brief

This file contains the material for reactive AND chemotactic scalars. It is
in diamond inheritance with MatListReactions and MatListChemotaxis,
which govern the actual doings

\level 3
*----------------------------------------------------------------------*/


#include "4C_mat_list_chemoreac.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | standard constructor                                     thon 06/15 |
 *----------------------------------------------------------------------*/
Mat::PAR::MatListChemoReac::MatListChemoReac(Teuchos::RCP<Core::Mat::PAR::Material> matdata)
    : MatList(matdata), MatListReactions(matdata), MatListChemotaxis(matdata)
{
}


Teuchos::RCP<Core::Mat::Material> Mat::PAR::MatListChemoReac::create_material()
{
  return Teuchos::rcp(new Mat::MatListChemoReac(this));
}


Mat::MatListChemoReacType Mat::MatListChemoReacType::instance_;


Core::Communication::ParObject* Mat::MatListChemoReacType::Create(const std::vector<char>& data)
{
  Mat::MatListChemoReac* MatListChemoReac = new Mat::MatListChemoReac();
  MatListChemoReac->Unpack(data);
  return MatListChemoReac;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           thon 06/15 |
 *----------------------------------------------------------------------*/
Mat::MatListChemoReac::MatListChemoReac()
    : MatList(), MatListChemotaxis(), MatListReactions(), paramsreachemo_(nullptr)
{
}


/*----------------------------------------------------------------------*
 | construct the material object given material paramete     thon 06/15 |
 *----------------------------------------------------------------------*/
Mat::MatListChemoReac::MatListChemoReac(Mat::PAR::MatListChemoReac* params)
    : MatList(params), MatListChemotaxis(params), MatListReactions(params), paramsreachemo_(params)
{
  // setup of material map
  if (paramsreachemo_->local_)
  {
    setup_mat_map();
  }
}


/*----------------------------------------------------------------------*
 | setup of material map                                     thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemoReac::setup_mat_map()
{
  // We just have to add the chemotactic/reactive materials, since the rest is already done in
  // Mat::MatList::setup_mat_map() called from the MatList constructor

  // here's the recursive creation of materials
  Mat::MatListReactions::setup_mat_map();
  Mat::MatListChemotaxis::setup_mat_map();

  return;
}


/*----------------------------------------------------------------------*
 | reset everything                                          thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemoReac::clear()
{
  paramsreachemo_ = nullptr;
  return;
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemoReac::Pack(Core::Communication::PackBuffer& data) const
{
  Core::Communication::PackBuffer::SizeMarker sm(data);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (paramsreachemo_ != nullptr)
    matid = paramsreachemo_->Id();  // in case we are in post-process mode

  add_to_pack(data, matid);

  // Pack base class material
  Mat::MatListReactions::Pack(data);
  Mat::MatListChemotaxis::Pack(data);
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void Mat::MatListChemoReac::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  clear();

  std::vector<char>::size_type position = 0;

  Core::Communication::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover paramsreac_
  int matid(-1);
  extract_from_pack(position, data, matid);
  paramsreachemo_ = nullptr;
  if (Global::Problem::Instance()->Materials() != Teuchos::null)
    if (Global::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = Global::Problem::Instance()->Materials()->GetReadFromProblem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramsreachemo_ = dynamic_cast<Mat::PAR::MatListChemoReac*>(mat);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  Mat::MatList::extract_from_pack(position, data, basedata);
  Mat::MatListReactions::Unpack(basedata);

  std::vector<char> basedata2(0);
  Mat::MatList::extract_from_pack(position, data, basedata2);
  Mat::MatListChemotaxis::Unpack(basedata2);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
