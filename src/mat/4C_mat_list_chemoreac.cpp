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
MAT::PAR::MatListChemoReac::MatListChemoReac(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : MatList(matdata), MatListReactions(matdata), MatListChemotaxis(matdata)
{
}


Teuchos::RCP<CORE::MAT::Material> MAT::PAR::MatListChemoReac::create_material()
{
  return Teuchos::rcp(new MAT::MatListChemoReac(this));
}


MAT::MatListChemoReacType MAT::MatListChemoReacType::instance_;


CORE::COMM::ParObject* MAT::MatListChemoReacType::Create(const std::vector<char>& data)
{
  MAT::MatListChemoReac* MatListChemoReac = new MAT::MatListChemoReac();
  MatListChemoReac->Unpack(data);
  return MatListChemoReac;
}


/*----------------------------------------------------------------------*
 | construct empty material object                           thon 06/15 |
 *----------------------------------------------------------------------*/
MAT::MatListChemoReac::MatListChemoReac()
    : MatList(), MatListChemotaxis(), MatListReactions(), paramsreachemo_(nullptr)
{
}


/*----------------------------------------------------------------------*
 | construct the material object given material paramete     thon 06/15 |
 *----------------------------------------------------------------------*/
MAT::MatListChemoReac::MatListChemoReac(MAT::PAR::MatListChemoReac* params)
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
void MAT::MatListChemoReac::setup_mat_map()
{
  // We just have to add the chemotactic/reactive materials, since the rest is already done in
  // MAT::MatList::setup_mat_map() called from the MatList constructor

  // here's the recursive creation of materials
  MAT::MatListReactions::setup_mat_map();
  MAT::MatListChemotaxis::setup_mat_map();

  return;
}


/*----------------------------------------------------------------------*
 | reset everything                                          thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemoReac::clear()
{
  paramsreachemo_ = nullptr;
  return;
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemoReac::Pack(CORE::COMM::PackBuffer& data) const
{
  CORE::COMM::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (paramsreachemo_ != nullptr)
    matid = paramsreachemo_->Id();  // in case we are in post-process mode

  AddtoPack(data, matid);

  // Pack base class material
  MAT::MatListReactions::Pack(data);
  MAT::MatListChemotaxis::Pack(data);
}


/*----------------------------------------------------------------------*
 | Unpack data from a char vector into this class            thon 06/15 |
 *----------------------------------------------------------------------*/
void MAT::MatListChemoReac::Unpack(const std::vector<char>& data)
{
  // make sure we have a pristine material
  clear();

  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover paramsreac_
  int matid(-1);
  ExtractfromPack(position, data, matid);
  paramsreachemo_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
      {
        // Note: We need to do a dynamic_cast here since Chemotaxis, Reaction, and Chemo-reaction
        // are in a diamond inheritance structure
        paramsreachemo_ = dynamic_cast<MAT::PAR::MatListChemoReac*>(mat);
      }
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  MAT::MatList::ExtractfromPack(position, data, basedata);
  MAT::MatListReactions::Unpack(basedata);

  std::vector<char> basedata2(0);
  MAT::MatList::ExtractfromPack(position, data, basedata2);
  MAT::MatListChemotaxis::Unpack(basedata2);

  // in the postprocessing mode, we do not unpack everything we have packed
  // -> position check cannot be done in this case
  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
