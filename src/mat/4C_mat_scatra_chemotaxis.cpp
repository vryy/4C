/*----------------------------------------------------------------------*/
/*! \file
 \brief
This file contains the base material for chemotactic scalars.

\level 3

*----------------------------------------------------------------------*/

#include "4C_mat_scatra_chemotaxis.hpp"

#include "4C_comm_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraChemotaxisMat::ScatraChemotaxisMat(Teuchos::RCP<CORE::MAT::PAR::Material> matdata)
    : Parameter(matdata),
      numscal_(matdata->Get<int>("NUMSCAL")),
      pair_(matdata->Get<std::vector<int>>("PAIR")),
      chemocoeff_(matdata->Get<double>("CHEMOCOEFF"))
{
  // Some checks for more safety
  if (numscal_ != (int)pair_.size())
    FOUR_C_THROW("number of materials %d does not fit to size of material vector %d", numscal_,
        pair_.size());

  // is there exactly one '1' (i.e. attractant) and at least one '-1' (i.e. chemotractant)?
  int numpos = 0;
  int numneg = 0;
  for (int i = 0; i < numscal_; i++)
  {
    if (pair_.at(i) > 1e-10)
      numpos++;
    else if (pair_.at(i) < -1e-10)
      numneg++;
  }
  if (numpos != 1 or numneg != 1)
    FOUR_C_THROW(
        "Each PAIR vector must contain exactly one '-1' (i.e. chemotractant) and exactly one '1' "
        "(i.e. attractant)!");

  return;
}


Teuchos::RCP<CORE::MAT::Material> MAT::PAR::ScatraChemotaxisMat::create_material()
{
  return Teuchos::rcp(new MAT::ScatraChemotaxisMat(this));
}


MAT::ScatraChemotaxisMatType MAT::ScatraChemotaxisMatType::instance_;


CORE::COMM::ParObject* MAT::ScatraChemotaxisMatType::Create(const std::vector<char>& data)
{
  MAT::ScatraChemotaxisMat* scatra_chemotaxis_mat = new MAT::ScatraChemotaxisMat();
  scatra_chemotaxis_mat->Unpack(data);
  return scatra_chemotaxis_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraChemotaxisMat::ScatraChemotaxisMat() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraChemotaxisMat::ScatraChemotaxisMat(MAT::PAR::ScatraChemotaxisMat* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraChemotaxisMat::Pack(CORE::COMM::PackBuffer& data) const
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraChemotaxisMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
  if (GLOBAL::Problem::Instance()->Materials() != Teuchos::null)
    if (GLOBAL::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = GLOBAL::Problem::Instance()->Materials()->GetReadFromProblem();
      CORE::MAT::PAR::Parameter* mat =
          GLOBAL::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraChemotaxisMat*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size())
    FOUR_C_THROW("Mismatch in size of data %d <-> %d", data.size(), position);
}

FOUR_C_NAMESPACE_CLOSE
