/*!----------------------------------------------------------------------
\file scatra_chemotaxis_mat.cpp

 \brief

This file contains the base material for chemotactic scalars.

<pre>
Maintainer: Moritz Thon
            thon@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-10364
</pre>
*----------------------------------------------------------------------*/

#include <vector>
#include "scatra_chemotaxis_mat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraChemotaxisMat::ScatraChemotaxisMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      numscal_(matdata->GetInt("NUMSCAL")),
      pair_(matdata->Get<std::vector<int>>("PAIR")),
      chemocoeff_(matdata->GetDouble("CHEMOCOEFF"))
{
  // Some checks for more safety
  if (numscal_ != (int)pair_->size())
    dserror("number of materials %d does not fit to size of material vector %d", numscal_,
        pair_->size());

  // is there exactly one '1' (i.e. attractant) and at least one '-1' (i.e. chemotractant)?
  int numpos = 0;
  int numneg = 0;
  for (int i = 0; i < numscal_; i++)
  {
    if (pair_->at(i) > 1e-10)
      numpos++;
    else if (pair_->at(i) < -1e-10)
      numneg++;
  }
  if (numpos != 1 or numneg != 1)
    dserror(
        "Each PAIR vector must contain exactly one '-1' (i.e. chemotractant) and exactly one '1' "
        "(i.e. attractant)!");

  return;
}


Teuchos::RCP<MAT::Material> MAT::PAR::ScatraChemotaxisMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraChemotaxisMat(this));
}


MAT::ScatraChemotaxisMatType MAT::ScatraChemotaxisMatType::instance_;


DRT::ParObject* MAT::ScatraChemotaxisMatType::Create(const std::vector<char>& data)
{
  MAT::ScatraChemotaxisMat* scatra_chemotaxis_mat = new MAT::ScatraChemotaxisMat();
  scatra_chemotaxis_mat->Unpack(data);
  return scatra_chemotaxis_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraChemotaxisMat::ScatraChemotaxisMat() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraChemotaxisMat::ScatraChemotaxisMat(MAT::PAR::ScatraChemotaxisMat* params)
    : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraChemotaxisMat::Pack(DRT::PackBuffer& data) const
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraChemotaxisMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror(
        "wrong instance type data. type = %d, UniqueParObjectId()=%d", type, UniqueParObjectId());
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
        params_ = static_cast<MAT::PAR::ScatraChemotaxisMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
