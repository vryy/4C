/*----------------------------------------------------------------------*/
/*! \file
\brief scatra_mat_aniso.cpp

\level 3

\maintainer Christoph Schmidt
*----------------------------------------------------------------------*/


#include <vector>
#include "scatra_mat_aniso.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatAniso::ScatraMatAniso(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(diff1)->PutScalar(matdata->GetDouble("DIFF1"));
  matparams_.at(diff2)->PutScalar(matdata->GetDouble("DIFF2"));
  matparams_.at(diff3)->PutScalar(matdata->GetDouble("DIFF3"));
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatAniso::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatAniso(this));
}


MAT::ScatraMatAnisoType MAT::ScatraMatAnisoType::instance_;

void MAT::PAR::ScatraMatAniso::OptParams(std::map<std::string, int>* pnames) {}

DRT::ParObject* MAT::ScatraMatAnisoType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatAniso* scatra_mat_aniso = new MAT::ScatraMatAniso();
  scatra_mat_aniso->Unpack(data);
  return scatra_mat_aniso;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatAniso::ScatraMatAniso() : params_(NULL) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatAniso::ScatraMatAniso(MAT::PAR::ScatraMatAniso* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatAniso::Pack(DRT::PackBuffer& data) const
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
void MAT::ScatraMatAniso::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::ScatraMatAniso*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
