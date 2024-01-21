/*----------------------------------------------------------------------*/
/*! \file
\brief scatra_mat_aniso.cpp

\level 3

*----------------------------------------------------------------------*/


#include "baci_mat_scatra_mat_aniso.H"

#include "baci_comm_utils.H"
#include "baci_global_data.H"
#include "baci_mat_par_bundle.H"

#include <vector>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatAniso::ScatraMatAniso(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(DRT::Problem::Instance()->GetCommunicators()->LocalComm()));
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

CORE::COMM::ParObject* MAT::ScatraMatAnisoType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatAniso* scatra_mat_aniso = new MAT::ScatraMatAniso();
  scatra_mat_aniso->Unpack(data);
  return scatra_mat_aniso;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatAniso::ScatraMatAniso() : params_(nullptr) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatAniso::ScatraMatAniso(MAT::PAR::ScatraMatAniso* params) : params_(params) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatAniso::Pack(CORE::COMM::PackBuffer& data) const
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
void MAT::ScatraMatAniso::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;

  CORE::COMM::ExtractAndAssertId(position, data, UniqueParObjectId());

  // matid and recover params_
  int matid;
  ExtractfromPack(position, data, matid);
  params_ = nullptr;
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

BACI_NAMESPACE_CLOSE
