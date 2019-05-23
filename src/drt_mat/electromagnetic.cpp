/*!----------------------------------------------------------------------
\brief Contains conductivity, permittivity and permeability of the medium for isotropic
       electromagetic field evolution.
       MAT 1 MAT_Electromagnetic CONDUCTIVITY 0.0 PERMITTIVITY 1.732 PERMEABILITY 1.732


\level 2

\maintainer Luca Berardocco

*/
/*----------------------------------------------------------------------*/

#include <vector>
#include "electromagnetic.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ElectromagneticMat::ElectromagneticMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(sigma_)->PutScalar(matdata->GetDouble("CONDUCTIVITY"));
  matparams_.at(epsilon_)->PutScalar(matdata->GetDouble("PERMITTIVITY"));
  matparams_.at(mu_)->PutScalar(matdata->GetDouble("PERMEABILITY"));

  return;
}

Teuchos::RCP<MAT::Material> MAT::PAR::ElectromagneticMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ElectromagneticMat(this));
}

MAT::ElectromagneticMatType MAT::ElectromagneticMatType::instance_;

void MAT::PAR::ElectromagneticMat::OptParams(std::map<std::string, int> *pnames)
{
  pnames->insert(std::pair<std::string, int>("CONDUCTIVITY", sigma_));
  pnames->insert(std::pair<std::string, int>("PERMITTIVITY", epsilon_));
  pnames->insert(std::pair<std::string, int>("PERMEABILITY", mu_));
}

DRT::ParObject *MAT::ElectromagneticMatType::Create(const std::vector<char> &data)
{
  MAT::ElectromagneticMat *soundprop = new MAT::ElectromagneticMat();
  soundprop->Unpack(data);
  return soundprop;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ElectromagneticMat::ElectromagneticMat() : params_(NULL) {}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::ElectromagneticMat::ElectromagneticMat(MAT::PAR::ElectromagneticMat *params) : params_(params)
{
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::ElectromagneticMat::Pack(DRT::PackBuffer &data) const
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

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::ElectromagneticMat::Unpack(const std::vector<char> &data)
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
      MAT::PAR::Parameter *mat =
          DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ElectromagneticMat *>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
