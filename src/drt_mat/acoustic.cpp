/*!----------------------------------------------------------------------
\file acoustic.cpp
\brief contains a density and a speed of sound to characterize all
       necessary properties for isotropic sound propagation in fluid
       like material
       example input line:
       MAT 1
<pre>
\level 2
\maintainer Luca Berardocco
            berardoccoo@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/

#include <vector>
#include "acoustic.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::AcousticMat::AcousticMat(Teuchos::RCP<MAT::PAR::Material> matdata) : Parameter(matdata)
{
  Epetra_Map dummy_map(1, 1, 0, *(DRT::Problem::Instance()->GetNPGroup()->LocalComm()));
  for (int i = first; i <= last; i++)
  {
    matparams_.push_back(Teuchos::rcp(new Epetra_Vector(dummy_map, true)));
  }
  matparams_.at(density)->PutScalar(matdata->GetDouble("DENSITY"));
  matparams_.at(c)->PutScalar(matdata->GetDouble("C"));

  return;
}


Teuchos::RCP<MAT::Material> MAT::PAR::AcousticMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AcousticMat(this));
}

MAT::AcousticMatType MAT::AcousticMatType::instance_;

void MAT::PAR::AcousticMat::OptParams(std::map<std::string, int>* pnames)
{
  pnames->insert(std::pair<std::string, int>("DENSITY", density));
  pnames->insert(std::pair<std::string, int>("C", c));
}

DRT::ParObject* MAT::AcousticMatType::Create(const std::vector<char>& data)
{
  MAT::AcousticMat* soundprop = new MAT::AcousticMat();
  soundprop->Unpack(data);
  return soundprop;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::AcousticMat::AcousticMat() : params_(NULL) {}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::AcousticMat::AcousticMat(MAT::PAR::AcousticMat* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::AcousticMat::Pack(DRT::PackBuffer& data) const
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
void MAT::AcousticMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::AcousticMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
