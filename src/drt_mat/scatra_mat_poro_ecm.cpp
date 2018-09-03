/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_poro_ecm.cpp

 \brief scatra material for transport within porous model with special implementations
        for ECM model

   \maintainer Andreas Rauch
               rauch@lnm.mw.tum.de
               http://www.lnm.mw.tum.de

\level 3
 *----------------------------------------------------------------------*/


#include <vector>
#include "scatra_mat_poro_ecm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatPoroECM::ScatraMatPoroECM(Teuchos::RCP<MAT::PAR::Material> matdata)
    : ScatraReactionMat(matdata), reacscale_(matdata->GetDouble("REACSCALE"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatPoroECM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatPoroECM(this));
}


MAT::ScatraMatPoroECMType MAT::ScatraMatPoroECMType::instance_;

DRT::ParObject* MAT::ScatraMatPoroECMType::Create(const std::vector<char>& data)
{
  MAT::ScatraMatPoroECM* scatra_mat = new MAT::ScatraMatPoroECM();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatPoroECM::ScatraMatPoroECM() : params_(NULL), reaccoeff_(0.0) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatPoroECM::ScatraMatPoroECM(MAT::PAR::ScatraMatPoroECM* params)
    : ScatraReactionMat(params), params_(params), reaccoeff_(0.0)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatPoroECM::Pack(DRT::PackBuffer& data) const
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

  // reaccoeff_
  AddtoPack(data, reaccoeff_);

  // add base class material
  ScatraReactionMat::Pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatPoroECM::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
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
        params_ = static_cast<MAT::PAR::ScatraMatPoroECM*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  // reaccoeff_
  ExtractfromPack(position, data, reaccoeff_);

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position, data, basedata);
  ScatraReactionMat::Unpack(basedata);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatPoroECM::ComputeReacCoeff(double chempot)
{
  reaccoeff_ = params_->reaccoeff_ * exp(params_->reacscale_ * chempot);
  return;
}
