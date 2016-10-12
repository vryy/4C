/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_multiporo.cpp

 \brief scatra material for transport within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/




#include <vector>
#include "scatra_mat_multiporo.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatMultiPoro::ScatraMatMultiPoro(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: ScatraMat(matdata),
  density_(matdata->GetDouble("DENSITY")),
  phaseID_(matdata->GetInt("PHASEID"))
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatMultiPoro::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatMultiPoro(this));
}


MAT::ScatraMatMultiPoroType MAT::ScatraMatMultiPoroType::instance_;

DRT::ParObject* MAT::ScatraMatMultiPoroType::Create( const std::vector<char> & data )
{
  MAT::ScatraMatMultiPoro* scatra_mat = new MAT::ScatraMatMultiPoro();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoro::ScatraMatMultiPoro()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatMultiPoro::ScatraMatMultiPoro(MAT::PAR::ScatraMatMultiPoro* params)
  :   ScatraMat(params),
      params_(params)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoro::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // matid
  int matid = -1;
  if (params_ != NULL)
    matid = params_->Id(); // in case we are in post-process mode
  AddtoPack(data, matid);

  // add base class material
  ScatraMat::Pack(data);

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::ScatraMatMultiPoro::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId())
    dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::ScatraMatMultiPoro*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  // extract base class material
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ScatraMat::Unpack(basedata);

}


