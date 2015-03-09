/*----------------------------------------------------------------------*/
/*!
 \file scatra_mat_poro_ecm.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#include <vector>
#include "scatra_mat_poro_ecm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_comm/comm_utils.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::ScatraMatPoroECM::ScatraMatPoroECM(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: ScatraMat(matdata)
{
}

Teuchos::RCP<MAT::Material> MAT::PAR::ScatraMatPoroECM::CreateMaterial()
{
  return Teuchos::rcp(new MAT::ScatraMatPoroECM(this));
}


MAT::ScatraMatPoroECMType MAT::ScatraMatPoroECMType::instance_;

DRT::ParObject* MAT::ScatraMatPoroECMType::Create( const std::vector<char> & data )
{
  MAT::ScatraMatPoroECM* scatra_mat = new MAT::ScatraMatPoroECM();
  scatra_mat->Unpack(data);
  return scatra_mat;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatPoroECM::ScatraMatPoroECM()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::ScatraMatPoroECM::ScatraMatPoroECM(MAT::PAR::ScatraMatPoroECM* params)
  :   ScatraMat(params),
      params_(params)
{
}



