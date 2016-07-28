/*!----------------------------------------------------------------------
\file meshFreeAMmat.cpp

\brief material specifications of the meshfree method implemented by katta

\level 2

\maintainer Alessandro Cattabiani
*----------------------------------------------------------------------*/


#include <vector>
#include "meshFreeAMmat.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::MeshFreeAMMat::MeshFreeAMMat(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  initMass_(matdata->GetDouble("INITMASS")),
  initRadius_(matdata->GetDouble("INITRADIUS")),
  initTemperature_(matdata->GetDouble("INITTEMPERATURE")),
  CPS_(matdata->GetDouble("CPS")),
  CPL_(matdata->GetDouble("CPL")),
  SL_latent_heat_max_(matdata->GetDouble("SL_LATENT_HEAT")),
  SL_transitionTemperature_(matdata->GetDouble("SL_TRANSITION_TEMPERATURE")),
  S_thermalExpansion_(matdata->GetDouble("S_THERMAL_EXPANSION")),
  L_thermalExpansion_(matdata->GetDouble("L_THERMAL_EXPANSION")),
  SL_thermalExpansion_(matdata->GetDouble("SL_THERMAL_EXPANSION"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::MeshFreeAMMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::MeshFreeAMMat(this));
}

MAT::MeshFreeAMMatType MAT::MeshFreeAMMatType::instance_;

DRT::ParObject* MAT::MeshFreeAMMatType::Create( const std::vector<char> & data )
{
  MAT::MeshFreeAMMat* MeshFreeAMMat = new MAT::MeshFreeAMMat();
  MeshFreeAMMat->Unpack(data);
  return MeshFreeAMMat;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MeshFreeAMMat::MeshFreeAMMat()
  : params_(NULL)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::MeshFreeAMMat::MeshFreeAMMat(MAT::PAR::MeshFreeAMMat* params)
  : params_(params)
{
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MeshFreeAMMat::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = -1;
  if (params_ != NULL) matid = params_->Id();  // in case we are in post-process mode
  AddtoPack(data,matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MAT::MeshFreeAMMat::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid and recover params_
  int matid;
  ExtractfromPack(position,data,matid);
  params_ = NULL;
  if (DRT::Problem::Instance()->Materials() != Teuchos::null)
    if (DRT::Problem::Instance()->Materials()->Num() != 0)
    {
      const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
      MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
      if (mat->Type() == MaterialType())
        params_ = static_cast<MAT::PAR::MeshFreeAMMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}
