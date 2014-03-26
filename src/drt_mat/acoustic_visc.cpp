/*!----------------------------------------------------------------------
\file acoustic_visc.cpp
\brief contains a density, a speed of sound and a viscsity to characterize
       all necessary properties for sound propagation in acoustic material
       example input line:
       MAT 1 MAT_AcousticVisc DENSITY 1000.0 C 1500.0 VISC  1.0
<pre>
Maintainer: Svenja Schoeder
</pre>
*/
/*----------------------------------------------------------------------*/

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "acoustic_visc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::AcousticViscMat::AcousticViscMat(
  Teuchos::RCP<MAT::PAR::Material> matdata
  )
: Parameter(matdata),
  c_(matdata->GetDouble("C")),
  density_(matdata->GetDouble("DENSITY")),
  visc_(matdata->GetDouble("VISC")),
  therm_(matdata->GetDouble("THERM"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::AcousticViscMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AcousticViscMat(this));
}

MAT::AcousticViscMatType MAT::AcousticViscMatType::instance_;

DRT::ParObject* MAT::AcousticViscMatType::Create( const std::vector<char> & data )
{
  MAT::AcousticViscMat* soundprop = new MAT::AcousticViscMat();
  soundprop->Unpack(data);
  return soundprop;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::AcousticViscMat::AcousticViscMat()
  : params_(NULL)
{
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::AcousticViscMat::AcousticViscMat(MAT::PAR::AcousticViscMat* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::AcousticViscMat::Pack(DRT::PackBuffer& data) const
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


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::AcousticViscMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::AcousticViscMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}

