/*!----------------------------------------------------------------------
\file acoustic_sol.cpp
\brief contains a density, a speed of sound and a solsity to characterize
       all necessary properties for sound propagation in acoustic material
       example input line:
       MAT 1 MAT_AcousticSol DENSITY 1000.0 C 1500.0 VISC 1.0
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
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "acoustic_sol.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_bundle.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::AcousticSolMat::AcousticSolMat(Teuchos::RCP<MAT::PAR::Material> matdata)
    : Parameter(matdata),
      c_(matdata->GetDouble("C")),
      density_(matdata->GetDouble("DENSITY")),
      visc_(matdata->GetDouble("VISC"))
{
}


Teuchos::RCP<MAT::Material> MAT::PAR::AcousticSolMat::CreateMaterial()
{
  return Teuchos::rcp(new MAT::AcousticSolMat(this));
}

MAT::AcousticSolMatType MAT::AcousticSolMatType::instance_;

DRT::ParObject* MAT::AcousticSolMatType::Create(const std::vector<char>& data)
{
  MAT::AcousticSolMat* soundprop = new MAT::AcousticSolMat();
  soundprop->Unpack(data);
  return soundprop;
}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::AcousticSolMat::AcousticSolMat() : params_(NULL) {}

/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
MAT::AcousticSolMat::AcousticSolMat(MAT::PAR::AcousticSolMat* params) : params_(params) {}


/*----------------------------------------------------------------------*
 |                                                                      |
 *----------------------------------------------------------------------*/
void MAT::AcousticSolMat::Pack(DRT::PackBuffer& data) const
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
void MAT::AcousticSolMat::Unpack(const std::vector<char>& data)
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
        params_ = static_cast<MAT::PAR::AcousticSolMat*>(mat);
      else
        dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(),
            MaterialType());
    }

  if (position != data.size()) dserror("Mismatch in size of data %d <-> %d", data.size(), position);
}
