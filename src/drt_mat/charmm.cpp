#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "sys/types.h"
#include "sys/stat.h"
#include <math.h>
#include "charmm.H"
#include "../drt_so3/so_hex8.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MAT::PAR::CHARMM::CHARMM( Teuchos::RCP<MAT::PAR::Material> matdata )
: Parameter(matdata),
  density_(matdata->GetDouble("DENS"))
{
}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/

MAT::CHARMM::CHARMM()
  : params_(NULL)
{
}


MAT::CHARMM::CHARMM(MAT::PAR::CHARMM* params)
  : params_(params)
{
}


/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::CHARMM::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // matid
  int matid = params_->Id();
  AddtoPack(data,matid);
}

/*----------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
void MAT::CHARMM::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matid
  int matid;
  ExtractfromPack(position,data,matid);
  // in post-process mode we do not have any instance of DRT::Problem
  if (DRT::Problem::NumInstances() > 0)
  {
    const int probinst = DRT::Problem::Instance()->Materials()->GetReadFromProblem();
    MAT::PAR::Parameter* mat = DRT::Problem::Instance(probinst)->Materials()->ParameterById(matid);
    if (mat->Type() == MaterialType())
      params_ = static_cast<MAT::PAR::CHARMM*>(mat);
    else
      dserror("Type of parameter material %d does not fit to calling type %d", mat->Type(), MaterialType());
  }
  else
  {
    params_ = NULL;
  }

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


/*----------------------------------------------------------------------*
 |  MD-FE energy based coupling schema material                 rm 03/08|
 *----------------------------------------------------------------------*/
void MAT::CHARMM::Evaluate( const LINALG::Matrix<6,1>& glstrain,
                            LINALG::Matrix<6,6>& cmat,
                            LINALG::Matrix<6,1>& stress,
                            const int ele_ID,
                            const int gp,
                            DRT::Container& data_,
                            const double time,
                            Epetra_SerialDenseMatrix* xrefe,
                            Epetra_SerialDenseMatrix* xcurr)
{
	cout << "Too be implemented....." << endl;
	exit(1);
	return;
}


/*----------------------------------------------------------------------*
 |  File based API to CHARMM                                    rm 03/08|
 *----------------------------------------------------------------------*/
Epetra_SerialDenseVector MAT::CHARMM::charmmfileapi ( const double STARTD,
                                    const double ENDD,
                                    const LINALG::SerialDenseVector direction,
                                    LINALG::SerialDenseVector& charmm_result)
{
	cout << "Too be implemented....." << endl;
	exit(1);
	return(0);
}


#endif

