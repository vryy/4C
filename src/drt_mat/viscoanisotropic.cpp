/*!----------------------------------------------------------------------
\file viscoanisotropic.cpp
\brief

<pre>
Maintainer: Moritz Frenzel & Thomas Kloeppel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>
*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <vector>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include "Epetra_SerialDenseSolver.h"
#include "viscoanisotropic.H"
#include "../drt_lib/linalg_serialdensevector.H"


extern struct _MATERIAL *mat;  ///< C-style material struct


/*----------------------------------------------------------------------*
 |  Constructor                                   (public)         05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoAnisotropic::ViscoAnisotropic()
  : matdata_(NULL)
{
}


/*----------------------------------------------------------------------*
 |  Copy-Constructor                             (public)          05/08|
 *----------------------------------------------------------------------*/
MAT::ViscoAnisotropic::ViscoAnisotropic(MATERIAL* matdata)
  : matdata_(matdata)
{
}


/*----------------------------------------------------------------------*
 |  Pack                                          (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // matdata
  int matdata = matdata_ - mat;   // pointer difference to reach 0-entry
  AddtoPack(data,matdata);
  //Pack history data
  int histsize;
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack                                        (public)         05/08|
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Unpack(const vector<char>& data)
{
  isinit_=true;
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // matdata
  int matdata;
  ExtractfromPack(position,data,matdata);
  matdata_ = &mat[matdata];     // unpack pointer to my specific matdata_
  // history data
  int twicehistsize;
  ExtractfromPack(position,data,twicehistsize);


  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Update()
{
  histstresslast_=histstresscurr_;
  artstresslast_=artstresscurr_;
  const LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> emptyvec;//6 stresses for 3D
  histstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  artstresscurr_=rcp(new vector<LINALG::FixedSizeSerialDenseMatrix<NUM_STRESS_3D,1> >);
  const int histsize=histstresslast_->size();
  histstresscurr_->resize(histsize);
  artstresscurr_->resize(histsize);
  for (int j=0; j<histsize; ++j)
  {
    histstresscurr_->at(j) = emptyvec;
    artstresscurr_->at(j) = emptyvec;
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void MAT::ViscoAnisotropic::Evaluate
(
  const Epetra_SerialDenseVector* glstrain,  ///<green lagrange strain
  const int gp,   ///< number of Gauss points
  Teuchos::ParameterList& params,  ///< parameter list for communication
  Epetra_SerialDenseMatrix* cmat,  ///< material stiffness matrix
  Epetra_SerialDenseVector* stress  ///< 2nd PK-stress
)
{
  return;
}


#endif

