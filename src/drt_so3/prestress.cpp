/*!----------------------------------------------------------------------
\file prestress.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET
#if defined(PRESTRESS) || defined(POSTSTRESS)

#include "prestress.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 07/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PreStress::PreStress(const int numnode,
                                    const int ngp,
                                    const bool istet4) :
ParObject(),
isinit_(false),
numnode_(numnode)
{
  // allocate history memory
  Fhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp,9));
  if (!istet4)
    invJhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp,9));
  else
    invJhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp,12));

  // init the deformation gradient history
  LINALG::Matrix<3,3> F(true);
  F(0,0) = F(1,1) = F(2,2) = 1.0;
  for (int i=0; i<NGP(); ++i) MatrixtoStorage(i,F,FHistory());
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::PreStress::PreStress(const DRT::ELEMENTS::PreStress& old) :
ParObject(old),
isinit_(old.isinit_),
numnode_(old.numnode_),
Fhist_(Teuchos::rcp(new Epetra_SerialDenseMatrix(old.FHistory()))),
invJhist_(Teuchos::rcp(new Epetra_SerialDenseMatrix(old.JHistory())))
{
  return;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PreStress::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // pack isinit_
  AddtoPack(data,isinit_);

  // pack numnode_
  AddtoPack(data,numnode_);

  // pack Fhist_
  AddtoPack(data,*Fhist_);

  // pack invJhist_
  AddtoPack(data,*invJhist_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::PreStress::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract isinit_
  ExtractfromPack(position,data,isinit_);

  // extract numnode_
  ExtractfromPack(position,data,numnode_);

  // extract Fhist_
  ExtractfromPack(position,data,*Fhist_);

  // extract invJhist_
  ExtractfromPack(position,data,*invJhist_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}







#endif  // #if defined(PRESTRESS) || defined(POSTSTRESS)
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
