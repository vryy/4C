/*!----------------------------------------------------------------------
\file inversedesign.cpp

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)

#include "inversedesign.H"
#include "../drt_lib/drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::InvDesign::InvDesign(const int numnod, const int ngp,  const bool istet4) :
ParObject(),
numnod_(numnod),
ngp_(ngp),
isinit_(false)
{
  // allocate history memory
  detJ_.resize(ngp);

  Fhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp,9));
  LINALG::FixedSizeSerialDenseMatrix<3,3> F(true); // set to zero
  F(0,0) = F(1,1) = F(2,2) = 1.0;
  for (int i=0; i<ngp; ++i) MatrixtoStorage(i,F,FHistory());

  if (!istet4)
    invJhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp,9));
  else
    invJhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp,12));
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::InvDesign::InvDesign(const DRT::ELEMENTS::InvDesign& old) :
ParObject(old),
numnod_(old.numnod_),
ngp_(old.ngp_),
isinit_(old.isinit_),
Fhist_(Teuchos::rcp(new Epetra_SerialDenseMatrix(old.FHistory()))),
invJhist_(Teuchos::rcp(new Epetra_SerialDenseMatrix(old.JHistory()))),
detJ_(old.detJ_)
{
  return;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // numnod_
  AddtoPack(data,numnod_);

  // ngp_
  AddtoPack(data,ngp_);

  // isinit_
  AddtoPack(data,isinit_);

  // Fhist_
  AddtoPack(data,*Fhist_);

  // invJhist_
  AddtoPack(data,*invJhist_);

  // detJ_
  AddtoPack(data,detJ_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // numnod_
  ExtractfromPack(position,data,numnod_);

  // ngp_
  ExtractfromPack(position,data,ngp_);

  // isinit_
  ExtractfromPack(position,data,isinit_);

  // Fhist_
  ExtractfromPack(position,data,*Fhist_);

  // invJhist_
  ExtractfromPack(position,data,*invJhist_);

  // detJ_
  ExtractfromPack(position,data,detJ_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}







#endif  // #if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#endif  // #ifdef CCADISCRET
