/*----------------------------------------------------------------------*/
/*! \file

\brief special element adaptions for inverse design
\maintainer Christoph Meier
\level 2

*----------------------------------------------------------------------*/

#include "inversedesign.H"
#include "../drt_lib/drt_dserror.H"


DRT::ELEMENTS::InvDesignType DRT::ELEMENTS::InvDesignType::instance_;


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::InvDesign::InvDesign(const int numnod, const int ngp, const bool istet4)
    : ParObject(), numnod_(numnod), ngp_(ngp), isinit_(false)
{
  // allocate history memory
  detJ_.resize(ngp);

  Fhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp, 9));
  LINALG::Matrix<3, 3> F(true);  // set to zero
  F(0, 0) = F(1, 1) = F(2, 2) = 1.0;
  for (int i = 0; i < ngp; ++i) MatrixtoStorage(i, F, FHistory());

  if (!istet4)
    invJhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp, 9));
  else
    invJhist_ = Teuchos::rcp(new Epetra_SerialDenseMatrix(ngp, 12));
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::InvDesign::InvDesign(const DRT::ELEMENTS::InvDesign& old)
    : ParObject(old),
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
void DRT::ELEMENTS::InvDesign::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm(data);
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data, type);

  // numnod_
  AddtoPack(data, numnod_);

  // ngp_
  AddtoPack(data, ngp_);

  // isinit_
  AddtoPack(data, isinit_);

  // Fhist_
  AddtoPack(data, *Fhist_);

  // invJhist_
  AddtoPack(data, *invJhist_);

  // detJ_
  AddtoPack(data, detJ_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 08/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::InvDesign::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position, data, type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // numnod_
  ExtractfromPack(position, data, numnod_);

  // ngp_
  ExtractfromPack(position, data, ngp_);

  // isinit_
  isinit_ = ExtractInt(position, data);

  // Fhist_
  ExtractfromPack(position, data, *Fhist_);

  // invJhist_
  ExtractfromPack(position, data, *invJhist_);

  // detJ_
  ExtractfromPack(position, data, detJ_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d", (int)data.size(), position);
  return;
}


int DRT::ELEMENTS::InvDesign::UniqueParObjectId() const
{
  return InvDesignType::Instance().UniqueParObjectId();
}
