/*!----------------------------------------------------------------------
\file meshtying_element.cpp
\brief

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "meshtying_element.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtElement::MtElement(int id, ElementType etype, int owner,
                              const DRT::Element::DiscretizationType& shape,
                              const int numnode,
                              const int* nodeids,
                              const bool isslave) :
MORTAR::MortarElement(id,etype,owner,shape,numnode,nodeids,isslave)
{
  // empty constructor
  
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtElement::MtElement(const CONTACT::MtElement& old) :
MORTAR::MortarElement(old)
{
  // empty copy-constructor
  
  return;
}

/*----------------------------------------------------------------------*
 |  clone-ctor (public)                                      mwgee 10/07|
 *----------------------------------------------------------------------*/
CONTACT::MtElement* CONTACT::MtElement::Clone() const
{
  CONTACT::MtElement* newele = new CONTACT::MtElement(*this);
  return newele;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 10/07|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const CONTACT::MtElement& element)
{
  element.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                   mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtElement::Print(ostream& os) const
{
  os << "Meshtying ";
  MORTAR::MortarElement::Print(os);

  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtElement::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  
  // add base class MORTAR::MortarElement
  vector<char> basedata(0);
  MORTAR::MortarElement::Pack(basedata);
  AddtoPack(data,basedata);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                           mwgee 10/07|
 *----------------------------------------------------------------------*/
void CONTACT::MtElement::Unpack(const vector<char>& data)
{
  int position = 0;
  
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  
  // extract base class MORTAR::MortarElement
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  MORTAR::MortarElement::Unpack(basedata);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}



/*----------------------------------------------------------------------*
 |  evaluate element (public)                                mwgee 10/07|
 *----------------------------------------------------------------------*/
int CONTACT::MtElement::Evaluate(ParameterList&            params,
                                 DRT::Discretization&      discretization,
                                 vector<int>&              lm,
                                 Epetra_SerialDenseMatrix& elemat1,
                                 Epetra_SerialDenseMatrix& elemat2,
                                 Epetra_SerialDenseVector& elevec1,
                                 Epetra_SerialDenseVector& elevec2,
                                 Epetra_SerialDenseVector& elevec3)
{
  dserror("CONTACT::MtElement::Evaluate not implemented!");
  return -1;
}

#endif  // #ifdef CCADISCRET
