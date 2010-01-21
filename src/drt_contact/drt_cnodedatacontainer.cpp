/*!----------------------------------------------------------------------
\file drt_cnodedatacontainer.cpp
\brief A class for additional data of a contact node

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
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_cnode.H"
#include "drt_cnodedatacontainer.H"
#include "../drt_lib/drt_dserror.H"
#include "drt_celement.H"
#include "contactdefines.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::CNodeDataContainer::CNodeDataContainer():
activeold_(false),
slip_(false),
grow_(1.0e12)
{
  for (int i=0;i<3;++i)
  {
    n()[i]=0.0;
    txi()[i]=0.0;
    teta()[i]=0.0;
  	lm()[i]=0.0;
    lmold()[i]=0.0;
    lmuzawa()[i]=0.0;
    jump()[i]=0.0;
    traction()[i]=0.0;
    tractionold()[i]=0.0;
  }

  Kappa() = 1.0;

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::CNodeDataContainer::CNodeDataContainer(const CONTACT::CNodeDataContainer& old):
activeold_(old.activeold_),
slip_(old.slip_),
drows_(old.drows_),
mrows_(old.mrows_),
drowsold_(old.drowsold_),
mrowsold_(old.mrowsold_),
mmodrows_(old.mmodrows_),
drowsPG_(old.drowsPG_),
mrowsPG_(old.mrowsPG_),
drowsoldPG_(old.drowsoldPG_),
mrowsoldPG_(old.mrowsoldPG_),
snodes_(old.snodes_),
mnodes_(old.mnodes_),
mnodesold_(old.mnodesold_),
grow_(old.grow_)

{
  for (int i=0;i<3;++i)
  {
    n()[i]=old.n_[i];
  	txi()[i]=old.txi_[i];
    teta()[i]=old.teta_[i];
    lm()[i]=old.lm_[i];
    lmold()[i]=old.lmold_[i];
    lmuzawa()[i]=old.lmuzawa_[i];
    jump()[i]=old.jump_[i];
    traction()[i]=old.traction_[i];
    tractionold()[i]=old.tractionold_[i];
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of CNodeDataContainer and
    return pointer to it (public)                              mgit 01/10|
 *----------------------------------------------------------------------*/
CONTACT::CNodeDataContainer* CONTACT::CNodeDataContainer::Clone() const
{
  CONTACT::CNodeDataContainer* newnodedc = new CONTACT::CNodeDataContainer(*this);
  return newnodedc;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CNodeDataContainer::Pack(vector<char>& data) const
{
  // add n_
	DRT::ParObject::AddtoPack(data,n_,3);
  // add txi_
	DRT::ParObject::AddtoPack(data,txi_,3);
  // add teta_
  DRT::ParObject::AddtoPack(data,teta_,3);
  // add lm_
	DRT::ParObject::AddtoPack(data,lm_,3);
  // add lmold_
  DRT::ParObject::AddtoPack(data,lmold_,3);
  // add lmuzawa_
	DRT::ParObject::AddtoPack(data,lmuzawa_,3);
  // add jump_
	DRT::ParObject::AddtoPack(data,jump_,3);
	// add activeold_
	DRT::ParObject::AddtoPack(data,activeold_);
  // add slip_
  DRT::ParObject::AddtoPack(data,slip_);
  // add traction_
  DRT::ParObject::AddtoPack(data,traction_,3);
  // add tractionold_
  DRT::ParObject::AddtoPack(data,tractionold_,3);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 01/10|
 *----------------------------------------------------------------------*/
void CONTACT::CNodeDataContainer::Unpack(int& position, const vector<char>& data)
{
  // n_
	DRT::ParObject:: ExtractfromPack(position,data,n_,3);
  // txi_
	DRT::ParObject::ExtractfromPack(position,data,txi_,3);
  // teta_
  DRT::ParObject::ExtractfromPack(position,data,teta_,3);
  // lm_
	DRT::ParObject::ExtractfromPack(position,data,lm_,3);
  // lmold_
  DRT::ParObject::ExtractfromPack(position,data,lmold_,3);
  // lmuzawa_
	DRT::ParObject::ExtractfromPack(position,data,lmuzawa_,3);
  // jump_
	DRT::ParObject::ExtractfromPack(position,data,jump_,3);
  // activeold_
	DRT::ParObject::ExtractfromPack(position,data,activeold_);
  // slip_
  DRT::ParObject::ExtractfromPack(position,data,slip_);
  // traction_
	DRT::ParObject::ExtractfromPack(position,data,traction_,3);
  // tractionold_
	DRT::ParObject::ExtractfromPack(position,data,tractionold_,3);

	return;
}



#endif  // #ifdef CCADISCRET
