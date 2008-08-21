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
DRT::ELEMENTS::InvDesign::InvDesign(const int numnod, const int ngp) :
ParObject(),
numnod_(numnod),
ngp_(ngp)
{
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 08/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::InvDesign::InvDesign(const DRT::ELEMENTS::InvDesign& old) :
ParObject(old),
numnod_(old.numnod_),
ngp_(old.ngp_)
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


  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}







#endif  // #if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#endif  // #ifdef CCADISCRET
