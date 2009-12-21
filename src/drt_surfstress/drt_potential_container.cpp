/*!----------------------------------------------------------------------
\file potential_container.cpp

\brief A container class which stores all necessary data for a potential element

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "../drt_surfstress/drt_potential_container.H"
#include "../drt_lib/drt_element.H"


/*----------------------------------------------------------------------*
 |  standard constructor                                       (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
POTENTIAL::PotentialElementContainer::PotentialElementContainer():
  ParObject(),
  gid_(-1),
  distype_(DRT::Element::dis_none),
  body_label_(-1)
{
  xyz_e_.Shape(1,1);
  XYZ_e_.Shape(1,1);
  lm_.resize(1,-1);
  return;
}



/*----------------------------------------------------------------------*
 |  standard constructor                                       (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
POTENTIAL::PotentialElementContainer::PotentialElementContainer(
    const int                                 gid,
    const DRT::Element::DiscretizationType    distype,
    const int                                 body_label,
    const Epetra_SerialDenseMatrix&           xyz_e,
    const Epetra_SerialDenseMatrix&           XYZ_e,
    const std::vector<int>&      	            lm):
      ParObject(),
  gid_(gid),
  distype_(distype),
  body_label_(body_label),
  xyz_e_(xyz_e),
  XYZ_e_(XYZ_e),
  lm_(lm)
{
  return;
}

    
    
/*----------------------------------------------------------------------*
 |  dtor (public)                                           u.may 12/09 |
 *----------------------------------------------------------------------*/
POTENTIAL::PotentialElementContainer::~PotentialElementContainer()
{
  return;
}


/*----------------------------------------------------------------------*
 |  copy constructor                                           (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
POTENTIAL::PotentialElementContainer::PotentialElementContainer(
	const PotentialElementContainer&  old):
	  ParObject(old),
  gid_(old.gid_),
  distype_(old.distype_),
  body_label_(old.body_label_),
  xyz_e_(old.xyz_e_),
  XYZ_e_(old.XYZ_e_),
  lm_(old.lm_)
{
 return;
}



/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
void POTENTIAL::PotentialElementContainer::Pack(vector<char>& data) const
{
  // TODO check if data.resize(0);

  // global id gid_
  AddtoPack(data, gid_);
  // distype  
  AddtoPack(data, distype_);
  // body_label_
  AddtoPack(data, body_label_);
  // xyz_e_
  AddtoPack(data, xyz_e_);
  // XYZ_e_
  AddtoPack(data, XYZ_e_);
  // lm_
  AddtoPack(data, lm_);

  return;
}



/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
void POTENTIAL::PotentialElementContainer::Unpack(
  const vector<char>&     data)
{
  int position = 0;
  // gid
  ExtractfromPack(position,data, gid_);
  // distype_
  ExtractfromPack(position,data, distype_);
  // body_label_
  ExtractfromPack(position,data, body_label_);
  // xyz_e_
  ExtractfromPack(position,data, xyz_e_);
  // XYZ_e_
  ExtractfromPack(position,data, XYZ_e_);
  // lm_
  ExtractfromPack(position,data, lm_);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}



/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
void POTENTIAL::PotentialElementContainer::Unpack(
	const vector<char>& 		data, 
	int& 			              position)
{
  // gid
  ExtractfromPack(position,data, gid_);
  // distype_
  ExtractfromPack(position,data, distype_);
  // body_label_
  ExtractfromPack(position,data, body_label_);
  // xyz_e_
  ExtractfromPack(position,data, xyz_e_);
  // XYZ_e_
  ExtractfromPack(position,data, XYZ_e_);
  // lm_
  ExtractfromPack(position,data, lm_);

  // if (position != (int)data.size())
  //  dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


#endif  // #ifdef CCADISCRET



