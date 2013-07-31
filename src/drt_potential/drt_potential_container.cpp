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

#include "../drt_potential/drt_potential_container.H"
#include "../drt_lib/drt_element.H"

POTENTIAL::PecType POTENTIAL::PecType::instance_;


DRT::ParObject* POTENTIAL::PecType::Create( const std::vector<char> & data )
{
  return NULL;
}


/*----------------------------------------------------------------------*
 |  standard constructor                                       (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
POTENTIAL::PotentialElementContainer::PotentialElementContainer():
  ParObject(),
  gid_(-1),
  distype_(DRT::Element::dis_none),
  body_label_(-1),
  beta_(-1.0)
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
    const double							  beta,
    const Epetra_SerialDenseMatrix&           xyz_e,
    const Epetra_SerialDenseMatrix&           XYZ_e,
    const std::vector<int>&      	          lm):
      ParObject(),
  gid_(gid),
  distype_(distype),
  body_label_(body_label),
  beta_(beta),
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
  beta_(old.beta_),
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
void POTENTIAL::PotentialElementContainer::Pack(DRT::PackBuffer& data) const
{
  // global id gid_
  AddtoPack(data, gid_);
  // distype
  AddtoPack(data, distype_);
  // body_label_
  AddtoPack(data, body_label_);
  // beta
  AddtoPack(data, beta_);
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
  const std::vector<char>&     data)
{
  dserror("this unpack method is not used");
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
void POTENTIAL::PotentialElementContainer::Unpack(
	const std::vector<char>& 			data,
	std::vector<char>::size_type&  position)
{
  // gid
  ExtractfromPack(position,data, gid_);
  // distype_
  distype_ = static_cast<DRT::Element::DiscretizationType>( ExtractInt(position,data) );
  // body_label_
  ExtractfromPack(position,data, body_label_);
  // beta_
  ExtractfromPack(position,data, beta_);
  // xyz_e_
  ExtractfromPack(position,data, xyz_e_);
  // XYZ_e_
  ExtractfromPack(position,data, XYZ_e_);
  // lm_
  ExtractfromPack(position,data, lm_);

  return;
}


/*----------------------------------------------------------------------*
 |  Print                                                      (public) |
 |                                                          u.may 12/09 |
 *----------------------------------------------------------------------*/
void POTENTIAL::PotentialElementContainer::Print()
{
  std::cout << "Print: Potential Element Container "           << std::endl;
  std::cout << "PEC: global Id = "             << gid_         << std::endl;
  std::cout << "PEC: discretization type = "   << distype_     << std::endl;
  std::cout << "PEC: label = "                 << body_label_  << std::endl;
  std::cout << "PEC: beta = "                  << beta_  	  << std::endl;
  std::cout << "PEC: number of nodes = "       << xyz_e_.M()   << std::endl;
  std::cout << std::endl;
  std::cout << "PEC: spatial configuration";
  xyz_e_.Print(std::cout); std::cout << std::endl;
  std::cout << "PEC: reference configuration";
  XYZ_e_.Print(std::cout); std::cout << std::endl;

  std::cout << "PEC: global dof ids " << std::endl;
  for(int i_lm = 0; i_lm < (int) lm_.size(); i_lm++)
    std::cout << "lm[" << i_lm << "] =      " << lm_[i_lm] << std::endl;
  std::cout << std::endl;
}




