/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast.cpp
\brief

<pre>
   Maintainer: Alexander Seitz
               seitz@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15271
</pre>
*/


/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast.H"

#include "../drt_lib/drt_linedefinition.H"
#include "../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_mat/plasticelasthyper.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Plast<so3_ele,distype>::So3_Plast(
  int id,
  int owner
  )
: so3_ele(id,owner),
  So3_Base(),
  intpoints_(distype),
  last_plastic_defgrd_inverse_(Teuchos::null),
  DalphaK_last_iter_(Teuchos::null),
  last_alpha_isotropic_(Teuchos::null),
  last_alpha_kinematic_(Teuchos::null),
  activity_state_(Teuchos::null),
  KbbInv_(Teuchos::null),
  Kbd_(Teuchos::null),
  fbeta_(Teuchos::null),
  KbbInvHill_(Teuchos::null),
  KbdHill_(Teuchos::null),
  fbetaHill_(Teuchos::null),
  mDLp_last_iter_(Teuchos::null),
  deltaAlphaI_(Teuchos::null)
{
  numgpt_ = intpoints_.NumPoints();
  kintype_ = geo_nonlinear;

  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Plast<so3_ele,distype>::So3_Plast(
  const DRT::ELEMENTS::So3_Plast<so3_ele,distype>& old
  )
: so3_ele(old),
  So3_Base(),
  intpoints_(distype)
{
  numgpt_ = intpoints_.NumPoints();
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to  seitz 07/13 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Clone() const
{
  DRT::ELEMENTS::So3_Plast< so3_ele, distype>* newelement
    = new DRT::ELEMENTS::So3_Plast< so3_ele, distype>(*this);

  return newelement;
}


template<class so3_ele, DRT::Element::DiscretizationType distype>
const int DRT::ELEMENTS::So3_Plast<so3_ele,distype>::VOIGT3X3SYM_[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
template<class so3_ele, DRT::Element::DiscretizationType distype>
const int DRT::ELEMENTS::So3_Plast<so3_ele,distype>::VOIGT3X3NONSYM_[3][3] = {{0,3,5},{6,1,4},{8,7,2}};


/*----------------------------------------------------------------------*
 | pack data (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Pack(
  DRT::PackBuffer& data
  ) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  so3_ele::AddtoPack(data,type);
  // data_
  so3_ele::AddtoPack(data,data_);
  // kintype_
  so3_ele::AddtoPack(data,kintype_);
  // detJ_
  so3_ele::AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  so3_ele::AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    so3_ele::AddtoPack(data,invJ_[i]);

  // add base class Element
  so3_ele::Pack(data);

  // plasticity stuff
  int histsize=0;
  if (last_plastic_defgrd_inverse_!=Teuchos::null)
      histsize=last_plastic_defgrd_inverse_->size();
  so3_ele::AddtoPack(data,histsize);
  bool Hill=(bool)(mDLp_last_iter_!=Teuchos::null);
  so3_ele::AddtoPack(data,(int)Hill);

  if (histsize!=0)
    for (int i=0; i<histsize; i++)
    {
      so3_ele::AddtoPack(data,last_plastic_defgrd_inverse_->at(i));
      so3_ele::AddtoPack(data,last_alpha_isotropic_->at(i));
      so3_ele::AddtoPack(data,last_alpha_kinematic_->at(i));
      so3_ele::AddtoPack(data,(int)activity_state_->at(i));
    }

  return;

}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                     seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Unpack(
  const std::vector<char>& data
  )
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  so3_ele::ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class element data_
  std::vector<char> tmp(0);
  so3_ele::ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);
  // kintype_
  kintype_ = static_cast<GenKinematicType>( so3_ele::ExtractInt(position,data) );
  // detJ_
  so3_ele::ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  so3_ele::ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<nsd_,nsd_>(true));
  for (int i=0; i<size; ++i)
    so3_ele::ExtractfromPack(position,data,invJ_[i]);

  // extract base class Element
  std::vector<char> basedata(0);
  so3_ele::ExtractfromPack(position,data,basedata);
  so3_ele::Unpack(basedata);

   int histsize=so3_ele::ExtractInt(position,data);

   bool Hill=(bool)so3_ele::ExtractInt(position,data);

   // initialize
   last_plastic_defgrd_inverse_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >);
   last_alpha_isotropic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<1,1> >);
   last_alpha_kinematic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >);
   activity_state_              = Teuchos::rcp( new std::vector<bool>);

   LINALG::Matrix<3,3> tmp33;
   LINALG::Matrix<1,1> tmp11;
    for (int i=0; i<histsize; i++)
    {
      so3_ele::ExtractfromPack(position,data,tmp33);
      last_plastic_defgrd_inverse_->push_back(tmp33);
      so3_ele::ExtractfromPack(position,data,tmp11);
      last_alpha_isotropic_->push_back(tmp11);
      so3_ele::ExtractfromPack(position,data,tmp33);
      last_alpha_kinematic_->push_back(tmp33);
      activity_state_->push_back((bool)so3_ele::ExtractInt(position,data));
    }

    // initialize yield function specific stuff
    if (Hill)
    {
      KbbInvHill_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<8,8> >);
      KbdHill_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<8,numdofperelement_> >);
      fbetaHill_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >);
      mDLp_last_iter_                  = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >);
      deltaAlphaI_                     = Teuchos::rcp( new std::vector<double>);
      LINALG::Matrix<8,8> tmp88(true);
      LINALG::Matrix<8,numdofperelement_> tmp8x(true);
      LINALG::Matrix<8,1> tmp81(true);
      KbbInvHill_                      ->resize(numgpt_,tmp88);
      KbdHill_                         ->resize(numgpt_,tmp8x);
      fbetaHill_                       ->resize(numgpt_,tmp81);
      mDLp_last_iter_                  ->resize(numgpt_,tmp81);
      deltaAlphaI_                     ->resize(numgpt_,0.);
      DalphaK_last_iter_ = Teuchos::null;
      KbbInv_            = Teuchos::null;
      Kbd_               = Teuchos::null;
      fbeta_             = Teuchos::null;
    }
    else // von Mises
    {
      KbbInv_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<5,5> >);
      Kbd_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<5,numdofperelement_> >);
      fbeta_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >);
      DalphaK_last_iter_           = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >);
      LINALG::Matrix<5,5> tmp55(true);
      LINALG::Matrix<5,numdofperelement_> tmp5x(true);
      LINALG::Matrix<5,1> tmp51(true);
      KbbInv_                      ->resize(numgpt_,tmp55);
      Kbd_                         ->resize(numgpt_,tmp5x);
      fbeta_                       ->resize(numgpt_,tmp51);
      DalphaK_last_iter_           ->resize(numgpt_,tmp51);
      mDLp_last_iter_        = Teuchos::null;
      deltaAlphaI_           = Teuchos::null;
      KbbInvHill_            = Teuchos::null;
      KbdHill_               = Teuchos::null;
      fbetaHill_             = Teuchos::null;
    }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;

}  // Unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                              seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Print(std::ostream& os) const
{
  os << "So3_Plast ";
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)             seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Plast<so3_ele,distype>::ReadElement(
  const std::string& eletype,
  const std::string& eledistype,
  DRT::INPUT::LineDefinition* linedef
  )
{
  so3_ele::ReadElement(eletype,eledistype,linedef);

  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    kintype_ = geo_linear;
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
    kintype_ = geo_nonlinear;
  else
    dserror("Reading of SO3_PLAST element failed! KINEM unknown");

  // allocate plastic history variables
  last_plastic_defgrd_inverse_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >);
  last_alpha_isotropic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<1,1> >);
  last_alpha_kinematic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >);
  activity_state_              = Teuchos::rcp( new std::vector<bool>);
  LINALG::Matrix<3,3> tmp33(true);
  LINALG::Matrix<3,3> id2(true);
  for (int i=0; i<3; i++) id2(i,i)=1.;
  LINALG::Matrix<1,1> tmp11(true);
  last_plastic_defgrd_inverse_ ->resize(numgpt_,id2);
  last_alpha_isotropic_        ->resize(numgpt_,tmp11);
  last_alpha_kinematic_        ->resize(numgpt_,tmp33);
  activity_state_              ->resize(numgpt_,false);
  // hill plasticity related stuff
  if (HaveHillPlasticity())
  {
    KbbInvHill_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<8,8> >);
    KbdHill_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<8,numdofperelement_> >);
    fbetaHill_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >);
    mDLp_last_iter_                  = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >);
    deltaAlphaI_                     = Teuchos::rcp( new std::vector<double>);
    LINALG::Matrix<8,8> tmp88(true);
    LINALG::Matrix<8,numdofperelement_> tmp8x(true);
    LINALG::Matrix<8,1> tmp81(true);
    KbbInvHill_                      ->resize(numgpt_,tmp88);
    KbdHill_                         ->resize(numgpt_,tmp8x);
    fbetaHill_                       ->resize(numgpt_,tmp81);
    mDLp_last_iter_                  ->resize(numgpt_,tmp81);
    deltaAlphaI_                     ->resize(numgpt_,0.);

    // don't need
    DalphaK_last_iter_ = Teuchos::null;
    KbbInv_            = Teuchos::null;
    Kbd_               = Teuchos::null;
    fbeta_             = Teuchos::null;
  }
  else //von-Mises
  {
    KbbInv_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<5,5> >);
    Kbd_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<5,numdofperelement_> >);
    fbeta_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >);
    DalphaK_last_iter_           = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >);
    LINALG::Matrix<5,5> tmp55(true);
    LINALG::Matrix<5,numdofperelement_> tmp5x(true);
    LINALG::Matrix<5,1> tmp51(true);
    KbbInv_                      ->resize(numgpt_,tmp55);
    Kbd_                         ->resize(numgpt_,tmp5x);
    fbeta_                       ->resize(numgpt_,tmp51);
    DalphaK_last_iter_           ->resize(numgpt_,tmp51);

    // don't need
    mDLp_last_iter_        = Teuchos::null;
    deltaAlphaI_           = Teuchos::null;
    KbbInvHill_            = Teuchos::null;
    KbdHill_               = Teuchos::null;
    fbetaHill_             = Teuchos::null;
  }

  return true;

}  // ReadElement()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<so3_ele,distype>::UniqueParObjectId() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    // cast the most specialised element
    // otherwise cast fails, because hex8fbar == hex8
    const DRT::ELEMENTS::So_hex8fbar* ele
      = dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);
    if(ele != NULL)
      return So_hex8fbarPlastType::Instance().UniqueParObjectId();
    else
      return So_hex8PlastType::Instance().UniqueParObjectId();
    break;
  }  // hex8
  case DRT::Element::tet4:
    return So_tet4PlastType::Instance().UniqueParObjectId();
  case DRT::Element::hex27:
    return So_hex27PlastType::Instance().UniqueParObjectId();
    break;
  default:
    dserror("unknown element type!");
    break;
  }
  // Intel compiler needs a return
  return -1;

} // UniqueParObjectId()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Plast<so3_ele,distype>::ElementType() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    // cast the most specialised element
    // caution: otherwise does not work, because hex8fbar == hex8
    const DRT::ELEMENTS::So_hex8fbar* ele
      = dynamic_cast<const DRT::ELEMENTS::So_hex8fbar*>(this);
    if(ele != NULL)
      return So_hex8fbarPlastType::Instance();
    else
      return So_hex8PlastType::Instance();
    break;
  }
  case DRT::Element::tet4:
    return So_tet4PlastType::Instance();
  case DRT::Element::hex27:
    return So_hex27PlastType::Instance();
    break;
  default:
    dserror("unknown element type!");
    break;
  }
  // Intel compiler needs a return
  return So_hex8PlastType::Instance();

};  // ElementType()


/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline DRT::Node** DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Nodes()
{
  return so3_ele::Nodes();
}


/*----------------------------------------------------------------------*
 | get the material from so3 (public)                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline Teuchos::RCP<MAT::Material> DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Material(
  ) const
{
  return so3_ele::Material();
}


/*----------------------------------------------------------------------*
 | get the node Ids from so3 (public)                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
inline int DRT::ELEMENTS::So3_Plast<so3_ele,distype>::Id() const
{
  return so3_ele::Id();
}


/*----------------------------------------------------------------------*
 | return names of visualization data (public)              seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::VisNames(std::map<std::string,int>& names)
{
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1; // scalar
  std::string plastic_zone = "plastic_zone";
  names[plastic_zone] = 1; // scalar
  so3_ele::VisNames(names);

  return;
}  // VisNames()

/*----------------------------------------------------------------------*
 | return visualization data (public)                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Plast<so3_ele,distype>::VisData(const std::string& name, std::vector<double>& data)
{
  if (name == "accumulatedstrain")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double tmp=0.;
    for (int gp=0; gp<numgpt_; gp++)
      tmp+= AccumulatedStrain(gp);
    data[0] = tmp/numgpt_;
  }
  if (name == "plastic_zone")
  {
    bool plastic_history=false;
    bool curr_active=false;
    if ((int)data.size()!=1) dserror("size mismatch");
    for (int gp=0; gp<numgpt_; gp++)
    {
      if (AccumulatedStrain(gp)!=0.) plastic_history=true;
      if (activity_state_->at(gp)==true) curr_active=true;
    }
    data[0] = plastic_history+curr_active;
  }

  return so3_ele::VisData(name, data);

}  // VisData()

/*----------------------------------------------------------------------*
 | return if there is hill plasticity                       seitz 09/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Plast<so3_ele,distype>::HaveHillPlasticity()
{
  // get plastic hyperelastic material
     MAT::PlasticElastHyper* plmat = NULL;
     if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
       plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
     else
       dserror("so3_ssn_plast elements only with PlasticElastHyper material");
     return (bool)plmat->HaveHillPlasticity();
}


/*----------------------------------------------------------------------*/
// include the file at the end of so3_ssn_plast.cpp
#include "so3_ssn_plast_fwd.hpp"
