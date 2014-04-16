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
  stab_s_(-1.),
  cpl_(-1.),
  last_plastic_defgrd_inverse_(Teuchos::null),
  last_alpha_isotropic_(Teuchos::null),
  last_alpha_kinematic_(Teuchos::null),
  activity_state_(Teuchos::null),
  KbbInv_(Teuchos::null),
  Kbd_(Teuchos::null),
  fbeta_(Teuchos::null),
  DalphaK_last_iter_(Teuchos::null),
  KbbInvHill_(Teuchos::null),
  KbdHill_(Teuchos::null),
  fbetaHill_(Teuchos::null),
  mDLp_last_iter_(Teuchos::null),
  deltaAlphaI_(Teuchos::null),
  KaaInv_(Teuchos::null),
  Kad_(Teuchos::null),
  feas_(Teuchos::null),
  Kba_(Teuchos::null),
  alpha_eas_(Teuchos::null),
  alpha_eas_last_timestep_(Teuchos::null),
  alpha_eas_delta_over_last_timestep_(Teuchos::null),
  eastype_(soh8p_easnone),
  neas_(0)
{
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
  So3_Base()
{
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

  // Gauss points and weights
  const int size2 = (int)xsi_.size();
  so3_ele::AddtoPack(data,size2);
  for (int i=0; i<size2;++i)
    so3_ele::AddtoPack(data,xsi_[i]);
  so3_ele::AddtoPack(data,wgt_);

  // parameters
  so3_ele::AddtoPack(data,stab_s_);
  so3_ele::AddtoPack(data,cpl_);

  // plasticity stuff
  int histsize=0;
  if (last_plastic_defgrd_inverse_!=Teuchos::null)
      histsize=last_plastic_defgrd_inverse_->size();
  so3_ele::AddtoPack(data,histsize);
  bool Hill=(bool)(mDLp_last_iter_!=Teuchos::null);
  so3_ele::AddtoPack(data,(int)Hill);

  // EAS element technology
  so3_ele::AddtoPack(data,(int)eastype_);
  so3_ele::AddtoPack(data,neas_);
  if (eastype_!=soh8p_easnone)
  {
    so3_ele::AddtoPack(data,(*alpha_eas_));
    so3_ele::AddtoPack(data,(*alpha_eas_last_timestep_));
    so3_ele::AddtoPack(data,(*alpha_eas_delta_over_last_timestep_));
  }

  // history at each Gauss point
  if (histsize!=0)
    for (int i=0; i<histsize; i++)
    {
      so3_ele::AddtoPack(data,last_plastic_defgrd_inverse_->at(i));
      so3_ele::AddtoPack(data,last_alpha_isotropic_->at(i));
      so3_ele::AddtoPack(data,last_alpha_kinematic_->at(i));
      so3_ele::AddtoPack(data,(int)activity_state_->at(i));
      if (Hill)
        so3_ele::AddtoPack(data,mDLp_last_iter_->at(i));
      else
        so3_ele::AddtoPack(data,DalphaK_last_iter_->at(i));
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

  // Gauss points and weights
  int size2 = so3_ele::ExtractInt(position,data);
  xsi_.resize(size2, LINALG::Matrix<nsd_,1>(true));
  for (int i=0; i<size2;++i)
    so3_ele::ExtractfromPack(position,data,xsi_[i]);
  so3_ele::ExtractfromPack(position,data,wgt_);
  numgpt_=wgt_.size();

  // paramters
  so3_ele::ExtractfromPack(position,data,stab_s_);
  so3_ele::ExtractfromPack(position,data,cpl_);

   int histsize=so3_ele::ExtractInt(position,data);

   bool Hill=(bool)so3_ele::ExtractInt(position,data);

   // EAS element technology
   eastype_=static_cast<EASType>(so3_ele::ExtractInt(position,data));
   so3_ele::ExtractfromPack(position,data,neas_);
   if ((int)eastype_!=neas_)
     dserror("mismatch in EAS");

   // no EAS
   if (eastype_==soh8p_easnone)
   {
     KaaInv_                             =Teuchos::null;
     Kad_                                =Teuchos::null;
     feas_                               =Teuchos::null;
     Kba_                                =Teuchos::null;
     alpha_eas_                          =Teuchos::null;
     alpha_eas_last_timestep_            =Teuchos::null;
     alpha_eas_delta_over_last_timestep_ =Teuchos::null;
   }
   else
   {
     KaaInv_                             =Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_,neas_));
     Kad_                                =Teuchos::rcp(new Epetra_SerialDenseMatrix(neas_,numdofperelement_));
     feas_                               =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
     Kba_                                =Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(histsize,Epetra_SerialDenseMatrix(5,neas_)));
     alpha_eas_                          =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
     alpha_eas_last_timestep_            =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
     alpha_eas_delta_over_last_timestep_ =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
   }

   // initialize plastic history data
   last_plastic_defgrd_inverse_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >(numgpt_));
   last_alpha_isotropic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<1,1> >(numgpt_));
   last_alpha_kinematic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >(numgpt_));
   activity_state_              = Teuchos::rcp( new std::vector<bool>(numgpt_));

   // initialize yield function specific stuff
   if (Hill)
   {
     LINALG::Matrix<8,8> tmp88(true);
     LINALG::Matrix<8,numdofperelement_> tmp8x(true);
     LINALG::Matrix<8,1> tmp81(true);
     KbbInvHill_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<8,8> >(numgpt_,tmp88));
     KbdHill_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<8,numdofperelement_> >(numgpt_,tmp8x));
     fbetaHill_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >(numgpt_,tmp81));
     mDLp_last_iter_                  = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >(numgpt_));
     deltaAlphaI_                     = Teuchos::rcp( new std::vector<double>(numgpt_,0.));

     // don't need
     DalphaK_last_iter_     = Teuchos::null;
     KbbInv_                = Teuchos::null;
     Kbd_                   = Teuchos::null;
     fbeta_                 = Teuchos::null;
   }
   else // von Mises
   {
     LINALG::Matrix<5,5> tmp55(true);
     LINALG::Matrix<5,numdofperelement_> tmp5x(true);
     LINALG::Matrix<5,1> tmp51(true);
     KbbInv_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<5,5> >(numgpt_,tmp55));
     Kbd_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<5,numdofperelement_> >(numgpt_,tmp5x));
     fbeta_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >(numgpt_,tmp51));
     DalphaK_last_iter_           = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >(numgpt_));

     // don't need
     mDLp_last_iter_        = Teuchos::null;
     deltaAlphaI_           = Teuchos::null;
     KbbInvHill_            = Teuchos::null;
     KbdHill_               = Teuchos::null;
     fbetaHill_             = Teuchos::null;
   }

   if (eastype_!=soh8p_easnone)
   {
     so3_ele::ExtractfromPack(position,data,(*alpha_eas_));
     so3_ele::ExtractfromPack(position,data,(*alpha_eas_last_timestep_));
     so3_ele::ExtractfromPack(position,data,(*alpha_eas_delta_over_last_timestep_));
   }

   for (int i=0; i<histsize; i++)
    {
      so3_ele::ExtractfromPack(position,data,(*last_plastic_defgrd_inverse_)[i]);
      so3_ele::ExtractfromPack(position,data,(*last_alpha_isotropic_)[i]);
      so3_ele::ExtractfromPack(position,data,(*last_alpha_kinematic_)[i]);
      (*activity_state_)[i]=(bool)(so3_ele::ExtractInt(position,data));
      if (Hill)
        so3_ele::ExtractfromPack(position,data,(*mDLp_last_iter_)[i]);
      else
        so3_ele::ExtractfromPack(position,data,(*DalphaK_last_iter_)[i]);
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

  // quadrature
  if (linedef->HaveNamed("NUMGP"))
  {
    if (ElementType()!=DRT::ELEMENTS::So_hex8PlastType::Instance())
      dserror("You may only choose the Gauss point number for SOLIDH8PLAST");

    int ngp =0;
    linedef->ExtractInt("NUMGP",ngp);

    switch(ngp)
    {
    case 8:
    {
      DRT::UTILS::GaussIntegration ip(distype,3);
      numgpt_=ip.NumPoints();
      xsi_.resize(numgpt_);
      wgt_.resize(numgpt_);
      for (int gp=0; gp<numgpt_; ++gp)
      {
       wgt_[gp]=ip.Weight(gp);
        const double* gpcoord = ip.Point(gp);
        for (int idim=0; idim<nsd_; idim++)
          xsi_[gp](idim) = gpcoord[idim];
      }
      break;
    }
    case 9:
    {
      DRT::UTILS::GaussIntegration ip(distype,3);
      numgpt_=ip.NumPoints()+1;
      xsi_.resize(numgpt_);
      wgt_.resize(numgpt_);
      for (int gp=0; gp<numgpt_-1; ++gp)
      {
        wgt_[gp]=5./9.;
        const double* gpcoord = ip.Point(gp);
        for (int idim=0; idim<nsd_; idim++)
          xsi_[gp](idim) = gpcoord[idim];
      }
      // 9th quadrature point at element center
      xsi_[numgpt_-1](0) = 0.;
      xsi_[numgpt_-1](1) = 0.;
      xsi_[numgpt_-1](2) = 0.;
      wgt_[numgpt_-1] = 32./9.;
      break;
    }
    case 27:
    {
      DRT::UTILS::GaussIntegration ip(distype,4);
      numgpt_=ip.NumPoints();
      xsi_.resize(numgpt_);
      wgt_.resize(numgpt_);
      for (int gp=0; gp<numgpt_; ++gp)
      {
       wgt_[gp]=ip.Weight(gp);
        const double* gpcoord = ip.Point(gp);
        for (int idim=0; idim<nsd_; idim++)
          xsi_[gp](idim) = gpcoord[idim];
      }
      break;
    }
    default:
      dserror("so3_plast doesn't know what to do with %i Gauss points",ngp);
      break;
    }
  }
  else // default integration
  {
    DRT::UTILS::GaussIntegration ip(distype,3);
    numgpt_=ip.NumPoints();
    xsi_.resize(numgpt_);
    wgt_.resize(numgpt_);
    for (int gp=0; gp<numgpt_; ++gp)
    {
     wgt_[gp]=ip.Weight(gp);
      const double* gpcoord = ip.Point(gp);
      for (int idim=0; idim<nsd_; idim++)
        xsi_[gp](idim) = gpcoord[idim];
    }
  }

  // EAS
  if (linedef->HaveNamed("EAS"))
  {
    if (distype != DRT::Element::hex8)
      dserror("EAS in so3 plast currently only for HEX8 elements");

    if (ElementType()==DRT::ELEMENTS::So_hex8fbarPlastType::Instance())
      dserror("no combination of Fbar and EAS ... And, by the way, how did you get here???");
    linedef->ExtractString("EAS",buffer);

    if (buffer == "none")
      eastype_=soh8p_easnone;
    else if (buffer == "mild")
      eastype_=soh8p_easmild;
    else if (buffer=="full")
      eastype_=soh8p_easfull;
    else
      dserror("unknown EAS type for so3_plast");

    if (HaveHillPlasticity() && eastype_!=soh8p_easnone)
      dserror("no EAS for Hill-plasticity yet");
  }
  else
    eastype_=soh8p_easnone;

  // initialize EAS data
  EasInit();

  LINALG::Matrix<3,3> tmp33(true);
  LINALG::Matrix<3,3> id2(true);
  for (int i=0; i<3; i++) id2(i,i)=1.;
  LINALG::Matrix<1,1> tmp11(true);
  // allocate plastic history variables
  last_plastic_defgrd_inverse_ = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >(numgpt_,id2));
  last_alpha_isotropic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<1,1> >(numgpt_,tmp11));
  last_alpha_kinematic_        = Teuchos::rcp( new std::vector<LINALG::Matrix<3,3> >(numgpt_,tmp33));
  activity_state_              = Teuchos::rcp( new std::vector<bool>(numgpt_,false));

  // hill plasticity related stuff
  if (HaveHillPlasticity())
  {
    LINALG::Matrix<8,8> tmp88(true);
    LINALG::Matrix<8,numdofperelement_> tmp8x(true);
    LINALG::Matrix<8,1> tmp81(true);
    KbbInvHill_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<8,8> >(numgpt_,tmp88));
    KbdHill_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<8,numdofperelement_> >(numgpt_,tmp8x));
    fbetaHill_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >(numgpt_,tmp81));
    mDLp_last_iter_                  = Teuchos::rcp( new std::vector<LINALG::Matrix<8,1> >(numgpt_,tmp81));
    deltaAlphaI_                     = Teuchos::rcp( new std::vector<double>(numgpt_,0.));

    // don't need
    DalphaK_last_iter_     = Teuchos::null;
    KbbInv_                = Teuchos::null;
    Kbd_                   = Teuchos::null;
    fbeta_                 = Teuchos::null;
  }
  else //von-Mises
  {
    LINALG::Matrix<5,5> tmp55(true);
    LINALG::Matrix<5,numdofperelement_> tmp5x(true);
    LINALG::Matrix<5,1> tmp51(true);
    KbbInv_                      = Teuchos::rcp( new std::vector<LINALG::Matrix<5,5> >(numgpt_,tmp55));
    Kbd_                         = Teuchos::rcp( new std::vector<LINALG::Matrix<5,numdofperelement_> >(numgpt_,tmp5x));
    fbeta_                       = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >(numgpt_,tmp51));
    DalphaK_last_iter_           = Teuchos::rcp( new std::vector<LINALG::Matrix<5,1> >(numgpt_,tmp51));

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

/*----------------------------------------------------------------------*
 | read relevant parameters from paramter list              seitz 01/14 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::ReadParameterList(Teuchos::RCP<Teuchos::ParameterList> plparams)
{

  cpl_=plparams->get<double>("SEMI_SMOOTH_CPL");
  stab_s_=plparams->get<double>("STABILIZATION_S");
  if (eastype_!=soh8p_easnone)
    plparams->get<int>("have_EAS")=1;

  return;
}


/*----------------------------------------------------------------------*
 | extrapolate stresses for hex8 elements (public)           seitz 12/13 |
 *----------------------------------------------------------------------*/
template<class so3_ele, DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<so3_ele,distype>::soh8_expol(
    LINALG::Matrix<numgpt_post,numstr_>& stresses,
    Epetra_MultiVector& expolstresses
    )
{
  if (distype!=DRT::Element::hex8)
    dserror("soh8_expol called from non-hex8 element");

  // reorder GP
  // as the GPs are ordered differently in the hex8 element and the GP repository
  LINALG::Matrix<numgpt_post,numstr_> stresses_ordered(false);
  for (int i=0; i<numstr_; ++i)
  {
    stresses_ordered(7,i) = stresses(0,i);
    stresses_ordered(6,i) = stresses(1,i);
    stresses_ordered(4,i) = stresses(2,i);
    stresses_ordered(5,i) = stresses(3,i);
    stresses_ordered(3,i) = stresses(4,i);
    stresses_ordered(2,i) = stresses(5,i);
    stresses_ordered(0,i) = stresses(6,i);
    stresses_ordered(1,i) = stresses(7,i);
  }

  // static variables, that are the same for every element
  static LINALG::Matrix<nen_,numgpt_post> expol;
  static bool isfilled;

  if (isfilled==false)
  {
    double sq3=sqrt(3.0);

      expol(0,0)=1.25+0.75*sq3;
      expol(0,1)=-0.25-0.25*sq3;
      expol(0,2)=-0.25+0.25*sq3;
      expol(0,3)=-0.25-0.25*sq3;
      expol(0,4)=-0.25-0.25*sq3;
      expol(0,5)=-0.25+0.25*sq3;
      expol(0,6)=1.25-0.75*sq3;
      expol(0,7)=-0.25+0.25*sq3;
      expol(1,1)=1.25+0.75*sq3;
      expol(1,2)=-0.25-0.25*sq3;
      expol(1,3)=-0.25+0.25*sq3;
      expol(1,4)=-0.25+0.25*sq3;
      expol(1,5)=-0.25-0.25*sq3;
      expol(1,6)=-0.25+0.25*sq3;
      expol(1,7)=1.25-0.75*sq3;
      expol(2,2)=1.25+0.75*sq3;
      expol(2,3)=-0.25-0.25*sq3;
      expol(2,4)=1.25-0.75*sq3;
      expol(2,5)=-0.25+0.25*sq3;
      expol(2,6)=-0.25-0.25*sq3;
      expol(2,7)=-0.25+0.25*sq3;
      expol(3,3)=1.25+0.75*sq3;
      expol(3,4)=-0.25+0.25*sq3;
      expol(3,5)=1.25-0.75*sq3;
      expol(3,6)=-0.25+0.25*sq3;
      expol(3,7)=-0.25-0.25*sq3;
      expol(4,4)=1.25+0.75*sq3;
      expol(4,5)=-0.25-0.25*sq3;
      expol(4,6)=-0.25+0.25*sq3;
      expol(4,7)=-0.25-0.25*sq3;
      expol(5,5)=1.25+0.75*sq3;
      expol(5,6)=-0.25-0.25*sq3;
      expol(5,7)=-0.25+0.25*sq3;
      expol(6,6)=1.25+0.75*sq3;
      expol(6,7)=-0.25-0.25*sq3;
      expol(7,7)=1.25+0.75*sq3;

      for (int i=0;i<NUMNOD_SOH8;++i)
      {
        for (int j=0;j<i;++j)
        {
          expol(i,j)=expol(j,i);
        }
      }

    isfilled = true;
  }

  LINALG::Matrix<nen_,numstr_> nodalstresses;
  nodalstresses.Multiply(expol, stresses_ordered);

  // "assembly" of extrapolated nodal stresses
  for (int i=0;i<nen_;++i)
  {
    int gid = so3_ele::NodeIds()[i];
    if (expolstresses.Map().MyGID(so3_ele::NodeIds()[i])) // rownode
    {
      int lid = expolstresses.Map().LID(gid);
      int myadjele = Nodes()[i]->NumElement();
      for (int j=0;j<6;j++)
        (*(expolstresses(j)))[lid] += nodalstresses(i,j)/myadjele;
    }
  }
}


/*----------------------------------------------------------------------*/
// include the file at the end of so3_ssn_plast.cpp
#include "so3_ssn_plast_fwd.hpp"
