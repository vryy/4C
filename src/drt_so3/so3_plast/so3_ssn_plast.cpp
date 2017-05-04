/*----------------------------------------------------------------------*/
/*!
\file so3_ssn_plast.cpp
\brief
\maintainer Alexander Seitz
\level 2

*/


/*----------------------------------------------------------------------*
 | headers                                                  seitz 07/13 |
 *----------------------------------------------------------------------*/
#include "so3_ssn_plast.H"

#include "../../drt_lib/drt_linedefinition.H"
#include "../../drt_fem_general/drt_utils_shapefunctions_service.H"
#include "../../drt_lib/drt_utils_factory.H"
#include "../../drt_mat/plasticelasthyper.H"
#include "../so_surface.H"
#include "../so_line.H"
#include "../../drt_inpar/inpar_tsi.H"
#include "../../linalg/linalg_serialdensevector.H"

// include this thermo-implementation to make sure, that the same Gauss-rule
// is used in the structural and the thermal part in a TSI problem
#include "../../drt_thermo/thermo_ele_impl_utils.H"
#include "../../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 | ctor (public)                                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Plast<distype>::So3_Plast(
  int id,
  int owner
  )
: So_base(id,owner),
  fbar_(false),
//  data_(Teuchos::null),
  KbbInv_(std::vector<LINALG::SerialDenseMatrix>(0)),
  Kbd_(std::vector<LINALG::SerialDenseMatrix>(0)),
  fbeta_(std::vector<LINALG::SerialDenseVector>(0)),
  dDp_last_iter_(std::vector<LINALG::SerialDenseVector>(0)),
  dDp_inc_(std::vector<LINALG::SerialDenseVector>(0)),
  plspintype_(plspin),
  KaaInv_(Teuchos::null),
  Kad_(Teuchos::null),
  KaT_(Teuchos::null),
  KdT_eas_(Teuchos::null),
  feas_(Teuchos::null),
  Kba_(Teuchos::null),
  alpha_eas_(Teuchos::null),
  alpha_eas_last_timestep_(Teuchos::null),
  alpha_eas_delta_over_last_timestep_(Teuchos::null),
  alpha_eas_inc_(Teuchos::null),
  eastype_(soh8p_easnone),
  neas_(0),
  tsi_(false),
  is_nitsche_contact_(false)
{
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Plast<distype>::So3_Plast(
  const DRT::ELEMENTS::So3_Plast<distype>& old
  ):
  So_base(old)
{
  return;
}


/*----------------------------------------------------------------------*
 | deep copy this instance of Solid3 and return pointer to  seitz 07/13 |
 | it (public)                                                          |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::Element* DRT::ELEMENTS::So3_Plast<distype>::Clone() const
{
  DRT::ELEMENTS::So3_Plast<distype>* newelement
    = new DRT::ELEMENTS::So3_Plast<distype>(*this);

  return newelement;
}


template<DRT::Element::DiscretizationType distype>
const int DRT::ELEMENTS::So3_Plast<distype>::VOIGT3X3SYM_[3][3] = {{0,3,5},{3,1,4},{5,4,2}};
template<DRT::Element::DiscretizationType distype>
const int DRT::ELEMENTS::So3_Plast<distype>::VOIGT3X3NONSYM_[3][3] = {{0,3,5},{6,1,4},{8,7,2}};

template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nen_,1> >
DRT::ELEMENTS::So3_Plast<distype>::shapefunct_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nen_> >
DRT::ELEMENTS::So3_Plast<distype>::deriv_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::invJ_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,double>
DRT::ELEMENTS::So3_Plast<distype>::detJ_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nen_> >
DRT::ELEMENTS::So3_Plast<distype>::N_XYZ_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::defgrd_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::defgrd_mod_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::rcg_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::delta_Lp_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::numstr_,DRT::ELEMENTS::So3_Plast<distype>::numdofperelement_> >
DRT::ELEMENTS::So3_Plast<distype>::bop_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::numstr_,1> >
DRT::ELEMENTS::So3_Plast<distype>::pk2_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::numstr_,DRT::ELEMENTS::So3_Plast<distype>::numstr_> >
DRT::ELEMENTS::So3_Plast<distype>::cmat_;

template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nen_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::xrefe_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nen_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::xcurr_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nen_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::xcurr_rate_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nen_,1> >
DRT::ELEMENTS::So3_Plast<distype>::etemp_;

template<DRT::Element::DiscretizationType distype>
std::pair<bool,double>
DRT::ELEMENTS::So3_Plast<distype>::detF_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,double>
DRT::ELEMENTS::So3_Plast<distype>::detF_0_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::inv_defgrd_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::inv_defgrd_0_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nen_> >
DRT::ELEMENTS::So3_Plast<distype>::N_XYZ_0_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::numstr_,1> >
DRT::ELEMENTS::So3_Plast<distype>::rcg_vec_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,double>
DRT::ELEMENTS::So3_Plast<distype>::f_bar_fac_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::numdofperelement_,1> >
DRT::ELEMENTS::So3_Plast<distype>::htensor_;

template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::numstr_,DRT::ELEMENTS::So3_Plast<distype>::numstr_> >
DRT::ELEMENTS::So3_Plast<distype>::T0invT_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::Matrix<DRT::ELEMENTS::So3_Plast<distype>::nsd_,DRT::ELEMENTS::So3_Plast<distype>::nsd_> >
DRT::ELEMENTS::So3_Plast<distype>::jac_0_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,double>
DRT::ELEMENTS::So3_Plast<distype>::det_jac_0_;
template<DRT::Element::DiscretizationType distype>
std::pair<bool,LINALG::SerialDenseMatrix>
DRT::ELEMENTS::So3_Plast<distype>::M_eas_;

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::NumVolume() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex18:
  case DRT::Element::hex27:
    return 0;
    break;
  default:
    dserror("unknown distpye for So3_Plast");
    break;
    return 0;
  }
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::NumSurface() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex18:
  case DRT::Element::hex27:
    return 6;
    break;
  default:
    dserror("unknown distpye for So3_Plast");
    break;
    return 0;
  }
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::NumLine() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  case DRT::Element::hex18:
  case DRT::Element::hex27:
    return 12;
    break;
  default:
    dserror("unknown distpye for So3_Plast");
    break;
    return 0;
  }
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So3_Plast<distype>::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So3_Plast<distype>::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So3_Plast<distype>::Volumes()
{
  std::vector<Teuchos::RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
 | pack data (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::Pack(
  DRT::PackBuffer& data
  ) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);

  // add base class Element
  So_base::Pack(data);

  // Gauss points and weights
  const int size2 = (int)xsi_.size();
  AddtoPack(data,size2);
  for (int i=0; i<size2;++i)
    AddtoPack(data,xsi_[i]);
  AddtoPack(data,wgt_);

  // parameters
  AddtoPack(data,(int)fbar_);

  // plastic spin type
  AddtoPack(data,(int)plspintype_);

  // tsi
  AddtoPack(data,(int)tsi_);
  if (tsi_)
  {
    AddtoPack(data,(int)KbT_->size());
    for (unsigned i=0; i<KbT_->size() ; i++)
    {
      AddtoPack(data,(*dFintdT_)[i]);
      AddtoPack(data,(*KbT_)[i]);
      AddtoPack(data,(*temp_last_)[i]);
    }
  }

  // EAS element technology
  AddtoPack(data,(int)eastype_);
  AddtoPack(data,neas_);
  if (eastype_!=soh8p_easnone)
  {
    AddtoPack(data,(*alpha_eas_));
    AddtoPack(data,(*alpha_eas_last_timestep_));
    AddtoPack(data,(*alpha_eas_delta_over_last_timestep_));
  }

  // history at each Gauss point
  int histsize=dDp_last_iter_.size();
  AddtoPack(data,histsize);
  if (histsize!=0)
    for (int i=0; i<histsize; i++)
      AddtoPack(data,dDp_last_iter_[i]);

  // nitsche contact
  AddtoPack(data,(int)is_nitsche_contact_);
  if (is_nitsche_contact_)
  {
    AddtoPack(data,cauchy_);
    AddtoPack(data,cauchy_deriv_);
    if (tsi_)
      AddtoPack(data,cauchy_deriv_T_);
  }

  return;
}  // Pack()


/*----------------------------------------------------------------------*
 | unpack data (public)                                     seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::Unpack(
  const std::vector<char>& data
  )
{
  std::vector<char>::size_type position = 0;

  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");

  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  So_base::Unpack(basedata);

  // Gauss points and weights
  int size2 = ExtractInt(position,data);
  xsi_.resize(size2, LINALG::Matrix<nsd_,1>(true));
  for (int i=0; i<size2;++i)
    ExtractfromPack(position,data,xsi_[i]);
  ExtractfromPack(position,data,wgt_);
  numgpt_=wgt_.size();

  // paramters
  fbar_=(bool)ExtractInt(position,data);

   // plastic spin type
   plspintype_=static_cast<PlSpinType>(ExtractInt(position,data));

   //tsi
   tsi_=(bool)ExtractInt(position,data);
   if (tsi_)
   {
     dFintdT_=Teuchos::rcp(new std::vector<LINALG::Matrix<numdofperelement_,1> >(numgpt_));
     KbT_=Teuchos::rcp(new std::vector<LINALG::SerialDenseVector>
       (numgpt_,LINALG::SerialDenseVector(plspintype_,true)));
     temp_last_=Teuchos::rcp(new std::vector<double>(numgpt_));
     int size=ExtractInt(position,data);
     for (int i=0; i<size; i++)
     {
       ExtractfromPack(position,data,(*dFintdT_)[i]);
       ExtractfromPack(position,data,(*KbT_)[i]);
       ExtractfromPack(position,data,(*temp_last_)[i]);
     }
   }

   // EAS element technology
   eastype_=static_cast<DRT::ELEMENTS::So3Plast_EASType>(ExtractInt(position,data));
   ExtractfromPack(position,data,neas_);

   // no EAS
   if (eastype_==soh8p_easnone)
   {
     KaaInv_                             =Teuchos::null;
     Kad_                                =Teuchos::null;
     KaT_                                =Teuchos::null;
     KdT_eas_                            =Teuchos::null;
     feas_                               =Teuchos::null;
     Kba_                                =Teuchos::null;
     alpha_eas_                          =Teuchos::null;
     alpha_eas_last_timestep_            =Teuchos::null;
     alpha_eas_delta_over_last_timestep_ =Teuchos::null;
     alpha_eas_inc_                      =Teuchos::null;
   }
   else
   {
     KaaInv_                             =Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,neas_,true));
     Kad_                                =Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,numdofperelement_,true));
     if (tsi_)
     {
       KaT_                      =Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,nen_,true));
       KdT_eas_                  =Teuchos::rcp(new LINALG::Matrix<numdofperelement_,nen_>);
     }
     feas_                               =Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
     Kba_                                =Teuchos::rcp(new std::vector<LINALG::SerialDenseMatrix>(numgpt_,LINALG::SerialDenseMatrix(plspintype_,neas_,true)));
     alpha_eas_                          =Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
     alpha_eas_last_timestep_            =Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
     alpha_eas_delta_over_last_timestep_ =Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
     alpha_eas_inc_                      =Teuchos::rcp(new LINALG::SerialDenseVector(neas_,true));
   }

   KbbInv_        .resize(numgpt_,LINALG::SerialDenseMatrix(plspintype_,plspintype_,true));
   Kbd_           .resize(numgpt_,LINALG::SerialDenseMatrix(plspintype_,numdofperelement_,true));
   fbeta_         .resize(numgpt_,LINALG::SerialDenseVector(plspintype_,true));
   dDp_last_iter_ .resize(numgpt_,LINALG::SerialDenseVector(plspintype_,true));
   dDp_inc_       .resize(numgpt_,LINALG::SerialDenseVector(plspintype_,true));

   if (eastype_!=soh8p_easnone)
   {
     ExtractfromPack(position,data,(*alpha_eas_));
     ExtractfromPack(position,data,(*alpha_eas_last_timestep_));
     ExtractfromPack(position,data,(*alpha_eas_delta_over_last_timestep_));
   }

   int size=ExtractInt(position,data);
   for (int i=0; i<size; i++)
     ExtractfromPack(position,data,dDp_last_iter_[i]);

   // Nitsche contact stuff
   is_nitsche_contact_=(bool)ExtractInt(position,data);
   if (is_nitsche_contact_)
   {
     ExtractfromPack(position,data,cauchy_);
     ExtractfromPack(position,data,cauchy_deriv_);
     if (tsi_)
       ExtractfromPack(position,data,cauchy_deriv_T_);
     else
       cauchy_deriv_T_.resize(0);
   }
   else
   {
     cauchy_      .resize(0);
     cauchy_deriv_.resize(0);
     cauchy_deriv_T_.resize(0);
   }

   if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;

}  // Unpack()


/*----------------------------------------------------------------------*
 | print this element (public)                              seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::Print(std::ostream& os) const
{
  os << "So3_Plast ";
  return;
}


/*----------------------------------------------------------------------*
 | read this element, get the material (public)             seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Plast<distype>::ReadElement(
  const std::string& eletype,
  const std::string& eledistype,
  DRT::INPUT::LineDefinition* linedef
  )
{
  std::string buffer;
  linedef->ExtractString("KINEM",buffer);

  // geometrically linear
  if (buffer == "linear")
  {
    dserror("no linear kinematics");
  }
  // geometrically non-linear with Total Lagrangean approach
  else if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::kinem_nonlinearTotLag;
    // everything ok
  }
  else
    dserror("Reading of SO3_PLAST element failed! KINEM unknown");

  // fbar
  if (linedef->HaveNamed("FBAR"))
  {
    std::string fb;
    linedef->ExtractString("FBAR",fb);
    if (fb=="yes")
      fbar_=true;
    else if (fb=="no")
      fbar_=false;
    else
      dserror("unknown fbar option (valid: yes/no)");
  }

  // quadrature
  if (linedef->HaveNamed("NUMGP"))
  {
    if (distype!=DRT::Element::hex8)
      dserror("You may only choose the Gauss point number for SOLIDH8PLAST");
    if (DRT::Problem::Instance()->ProblemType() == prb_tsi)
      dserror("You may not choose the Gauss point number in TSI problems");

    int ngp =0;
    linedef->ExtractInt("NUMGP",ngp);

    switch(ngp)
    {
    case 8:
    {
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::UTILS::intrule_hex_8point);
      numgpt_=intpoints.IP().nquad;
      xsi_.resize(numgpt_);
      wgt_.resize(numgpt_);
      for (int gp=0; gp<numgpt_; ++gp)
      {
        wgt_[gp]=(intpoints.IP().qwgt)[gp];
        const double* gpcoord = (intpoints.IP().qxg)[gp];
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
      DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::UTILS::intrule_hex_27point);
      numgpt_=intpoints.IP().nquad;
      xsi_.resize(numgpt_);
      wgt_.resize(numgpt_);
      for (int gp=0; gp<numgpt_; ++gp)
      {
        wgt_[gp]=(intpoints.IP().qwgt)[gp];
        const double* gpcoord = (intpoints.IP().qxg)[gp];
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
    DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(THR::DisTypeToOptGaussRule<distype>::rule);
    numgpt_=intpoints.IP().nquad;
    xsi_.resize(numgpt_);
    wgt_.resize(numgpt_);
    for (int gp=0; gp<numgpt_; ++gp)
    {
      wgt_[gp]=(intpoints.IP().qwgt)[gp];
      const double* gpcoord = (intpoints.IP().qxg)[gp];
      for (int idim=0; idim<nsd_; idim++)
        xsi_[gp](idim) = gpcoord[idim];
    }
  }

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);

  SetMaterial(material);

  Teuchos::RCP<MAT::So3Material> so3mat = SolidMaterial();
  so3mat->Setup(numgpt_, linedef);
  so3mat->ValidKinematics(INPAR::STR::kinem_nonlinearTotLag);
  if (so3mat->MaterialType()!=INPAR::MAT::m_plelasthyper)
    std::cout << "*** warning *** so3plast used w/o PlasticElastHyper material. Better use standard solid element!\n";
  if (HavePlasticSpin())
    plspintype_=plspin;
  else
    plspintype_=zerospin;

  // EAS
  if (linedef->HaveNamed("EAS"))
  {
    if (distype != DRT::Element::hex8)
      dserror("EAS in so3 plast currently only for HEX8 elements");

    linedef->ExtractString("EAS",buffer);

    if (buffer == "none")
      eastype_=soh8p_easnone;
    else if (buffer == "mild")
      eastype_=soh8p_easmild;
    else if (buffer=="full")
      eastype_=soh8p_easfull;
    else
      dserror("unknown EAS type for so3_plast");

    if (fbar_ && eastype_!=soh8p_easnone)
      dserror("no combination of Fbar and EAS");
  }
  else
    eastype_=soh8p_easnone;

  // initialize EAS data
  EasInit();

  // plasticity related stuff
  KbbInv_        .resize(numgpt_,LINALG::SerialDenseMatrix(plspintype_,plspintype_,true));
  Kbd_           .resize(numgpt_,LINALG::SerialDenseMatrix(plspintype_,numdofperelement_,true));
  fbeta_         .resize(numgpt_,LINALG::SerialDenseVector(plspintype_,true));
  dDp_last_iter_ .resize(numgpt_,LINALG::SerialDenseVector(plspintype_,true));
  dDp_inc_       .resize(numgpt_,LINALG::SerialDenseVector(plspintype_,true));

  Teuchos::ParameterList plparams   =
      DRT::Problem::Instance()->SemiSmoothPlastParams();
  plparams.set<PROBLEM_TYP>("PROBLEM_TYP",DRT::Problem::Instance()->ProblemType());
  ReadParameterList(Teuchos::rcpFromRef<Teuchos::ParameterList>(plparams));


  return true;

}  // ReadElement()

/*----------------------------------------------------------------------*
 | get the nodes from so3 (public)                          seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::UniqueParObjectId() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return So_hex8PlastType::Instance().UniqueParObjectId();
    break;
  }  // hex8
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
template<DRT::Element::DiscretizationType distype>
DRT::ElementType& DRT::ELEMENTS::So3_Plast<distype>::ElementType() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return So_hex8PlastType::Instance();
    break;
  }
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
 | return names of visualization data (public)              seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::VisNames(std::map<std::string,int>& names)
{
  DRT::Element::VisNames(names);
  SolidMaterial()->VisNames(names);

  return;
}  // VisNames()

/*----------------------------------------------------------------------*
 | return visualization data (public)                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Plast<distype>::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name,data))
    return true;

  return SolidMaterial()->VisData(name,data,numgpt_,Id());

}  // VisData()

/*----------------------------------------------------------------------*
 | read relevant parameters from paramter list              seitz 01/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::ReadParameterList(Teuchos::RCP<Teuchos::ParameterList> plparams)
{
  double cpl=plparams->get<double>("SEMI_SMOOTH_CPL");
  double s=plparams->get<double>("STABILIZATION_S");
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    static_cast<MAT::PlasticElastHyper*>(Material().get())->GetParams(s,cpl);

  PROBLEM_TYP probtype = plparams->get<PROBLEM_TYP>("PROBLEM_TYP");
  if (probtype == prb_tsi)
    tsi_=true;
  else
    tsi_=false;
  if (tsi_)
  {
    // get plastic hyperelastic material
    MAT::PlasticElastHyper* plmat = NULL;
    if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
      plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
    else
      dserror("so3_ssn_plast elements only with PlasticElastHyper material");

    // get dissipation mode
    INPAR::TSI::DissipationMode mode =
        DRT::INPUT::IntegralValue<INPAR::TSI::DissipationMode>(*plparams,"DISSIPATION_MODE");

    // prepare material for tsi
    plmat->SetupTSI(numgpt_,numdofperelement_,(eastype_!=soh8p_easnone),mode);

    // setup element data
    dFintdT_=Teuchos::rcp(new std::vector<LINALG::Matrix<numdofperelement_,1> >(numgpt_));
    temp_last_=Teuchos::rcp(new std::vector<double>(numgpt_,plmat->InitTemp()));
    KbT_=Teuchos::rcp(new std::vector<LINALG::SerialDenseVector>
    (numgpt_,LINALG::SerialDenseVector(plspintype_,true)));

    if (eastype_!=soh8p_easnone)
    {
      KaT_     = Teuchos::rcp(new LINALG::SerialDenseMatrix(neas_,nen_,true));
      KdT_eas_ = Teuchos::rcp(new LINALG::Matrix<numdofperelement_,nen_>);
    }
    else
    {
      KaT_     = Teuchos::null;
      KdT_eas_ = Teuchos::null;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | extrapolate stresses for hex8 elements (public)           seitz 12/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::soh8_expol(
    LINALG::Matrix<numgpt_post,numstr_>& stresses,
    Epetra_MultiVector& expolstresses
    )
{
  if (distype!=DRT::Element::hex8)
    dserror("soh8_expol called from non-hex8 element");

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
  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i=0; i<nen_; ++i)
  {
    const int lid = expolstresses.Map().LID(NodeIds()[i]);
    if (lid >= 0) // rownode
    {
      const double invmyadjele = 1.0/Nodes()[i]->NumElement();
      for (int j=0; j<MAT::NUM_STRESS_3D; ++j)
        (*(expolstresses(j)))[lid] += nodalstresses(i,j)*invmyadjele;
    }
  }
  return;
}


/*----------------------------------------------------------------------*
 | Have plastic spin                                        seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
bool DRT::ELEMENTS::So3_Plast<distype>::HavePlasticSpin()
{
  // get plastic hyperelastic material
  MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());

  if (plmat!=NULL)
    return plmat->HavePlasticSpin();

  return false;
}

int DRT::ELEMENTS::PlastEasTypeToNumEasV(DRT::ELEMENTS::So3Plast_EASType et)
{
  switch (et)
  {
  case soh8p_easnone: return PlastEasTypeToNumEas<soh8p_easnone>::neas; break;
  case soh8p_easmild: return PlastEasTypeToNumEas<soh8p_easmild>::neas; break;
  case soh8p_easfull: return PlastEasTypeToNumEas<soh8p_easfull>::neas; break;
  case soh8p_eassosh8: return PlastEasTypeToNumEas<soh8p_eassosh8>::neas; break;
  case soh18p_eassosh18: return PlastEasTypeToNumEas<soh18p_eassosh18>::neas; break;
  default: dserror("EAS type not implemented");
  }
  return -1;
}


#include "so3_ssn_plast_fwd.hpp"
