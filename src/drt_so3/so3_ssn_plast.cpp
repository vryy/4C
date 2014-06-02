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
#include "so_surface.H"
#include "so_line.H"


/*----------------------------------------------------------------------*
 | ctor (public)                                            seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Plast<distype>::So3_Plast(
  int id,
  int owner
  )
: DRT::Element(id,owner),
  KbbInv_(std::vector<Epetra_SerialDenseMatrix>(0)),
  Kbd_(std::vector<Epetra_SerialDenseMatrix>(0)),
  fbeta_(std::vector<Epetra_SerialDenseVector>(0)),
  dDp_last_iter_(std::vector<Epetra_SerialDenseVector>(0)),
  plspintype_(plspin),
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
  return;
}


/*----------------------------------------------------------------------*
 | copy-ctor (public)                                       seitz 07/13 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::So3_Plast<distype>::So3_Plast(
  const DRT::ELEMENTS::So3_Plast<distype>& old
  ):
   DRT::Element(old)
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


/*----------------------------------------------------------------------*
 |                                                          seitz 05/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::So3_Plast<distype>::NumVolume() const
{
  switch(distype)
  {
  case DRT::Element::hex8:
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
  Element::Pack(data);

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

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
  Element::Unpack(basedata);

  // kintype_
//  kintype_ = static_cast<GenKinematicType>( ExtractInt(position,data) );
  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<nsd_,nsd_>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

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

   // EAS element technology
   eastype_=static_cast<EASType>(ExtractInt(position,data));
   ExtractfromPack(position,data,neas_);
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
     Kba_                                =Teuchos::rcp(new std::vector<Epetra_SerialDenseMatrix>(numgpt_,Epetra_SerialDenseMatrix(plspintype_,neas_)));
     alpha_eas_                          =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
     alpha_eas_last_timestep_            =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
     alpha_eas_delta_over_last_timestep_ =Teuchos::rcp(new Epetra_SerialDenseVector(neas_));
   }

     KbbInv_        .resize(numgpt_,Epetra_SerialDenseMatrix(plspintype_,plspintype_));
     Kbd_           .resize(numgpt_,Epetra_SerialDenseMatrix(plspintype_,numdofperelement_));
     fbeta_         .resize(numgpt_,Epetra_SerialDenseVector(plspintype_));
     dDp_last_iter_ .resize(numgpt_,Epetra_SerialDenseVector(plspintype_));

   if (eastype_!=soh8p_easnone)
   {
     ExtractfromPack(position,data,(*alpha_eas_));
     ExtractfromPack(position,data,(*alpha_eas_last_timestep_));
     ExtractfromPack(position,data,(*alpha_eas_delta_over_last_timestep_));
   }

   size=ExtractInt(position,data);
   for (int i=0; i<size; i++)
     ExtractfromPack(position,data,dDp_last_iter_[i]);

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

  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);

  SetMaterial(material);

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->Setup(numgpt_, linedef);
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
  KbbInv_        .resize(numgpt_,Epetra_SerialDenseMatrix(plspintype_,plspintype_));
  Kbd_           .resize(numgpt_,Epetra_SerialDenseMatrix(plspintype_,numdofperelement_));
  fbeta_         .resize(numgpt_,Epetra_SerialDenseVector(plspintype_));
  dDp_last_iter_ .resize(numgpt_,Epetra_SerialDenseVector(plspintype_));

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
  std::string accumulatedstrain = "accumulatedstrain";
  names[accumulatedstrain] = 1; // scalar
  std::string plastic_zone = "plastic_zone";
  names[plastic_zone] = 1; // scalar

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

  // get plastic hyperelastic material
     MAT::PlasticElastHyper* plmat = NULL;
  if (Material()->MaterialType()==INPAR::MAT::m_plelasthyper)
    plmat= static_cast<MAT::PlasticElastHyper*>(Material().get());
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  if (name == "accumulatedstrain")
  {
    if ((int)data.size()!=1) dserror("size mismatch");
    double tmp=0.;
    for (int gp=0; gp<numgpt_; gp++)
      tmp+= plmat->AccumulatedStrain(gp);
    data[0] = tmp/numgpt_;
  }
  if (name == "plastic_zone")
  {
    bool plastic_history=false;
    bool curr_active=false;
    if ((int)data.size()!=1) dserror("size mismatch");
    for (int gp=0; gp<numgpt_; gp++)
    {
      if (plmat->AccumulatedStrain(gp)!=0.) plastic_history=true;
      if (plmat->Active(gp))                curr_active=true;
    }
    data[0] = plastic_history+curr_active;
  }

  return false;

}  // VisData()

/*----------------------------------------------------------------------*
 | read relevant parameters from paramter list              seitz 01/14 |
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::So3_Plast<distype>::ReadParameterList(Teuchos::RCP<Teuchos::ParameterList> plparams)
{
  double cpl=plparams->get<double>("SEMI_SMOOTH_CPL");
  double s=plparams->get<double>("STABILIZATION_S");
  static_cast<MAT::PlasticElastHyper*>(Material().get())->GetParams(s,cpl);
  if (eastype_!=soh8p_easnone)
    plparams->get<int>("have_EAS")=1;

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

  // reorder GP
  // as the GPs are ordered differently in the hex8 element and the GP repository
  LINALG::Matrix<numgpt_post,numstr_> stresses_ordered(false);
  for (int i=0; i<numstr_; ++i)
  {
    stresses_ordered(6,i) = stresses(0,i);
    stresses_ordered(7,i) = stresses(1,i);
    stresses_ordered(5,i) = stresses(2,i);
    stresses_ordered(4,i) = stresses(3,i);
    stresses_ordered(2,i) = stresses(4,i);
    stresses_ordered(3,i) = stresses(5,i);
    stresses_ordered(1,i) = stresses(6,i);
    stresses_ordered(0,i) = stresses(7,i);
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
    int gid = NodeIds()[i];
    if (expolstresses.Map().MyGID(NodeIds()[i])) // rownode
    {
      int lid = expolstresses.Map().LID(gid);
      int myadjele = Nodes()[i]->NumElement();
      for (int j=0;j<6;j++)
        (*(expolstresses(j)))[lid] += nodalstresses(i,j)/myadjele;
    }
  }
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
  else
    dserror("so3_ssn_plast elements only with PlasticElastHyper material");

  return plmat->HavePlasticSpin();
}

#include "so3_ssn_plast_fwd.hpp"
