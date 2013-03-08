/*!----------------------------------------------------------------------
\file so_weg6.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/

#include "so_weg6.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"

DRT::ELEMENTS::So_weg6Type DRT::ELEMENTS::So_weg6Type::instance_;


DRT::ParObject* DRT::ELEMENTS::So_weg6Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_weg6* object =
    new DRT::ELEMENTS::So_weg6(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDW6" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_weg6(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_weg6Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_weg6(id,owner));
  return ele;
}


void DRT::ELEMENTS::So_weg6Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_weg6Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_weg6Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDW6"];

  defs["WEDGE6"]
    .AddIntVector("WEDGE6",6)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddOptionalNamedDoubleVector("RAD",3)
    .AddOptionalNamedDoubleVector("AXI",3)
    .AddOptionalNamedDoubleVector("CIR",3)
    .AddOptionalNamedDoubleVector("FIBER1",3)
    .AddOptionalNamedDoubleVector("FIBER2",3)
    .AddOptionalNamedDoubleVector("FIBER3",3)
    .AddOptionalNamedDouble("HU")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::So_weg6(int id, int owner) :
DRT::Element(id,owner),
data_(),
pstype_(INPAR::STR::prestress_none),
pstime_(0.0),
time_(0.0)
{
  kintype_ = sow6_nonlinear;
  invJ_.resize(NUMGPT_WEG6);
  detJ_.resize(NUMGPT_WEG6);
  for (int i=0; i<NUMGPT_WEG6; ++i)
  {
    detJ_[i]= 0.0;
    invJ_[i]=LINALG::Matrix<NUMDIM_WEG6,NUMDIM_WEG6>(true);
  }

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params!=Teuchos::null)
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
    pstime_ = sdyn.get<double>("PRESTRESSTIME");
  }

  if (pstype_==INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_WEG6,NUMGPT_WEG6));

  if (pstype_==INPAR::STR::prestress_id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_WEG6,NUMGPT_WEG6));

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::So_weg6(const DRT::ELEMENTS::So_weg6& old) :
DRT::Element(old),
kintype_(old.kintype_),
data_(old.data_),
detJ_(old.detJ_),
pstype_(old.pstype_),
pstime_(old.pstime_),
time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (unsigned int i=0; i<invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }

  if (pstype_==INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  if (pstype_==INPAR::STR::prestress_id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(*(old.invdesign_)));

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_weg6::Clone() const
{
  DRT::ELEMENTS::So_weg6* newelement = new DRT::ELEMENTS::So_weg6(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_weg6::Shape() const
{
  return wedge6;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  // kintype_
  AddtoPack(data,kintype_);
  // data_
  AddtoPack(data,data_);

  // prestress_
  AddtoPack(data,pstype_);
  AddtoPack(data,pstime_);
  AddtoPack(data,time_);
  if (pstype_==INPAR::STR::prestress_mulf)
  {
    DRT::ParObject::AddtoPack(data,*prestress_);
  }

  // invdesign_
  if (pstype_==INPAR::STR::prestress_id)
  {
    DRT::ParObject::AddtoPack(data,*invdesign_);
  }

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const unsigned int size = invJ_.size();
  AddtoPack(data,size);
  for (unsigned int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::Unpack(const std::vector<char>& data)
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
  kintype_ = static_cast<KinematicType>( ExtractInt(position,data) );
  // data_
  std::vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  // prestress_
  pstype_ = static_cast<INPAR::STR::PreStress>( ExtractInt(position,data) );
  ExtractfromPack(position,data,pstime_);
  ExtractfromPack(position,data,time_);
  if (pstype_==INPAR::STR::prestress_mulf)
  {
    std::vector<char> tmpprestress(0);
    ExtractfromPack(position,data,tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_WEG6,NUMGPT_WEG6));
    prestress_->Unpack(tmpprestress);
  }

  // invdesign_
  if (pstype_==INPAR::STR::prestress_id)
  {
    std::vector<char> tmpinvdesign(0);
    ExtractfromPack(position,data,tmpinvdesign);
    if (invdesign_ == Teuchos::null)
      invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_WEG6,NUMGPT_WEG6));
    invdesign_->Unpack(tmpinvdesign);
  }

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size);
  for (int i=0; i<size; ++i)
  {
    invJ_[i] = LINALG::Matrix<NUMDIM_WEG6,NUMDIM_WEG6>(true);
    ExtractfromPack(position,data,invJ_[i]);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::~So_weg6()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::Print(ostream& os) const
{
  os << "So_weg6 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes     maf 02/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::soweg6_expol
(
    LINALG::Matrix<NUMGPT_WEG6,NUMSTR_WEG6>& stresses,
    Epetra_MultiVector& expolstresses
)
{
  static LINALG::Matrix<NUMNOD_WEG6,NUMGPT_WEG6> expol;
  static bool isfilled;

  if (isfilled==false)
  {
   expol(0,0)=  -0.61004233964073;
   expol(0,1)=   0.12200846792815;
   expol(0,2)=   0.12200846792815;
   expol(0,3)=   2.27670900630740;
   expol(0,4)=  -0.45534180126148;
   expol(0,5)=  -0.45534180126148;
   expol(1,1)=  -0.61004233964073;
   expol(1,2)=   0.12200846792815;
   expol(1,3)=  -0.45534180126148;
   expol(1,4)=   2.27670900630740;
   expol(1,5)=  -0.45534180126148;
   expol(2,2)=  -0.61004233964073;
   expol(2,3)=  -0.45534180126148;
   expol(2,4)=  -0.45534180126148;
   expol(2,5)=   2.27670900630740;
   expol(3,3)=  -0.61004233964073;
   expol(3,4)=   0.12200846792815;
   expol(3,5)=   0.12200846792815;
   expol(4,4)=  -0.61004233964073;
   expol(4,5)=   0.12200846792815;
   expol(5,5)=  -0.61004233964073;
   for (int i=0;i<NUMNOD_WEG6;++i)
    {
      for (int j=0;j<i;++j)
      {
        expol(i,j)=expol(j,i);
      }
    }
   isfilled = true;
  }

  LINALG::Matrix<NUMNOD_WEG6,NUMSTR_WEG6> nodalstresses;
  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i=0;i<NUMNOD_WEG6;++i)
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
std::vector<double> DRT::ELEMENTS::So_weg6::ElementCenterRefeCoords()
{

  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_WEG6,NUMDIM_WEG6> xrefe;  // material coord. of element
  for (int i=0; i<NUMNOD_WEG6; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_WEG6,1> funct;
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  LINALG::Matrix<1,NUMDIM_WEG6> midpoint;
  //midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0,0);
  centercoords[1] = midpoint(0,1);
  centercoords[2] = midpoint(0,2);
  return centercoords;
}
/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::VisNames(std::map<string,int>& names)
{

  if ((Material()->MaterialType() == INPAR::MAT::m_artwallremod) ||
	  (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular))
  {
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
  }
  if (Material()->MaterialType() == INPAR::MAT::m_elasthyper)
  {
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(Material().get());
    if (elahy->AnisotropicPrincipal() or elahy->AnisotropicModified())
    {
      std::vector<LINALG::Matrix<3,1> > fibervecs;
      elahy->GetFiberVecs(fibervecs);
      int vissize = fibervecs.size();
      string fiber;
      for (int i = 0; i < vissize; i++)
      {
        std::ostringstream s;
        s << "Fiber" << i+1;
        fiber = s.str();
        names[fiber] = 3; // 3-dim vector
      }
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular)
  {
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3;
    fiber = "Fiber3";
    names[fiber] = 3;
    fiber = "Fiber4";
    names[fiber] = 3;
  }
  if (Material()->MaterialType() == INPAR::MAT::m_growth)
  {
    string fiber = "Theta";
    names[fiber] = 1;
    fiber = "Mandel";
    names[fiber] = 1;
    MAT::Growth* grow = static_cast <MAT::Growth*>(Material().get());
    if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular)
    {
      fiber = "Fiber1";
      names[fiber] = 3; // 3-dim vector
      fiber = "Fiber2";
      names[fiber] = 3; // 3-dim vector
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_humphreycardiovascular)
    {
      fiber = "Fiber1";
      names[fiber] = 3; // 3-dim vector
      fiber = "Fiber2";
      names[fiber] = 3;
      fiber = "Fiber3";
      names[fiber] = 3;
      fiber = "Fiber4";
      names[fiber] = 3;
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_elasthyper)
    {
      MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(grow->Matelastic().get());
      if (elahy->AnisotropicPrincipal() or elahy->AnisotropicModified())
      {
        std::vector<LINALG::Matrix<3,1> > fibervecs;
        elahy->GetFiberVecs(fibervecs);
        int vissize = fibervecs.size();
        string fiber;
        for (int i = 0; i < vissize; i++)
        {
          std::ostringstream s;
          s << "Fiber" << i+1;
          fiber = s.str();
          names[fiber] = 3; // 3-dim vector
        }
      }
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture)
  {
    string fiber = "MassStress";
    names[fiber] = 3;
    fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
    fiber = "referentialMassDensity";
    names[fiber] = 1;
    fiber = "CollagenMassDensity";
    names[fiber] = 3;
    fiber = "Prestretch";
    names[fiber] = 3;
    fiber = "Homstress";
    names[fiber] = 3;
  }
  if (Material()->MaterialType() == INPAR::MAT::m_aaaneohooke_stopro)
  {
    string fiber = "beta";
    names[fiber] = 1; // scalar
    fiber = "youngs";
    names[fiber] = 1; // scalar
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 07/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_weg6::VisData(const string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if(DRT::Element::VisData(name,data))
    return true;

  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod)
  {
    MAT::ArtWallRemod* art = static_cast <MAT::ArtWallRemod*>(Material().get());
    std::vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
    std::vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      data[0] = a1[0];
      data[1] = a1[1];
      data[2] = a1[2];
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      data[0] = a2[0];
      data[1] = a2[1];
      data[2] = a2[2];
    }
    else if (name == "Owner")
    {
      if ((int)data.size()<1)
        dserror("Size mismatch");
      data[0] = Owner();
    }
    else
    {
      cout << name << endl;
      dserror("Unknown VisData!");
    }

  }
  if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular)
  {
    MAT::HolzapfelCardio* art = static_cast <MAT::HolzapfelCardio*>(Material().get());
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
      data[0] = a1[0];
      data[1] = a1[1];
      data[2] = a1[2];
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
      data[0] = a2[0];
      data[1] = a2[1];
      data[2] = a2[2];
    }
    else
    {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_elasthyper)
  {
    MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(Material().get());
    if (elahy->AnisotropicPrincipal() or elahy->AnisotropicModified())
    {
      std::vector<LINALG::Matrix<3,1> > fibervecs;
      elahy->GetFiberVecs(fibervecs);
      int vissize = fibervecs.size();
      for (int i = 0; i < vissize; i++)
      {
        std::ostringstream s;
        s << "Fiber" << i+1;
        string fiber;
        fiber = s.str();
        if (name == fiber)
        {
          if ((int)data.size()!=3)
            dserror("size mismatch");
          data[0] = fibervecs.at(i)(0);
          data[1] = fibervecs.at(i)(1);
          data[2] = fibervecs.at(i)(2);
        }
      }
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular)
  {
    MAT::HumphreyCardio* art = static_cast <MAT::HumphreyCardio*>(Material().get());
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
      data[0] = a1[0];
      data[1] = a1[1];
      data[2] = a1[2];
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
      data[0] = a2[0];
      data[1] = a2[1];
      data[2] = a2[2];
    }
    else if (name == "Fiber3")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a3 = art->Geta3()->at(0);  // get a3 of first gp
      data[0] = a3[0];
      data[1] = a3[1];
      data[2] = a3[2];
    }
    else if (name == "Fiber4")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a4 = art->Geta4()->at(0);  // get a4 of first gp
      data[0] = a4[0];
      data[1] = a4[1];
      data[2] = a4[2];
    }
    else
    {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_growth)
  {
    MAT::Growth* grow = static_cast <MAT::Growth*>(Material().get());
    if (name == "Theta")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp += grow->Gettheta()->at(iter);
      data[0] = temp/NUMGPT_WEG6;
    }
    else if (name == "Mandel")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp += grow->Getmandel()->at(iter);
      data[0] = temp/NUMGPT_WEG6;
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular)
    {
      MAT::HolzapfelCardio* art = static_cast <MAT::HolzapfelCardio*>(grow->Matelastic().get());
      if (name == "Fiber1")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
        data[0] = a1[0];
        data[1] = a1[1];
        data[2] = a1[2];
      }
      else if (name == "Fiber2")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
        data[0] = a2[0];
        data[1] = a2[1];
        data[2] = a2[2];
      }
      else
      {
        return false;
      }
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_humphreycardiovascular)
    {
      MAT::HumphreyCardio* art = static_cast <MAT::HumphreyCardio*>(grow->Matelastic().get());
      if (name == "Fiber1")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
        data[0] = a1[0];
        data[1] = a1[1];
        data[2] = a1[2];
      }
      else if (name == "Fiber2")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
        data[0] = a2[0];
        data[1] = a2[1];
        data[2] = a2[2];
      }
      else if (name == "Fiber3")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a3 = art->Geta3()->at(0);  // get a3 of first gp
        data[0] = a3[0];
        data[1] = a3[1];
        data[2] = a3[2];
      }
      else if (name == "Fiber4")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a4 = art->Geta4()->at(0);  // get a4 of first gp
        data[0] = a4[0];
        data[1] = a4[1];
        data[2] = a4[2];
      }
      else
      {
        return false;
      }
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_elasthyper)
    {
      MAT::ElastHyper* elahy = static_cast <MAT::ElastHyper*>(grow->Matelastic().get());
      if (elahy->AnisotropicPrincipal() or elahy->AnisotropicModified())
      {
        std::vector<LINALG::Matrix<3,1> > fibervecs;
        elahy->GetFiberVecs(fibervecs);
        int vissize = fibervecs.size();
        for (int i = 0; i < vissize; i++)
        {
          std::ostringstream s;
          s << "Fiber" << i+1;
          string fiber;
          fiber = s.str();
          if (name == fiber)
          {
            if ((int)data.size()!=3)
              dserror("size mismatch");
            data[0] = fibervecs.at(i)(0);
            data[1] = fibervecs.at(i)(1);
            data[2] = fibervecs.at(i)(2);
          }
        }
      }
    }
    else
    {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture)
  {
    MAT::ConstraintMixture* cons = static_cast <MAT::ConstraintMixture*>(Material().get());
    if (name == "MassStress")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp.Update(1.0,cons->GetVis(iter),1.0);
      data[0] = temp(0)/NUMGPT_WEG6;
      data[1] = temp(1)/NUMGPT_WEG6;
      data[2] = temp(2)/NUMGPT_WEG6;
    }
    else if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> a1 = cons->Geta1()->at(0);  // get a1 of first gp
      data[0] = a1(0);
      data[1] = a1(1);
      data[2] = a1(2);
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> a2 = cons->Geta2()->at(0);  // get a2 of first gp
      data[0] = a2(0);
      data[1] = a2(1);
      data[2] = a2(2);
    }
    else if (name == "referentialMassDensity")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp += cons->GetMassDensity(iter);
      data[0] = temp/NUMGPT_WEG6;
    }
    else if (name == "CollagenMassDensity")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp.Update(1.0,cons->GetMassDensityCollagen(iter),1.0);
      data[0] = temp(0)/NUMGPT_WEG6;
      data[1] = temp(1)/NUMGPT_WEG6;
      data[2] = temp(2)/NUMGPT_WEG6;
    }
    else if (name == "Prestretch")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp.Update(1.0,cons->GetPrestretch(iter),1.0);
      data[0] = temp(0)/NUMGPT_WEG6;
      data[1] = temp(1)/NUMGPT_WEG6;
      data[2] = temp(2)/NUMGPT_WEG6;
    }
    else if (name == "Homstress")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_WEG6; iter++)
        temp.Update(1.0,cons->GetHomstress(iter),1.0);
      data[0] = temp(0)/NUMGPT_WEG6;
      data[1] = temp(1)/NUMGPT_WEG6;
      data[2] = temp(2)/NUMGPT_WEG6;
    }
    else
    {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_aaaneohooke_stopro)
  {
    MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(Material().get());
    if (name=="beta")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      data[0] = aaa_stopro->Beta();
    }
    else if (name=="youngs")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      data[0] = aaa_stopro->Youngs();
    }
    else
    {
      return false;
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_weg6::Volumes()
{
  std::vector<Teuchos::RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_weg6::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_weg6::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
}


