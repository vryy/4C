/*!----------------------------------------------------------------------**##
\file so_tet4.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
writen by : Alexander Volf
			alexander.volf@mytum.de
</pre>

*----------------------------------------------------------------------*/

#include "so_tet4.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"


DRT::ELEMENTS::So_tet4Type DRT::ELEMENTS::So_tet4Type::instance_;


DRT::ParObject* DRT::ELEMENTS::So_tet4Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_tet4* object = new DRT::ELEMENTS::So_tet4(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDT4" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_tet4(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet4Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_tet4(id,owner));
  return ele;
}


void DRT::ELEMENTS::So_tet4Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_tet4Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_tet4Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT4"];

  defs["TET4"]
    .AddIntVector("TET4",4)
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


/*----------------------------------------------------------------------***
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4::So_tet4(int id, int owner) :
DRT::Element(id,owner),
//material_(0),
//data_(),
V_(-1.0),
pstype_(INPAR::STR::prestress_none),
pstime_(0.0),
time_(0.0)
{

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params!=Teuchos::null)
  {
    const ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
    pstime_ = sdyn.get<double>("PRESTRESSTIME");
  }

  if (pstype_==INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOTET4,NUMGPT_SOTET4,true));

  if (pstype_==INPAR::STR::prestress_id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOTET4,NUMGPT_SOTET4,true));

  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4::So_tet4(const DRT::ELEMENTS::So_tet4& old) :
DRT::Element(old),
//material_(old.material_),
//data_(old.data_),
V_(old.V_),
pstype_(old.pstype_),
pstime_(old.pstime_),
time_(old.time_)
{

  if (pstype_==INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));

  if (pstype_==INPAR::STR::prestress_id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(*(old.invdesign_)));

  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_tet4::Clone() const
{
  DRT::ELEMENTS::So_tet4* newelement = new DRT::ELEMENTS::So_tet4(*this);
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_tet4::Shape() const
{
  return tet4;
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);
  // ngp_
  //AddtoPack(data,ngp_,3*sizeof(int));
  // material_
  //AddtoPack(data,material_);
  // kintype_
  AddtoPack(data,kintype_);

  //vector<char> tmp(0);
  //data_.Pack(tmp);
  //AddtoPack(data,tmp);

  // V_
  AddtoPack(data,V_);

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

  return;
}


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // ngp_
  //ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // material_
  //ExtractfromPack(position,data,material_);
  // kintype_
  kintype_ = static_cast<KinematicType>( ExtractInt(position,data) );
  // data_
  //vector<char> tmp(0);
  //ExtractfromPack(position,data,tmp);
  //data_.Unpack(tmp);
  // V_
  ExtractfromPack(position,data,V_);

  // prestress_
  pstype_ = static_cast<INPAR::STR::PreStress>( ExtractInt(position,data) );
  ExtractfromPack(position,data,pstime_);
  ExtractfromPack(position,data,time_);
  if (pstype_==INPAR::STR::prestress_mulf)
  {
    vector<char> tmpprestress(0);
    ExtractfromPack(position,data,tmpprestress);
    if (prestress_ == Teuchos::null)
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOTET4,NUMGPT_SOTET4,true));
    prestress_->Unpack(tmpprestress);
  }

  // invdesign_
  if (pstype_==INPAR::STR::prestress_id)
  {
    vector<char> tmpinvdesign(0);
    ExtractfromPack(position,data,tmpinvdesign);
    if (invdesign_ == Teuchos::null)
      invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOTET4,NUMGPT_SOTET4,true));
    invdesign_->Unpack(tmpinvdesign);
  }

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 03/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_expol
(
    LINALG::Matrix<NUMGPT_SOTET4,NUMSTR_SOTET4>& stresses,
    Epetra_MultiVector& expolstresses
)
{
  static LINALG::Matrix<NUMNOD_SOTET4, NUMGPT_SOTET4> expol;
  static bool isfilled;

  if (isfilled==false)
  {
    expol(0,0)=1.0;
    expol(1,0)=1.0;
    expol(2,0)=1.0;
    expol(3,0)=1.0;

    isfilled=true;
  }

  LINALG::Matrix<NUMNOD_SOTET4,NUMSTR_SOTET4> nodalstresses;
  nodalstresses.Multiply(expol,stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i=0;i<NUMNOD_SOTET4;++i)
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
  return;
}




/*----------------------------------------------------------------------***
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet4::~So_tet4()
{
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::Print(ostream& os) const
{
  os << "So_tet4 ";
  Element::Print(os);
  cout << endl;
  //cout << data_;
  return;
}

  /*====================================================================*/
  /* 4-node tetrahedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (ksi1, ksi2, ksi3, ksi4) of nodes
   * of a common tetrahedron [-1,1]x[-1,1]x[-1,1]
   *  4-node hexahedron: node 0,1,...,3
   *
   * -----------------------
   *- this is the numbering used in GiD & EXODUS!!
   *      3-
   *      |\ ---
   *      |  \    ---
   *      |    \      ---
   *      |      \        -2
   *      |        \       /\
   *      |          \   /   \
   *      |            X      \
   *      |          /   \     \
   *      |        /       \    \
   *      |      /           \   \
   *      |    /               \  \
   *      |  /                   \ \
   *      |/                       \\
   *      0--------------------------1
   */
  /*====================================================================*/

/*----------------------------------------------------------------------***
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_tet4::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}


 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_tet4::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

vector<double> DRT::ELEMENTS::So_tet4::ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_SOTET4,NUMDIM_SOTET4> xrefe;  // material coord. of element
  for (int i=0; i<NUMNOD_SOTET4; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOTET4,1> funct;
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  LINALG::Matrix<1,NUMDIM_SOTET4> midpoint;
  //midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  vector<double> centercoords(3);
  centercoords[0] = midpoint(0,0);
  centercoords[1] = midpoint(0,1);
  centercoords[2] = midpoint(0,2);
  return centercoords;
}

/*----------------------------------------------------------------------***++
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_tet4::Lines()
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
 |  Return names of visualization data (public)                 st 01/10|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::VisNames(std::map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);
  if ((Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular))
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
        ostringstream s;
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
    } else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
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
          ostringstream s;
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
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                          st 01/10|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet4::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name,data))
    return true;

  if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular){
    MAT::HolzapfelCardio* art = static_cast <MAT::HolzapfelCardio*>(Material().get());
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
      data[0] = a1[0];
      data[1] = a1[1];
      data[2] = a1[2];
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
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
        ostringstream s;
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
  if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
    MAT::HumphreyCardio* art = static_cast <MAT::HumphreyCardio*>(Material().get());
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
      data[0] = a1[0];
      data[1] = a1[1];
      data[2] = a1[2];
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
      data[0] = a2[0];
      data[1] = a2[1];
      data[2] = a2[2];
    }
    else if (name == "Fiber3")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      vector<double> a3 = art->Geta3()->at(0);  // get a3 of first gp
      data[0] = a3[0];
      data[1] = a3[1];
      data[2] = a3[2];
    }
    else if (name == "Fiber4")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      vector<double> a4 = art->Geta4()->at(0);  // get a4 of first gp
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
      for (int iter=0; iter<NUMGPT_SOTET4; iter++)
        temp += grow->Gettheta()->at(iter);
      data[0] = temp/NUMGPT_SOTET4;
    }
    else if (name == "Mandel")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<NUMGPT_SOTET4; iter++)
        temp += grow->Getmandel()->at(iter);
      data[0] = temp/NUMGPT_SOTET4;
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular)
    {
      MAT::HolzapfelCardio* art = static_cast <MAT::HolzapfelCardio*>(grow->Matelastic().get());
      if (name == "Fiber1")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
        data[0] = a1[0];
        data[1] = a1[1];
        data[2] = a1[2];
      }
      else if (name == "Fiber2")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
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
        vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
        data[0] = a1[0];
        data[1] = a1[1];
        data[2] = a1[2];
      }
      else if (name == "Fiber2")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
        data[0] = a2[0];
        data[1] = a2[1];
        data[2] = a2[2];
      }
      else if (name == "Fiber3")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        vector<double> a3 = art->Geta3()->at(0);  // get a3 of first gp
        data[0] = a3[0];
        data[1] = a3[1];
        data[2] = a3[2];
      }
      else if (name == "Fiber4")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        vector<double> a4 = art->Geta4()->at(0);  // get a4 of first gp
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
          ostringstream s;
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
  if (Material()->MaterialType() == INPAR::MAT::m_constraintmixture){
    MAT::ConstraintMixture* cons = static_cast <MAT::ConstraintMixture*>(Material().get());
    if (name == "MassStress")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOTET4; iter++)
        temp.Update(1.0,cons->GetVis(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOTET4;
      data[1] = temp(1)/NUMGPT_SOTET4;
      data[2] = temp(2)/NUMGPT_SOTET4;
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
      for (int iter=0; iter<NUMGPT_SOTET4; iter++)
        temp += cons->GetMassDensity(iter);
      data[0] = temp/NUMGPT_SOTET4;
    }
    else if (name == "CollagenMassDensity")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOTET4; iter++)
        temp.Update(1.0,cons->GetMassDensityCollagen(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOTET4;
      data[1] = temp(1)/NUMGPT_SOTET4;
      data[2] = temp(2)/NUMGPT_SOTET4;
    }
    else
    {
      return false;
    }
  }

  return true;
}

