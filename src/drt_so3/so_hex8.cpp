/*!----------------------------------------------------------------------
\file so_hex8.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/

#include "so_hex8.H"
#include "so_hex8fbar.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/contchainnetw.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/structporo.H"
#include "../drt_mat/plasticlinelast.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/damage.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>

// inverse design object
#include "inversedesign.H"
#include "prestress.H"


DRT::ELEMENTS::So_hex8Type DRT::ELEMENTS::So_hex8Type::instance_;
namespace {
  const std::string name = DRT::ELEMENTS::So_hex8Type::Instance().Name();
}

DRT::ParObject* DRT::ELEMENTS::So_hex8Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_hex8* object = new DRT::ELEMENTS::So_hex8(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if (eletype == "SOLIDH8")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex8(id,owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_hex8Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_hex8(id,owner));
  return ele;
}


void DRT::ELEMENTS::So_hex8Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_hex8Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_hex8Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDH8"];

  defs["HEX8"]
    .AddIntVector("HEX8",8)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddNamedString("EAS")
    .AddOptionalNamedDoubleVector("RAD",3)
    .AddOptionalNamedDoubleVector("AXI",3)
    .AddOptionalNamedDoubleVector("CIR",3)
    .AddOptionalNamedDoubleVector("FIBER1",3)
    .AddOptionalNamedDoubleVector("FIBER2",3)
    .AddOptionalNamedDoubleVector("FIBER3",3)
    .AddOptionalNamedDouble("STRENGTH")
    .AddOptionalNamedDouble("HU")
    ;
}


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(int id, int owner) :
DRT::Element(id,owner),
data_(),
pstype_(INPAR::STR::prestress_none),
pstime_(0.0),
time_(0.0)
{
  kintype_ = soh8_nonlinear;
  eastype_ = soh8_easnone;
  neas_ = 0;
  invJ_.resize(NUMGPT_SOH8, LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>(true));
  detJ_.resize(NUMGPT_SOH8, 0.0);

  Teuchos::RCP<const Teuchos::ParameterList> params = DRT::Problem::Instance()->getParameterList();
  if (params!=Teuchos::null)
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    pstype_ = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
    pstime_ = sdyn.get<double>("PRESTRESSTIME");
  }
  if (pstype_==INPAR::STR::prestress_mulf)
    prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8,NUMGPT_SOH8));

  if (pstype_==INPAR::STR::prestress_id)
    invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOH8,NUMGPT_SOH8));

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::So_hex8(const DRT::ELEMENTS::So_hex8& old) :
DRT::Element(old),
kintype_(old.kintype_),
eastype_(old.eastype_),
neas_(old.neas_),
data_(old.data_),
detJ_(old.detJ_),
pstype_(old.pstype_),
pstime_(old.pstime_),
time_(old.time_)
{
  invJ_.resize(old.invJ_.size());
  for (int i=0; i<(int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH8 x NUMDIM_SOH8?
    //invJ_[i].Shape(old.invJ_[i].M(),old.invJ_[i].N());
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
DRT::Element* DRT::ELEMENTS::So_hex8::Clone() const
{
  DRT::ELEMENTS::So_hex8* newelement = new DRT::ELEMENTS::So_hex8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_hex8::Shape() const
{
  return hex8;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Pack(DRT::PackBuffer& data) const
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
  // eastype_
  AddtoPack(data,eastype_);
  // neas_
  AddtoPack(data,neas_);
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
  else if (pstype_==INPAR::STR::prestress_id)
  {
    DRT::ParObject::AddtoPack(data,*invdesign_);
  }

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Unpack(const std::vector<char>& data)
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
  // eastype_
  eastype_ = static_cast<EASType>( ExtractInt(position,data) );
  // neas_
  ExtractfromPack(position,data,neas_);
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
    {
      int numgpt = NUMGPT_SOH8;
      // see whether I am actually a So_hex8fbar element
      DRT::ELEMENTS::So_hex8fbar* me = dynamic_cast<DRT::ELEMENTS::So_hex8fbar*>(this);
      if (me) numgpt += 1; // one more history entry for centroid data in hex8fbar
      prestress_ = Teuchos::rcp(new DRT::ELEMENTS::PreStress(NUMNOD_SOH8,numgpt));
    }
    prestress_->Unpack(tmpprestress);
  }

  // invdesign_
  else if (pstype_==INPAR::STR::prestress_id)
  {
    std::vector<char> tmpinvdesign(0);
    ExtractfromPack(position,data,tmpinvdesign);
    if (invdesign_ == Teuchos::null)
      invdesign_ = Teuchos::rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_SOH8,NUMGPT_SOH8));
    invdesign_->Unpack(tmpinvdesign);
  }

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_hex8::~So_hex8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  Get number of degrees of freedom (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex8::NumDofPerNode(const unsigned nds,const DRT::Node& node) const
{
  if (nds==1)
  {
    // what's the current problem type?
    PROBLEM_TYP probtype = DRT::Problem::Instance()->ProblemType();
    switch (probtype)
    {
    // in case of TSI/TFSI_Aero: nds=1 (second discretization) temperature: 1Dof/Node
    case prb_tsi:
    case prb_tfsi_aero:
    case prb_thermo_fsi:
      return 1;
      break;
    default:
      dserror("Unknown problem type with a second dofset!");
    break;
    }
  }
  // structure, 3D: 3 Dofs per node
  //    return 3;
  return NumDofPerNode(node);
};


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::Print(ostream& os) const
{
  os << "So_hex8 ";
  Element::Print(os);
  //cout << endl;
  //cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      lw 02/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_expol
(
    LINALG::Matrix<NUMGPT_SOH8,NUMSTR_SOH8>& stresses,
    Epetra_MultiVector& expolstresses
)
{
  // static variables, that are the same for every element
  static LINALG::Matrix<NUMNOD_SOH8,NUMGPT_SOH8> expol;
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

  LINALG::Matrix<NUMNOD_SOH8,NUMSTR_SOH8> nodalstresses;
  nodalstresses.Multiply(expol, stresses);

  // "assembly" of extrapolated nodal stresses
  for (int i=0;i<NUMNOD_SOH8;++i)
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

  /*====================================================================*/
  /* 8-node hexhedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (r,s,t) of nodes
   * of biunit cube [-1,1]x[-1,1]x[-1,1]
   *  8-node hexahedron: node 0,1,...,7
   *                      t
   *                      |
   *             4========|================7
   *           //|        |               /||
   *          // |        |              //||
   *         //  |        |             // ||
   *        //   |        |            //  ||
   *       //    |        |           //   ||
   *      //     |        |          //    ||
   *     //      |        |         //     ||
   *     5=========================6       ||
   *    ||       |        |        ||      ||
   *    ||       |        o--------||---------s
   *    ||       |       /         ||      ||
   *    ||       0------/----------||------3
   *    ||      /      /           ||     //
   *    ||     /      /            ||    //
   *    ||    /      /             ||   //
   *    ||   /      /              ||  //
   *    ||  /      /               || //
   *    || /      r                ||//
   *    ||/                        ||/
   *     1=========================2
   *
   */
  /*====================================================================*/

/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Volumes()
{
  std::vector<Teuchos::RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Surfaces()
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
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_hex8::Lines()
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
 |  get location of element center                              jb 08/11|
 *----------------------------------------------------------------------*/
std::vector<double> DRT::ELEMENTS::So_hex8::ElementCenterRefeCoords()
{
  // update element geometry
  DRT::Node** nodes = Nodes();
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  for (int i=0; i<NUMNOD_SOH8; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  const DRT::Element::DiscretizationType distype = Shape();
  LINALG::Matrix<NUMNOD_SOH8,1> funct;
  // Element midpoint at r=s=t=0.0
  DRT::UTILS::shape_function_3D(funct, 0.0, 0.0, 0.0, distype);
  LINALG::Matrix<1,NUMDIM_SOH8> midpoint;
  //midpoint.Multiply('T','N',1.0,funct,xrefe,0.0);
  midpoint.MultiplyTN(funct, xrefe);
  std::vector<double> centercoords(3);
  centercoords[0] = midpoint(0,0);
  centercoords[1] = midpoint(0,1);
  centercoords[2] = midpoint(0,2);
  return centercoords;
}
/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 01/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::VisNames(std::map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);

  if (Material()->MaterialType() == INPAR::MAT::m_contchainnetw)
  {
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber3";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber4";
    names[fiber] = 3; // 3-dim vector
    fiber = "FiberCell1";
    names[fiber] = 3; // 3-dim vector
    fiber = "FiberCell2";
    names[fiber] = 3; // 3-dim vector
    fiber = "FiberCell3";
    names[fiber] = 3; // 3-dim vector
    fiber = "l1";
    names[fiber] = 1;
    fiber = "l2";
    names[fiber] = 1;
    fiber = "l3";
    names[fiber] = 1;
//    fiber = "l1_0";
//    names[fiber] = 1;
//    fiber = "l2_0";
//    names[fiber] = 1;
//    fiber = "l3_0";
//    names[fiber] = 1;
  }
  if ((Material()->MaterialType() == INPAR::MAT::m_artwallremod) ||
      (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic)||
      (Material()->MaterialType() == INPAR::MAT::m_anisotropic_balzani)||
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
  if (Material()->MaterialType() == INPAR::MAT::m_structporo)
  {
    string porosity = "porosity";
    names[porosity] = 1; // scalar
  }
  if ( (Material()->MaterialType() == INPAR::MAT::m_pllinelast) ||
       (Material()->MaterialType() == INPAR::MAT::m_thermopllinelast) ||
       (Material()->MaterialType() == INPAR::MAT::m_elpldamage)
       )
  {
    string accumulatedstrain = "accumulatedstrain";
    names[accumulatedstrain] = 1; // scalar
  }


  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex8::VisData(const string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name,data))
    return true;

  if (Material()->MaterialType() == INPAR::MAT::m_contchainnetw)
  {
    RCP<MAT::Material> mat = Material();
    MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
    if (!chain->Initialized())
    {
      data[0] = 0.0; data[1] = 0.0; data[2] = 0.0;
    }
    else
    {
      RCP<std::vector<std::vector<double> > > gplis = chain->Getli();
      RCP<std::vector<std::vector<double> > > gpli0s = chain->Getli0();
      RCP<std::vector<LINALG::Matrix<3,3> > > gpnis = chain->Getni();

      std::vector<double> centerli (3,0.0);
      std::vector<double> centerli_0 (3,0.0);
      for (int i = 0; i < (int)gplis->size(); ++i)
      {
        LINALG::Matrix<3,1> loc(&(gplis->at(i)[0]));
        //Epetra_SerialDenseVector loc(CV,&(gplis->at(i)[0]),3);
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,gpnis->at(i),loc,0.0);
        glo.Multiply(gpnis->at(i),loc);
        // Unfortunately gpnis is a vector of Epetras, to change this
        // I must begin at a deeper level...
        centerli[0] += glo(0);
        centerli[1] += glo(1);
        centerli[2] += glo(2);

//        centerli[0] += gplis->at(i)[0];
//        centerli[1] += gplis->at(i)[1];
//        centerli[2] += gplis->at(i)[2];
//
        centerli_0[0] += gplis->at(i)[0];
        centerli_0[1] += gplis->at(i)[1];
        centerli_0[2] += gplis->at(i)[2];
      }
      centerli[0] /= gplis->size();
      centerli[1] /= gplis->size();
      centerli[2] /= gplis->size();

      centerli_0[0] /= gplis->size();
      centerli_0[1] /= gplis->size();
      centerli_0[2] /= gplis->size();

      // just the unit cell of the first gp
      int gp = 0;
//      Epetra_SerialDenseVector loc(CV,&(gplis->at(gp)[0]),3);
//      Epetra_SerialDenseVector glo(3);
//      glo.Multiply('N','N',1.0,gpnis->at(gp),loc,0.0);
      LINALG::Matrix<3,3> T(gpnis->at(gp).A(),true);
      std::vector<double> gpli =  chain->Getli()->at(gp);

      if (name == "Fiber1")
      {
        if ((int)data.size()!=3) dserror("size mismatch");
        data[0] = centerli[0];
        data[1] = -centerli[1];
        data[2] = -centerli[2];
      }
      else if (name == "Fiber2")
      {
        data[0] = centerli[0];
        data[1] = centerli[1];
        data[2] = -centerli[2];
      }
      else if (name == "Fiber3")
      {
        data[0] = centerli[0];
        data[1] = centerli[1];
        data[2] = centerli[2];
      }
      else if (name == "Fiber4")
      {
        data[0] = -centerli[0];
        data[1] = -centerli[1];
        data[2] = centerli[2];
      }
      else if (name == "FiberCell1")
      {
        LINALG::Matrix<3,1> e(true);
        e(0) = gpli[0];
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,T,e,0.0);
        glo.Multiply(T, e);
        data[0] = glo(0);
        data[1] = glo(1);
        data[2] = glo(2);
      }
      else if (name == "FiberCell2")
      {
        LINALG::Matrix<3,1> e(true);
        e(1) = gpli[1];
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,T,e,0.0);
        glo.Multiply(T, e);
        data[0] = glo(0);
        data[1] = glo(1);
        data[2] = glo(2);
      }
      else if (name == "FiberCell3")
      {
        LINALG::Matrix<3,1> e(true);
        e(2) = gpli[2];
        LINALG::Matrix<3,1> glo;
        //glo.Multiply('N','N',1.0,T,e,0.0);
        glo.Multiply(T, e);
        data[0] = glo(0);
        data[1] = glo(1);
        data[2] = glo(2);
      }
      else if (name == "l1")
      {
        data[0] = centerli_0[0];
      }
      else if (name == "l2")
      {
        data[0] = centerli_0[1];
      }
      else if (name == "l3")
      {
        data[0] = centerli_0[2];
//      } else if (name == "l1_0"){
//        data[0] = centerli_0[0];
//      } else if (name == "l2_0"){
//        data[0] = centerli_0[1];
//      } else if (name == "l3_0"){
//        data[0] = centerli_0[2];
      }
      else
      {
        return false;
      }
    }
  }
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
    else
    {
     return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic)
  {
    MAT::ViscoAnisotropic* art = static_cast <MAT::ViscoAnisotropic*>(Material().get());
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
    else
    {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_anisotropic_balzani)
  {
    MAT::AnisotropicBalzani* balz = static_cast <MAT::AnisotropicBalzani*>(Material().get());
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      data[0] = balz->Geta1().at(0);
      data[1] = balz->Geta1().at(1);
      data[2] = balz->Geta1().at(2);
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      data[0] = balz->Geta2().at(0);
      data[1] = balz->Geta2().at(1);
      data[2] = balz->Geta2().at(2);
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
  if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular)
  {
    MAT::HumphreyCardio* art = static_cast <MAT::HumphreyCardio*>(Material().get());
    if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a1 = art->Geta1()->at(0); // get a1 of first gp
      data[0] = a1[0];
      data[1] = a1[1];
      data[2] = a1[2];
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a2 = art->Geta2()->at(0); // get a2 of first gp
      data[0] = a2[0];
      data[1] = a2[1];
      data[2] = a2[2];
    }
    else if (name == "Fiber3")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a3 = art->Geta3()->at(0); // get a3 of first gp
      data[0] = a3[0];
      data[1] = a3[1];
      data[2] = a3[2];
    }
    else if (name == "Fiber4")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      std::vector<double> a4 = art->Geta4()->at(0); // get a4 of first gp
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
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp += grow->Gettheta()->at(iter);
      data[0] = temp/NUMGPT_SOH8;
    }
    else if (name == "Mandel")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp += grow->Getmandel()->at(iter);
      data[0] = temp/NUMGPT_SOH8;
    }
    else if (grow->Matelastic()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular)
    {
      MAT::HolzapfelCardio* art = static_cast <MAT::HolzapfelCardio*>(grow->Matelastic().get());
      if (name == "Fiber1")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a1 = art->Geta1()->at(0); // get a1 of first gp
        data[0] = a1[0];
        data[1] = a1[1];
        data[2] = a1[2];
      }
      else if (name == "Fiber2")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a2 = art->Geta2()->at(0); // get a2 of first gp
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
        std::vector<double> a1 = art->Geta1()->at(0); // get a1 of first gp
        data[0] = a1[0];
        data[1] = a1[1];
        data[2] = a1[2];
      }
      else if (name == "Fiber2")
      {
        if ((int)data.size()!=3) dserror("size mismatch");
        std::vector<double> a2 = art->Geta2()->at(0); // get a2 of first gp
        data[0] = a2[0];
        data[1] = a2[1];
        data[2] = a2[2];
      }
      else if (name == "Fiber3")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a3 = art->Geta3()->at(0); // get a3 of first gp
        data[0] = a3[0];
        data[1] = a3[1];
        data[2] = a3[2];
      }
      else if (name == "Fiber4")
      {
        if ((int)data.size()!=3)
          dserror("size mismatch");
        std::vector<double> a4 = art->Geta4()->at(0); // get a4 of first gp
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
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,cons->GetVis(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
      data[1] = temp(1)/NUMGPT_SOH8;
      data[2] = temp(2)/NUMGPT_SOH8;
    }
    else if (name == "Fiber1")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> a1 = cons->Geta1()->at(0); // get a1 of first gp
      data[0] = a1(0);
      data[1] = a1(1);
      data[2] = a1(2);
    }
    else if (name == "Fiber2")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> a2 = cons->Geta2()->at(0); // get a2 of first gp
      data[0] = a2(0);
      data[1] = a2(1);
      data[2] = a2(2);
    }
    else if (name == "referentialMassDensity")
    {
      if ((int)data.size()!=1)
        dserror("size mismatch");
      double temp = 0.0;
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp += cons->GetMassDensity(iter);
      data[0] = temp/NUMGPT_SOH8;
    }
    else if (name == "CollagenMassDensity")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,cons->GetMassDensityCollagen(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
      data[1] = temp(1)/NUMGPT_SOH8;
      data[2] = temp(2)/NUMGPT_SOH8;
    }
    else if (name == "Prestretch")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,cons->GetPrestretch(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
      data[1] = temp(1)/NUMGPT_SOH8;
      data[2] = temp(2)/NUMGPT_SOH8;
    }
    else if (name == "Homstress")
    {
      if ((int)data.size()!=3)
        dserror("size mismatch");
      LINALG::Matrix<3,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,cons->GetHomstress(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
      data[1] = temp(1)/NUMGPT_SOH8;
      data[2] = temp(2)/NUMGPT_SOH8;
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
      if ((int)data.size()!=1) dserror("size mismatch");
      data[0] = aaa_stopro->Beta();
    }
    else if (name=="youngs")
    {
      if ((int)data.size()!=1) dserror("size mismatch");
      data[0] = aaa_stopro->Youngs();
    }
    else
    {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_structporo)
  {
    MAT::StructPoro* structporo = static_cast <MAT::StructPoro*>(Material().get());
      if (name=="porosity")
      {
        if ((int)data.size()!=1) dserror("size mismatch");
          data[0] = structporo->PorosityAv();
      }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_pllinelast)
  {
    MAT::PlasticLinElast* pllinelast = static_cast <MAT::PlasticLinElast*>(Material().get());
    if (name == "accumulatedstrain")
    {
      if ((int)data.size()!=1) dserror("size mismatch");
      LINALG::Matrix<1,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,pllinelast->AccumulatedStrain(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_thermopllinelast)
  {
    MAT::ThermoPlasticLinElast* thrpllinelast = static_cast <MAT::ThermoPlasticLinElast*>(Material().get());
    if (name == "accumulatedstrain")
    {
      if ((int)data.size()!=1) dserror("size mismatch");
      LINALG::Matrix<1,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,thrpllinelast->AccumulatedStrain(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_elpldamage)
  {
    MAT::Damage* damage = static_cast <MAT::Damage*>(Material().get());
    if (name == "accumulatedstrain")
    {
      if ((int)data.size()!=1) dserror("size mismatch");
      LINALG::Matrix<1,1> temp(true);
      for (int iter=0; iter<NUMGPT_SOH8; iter++)
        temp.Update(1.0,damage->AccumulatedStrain(iter),1.0);
      data[0] = temp(0)/NUMGPT_SOH8;
    }
  }

  return true;
}


