/*!----------------------------------------------------------------------**##
\file so_tet10.cpp
\brief

<pre>
Maintainer: Jonas Biehler
            biehler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15276
</pre>

*----------------------------------------------------------------------*/

#include "so_tet10.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_fem_general/drt_utils_integration.H"


DRT::ELEMENTS::So_tet10Type DRT::ELEMENTS::So_tet10Type::instance_;


DRT::ParObject* DRT::ELEMENTS::So_tet10Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_tet10* object = new DRT::ELEMENTS::So_tet10(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="SOLIDT10" )
  {
    Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_tet10(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_tet10Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = rcp(new DRT::ELEMENTS::So_tet10(id,owner));
  return ele;
}


void DRT::ELEMENTS::So_tet10Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_tet10Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_tet10Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDT10"];

  defs["TET10"]
    .AddIntVector("TET10",10)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddOptionalNamedDoubleVector("RAD",3)
    .AddOptionalNamedDoubleVector("AXI",3)
    .AddOptionalNamedDoubleVector("CIR",3)
    .AddOptionalNamedDoubleVector("FIBER1",3)
    .AddOptionalNamedDoubleVector("FIBER2",3)
    .AddOptionalNamedDoubleVector("FIBER3",3)
    ;
}


/*----------------------------------------------------------------------***
 |  ctor (public)                                              			|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet10::So_tet10(int id, int owner) :
DRT::Element(id,owner),
data_()
{
  kintype_ = so_tet10_nonlinear;
  invJ_.resize(NUMGPT_SOTET10, LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10>(true));
  detJ_.resize(NUMGPT_SOTET10, 0.0);
  invJ_mass_.resize(NUMGPT_MASS_SOTET10, LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10>(true));
  detJ_mass_.resize(NUMGPT_MASS_SOTET10, 0.0);
  return;
}

/*----------------------------------------------------------------------***
 |  copy-ctor (public)                                         			|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet10::So_tet10(const DRT::ELEMENTS::So_tet10& old) :
DRT::Element(old),
kintype_(old.kintype_),
data_(old.data_),
detJ_(old.detJ_),
detJ_mass_(old.detJ_mass_)
//try out later detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (int i=0; i<(int)invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }

  invJ_mass_.resize(old.invJ_mass_.size());
  for (int i=0; i<(int)invJ_mass_.size(); ++i)
    {
      invJ_mass_[i] = old.invJ_mass_[i];
    }

  return;
}

/*----------------------------------------------------------------------***
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            		    |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_tet10::Clone() const
{
  DRT::ELEMENTS::So_tet10* newelement = new DRT::ELEMENTS::So_tet10(*this);
  return newelement;
}

/*----------------------------------------------------------------------***
 |                                                             (public) |
 |                                                            		    |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_tet10::Shape() const
{
  return tet10;
}

/*----------------------------------------------------------------------***
 |  Pack data                                                  (public) |
 |                                                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Pack(DRT::PackBuffer& data) const
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
  //data
  AddtoPack(data,data_);
  // detJ_
  AddtoPack(data,detJ_);
  AddtoPack(data,detJ_mass_);

  //invJ
  const int size = (int)invJ_.size();
  AddtoPack(data,size);
  for (int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  const int size_mass = (int)invJ_mass_.size();
    AddtoPack(data,size_mass);
    for (int i=0; i<size_mass; ++i)
      AddtoPack(data,invJ_mass_[i]);

  return;
}


/*----------------------------------------------------------------------***
 |  Unpack data                                                (public) |
 |                                                            			|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Unpack(const vector<char>& data)
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
  // kintype_
  kintype_ = static_cast<KinematicType>( ExtractInt(position,data) );
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  // detJ_
  ExtractfromPack(position,data,detJ_);
  ExtractfromPack(position,data,detJ_mass_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  int size_mass = 0;
    ExtractfromPack(position,data,size_mass);
    invJ_mass_.resize(size_mass, LINALG::Matrix<NUMDIM_SOTET10,NUMDIM_SOTET10>(true));
    for (int i=0; i<size_mass; ++i)
      ExtractfromPack(position,data,invJ_mass_[i]);

    if (position != data.size())
      dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}



/*----------------------------------------------------------------------***
 |  dtor (public)                                              			|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_tet10::~So_tet10()
{
  return;
}


/*----------------------------------------------------------------------***
 |  print this element (public)                             	 		|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::Print(ostream& os) const
{
  os << "So_tet10 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}



/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      popp 12/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_expol
(
    LINALG::Matrix<NUMGPT_SOTET10,NUMSTR_SOTET10>& stresses,
    Epetra_MultiVector& expolstresses
)
{
  static LINALG::Matrix<NUMNOD_SOTET10,NUMGPT_SOTET10> expol;
  static bool isfilled;
  LINALG::Matrix<NUMNOD_SOTET10,NUMSTR_SOTET10> nodalstresses;

  if (isfilled==true)
  {
    nodalstresses.Multiply(expol, stresses);
  }
  else
  {
    // get gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(DRT::UTILS::intrule_tet_4point);

    // loop over all nodes
    for (int ip=0; ip<NUMNOD_SOTET10; ++ip)
    {
      // gaussian coordinates
      const double e1 = intpoints.qxg[ip][0];
      const double e2 = intpoints.qxg[ip][1];
      const double e3 = intpoints.qxg[ip][2];

      // coordinates of node in the fictitious GP element
      double e1expol;
      double e2expol;
      double e3expol;

      if (e1!=0) e1expol = 1/e1;
      else       e1expol = 0;
      if (e2!=0) e2expol = 1/e2;
      else       e2expol = 0;
      if (e3!=0) e3expol = 1/e3;
      else       e3expol = 0;

      // shape functions for the extrapolated coordinates
      // (yes, we REALLY mean DiscretizationType tet10 here!)
      LINALG::Matrix<NUMNOD_SOTET10,1> funct;
      DRT::UTILS::shape_function_3D(funct,e1expol,e2expol,e3expol,tet10);

      // extrapolation matrix
      for (int i=0;i<NUMGPT_SOTET10;++i) expol(ip,i) = funct(i);
    }

    // do extrapolation
    nodalstresses.Multiply(expol, stresses);

    isfilled = true;
  }

  // "assembly" of extrapolated nodal stresses
  for (int i=0;i<NUMNOD_SOTET10;++i)
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
  /* 10-node tetrahedra node topology*/
  /*--------------------------------------------------------------------*/
  /* parameter coordinates (ksi1, ksi2, ksi3, ksi4) of nodes
   * of a common tetrahedron [-1,1]x[-1,1]x[-1,1]
   *  10-node hexahedron: node 0,1,...,9
   *
   * -----------------------
   *- this is the numbering used in GiD & EXODUS!!
   *      3-
   *      |\ ---
   *      |  \    --9
   *      |    \      ---
   *      |      \        -2
   *      |        \       /\
   *      |          \   /   \
   *      7            8      \
   *      |          /   \     \
   *      |        6       \    5
   *      |      /           \   \
   *      |    /               \  \
   *      |  /                   \ \
   *      |/                       \\
   *      0------------4-------------1
   */
  /*====================================================================*/

/*----------------------------------------------------------------------***
 |  get vector of volumes (length 1) (public)                  			|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_tet10::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}


 /*----------------------------------------------------------------------**#
 |  get vector of surfaces (public)                             		|
 |  surface normals always point outward                                |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_tet10::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------***++
 |  get vector of lines (public)                               			|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_tet10::Lines()
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
void DRT::ELEMENTS::So_tet10::VisNames(map<string,int>& names)
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

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                          st 01/10|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_tet10::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name,data))
    return true;

  if (Material()->MaterialType() == INPAR::MAT::m_holzapfelcardiovascular){
    MAT::HolzapfelCardio* art = static_cast <MAT::HolzapfelCardio*>(Material().get());
    vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
    vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
    if (name == "Fiber1"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a1[0]; data[1] = a1[1]; data[2] = a1[2];
    } else if (name == "Fiber2"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a2[0]; data[1] = a2[1]; data[2] = a2[2];
    } else {
      return false;
    }
  }
  if (Material()->MaterialType() == INPAR::MAT::m_humphreycardiovascular){
    MAT::HumphreyCardio* art = static_cast <MAT::HumphreyCardio*>(Material().get());
    vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
    vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
    vector<double> a3 = art->Geta3()->at(0);  // get a3 of first gp
    vector<double> a4 = art->Geta4()->at(0);  // get a4 of first gp
    if (name == "Fiber1"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a1[0]; data[1] = a1[1]; data[2] = a1[2];
    } else if (name == "Fiber2"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a2[0]; data[1] = a2[1]; data[2] = a2[2];
    } else if (name == "Fiber3"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a3[0]; data[1] = a3[1]; data[2] = a3[2];
    } else if (name == "Fiber4"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a4[0]; data[1] = a4[1]; data[2] = a4[2];
    } else {
      return false;
    }
  }

  return true;
}


