/*!----------------------------------------------------------------------
\file so_sh18.cpp
\brief

<pre>
Maintainer: Alexander Seitz
            seitz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>

*----------------------------------------------------------------------*/

#include "so_sh18.H"
#include "so_surface.H"
#include "so_line.H"
#include "../drt_lib/drt_utils_factory.H"
#include "../drt_lib/drt_utils_nullspace.H"
#include "../drt_mat/so3_material.H"

DRT::ELEMENTS::So_sh18Type DRT::ELEMENTS::So_sh18Type::instance_;
namespace {
  const std::string name = DRT::ELEMENTS::So_sh18Type::Instance().Name();
}

DRT::ParObject* DRT::ELEMENTS::So_sh18Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::So_sh18* object = new DRT::ELEMENTS::So_sh18(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh18Type::Create( const std::string eletype,
                                                            const std::string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if (eletype == "SOLIDSH18")
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh18(id,owner));
    return ele;
  }

  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::So_sh18Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::So_sh18(id,owner));
  return ele;
}


void DRT::ELEMENTS::So_sh18Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  numdf = 3;
  dimns = 6;
  nv = 3;
}

void DRT::ELEMENTS::So_sh18Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeStructure3DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::So_sh18Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["SOLIDSH18"];

  defs["HEX18"]
    .AddIntVector("HEX18",18)
    .AddNamedInt("MAT")
    .AddNamedString("KINEM")
    .AddNamedString("TSL")
    .AddNamedString("MEL")
    .AddNamedString("CTL")
    .AddNamedString("VOL")
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
 |  ctor (public)                                           seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh18::So_sh18(int id, int owner) :
DRT::Element(id,owner)
{
  invJ_.resize(NUMGPT_SOH18, LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18>(true));
  detJ_.resize(NUMGPT_SOH18, 0.0);

  // setup integration rule
  xsi_.resize(NUMGPT_SOH18);
  wgt_.resize(NUMGPT_SOH18);

  // in plane
  DRT::UTILS::IntPointsAndWeights<2> ip_p(DRT::UTILS::intrule_quad_9point);
  // director
  DRT::UTILS::IntPointsAndWeights<1> ip_d(DRT::UTILS::intrule_line_2point);

  for (int d=0; d<ip_d.IP().nquad; ++d)
    for (int p=0; p<ip_p.IP().nquad; ++p)
    {
      wgt_[p+d*ip_p.IP().nquad] = ip_d.IP().qwgt[d]*
                                  ip_p.IP().qwgt[p];
      xsi_.at(p+d*ip_p.IP().nquad)(0)=ip_p.IP().qxg[p][0];
      xsi_.at(p+d*ip_p.IP().nquad)(1)=ip_p.IP().qxg[p][1];
      xsi_.at(p+d*ip_p.IP().nquad)(2)=ip_d.IP().qxg[d][0];
    }

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                      seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh18::So_sh18(const DRT::ELEMENTS::So_sh18& old) :
DRT::Element(old),
detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (int i=0; i<(int)invJ_.size(); ++i)
  {
    // can this size be anything but NUMDIM_SOH8 x NUMDIM_SOH8?
    //invJ_[i].Shape(old.invJ_[i].M(),old.invJ_[i].N());
    invJ_[i] = old.invJ_[i];
  }

  // setup integration rule
  xsi_.resize(NUMGPT_SOH18);
  wgt_.resize(NUMGPT_SOH18);

  // in plane
  DRT::UTILS::IntPointsAndWeights<2> ip_p(DRT::UTILS::intrule_quad_9point);
  // director
  DRT::UTILS::IntPointsAndWeights<1> ip_d(DRT::UTILS::intrule_line_2point);

  for (int d=0; d<ip_d.IP().nquad; ++d)
    for (int p=0; p<ip_p.IP().nquad; ++p)
    {
      wgt_[p+d*ip_p.IP().nquad] = ip_d.IP().qwgt[d]*
                                  ip_p.IP().qwgt[p];
      xsi_.at(p+d*ip_p.IP().nquad)(0)=ip_p.IP().qxg[p][0];
      xsi_.at(p+d*ip_p.IP().nquad)(1)=ip_p.IP().qxg[p][1];
      xsi_.at(p+d*ip_p.IP().nquad)(2)=ip_d.IP().qxg[d][0];
    }

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_sh18::Clone() const
{
  DRT::ELEMENTS::So_sh18* newelement = new DRT::ELEMENTS::So_sh18(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_sh18::Shape() const
{
  return hex18;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Pack(DRT::PackBuffer& data) const
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

  // element technology bools
  AddtoPack(data,(int)dsg_shear_);
  AddtoPack(data,(int)dsg_membrane_);
  AddtoPack(data,(int)dsg_ctl_);
  AddtoPack(data,(int)eas_);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Unpack(const std::vector<char>& data)
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

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size = 0;
  ExtractfromPack(position,data,size);
  invJ_.resize(size, LINALG::Matrix<NUMDIM_SOH18,NUMDIM_SOH18>(true));
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  // element technology bools
  dsg_shear_    = ExtractInt(position,data);
  dsg_membrane_ = ExtractInt(position,data);
  dsg_ctl_      = ExtractInt(position,data);
  eas_          = ExtractInt(position,data);
  SetupDSG();

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            seitz 11/14 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_sh18::~So_sh18()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                             seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::Print(std::ostream& os) const
{
  os << "So_sh18 ";
  Element::Print(os);
  return;
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)               seitz 11/14 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_sh18::Volumes()
{
  std::vector<Teuchos::RCP<Element> > volumes(1);
  volumes[0]= Teuchos::rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          seitz 11/14 |
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_sh18::Surfaces()
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
 |  get vector of lines (public)                            seitz 11/14 |
 *----------------------------------------------------------------------*/
std::vector<Teuchos::RCP<DRT::Element> > DRT::ELEMENTS::So_sh18::Lines()
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
 |  Return names of visualization data (public)             seitz 11/14 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_sh18::VisNames(std::map<std::string,int>& names)
{
  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
  so3mat->VisNames(names);

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                      seitz 11/14 |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_sh18::VisData(const std::string& name, std::vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if (DRT::Element::VisData(name,data))
    return true;

  Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());

  return so3mat->VisData(name, data, NUMGPT_SOH18, this->Id());
}


