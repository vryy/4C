/*!----------------------------------------------------------------------
\file bele3_4.cpp
\brief

<pre>
 * Maintainer: Benedikt Schott
 *             schott@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15241
</pre>

*----------------------------------------------------------------------*/

#include "bele3_4.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"


DRT::ELEMENTS::Bele3_4Type DRT::ELEMENTS::Bele3_4Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Bele3_4Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Bele3_4* object = new DRT::ELEMENTS::Bele3_4(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Bele3_4Type::Create( const std::string eletype,
                                                             const std::string eledistype,
                                                             const int id,
                                                             const int owner )
{
  if ( eletype=="BELE3_4" )
  {
    Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Bele3_4(id,owner));
    return ele;
  }
  return Teuchos::null;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Bele3_4Type::Create( const int id, const int owner )
{
  Teuchos::RCP<DRT::Element> ele = Teuchos::rcp(new DRT::ELEMENTS::Bele3_4(id,owner));
  return ele;
}


void DRT::ELEMENTS::Bele3_4Type::NodalBlockInformation( DRT::Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
}

void DRT::ELEMENTS::Bele3_4Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3_4::Bele3_4(int id, int owner) :
DRT::Element(id,owner)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3_4::Bele3_4(const DRT::ELEMENTS::Bele3_4& old) :
DRT::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Bele3_4::Clone() const
{
  DRT::ELEMENTS::Bele3_4* newelement = new DRT::ELEMENTS::Bele3_4(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Bele3_4::Shape() const
{
  switch (NumNode())
  {
  case  3: return tri3;
  case  4: return quad4;
  case  6: return tri6;
  case  8: return quad8;
  case  9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3_4::Pack(DRT::PackBuffer& data) const
{
  DRT::PackBuffer::SizeMarker sm( data );
  sm.Insert();

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  Element::Pack(data);

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3_4::Unpack(const std::vector<char>& data)
{
  std::vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  std::vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",data.size(),position);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Bele3_4::~Bele3_4()
{
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Bele3_4::Print(ostream& os) const
{
  os << "Bele3_4 " << DRT::DistypeToString(Shape());
  Element::Print(os);
  return;
}



DRT::UTILS::GaussRule2D DRT::ELEMENTS::Bele3_4::getOptimalGaussrule(const DRT::Element::DiscretizationType& distype) const
{
  DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
    switch (distype)
    {
    case DRT::Element::quad4:
        rule = DRT::UTILS::intrule_quad_4point;
        break;
    case DRT::Element::quad8: case DRT::Element::quad9:
        rule = DRT::UTILS::intrule_quad_9point;
        break;
    case DRT::Element::tri3:
        rule = DRT::UTILS::intrule_tri_3point;
        break;
    case DRT::Element::tri6:
        rule = DRT::UTILS::intrule_tri_6point;
        break;
    default:
        dserror("unknown number of nodes for gaussrule initialization");
  }
  return rule;
}



