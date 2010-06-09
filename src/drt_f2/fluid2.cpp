/*!----------------------------------------------------------------------
\file fluid2.cpp
\brief

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;


DRT::ELEMENTS::Fluid2Type DRT::ELEMENTS::Fluid2Type::instance_;


DRT::ParObject* DRT::ELEMENTS::Fluid2Type::Create( const std::vector<char> & data )
{
  DRT::ELEMENTS::Fluid2* object = new DRT::ELEMENTS::Fluid2(-1,-1);
  object->Unpack(data);
  return object;
}


Teuchos::RCP<DRT::Element> DRT::ELEMENTS::Fluid2Type::Create( const string eletype,
                                                            const string eledistype,
                                                            const int id,
                                                            const int owner )
{
  if ( eletype=="FLUID2" )
  {
    if ( eledistype!="NURBS4" and eledistype!="NURBS9" )
      return rcp(new DRT::ELEMENTS::Fluid2(id,owner));
  }
  return Teuchos::null;
}


void DRT::ELEMENTS::Fluid2Type::NodalBlockInformation( Element * dwele, int & numdf, int & dimns, int & nv, int & np )
{
  nv = 2;
  np = 1;
  numdf = 3;
  dimns = 3;
}


void DRT::ELEMENTS::Fluid2Type::ComputeNullSpace( DRT::Discretization & dis, std::vector<double> & ns, const double * x0, int numdf, int dimns )
{
  DRT::UTILS::ComputeFluid2DNullSpace( dis, ns, x0, numdf, dimns );
}

void DRT::ELEMENTS::Fluid2Type::SetupElementDefinition( std::map<std::string,std::map<std::string,DRT::INPUT::LineDefinition> > & definitions )
{
  std::map<std::string,DRT::INPUT::LineDefinition>& defs = definitions["FLUID2"];

  defs["QUAD4"]
    .AddIntVector("QUAD4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["QUAD8"]
    .AddIntVector("QUAD8",8)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["QUAD9"]
    .AddIntVector("QUAD9",9)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TRI3"]
    .AddIntVector("TRI3",3)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["TRI6"]
    .AddIntVector("TRI6",6)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS4"]
    .AddIntVector("NURBS4",4)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["NURBS9"]
    .AddIntVector("NURBS9",9)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;

  defs["THQ9"]
    .AddIntVector("THQ9",9)
    .AddNamedInt("MAT")
    .AddNamedString("NA")
    ;
}


map<string,DRT::ELEMENTS::Fluid2::StabilisationAction> DRT::ELEMENTS::Fluid2::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::Fluid2(int id, int owner) :
DRT::Element(id,owner),
is_ale_(false),
dismode_(dismod_equalorder),
data_()
{
  gaussrule_ = intrule2D_undefined;

  saccn_ .Shape(0,0);
  svelnp_.Shape(0,0);
  sveln_ .Shape(0,0);

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::Fluid2(const DRT::ELEMENTS::Fluid2& old) :
DRT::Element(old           ),
gaussrule_  (old.gaussrule_),
is_ale_     (old.is_ale_   ),
dismode_    (old.dismode_  ),
data_       (old.data_     ),
saccn_      (old.saccn_    ),
svelnp_     (old.svelnp_   ),
sveln_      (old.sveln_    )
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid2 and return pointer to it (public) |
 |                                                          gammi 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid2::Clone() const
{
  DRT::ELEMENTS::Fluid2* newelement = new DRT::ELEMENTS::Fluid2(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Fluid2::Shape() const
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
 |  Pack data                                                  (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // Gaussrule
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  // is_ale_
  AddtoPack(data,is_ale_);

  // history variables
  AddtoPack(data,saccn_.M());
  AddtoPack(data,saccn_.N());

  int size = saccn_.N()*saccn_.M()*sizeof(double);
  AddtoPack(data,saccn_ .A(),size);
  AddtoPack(data,svelnp_.A(),size);
  AddtoPack(data,sveln_ .A(),size);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  // discretization mode
  AddtoPack(data,dismode_);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2::Unpack(const vector<char>& data)
{
  vector<char>::size_type position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  dsassert(type == UniqueParObjectId(), "wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // Gaussrule
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = GaussRule2D(gausrule_integer); //explicit conversion from integer to enum
  // is_ale_
  ExtractfromPack(position,data,is_ale_);

  // history variables (subscale velocities, accelerations and pressure)
  {
    int firstdim;
    int secondim;

    ExtractfromPack(position,data,firstdim);
    ExtractfromPack(position,data,secondim);

    saccn_ .Shape(firstdim,secondim);
    svelnp_.Shape(firstdim,secondim);
    sveln_ .Shape(firstdim,secondim);

    int size = firstdim*secondim*sizeof(double);

    ExtractfromPack(position,data,&(saccn_ .A()[0]),size);
    ExtractfromPack(position,data,&(svelnp_.A()[0]),size);
    ExtractfromPack(position,data,&(sveln_ .A()[0]),size);
  }

  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  ExtractfromPack(position,data,dismode_);

  if (position != data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 11/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::~Fluid2()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 11/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2::Print(ostream& os) const
{
  os << "Fluid2 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> >  DRT::ELEMENTS::Fluid2::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Fluid2Line,Fluid2>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of Surfaces (length 1) (public)               gammi 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> >  DRT::ELEMENTS::Fluid2::Surfaces()
{
  vector<RCP<Element> > surfaces(1);
  surfaces[0]= rcp(this, false);
  return surfaces;
}

int DRT::ELEMENTS::Fluid2::NumDofPerNode(const DRT::Node& node) const
{
	// non-equal order elements
	// distinguish type of node
	switch(Shape())
	{
	case tri3:
	case quad4:
	case tri6:
	case quad8:
		return 3;	// standard
		break;
	case quad9:
	{
		if(dismode_ == dismod_equalorder)	return 3;
		else if(dismode_ == dismod_taylorhood)	// Taylor Hood
		{
			if(node.Id() == (NodeIds())[4] ||
					node.Id() == (NodeIds())[5] ||
					node.Id() == (NodeIds())[6] ||
					node.Id() == (NodeIds())[7] ||
					node.Id() == (NodeIds())[8])
				 return 2; // no "corner"-node, but edge-node
			else return 3;
		}
		break;
	}
	case nurbs4:	// 4 control point first order nurbs surface element
	case nurbs9:	// 9 control point second order nurbs surface element
		return 3;
		break;
	default:
		dserror("Shape of Element not supported. (use tri or quad elements!)");
	}

	return 3;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
