/*!----------------------------------------------------------------------
\file wall1.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 01/08/|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::Wall1(int id, int owner) :
DRT::Element(id,element_wall1,owner),
data_(),
material_(0),
thickness_(0.0),
gaussrule_(DRT::UTILS::intrule2D_undefined),
wtype_(plane_none),
stresstype_(w1_none),
iseas_(false),
eastype_(eas_vague)

{
//  tsi_couptyp_ = tsi_coup_none;
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 01/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::Wall1(const DRT::ELEMENTS::Wall1& old) :
DRT::Element(old),
data_(old.data_),
material_(old.material_),
thickness_(old.thickness_),
gaussrule_(old.gaussrule_),
wtype_(old.wtype_),
stresstype_(old.stresstype_),
iseas_(old.iseas_),
eastype_(old.eas_vague)
// tsi_couptyp_(old.tsi_couptyp_)

{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Wall1 and return pointer to it (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Wall1::Clone() const
{
  DRT::ELEMENTS::Wall1* newelement = new DRT::ELEMENTS::Wall1(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          mgit 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Wall1::Shape() const
{
  switch (NumNode())
  {
  case 4: return quad4;
  case 8: return quad8;
  case 9: return quad9;
  case 3: return tri3;
  case 6: return tri6;

  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // material_
  AddtoPack(data,material_);
  //thickness
  AddtoPack(data,thickness_);
  // plane strain or plane stress information
  AddtoPack(data,wtype_);
  // gaussrule_
  AddtoPack(data,gaussrule_); //implicit conversion from enum to integer
  // stresstype
  AddtoPack(data,stresstype_);
  // eas
  AddtoPack(data,iseas_);
  // eas type
  AddtoPack(data,eastype_);
//  //tsi
//  AddtoPack(data,tsi_couptyp_);
  //data
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Element
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Element::Unpack(basedata);
  // material_
  ExtractfromPack(position,data,material_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // plane strain or plane stress information_
  ExtractfromPack(position,data,wtype_);
  // gaussrule_
  int gausrule_integer;
  ExtractfromPack(position,data,gausrule_integer);
  gaussrule_ = DRT::UTILS::GaussRule2D(gausrule_integer); //explicit conversion from integer to enum
  // stresstype_
  ExtractfromPack(position,data,stresstype_);
  // iseas_
  ExtractfromPack(position,data,iseas_);
  // eastype_
  ExtractfromPack(position,data,eastype_);
//  // tsi_couptype
//  ExtractfromPack(position,data,tsi_couptyp_);
  //data
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1::~Wall1()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Print(ostream& os) const
{
  os << "Wall1 ";
  Element::Print(os);
  os << " gaussrule_: " << gaussrule_ << " ";
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Wall1Register (public)              mgit 03/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Wall1::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Wall1Register(Type()));
}



/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mgit 07/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> >  DRT::ELEMENTS::Wall1::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Wall1Line,Wall1>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mgit 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> >  DRT::ELEMENTS::Wall1::Surfaces()
{
  vector<RCP<Element> > surfaces(1);
  surfaces[0]= rcp(this, false);
  return surfaces;
}


/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes      popp 08/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_expol
(
    Epetra_SerialDenseMatrix& stresses,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2
)
{
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
  const DiscretizationType dt = this->Shape();
  const int numgp = intpoints.nquad;
  const int numnode = NumNode();
  Epetra_SerialDenseVector funct(numnode);

  static Epetra_SerialDenseMatrix expol(numnode,numgp);
  static bool isfilled;

  if (isfilled==false)
  {
  	// tri3, quad4, tri6, quad8 and quad9
    if (dt==tri3 or dt==tri6 or dt==quad4 or dt==quad8 or dt==quad9)
    {
    	// loop over gaussian points
      for (int ip=0; ip<intpoints.nquad; ++ip)
      {
        // gaussian coordinates
        const double e1 = intpoints.qxg[ip][0];
        const double e2 = intpoints.qxg[ip][1];

        // coordinates of the extrapolated points
        double e1expol;
        double e2expol;

        if (e1!=0)
        {
        	e1expol = 1/e1;
        }	
        else
        {
        	e1expol = 0;	
        }

        if (e2!=0)
        {
        e2expol = 1/e2;
        }	
        else
        {
        e2expol = 0;	
        }

        // shape functions for the extrapolated coordinates
        DRT::UTILS::shape_function_2D(funct,e1expol,e2expol,dt);

        // extrapolation matrix
        for(int i=0;i<numnode;++i)
        {
        	expol(ip,i)=funct(i);
        }
      }
      isfilled = true;
    }

    else dserror("extrapolation not yet implemented for this element type");
  }

  Epetra_SerialDenseMatrix nodalstresses(numnode,Wall1::numstr_);
  nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);

  // distribute nodal stresses to elevectors for assembling
  for (int i=0;i<numnode;++i)
  {
    elevec1(2*i)=nodalstresses(i,0);
    elevec1(2*i+1)=nodalstresses(i,1);
  }
  for (int i=0;i<numnode;++i)
  {
    elevec2(2*i)=nodalstresses(i,2);
    elevec2(2*i+1)=nodalstresses(i,3);
  }

}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1Register::Wall1Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1Register::Wall1Register(
                               const DRT::ELEMENTS::Wall1Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1Register* DRT::ELEMENTS::Wall1Register::Clone() const
{
  return new DRT::ELEMENTS::Wall1Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            mgit 03/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Register::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // base class ElementRegister
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  ElementRegister::Unpack(basedata);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mgit 03/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Wall1Register::~Wall1Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Register::Print(ostream& os) const
{
  os << "Wall1Register ";
  ElementRegister::Print(os);
  return;
}


int DRT::ELEMENTS::Wall1Register::Initialize(DRT::Discretization& dis)
{
  return 0;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_WALL1
