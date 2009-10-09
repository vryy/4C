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

using namespace DRT::UTILS;

map<string,DRT::ELEMENTS::Fluid2::StabilisationAction> DRT::ELEMENTS::Fluid2::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2::Fluid2(int id, int owner) :
DRT::Element(id,element_fluid2,owner),
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
  int position = 0;
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

  if (position != (int)data.size())
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
 |  allocate and return Fluid2Register (public)              gammi 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Fluid2::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Fluid2Register(Type()));
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

/*----------------------------------------------------------------------*
 |  calculate volume (=surface) of current element (public)     tw 10/09|
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::Fluid2::CalculateSurface(ParameterList& params, DRT::Discretization& discretization, vector<int>& lm) const
{
  // the number of nodes
  const int numnode = NumNode();

  const DRT::Node*const* nodes = Nodes();
  if(!nodes)  dserror("no nodes? make sure that Discretization::FillComplete is called!");

  // --------------------------------------------------
  // create matrix objects for nodal values
  Epetra_SerialDenseMatrix  edispnp(2,numnode);


  if(is_ale_)
  {
    // get most recent displacements
    RCP<const Epetra_Vector> dispnp = discretization.GetState("dispnp");

    if (dispnp==null)
    {
      dserror("Cannot get state vector 'dispnp'");
    }

    vector<double> mydispnp(lm.size());
    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);

    // extract velocity part from "mygridvelaf" and get
    // set element displacements
    for (int i=0;i<numnode;++i)
    {
      int fi    =3*i;
      int fip   =fi+1;
      edispnp(0,i)    = mydispnp   [fi  ];
      edispnp(1,i)    = mydispnp   [fip ];
    }
  }

  Epetra_SerialDenseMatrix xyze(2,numnode);

  // get node coordinates
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
  }

  // add displacement, when fluid nodes move in the ALE case
  if (is_ale_)
  {
    for (int inode=0; inode<numnode; inode++)
    {
      xyze(0,inode) += edispnp(0,inode);
      xyze(1,inode) += edispnp(1,inode);
    }
  }

  switch(Shape())
  {
  case DRT::Element::tri3:
  case DRT::Element::tri6:
  {
    /*  2
     *   +
     *  /|\.
     *   |   .
     * b |     .
     *   |       .
     *   +-------->+
     *  0    a      1
     */

    double ax = 0.0; double ay = 0.0;
    double bx = 0.0; double by = 0.0;
    ax = xyze(0,1) - xyze(0,0); ay = xyze(1,1) - xyze(1,0);
    bx = xyze(0,2) - xyze(0,0); by = xyze(1,2) - xyze(1,0);

    // calculate determinant
    double det = ax*by - ay*bx;
    if(det < 0.0) dserror("error in FLUID2::CalculateSurface. Surface < 0 ???");

    return det/2.0; // area of tri element
  }
  break;
  case DRT::Element::quad4:
  case DRT::Element::quad8:
  case DRT::Element::quad9:
  {
    /*  3     c     2
     *   +<--------+
     *  /|\.       |
     *   |   .     |
     * b |     .   | d
     *   |       .\|/
     *   +-------->+
     *  0    a      1
     */

    double ax = 0.0; double ay = 0.0;
    double bx = 0.0; double by = 0.0;
    double cx = 0.0; double cy = 0.0;
    double dx = 0.0; double dy = 0.0;
    ax = xyze(0,1) - xyze(0,0); ay = xyze(1,1) - xyze(1,0);
    bx = xyze(0,3) - xyze(0,0); by = xyze(1,3) - xyze(1,0);
    cx = xyze(0,3) - xyze(0,2); cy = xyze(1,3) - xyze(1,2);
    dx = xyze(0,1) - xyze(0,2); dy = xyze(1,1) - xyze(1,2);

    // calculate vol triangle (a x b)
    double tri1 = 0.5 * (ax*by - ay*bx);

    // calculate vol triangle (c x d)
    double tri2 = 0.5 * (cx*dy - cy*dx);

    double ret = tri1 + tri2;
    if(ret < 0.0) dserror("error in FLUID2::CalculateSurface. Surface < 0 ???");

    return ret;
  }
  break;
  default:
    dserror("CalculateSurface for current element type not implemented. (only tri and quad supported)");
  }

  return -1.0;
}

//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Register::Fluid2Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Register::Fluid2Register(
                               const DRT::ELEMENTS::Fluid2Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Register* DRT::ELEMENTS::Fluid2Register::Clone() const
{
  return new DRT::ELEMENTS::Fluid2Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Register::Pack(vector<char>& data) const
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
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid2Register::~Fluid2Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid2Register::Print(ostream& os) const
{
  os << "Fluid2Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
