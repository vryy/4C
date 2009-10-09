/*!----------------------------------------------------------------------
\file fluid3.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*/
// map to convert strings to actions (stabilization)
/*----------------------------------------------------------------------*/
map<string,DRT::ELEMENTS::Fluid3::StabilisationAction> DRT::ELEMENTS::Fluid3::stabstrtoact_;

/*----------------------------------------------------------------------*
 |  ctor (public)                                            gammi 02/08|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::Fluid3(int id, int owner) :
DRT::Element(id,element_fluid3,owner),
is_ale_(false),
data_()
{
    gaussrule_ = intrule3D_undefined;

    Cs_delta_sq_=0;

    saccn_ .Shape(0,0);
    svelnp_.Shape(0,0);
    sveln_ .Shape(0,0);

    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::Fluid3(const DRT::ELEMENTS::Fluid3& old) :
DRT::Element(old             ),
gaussrule_  (old.gaussrule_  ),
is_ale_     (old.is_ale_     ),
data_       (old.data_       ),
Cs_delta_sq_(old.Cs_delta_sq_),
saccn_      (old.saccn_      ),
svelnp_     (old.svelnp_     ),
sveln_      (old.sveln_      )
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Fluid3 and return pointer to it (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Fluid3::Clone() const
{
  DRT::ELEMENTS::Fluid3* newelement = new DRT::ELEMENTS::Fluid3(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Fluid3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  5: return pyramid5;
  case  6: return wedge6;
  case  8: return hex8;
  case 10: return tet10;
  case 15: return wedge15;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Pack(vector<char>& data) const
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
  // Cs_delta_sq_, the Smagorinsky constant for the dynamic Smagorinsky model
  AddtoPack(data,Cs_delta_sq_);

  // history variables
  AddtoPack(data,saccn_.M());
  AddtoPack(data,saccn_.N());

  int size = saccn_.M()*saccn_.N()*sizeof(double);

  AddtoPack(data,saccn_ .A(),size);
  AddtoPack(data,svelnp_.A(),size);
  AddtoPack(data,sveln_ .A(),size);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                          gammi 02/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Unpack(const vector<char>& data)
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
  gaussrule_ = GaussRule3D(gausrule_integer); //explicit conversion from integer to enum
  // is_ale_
  ExtractfromPack(position,data,is_ale_);
  // extract Cs_delta_sq_, the Smagorinsky constant for the dynamic
  // Smagorinsky model
  ExtractfromPack(position,data,Cs_delta_sq_);

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

  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            gammi 02/08|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3::~Fluid3()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              gammi 02/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3::Print(ostream& os) const
{
  os << "Fluid3 ";
  Element::Print(os);
  //cout << endl;
  cout << data_;
  return;
}


/*----------------------------------------------------------------------*
 |  allocate and return Fluid3Register (public)              mwgee 02/08|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Fluid3::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Fluid3Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines              (public)                  gjb 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Fluid3Line,Fluid3>(DRT::UTILS::buildLines,this);
}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                            gjb 05/08|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<Fluid3Surface,Fluid3>(DRT::UTILS::buildSurfaces,this);
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Fluid3::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}

/*----------------------------------------------------------------------*
 |  calculate volume of current element (public)                tw 10/09|
 *----------------------------------------------------------------------*/
double DRT::ELEMENTS::Fluid3::CalculateVolume(ParameterList& params, DRT::Discretization& discretization,vector<int>& lm) const
{
  // the number of nodes
  const int numnode = NumNode();

  const DRT::Node*const* nodes = Nodes();
  if(!nodes)  dserror("no nodes? make sure that Discretization::FillComplete is called!");

  // --------------------------------------------------
  // create matrix objects for nodal values
  Epetra_SerialDenseMatrix  edispnp(3,numnode);


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
      int fi    =4*i;
      int fip   =fi+1;
      int fipp  =fip+1;
      edispnp(0,i)    = mydispnp   [fi  ];
      edispnp(1,i)    = mydispnp   [fip ];
      edispnp(2,i)    = mydispnp   [fipp];
    }
  }

  Epetra_SerialDenseMatrix xyze(3,numnode);

  // get node coordinates
  for (int inode=0; inode<numnode; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze(0,inode) = x[0];
    xyze(1,inode) = x[1];
    xyze(2,inode) = x[2];
  }

  // add displacement, when fluid nodes move in the ALE case
  if (is_ale_)
  {
    for (int inode=0; inode<numnode; inode++)
    {
      xyze(0,inode) += edispnp(0,inode);
      xyze(1,inode) += edispnp(1,inode);
      xyze(2,inode) += edispnp(2,inode);
    }
  }

  switch(Shape())
  {
  case DRT::Element::tet4:
  case DRT::Element::tet10:
  {
    /*
     *        3 +
     *         /|\.
     *          |   .
     *        c |     .
     *          |       .
     *        0 +-------->+
     *      a  /    b     2
     *        /
     *      |/_ 1
     */
    // TODO: not tested
    double a1 = 0.0; double a2 = 0.0; double a3 = 0.0;
    double b1 = 0.0; double b2 = 0.0; double b3 = 0.0;
    double c1 = 0.0; double c2 = 0.0; double c3 = 0.0;
    a1 = xyze(0,1) - xyze(0,0); a2 = xyze(1,1) - xyze(1,0); a3 = xyze(2,1) - xyze(2,0);
    b1 = xyze(0,2) - xyze(0,0); b2 = xyze(1,2) - xyze(1,0); b3 = xyze(2,2) - xyze(2,0);
    c1 = xyze(0,3) - xyze(0,0); c2 = xyze(1,3) - xyze(1,0); c3 = xyze(2,3) - xyze(2,0);

    // V = 1/6 * |(a x b) * c|
    double ret = (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    return ret; // volume of tet
  }
  break;
  case DRT::Element::wedge6:
  case DRT::Element::wedge15:
  {
    // TODO: not tested
    double a1 = 0.0; double a2 = 0.0; double a3 = 0.0;
    double b1 = 0.0; double b2 = 0.0; double b3 = 0.0;
    double c1 = 0.0; double c2 = 0.0; double c3 = 0.0;
    // ~ tet1
    a1 = xyze(0,1) - xyze(0,0); a2 = xyze(1,1) - xyze(1,0); a3 = xyze(2,1) - xyze(2,0);
    b1 = xyze(0,2) - xyze(0,0); b2 = xyze(1,2) - xyze(1,0); b3 = xyze(2,2) - xyze(2,0);
    c1 = xyze(0,3) - xyze(0,0); c2 = xyze(1,3) - xyze(1,0); c3 = xyze(2,3) - xyze(2,0);
    double ret = (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    // ~ tet2
    a1 = xyze(0,3) - xyze(0,4); a2 = xyze(1,3) - xyze(1,4); a3 = xyze(2,3) - xyze(2,4);
    b1 = xyze(0,5) - xyze(0,4); b2 = xyze(1,5) - xyze(1,4); b3 = xyze(2,5) - xyze(2,4);
    c1 = xyze(0,1) - xyze(0,4); c2 = xyze(1,1) - xyze(1,4); c3 = xyze(2,1) - xyze(2,4);
    ret += (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    // ~ tet3
    a1 = xyze(0,3) - xyze(0,5); a2 = xyze(1,3) - xyze(1,5); a3 = xyze(2,3) - xyze(2,5);
    b1 = xyze(0,2) - xyze(0,5); b2 = xyze(1,2) - xyze(1,5); b3 = xyze(2,2) - xyze(2,5);
    c1 = xyze(0,1) - xyze(0,5); c2 = xyze(1,1) - xyze(1,5); c3 = xyze(2,1) - xyze(2,5);
    ret += (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    return ret;
  }
  break;
  case DRT::Element::hex8:
  case DRT::Element::hex20:
  case DRT::Element::hex27:
  {
    // calculate volume of the 6 tetrahedal elements that "form" the hexahedral element
    double a1 = 0.0; double a2 = 0.0; double a3 = 0.0;
    double b1 = 0.0; double b2 = 0.0; double b3 = 0.0;
    double c1 = 0.0; double c2 = 0.0; double c3 = 0.0;

    // ~ wedge 1
    a1 = xyze(0,1) - xyze(0,0); a2 = xyze(1,1) - xyze(1,0); a3 = xyze(2,1) - xyze(2,0);
    b1 = xyze(0,3) - xyze(0,0); b2 = xyze(1,3) - xyze(1,0); b3 = xyze(2,3) - xyze(2,0);
    c1 = xyze(0,4) - xyze(0,0); c2 = xyze(1,4) - xyze(1,0); c3 = xyze(2,4) - xyze(2,0);

    // V = 1/6 * |(a x b) * c|
    double ret = (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    a1 = xyze(0,4) - xyze(0,7); a2 = xyze(1,4) - xyze(1,7); a3 = xyze(2,4) - xyze(2,7);
    b1 = xyze(0,3) - xyze(0,7); b2 = xyze(1,3) - xyze(1,7); b3 = xyze(2,3) - xyze(2,7);
    c1 = xyze(0,5) - xyze(0,7); c2 = xyze(1,5) - xyze(1,7); c3 = xyze(2,5) - xyze(2,7);

    ret += (double)  1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    a1 = xyze(0,5) - xyze(0,1); a2 = xyze(1,5) - xyze(1,1); a3 = xyze(2,5) - xyze(2,1);
    b1 = xyze(0,3) - xyze(0,1); b2 = xyze(1,3) - xyze(1,1); b3 = xyze(2,3) - xyze(2,1);
    c1 = xyze(0,4) - xyze(0,1); c2 = xyze(1,4) - xyze(1,1); c3 = xyze(2,4) - xyze(2,1);

    ret += (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    // ~ wedge 2
    a1 = xyze(0,2) - xyze(0,1); a2 = xyze(1,2) - xyze(1,1); a3 = xyze(2,2) - xyze(2,1);
    b1 = xyze(0,3) - xyze(0,1); b2 = xyze(1,3) - xyze(1,1); b3 = xyze(2,3) - xyze(2,1);
    c1 = xyze(0,5) - xyze(0,1); c2 = xyze(1,5) - xyze(1,1); c3 = xyze(2,5) - xyze(2,1);
    ret += (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    a1 = xyze(0,5) - xyze(0,3); a2 = xyze(1,5) - xyze(1,3); a3 = xyze(2,5) - xyze(2,3);
    b1 = xyze(0,2) - xyze(0,3); b2 = xyze(1,2) - xyze(1,3); b3 = xyze(2,2) - xyze(2,3);
    c1 = xyze(0,7) - xyze(0,3); c2 = xyze(1,7) - xyze(1,3); c3 = xyze(2,7) - xyze(2,3);
    ret += (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    a1 = xyze(0,5) - xyze(0,6); a2 = xyze(1,5) - xyze(1,6); a3 = xyze(2,5) - xyze(2,6);
    b1 = xyze(0,7) - xyze(0,6); b2 = xyze(1,7) - xyze(1,6); b3 = xyze(2,7) - xyze(2,6);
    c1 = xyze(0,2) - xyze(0,6); c2 = xyze(1,2) - xyze(1,6); c3 = xyze(2,2) - xyze(2,6);
    ret += (double) 1/6 * abs((a2*b3-a3*b2)*c1 + (a3*b1-a1*b3)*c2 + (a1*b2-a2*b1)*c3);

    return ret;
  }
  break;
  default:
    dserror("CalculateVolume for current element type not implemented.");
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
DRT::ELEMENTS::Fluid3Register::Fluid3Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Register::Fluid3Register(
                               const DRT::ELEMENTS::Fluid3Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Fluid3Register* DRT::ELEMENTS::Fluid3Register::Clone() const
{
  return new DRT::ELEMENTS::Fluid3Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Register::Pack(vector<char>& data) const
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
void DRT::ELEMENTS::Fluid3Register::Unpack(const vector<char>& data)
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
DRT::ELEMENTS::Fluid3Register::~Fluid3Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Fluid3Register::Print(ostream& os) const
{
  os << "Fluid3Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
