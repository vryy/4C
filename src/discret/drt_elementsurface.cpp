/*!----------------------------------------------------------------------
\file drt_elementsurface.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_elementsurface.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "linalg_utils.H"
#include "drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ElementSurface::ElementSurface(int id, int owner,
                                    int nnode, const int* nodeids,
                                    DRT::Node** nodes) :
DRT::Element(id,element_elementsurface,owner)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 01/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ElementSurface::ElementSurface(const DRT::ElementSurface& old) :
DRT::Element(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ElementSurface::Clone() const
{
  DRT::ElementSurface* newelement = new DRT::ElementSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
const char* DRT::ElementSurface::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // pack base class
  int         basesize = 0;
  const char* basedata = Element::Pack(basesize);
  
  // calculate size of vector data
  size = 
  sizeint +                    // size itself
  sizeint +                    // type of this instance of ParObject, see top of ParObject.H
  basesize +                   // base class data
  0;            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = UniqueParObjectId();
  AddtoPack(position,data,type);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // continue to add stuff here  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 01/07 |
 *----------------------------------------------------------------------*/
bool DRT::ElementSurface::Unpack(const char* data)
{
  //const int sizeint    = sizeof(int);
  //const int sizeforcetype = sizeof(enum ForceType);
  //const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  int position = 0;

  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("Wrong instance type in data");
  // extract base class
  int basesize = SizePack(&data[position]);
  Element::Unpack(&data[position]);
  position += basesize;
  // extract more stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::ElementSurface::~ElementSurface()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ElementSurface::Print(ostream& os) const
{
  os << "ElementSurface ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return ElementSurfaceRegister (public)      mwgee 01/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ElementSurface::ElementRegister() const
{
  return null;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  mwgee 01/07|
 *----------------------------------------------------------------------*/
int DRT::ElementSurface::EvaluateNeumann(ParameterList& params, 
                                         DRT::Discretization&      discretization,
                                         DRT::Condition&           condition,
                                         vector<int>&              lm,
                                         Epetra_SerialDenseVector& elevec1)
{
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::Utils::ExtractMyValues(*disp,mydisp,lm);
  
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // no. of nodes on this surface
  const int nnode = NumNode();
  if (nnode==3 || nnode==6) dserror("triangles not yet impl.");

  // decide how many gaussian points to use
  // here, 8 or 9 noded surface uses 3x3 integration
  // 4 noded surface uses 2x2 integration
  // 3 or 6 noded surface uses 1 or 3 gaussian points
  double xgpr[3];
  double wgtr[3];
  double xgps[3];
  double wgts[3];
  int    ngp;
  GaussianPoints(nnode,ngp,xgpr,wgtr,xgps,wgts);
  
  // find out how many degrees of freedom we have per node
  int numdf = 1000;
  for (int i=0; i<NumNode(); ++i)
    numdf = min(numdf,Nodes()[i]->Dof().NumDof());
  
    
  
  
  return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::ElementSurface::GaussianPoints(const int nnode,
                                         int&   ngp,
                                         double xgpr[],
                                         double wgtr[],
                                         double xgps[],
                                         double wgts[])
{
  const double invsqrtthree = 1./sqrt(3.);
  const double sqrtthreeinvfive = sqrt(3./5.);
  const double wgt  = 5.0/9.0;
  const double wgt0 = 8.0/9.0;
  
  switch (nnode)
  {
  case 4:
    ngp = 2;
    xgpr[0] = -invsqrtthree;
    xgpr[1] =  invsqrtthree;
    xgpr[2] =  0.0;
    wgtr[0] =  1.0;
    wgtr[1] =  1.0;
    wgtr[2] =  0.0;
    xgps[0] = -invsqrtthree;
    xgps[1] =  invsqrtthree;
    xgps[2] =  0.0;
    wgts[0] =  1.0;
    wgts[1] =  1.0;
    wgts[2] =  0.0;
  break;
  case 8:
  case 9:
    ngp = 3;
    xgpr[0] = -sqrtthreeinvfive;
    xgpr[1] =  0.0;
    xgpr[2] =  sqrtthreeinvfive;
    wgtr[0] =  wgt;
    wgtr[1] =  wgt0;
    wgtr[2] =  wgt;
    xgps[0] = -sqrtthreeinvfive;
    xgps[1] =  0.0;
    xgps[2] =  sqrtthreeinvfive;
    wgts[0] =  wgt;
    wgts[1] =  wgt0;
    wgts[2] =  wgt;
  break;
  case 3:
  {
    ngp = 1;
    const double third = 1.0/3.0;
    xgpr[0] =  third;
    xgpr[1] =  0.0;
    xgpr[2] =  0.0;
    xgps[0] =  third;
    xgps[1] =  0.0;
    xgps[2] =  0.0;
    wgtr[0] =  0.5;
    wgtr[1] =  0.0;
    wgtr[2] =  0.0;
    wgts[0] =  0.5;
    wgts[1] =  0.0;
    wgts[2] =  0.0;
  }
  break;
  case 6:
  {
    ngp = 3;
    const double wgt = 1.0/6.0;
    xgpr[0] =  0.5;
    xgpr[1] =  0.5;
    xgpr[2] =  0.0;
    xgps[0] =  0.0;
    xgps[1] =  0.5;
    xgps[2] =  0.5;
    wgtr[0] =  wgt;
    wgtr[1] =  wgt;
    wgtr[2] =  wgt;
    wgts[0] =  wgt;
    wgts[1] =  wgt;
    wgts[2] =  wgt;
  }
  break;
  default:
    dserror("Unknown number of nodes");
  break;
  }
  
    
  return;
}


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
