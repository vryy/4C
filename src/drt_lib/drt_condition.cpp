/*!----------------------------------------------------------------------
\file drt_condition.cpp

\brief A condition of any kind

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_condition.H"
#include "drt_element.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(const int id, const ConditionType type,
                          const bool buildgeometry, const GeometryType gtype) :
Container(),
id_(id),
buildgeometry_(buildgeometry),
type_(type),
gtype_(gtype),
comm_(null)
{ 
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition() :
Container(),
id_(-1),
buildgeometry_(false),
type_(none),
gtype_(NoGeom),
comm_(null)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::Condition(const DRT::Condition& old) :
Container(old),
id_(old.id_),
buildgeometry_(old.buildgeometry_),
type_(old.type_),
gtype_(old.gtype_),
comm_(old.comm_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Condition::~Condition()
{
  return;
}


/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Condition& cond)
{
  cond.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Condition::Print(ostream& os) const
{
  os << "Condition " << Id() << " ";
  if      (Type()==PointDirichlet)              os << "Point Dirichlet boundary condition: ";
  else if (Type()==LineDirichlet)               os << "Line Dirichlet boundary condition: ";
  else if (Type()==SurfaceDirichlet)            os << "Surface Dirichlet boundary condition: ";
  else if (Type()==VolumeDirichlet)             os << "Volume Dirichlet boundary condition: ";
  else if (Type()==PointNeumann)                os << "Point Neumann boundary condition: ";
  else if (Type()==LineNeumann)                 os << "Line Neumann boundary condition: ";
  else if (Type()==SurfaceNeumann)              os << "Surface Neumann boundary condition: ";
  else if (Type()==VolumeNeumann)               os << "Volume Neumann boundary condition: ";
  else if (Type()==Contact)                     os << "Contact boundary condition: ";
  else if (Type()==LineIsothermalNoslip)        os << "Isothermal no-slip wall boundary condition";
  else if (Type()==LineSubsonicInflow)          os << "Subsonic inflow boundary condition:";
  else if (Type()==LineSubsonicOutflow)         os << "Subsonic outflow boundary condition:";
  else if (Type()==XFEMCoupling)                os << "XFEM Coupling condition:";
  else if (Type()==LineLIFTDRAG)                os << "Line LIFTDRAG condition:";
  else if (Type()==SurfLIFTDRAG)                os << "Surf LIFTDRAG condition:";
  else if (Type()==VolumeConstraint_3D)            os << "Volume constraint surface boundary condition:";
 
  else dserror("no output string for condition defined in DRT::Condition::Print");
  Container::Print(os);
  if ((int)geometry_.size())
  {
    cout << endl;
    cout << "Elements of this condition:\n";
    map<int,RefCountPtr<DRT::Element> >::const_iterator curr;
    for (curr=geometry_.begin(); curr!=geometry_.end(); ++curr)
      cout << "      " << *(curr->second) << endl;
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Condition::Pack(vector<char>& data) const
{
  data.resize(0);
  
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class container
  vector<char> basedata;
  Container::Pack(basedata);
  AddtoPack(data,basedata);
  // id_
  AddtoPack(data,id_);
  // buildgeometry_
  AddtoPack(data,buildgeometry_);
  // type_
  AddtoPack(data,type_);
  // gtype_
  AddtoPack(data,gtype_);
  
  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gee 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Condition::Unpack(const vector<char>& data)
{
  int position = 0;
  // extract type
  int type = 0;
  ExtractfromPack(position,data,type);
  if (type != UniqueParObjectId()) dserror("wrong instance type data");
  // extract base class Container
  vector<char> basedata(0);
  ExtractfromPack(position,data,basedata);
  Container::Unpack(basedata);
  // id_
  ExtractfromPack(position,data,id_);
  // buildgeometry_
  ExtractfromPack(position,data,buildgeometry_);
  // type_
  ExtractfromPack(position,data,type_);
  // gtype_
  ExtractfromPack(position,data,gtype_);
  
  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
} 




#endif  // #ifdef CCADISCRET
