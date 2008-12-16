/*!----------------------------------------------------------------------
\file drt_condition.cpp

<pre>
-------------------------------------------------------------------------
                 BACI finite element library subsystem
            Copyright (2008) Technical University of Munich
              
Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed, 
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library may solemnly used in conjunction with the BACI contact library
for purposes described in the above mentioned contract.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de) 
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de                   

-------------------------------------------------------------------------
</pre>

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
  else if (Type()==PointLocsys)                 os << "Point local coordinate system condition: ";
  else if (Type()==LineLocsys)                  os << "Line local coordinate system condition: ";
  else if (Type()==SurfaceLocsys)               os << "Surface local coordinate system condition: ";
  else if (Type()==VolumeLocsys)                os << "Volume local coordinate system condition: ";
  else if (Type()==FSICoupling)                 os << "FSI Coupling condition:";
  else if (Type()==XFEMCoupling)                os << "XFEM Coupling condition:";
  else if (Type()==LineLIFTDRAG)                os << "Line LIFTDRAG condition:";
  else if (Type()==SurfLIFTDRAG)                os << "Surf LIFTDRAG condition:";
  else if (Type()==SurfaceTension)              os << "Surface tension condition:";
  else if (Type()==Surfactant)                  os << "Surfactant condition:";
  else if (Type()==MicroBoundary)               os << "Microscale boundary condition:";
  else if (Type()==VolumeConstraint_3D)         os << "Volume constraint surface boundary condition:";
  else if (Type()==AreaConstraint_3D)           os << "Area constraint surface boundary condition:";
  else if (Type()==AreaConstraint_2D)           os << "Area constraint surface boundary condition:";
  else if (Type()==VolumeMonitor_3D)            os << "Volume monitor condition";
  else if (Type()==AreaMonitor_3D)              os << "Area monitor condition";
  else if (Type()==AreaMonitor_2D)              os << "Area monitor condition";
  else if (Type()==ImpedanceCond)               os << "Impedance boundary condition";
  else if (Type()==MPC_NodeOnPlane_3D)          os << "Multipoint constraint on a plane";
  else if (Type()==MPC_NodeOnLine_2D)           os << "Multipoint constraint on a line";
  else if (Type()==LJ_Potential)                os << "Lennard-Jones potential on a surface";
  else if (Type()==LineWeakDirichlet)           os << "line weak Dirichlet condition";
  else if (Type()==LinePeriodic)                os << "line periodic boundary condition";
  else if (Type()==SurfacePeriodic)             os << "surface periodic boundary condition";
  else if (Type()==Brownian_Motion)             os << "stochastical surface condition (Brownian Motion)";
  else if (Type()==FilamentNumber)              os << "line condition for polymer networks";
  else if (Type()==ForceSensor)                 os << "marking points in a system where force sensors are applied";
  else if (Type()==FlowRateThroughSurface_3D)   os << "Monitor flow rate through a interface";
  else if (Type()==ImpulsRateThroughSurface_3D) os << "Monitor impuls rate through a interface";
  
  
  else dserror("no output string for condition defined in DRT::Condition::Print");
  Container::Print(os);
  if ((int)geometry_.size())
  {
    os << endl;
    os << "Elements of this condition:\n";
    map<int,RefCountPtr<DRT::Element> >::const_iterator curr;
    for (curr=geometry_.begin(); curr!=geometry_.end(); ++curr)
      os << "      " << *(curr->second) << endl;
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


/*----------------------------------------------------------------------*
 |                                                             (public) |
 |  Adjust Ids of elements in order to obtain unique Ids within one     |
 |  condition type                                                      |
 |                                                             lw 12/07 |
 *----------------------------------------------------------------------*/
void DRT::Condition::AdjustId(const int shift)
{
  map<int,RefCountPtr<DRT::Element> > geometry;
  map<int,RefCountPtr<DRT::Element> >::iterator iter;

  for (iter=geometry_.begin();iter!=geometry_.end();++iter)
  {
    iter->second->SetId(iter->first+shift);
    geometry[iter->first+shift]=geometry_[iter->first];
  }

  swap(geometry_, geometry);

  return;
}


#endif  // #ifdef CCADISCRET
