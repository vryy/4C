/*!----------------------------------------------------------------------
\file shell8.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "shell8.H"
#include "drt_discret.H"
#include "drt_utils.H"
#include "drt_dserror.H"




/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::Shell8(int id, int owner) :
DRT::Element(id,element_shell8,owner),
forcetype_(s8_none),
thickness_(0.0),
ngptri_(0),
nhyb_(0),
ans_(0),
sdc_(1.0),
material_(0),
data_()
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  eas_[0] = eas_[1] = eas_[2] = eas_[3] = eas_[4] = 0;
  lines_.resize(0);
  lineptrs_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::Shell8(const DRT::Elements::Shell8& old) :
DRT::Element(old),
forcetype_(old.forcetype_),
thickness_(old.thickness_),
ngptri_(old.ngptri_),
nhyb_(old.nhyb_),
ans_(old.ans_),
sdc_(old.sdc_),
material_(old.material_),
data_(old.data_),
lines_(old.lines_),
lineptrs_(old.lineptrs_)
{
  for (int i=0; i<3; ++i) ngp_[i] = old.ngp_[i];
  for (int i=0; i<5; ++i) eas_[i] = old.eas_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Shell8 and return pointer to it (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::Shell8::Clone() const
{
  DRT::Elements::Shell8* newelement = new DRT::Elements::Shell8(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Elements::Shell8::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizeforcetype = sizeof(enum ForceType);
  const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // pack base class
  int         basesize = 0;
  const char* basedata = Element::Pack(basesize);
  
  // pack data_
  int         consize = 0;
  const char* condata = data_.Pack(consize);
   
  // calculate size of vector data
  size = 
  sizeint +                    // size itself
  sizeint +                    // type of this instance of ParObject, see top of ParObject.H
  basesize +                   // base class data
  sizeforcetype +              // forcetype_
  sizedouble +                 // thickness_
  3*sizeint +                  // ngp_
  sizeint +                    // ngptri_  
  sizeint +                    // nhyb_
  5*sizeint +                  // eas_
  sizeint +                    // ans_
  sizedouble +                 // sdc_
  sizeint +                    // material_
  consize +                    // data_
  0;            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Shell8;
  AddtoPack(position,data,type);
  // add base class
  AddtoPack(position,data,basedata,basesize);
  delete [] basedata;
  // forcetype_
  AddtoPack(position,data,&forcetype_,sizeof(enum ForceType));
  // thickness_
  AddtoPack(position,data,thickness_);
  // ngp_
  AddtoPack(position,data,ngp_,3*sizeint);
  // ngptri_
  AddtoPack(position,data,ngptri_);
  // nhyb_
  AddtoPack(position,data,nhyb_);
  // eas_
  AddtoPack(position,data,eas_,5*sizeint);
  // ans_
  AddtoPack(position,data,ans_);
  // sdc_
  AddtoPack(position,data,sdc_);
  // material_
  AddtoPack(position,data,material_);
  // data_
  AddtoPack(position,data,condata,consize);
  delete [] condata;
  // continue to add stuff here  
    
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}

/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Elements::Shell8::Unpack(const char* data)
{
  const int sizeint    = sizeof(int);
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
  if (type != ParObject_Shell8) dserror("Wrong instance type in data");
  // extract base class
  int basesize = SizePack(&data[position]);
  Element::Unpack(&data[position]);
  position += basesize;
  // forcetype_
  ExtractfromPack(position,data,forcetype_);
  // thickness_
  ExtractfromPack(position,data,thickness_);
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeint);
  // ngptri_
  ExtractfromPack(position,data,ngptri_);
  // nhyb_
  ExtractfromPack(position,data,nhyb_);
  // eas_
  ExtractfromPack(position,data,eas_,5*sizeint);
  // ans_
  ExtractfromPack(position,data,ans_);
  // sdc_
  ExtractfromPack(position,data,sdc_);
  // material_
  ExtractfromPack(position,data,material_);
  // data_
  int consize = SizePack(&data[position]);
  data_.Unpack(&data[position]);
  position += consize;
  
  // extract more stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8::~Shell8()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8::Print(ostream& os) const
{
  os << "Shell8 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return Shell8Register (public)              mwgee 12/06|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::Elements::Shell8::ElementRegister() const
{
  return rcp(new DRT::Elements::Shell8Register(Type()));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Shell8::Lines()
{
  const int nline = NumLine();
  const int numnode = NumNode();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  int nodeids[100];
  DRT::Node* nodes[100];
  if (nline==4)
  {
    if (numnode==4)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      lines_[0] = 
        rcp(new DRT::Elements::Shell8Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[1] = 
        rcp(new DRT::Elements::Shell8Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      lines_[2] = 
        rcp(new DRT::Elements::Shell8Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      lines_[3] = 
        rcp(new DRT::Elements::Shell8Line(3,Owner(),2,nodeids,nodes,this,3));
      lineptrs_[3] = lines_[3].get();
    }
    else if (numnode==9)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[4];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[4];
      lines_[0] = 
        rcp(new DRT::Elements::Shell8Line(0,Owner(),3,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[5];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[5];
      lines_[1] = 
        rcp(new DRT::Elements::Shell8Line(1,Owner(),3,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodeids[2] = NodeIds()[6];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      nodes[2] = Nodes()[6];
      lines_[2] = 
        rcp(new DRT::Elements::Shell8Line(2,Owner(),3,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[7];
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[7];
      lines_[3] = 
        rcp(new DRT::Elements::Shell8Line(3,Owner(),3,nodeids,nodes,this,3));
      lineptrs_[3] = lines_[3].get();
    }
    else dserror("Number of nodes not supported");
  }
  else if (nline==3)
  {
    if (numnode==3)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      lines_[0] = 
        rcp(new DRT::Elements::Shell8Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[1] = 
        rcp(new DRT::Elements::Shell8Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      lines_[2] = 
        rcp(new DRT::Elements::Shell8Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();
    }
    else if (numnode==6)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[3];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[3];
      lines_[0] = 
        rcp(new DRT::Elements::Shell8Line(0,Owner(),3,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[4];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[4];
      lines_[1] = 
        rcp(new DRT::Elements::Shell8Line(1,Owner(),3,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[5];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[5];
      lines_[2] = 
        rcp(new DRT::Elements::Shell8Line(2,Owner(),3,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();
    }
    else dserror("Number of nodes not supported");
  }
  else dserror("Number of lines not supported");
  return (DRT::Element**)(&(lineptrs_[0]));
}

/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          mwgee 01/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::Shell8::Surfaces()
{
  surfaces_.resize(1);
  surfaces_[0] = this;
  return &surfaces_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register::Shell8Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register::Shell8Register(
                               const DRT::Elements::Shell8Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register* DRT::Elements::Shell8Register::Clone() const
{
  return new DRT::Elements::Shell8Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Elements::Shell8Register::Pack(int& size) const
{
  const int sizeint    = sizeof(int);

  // pack base class
  int         basesize = 0;
  const char* basedata = ElementRegister::Pack(basesize);
   
  // calculate size of vector data
  size = 
  sizeint +                     // size itself
  sizeint +                     // type of this instance of ParObject, see top of ParObject.H
  basesize +                    // base class data
  0;                            // continue to add data here...
  
  char* data = new char[size];
  
  // pack stuff into vector
  int position = 0;
  
  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Shell8Register;
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
 |  Unpack data                                                (public) |
 |                                                            gee 12/06 |
 *----------------------------------------------------------------------*/
bool DRT::Elements::Shell8Register::Unpack(const char* data)
{
  int position = 0;
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_Shell8Register) dserror("Wrong instance type in data");
  // extract base class
  int basesize = SizePack(&data[position]);
  ElementRegister::Unpack(&data[position]);
  position += basesize;
  // extract more stuff here

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return true;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 12/06|
 *----------------------------------------------------------------------*/
DRT::Elements::Shell8Register::~Shell8Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                           mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Elements::Shell8Register::Print(ostream& os) const
{
  os << "Shell8Register ";
  ElementRegister::Print(os);
  return;
}





#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
