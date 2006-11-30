/*!----------------------------------------------------------------------
\file node.cpp
\brief A virtual class for a node

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_node.H"
#include "drt_dserror.H"



/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::Node(int id, const double* coords, const int owner) :
ParObject(),
id_(id),
owner_(owner),
dofset_(),
dentitytype_(on_none),
dentityid_(-1)
{
  for (int i=0; i<3; ++i) x_[i] = coords[i];
  element_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::Node(const DRT::Node& old) :
ParObject(old),
id_(old.id_),
owner_(old.owner_),
dofset_(old.dofset_),
element_(old.element_),
dentitytype_(old.dentitytype_),
dentityid_(old.dentityid_),
condition_(old.condition_) // scary: this is a RefCountPtr !
{
  for (int i=0; i<3; ++i) x_[i] = old.x_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 11/06|
 *----------------------------------------------------------------------*/
DRT::Node::~Node()
{
  return;
}


/*----------------------------------------------------------------------*
 |  Deep copy this instance of Node and return pointer to it (public)   |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
DRT::Node* DRT::Node::Clone() const
{
  DRT::Node* newnode = new DRT::Node(*this);
  return newnode;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 11/06|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const DRT::Node& node)
{
  node.Print(os); 
  return os;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
void DRT::Node::Print(ostream& os) const
{
  // Print id and coordinates
  os << "Node " << setw(12) << Id()
     << " Owner " << setw(4) << Owner()
     << " Coords " 
     << setw(12) << X()[0] << " " 
     << setw(12) << X()[1] << " " 
     << setw(12) << X()[2] << " ";
  // print dofs if there are any
  if (Dof().NumDof())
  {
    os << Dof();
  }
  
  // Print design entity if there is any
  if (dentitytype_ != on_none)
  {
    if      (dentitytype_==on_dnode) os << "on DNODE " << dentityid_ << " ";
    else if (dentitytype_==on_dline) os << "on DLINE " << dentityid_ << " ";
    else if (dentitytype_==on_dsurface) os << "on DSURF " << dentityid_ << " ";
    else if (on_dsurface==on_dvolume) os << "on DVOL " << dentityid_ << " ";
    else dserror("Unknown type of design entity");
  }
  // Print conditions if there are any
  int numcond = condition_.size();
  if (numcond)
  {
    os << endl; // end the previous line
    map<string,RefCountPtr<Condition> >::const_iterator curr;
    for (curr=condition_.begin(); curr != condition_.end(); ++curr)
    {
      os << "Condition : " << curr->first << endl;
      os << *(curr->second);
    }
  }
  return;
}

/*----------------------------------------------------------------------*
 |  Pack data from this element into vector of length size     (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const char* DRT::Node::Pack(int& size) const
{
  const int sizeint    = sizeof(int);
  const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  // get data and size of dofset_
  int dofsetsize=0;
  const char* dofsetpack = dofset_.Pack(dofsetsize);

  // get the size of all conditions
  int condsize=0;
  map<string,RefCountPtr<Condition> >::const_iterator curr;
  for (curr=condition_.begin(); curr != condition_.end(); ++curr)
  {
    condsize += SizeString(curr->first);
    int tmp=0;
    const char* condpack = curr->second->Pack(tmp);
    condsize += tmp;
    delete [] condpack;
  }

  size = 
  sizeint                +   // holds size itself
  sizeint                +   // type of this instance of ParObject, see top of ParObject.H
  sizeint                +   // holds id
  sizeint                +   // owner_
  sizedouble*3           +   // holds x_
  dofsetsize             +   // dofset_
  sizeof(OnDesignEntity) +   // dentitytype_
  sizeint                +   // dentityid_
  sizeint                +   // no. objects in condition_
  condsize               +   // condition_
  0;                         // continue to add data here...

  char* data = new char[size];

  // pack stuff into vector
  int position = 0;

  // add size
  AddtoPack(position,data,size);
  // ParObject type
  int type = ParObject_Node;
  AddtoPack(position,data,type);
  // add id_
  int id = Id();
  AddtoPack(position,data,id);
  // add owner_
  int owner = Owner();
  AddtoPack(position,data,owner);
  // add x_
  AddtoPack(position,data,x_,3*sizedouble);
  // dofset_
  AddtoPack(position,data,dofsetpack,dofsetsize);
  delete [] dofsetpack;
  // dentitytype_
  AddtoPack(position,data,dentitytype_);
  // dentityid_
  AddtoPack(position,data,dentityid_);
  // condition_
  int num = condition_.size(); // no. of objects
  AddtoPack(position,data,num);
  for (curr=condition_.begin(); curr != condition_.end(); ++curr)
  {
    AddStringtoPack(position,data,curr->first);
    int tmp=0;
    const char* condpack = curr->second->Pack(tmp);
    AddtoPack(position,data,condpack,tmp);
    delete [] condpack;
  }
  // continue to add stuff here
  
  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);

  return data;
}


/*----------------------------------------------------------------------*
 |  Unpack data into this element                              (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
bool DRT::Node::Unpack(const char* data)
{
  //const int sizeint    = sizeof(int);
  const int sizedouble = sizeof(double);
  //const int sizechar   = sizeof(char);

  int position = 0;
  
  // extract size
  int size = 0;
  ExtractfromPack(position,data,size);
  // ParObject instance type
  int type=0;
  ExtractfromPack(position,data,type);
  if (type != ParObject_Node) dserror("Wrong instance type in data");
  // extract id_
  ExtractfromPack(position,data,id_);
  // extract id_
  ExtractfromPack(position,data,owner_);
  // extract x_
  ExtractfromPack(position,data,x_,3*sizedouble);  
  // dofset_
  dofset_.Unpack(&data[position]);
  int dofsetsize = dofset_.SizePack(&data[position]);
  position += dofsetsize;
  // dentitytype_
  ExtractfromPack(position,data,dentitytype_);
  // dentityid_
  ExtractfromPack(position,data,dentityid_);
  // condition_
  int num=0;
  ExtractfromPack(position,data,num);
  for (int i=0; i<num; ++i)
  {
    string name;
    ExtractStringfromPack(position,data,name);
    RefCountPtr<Condition> cond = rcp(new Condition());
    cond->Unpack(&data[position]);
    int condsize = cond->SizePack(&data[position]);
    position += condsize;
    SetCondition(name,cond);
  }

  if (position != size)
    dserror("Mismatch in size of data %d <-> %d",size,position);
  return true;
}


/*----------------------------------------------------------------------*
 |  Get a condition of a certain name                          (public) |
 |                                                            gee 11/06 |
 *----------------------------------------------------------------------*/
const DRT::Condition* DRT::Node::GetCondition(const string& name) const
{
  map<string,RefCountPtr<Condition> >::const_iterator curr = condition_.find(name);
  if (curr != condition_.end()) return curr->second.get();
  else                          return NULL;
  return NULL;
}










#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
