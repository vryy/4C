#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "ale3.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


DRT::Elements::Ale3::Ale3(int id, int owner)
  : DRT::Element(id,element_ale3,owner),
    material_(0),
    data_()
{
  ngp_[0] = ngp_[1] = ngp_[2] = 0;
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
}


DRT::Elements::Ale3::Ale3(const DRT::Elements::Ale3& old)
  : DRT::Element(old),
    material_(old.material_),
    data_(old.data_),
    surfaces_(old.surfaces_),
    surfaceptrs_(old.surfaceptrs_)
{
  for (int i=0; i<3; ++i)
    ngp_[i] = old.ngp_[i];
}


DRT::Element* DRT::Elements::Ale3::Clone() const
{
  DRT::Elements::Ale3* newelement = new DRT::Elements::Ale3(*this);
  return newelement;
}


DRT::Element::DiscretizationType DRT::Elements::Ale3::Shape() const
{
  switch (NumNode())
  {
  case  4: return tet4;
  case  8: return hex8;
  case 10: return tet10;
  case 20: return hex20;
  case 27: return hex27;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}


void DRT::Elements::Ale3::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // ngp_
  AddtoPack(data,ngp_,3*sizeof(int));
  // material_
  AddtoPack(data,material_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);
}


void DRT::Elements::Ale3::Unpack(const vector<char>& data)
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
  // ngp_
  ExtractfromPack(position,data,ngp_,3*sizeof(int));
  // material_
  ExtractfromPack(position,data,material_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


DRT::Elements::Ale3::~Ale3()
{
}


void DRT::Elements::Ale3::Print(ostream& os) const
{
  os << "Ale3 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


RefCountPtr<DRT::ElementRegister> DRT::Elements::Ale3::ElementRegister() const
{
  return rcp(new DRT::Elements::Ale3Register(Type()));
}


DRT::Element** DRT::Elements::Ale3::Surfaces()
{

  const int nsurf = NumSurface();
  const int numnode = NumNode();
  surfaces_.resize(nsurf);
  surfaceptrs_.resize(nsurf);
  int nodeids[100];
  DRT::Node* nodes[100];

  if (nsurf==4)
  {
    if (numnode==4)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[2];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[2];
      surfaces_[0] =
        rcp(new DRT::Elements::Ale3Surface(0,Owner(),3,nodeids,nodes,this,0));
      surfaceptrs_[0] = surfaces_[0].get();

      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[3];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[3];
      surfaces_[1] =
        rcp(new DRT::Elements::Ale3Surface(1,Owner(),3,nodeids,nodes,this,1));
      surfaceptrs_[1] = surfaces_[1].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      nodeids[2] = NodeIds()[3];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[0];
      nodes[2] = Nodes()[3];
      surfaces_[2] =
        rcp(new DRT::Elements::Ale3Surface(2,Owner(),3,nodeids,nodes,this,2));
      surfaceptrs_[2] = surfaces_[2].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[3];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[3];
      surfaces_[3] =
        rcp(new DRT::Elements::Ale3Surface(3,Owner(),3,nodeids,nodes,this,3));
      surfaceptrs_[3] = surfaces_[3].get();

    }
    else if (numnode==10)
    {
      dserror("TET10 surfaces not implemented.");
    }
    else dserror("Number of nodes not supported");
  }
  else if (nsurf==6)
  {
    if (numnode==8)
    {
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[3];
      nodeids[2] = NodeIds()[2];
      nodeids[3] = NodeIds()[1];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[3];
      nodes[2] = Nodes()[2];
      nodes[3] = Nodes()[1];
      surfaces_[0] =
        rcp(new DRT::Elements::Ale3Surface(0,Owner(),4,nodeids,nodes,this,0));
      surfaceptrs_[0] = surfaces_[0].get();

      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      nodeids[2] = NodeIds()[5];
      nodeids[3] = NodeIds()[4];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      nodes[2] = Nodes()[5];
      nodes[3] = Nodes()[4];
      surfaces_[1] =
        rcp(new DRT::Elements::Ale3Surface(1,Owner(),4,nodeids,nodes,this,1));
      surfaceptrs_[1] = surfaces_[1].get();

      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[4];
      nodeids[2] = NodeIds()[7];
      nodeids[3] = NodeIds()[3];
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[4];
      nodes[2] = Nodes()[7];
      nodes[3] = Nodes()[3];
      surfaces_[2] =
        rcp(new DRT::Elements::Ale3Surface(2,Owner(),4,nodeids,nodes,this,2));
      surfaceptrs_[2] = surfaces_[2].get();

      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      nodeids[2] = NodeIds()[7];
      nodeids[3] = NodeIds()[6];
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      nodes[2] = Nodes()[7];
      nodes[3] = Nodes()[6];
      surfaces_[3] =
        rcp(new DRT::Elements::Ale3Surface(3,Owner(),4,nodeids,nodes,this,3));
      surfaceptrs_[3] = surfaces_[3].get();

      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      nodeids[2] = NodeIds()[6];
      nodeids[3] = NodeIds()[5];
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      nodes[2] = Nodes()[6];
      nodes[3] = Nodes()[5];
      surfaces_[4] =
        rcp(new DRT::Elements::Ale3Surface(4,Owner(),4,nodeids,nodes,this,4));
      surfaceptrs_[4] = surfaces_[4].get();

      nodeids[0] = NodeIds()[4];
      nodeids[1] = NodeIds()[5];
      nodeids[2] = NodeIds()[6];
      nodeids[3] = NodeIds()[7];
      nodes[0] = Nodes()[4];
      nodes[1] = Nodes()[5];
      nodes[2] = Nodes()[6];
      nodes[3] = Nodes()[7];
      surfaces_[5] =
        rcp(new DRT::Elements::Ale3Surface(5,Owner(),4,nodeids,nodes,this,5));
      surfaceptrs_[5] = surfaces_[5].get();
    }
    else if (numnode==20)
    {
      dserror("hex20 surfaces not implemented.");
    }
    else if (numnode==27)
    {
      dserror("hex27 surfaces not implemented.");
    }
    else dserror("Number of nodes not supported");
  }
  else dserror("Number of lines not supported");

  return (DRT::Element**)(&(surfaceptrs_[0]));
}


DRT::Element** DRT::Elements::Ale3::Volumes()
{
  volume_.resize(1);
  volume_[0] = this; //points to Ale3 element itself
  return &volume_[0];
}


void DRT::Elements::Ale3::SetGaussPoints()
{
  switch (Shape())
  {
  case hex8:
    ngp_[0] = 2; ngp_[1] = 2; ngp_[2] = 2;
    break;
  case hex20:
  case hex27:
    ngp_[0] = 3; ngp_[1] = 3; ngp_[2] = 3;
    break;
  default:
    dserror("unsupported shape %d for gauss point selection", Shape());
  }
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


DRT::Elements::Ale3Register::Ale3Register(DRT::Element::ElementType etype)
  : ElementRegister(etype)
{
}


DRT::Elements::Ale3Register::Ale3Register(const DRT::Elements::Ale3Register& old)
  : ElementRegister(old)
{
}


DRT::Elements::Ale3Register* DRT::Elements::Ale3Register::Clone() const
{
  return new DRT::Elements::Ale3Register(*this);
}


void DRT::Elements::Ale3Register::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class ElementRegister
  vector<char> basedata(0);
  ElementRegister::Pack(basedata);
  AddtoPack(data,basedata);
}


void DRT::Elements::Ale3Register::Unpack(const vector<char>& data)
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
}


DRT::Elements::Ale3Register::~Ale3Register()
{
}


void DRT::Elements::Ale3Register::Print(ostream& os) const
{
  os << "Ale3Register ";
  ElementRegister::Print(os);
}


#endif
#endif
#endif
