#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "ale2.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"


DRT::Elements::Ale2::Ale2(int id, int owner)
  : DRT::Element(id,element_ale2,owner),
    material_(0),
    data_()
{
  ngp_[0] = ngp_[1] = 0;
  lines_.resize(0);
  lines_.resize(0);
}


DRT::Elements::Ale2::Ale2(const DRT::Elements::Ale2& old)
  : DRT::Element(old),
    material_(old.material_),
    data_(old.data_),
    lines_(old.lines_),
    lineptrs_(old.lineptrs_)
{
  for (int i=0; i<2; ++i)
    ngp_[i] = old.ngp_[i];
}


DRT::Element* DRT::Elements::Ale2::Clone() const
{
  DRT::Elements::Ale2* newelement = new DRT::Elements::Ale2(*this);
  return newelement;
}


DRT::Element::DiscretizationType DRT::Elements::Ale2::Shape() const
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


void DRT::Elements::Ale2::Pack(vector<char>& data) const
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
  AddtoPack(data,ngp_,2*sizeof(int));
  // material_
  AddtoPack(data,material_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);
}


void DRT::Elements::Ale2::Unpack(const vector<char>& data)
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
  ExtractfromPack(position,data,ngp_,2*sizeof(int));
  // material_
  ExtractfromPack(position,data,material_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
}


DRT::Elements::Ale2::~Ale2()
{
}


void DRT::Elements::Ale2::Print(ostream& os) const
{
  os << "Ale2 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}


RefCountPtr<DRT::ElementRegister> DRT::Elements::Ale2::ElementRegister() const
{
  return rcp(new DRT::Elements::Ale2Register(Type()));
}


DRT::Element** DRT::Elements::Ale2::Lines()
{

  const int nline   = NumLine();
  const int numnode = NumNode();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  int nodeids[100];
  DRT::Node* nodes[100];

  /* triangle elements */
  if (nline==3)
  {
    /* mind: The displayed element is not the reference element!
     *
     *               2
     *                X
     *                |\
     *                | \       edge1
     *   edge2       6o  o5    (line1)
     *  (line2)       |   \
     *                |    \
     *                X--o--X
     *               0   4   1
     *
     *             edge0
     *            (line0)
     *
     *
     *      X nodes for tri3
     *      o addotional nodes for tri6
     *                                                                       */
    /* linear triangles*/
    if (numnode==3)
    {
      /* first edge */
      // set node id's
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      // get nodes
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      // create the line and get the line pointer
      lines_[0] =
        rcp(new DRT::Elements::Ale2Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      /* second edge */
      // set node id's
      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      // get nodes
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      // create the line and get the line pointer
      lines_[1] =
        rcp(new DRT::Elements::Ale2Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      /* third edge */
      // set node id's
      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[0];
      // get nodes
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[0];
      // create the line and get the line pointer
      lines_[2] =
        rcp(new DRT::Elements::Ale2Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();
    }
    /* quadratic triangles*/
    else if (numnode==6)
    {
      dserror("TRI6 lines not implemented.");
    }
  }
  /* quad elements*/
  else if (nline==4)
  {
    /* mind: The displayed element is not the reference element!
     *
     *
     *                  edge2
     *		       (line2)
     *
     *              3     6     2
     *               X----o----X
     *               |         |
     *               |         |
     *   edge3      7o    O    o5      edge1
     *  (line3)      |    8    |      (line1)
     *               |         |
     *               X----o----X
     *              0     4     1
     *
     *                  edge0
     *                 (line0)
     *
     *
     *      X nodes for quad
     *      o addotional nodes for quad8
     *      O addotional nodes for full quadratic quad9
     *
     *                                                                 */
    if (numnode==4)
    {
      /* first edge */
      // set node id's
      nodeids[0] = NodeIds()[0];
      nodeids[1] = NodeIds()[1];
      // get nodes
      nodes[0] = Nodes()[0];
      nodes[1] = Nodes()[1];
      // create the line and get the line pointer
      lines_[0] =
        rcp(new DRT::Elements::Ale2Line(0,Owner(),2,nodeids,nodes,this,0));
      lineptrs_[0] = lines_[0].get();

      /* second edge */
      // set node id's
      nodeids[0] = NodeIds()[1];
      nodeids[1] = NodeIds()[2];
      // get nodes
      nodes[0] = Nodes()[1];
      nodes[1] = Nodes()[2];
      // create the line and get the line pointer
      lines_[1] =
        rcp(new DRT::Elements::Ale2Line(1,Owner(),2,nodeids,nodes,this,1));
      lineptrs_[1] = lines_[1].get();

      /* third edge */
      // set node id's
      nodeids[0] = NodeIds()[2];
      nodeids[1] = NodeIds()[3];
      // get nodes
      nodes[0] = Nodes()[2];
      nodes[1] = Nodes()[3];
      // create the line and get the line pointer
      lines_[2] =
        rcp(new DRT::Elements::Ale2Line(2,Owner(),2,nodeids,nodes,this,2));
      lineptrs_[2] = lines_[2].get();

      /* fourth edge */
      // set node id's
      nodeids[0] = NodeIds()[3];
      nodeids[1] = NodeIds()[0];
      // get nodes
      nodes[0] = Nodes()[3];
      nodes[1] = Nodes()[0];
      // create the line and get the line pointer
      lines_[3] =
        rcp(new DRT::Elements::Ale2Line(3,Owner(),2,nodeids,nodes,this,3));
      lineptrs_[3] = lines_[3].get();

    }
    else if (numnode==8)
    {
      dserror("quad8 lines not implemented.");
    }
    else if (numnode==9)
    {
      dserror("quad9 lines not implemented.");
    }
    else dserror("Number of nodes not supported");
  }
  else dserror("Number of lines not supported");

  return (DRT::Element**)(&(lineptrs_[0]));
}


DRT::Element** DRT::Elements::Ale2::Surfaces()
{
  surface_.resize(1);
  surface_[0] = this; //points to Ale2 element itself
  return &surface_[0];
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


DRT::Elements::Ale2Register::Ale2Register(DRT::Element::ElementType etype)
  : ElementRegister(etype)
{
}


DRT::Elements::Ale2Register::Ale2Register(const DRT::Elements::Ale2Register& old)
  : ElementRegister(old)
{
}


DRT::Elements::Ale2Register* DRT::Elements::Ale2Register::Clone() const
{
  return new DRT::Elements::Ale2Register(*this);
}


void DRT::Elements::Ale2Register::Pack(vector<char>& data) const
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


void DRT::Elements::Ale2Register::Unpack(const vector<char>& data)
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


DRT::Elements::Ale2Register::~Ale2Register()
{
}


void DRT::Elements::Ale2Register::Print(ostream& os) const
{
  os << "Ale2Register ";
  ElementRegister::Print(os);
}


#endif
#endif
#endif
