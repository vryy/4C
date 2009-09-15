/*!
\file scatra_element.cpp
\brief

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#include "scatra_element.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_mat/matlist.H"
#include "../drt_lib/drt_globalproblem.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::Transport(int id, int owner) :
DRT::Element(id,element_transport,owner),
data_(),
numdofpernode_(-1),
distype_(dis_none)
{
    return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::Transport(const DRT::ELEMENTS::Transport& old) :
DRT::Element(old),
data_(old.data_),
numdofpernode_(old.numdofpernode_),
distype_(old.distype_)
{
    return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Transport and return pointer to it (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::Transport::Clone() const
{
  DRT::ELEMENTS::Transport* newelement = new DRT::ELEMENTS::Transport(*this);
  return newelement;
}


/*----------------------------------------------------------------------*
 |  create material class (public)                            gjb 07/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::SetMaterial(int matnum)
{
  // the standard part:
  //mat_ = MAT::Material::Factory(matnum);  // not allowed since mat_ is private
  DRT::Element::SetMaterial(matnum);

  // the special part:
  // now the element knows its material, and we can use it to determine numdofpernode
  RefCountPtr<MAT::Material> mat = Material();
  if(mat->MaterialType() == INPAR::MAT::m_scatra or
     mat->MaterialType() == INPAR::MAT::m_mixfrac_scatra or
     mat->MaterialType() == INPAR::MAT::m_sutherland_scatra or
     mat->MaterialType() == INPAR::MAT::m_arrhenius_pv_scatra)
  {
    numdofpernode_=1; // we only have a single scalar
  }
  else if (mat->MaterialType() == INPAR::MAT::m_matlist) // we have a system of scalars
  {
    const MAT::MatList* actmat = static_cast<const MAT::MatList*>(mat.get());
    numdofpernode_=actmat->NumMat();
  }
  else
    dserror("condif material expected but got type %d", mat->MaterialType());

  // for problem type ELCH we have one additional degree of freedom per node
  // for the electric potential
  if (DRT::Problem::Instance()->ProblemType()=="elch")
  {
    numdofpernode_ += 1;
    dsassert(numdofpernode_>2,"numdofpernode_ is not > 2 for ELCH problem");
  }

  return;
}


/*----------------------------------------------------------------------*
 |  Return the shape of a Transport element                      (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::Transport::Shape() const
{
  return distype_;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // numdofpernode
  AddtoPack(data,numdofpernode_);
  // distype
  AddtoPack(data,distype_);

  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Unpack(const vector<char>& data)
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
  // numdofpernode
  ExtractfromPack(position,data,numdofpernode_);
  // distype
  ExtractfromPack(position,data,distype_);

  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 05/08 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Transport::~Transport()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 05/08 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::Print(ostream& os) const
{
  os << "Transport element";
  Element::Print(os);
  cout << endl;
  cout << "DiscretizationType:  "<<distype_<<endl;
  cout << endl;
  cout << "Number DOF per Node: "<<numdofpernode_<<endl;
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  allocate and return register element (public)             gjb 05/08 |
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::Transport::ElementRegister() const
{
  //Assuming that this element do not need initialization, we return a
  //dummy base class here.
  return rcp(new DRT::ElementRegister(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of lines            (public)                  g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Transport::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  if (NumLine() > 1) // 3D and 2D
    return DRT::UTILS::ElementBoundaryFactory<TransportBoundary,Transport>(DRT::UTILS::buildLines,this);
  else
  {
    // 1D (we return the element itself)
    vector<RCP<Element> > lines(1);
    lines[0]= rcp(this, false);
    return lines;
  }

}


/*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                          g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Transport::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  if (NumSurface() > 1) // 3D
    return DRT::UTILS::ElementBoundaryFactory<TransportBoundary,Transport>(DRT::UTILS::buildSurfaces,this);
  else if (NumSurface() == 1)
  {
    // 2D (we return the element itself)
    vector<RCP<Element> > surfaces(1);
    surfaces[0]= rcp(this, false);
    return surfaces;
  }
  else
  {
    // 1D
    dserror("Surfaces() for 1D-Transport element not implemented");
    return DRT::Element::Surfaces();
  }
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                g.bau 03/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::Transport::Volumes()
{
  if (NumVolume() == 1)
  {
    vector<RCP<Element> > volumes(1);
    volumes[0]= rcp(this, false);
    return volumes;
  }
  else
  {
    dserror("Surfaces() for 1D-/2D-Transport element not implemented");
    return DRT::Element::Volumes();
  }
}


/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                gjb 01/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Transport::VisNames(map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);

  // see whether we have additional data for visualization in our container
  for (int k = 0 ;k<numdofpernode_; k++)
  {
    ostringstream temp;
    temp << k;

    // element Peclet number
    string name = "Pe_"+temp.str();
    const vector<double>* Pe = data_.Get<vector<double> >(name);
    if (Pe) names.insert(pair<string,int>(name,1));

    // element Peclet number (only migration term)
    name = "Pe_mig_"+temp.str();
    const vector<double>* Pe_mig = data_.Get<vector<double> >(name);
    if (Pe_mig) names.insert(pair<string,int>(name,1));

    //characteristic element length
    name = "hk_"+temp.str();
    const vector<double>* hk = data_.Get<vector<double> >(name);
    if (hk) names.insert(pair<string,int>(name,1));

    // Stabilization parameter at element center
    name = "tau_"+temp.str();
    const vector<double>* tau = data_.Get<vector<double> >(name);
    if (tau) names.insert(pair<string,int>(name,1));

  } // loop over transported scalars

  return;
}


/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         gjb 01/09|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport ::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  if(DRT::Element::VisData(name,data))
    return true;

  for (int k = 0 ;k<numdofpernode_; k++)
  {
    ostringstream temp;
    temp << k;
    if (   (name == "Pe_"+temp.str()    )
        || (name == "Pe_mig_"+temp.str())
        || (name == "hk_"+temp.str()    )
        || (name == "tau_"+temp.str()   )
    )
    {
      if ((int)data.size()!=1) dserror("size mismatch");
      const double value = data_.GetDouble(name);
      data[0] = value;
      return true;
    }
  } // loop over transported scalars

  return false;
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================


/*----------------------------------------------------------------------*
 |  ctor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::TransportBoundary(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::Transport* parent,
                              const int lbeleid) :
DRT::Element(id,element_transportboundary,owner),
parent_(parent),
lbeleid_(lbeleid)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::TransportBoundary(const DRT::ELEMENTS::TransportBoundary& old) :
DRT::Element(old),
parent_(old.parent_),
lbeleid_(old.lbeleid_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it     (public) gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::TransportBoundary::Clone() const
{
  DRT::ELEMENTS::TransportBoundary* newelement = new DRT::ELEMENTS::TransportBoundary(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |  Return shape of this element                    (public)  gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::TransportBoundary::Shape() const
{
  switch (NumNode())
  {
  case 2: return line2;
  case 3:
    if ((parent_->Shape() == quad8) || (parent_->Shape() == quad9))
      return line3;
    else
      return tri3;
  case 4: return quad4;
  case 6: return tri6;
  case 8: return quad8;
  case 9: return quad9;
  default:
    dserror("unexpected number of nodes %d", NumNode());
  }
  return dis_none;
}

/*----------------------------------------------------------------------*
 |  Pack data (public)                                        gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("This TransportBoundary element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data (public)                                      gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Unpack(const vector<char>& data)
{
  dserror("This TransportBoundary element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                             gjb 01/09 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::TransportBoundary::~TransportBoundary()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                               gjb 01/09 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::Print(ostream& os) const
{
  os << "TransportBoundary ";
  Element::Print(os);
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::TransportBoundary::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  dserror("Lines of TransportBoundary not implemented");
  vector<RCP<DRT::Element> > lines(0);
  return lines;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                              gjb 01/09 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::TransportBoundary::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new surface elements:
  dserror("Surfaces of TransportBoundary not implemented");
  vector<RCP<DRT::Element> > surfaces(0);
  return surfaces;
}


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
