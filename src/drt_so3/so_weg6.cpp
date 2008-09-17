/*!----------------------------------------------------------------------
\file so_weg6.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_weg6.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/artwallremod.H"

using namespace DRT::UTILS;

// inverse design object
#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
#include "inversedesign.H"
#endif

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::So_weg6(int id, int owner) :
DRT::Element(id,element_so_weg6,owner),
data_()
{
  kintype_ = sow6_totlag;
  invJ_.resize(NUMGPT_WEG6);
  detJ_.resize(NUMGPT_WEG6);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  prestress_ = rcp(new DRT::ELEMENTS::PreStress(NUMNOD_WEG6,NUMGPT_WEG6));
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  invdesign_ = rcp(new DRT::ELEMENTS::InvDesign(NUMNOD_WEG6,NUMGPT_WEG6));
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::So_weg6(const DRT::ELEMENTS::So_weg6& old) :
DRT::Element(old),
kintype_(old.kintype_),
data_(old.data_),
detJ_(old.detJ_)
{
  invJ_.resize(old.invJ_.size());
  for (unsigned int i=0; i<invJ_.size(); ++i)
  {
    invJ_[i] = old.invJ_[i];
  }

#if defined(PRESTRESS) || defined(POSTSTRESS)
  prestress_ = rcp(new DRT::ELEMENTS::PreStress(*(old.prestress_)));
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  invdesign_ = rcp(new DRT::ELEMENTS::InvDesign(*(old.invdesign_)));
#endif

  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance of Solid3 and return pointer to it (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::So_weg6::Clone() const
{
  DRT::ELEMENTS::So_weg6* newelement = new DRT::ELEMENTS::So_weg6(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::So_weg6::Shape() const
{
  return wedge6;
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::Pack(vector<char>& data) const
{
  data.resize(0);

  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  // kintype_
  AddtoPack(data,kintype_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  // prestress_
  vector<char> tmpprestress(0);
  prestress_->Pack(tmpprestress);
  AddtoPack(data,tmpprestress);
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  // invdesign_
  vector<char> tmpinvdesign(0);
  invdesign_->Pack(tmpinvdesign);
  AddtoPack(data,tmpinvdesign);
#endif

  // detJ_
  AddtoPack(data,detJ_);

  // invJ_
  const unsigned int size = invJ_.size();
  AddtoPack(data,size);
  for (unsigned int i=0; i<size; ++i)
    AddtoPack(data,invJ_[i]);

  return;
}


/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::Unpack(const vector<char>& data)
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
  // kintype_
  ExtractfromPack(position,data,kintype_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

#if defined(PRESTRESS) || defined(POSTSTRESS)
  // prestress_
  vector<char> tmpprestress(0);
  ExtractfromPack(position,data,tmpprestress);
  prestress_->Unpack(tmpprestress);
#endif

#if defined(INVERSEDESIGNCREATE) || defined(INVERSEDESIGNUSE)
  // invdesign_
  vector<char> tmpinvdesign(0);
  ExtractfromPack(position,data,tmpinvdesign);
  invdesign_->Unpack(tmpinvdesign);
#endif

  // detJ_
  ExtractfromPack(position,data,detJ_);
  // invJ_
  int size;
  ExtractfromPack(position,data,size);
  invJ_.resize(size);
  for (int i=0; i<size; ++i)
    ExtractfromPack(position,data,invJ_[i]);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);
  return;
}


/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::~So_weg6()
{
  return;
}


/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::Print(ostream& os) const
{
  os << "So_weg6 ";
  Element::Print(os);
  cout << endl;
  cout << data_;
  return;
}

/*----------------------------------------------------------------------*
 |  extrapolation of quantities at the GPs to the nodes     maf 02/08   |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::soweg6_expol(LINALG::FixedSizeSerialDenseMatrix<NUMGPT_WEG6,NUMSTR_WEG6>& stresses,
                                          LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMSTR_WEG6>& nodalstresses)
{
  static LINALG::FixedSizeSerialDenseMatrix<NUMNOD_WEG6,NUMGPT_WEG6> expol;
  static bool isfilled;

  if (isfilled==true)
  {
    nodalstresses.Multiply(expol,stresses);
  }
  else
  {
   expol(0,0)=  -0.61004233964073;
   expol(0,1)=   0.12200846792815;
   expol(0,2)=   0.12200846792815;
   expol(0,3)=   2.27670900630740;
   expol(0,4)=  -0.45534180126148;
   expol(0,5)=  -0.45534180126148;
   expol(1,1)=  -0.61004233964073;
   expol(1,2)=   0.12200846792815;
   expol(1,3)=  -0.45534180126148;
   expol(1,4)=   2.27670900630740;
   expol(1,5)=  -0.45534180126148;
   expol(2,2)=  -0.61004233964073;
   expol(2,3)=  -0.45534180126148;
   expol(2,4)=  -0.45534180126148;
   expol(2,5)=   2.27670900630740;
   expol(3,3)=  -0.61004233964073;
   expol(3,4)=   0.12200846792815;
   expol(3,5)=   0.12200846792815;
   expol(4,4)=  -0.61004233964073;
   expol(4,5)=   0.12200846792815;
   expol(5,5)=  -0.61004233964073;
   for (int i=0;i<NUMNOD_WEG6;++i)
    {
      for (int j=0;j<i;++j)
      {
        expol(i,j)=expol(j,i);
      }
    }

    nodalstresses.Multiply(expol,stresses);
  }
}

/*----------------------------------------------------------------------*
 |  Return names of visualization data (public)                maf 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::VisNames(map<string,int>& names)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisNames(names);
  if (Material()->MaterialType() == m_artwallremod){
    string fiber = "Fiber1";
    names[fiber] = 3; // 3-dim vector
    fiber = "Fiber2";
    names[fiber] = 3; // 3-dim vector
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Return visualization data (public)                         maf 07/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::VisData(const string& name, vector<double>& data)
{
  // Put the owner of this element into the file (use base class method for this)
  DRT::Element::VisData(name,data);
  if (Material()->MaterialType() == m_artwallremod){
    MAT::ArtWallRemod* art = static_cast <MAT::ArtWallRemod*>(Material().get());
    vector<double> a1 = art->Geta1()->at(0);  // get a1 of first gp
    vector<double> a2 = art->Geta2()->at(0);  // get a2 of first gp
    if (name == "Fiber1"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a1[0]; data[1] = a1[1]; data[2] = a1[2];
    } else if (name == "Fiber2"){
      if ((int)data.size()!=3) dserror("size mismatch");
      data[0] = a2[0]; data[1] = a2[1]; data[2] = a2[2];
    } else if (name == "Owner"){
      if ((int)data.size()<1) dserror("Size mismatch");
      data[0] = Owner();
    } else {
      cout << name << endl;
      dserror("Unknown VisData!");
    }
 }


  return;
}


/*----------------------------------------------------------------------*
 |  allocate and return So_weg6Register (public)                maf 04/07|
 *----------------------------------------------------------------------*/
RefCountPtr<DRT::ElementRegister> DRT::ELEMENTS::So_weg6::ElementRegister() const
{
  return rcp(new DRT::ELEMENTS::Sow6Register(Type()));
}


/*----------------------------------------------------------------------*
 |  get vector of volumes (length 1) (public)                  maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_weg6::Volumes()
{
  vector<RCP<Element> > volumes(1);
  volumes[0]= rcp(this, false);
  return volumes;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_weg6::Surfaces()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralSurface,DRT::Element>(DRT::UTILS::buildSurfaces,this);
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
vector<RCP<DRT::Element> > DRT::ELEMENTS::So_weg6::Lines()
{
  // do NOT store line or surface elements inside the parent element
  // after their creation.
  // Reason: if a Redistribute() is performed on the discretization,
  // stored node ids and node pointers owned by these boundary elements might
  // have become illegal and you will get a nice segmentation fault ;-)

  // so we have to allocate new line elements:
  return DRT::UTILS::ElementBoundaryFactory<StructuralLine,DRT::Element>(DRT::UTILS::buildLines,this);
}


//=======================================================================
//=======================================================================
//=======================================================================
//=======================================================================

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Register::Sow6Register(DRT::Element::ElementType etype) :
ElementRegister(etype)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Register::Sow6Register(
                               const DRT::ELEMENTS::Sow6Register& old) :
ElementRegister(old)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Register* DRT::ELEMENTS::Sow6Register::Clone() const
{
  return new DRT::ELEMENTS::Sow6Register(*this);
}

/*----------------------------------------------------------------------*
 |  Pack data                                                  (public) |
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sow6Register::Pack(vector<char>& data) const
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
 |                                                            maf 04/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sow6Register::Unpack(const vector<char>& data)
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
 |  dtor (public)                                              maf 04/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::Sow6Register::~Sow6Register()
{
  return;
}

/*----------------------------------------------------------------------*
 |  print (public)                                             maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Sow6Register::Print(ostream& os) const
{
  os << "Sow6Register ";
  ElementRegister::Print(os);
  return;
}

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
