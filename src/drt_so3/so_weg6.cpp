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
#include "../drt_lib/drt_utils_fem_shapefunctions.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::So_weg6::So_weg6(int id, int owner) :
DRT::Element(id,element_so_weg6,owner),
data_()
{
  volume_.resize(0);
  surfaces_.resize(0);
  surfaceptrs_.resize(0);
  lines_.resize(0);
  lineptrs_.resize(0);
  kintype_ = sow6_totlag;
  donerewinding_ = false;
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
volume_(old.volume_),
surfaces_(old.surfaces_),
surfaceptrs_(old.surfaceptrs_),
lines_(old.lines_),
lineptrs_(old.lineptrs_),
donerewinding_(old.donerewinding_)
{
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
  // rewind flags
  AddtoPack(data,donerewinding_);
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

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
  // rewinding flags
  ExtractfromPack(position,data,donerewinding_);
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

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
void DRT::ELEMENTS::So_weg6::soweg6_expol(Epetra_SerialDenseMatrix& stresses,
                                          Epetra_SerialDenseMatrix& nodalstresses)
{
  static Epetra_SerialDenseMatrix expol(NUMNOD_WEG6,NUMGPT_WEG6);
  static bool isfilled;

  if (isfilled==true)
  {
    nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);
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

    nodalstresses.Multiply('N','N',1.0,expol,stresses,0.0);
  }
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
DRT::Element** DRT::ELEMENTS::So_weg6::Volumes()
{
  volume_.resize(1);
  return 0;
}

 /*----------------------------------------------------------------------*
 |  get vector of surfaces (public)                             maf 04/07|
 |  surface normals always point outward                                 |
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::So_weg6::Surfaces()
{
  const int nsurf = NumSurface();
  surfaces_.resize(nsurf);
  surfaceptrs_.resize(nsurf);
  // first the 3 quad surfaces (#0..2)
  for (int qisurf = 0; qisurf < 3; ++qisurf) 
  {
        const int nnode_surf = 4;
        const int surfid = qisurf;
        int nodeids[nnode_surf];
        DRT::Node* nodes[nnode_surf];
        for (int qinode = 0; qinode < nnode_surf; ++qinode) {
          nodeids[qinode] = NodeIds()[eleNodeNumbering_wedge15_quadsurfaces[surfid][qinode]];
          nodes[qinode] = Nodes()[eleNodeNumbering_wedge15_quadsurfaces[surfid][qinode]];
        }
        surfaces_[qisurf] = rcp(new DRT::ELEMENTS::Sow6Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
        surfaceptrs_[qisurf] = surfaces_[qisurf].get();
  };
  // then the tri's...
  for (int tisurf = 0; tisurf < 2; ++tisurf) 
  {
      const int nnode_surf = 3;
      const int surfid = tisurf + 3;
      int nodeids[nnode_surf];
      DRT::Node* nodes[nnode_surf];
      for (int tinode = 0; tinode < nnode_surf; ++tinode) {
        nodeids[tinode] = NodeIds()[eleNodeNumbering_wedge15_trisurfaces[tisurf][tinode]];
        nodes[tinode] = Nodes()[eleNodeNumbering_wedge15_trisurfaces[tisurf][tinode]];
      }
      surfaces_[surfid] = rcp(new DRT::ELEMENTS::Sow6Surface(surfid,Owner(),nnode_surf,nodeids,nodes,this,surfid));
      surfaceptrs_[surfid] = surfaces_[surfid].get();
  };
  return (DRT::Element**)(&(surfaceptrs_[0]));
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                               maf 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::So_weg6::Lines()
{
  dserror("So_weg6 lines not yet implemented");
  const int nline = NumLine();
  lines_.resize(nline);
  lineptrs_.resize(nline);
  return (DRT::Element**)(&(lineptrs_[0]));
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
