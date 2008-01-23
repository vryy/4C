/*!----------------------------------------------------------------------
\file so_disp_surface.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOH8
#ifdef CCADISCRET

#include "so_disp.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::SoDispSurface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::ELEMENTS::SoDisp* parent,
                              const int lsurface) :
DRT::Element(id,element_sodispsurface,owner),
parent_(parent),
lsurface_(lsurface),
data_()
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}


/*----------------------------------------------------------------------*
 |  ctor (public)  for sending SoDispSurfaces                 umay 10/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::SoDispSurface(int id, int owner) :
DRT::Element(id,element_sodispsurface,owner),
parent_(),
lsurface_(),
data_()
{

  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::SoDispSurface(const DRT::ELEMENTS::SoDispSurface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_),
data_(old.data_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::ELEMENTS::SoDispSurface::Clone() const
{
  DRT::ELEMENTS::SoDispSurface* newelement = new DRT::ELEMENTS::SoDispSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::ELEMENTS::SoDispSurface::Shape() const
{
    switch (NumNode())
    {
    case 3: return tri3;
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
 |  Pack data                                                  (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispSurface::Pack(vector<char>& data) const
{
  data.resize(0);
  // pack type of this instance of ParObject
  int type = UniqueParObjectId();
  AddtoPack(data,type);
  // add base class Element
  vector<char> basedata(0);
  Element::Pack(basedata);
  AddtoPack(data,basedata);
  
  // data_
  vector<char> tmp(0);
  data_.Pack(tmp);
  AddtoPack(data,tmp);

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispSurface::Unpack(const vector<char>& data)
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
  // data_
  vector<char> tmp(0);
  ExtractfromPack(position,data,tmp);
  data_.Unpack(tmp);

  if (position != (int)data.size())
    dserror("Mismatch in size of data %d <-> %d",(int)data.size(),position);

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::SoDispSurface::~SoDispSurface()
{
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::ELEMENTS::SoDispSurface::Lines()
{
    const DiscretizationType distype = Shape();
    const int nline   = NumLine();
    lines_.resize(nline);
    lineptrs_.resize(nline);

    switch (distype)
    {
    case tri3:
        CreateLinesTri(3, 2);
        break;
    case tri6:
        CreateLinesTri(3, 3);
        break;
    case quad4:
        CreateLinesQuad(4, 2);
        break;
    case quad8:
        CreateLinesQuad(4, 3);
        break;
    case quad9:
        CreateLinesQuad(4, 3);
        break;
    default:
        dserror("distype not supported");
    }

    return (DRT::Element**)(&(lineptrs_[0]));
}

void DRT::ELEMENTS::SoDispSurface::CreateLinesTri(const int& nline,
                                                  const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];
        
        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[DRT::UTILS::eleNodeNumbering_tri6_lines[iline][inode]];
             nodes[inode]   = Nodes()[  DRT::UTILS::eleNodeNumbering_tri6_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::ELEMENTS::SoDispLine(iline,Owner(),nnode,nodeids,nodes,this,NULL,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}        

void DRT::ELEMENTS::SoDispSurface::CreateLinesQuad(const int& nline,
                                                   const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];
        
        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[DRT::UTILS::eleNodeNumbering_quad9_lines[iline][inode]];
             nodes[inode]   = Nodes()[  DRT::UTILS::eleNodeNumbering_quad9_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::ELEMENTS::SoDispLine(iline,Owner(),nnode,nodeids,nodes,this,NULL,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}    

/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDispSurface::Print(ostream& os) const
{
  os << "SoDispSurface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
