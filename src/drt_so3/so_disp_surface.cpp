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
#ifdef TRILINOS_PACKAGE

#include "so_disp.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 04/07|
 |  id             (in)  this element's global id                       |
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispSurface::SoDispSurface(int id, int owner,
                              int nnode, const int* nodeids,
                              DRT::Node** nodes,
                              DRT::Elements::SoDisp* parent,
                              const int lsurface) :
DRT::Element(id,element_sodispsurface,owner),
parent_(parent),
lsurface_(lsurface)
{
  SetNodeIds(nnode,nodeids);
  BuildNodalPointers(nodes);
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispSurface::SoDispSurface(const DRT::Elements::SoDispSurface& old) :
DRT::Element(old),
parent_(old.parent_),
lsurface_(old.lsurface_)
{
  return;
}

/*----------------------------------------------------------------------*
 |  Deep copy this instance return pointer to it               (public) |
 |                                                            maf 01/07 |
 *----------------------------------------------------------------------*/
DRT::Element* DRT::Elements::SoDispSurface::Clone() const
{
  DRT::Elements::SoDispSurface* newelement = new DRT::Elements::SoDispSurface(*this);
  return newelement;
}

/*----------------------------------------------------------------------*
 |                                                             (public) |
 |                                                          u.kue 03/07 |
 *----------------------------------------------------------------------*/
DRT::Element::DiscretizationType DRT::Elements::SoDispSurface::Shape() const
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
void DRT::Elements::SoDispSurface::Pack(vector<char>& data) const
{
  data.resize(0);
  dserror("this SoDispSurface element does not support communication");

  return;
}

/*----------------------------------------------------------------------*
 |  Unpack data                                                (public) |
 |                                                            maf 02/07 |
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispSurface::Unpack(const vector<char>& data)
{
  dserror("this SoDispSurface element does not support communication");
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/07|
 *----------------------------------------------------------------------*/
DRT::Elements::SoDispSurface::~SoDispSurface()
{
  return;
}

/*----------------------------------------------------------------------*
 |  get vector of lines (public)                             gammi 04/07|
 *----------------------------------------------------------------------*/
DRT::Element** DRT::Elements::SoDispSurface::Lines()
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

void DRT::Elements::SoDispSurface::CreateLinesTri(const int& nline,
                                                  const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];
        
        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[DRT::Utils::eleNodeNumbering_tri6_lines[iline][inode]];
             nodes[inode]   = Nodes()[  DRT::Utils::eleNodeNumbering_tri6_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::SoDispLine(iline,Owner(),nnode,nodeids,nodes,this,NULL,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}        

void DRT::Elements::SoDispSurface::CreateLinesQuad(const int& nline,
                                                   const int& nnode)
{
    for(int iline=0;iline<nline;iline++)
    {
        int nodeids[nnode];
        DRT::Node* nodes[nnode];
        
        for (int inode=0;inode<nnode;inode++)
        {
             nodeids[inode] = NodeIds()[DRT::Utils::eleNodeNumbering_quad9_lines[iline][inode]];
             nodes[inode]   = Nodes()[  DRT::Utils::eleNodeNumbering_quad9_lines[iline][inode]];
        }
        lines_[iline] = rcp(new DRT::Elements::SoDispLine(iline,Owner(),nnode,nodeids,nodes,this,NULL,iline));
        lineptrs_[iline] = lines_[iline].get();
    }
}    

/*----------------------------------------------------------------------*
 |  print this element (public)                                maf 01/07|
 *----------------------------------------------------------------------*/
void DRT::Elements::SoDispSurface::Print(ostream& os) const
{
  os << "SoDispSurface ";
  Element::Print(os);
  return;
}



#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif // #ifdef D_SOH8
