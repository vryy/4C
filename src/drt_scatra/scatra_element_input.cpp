/*----------------------------------------------------------------------*/
/*!
\file scatra_element_input.cpp
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
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport::ReadElement(const std::string& eletype,
                                           const std::string& distype,
                                           DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  SetDisType(DRT::StringToDistype(distype));

  return true;
}

#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                                gjb 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Transport::ReadElement()
{
    map<string, DiscretizationType> gid2distype;
    gid2distype["HEX8"]     = hex8;
    gid2distype["HEX20"]    = hex20;
    gid2distype["HEX27"]    = hex27;
    gid2distype["TET4"]     = tet4;
    gid2distype["TET10"]    = tet10;
    gid2distype["WEDGE6"]   = wedge6;
    gid2distype["WEDGE15"]  = wedge15;
    gid2distype["PYRAMID5"] = pyramid5;
    gid2distype["QUAD4"]    = quad4;
    gid2distype["QUAD8"]    = quad8;
    gid2distype["QUAD9"]    = quad9;
    gid2distype["TRI3"]     = tri3;
    gid2distype["TRI6"]     = tri6;
    gid2distype["NURBS4"]   = nurbs4;
    gid2distype["NURBS9"]   = nurbs9;
    gid2distype["LINE2"]    = line2;
    gid2distype["LINE3"]    = line3;

    // read element's nodes
    int   ierr = 0;
    int   nnode = 0;
    int   nodes[27];
    DiscretizationType distype = dis_none;

    map<string, DiscretizationType>::const_iterator iter;
    for( iter = gid2distype.begin(); iter != gid2distype.end(); iter++ )
    {
        const string eletext = iter->first;
        frchk(eletext.c_str(), &ierr);
        if (ierr == 1)
        {
            distype = gid2distype[eletext];
            // set the discretization type
            SetDisType(distype);
            nnode = DRT::UTILS::getNumberOfElementNodes(distype);
            frint_n(eletext.c_str(), nodes, nnode, &ierr);
            dsassert(ierr==1, "Reading of ELEMENT Topology failed\n");
            break;
        }
    }

    // reduce node numbers by one
    for (int i=0; i<nnode; ++i) nodes[i]--;

    // apply the node ids
    SetNodeIds(nnode,nodes);

    // read number of material model
    int material = 0;
    frint("MAT",&material,&ierr);
    dsassert(ierr==1, "Reading of material for Transport element failed\n");
    dsassert(material!=0, "No material defined for Transport element\n");
    // set material (and numdofpernode_) of transport element
    SetMaterial(material);


  return true;

} // Transport::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
#endif  // #if defined(D_FLUID2) || defined(D_FLUID3)
