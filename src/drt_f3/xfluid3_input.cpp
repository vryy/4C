/*!----------------------------------------------------------------------
\file xfluid3_input.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

extern "C"
{
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
}
#include "xfluid3.H"
#include "../drt_lib/drt_utils.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  read element input (public)                              g.bau 03/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::XFluid3::ReadElement()
{
    typedef map<string, DiscretizationType> Gid2DisType;
    Gid2DisType gid2distype;
    gid2distype["HEX8"]  = hex8;
    gid2distype["HEX20"] = hex20;
    gid2distype["HEX27"] = hex27;
    gid2distype["TET4"]  = tet4;
    gid2distype["TET10"] = tet10;
    gid2distype["WEDGE6"] = wedge6;
    gid2distype["WEDGE15"] = wedge15;
    gid2distype["PYRAMID5"] = pyramid5;

    typedef map<DiscretizationType, int> DisType2NumNodes;
    DisType2NumNodes distype2NumNodes;
    distype2NumNodes[hex8]  = 8;
    distype2NumNodes[hex20] = 20;
    distype2NumNodes[hex27] = 27;
    distype2NumNodes[tet4]  = 4;
    distype2NumNodes[tet10] = 10;
    distype2NumNodes[wedge6] = 6;
    distype2NumNodes[wedge15] = 15;
    distype2NumNodes[pyramid5] = 5;

    // read element's nodes
    int   ierr = 0;
    int   nnode = 0;
    int   nodes[27];
    DiscretizationType distype = dis_none;

    Gid2DisType::const_iterator iter;
    for( iter = gid2distype.begin(); iter != gid2distype.end(); iter++ )
    {
        const string eletext = iter->first;
        frchk(eletext.c_str(), &ierr);
        if (ierr == 1)
        {
            distype = gid2distype[eletext];
            nnode = distype2NumNodes[distype];
            frint_n(eletext.c_str(), nodes, nnode, &ierr);
            dsassert(ierr==1, "Reading of ELEMENT Topology failed\n");
            break;
        }
    }

    // reduce node numbers by one
    for (int i=0; i<nnode; ++i) nodes[i]--;

    SetNodeIds(nnode,nodes);

    // read number of material model
    int material = 0;
    frint("MAT",&material,&ierr);
    dsassert(ierr==1, "Reading of material for FLUID3 element failed\n");
    dsassert(material!=0, "No material defined for FLUID3 element\n");
    SetMaterial(material);
    
    // Initialize winding flags
    donerewinding_ = false;
    
    // read net algo
    char  buffer[50];
    frchar("NA",buffer,&ierr);
    if (ierr==1)
    {
        if (strncmp(buffer,"ale",3)==0 ||
            strncmp(buffer,"ALE",3)==0 ||
            strncmp(buffer,"Ale",3)==0 )
        {
            is_ale_ = true;
        }
        else if (strncmp(buffer,"euler",5)==0 ||
                 strncmp(buffer,"EULER",5)==0 ||
                 strncmp(buffer,"Euler",5)==0 )
            is_ale_ = false;
        else
            dserror("Reading of FLUID3 element failed: Euler/Ale\n");
    }
    else
        dserror("Reading of FLUID3 element net algorithm failed: NA\n");

  return true;

} // XFluid3::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
