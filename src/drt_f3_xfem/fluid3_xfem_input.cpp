/*!----------------------------------------------------------------------
\file fluid3_xfem_input.cpp
\brief

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3_XFEM
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

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
#include "fluid3_xfem.H"

/*----------------------------------------------------------------------*
 |  read element input                                       a.ger 06/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::XFluid3::ReadElement()
{
    typedef map<string, DiscretizationType> Gid2DisType;
    Gid2DisType gid2distype;
    gid2distype["HEX8"]  = hex8;
    gid2distype["HEX20"] = hex20;
    gid2distype["HEX27"] = hex27;
    gid2distype["TET4"]  = tet4;
    gid2distype["TET10"] = tet10;

    typedef map<DiscretizationType, int> DisType2NumNodes;
    DisType2NumNodes distype2NumNodes;
    distype2NumNodes[hex8]  = 8;
    distype2NumNodes[hex20] = 20;
    distype2NumNodes[hex27] = 27;
    distype2NumNodes[tet4]  = 4;
    distype2NumNodes[tet10] = 10;
    
//    typedef vector<int> Volume2Surface(9);
//    typedef vector<Volume2Surface> Volume2Surfaces(6);
//    Volume2Surfaces vol2surf;
//    vol2surf[0][0] = 0;
//    vol2surf[0][1] = 3;
//    vol2surf[0][2] = 2;
//    vol2surf[0][3] = 1;
//    vol2surf[0][4] = 11;
//    vol2surf[0][5] = 10;
//    vol2surf[0][6] = 9;
//    vol2surf[0][7] = 8;
//    vol2surf[0][8] = 20;
    
    vector<int> gid2baciNodeNumbering(27);
    gid2baciNodeNumbering[0] =   1;
    gid2baciNodeNumbering[1] =   2;
    gid2baciNodeNumbering[2] =   3;
    gid2baciNodeNumbering[3] =   0;
    gid2baciNodeNumbering[4] =   5;
    gid2baciNodeNumbering[5] =   6;
    gid2baciNodeNumbering[6] =   7;
    gid2baciNodeNumbering[7] =   4;
    gid2baciNodeNumbering[8] =   9;
    gid2baciNodeNumbering[9] =  10;
    gid2baciNodeNumbering[10] = 11;
    gid2baciNodeNumbering[11] =  8;
    gid2baciNodeNumbering[12] = 17;
    gid2baciNodeNumbering[13] = 18;
    gid2baciNodeNumbering[14] = 19;
    gid2baciNodeNumbering[15] = 16;
    gid2baciNodeNumbering[16] = 13;
    gid2baciNodeNumbering[17] = 14;
    gid2baciNodeNumbering[18] = 15;
    gid2baciNodeNumbering[19] = 12;       
    gid2baciNodeNumbering[20] = 20;
    gid2baciNodeNumbering[21] = 25;
    gid2baciNodeNumbering[22] = 23;
    gid2baciNodeNumbering[23] = 24;
    gid2baciNodeNumbering[24] = 22;
    gid2baciNodeNumbering[25] = 21;
    gid2baciNodeNumbering[26] = 26;
    
    int   ierr = 0;
    int   nnode = 0;
    int   nodes_read[27];
    
    Gid2DisType::iterator iter;
    for( iter = gid2distype.begin(); iter != gid2distype.end(); iter++ ) 
    {
        const string eletext = iter->first;
        frchk(eletext.c_str(), &ierr);
        if (ierr == 1)
        {
            DiscretizationType distype = gid2distype[eletext];
            nnode = distype2NumNodes[distype];
            frint_n(eletext.c_str(), nodes_read, nnode, &ierr);
            dsassert(ierr==1, "Reading of ELEMENT Topology failed\n");
            break;
        }
    }
  
    // reduce node numbers by one
    for (int i=0; i<nnode; ++i) nodes_read[i]--;
    
    // map gid numbering to baci numbering
    int   nodes[27];
    for (int i=0; i<nnode; ++i) nodes[i] = nodes_read[gid2baciNodeNumbering[i]];
    //for (int i=0; i<nnode; ++i) nodes[i] = nodes_read[i];
    
    SetNodeIds(nnode,nodes);

    // read number of material model
    material_ = 0;
    frint("MAT",&material_,&ierr);
    dsassert(ierr==1, "Reading of material for XFLUID3 element failed\n");
    dsassert(material_!=0, "No material defined for XFLUID3 element\n");

    // read/set gaussian rule
    gaussrule_ = get_optimal_gaussrule(Shape());
    
    return true;
} // Fluid3::ReadElement()


#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3_XFEM
