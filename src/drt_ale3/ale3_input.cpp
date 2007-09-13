#ifdef D_ALE
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#ifdef PARALLEL
#include <mpi.h>
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

#include "ale3.H"
#include "../drt_lib/drt_utils.H"

using namespace DRT::Utils;


bool DRT::Elements::Ale3::ReadElement()
{
  cout << "input" <<endl;
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
    DiscretizationType distype;

    Gid2DisType::iterator iter;
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
  if (ierr!=1) dserror("Reading of ALE3 element failed\n");
  if (material==0) dserror("No material defined for ALE3 element\n");
  SetMaterial(material);

    // read gaussian points and set gaussrule (not used in the moment). Instead method GetOptimalGaussrule is used during
    // the evaluation
#if 0
    char  buffer[50];
    int ngp[3];
    switch (distype)
    {
    case hex8: case hex20: case hex27:
    {
        frint_n("GP",ngp,3,&ierr);
        dsassert(ierr==1, "Reading of FLUID3 element failed: GP\n");
        switch (ngp[0])
        {
        case 1:
            gaussrule_ = intrule_hex_1point;
            break;
        case 2:
	    gaussrule_ = intrule_hex_8point;
            break;
        case 3:
            gaussrule_ = intrule_hex_27point;
            break;
        default:
            dserror("Reading of FLUID3 element failed: Gaussrule for hexaeder not supported!\n");
        }
        break;
    }
    case tet4: case tet10:
    {
        frint("GP_TET",&ngp[0],&ierr);
        dsassert(ierr==1, "Reading of FLUID3 element failed: GP_TET\n");

        frchar("GP_ALT",buffer,&ierr);
        dsassert(ierr==1, "Reading of FLUID3 element failed: GP_ALT\n");

        switch(ngp[0])
        {
        case 1:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = intrule_tet_1point;
            else
                dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
            break;
        case 4:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = intrule_tet_4point;
            else if (strncmp(buffer,"gaussrad",8)==0)
                gaussrule_ = intrule_tet_4point_alternative;
            else
                dserror("Reading of FLUID3 element failed: GP_ALT\n");
            break;
        case 10:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = intrule_tet_10point;
            else
                dserror("Reading of FLUID3 element failed: GP_ALT: gauss-radau not possible!\n");
            break;
        default:
            dserror("Reading of FLUID3 element failed: Gaussrule for tetraeder not supported!\n");
        }
        break;
    } // end reading gaussian points for tetrahedral elements
    case wedge6: case wedge15:
    {
      frint("GP_WEDGE",&ngp[0],&ierr);
        dsassert(ierr==1, "Reading of FLUID3 element failed: GP_WEDGE\n");
      switch (ngp[0])
        {
        case 1:
            gaussrule_ = intrule_wedge_1point;
            break;
        case 6:
            gaussrule_ = intrule_wedge_6point;
            break;
        case 9:
            gaussrule_ = intrule_wedge_9point;
            break;
        default:
            dserror("Reading of FLUID3 element failed: Gaussrule for wedge  not supported!\n");
        }
        break;
    }

    case pyramid5:
    {
      frint("GP_PYRAMID",&ngp[0],&ierr);
        dsassert(ierr==1, "Reading of FLUID3 element failed: GP_PYRAMID\n");
      switch (ngp[0])
        {
        case 1:
            gaussrule_ = intrule_pyramid_1point;
            break;
        case 8:
            gaussrule_ = intrule_pyramid_8point;
            break;
        default:
            dserror("Reading of FLUID3 element failed: Gaussrule for pyramid  not supported!\n");
        }
        break;
    }
    default:
        dserror("Reading of FLUID3 element failed: integration points\n");
    } // end switch distype
#endif

  return true;
}


#endif
#endif
#endif
