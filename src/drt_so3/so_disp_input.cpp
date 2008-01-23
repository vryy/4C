/*!----------------------------------------------------------------------
\file so_disp_input.cpp
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
#include "so_disp.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::Elements::SoDisp::ReadElement()
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
    
    const bool allowed_element = (distype == hex27) || (distype == hex20) || (distype == tet10) || (distype == wedge15);
    // The intention of this element is to help debugging the xfem intersection routines.
    // The element is purely displacement based and has not been tested for correct compution
    // It serves solely as higher order geometry input data to XFEm problems
    //dsassert(allowed_element, "Only quadratic order for DISP (displacement based) element.");

    // reduce node numbers by one
    for (int i=0; i<nnode; ++i) nodes[i]--;

    SetNodeIds(nnode,nodes);

    // read number of material model
    int material = 0;
    frint("MAT",&material,&ierr);
    dsassert(ierr==1, "Reading of material for SOLID3 element failed\n");
    dsassert(material!=0, "No material defined for SOLID3 element\n");
    SetMaterial(material);


    // read gaussian points and set gaussrule
    char  buffer[50];
    int ngp[3];
    switch (distype)
    {
    case hex8: case hex20: case hex27:
    {
        frint_n("GP",ngp,3,&ierr);
        dsassert(ierr==1, "Reading of SOLID3 element failed: GP\n");
        switch (ngp[0])
        {
        case 1:  
            gaussrule_ = DRT::UTILS::intrule_hex_1point; 
            break; 
        case 2:  
            gaussrule_ = DRT::UTILS::intrule_hex_8point; 
            break;
        case 3:  
            gaussrule_ = DRT::UTILS::intrule_hex_27point; 
            break;
        default:
            dserror("Reading of SOLID3 element failed: Gaussrule for hexaeder not supported!\n");
        }
        break;
    }
    case tet4: case tet10:
    {
//        gaussrule_ = DRT::UTILS::intrule_tet_4point;
        frint("GP_TET",&ngp[0],&ierr);
        dsassert(ierr==1, "Reading of SOLID3 element failed: GP_TET\n");

        frchar("GP_ALT",buffer,&ierr);
        dsassert(ierr==1, "Reading of SOLID3 element failed: GP_ALT\n");

        switch(ngp[0])
        {
        case 1:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = DRT::UTILS::intrule_tet_1point;
            else
                dserror("Reading of SOLID3 element failed: GP_ALT: gauss-radau not possible!\n");
            break;
        case 4:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = DRT::UTILS::intrule_tet_4point;
            else if (strncmp(buffer,"gaussrad",8)==0)
                gaussrule_ = DRT::UTILS::intrule_tet_4point_alternative;
            else
                dserror("Reading of SOLID3 element failed: GP_ALT\n");
            break;
        case 10:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = DRT::UTILS::intrule_tet_5point;
            else
                dserror("Reading of SOLID3 element failed: GP_ALT: gauss-radau not possible!\n");
            break;
        default:
            dserror("Reading of SOLID3 element failed: Gaussrule for tetraeder not supported!\n");
        }
        break;
    } // end reading gaussian points for tetrahedral elements
    default:
        dserror("Reading of SOLID3 element failed: integration points\n");
    } // end switch distype

    // read kinematic type
    frchar("KINEM",buffer,&ierr);
    if (ierr)
    {
     // geometrically linear
     if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = sodisp_geolin;
     // geometrically non-linear with Total Lagrangean approach
     else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = sodisp_totlag;
     // geometrically non-linear with Updated Lagrangean approach
     else if (strncmp(buffer,"Updlag",6)==0)
     {
         kintype_ = sodisp_updlag;
         dserror("Updated Lagrange for SOLID3 is not implemented!");
     }
     else dserror("Reading of SOLID3 element failed");
    }
    

    // read stress evaluation/output type
    frchar("STRESS",buffer,&ierr);
    if (ierr!=1) dserror("reading of SODISP stress failed");
    if (strncmp(buffer,"none",4)==0)  stresstype_= sodisp_stress_none;
    if (strncmp(buffer,"Gpxyz",5)==0) stresstype_= sodisp_stress_gpxyz;
    if (strncmp(buffer,"Gprst",5)==0) stresstype_= sodisp_stress_gprst;
    if (strncmp(buffer,"Gp123",5)==0) stresstype_= sodisp_stress_gp123;
    if (strncmp(buffer,"Ndxyz",5)==0) stresstype_= sodisp_stress_ndxyz;
    if (strncmp(buffer,"Ndrst",5)==0) stresstype_= sodisp_stress_ndrst;
    if (strncmp(buffer,"Nd123",5)==0) stresstype_= sodisp_stress_nd123;
    // set default: no stresses
    else stresstype_= sodisp_stress_none;

    numnod_disp_ = NumNode();      // number of nodes
    numdof_disp_ = NumNode() * NODDOF_DISP;     // total dofs per element
    const DRT::UTILS::IntegrationPoints3D  intpoints = DRT::UTILS::getIntegrationPoints3D(gaussrule_);
    numgpt_disp_ = intpoints.nquad;      // total gauss points per element


  return true;

} // SoDisp::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
