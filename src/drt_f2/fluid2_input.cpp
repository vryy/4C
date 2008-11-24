/*!----------------------------------------------------------------------
\file fluid2_input.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID2
#ifdef CCADISCRET

#include "fluid2.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  read element input (public)                              g.bau 03/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Fluid2::ReadElement()
{
    typedef map<string, DiscretizationType> Gid2DisType;
    Gid2DisType gid2distype;
    gid2distype["QUAD4"]    = quad4;
    gid2distype["QUAD8"]    = quad8;
    gid2distype["QUAD9"]    = quad9;
    gid2distype["TRI3"]     = tri3;
    gid2distype["TRI6"]     = tri6;
    gid2distype["NURBS4"]   = nurbs4;
    gid2distype["NURBS9"]   = nurbs9;

    // read element's nodes
    int   ierr = 0;
    int   nnode = 0;
    int   nodes[9];
    DiscretizationType distype = dis_none;

    Gid2DisType::const_iterator iter;
    for( iter = gid2distype.begin(); iter != gid2distype.end(); iter++ )
    {
        const string eletext = iter->first;
        frchk(eletext.c_str(), &ierr);
        if (ierr == 1)
        {
            distype = gid2distype[eletext];
            nnode = DRT::UTILS::getNumberOfElementNodes(distype);
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
    dsassert(ierr==1, "Reading of material for FLUID2 element failed\n");
    dsassert(material!=0, "No material defined for FLUID2 element\n");
    SetMaterial(material);


    // read gaussian points and set gaussrule
    char  buffer[50];
    int ngp[2];
    switch (distype)
    {
    case quad4: case quad8: case quad9: case nurbs4: case nurbs9:
    {
        frint_n("GP",ngp,2,&ierr);
        dsassert(ierr==1, "Reading of FLUID2 element failed: GP\n");
        switch (ngp[0])
        {
        case 1:
            gaussrule_ = intrule_quad_1point;
            break;
        case 2:
            gaussrule_ = intrule_quad_4point;
            break;
        case 3:
            gaussrule_ = intrule_quad_9point;
            break;
        default:
            dserror("Reading of FLUID2 element failed: Gaussrule for quad not supported!\n");
        }
        break;
    }
    case tri3: case tri6:
    {
        frint("GP_TRI",&ngp[0],&ierr);
        dsassert(ierr==1, "Reading of FLUID2 element failed: GP_TRI\n");

        frchar("GP_ALT",buffer,&ierr);
        dsassert(ierr==1, "Reading of FLUID2 element failed: GP_ALT\n");

        switch(ngp[0])
        {
        case 1:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = intrule_tri_1point;
            else
                dserror("Reading of FLUID2 element failed: GP_ALT: gauss-radau not possible!\n");
            break;
        case 3:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = intrule_tri_3point;
            else
                dserror("Reading of FLUID2 element failed: GP_ALT\n");
            break;
        case 6:
            if (strncmp(buffer,"standard",8)==0)
                gaussrule_ = intrule_tri_6point;
            else
                dserror("Reading of FLUID2 element failed: GP_ALT: gauss-radau not possible!\n");
            break;
        default:
            dserror("Reading of FLUID2 element failed: Gaussrule for triangle not supported!\n");
        }
        break;
    } // end reading gaussian points for triangular elements
    default:
        dserror("Reading of FLUID2 element failed: integration points\n");
    } // end switch distype

    // read net algo
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
            dserror("Reading of FLUID2 element failed: Euler/Ale\n");
    }
    else
        dserror("Reading of FLUID2 element net algorithm failed: NA\n");


  // input of ale and free surface related stuff is not supported
  // at the moment. TO DO!


  return true;

} // Fluid2::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
