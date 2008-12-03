/*!----------------------------------------------------------------------
\file xfluid3_input.cpp
\brief input stuff

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "xfluid3.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |  read element input (public)                              g.bau 03/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::XFluid3::ReadElement()
{
  if (genprob.ndim!=3)
    dserror("Not a 3d problem. Panic.");

    typedef map<string, DiscretizationType> Gid2DisType;
    Gid2DisType gid2distype;
    gid2distype.insert(make_pair("HEX8"     , hex8));
    gid2distype.insert(make_pair("HEX20"    , hex20));
    gid2distype.insert(make_pair("HEX27"    , hex27));
    gid2distype.insert(make_pair("TET4"     , tet4));
    gid2distype.insert(make_pair("TET10"    , tet10));
    gid2distype.insert(make_pair("WEDGE6"   , wedge6));
    gid2distype.insert(make_pair("WEDGE15"  , wedge15));
    gid2distype.insert(make_pair("PYRAMID5" , pyramid5));
    gid2distype.insert(make_pair("NURBS8"   , nurbs8));
    gid2distype.insert(make_pair("NURBS27"  , nurbs27));

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
    dsassert(ierr==1, "Reading of material for FLUID3 element failed\n");
    dsassert(material!=0, "No material defined for FLUID3 element\n");
    SetMaterial(material);


    // read gaussian points and set gaussrule
    // removed
    char  buffer[50];
    
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
            dserror("Reading of FLUID3 element failed: Euler/Ale\n");
    }
    else
        dserror("Reading of FLUID3 element net algorithm failed: NA\n");


  // input of ale and free surface related stuff is not supported
  // at the moment. TO DO!
  // see src/fluid3/f3_inpele.c for missing details


  return true;

} // XFluid3::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
