/*!----------------------------------------------------------------------
\file combust3_input.cpp
\brief

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>

*----------------------------------------------------------------------*/
#ifdef D_FLUID3
#ifdef CCADISCRET

#include "combust3.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Combust3::ReadElement(const std::string& eletype,
                                          const std::string& distype,
                                          DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  return true;
}


#if 0
/*----------------------------------------------------------------------*
 |  read element input (public)                              g.bau 03/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Combust3::ReadElement()
{
    typedef map<string, DiscretizationType> Gid2DisType;
    Gid2DisType gid2distype;
    gid2distype["HEX8"]     = hex8;
    gid2distype["HEX20"]    = hex20;
    gid2distype["HEX27"]    = hex27;
    gid2distype["TET4"]     = tet4;
    gid2distype["TET10"]    = tet10;
    gid2distype["WEDGE6"]   = wedge6;
    gid2distype["WEDGE15"]  = wedge15;
    gid2distype["PYRAMID5"] = pyramid5;

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

    // read net algo
    char  buffer[50];
    frchar("NA",buffer,&ierr);
    if (ierr==1)
    {
        if (strncmp(buffer,"ale",3)==0 ||
            strncmp(buffer,"ALE",3)==0 ||
            strncmp(buffer,"Ale",3)==0 );
        else if (strncmp(buffer,"euler",5)==0 ||
                 strncmp(buffer,"EULER",5)==0 ||
                 strncmp(buffer,"Euler",5)==0 );
        else
            dserror("Reading of Combust3 element failed: Euler/Ale\n");
    }
    else
        dserror("Reading of Combust3 element net algorithm failed: NA\n");


  // input of ale and free surface related stuff is not supported
  // at the moment. TO DO!
  // see src/fluid3/f3_inpele.c for missing details


  return true;

} // XFluid3::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID3
