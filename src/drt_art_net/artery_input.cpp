/*!----------------------------------------------------------------------
\file artery_input.cpp
\brief

<pre>
Maintainer: Mahmoud Ismail
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15268
</pre>

*----------------------------------------------------------------------*/
#ifdef D_ARTNET
#ifdef CCADISCRET

#include "artery.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"

using namespace DRT::UTILS;

/*----------------------------------------------------------------------*
 |  read element input (public)                             ismail 03/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Artery::ReadElement()
{

  typedef map<string, DiscretizationType> Gid2DisType;
    Gid2DisType gid2distype;
    gid2distype["LINE2"]    = line2;

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
    dsassert(ierr==1, "Reading of material for ART element failed\n");
    dsassert(material!=0, "No material defined for ART  element\n");
    SetMaterial(material);


    // read gaussian points and set gaussrule
    char  buffer[50];
    int ngp[1];

    switch (distype)
    {
    case line2:
    {
      frint_n("GP",ngp,1,&ierr);
        dsassert(ierr==1, "Reading of ART element failed: GP\n");
        switch (ngp[0])
        {
        case 1:
           gaussrule_ = intrule_line_1point;
            break;
        case 2:
            gaussrule_ = intrule_line_2point;
            break;
        case 3:
            gaussrule_ = intrule_line_3point;
            break;
        case 4:
            gaussrule_ = intrule_line_4point;
            break;
        case 5:
            gaussrule_ = intrule_line_5point;
            break;
        case 6:
            gaussrule_ = intrule_line_6point;
            break;
        case 7:
            gaussrule_ = intrule_line_7point;
            break;
        case 8:
            gaussrule_ = intrule_line_8point;
            break;
        case 9:
            gaussrule_ = intrule_line_9point;
            break;
        case 10:
            gaussrule_ = intrule_line_10point;
            break;
        default:
            dserror("Reading of ART element failed: Gaussrule for line not supported!\n");
        }
        break;
    }
    default:
        dserror("Reading of ART element failed: integration points\n");
    } // end switch distype

    // read net algo
    frchar("TYPE",buffer,&ierr);
    if (ierr==1)
    {
        if (strncmp(buffer,"nonlin",6)==0 ||
                 strncmp(buffer,"Nonlin",6)==0 ||
                 strncmp(buffer,"NonLin",6)==0 )
        {
          //cout<<"Using nonlinear artery elements"<<endl;
        }
        else
            dserror("Reading of ART element failed: Lin/NonLin\n");
    }
    else
        dserror("Reading of ART element net algorithm failed: NA\n");


  // input of ale and free surface related stuff is not supported
  // at the moment. TO DO!


  return true;

} // Artery::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_ARTNET
