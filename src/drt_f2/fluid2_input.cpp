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
#include "../drt_lib/drt_linedefinition.H"

using namespace DRT::UTILS;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Fluid2::ReadElement(const std::string& eletype,
                                        const std::string& distype,
                                        DRT::INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT",material);
  SetMaterial(material);

  if (distype=="THQ9")
    dismode_ = dismod_taylorhood;
  else
    dismode_ = dismod_equalorder;

  // The distype is determined by the number of node. This is the wrong way
  // around, but we do not need to switch on the distype string.
  DiscretizationType shape = Shape();

  gaussrule_ = getOptimalGaussrule(shape);

  std::string na;
  linedef->ExtractString("NA",na);
  if (na=="ale" or na=="ALE" or na=="Ale")
  {
    is_ale_ = true;
  }
  else if (na=="euler" or na=="EULER" or na=="Euler")
    is_ale_ = false;
  else
    dserror("Reading of FLUID3 element failed: Euler/Ale");

  return true;
}


#if 0
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
    gid2distype["THQ9"]	    = quad9;		// Taylor-Hood (Q2Q1)

    // discretization mode (equal-order, nonequal order)
    typedef map<string, DiscretizationMode> Gid2DisMode;
    Gid2DisMode gid2dismode;
    gid2dismode["QUAD4"]	= dismod_equalorder;
    gid2dismode["QUAD8"]	= dismod_equalorder;
    gid2dismode["QUAD9"]	= dismod_equalorder;
    gid2dismode["TRI3"] 	= dismod_equalorder;
    gid2dismode["TRI6"] 	= dismod_equalorder;
    gid2dismode["NURBS4"]	= dismod_equalorder;
    gid2dismode["NURBS9"]	= dismod_equalorder;
    gid2dismode["THQ9"]	 	= dismod_taylorhood;

    // read element's nodes
    int   ierr = 0;
    int   nnode = 0;
    int   nodes[9];
    DiscretizationType distype = dis_none;
    dismode_ = dismod_equalorder;		// standard

    Gid2DisType::const_iterator iter;
    for( iter = gid2distype.begin(); iter != gid2distype.end(); iter++ )
    {
        const string eletext = iter->first;
        frchk(eletext.c_str(), &ierr);
        if (ierr == 1)
        {
            distype = gid2distype[eletext];
            dismode_ = gid2dismode[eletext];	// equal or nonequal order element
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

    // set discretization mode for non-equal or equal order elements
    /*switch(distype)
    {
    case quad4: case tri3: case tri6: case nurbs9: case nurbs4: case quad8:
    {
    	dismode_=dismod_equal;
    	break;
    }
    case quad9:
    {
    	dismode_ = dismod_nonequal;
    	break;
    }
    default:
    	dserror("distype not recognized? TO DO: more accurate distinction of discretization types necessary. => Fluid2::DiscretizationMode");
    	break;
    }*/


  // input of ale and free surface related stuff is not supported
  // at the moment. TO DO!

  return true;

} // Fluid2::ReadElement()
#endif

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_FLUID2
