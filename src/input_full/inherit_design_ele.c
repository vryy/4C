/*!----------------------------------------------------------------------
\file
\brief Inherit stabilisation condition to elements

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6572
</pre>

*----------------------------------------------------------------------*/
#ifndef CCADISCRET
/*!
\addtogroup INPUT
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"
#include "../fluid2_is/fluid2_is.h"
#include "../fluid3_is/fluid3_is.h"


#ifdef D_FLUID
/*!----------------------------------------------------------------------
\brief inherit stabilisation condition from design to elements

<pre>                                                       chfoe 01/04
The stabilisation parameters have been read within a condition assigned
to the dsurf (2D) or dvol (3D). Now pointers within the elements are set
to the stabilisation constants of the corresponding dsurface or dvolume.

This routine has to be called for discretisations of the fluid field only!
</pre>

\warning There is nothing special to this routine
\return void
\sa

*----------------------------------------------------------------------*/
void inherit_design_ele(DISCRET *actdis)
{

INT 	 i;
ELEMENT	*actele;

#ifdef D_FLUID2
FLUID2	*fluid2;
DSURF	*actdsurf;
#endif

#ifdef D_FLUID3
FLUID3	*fluid3;
DVOL	*actdvol;
#endif

#ifdef DEBUG
dstrc_enter("inherit_design_ele");
#endif
/*----------------------------------------------------------------------*/
for (i=0; i<actdis->numele; i++)
{
  actele = &(actdis->element[i]);
  switch (actele->eltyp)
  {
#ifdef D_FLUID2
    case el_fluid2:
      /*------------------------------------------ set some pointers ---*/
      fluid2 = actele->e.f2;
      actdsurf = actele->g.gsurf->dsurf;

      /*------------------- get type of stabilisation to the element ---*/
      fluid2->stab_type = actdsurf->stab_type;

      /*--- assign pointer to corresponding stabilisation parameters ---*/
      switch (fluid2->stab_type)
      {
        case stab_gls:
          dsassert(actdsurf->stabi.gls!=NULL,"no stabilisation at DSURF!\n");
          fluid2->stabi.gls = actdsurf->stabi.gls;
        break;
        case stab_usfem:
#ifdef D_FLUID2_TDS
	  case stab_tds:
        break;
#endif
        case stab_prespro:
         /* nothing needs to be done at the moment (no parameters!) */
        break;
        default:
          dserror("Unknown stabilisation");
      }
    break;
    /*---------------------------- do nothin in the following cases: ---*/
    case el_fluid2_tu:        /* turbulence elements                    */
    break;
#endif
#ifdef D_FLUID2_PRO
    case el_fluid2_pro:       /* projection method                      */
    break;
#endif
#ifdef D_FLUID3_PRO
    case el_fluid3_pro:       /* projection method                      */
    break;
#endif
#ifdef D_FLUID3
    case el_fluid3:
      /*------------------------------------------ set some pointers ---*/
      fluid3 = actele->e.f3;
      actdvol = actele->g.gvol->dvol;

      /*------------------- get type of stabilisation to the element ---*/
      fluid3->stab_type = actdvol->stab_type;

      /*--- assign pointer to corresponding stabilisation parameters ---*/
      switch (fluid3->stab_type)
      {
        case stab_gls:
        dsassert(actdvol->stabi.gls!=NULL,"no stabilisation at DVOL!\n");
        fluid3->stabi.gls = actdvol->stabi.gls;
        break;
        case stab_usfem:
         /* nothing needs to be done at the moment (no parameters!) */
        break;
#ifdef D_FLUID3_TDS
	  case stab_tds:
        break;
#endif
        default:
        dserror("Unknown stabilisation for fluid3!");
      }
#endif
    break;

#ifdef D_FLUID3_F
    case el_fluid3_fast:
      /*------------------------------------------ set some pointers ---*/
      fluid3 = actele->e.f3;
      actdvol = actele->g.gvol->dvol;

      /*------------------- get type of stabilisation to the element ---*/
      fluid3->stab_type = actdvol->stab_type;

      /*--- assign pointer to corresponding stabilisation parameters ---*/
      switch (fluid3->stab_type)
      {
        case stab_gls:
        dsassert(actdvol->stabi.gls!=NULL,"no stabilisation at DVOL!\n");
        fluid3->stabi.gls = actdvol->stabi.gls;
        break;
        case stab_usfem:
         /* nothing needs to be done at the moment (no parameters!) */
        break;
#ifdef D_FLUID3_TDS
	  case stab_tds:
        break;
#endif
        default:
        dserror("Unknown stabilisation for fluid3f!");
      }
    break;
#endif

#ifdef D_FLUID2_IS
  case el_fluid2_is:
  {
    FLUID2_IS* f2is;
    f2is = actele->e.f2is;
    actdsurf = actele->g.gsurf->dsurf;

    /*------------------- get type of stabilisation to the element ---*/
    f2is->stab_type = actdsurf->stab_type;

    /*--- assign pointer to corresponding stabilisation parameters ---*/
    switch (f2is->stab_type)
    {
    case stab_gls:
      dsassert(actdsurf->stabi.gls!=NULL,"no stabilisation at DSURF!\n");
      f2is->stabi.gls = actdsurf->stabi.gls;
      break;
    case stab_usfem:
#ifdef D_FLUID2_TDS
	  case stab_tds:
        break;
#endif
    case stab_prespro:
      /* nothing needs to be done at the moment (no parameters!) */
      break;
    default:
      dserror("Unknown stabilisation");
    }
    break;
  }
#endif

#ifdef D_FLUID3_IS
  case el_fluid3_is:
  {
    FLUID3_IS* f3is;
    f3is = actele->e.f3is;
    actdvol = actele->g.gvol->dvol;

    /*------------------- get type of stabilisation to the element ---*/
    f3is->stab_type = actdvol->stab_type;

    /*--- assign pointer to corresponding stabilisation parameters ---*/
    switch (f3is->stab_type)
    {
    case stab_gls:
      dsassert(actdvol->stabi.gls!=NULL,"no stabilisation at DSURF!\n");
      f3is->stabi.gls = actdvol->stabi.gls;
      break;
    case stab_usfem:
#ifdef D_FLUID3_TDS
	  case stab_tds:
        break;
#endif
    case stab_prespro:
      /* nothing needs to be done at the moment (no parameters!) */
      break;
    default:
      dserror("Unknown stabilisation");
    }
    break;
  }
#endif

    default: dserror("Unknown element (eltyp) in fluid field!");
  }
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of inherit_design_ele */

#endif
/*! @} (documentation module close)*/
#endif
