/*!----------------------------------------------------------------------
\file
\brief Inherit stabilisation condition to elements

<pre>
Maintainer: Christiane Foerster
            foerster@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/foerster/
            0711 - 685-6127
</pre>

*----------------------------------------------------------------------*/
/*!
\addtogroup INPUT
*//*! @{ (documentation module open)*/
#include "../headers/standardtypes.h"
#include "../fluid2/fluid2.h"
#include "../fluid3/fluid3.h"


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
    case el_fluid2:
#ifdef D_FLUID2
      /*------------------------------------------ set some pointers ---*/
      fluid2 = actele->e.f2;
      actdsurf = actele->g.gsurf->dsurf;

      /*------------------- get type of stabilisation to the element ---*/
      fluid2->stab_type = actdsurf->stab_type;

      /*--- assign pointer to corresponding stabilisation parameters ---*/
      if (fluid2->stab_type == stab_gls)
      {
        dsassert(actdsurf->stabi.gls!=NULL,"no stabilisation at DSURF!\n");
        fluid2->stabi.gls = actdsurf->stabi.gls;
      }
      else if (fluid2->stab_type == stab_prespro)
        ;   /* nothing needs to be done at the moment (no parameters!) */
      else
        dserror("Unknown stabilisation");
    break;
    /*---------------------------- do nothin in the following cases: ---*/
    case el_fluid2_pro:         /* projection method                    */
    case el_fluid2_tu:          /* turbulence elements                  */
#endif
    break;
    case el_fluid2_xfem:
#ifdef D_XFEM
      /*------------------------------------------ set some pointers ---*/
      fluid2 = actele->e.f2;
      actdsurf = actele->g.gsurf->dsurf;

      /*------------------- get type of stabilisation to the element ---*/
      fluid2->stab_type = actdsurf->stab_type;

      /*--- assign pointer to corresponding stabilisation parameters ---*/
      if (fluid2->stab_type == stab_gls)
        fluid2->stabi.gls = actdsurf->stabi.gls;
      else if (fluid2->stab_type == stab_prespro)
        ;   /* nothing needs to be done at the moment (no parameters!) */
      else
        dserror("Unknown stabilisation");
    break;
#endif
    break;
    case el_fluid3:
#ifdef D_FLUID3
      /*------------------------------------------ set some pointers ---*/
      fluid3 = actele->e.f3;
      actdvol = actele->g.gvol->dvol;

      /*------------------- get type of stabilisation to the element ---*/
      fluid3->stab_type = actdvol->stab_type;

      /*--- assign pointer to corresponding stabilisation parameters ---*/
      if (fluid3->stab_type == stab_gls)
      {
        dsassert(actdvol->stabi.gls!=NULL,"no stabilisation at DVOL!\n");
        fluid3->stabi.gls = actdvol->stabi.gls;
      }
      else
        dserror("Other than gls stabilisation not yet implemented!");
#endif
    break;
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

