/*!---------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/


#ifndef VISUAL_PROTOTYPES_H
#define VISUAL_PROTOTYPES_H

/*----------------------------------------------------------------------*
 | global_visual.c                                       genk  07/02    |
 *----------------------------------------------------------------------*/
/*!---------------------------------------------------------------------
\brief call of Visualisation tools

<pre>                                                         genk 07/02

This routine checks the type of problem and based on the program  options
a visualisation tool is called.
At the moment implemented:
VISUAL2

</pre>
\return void

------------------------------------------------------------------------*/
void ntavisual(void);


/*----------------------------------------------------------------------*
 | ccarat_visual2.c                                        genk 07/02   |
 *----------------------------------------------------------------------*/
void vis2caf(
    INT      numf,
    INT      numa,
    INT      nums);

void v2movie(void);

void v2cell(
    FIELD   *actfield);

void v2_init(
    char     *titl,
    INT      *iopt,
    INT      *cmncol,
    char     *cmfile,
    INT      *cmunit,
    INT      *xypix,
    float    *xymin,
    float    *xymax,
    INT      *nkeys,
    INT      *ikeys,
    INT      *tkeys,
    INT      *fkeys,
    float   **flims,
    INT      *mnode,
    INT      *mptri,
    INT      *mpptri,
    INT      *mface,
    INT      *mpface,
    INT      *medge,
    INT      *mpedge);




/*----------------------------------------------------------------------*
 | ccarat_visual3.c                                        genk 01/04   |
 *----------------------------------------------------------------------*/
void vis3caf(
    INT     numff,
    INT     numaf,
    INT     numsf);



/*----------------------------------------------------------------------*
 | visual_readflaviares.c                                  genk 01/04   |
 *----------------------------------------------------------------------*/
void visual_readflaviares(
    FIELD   *actfield,
    INT     *ntsteps,
    ARRAY   *time_a,
    ARRAY   *step_a,
    INT     *FIRSTSTEP,
    INT     *LASTSTEP,
    INT     *DSTEP);



#endif
