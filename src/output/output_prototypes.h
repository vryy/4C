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


#ifndef OUTPUT_PROTOTYPES_H
#define OUTPUT_PROTOTYPES_H

/*----------------------------------------------------------------------*
 |  out_global.c                                         m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_general(void);
void out_sol(FIELD *actfield, PARTITION *actpart, INTRA *actintra, 
             INT step, INT place);

void out_fluidmf(FIELD *fluidfield);
void out_fsi(FIELD *fluidfield);
void out_fluidtu(FIELD *actfield, INTRA *actintra, INT step, INT place);
/*----------------------------------------------------------------------*
 |  out_gid_sol.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_sol_init(void);
void out_gid_domains(FIELD *actfield);
void out_gid_sol(char string[], FIELD *actfield, INTRA  *actintra, INT step,
                 INT place, DOUBLE time);
void out_fsi(FIELD *fluidfield);
void out_gid_sol_fsi(FIELD *fluidfield, FIELD *structfield);
/*----------------------------------------------------------------------*
 |  out_gid_soldyn.c                                     m.gee 5/03     |
 *----------------------------------------------------------------------*/
void out_gid_soldyn(char string[], FIELD *actfield, INTRA  *actintra, INT step,
                   INT place, DOUBLE totaltime);
/*----------------------------------------------------------------------*
 |  out_gid_msh.c                                        m.gee 12/01    |
 *----------------------------------------------------------------------*/
void out_gid_msh(void);
void out_gid_allcoords(FILE *out);
/*----------------------------------------------------------------------*
 |  out_monitor.c                                         genk 01/03    |
 *----------------------------------------------------------------------*/
void out_monitor(FIELD *actfield, INT numf,DOUBLE time,INT init);
void out_area(ARRAY totarea_a, DOUBLE time, INT itnum, INT init);

/*----------------------------------------------------------------------*
 |  out_checkfilesize.c                                   genk 08/03    |
 *----------------------------------------------------------------------*/
void out_checkfilesize(INT opt); 
/*----------------------------------------------------------------------*
 |  out_plt.c                                            chfoe 01/04    |
 *----------------------------------------------------------------------*/
void plot_liftdrag(DOUBLE time, DOUBLE *liftdrag);
void plot_lte(	DOUBLE  time, 
                INT     step, 
                DOUBLE  norm, 
                DOUBLE  dt, 
                INT     itnum);
void plot_ale_quality(FIELD *field,INT step, INTRA *actintra, 
                      PARTITION *actpart);

/*----------------------------------------------------------------------*
 |  out_gid_solfsi.c                                       mn 05/03     |
 *----------------------------------------------------------------------*/       
void out_gid_sol_fsi(
    FIELD     *fluidfield,
    FIELD     *structfield);


#endif

