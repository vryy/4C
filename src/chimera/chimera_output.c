#ifdef D_CHIMERA
#include "../headers/standardtypes.h"
/*#include "stdio.h"*/
/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
  *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;
/*----------------------------------------------------------
Write output to binary files --- may be plotted with pager
----------------------------------------------------------*/
void chimera_output (
  INT              numdis
  )
{
    INT  *ien,rr,mm,binaer=0;
    double *xy,*data,sdc;
    FILE  *outputfile;

#ifdef DEBUG
    dstrc_enter("chimera_output");
#endif
    for (binaer=0;binaer<2;binaer++){
	if (field[genprob.numff].dis[numdis].element[0].distyp == quad4)
	{
/*---------------------------------------- Tempor"are Zwischenspeicher */
	    ien=(int*)malloc(4*field[0].dis[numdis].numele*sizeof(int));
	    for (rr=0;rr<field[0].dis[numdis].numele;rr++)
	    {
		for (mm=0;mm<4;mm++)
		{
		    ien[4*rr+mm]=1+field[0].dis[numdis].element[rr].node[mm][0].Id_loc;
		}
	    }

	    xy=(double*)malloc(2*field[0].dis[numdis].numnp*sizeof(double));
	    for (rr=0;rr<field[0].dis[numdis].numnp;rr++)
	    {
		for (mm=0;mm<2;mm++)
		{
		    xy[2*rr+mm]=field[0].dis[numdis].node[rr].x[mm];
		}
	    }

	    data=(double*)malloc(3*field[0].dis[numdis].numnp*sizeof(double));
	    for (rr=0;rr<field[0].dis[numdis].numnp;rr++)
	    {
		for (mm=0;mm<3;mm++)
		{
		    data[3*rr+mm]=field[0].dis[numdis].node[rr].sol_increment.a.da[1][mm];
		}
		if (field[0].dis[numdis].node[rr].gnode[0].Knotentyp==Loch_Standard)
		{
		    for (mm=0;mm<3;mm++)
		    {
			data[3*rr+mm]=0;
		    }
		}
	    }
	    if (binaer==0) {
/*---------------------------------------------------------- write ien */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien0","w""b");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien1","w""b");
		}
		fwrite(&ien[0],sizeof(int),4*field[0].dis[numdis].numele,outputfile);
		fclose(outputfile);
/*-------------------------------------------------- write coordinates */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy0","wb");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy1","wb");
		}
		fwrite(&xy[0],sizeof(double),2*field[0].dis[numdis].numnp,outputfile);
		fclose(outputfile);
/*------------------------------------------------------- write values */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_data0","w""b");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_data1","w""b");
		}
		fwrite(&data[0],sizeof(double),3*field[0].dis[numdis].numnp,outputfile);
		fclose(outputfile);
	    }
	    else
	    {
/*---------------------------------------------------------- write ien */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien0asci","w");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien1asci","w");
		}
		for (rr=0;rr<field[0].dis[numdis].numele;rr++){
		    for (mm=0;mm<4;mm++){
			fprintf(outputfile," %14d   ",ien[4*rr+mm]);
		    }
		    fprintf(outputfile,"\n");
		}
		fclose(outputfile);
/*-------------------------------------------------- write coordinates */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy0asci","w");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy1asci","w");
		}
		for (rr=0;rr<field[0].dis[numdis].numnp;rr++){
		    for (mm=0;mm<2;mm++){
			fprintf(outputfile," %14g   ",xy[2*rr+mm]);
		    }
		    fprintf(outputfile,"\n");
		}
		fclose(outputfile);
/*------------------------------------------------------- write values */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_data0asci","w");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_data1asci","w");
		}
		for (rr=0;rr<field[0].dis[numdis].numnp;rr++){
		    for (mm=0;mm<3;mm++){
			fprintf(outputfile," %14g   ",data[3*rr+mm]);
		    }
		    fprintf(outputfile,"\n");
		}
		fclose(outputfile);
	    }
	}
	if (field[genprob.numff].dis[numdis].element[0].distyp == tri3)
	{
/*---------------------------------------- Tempor"are Zwischenspeicher */
	    ien=(int*)malloc(3*field[0].dis[numdis].numele*sizeof(int));
	    for (rr=0;rr<field[0].dis[numdis].numele;rr++)
	    {
		for (mm=0;mm<3;mm++)
		{
		    ien[3*rr+mm]=1+field[0].dis[numdis].element[rr].node[mm][0].Id_loc;
		}
	    }

	    xy=(double*)malloc(2*field[0].dis[numdis].numnp*sizeof(double));
	    for (rr=0;rr<field[0].dis[numdis].numnp;rr++)
	    {
		for (mm=0;mm<2;mm++)
		{
		    xy[2*rr+mm]=field[0].dis[numdis].node[rr].x[mm];
		}
	    }

	    data=(double*)malloc(3*field[0].dis[numdis].numnp*sizeof(double));
	    for (rr=0;rr<field[0].dis[numdis].numnp;rr++)
	    {
		for (mm=0;mm<3;mm++)
		{
		    data[3*rr+mm]=field[0].dis[numdis].node[rr].sol_increment.a.da[1][mm];
		}
		if (field[0].dis[numdis].node[rr].gnode[0].Knotentyp==Loch_Standard)
		{
		    for (mm=0;mm<3;mm++)
		    {
			data[3*rr+mm]=0;
		    }
		}
	    }
	    if (binaer==0) {
/*---------------------------------------------------------- write ien */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien0","w""b");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien1","w""b");
		}
		fwrite(&ien[0],sizeof(int),4*field[0].dis[numdis].numele,outputfile);
		fclose(outputfile);
/*-------------------------------------------------- write coordinates */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy0","wb");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy1","wb");
		}
		fwrite(&xy[0],sizeof(double),2*field[0].dis[numdis].numnp,outputfile);
		fclose(outputfile);
/*------------------------------------------------------- write values */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_data0","w""b");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_data1","w""b");
		}
		fwrite(&data[0],sizeof(double),3*field[0].dis[numdis].numnp,outputfile);
		fclose(outputfile);
	    }
	    else
	    {
/*---------------------------------------------------------- write ien */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien0asci","w");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_ien1asci","w");
		}
		for (rr=0;rr<field[0].dis[numdis].numele;rr++){
		    for (mm=0;mm<3;mm++){
			fprintf(outputfile," %14d   ",ien[3*rr+mm]);
		    }
		    fprintf(outputfile,"\n");
		}
		fclose(outputfile);
/*-------------------------------------------------- write coordinates */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy0asci","w");
		}
		else
		{
		    outputfile=fopen("./src/chimera/to_pager/chimera_xy1asci","w");
		}
		for (rr=0;rr<field[0].dis[numdis].numnp;rr++){
		    for (mm=0;mm<2;mm++){
			fprintf(outputfile," %14g   ",xy[2*rr+mm]);
		    }
		    fprintf(outputfile,"\n");
		}
		fclose(outputfile);
/*------------------------------------------------------- write values */
		if (numdis == 0){
		    outputfile=fopen("./src/chimera/to_pager/chimera_data0asci","w");
		}
		else
		{
                  outputfile=fopen("./src/chimera/to_pager/chimera_data1asci","w");
		}
		for (rr=0;rr<field[0].dis[numdis].numnp;rr++){
		    for (mm=0;mm<3;mm++){
			fprintf(outputfile," %14g   ",data[3*rr+mm]);
		    }
		    fprintf(outputfile,"\n");
		}
		fclose(outputfile);
	    }
	}
    }
    free(ien);
    free(xy);
    free(data);
#ifdef DEBUG
    dstrc_exit();
#endif

    return;
} /* end of chimera_output */
#endif
