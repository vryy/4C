/*!----------------------------------------------------------------------
\file
\brief service routines for multi-level fluid

<pre>
Maintainer:  name
             e-mail
             homepage
             telephone number


</pre>

------------------------------------------------------------------------*/
#ifdef D_FLUID
#include "../headers/standardtypes.h"
#include "fluid_prototypes.h"

/*!---------------------------------------------------------------------                                         
\brief initialization of (sub-)submesh for multi-level fluid

<pre>                                                       gravem 07/03

In this routine, the necessary variables and arrays for the creation of
submeshes and sub-submeshes are initialized (only quadrilaterals and
hexaeders so far).
			     
</pre>   
\param *actfield   FIELD           (i)  actual field  
\param *fdyn	   FLUID_DYNAMIC   (i)  
\return void 

------------------------------------------------------------------------*/
void fluid_ml_init(FIELD         *actfield,  
                   FLUID_DYNAMIC *fdyn)      
{
INT              numsd;        /* number of spatial dimensions          */
INT              numen;        /* number of element nodes               */
INT              ndum;         /* dummy variable                        */
INT              xele,yele,zele; /* number of elements in coord. dir.   */
INT              ngpr,ngps,ngpt; /*  number of gauss p. in coord. dir.  */
INT              smorder;      /* interpolation order for submesh       */
INT              ssmorder;     /* interpolation order for sub-submesh   */
FLUID_DYN_ML    *mlvar;        /* pointer to fluid_dyn_ml               */
FLUID_ML_SMESH  *submesh;      /* pointer to fluid_ml_smesh             */
FLUID_ML_SMESH  *ssmesh;       /* pointer to fluid_ml_smesh             */

#ifdef DEBUG 
dstrc_enter("fluid_ml_init");
#endif

/*---------------------------------------------------- set some values */
numsd   = fdyn->numdf-1;
numen   = actfield->dis[0].element[0].numnp;
mlvar   = &(fdyn->mlvar);
submesh = &(mlvar->submesh);
if (mlvar->smsgvi>2) ssmesh = &(mlvar->ssmesh);

/*--------- number of velocity, pressure and overall bubble functions */
mlvar->nvbub  = numen; 
mlvar->npbub  = numsd*numen;
mlvar->nelbub = mlvar->nvbub + mlvar->npbub + numsd;
   
/*----- submesh values (only linear/quadratic rectangle/hexaeder yet) */
math_intextract(mlvar->smelenum,&ndum,&xele,&yele,&zele);
math_intextract(mlvar->smnumgp,&ndum,&ngpr,&ngps,&ngpt);
smorder = mlvar->smorder;
submesh->ntyp = 1;
if (smorder>1) 
{
  if (numsd>2) submesh->typ=hex20;
  else submesh->typ=quad9;
}    
else 
{
  if (numsd>2) submesh->typ=hex8;
  else submesh->typ=quad4;
}  
/*---------------------------------------------- basic values for 2-D */
submesh->ngpr   = ngpr;   
submesh->ngps   = ngps;   
submesh->ngpt   = ngpt;   
submesh->numen  = (smorder+1)*(smorder+1); 
submesh->numele = xele*yele; 
submesh->numnp  = (smorder*xele+1)*(smorder*yele+1); 
submesh->numeq  = (smorder*xele-1)*(smorder*yele-1); 
/*--------------------------------------------- modifications for 3-D */
if (numsd>2)
{
  submesh->numen  = submesh->numen*(smorder+1); 
  submesh->numele = submesh->numele*zele;
  submesh->numnp  = submesh->numnp*(smorder*zele+1); 
  submesh->numeq  = submesh->numeq*(smorder*zele-1); 
}   
   
/*---------------------------------------------------- submesh arrays */
amdef("smxyzpd",&(submesh->xyzpd),numsd,submesh->numnp,"DA");
amzero(&(submesh->xyzpd));

amdef("smid",&(submesh->id),submesh->numnp,1,"IV");
amzero(&(submesh->id));

amdef("smien",&(submesh->ien),submesh->numele,submesh->numen,"IA");
amzero(&(submesh->ien));

amdef("smmat",&(submesh->mat),submesh->numeq,submesh->numeq,"DA");
amzero(&(submesh->mat));

amdef("smrhs",&(submesh->rhs),submesh->numeq,mlvar->nelbub,"DA");
amzero(&(submesh->rhs));

amdef("smipiv",&(submesh->ipiv),submesh->numeq,1,"IV");
amzero(&(submesh->ipiv));

if (mlvar->smsgvi>2)
{
/*- sub-submesh values (only linear/quadratic rectangle/hexaeder yet) */
  math_intextract(mlvar->ssmelenum,&ndum,&xele,&yele,&zele);
  math_intextract(mlvar->ssmnumgp,&ndum,&ngpr,&ngps,&ngpt);
  ssmorder = mlvar->ssmorder; 
  ssmesh->ntyp = 1;
  if (smorder>1) 
  {
    if (numsd>2) ssmesh->typ=hex20;
    else ssmesh->typ=quad9;
  }    
  else 
  {
    if (numsd>2) ssmesh->typ=hex8;
    else ssmesh->typ=quad4;
  }  
/*---------------------------------------------- basic values for 2-D */
  ssmesh->ngpr   = ngpr;   
  ssmesh->ngps   = ngps;   
  ssmesh->ngpt   = ngpt;   
  ssmesh->numen  = (ssmorder+1)*(ssmorder+1); 
  ssmesh->numele = xele*yele; 
  ssmesh->numnp  = (ssmorder*xele+1)*(ssmorder*yele+1); 
  ssmesh->numeq  = (ssmorder*xele-1)*(ssmorder*yele-1); 
  if (numsd>2)
  {
/*--------------------------------------------- modifications for 3-D */
    ssmesh->numen  = ssmesh->numen*(ssmorder+1); 
    ssmesh->numele = ssmesh->numele*zele;
    ssmesh->numnp  = ssmesh->numnp*(ssmorder*zele+1); 
    ssmesh->numeq  = ssmesh->numeq*(ssmorder*zele-1); 
  }   
   
/*------------------------------------------------ sub-submesh arrays */
  amdef("ssxyzpd",&(ssmesh->xyzpd),numsd,ssmesh->numnp,"DA");
  amzero(&(ssmesh->xyzpd));

  amdef("ssid",&(ssmesh->id),ssmesh->numnp,1,"IV");
  amzero(&(ssmesh->id));

  amdef("ssien",&(ssmesh->ien),ssmesh->numele,ssmesh->numen,"IA");
  amzero(&(ssmesh->ien));

  amdef("ssmat",&(ssmesh->mat),ssmesh->numeq,ssmesh->numeq,"DA");
  amzero(&(ssmesh->mat));

  amdef("ssrhs",&(ssmesh->rhs),ssmesh->numeq,1,"DV");
  amzero(&(ssmesh->rhs));

  amdef("ssipiv",&(ssmesh->ipiv),ssmesh->numeq,1,"IV");
  amzero(&(ssmesh->ipiv));
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_ml_init */ 
   
/*!---------------------------------------------------------------------                                         
\brief copy submesh solution at (n+1) to (n) for multi-level fluid

<pre>                                                       gravem 07/03

In this routine, the submesh solution obtained at time step (n+1) is
copied to the place (n) in the solution history.
			     
</pre>   
\param *actpart    PARTITION       (i)  
\param *fdyn	   FLUID_DYNAMIC   (i)  
\return void 

------------------------------------------------------------------------*/
void fluid_smcopy(PARTITION       *actpart,      
                  FLUID_DYNAMIC   *fdyn)      
{
INT           i;
INT           numrhs,numeq;
ELEMENT      *actele;
FLUID_DYN_ML *mlvar;            

#ifdef DEBUG 
dstrc_enter("fluid_smcopy");
#endif

/*---------------------------- check for quasi-static bubble assumption */
mlvar = &(fdyn->mlvar);
if (mlvar->quastabub!=0) goto end;
/*----------------------------- set number of equations and rhs vectors */
numrhs  = mlvar->nelbub;
numeq   = mlvar->submesh.numeq;

for (i=0; i<actpart->pdis[0].numele; i++)
{
/*--------------------------------------- set pointer to active element */
   actele = actpart->pdis[0].element[i];
/*--------------------------------- go to specific element copy routine */
   switch(actele->eltyp)
   {
     case el_fluid2: 
       f2_smcopy2(actele,numeq,numrhs);
     break;
     case el_fluid3: 
       f3_smcopy2(actele,numeq,numrhs);
     break; 
     default:
      dserror("Type of element unknown");
   }
}

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_smcopy */ 
   
/*!---------------------------------------------------------------------                                         
\brief print out global submesh matrix and rhs vectors

<pre>                                                       gravem 07/03

In this routine, the global submesh matrix and rhs vectors are written 
in the file 'gmr.dat' for control purposes.
			     
</pre>   
\param  **smmat      DOUBLE  (i)    submesh global matrix  
\param  **semat      DOUBLE  (i)    submesh global rhs vectors  
\param    numeq      INT     (i)    number of submesh equations  
\param    numrhs     INT     (i)    number of submesh rhs vectors  
\return void 

------------------------------------------------------------------------*/
void fluid_prgmr(DOUBLE    **smmat,
                 DOUBLE    **smrhs,
		 INT         numeq,
                 INT         numrhs)      
{
INT           i,j;
FILE         *matrhs;

#ifdef DEBUG 
dstrc_enter("fluid_prgmr");
#endif

/*----------------------------------------------------------- open file */
matrhs = fopen("gmr.dat","w");

/*--- write numbering in header: first rhs vectors, then matrix columns */
for (i=0;i<numrhs;i++)
{
  fprintf(matrhs,"%12d",i);
}  
for (i=0;i<numeq;i++)
{
  fprintf(matrhs,"%12d",i);
}  

fprintf(matrhs,"\n");

/*----------- write out values: first rhs vectors, then matrix columns */
for (i=0;i<numeq;i++)
{
  for (j=0;j<numrhs;j++)
  {
    fprintf(matrhs,"%12.3e",smrhs[i][j]);
  }  
  for (j=0;j<numeq;j++)
  {
    fprintf(matrhs,"%12.3e",smmat[i][j]);
  }
  fprintf(matrhs,"\n");
}    

/*---------------------------------------------------------- close file */
fclose(matrhs);

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_prgmr */ 

/*!---------------------------------------------------------------------                                         
\brief assembling routine for (sub-)submesh matrices

<pre>                                                       gravem 07/03

In this routine, the addition of (sub-)submesh element matrices to the 
global (sub-)submesh matrix is performed
			     
</pre>   
\param  **smat       DOUBLE  (i/o)  (sub-)submesh global matrix  
\param  **semat      DOUBLE  (i)    (sub-)submesh element matrix  
\param    numen      INT     (i)    number of (sub-)submesh element nodes  
\param   *slme       INT     (i)    (sub-)submesh element location vector   
\param    fac        DOUBLE  (i)    factor   
\return void 

------------------------------------------------------------------------*/
void fluid_add_smat(DOUBLE  **smat,
		    DOUBLE  **semat,
		    INT       numen,
		    INT      *slme,  
                    DOUBLE    fac)      
{
INT  icol,lcc,ic,irow,lcr,jr;

#ifdef DEBUG 
dstrc_enter("fluid_add_smat");
#endif

/*--------------------------------------------- loop over element nodes */
icol=-1;
for (lcc=0;lcc<numen;lcc++)
{
  icol++;
  ic=slme[lcc];
  if (ic==-1) continue;
  irow=-1;
/*--------------------------------------------- loop over element nodes */
  for (lcr=0;lcr<numen;lcr++)
  {
    irow++;
    jr=slme[lcr];
    if (jr==-1) continue;
    smat[jr][ic] += semat[irow][icol]*fac;
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_add_smat */ 
   
/*!---------------------------------------------------------------------                                         
\brief assembling routine for submesh rhs vectors

<pre>                                                       gravem 07/03

In this routine, the addition of submesh element rhs vectors to the 
global submesh rhs vectors is performed. There are (numrhs) rhs vectors.
			     
</pre>   
\param  **smrhs      DOUBLE  (i/o)  submesh global rhs vectors  
\param  **smerhs     DOUBLE  (i)    submesh element rhs vectors  
\param    numrhs     INT     (i)    number of rhs vectors  
\param    numen      INT     (i)    number of submesh element nodes  
\param   *smlme      INT     (i)    submesh element location vector   
\return void 

------------------------------------------------------------------------*/
void fluid_add_smrhs(DOUBLE  **smrhs,
		     DOUBLE  **smerhs,
                     INT       numrhs,
		     INT       numen,
		     INT      *smlme)	  
{
INT  irhs,irow,inod,nn;

#ifdef DEBUG 
dstrc_enter("fluid_add_smrhs");
#endif

/*--------------------------------------------- loop over number of rhs */
for (irhs=0;irhs<numrhs;irhs++)
{
  irow=-1;
/*--------------------------------------------- loop over element nodes */
  for (inod=0;inod<numen;inod++)
  {
    irow++;
    nn=smlme[inod];
    if (nn==-1) continue;
    smrhs[nn][irhs] += smerhs[irow][irhs];
  }
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_add_smrhs */ 
   
/*!---------------------------------------------------------------------                                         
\brief assembling routine for sub-submesh rhs vector

<pre>                                                       gravem 07/03

In this routine, the addition of sub-submesh element rhs vectors to the 
global sub-submesh rhs vector is performed.
			     
</pre>   
\param  **ssrhs      DOUBLE  (i/o)  sub-submesh global rhs vector  
\param  **sserhs     DOUBLE  (i)    sub-submesh element rhs vector  
\param    numen      INT     (i)    number of sub-submesh element nodes  
\param   *sslme      INT     (i)    sub-submesh element location vector   
\return void 

------------------------------------------------------------------------*/
void fluid_add_ssrhs(DOUBLE   *ssrhs,
		     DOUBLE   *sserhs,
		     INT       numen,
		     INT      *sslme)	  
{
INT  inod,nn;

#ifdef DEBUG 
dstrc_enter("fluid_add_ssrhs");
#endif

/*--------------------------------------------- loop over element nodes */
for (inod=0;inod<numen;inod++)
{
  nn=sslme[inod];
  if (nn==-1) continue;
  ssrhs[nn] += sserhs[inod];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_add_ssrhs */ 
   
/*!---------------------------------------------------------------------                                         
\brief integration routine for submesh lhs (only diffusive part)

<pre>                                                       gravem 07/03

In this routine, the addition of the submesh lhs element integral (only 
diffusive part) to the global submesh lhs integral is performed.
			     
</pre>   
\param   *smidiff    DOUBLE  (i/o)  submesh global lhs integral 
\param  **smiediff   DOUBLE  (i)    submesh element lhs integral  
\param    numen      INT     (i)    number of submesh element nodes  
\param   *smlme      INT     (i)    submesh element location vector   
\return void 

------------------------------------------------------------------------*/
void fluid_add_intlhs(DOUBLE   *smidiff,
		      DOUBLE  **smiediff,
		      INT       numen,
		      INT      *smlme)      
{
INT  inod,nn;

#ifdef DEBUG 
dstrc_enter("fluid_add_intlhs");
#endif

/*--------------------------------------------- loop over element nodes */
for (inod=0;inod<numen;inod++)
{
  nn=smlme[inod];
  if (nn==-1) continue;
  *smidiff += smiediff[inod][inod];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_add_intlhs */ 
   
/*!---------------------------------------------------------------------                                         
\brief integration routine for submesh rhs 

<pre>                                                       gravem 07/03

In this routine, the addition of the submesh lhs element integral (only 
diffusive part) to the global submesh lhs integral is performed.
			     
</pre>   
\param   *smirhs     DOUBLE  (i/o)  submesh global rhs integral 
\param  **smierhs    DOUBLE  (i)    submesh element rhs integral  
\param    numen      INT     (i)    number of submesh element nodes  
\param   *smlme      INT     (i)    submesh element location vector   
\return void 

------------------------------------------------------------------------*/
void fluid_add_intrhs(DOUBLE   *smirhs,
		      DOUBLE   *smierhs,
		      INT       numen,
		      INT      *smlme)      
{
INT  inod,nn;

#ifdef DEBUG 
dstrc_enter("fluid_add_intrhs");
#endif

/*--------------------------------------------- loop over element nodes */
for (inod=0;inod<numen;inod++)
{
  nn=smlme[inod];
  if (nn==-1) continue;
  *smirhs += smierhs[inod];
}

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_add_intrhs */ 
   
/*!---------------------------------------------------------------------                                         
\brief evaluate bubble functions at integration point of submesh element 

<pre>                                                       gravem 07/03

In this routine, the values of the bubble functions at the actual 
integration point of the submesh element are determined.
			     
</pre>   
\param   *bubint     DOUBLE  (o)    bubble functions at integration p. 
\param   *smfunct    DOUBLE  (i)    submesh element shape functions  
\param  **ebub       DOUBLE  (i)    bubble functions at sm element nodes  
\param    smiel      INT     (i)    number of submesh element nodes   
\param    nbub       INT     (i)    number of bubble functions   
\return void 

------------------------------------------------------------------------*/
void fluid_bubint(DOUBLE  *bubint,     
                  DOUBLE  *smfunct,    
	          DOUBLE **ebub,
		  INT      smiel,     
	          INT      nbub) 
{
INT     i,j;

#ifdef DEBUG 
dstrc_enter("fluid_bubint");
#endif

for (i=0;i<nbub;i++) /* loop over number of bubble functions */
{
   bubint[i]=ZERO; 
   for (j=0;j<smiel;j++) /* loop over all nodes of the submesh element */
   {
      bubint[i] += smfunct[j]*ebub[j][i];
   } /* end loop over j */
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of fluid_bubint */

/*!---------------------------------------------------------------------                                         
\brief evaluate pressure bubble functions at int. p. of sm element 

<pre>                                                       gravem 07/03

In this routine, the values of the pressure bubble functions at the 
actual integration point of the submesh element are determined.
			     
</pre>   
\param   *pbubint    DOUBLE  (o)    pressure bubble funct. at int. p. 
\param   *smfunct    DOUBLE  (i)    submesh element shape functions  
\param  **epbub      DOUBLE  (i)    press. bubble funct. at sm ele. nod.  
\param    smiel      INT     (i)    number of submesh element nodes   
\param    iel        INT     (i)    number of large-scale element nodes   
\param    nsd        INT     (i)    number of spatial dimensions   
\return void 

------------------------------------------------------------------------*/
void fluid_pbubint(DOUBLE **pbubint,     
                   DOUBLE  *smfunct,    
	           DOUBLE **epbub,
		   INT      smiel,     
		   INT      iel,     
		   INT      nsd) 
{
INT     i,j,k,ibub;

#ifdef DEBUG 
dstrc_enter("fluid_pbubint");
#endif

ibub=0;
for (i=0;i<iel;i++) /* loop over number of large-scale element nodes */
{
  for (j=0;j<nsd;j++) /* loop over number of spatial dimensions */
  {
    pbubint[j][i]=ZERO; 
    for (k=0;k<smiel;k++) /* loop over all nodes of the submesh element */
    {
      pbubint[j][i] += smfunct[k]*epbub[k][ibub];
    } /* end loop over j */
    ibub++;
  }
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of fluid_pbubint */

/*!---------------------------------------------------------------------                                         
\brief evaluate bubble function derivatives at int. p. of sm element 

<pre>                                                       gravem 07/03

In this routine, the values of the bubble function derivatives at the 
actual integration point of the submesh element are determined.
			     
</pre>   
\param  **bubderxy   DOUBLE  (o)    bubble funct. derivatives at int. p. 
\param  **smderxy    DOUBLE  (i)    sm element shape funct. derivatives  
\param  **ebub       DOUBLE  (i)    bubble functions at sm element nodes  
\param    smiel      INT     (i)    number of submesh element nodes   
\param    nbub       INT     (i)    number of bubble functions   
\param    nder       INT     (i)    number of derivatives   
\return void 

------------------------------------------------------------------------*/
void fluid_bubder(DOUBLE **bubderxy,     
                  DOUBLE **smderxy,    
	          DOUBLE **ebub,
		  INT      smiel,     
	          INT      nbub,
		  INT      nder) 
{
INT     i,j,k;

#ifdef DEBUG 
dstrc_enter("fluid_bubder");
#endif

for (i=0;i<nbub;i++) /* loop over number of bubble functions */
{
  for (j=0;j<nder;j++) /* loop over number of derivatives */
  {
    bubderxy[j][i]=ZERO; 
    for (k=0;k<smiel;k++) /* loop over all nodes of the submesh element */
    {
      bubderxy[j][i] += smderxy[j][k]*ebub[k][i];
    } /* end loop over j */
  }  
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of fluid_bubder */

/*!---------------------------------------------------------------------                                         
\brief evaluate pressure bubble funct. deriv. at int. p. of sm element 

<pre>                                                       gravem 07/03

In this routine, the values of the pressure bubble function derivatives 
at the actual integration point of the submesh element are determined.
			     
</pre>   
\param  **pbubderxy  DOUBLE  (o)    press. bub. funct. deriv. at int. p. 
\param  **smderxy    DOUBLE  (i)    sm element shape funct. derivatives  
\param  **epbub      DOUBLE  (i)    press. bubble funct. at sm ele. nod.  
\param    smiel      INT     (i)    number of submesh element nodes   
\param    iel        INT     (i)    number of large-scale element nodes   
\param    nsd        INT     (i)    number of spatial dimensions   
\param    nder       INT     (i)    number of derivatives   
\return void 

------------------------------------------------------------------------*/
void fluid_pbubder(DOUBLE ***pbubderxy,     
                   DOUBLE  **smderxy,    
	           DOUBLE  **epbub,
		   INT       smiel,     
		   INT       iel,     
	           INT       nsd,
		   INT       nder) 
{
INT     i,j,k,l,ibub;

#ifdef DEBUG 
dstrc_enter("fluid_pbubder");
#endif

ibub=0;
for (i=0;i<iel;i++) /* loop over number of large-scale element nodes */
{
  for (j=0;j<nsd;j++) /* loop over number of spatial dimensions */
  {
    for (k=0;k<nder;k++) /* loop over number of derivatives */
    {
      pbubderxy[k][j][i]=ZERO; 
      for (l=0;l<smiel;l++) /* loop over all nodes of the sm element */
      {
        pbubderxy[k][j][i] += smderxy[k][l]*epbub[l][ibub];
      } /* end loop over j */
    }  
    ibub++;
  }  
} /* end loop over i */

/*---------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return; 
} /* end of fluid_pbubder */

/*!---------------------------------------------------------------------                                         
\brief evaluate normalized rhs forces for fluid

<pre>                                                       gravem 07/03

In this routine, normalized rhs forces for fluid are calculated.

</pre>
\param **erhs      DOUBLE	   (i/o)  element rhs
\param  *funct     DOUBLE	   (i)    natural shape functions
\param   fac 	   DOUBLE	   (i)    weighting factor	      
\param   iel	   INT  	   (i)	  number of element nodes
\return void                                                                       

------------------------------------------------------------------------*/
void fluid_calnofo(DOUBLE       *erhs,  
		   DOUBLE  	*funct,  
		   DOUBLE  	 fac,	 
		   INT		 iel)
{
INT    irow;

#ifdef DEBUG 
dstrc_enter("fluid_calnofo");
#endif

for (irow=0;irow<iel;irow++)
{
  erhs[irow] += funct[irow]*fac;
} /* end loop over irow */

/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif

return;
} /* end of fluid_calnofo */

void fluid_mlcaldirich(
                       ELEMENT   *actele,  
		       DOUBLE    *dforces, 
                       DOUBLE   **estif,   
		       INT       *hasdirich
		      )     
{

INT         i,j;
INT         dof;
INT         numdf;                      /* number of fluid dofs         */
INT         nd=0;                      
DOUBLE      dirich[MAXDOFPERELE];       /* dirichlet values of act. ele */
INT         dirich_onoff[MAXDOFPERELE]; /* dirichlet flags of act. ele  */ 
GNODE      *actgnode;	                /* actual GNODE                 */
NODE       *actnode;	                /* actual NODE                  */

#ifdef DEBUG 
dstrc_enter("fluid_mlcaldirich");
#endif  

/*------------------------- check if there are any dirichlet conditions *
                                          for the nodes of this element */
for (i=0; i<actele->numnp; i++)
{
   actgnode = actele->node[i]->gnode;   
   if (actgnode->dirich==NULL) 
      continue;
   else
      *hasdirich=1;
      break;
} /* end loop over nodes */					  

if (*hasdirich==0) /* --> no nodes with DBC for this element */
   goto end;

/*---------------------------------- set number of dofs on this element */
for (i=0; i<actele->numnp; i++) nd += actele->node[i]->numdf;

/*---------------------------- init the vectors dirich and dirich_onoff */
for (i=0; i<nd; i++)
{
   dirich[i] = 0.0;
   dirich_onoff[i] = 0;
}

/*-------------------------------- fill vectors dirich and dirich_onoff */
/*                               dirichlet values at (n+1) were already */
/*                           written to the nodes (sol_increment[3][j]) */
for (i=0; i<actele->numnp; i++) /* loop nodes */
{
   numdf    = actele->node[i]->numdf;
   actnode  = actele->node[i];   
   actgnode = actnode->gnode;
   for (j=0; j<numdf; j++) /* loop dofs */
   {
      if (actgnode->dirich==NULL) continue;
      dirich_onoff[i*numdf+j] = actgnode->dirich->dirich_onoff.a.iv[j];
      dirich[i*numdf+j] = actnode->sol_increment.a.da[3][j];
   } /* end loop over dofs */
} /* end loop over nodes */
/*----------------------------------------- loop rows of element matrix */
for (i=0; i<nd; i++)
{
   /*------------------------------------- do nothing for supported row */
   if (dirich_onoff[i]!=0) continue;
   /*---------------------------------- loop columns of unsupported row */
   for (j=0; j<nd; j++)
   {
      /*---------------------------- do nothing for unsupported columns */
      if (dirich_onoff[j]==0) continue;
      dforces[i] -= estif[i][j] * dirich[j];
   }/* loop j over columns */
}/* loop i over rows */

end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of fluid_mlcaldirich*/ 

#endif








