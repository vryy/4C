/*!----------------------------------------------------------------------
\file
\brief contains the routines
 - 's9_surf' which calculates the entries of the load vector if the load 
             is applied on the surface of the shell 
 - 's9_surf_P' surface loads -> Pointloads
 - 's9_surf_onoff' switches the <neum_onoff> to 1 if the differential vector 
                   component is to be loaded  -> Point loads!  

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*! 
\addtogroup SHELL9 
*//*! @{ (documentation module open)*/


/*!----------------------------------------------------------------------
\brief calculates loadvector if load is applied on the shell surface 
                  -> Line and Surface Loads!                                        

<pre>                                                              sh 012/02
This routine calcuates the loadvector if the load is applied on the surface
of a multilayerd shell element (shell9)  -> Line and Surface Loads!
</pre>
\param **eload   DOUBLE  (i/o) the modified loadvector of the element
\param   nsurf   INT     (i)   1=MID; 2=TOP; 3=BOT
\param   numklay INT     (i)   number of kinematic layers to this element
\param   iel     INT     (i)   number of nodes to this element

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: s9eleload(), s9loadGP()

*-----------------------------------------------------------------------*/
void s9_surf(DOUBLE **eload, INT nsurf, INT numklay, INT iel)
{
INT     i,j,l;
INT     midlay;   /* middle layer if uneven number of numklay */
INT     fact;     /* factor depending on the norm of a3 */
INT     mod;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_surf");
#endif
/*--------- factor due to norm of a3 -----------------------------------*/
/*fact = 1; /* |a3_L| = h_L */
/*fact = 2; /* |a3_L| = 0.5 * h_L    -> A3_IST_EINHALB*/
fact = (1.0/A3FAC_SHELL9);
/*---------- NSURF = MID -----------------------------------------------*/
if (nsurf == 1) goto end;
else /*NSURF != MID */
{
   mod = numklay % 2;    /*divide modulo 2*/
}
/*---------- NSURF = TOP -----------------------------------------------*/
if (nsurf == 2)
{
  if (mod == 0)      /* even number of kinematic layers */
  {
    for (l=(numklay/2)+1; l<=numklay; l++) /*upper half of kinlay*/
    {
      for (i=0; i<iel; i++)
      {
         for (j=0; j<3; j++)
         {
            eload[j+3*l][i] = eload[j][i] * fact;
         }
      }/*nodes*/
    }/*layers*/
  }/*end even*/
  else                   /* uneven number of kinematic layers*/
  {
    midlay = (numklay+1)/2;
    for (i=0; i<iel; i++)
    {
       for (j=0; j<3; j++) /*only half of the load on the middle layer*/
       {
          eload[j+3*midlay][i] = 0.5*eload[j][i] * fact;
       }
    }
    for (l=midlay+1; l<=numklay; l++) /*upper half of kinlay*/
    {
      for (i=0; i<iel; i++)
      {
         for (j=0; j<3; j++)
         {
            eload[j+3*l][i] = eload[j][i] * fact;
         }
      }/*nodes*/
    }/*layers*/
  }/*end uneven*/
}
/*---------- NSURF = BOT -----------------------------------------------*/
if (nsurf == 3)
{
  if (mod == 0)      /* even number of kinematic layers */
  {
    for (l=1; l<=numklay/2; l++) /*lower half of kinlay*/
    {
      for (i=0; i<iel; i++)
      {
         for (j=0; j<3; j++)
         {
            eload[j+3*l][i] = -eload[j][i] * fact;
         }
      }/*nodes*/
    }/*layers*/
  }/*end even*/
  else                   /* uneven number of kinematic layers*/
  {
    midlay = (numklay+1)/2;
    for (i=0; i<iel; i++)
    {
       for (j=0; j<3; j++) /*only half of the load on the middle layer*/
       {
          eload[j+3*midlay][i] = -0.5*eload[j][i] * fact;
       }
    }
    for (l=1; l<midlay; l++) /*lower half of kinlay*/
    {
      for (i=0; i<iel; i++)
      {
         for (j=0; j<3; j++)
         {
            eload[j+3*l][i] = -eload[j][i] * fact;
         }
      }/*nodes*/
    }/*layers*/
  }/*end uneven*/
}
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9_surf */




/*!----------------------------------------------------------------------
\brief calculates loadvector if load is applied on the shell surface 
                  -> Point loads!                                         

<pre>                                                              sh 012/02
This routine calcuates the loadvector if the load is applied on the surface
of a multilayerd shell element (shell9)  -> Point loads!
</pre>
\param  *pload   DOUBLE  (i/o) the modified loadvector of the node
\param   nsurf   INT     (i)   1=MID; 2=TOP; 3=BOT
\param   numklay INT     (i)   number of kinematic layers to this node

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: rhs_point_neum()

*-----------------------------------------------------------------------*/
void s9_surf_P(DOUBLE *pload, INT nsurf, INT numklay)
{
INT     i,j,l;
INT     midlay;   /*middle layer if uneven number of numklay */
INT     fact;     /* factor depending on the norm of a3 */
INT     mod;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_surf_P");
#endif
/*--------- factor due to norm of a3 -----------------------------------*/
/*fact = 1; /* |a3_L| = h_L */
/*fact = 2; /* |a3_L| = 0.5 * h_L    -> A3_IST_EINHALB*/
fact = (1.0/A3FAC_SHELL9);
/*---------- NSURF = MID -----------------------------------------------*/
if (nsurf == 1) goto end;
else /*NSURF != MID */
{
   mod = numklay % 2;    /*divide modulo 2*/
}
/*---------- NSURF = TOP -----------------------------------------------*/
if (nsurf == 2)
{
  if (mod == 0)      /* even number of kinematic layers */
  {
    for (l=(numklay/2)+1; l<=numklay; l++) /*upper half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = pload[j] * fact;
      }
    }/*layers*/
  }/*end even*/
  else                   /* uneven number of kinematic layers*/
  {
    midlay = (numklay+1)/2;
    for (j=0; j<3; j++) /*only half of the load on the middle layer*/
    {
       pload[j+3*midlay] = 0.5*pload[j] * fact;
    }
    for (l=midlay+1; l<=numklay; l++) /*upper half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = pload[j] * fact;
      }
    }/*layers*/
  }/*end uneven*/
}
/*---------- NSURF = BOT -----------------------------------------------*/
if (nsurf == 3)
{
  if (mod == 0)      /* even number of kinematic layers */
  {
    for (l=1; l<=numklay/2; l++) /*lower half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = -pload[j] * fact;
      }
    }/*layers*/
  }/*end even*/
  else                   /* uneven number of kinematic layers*/
  {
    midlay = (numklay+1)/2;
    for (j=0; j<3; j++) /*only half of the load on the middle layer*/
    {
       pload[j+3*midlay] = -0.5*pload[j] * fact;
    }
    for (l=1; l<midlay; l++) /*lower half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = -pload[j] * fact;
      }
    }/*layers*/
  }/*end uneven*/
}
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9_surf_P */


/*!----------------------------------------------------------------------
\brief switches the <neum_onoff> to 1 if the differential vector component
       is to be loaded  -> Point loads!  
       (is similar to s9_surf_P but the array to be modified is a integer)                                         

<pre>                                                              sh 012/02
This routine switches the <neum_onoff> to 1 if the differential vector 
component is to be loaded -> shell element (shell9)  -> Point loads!
</pre>
\param  *pload   INT   (i/o) the modified loadvector of the node
\param   nsurf   INT     (i)   1=MID; 2=TOP; 3=BOT
\param   numklay INT     (i)   number of kinematic layers to this node

\warning There is nothing special to this routine
\return void                                               
\sa calling: ---; called by: rhs_point_neum()

*-----------------------------------------------------------------------*/
void s9_surf_onoff(INT *pload, INT nsurf, INT numklay)
{
INT     i,j,l;
INT     midlay;   /*middle layer if uneven number of numklay */
INT     mod;
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_enter("s9_surf_onoff");
#endif
/*---------- NSURF = MID -----------------------------------------------*/
if (nsurf == 1) goto end;
else /*NSURF != MID */
{
    mod = numklay % 2;    /*divide modulo 2*/
}
/*---------- NSURF = TOP -----------------------------------------------*/
if (nsurf == 2)
{
  if (mod == 0)      /* even number of kinematic layers */
  {
    for (l=(numklay/2)+1; l<=numklay; l++) /*upper half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = pload[j];
      }
    }/*layers*/
  }/*end even*/
  else                   /* uneven number of kinematic layers*/
  {
    midlay = (numklay+1)/2;
    for (l=midlay; l<=numklay; l++) /*upper half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = pload[j];
      }
    }/*layers*/
  }/*end uneven*/
}
/*---------- NSURF = BOT -----------------------------------------------*/
if (nsurf == 3)
{
  if (mod == 0)      /* even number of kinematic layers */
  {
    for (l=1; l<=numklay/2; l++) /*lower half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = pload[j];
      }
    }/*layers*/
  }/*end even*/
  else                   /* uneven number of kinematic layers*/
  {
    midlay = (numklay+1)/2;
    for (l=1; l<=midlay; l++) /*lower half of kinlay*/
    {
      for (j=0; j<3; j++)
      {
         pload[j+3*l] = pload[j];
      }
    }/*layers*/
  }/*end uneven*/
}
/*----------------------------------------------------------------------*/
end:
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return; 
} /* end of s9_surf_onoff */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/
