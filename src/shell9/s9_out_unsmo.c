/*!----------------------------------------------------------------------
\file
\brief calculates and writes output data for postprocessing with gid (unsmoothed way)
contains the routines:

The routines in here are necessary to show the stresses in a unsmoothed manner,
therefor the coordinates have to be written for each element, even for each layer
again, so that the stresses can be written to the different layers without
smoothing. -> Some coordinates may be written up to 8 times (4 adjacent Elements
at one node and one upper and one lower layer)

- s9_out_gid_allcoords_unsmo: calculates the coordinates for visualization as volume element and
                              writes them to the file -> .flavia.msh
- s9_out_gid_eletop_unsmo:    writes the element topology of the calculated 3D elements to the
                              file -> .flavia.msh
- s9_out_gid_sol_dis_unsmo:   calculates the displacements at the nodal points of the 3d element an
                              writes them to the file -> .flavia.res
- s9_out_gid_sol_str_unsmo:   writes the physical stresses, extrapolated from GP to nodal points
                              to the flavia.res file


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0711 - 685-6120
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL9
#include "../headers/standardtypes.h"
#include "shell9.h"

/*!
\addtogroup SHELL9
*//*! @{ (documentation module open)*/

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;
/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*!----------------------------------------------------------------------
\brief calcualtes and writes coordinates to the flavia.msh file
       for a unsmoothed visualisation in gid

<pre>                                                              sh 07/03
This routine calculates the coordinates of the nodal keypoints in the
different layers and writes them to the flavia.msh file
</pre>
\param  *out     FILE    (i)   file to be written on (flavia.msh)
\param   type    INT     (i)   element type (4/8/90-noded)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_gid_msh()

*-----------------------------------------------------------------------*/
void s9_out_gid_allcoords_unsmo(FILE *out,INT type)
{
INT           i,j,k,l,mlay,klay,jlay;
DOUBLE        zeta,e3;
DOUBLE        thick;              /*total thicknes of shell element*/
DOUBLE        klayhgt;            /*thickness of aktual layer in per cent of total thickness*/
DOUBLE        hl;                 /*thickness of aktual layer*/
DOUBLE        a3ref_l[3];         /*director a3 of aktual layer in ref. config*/
INT           numele;             /*number of elements*/
INT           numnp;              /*number of nodes to actele*/
ELEMENT      *actele;
NODE         *actnode;
INT           ele_ID;
INT           node_ID;
INT           print_ID;
INT           num_klay,num_mlay;
DOUBLE       *mlayhgt;
DOUBLE        sum_hgt,sum_hgt_mid;
INT           tot_lay,sum_lay,sum_lay_old;
DOUBLE        x[3],x_u[3],x_o[3]; /*coordinats*/
FIELD        *actfield;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_allcoords_unsmo");
#endif
/*---------------------------------------------------- loop over fields */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   numele = actfield->dis[0].numele;

   /*-------------------------------------- loop over all elements */
   for (j=0; j<numele; j++)
   {
      actele = &(actfield->dis[0].element[j]);
      if (actele->eltyp != el_shell9) continue;
      ele_ID = j;

      thick    = actele->e.s9->thick;
      num_klay = actele->e.s9->num_klay;
      numnp    = actele->numnp;

      /*get total number of layers to the element
        NOTE: it is assumed that there are only shell9 elements in the calculation and all of the
              have the same layer structure */
      tot_lay = 0;
      for (k=0; k<num_klay; k++) tot_lay += actele->e.s9->kinlay[k].num_mlay;

      /*----------------------------- loop over all nodes of actele */
      for (k=0; k<numnp; k++)
      {
         actnode = actele->node[k];
         node_ID  = k+1;
         /*get director of node:
           NOTE: in the undeformed geometrie, the directors of the different kinematic
                 layers are equal to the director in the total reference layer*/
         a3ref_l[0]= actele->e.s9->a3ref.a.da[0][k];
         a3ref_l[1]= actele->e.s9->a3ref.a.da[1][k];
         a3ref_l[2]= actele->e.s9->a3ref.a.da[2][k];

         /*switch the elementtype */
         if (type == 4) /* quad4 -> Hex8 */
         {
            /*init sum_lay: sum of all layers, kinematic and material*/
            sum_lay = 0;
            sum_lay_old = 0;

            /*calculate the coordinates of the bottom surface of the first kinematic layer*/
            klay = 0;
            /*initialize the coordinates to values of middle surface*/
            for (l=0; l<3; l++) x[l] = actnode->x[l];
            /*continuity matrix for aktual kin layer at bot*/
            for (jlay=0; jlay<num_klay; jlay++)
            {
               klayhgt = actele->e.s9->klayhgt[jlay];
               hl      = (klayhgt/100)*thick;
              /* hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
               hl      = A3FAC_SHELL9*hl;
               e3 = -1.0; /*bottom*/
               zeta = s9con(e3,num_klay,klay,jlay,1.0);
               /*calculate coordinates*/
               for (l=0; l<3; l++) x[l] = x[l] + zeta*hl*a3ref_l[l];
            }
            /* write coordinates of first kinematic layer at bottom */
            print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID;
            fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                         print_ID, x[0],x[1],x[2]);
            /*save coordinates for interpolation for the material layers*/
            for (l=0; l<3; l++) x_u[l] = x[l];

            /*calculate the coordinates of the top surface of the each kin layer*/
            for (klay=0; klay<num_klay; klay++)
            {
               /*initialize the coordinates to values of middle surface*/
               for (l=0; l<3; l++) x[l] = actnode->x[l];
               /*continuity matrix for aktual kin layer at top*/
               for (jlay=0; jlay<num_klay; jlay++)
               {
                  klayhgt = actele->e.s9->klayhgt[jlay];
                  hl      = (klayhgt/100)*thick;
                  /*hl      = 0.5*hl; */  /*norm of director is half of kinematic layer hight*/
                  hl      = A3FAC_SHELL9*hl;
                  e3 = +1.0; /*top*/
                  zeta = s9con(e3,num_klay,klay,jlay,1.0);
                  /*calculate coordinates*/
                  for (l=0; l<3; l++) x[l] = x[l] + zeta*hl*a3ref_l[l];
               }
               /*calculate the coordinates of the material layers of aktual kinematic layer*/
               num_mlay = actele->e.s9->kinlay[klay].num_mlay;
               mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
               sum_lay += num_mlay;
               /*save coordinates for interpolation for the material layers*/
               for (l=0; l<3; l++) x_o[l] = x[l];
               /*interpolate for the material layers within the actual kinematic layer */
               sum_hgt  = 0.0; /*initialize the total hgt*/
               for (mlay=1; mlay<num_mlay; mlay++) /*only if more than one material layer to this kinematic layer*/
               {
                  sum_hgt += mlayhgt[mlay-1];
                  for (l=0; l<3; l++) x[l] = x_u[l] + (sum_hgt/100.)*(x_o[l]-x_u[l]);
                  /* write the bottom coordinates of this material layer -> 2 times as node belongs to two layers */
                  print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + ((sum_lay_old+mlay)*2 - 1) * numnp;
                  fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                               print_ID, x[0],x[1],x[2]);
                  print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + ((sum_lay_old+mlay)*2 )    * numnp;
                  fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                               print_ID, x[0],x[1],x[2]);
               }

               if (klay == num_klay-1) /*last kinematic layer*/
               {
                 /* write the top coordinates of the last kinematic layer */
                 print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + ((tot_lay-1)*2 + 1) * numnp;
                 fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                              print_ID, x_o[0],x_o[1],x_o[2]);
               }
               else
               {
                 /* write the top coordinates of this kinematic layer -> 2 times as node belongs to two layers */
                 print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + (sum_lay*2-1)*numnp;
                 fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                              print_ID, x_o[0],x_o[1],x_o[2]);
                 print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + (sum_lay*2)*numnp;
                 fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                              print_ID, x_o[0],x_o[1],x_o[2]);


                 for (l=0; l<3; l++) x_u[l] = x_o[l];
                 sum_lay_old = sum_lay;
               }
            }

         }

         else if (type == 8 || type == 9) /* quad8 -> Hex20 / quad9 -> Hex27*/
         {
            /*init sum_lay: sum of all layers, kinematic and material*/
            sum_lay = 0;
            sum_lay_old = 0;

            /*calculate the coordinates of the bottom surface of the first kinematic layer*/
            klay = 0;
            /*initialize the coordinates to values of middle surface*/
            for (l=0; l<3; l++) x[l] = actnode->x[l];
            /*continuity matrix for aktual kin layer at bot*/
            for (jlay=0; jlay<num_klay; jlay++)
            {
               klayhgt = actele->e.s9->klayhgt[jlay];
               hl      = (klayhgt/100)*thick;
               /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
               hl      = A3FAC_SHELL9*hl;
               e3 = -1.0; /*bottom*/
               zeta = s9con(e3,num_klay,klay,jlay,1.0);
               /*calculate coordinates*/
               for (l=0; l<3; l++) x[l] = x[l] + zeta*hl*a3ref_l[l];
            }
            /* write coordinates of first kinematic layer at bottom */
            print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID;
            fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                         print_ID, x[0],x[1],x[2]);
            /*save coordinates for interpolation for the material layers*/
            for (l=0; l<3; l++) x_u[l] = x[l];

            /*calculate the coordinates of the top surface of the each kin layer*/
            for (klay=0; klay<num_klay; klay++)
            {
               /*initialize the coordinates to values of middle surface*/
               for (l=0; l<3; l++) x[l] = actnode->x[l];
               /*continuity matrix for aktual kin layer at top*/
               for (jlay=0; jlay<num_klay; jlay++)
               {
                  klayhgt = actele->e.s9->klayhgt[jlay];
                  hl      = (klayhgt/100)*thick;
                  /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
                  hl      = A3FAC_SHELL9*hl;   /*norm of director is half of kinematic layer hight*/
                  e3 = +1.0; /*top*/
                  zeta = s9con(e3,num_klay,klay,jlay,1.0);
                  /*calculate coordinates*/
                  for (l=0; l<3; l++) x[l] = x[l] + zeta*hl*a3ref_l[l];
               }
               /*calculate the coordinates of the material layers of aktual kinematic layer*/
               num_mlay = actele->e.s9->kinlay[klay].num_mlay;
               mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
               sum_lay += num_mlay;
               /*save coordinates for interpolation for the material layers*/
               for (l=0; l<3; l++) x_o[l] = x[l];
               /*interpolate for the material layers within the actual kinematic layer */
               sum_hgt     = 0.0; /*initialize the total hgt*/
               sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/
               for (mlay=0; mlay<num_mlay; mlay++)
               {
                  if (mlay > 0) /*if more than one material layer to this kinematic layer*/
                  {
                     sum_hgt += mlayhgt[mlay-1];
                     for (l=0; l<3; l++) x[l] = x_u[l] + (sum_hgt/100.)*(x_o[l]-x_u[l]);
                     /* write the bottom coordinates of this material layer -> 2 times as node belongs to two layers */
                     print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + ((sum_lay_old+mlay)*3 - 1) * 9;
                     fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                                  print_ID, x[0],x[1],x[2]);
                     print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + ((sum_lay_old+mlay)*3 )    * 9;
                     fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                                  print_ID, x[0],x[1],x[2]);
                  }

                  /*a midnode has to be calculated */
                  sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
                  for (l=0; l<3; l++) x[l] = x_u[l] + (sum_hgt_mid/100.)*(x_o[l]-x_u[l]);
                  /* write the bottom coordinates of this material layer  */
                  print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + ((sum_lay_old+mlay)*3 + 1)    * 9;
                  fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                               print_ID, x[0],x[1],x[2]);
               }

               if (klay == num_klay-1) /*last kinematic layer*/
               {
                   /* write the top coordinates of this kinematic layer */
                   print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (tot_lay * 3 - 1)  * 9;
                   fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                                print_ID, x_o[0],x_o[1],x_o[2]);
               }
               else
               {
                   /* write the top coordinates of this kinematic layer -> 2 times as node belongs to two layers */
                   print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (sum_lay * 3 - 1) * 9;
                   fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                                print_ID, x_o[0],x_o[1],x_o[2]);
                   print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (sum_lay * 3)     * 9;
                   fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                                print_ID, x_o[0],x_o[1],x_o[2]);

                   for (l=0; l<3; l++) x_u[l] = x_o[l];
                   sum_lay_old = sum_lay;
               }
            }
         }
         else  dserror("wrong element type for gid mesh -> unsmo !!!");
      } /*end loop over all nodes to actele */
   }/*end loop over all elements */
}/* end of (i=0; i<genprob.numfld; i++)*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_allcoords_unsmo */

/*!----------------------------------------------------------------------
\brief writes the element topology to the flavia.msh file
       for a unsmoothed visualisation in gid

<pre>                                                              sh 07/03
This routine writes the topology of the calculated 3-D elements to the
flavia.msh file
</pre>
\param  *out     FILE    (i)   file to be written on (flavia.msh)
\param   type    INT     (i)   element type (4/8/90-noded)
\param  *actele  ELEMENT (i)   actual element

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_gid_msh()

*-----------------------------------------------------------------------*/
void s9_out_gid_eletop_unsmo(FILE *out,INT type,INT j)
{
INT           i,k,lay;
INT           numnp;            /*number of nodal points*/
INT           numele;           /*number of elements*/
INT           num_klay;
INT           tot_lay;          /*sum of all layers: kinematic and material*/
INT           node_ID;          /*nodal ID*/
INT           ele_ID;           /*element ID*/
INT           print_ID;
FIELD        *actfield;
ELEMENT      *actele;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_eletop_unsmo");
#endif
/*---------------------------------------------------- loop over fields */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   numele = actfield->dis[0].numele;

   actele = &(actfield->dis[0].element[j]);
   if (actele->eltyp != el_shell9) continue;
   ele_ID = j;

   num_klay = actele->e.s9->num_klay;
   numnp    = actele->numnp;

   /*get total number of layers to the element
     NOTE: it is assumed that there are only shell9 elements in the calculation and all of the
           have the same layer structure */
   tot_lay = 0;
   for (k=0; k<num_klay; k++) tot_lay += actele->e.s9->kinlay[k].num_mlay;

   /*switch the elementtype */
   if (type == 4) /*quad4 -> 8-noded Hexahedra element */
   {
      /*--------------------------------- loop over all layers */
      for (lay=0; lay<tot_lay; lay++)
      {
          fprintf(out," %6d ",ele_ID+lay*numele+1);
          /*untere Ebene*/
          for (k=0; k<numnp; k++)
          {
             node_ID  = k+1;
             print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + lay * (numnp * 2);
             fprintf(out,"%6d ",print_ID);
          }
          /*obere Ebene*/
          for (k=0; k<numnp; k++)
          {
             node_ID = k+1;
             print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + lay * (numnp * 2) + numnp;
             fprintf(out,"%6d ",print_ID);
          }
          fprintf(out,"\n");
      }
   }

  else if (type == 8 || type == 9) /*quad8 -> 20-noded Hexahedra element / quad9 -> Hex27 */
  {
     /*--------------------------------- loop over all layers */
     for (lay=0; lay<tot_lay; lay++)
     {
         fprintf(out," %6d ",ele_ID+lay*numele+1);
         /*untere Ebene -> Eckknoten*/
         for (k=0; k<4; k++)
         {
            node_ID = k+1;
            print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3) * 9;
            fprintf(out,"%6d ",print_ID);
         }
         /*obere Ebene -> Eckknoten*/
         for (k=0; k<4; k++)
         {
            node_ID = k+1;
            print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+2) * 9;
            fprintf(out,"%6d ",print_ID);
         }
         /*untere Ebene -> Seitenmittelknoten*/
         for (k=4; k<8; k++)
         {
            node_ID = k+1;
            print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3) * 9;
            fprintf(out,"%6d ",print_ID);
         }
         /*mittlere Ebene -> Eckknoten*/
         for (k=0; k<4; k++)
         {
            node_ID = k+1;
            print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+1) * 9;
            fprintf(out,"%6d ",print_ID);
         }
         /*obere Ebene -> Seitenmittelknoten*/
         for (k=4; k<8; k++)
         {
            node_ID = k+1;
            print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+2) * 9;
            fprintf(out,"%6d ",print_ID);
         }

         /* write additional node numbers for quad9 */
         if (type == 9)
         {
            /*untere Ebene -> Mittelknoten*/
            for (k=8; k<9; k++)
            {
               node_ID = k+1;
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3) * 9;
               fprintf(out,"%6d ",print_ID);
            }
            /*mittlere Ebene -> Seitenmittelknoten*/
            for (k=4; k<8; k++)
            {
               node_ID = k+1;
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+1) * 9;
               fprintf(out,"%6d ",print_ID);
            }
            /*obere Ebene -> Mittelknoten*/
            for (k=8; k<9; k++)
            {
               node_ID = k+1;
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+2) * 9;
               fprintf(out,"%6d ",print_ID);
            }
            /*mittlere Ebene -> Mittelknoten*/
            for (k=8; k<9; k++)
            {
               node_ID = k+1;
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+1) * 9;
               fprintf(out,"%6d ",print_ID);
            }
         }

         fprintf(out,"\n");
     }
  }

   else  dserror("wrong element type for gid mesh -> unsmo !!!");
}/* end of (i=0; i<genprob.numfld; i++)*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_eletop_unsmo */

/*!----------------------------------------------------------------------
\brief calculates and writes the nodal displacements to the flavia.res file
       for a unsmoothed visualisation in gid

<pre>                                                              sh 07/03
This routine calculates the nodal displacement from the different vector
component in the different layers and writes them to the flavia.res file
</pre>
\param  *out      FILE    (i)   file to be written on (flavia.res)
\param  *actfield FIELD   (i)   actual field

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_gid_sol()

*-----------------------------------------------------------------------*/
void s9_out_gid_sol_dis_unsmo(FILE *out,FIELD *actfield,INT place)
{
INT           j,k,l,mlay,klay,jlay;
INT           numele;                   /*number of elements*/
INT           numnp;                    /*number of nodal points*/
DOUBLE        thick;
INT           num_klay,num_mlay;
DOUBLE       *mlayhgt;
DOUBLE        sum_hgt,sum_hgt_mid;
INT           tot_lay,sum_lay,sum_lay_old;
INT           node_ID;                  /*nodal ID*/
INT           ele_ID;
INT           print_ID;
DOUBLE        dis[3],dis_u[3],dis_o[3]; /*displacement x,y,z*/
INT           type;
DOUBLE        e3,zeta;
NODE         *actnode;
ELEMENT      *actele;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_sol_dis_unsmo");
#endif
/*----------------------------------------------------------------------*/
fprintf(out,"VALUES\n");

numele = actfield->dis[0].numele;
/*-------------------------------------- loop over all elements */
for (j=0; j<numele; j++)
{
   actele = &(actfield->dis[0].element[j]);
   if (actele->eltyp != el_shell9) continue;
   ele_ID = j;

   thick    = actele->e.s9->thick;
   num_klay = actele->e.s9->num_klay;
   numnp    = actele->numnp;
   type     = numnp;

   /*get total number of layers to the element
     NOTE: it is assumed that there are only shell9 elements in the calculation and all of the
           have the same layer structure */
   tot_lay = 0;
   for (k=0; k<num_klay; k++) tot_lay += actele->e.s9->kinlay[k].num_mlay;

   /*----------------------------- loop over all nodes of actele */
   for (k=0; k<numnp; k++)
   {
      actnode = actele->node[k];
      node_ID  = k+1;

      /*switch the elementtype */
      if (type == 4) /* quad4 -> Hex8 */
      {
         /*init sum_lay: sum of all layers, kinematic and material*/
         sum_lay = 0;
         sum_lay_old = 0;

         /*calculate the displacements of the bottom surface of the first kinematic layer*/
         klay = 0;
         /*initialize the displacements of this layer to zero */
         for(l=0; l<3; l++) dis[l]=0.0;
         /*calculate the relative displacements of this layer */
         for (jlay=0; jlay<num_klay; jlay++)
         {
            e3   = -1.0; /*bottom*/
            zeta = s9con(e3,num_klay,klay,jlay,1.0);
            for (l=0; l<3; l++) dis[l] = dis[l] + zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
         }
         /*add the displacement of reference layer to relative displacement*/
         for (l=0; l<3; l++) dis[l] = dis[l] + actnode->sol.a.da[place][l];

         /* write the displacements of first kinematic layer at bottom */
         print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID;
         fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                      print_ID, dis[0],dis[1],dis[2]);
         /*save coordinates for interpolation for the material layers*/
         for (l=0; l<3; l++) dis_u[l] = dis[l];

         /*calculate the displacements of the top surface of the each kin layer*/
         for (klay=0; klay<num_klay; klay++)
         {
            /*initialize the displacements of this layer to zero */
            for(l=0; l<3; l++) dis[l]=0.0;
            /*calculate the relative displacements of this layer */
            for (jlay=0; jlay<num_klay; jlay++)
            {
               e3   = +1.0; /*top*/
               zeta = s9con(e3,num_klay,klay,jlay,1.0);
               for (l=0; l<3; l++) dis[l] = dis[l] + zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
            }
            /*add the displacement of reference layer to relative displacement*/
            for (l=0; l<3; l++) dis[l] = dis[l] + actnode->sol.a.da[place][l];

            /*calculate the displacements of the material layers of aktual kinematic layer*/
            num_mlay = actele->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
            sum_lay += num_mlay;
            /*save displacements for interpolation for the material layers*/
            for (l=0; l<3; l++) dis_o[l] = dis[l];
            /*interpolate for the material layers within the actual kinematic layer */
            sum_hgt  = 0.0; /*initialize the total hgt*/
            for (mlay=1; mlay<num_mlay; mlay++) /*only if more than one material layer to this kinematic layer*/
            {
               sum_hgt += mlayhgt[mlay-1];
               for (l=0; l<3; l++) dis[l] = dis_u[l] + (sum_hgt/100.)*(dis_o[l]-dis_u[l]);
               /* write the bottom displacement of this material layer -> 2 times as node belongs to two layers  */
               print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + ((sum_lay_old+mlay)*2 - 1) * numnp;
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            print_ID, dis[0],dis[1],dis[2]);
               print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + ((sum_lay_old+mlay)*2 )    * numnp;
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            print_ID, dis[0],dis[1],dis[2]);

            }

            if (klay == num_klay-1) /*last kinematic layer*/
            {
              /* write the top displacement of the last kinematic layer */
              print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + ((tot_lay-1)*2 + 1) * numnp;
              fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                           print_ID, dis_o[0],dis_o[1],dis_o[2]);
            }
            else
            {
              /* write the top displacement of this kinematic layer -> 2 times as node belongs to two layers */
              print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + (sum_lay*2-1)*numnp;
              fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                           print_ID, dis_o[0],dis_o[1],dis_o[2]);
              print_ID = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + (sum_lay*2)*numnp;
              fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                           print_ID, dis_o[0],dis_o[1],dis_o[2]);


              for (l=0; l<3; l++) dis_u[l] = dis_o[l];
              sum_lay_old = sum_lay;
            }
         }
      }


     else if (type == 8 || type == 9) /* quad8 -> Hex20 / quad9 -> Hex27 */
     {
        /*init sum_lay: sum of all layers, kinematic and material*/
        sum_lay = 0;
        sum_lay_old = 0;

        /*calculate the displacements of the bottom surface of the first kinematic layer*/
        klay = 0;
        /*initialize the displacements of this layer to zero */
        for(l=0; l<3; l++) dis[l]=0.0;
        /*calculate the relative displacements of this layer */
        for (jlay=0; jlay<num_klay; jlay++)
        {
           e3   = -1.0; /*bottom*/
           zeta = s9con(e3,num_klay,klay,jlay,1.0);
           for (l=0; l<3; l++) dis[l] = dis[l] + zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
        }
        /*add the displacement of reference layer to relative displacement*/
        for (l=0; l<3; l++) dis[l] = dis[l] + actnode->sol.a.da[place][l];

        /* write the displacements of first kinematic layer at bottom */
        print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID;
        fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                     print_ID, dis[0],dis[1],dis[2]);
        /*save displacements for interpolation for the material layers*/
        for (l=0; l<3; l++) dis_u[l] = dis[l];

        /*calculate the displacements of the top surface of the each kin layer*/
        for (klay=0; klay<num_klay; klay++)
        {
           /*initialize the displacements of this layer to zero */
           for(l=0; l<3; l++) dis[l]=0.0;
           /*calculate the relative displacements of this layer */
           for (jlay=0; jlay<num_klay; jlay++)
           {
              e3   = +1.0; /*top*/
              zeta = s9con(e3,num_klay,klay,jlay,1.0);
              for (l=0; l<3; l++) dis[l] = dis[l] + zeta*actnode->sol.a.da[place][l+3*(jlay+1)];
           }
           /*add the displacement of reference layer to relative displacement*/
           for (l=0; l<3; l++) dis[l] = dis[l] + actnode->sol.a.da[place][l];

           /*calculate the displacements of the material layers of aktual kinematic layer*/
           num_mlay = actele->e.s9->kinlay[klay].num_mlay;
           mlayhgt  = actele->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
           sum_lay += num_mlay;
           /*save displacements for interpolation for the material layers*/
           for (l=0; l<3; l++) dis_o[l] = dis[l];
           /*interpolate for the material layers within the actual kinematic layer */
           sum_hgt     = 0.0; /*initialize the total hgt*/
           sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/
           for (mlay=0; mlay<num_mlay; mlay++)
           {
              if (mlay > 0) /*if more than one material layer to this kinematic layer*/
              {
                 sum_hgt += mlayhgt[mlay-1];
                 for (l=0; l<3; l++) dis[l] = dis_u[l] + (sum_hgt/100.)*(dis_o[l]-dis_u[l]);
                 /* write the bottom displacement of this material layer -> 2 times as node belongs to two layers   */
                 print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + ((sum_lay_old+mlay)*3 - 1) * 9;
                 fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                              print_ID, dis[0],dis[1],dis[2]);
                 print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + ((sum_lay_old+mlay)*3 )    * 9;
                 fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                              print_ID, dis[0],dis[1],dis[2]);
              }

              /*the displacement at midnode has to be calculated*/
              sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
              for (l=0; l<3; l++) dis[l] = dis_u[l] + (sum_hgt_mid/100.)*(dis_o[l]-dis_u[l]);
              /* write the bottom displacement of this material layer  */
              print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + ((sum_lay_old+mlay)*3 + 1)    * 9;
              fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                           print_ID, dis[0],dis[1],dis[2]);
           }

           if (klay == num_klay-1) /*last kinematic layer*/
           {
               /* write the top displacement of this kinematic layer */
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (tot_lay * 3 - 1)  * 9;
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            print_ID, dis_o[0],dis_o[1],dis_o[2]);
           }
           else
           {
               /* write the top displacement of this kinematic layer -> 2 times as node belongs to two layers */
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (sum_lay * 3 - 1) * 9;
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            print_ID, dis_o[0],dis_o[1],dis_o[2]);
               print_ID = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (sum_lay * 3)     * 9;
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            print_ID, dis_o[0],dis_o[1],dis_o[2]);

               for (l=0; l<3; l++) dis_u[l] = dis_o[l];
               sum_lay_old = sum_lay;
           }

         }
      }

      else  dserror("wrong element type for gid out sol (s9) -> unsmo !!!");
   }
}

fprintf(out,"END VALUES\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_sol_dis_unsmo */


/*!----------------------------------------------------------------------
\brief writes the physical stresses, extrapolated from GP to nodal points
       to the flavia.res file
       for a unsmoothed visualisation in gid

<pre>                                                              sh 07/03
This routine writes the physical stresses, extrapolated from GP to nodal
points to the flavia.res file. The stresses are averaged between the
different elements, so be careful with the interpretation!
The only problem is to write in the right order
</pre>
\param  *out      FILE    (i)   file to be written on (flavia.res)
\param  *actfield FIELD   (i)   actual field

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_gid_sol()

*-----------------------------------------------------------------------*/
void s9_out_gid_sol_str_unsmo(FILE *out,FIELD *actfield,INT place)
{
INT           j,k,lay;
INT           numnp;            /*number of nodal points in actfield */
INT           numele;           /*number of elements to actnode -> in plane*/
INT           num_node;         /*number of nodes to an element*/
DOUBLE      **stress;
INT           node_ID;                  /*nodal ID*/
INT           print_ID;
INT           ele_ID;
INT           stress_ID;
INT           stress_ID_o;
INT           stress_ID_u;
INT           num_klay;
INT           tot_lay;
ELEMENT      *actele;


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_sol_str_unsmo");
#endif
/*----------------------------------------------------------------------*/

numele = actfield->dis[0].numele;

/*-------------------------------------- loop over all elements */
for (j=0; j<numele; j++)
{
   actele = &(actfield->dis[0].element[j]);
   if (actele->eltyp != el_shell9) continue;
   ele_ID = j;

   num_klay = actele->e.s9->num_klay;
   numnp    = actele->numnp;
   num_node = numnp;
   stress   = actele->e.s9->stresses.a.d3[place];

   /*get total number of layers to the element
     NOTE: it is assumed that there are only shell9 elements in the calculation and all of the
           have the same layer structure */
   tot_lay = 0;
   for (k=0; k<num_klay; k++) tot_lay += actele->e.s9->kinlay[k].num_mlay;

   /*switch the type */
   if (num_node == 4) /* -> quad4 */
   {
      /*--------------------------------- loop over all layers */
      for (lay=0; lay<tot_lay; lay++)
      {
         /*untere Ebene*/
         for (k=0; k<numnp; k++)
         {
            node_ID   = k+1;
            print_ID  = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + lay * (numnp * 2);
            stress_ID =                                      node_ID + lay * (numnp * 2) - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                 stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
         }
         /*obere Ebene*/
         for (k=0; k<numnp; k++)
         {
            node_ID   = k+1;
            print_ID  = ( ele_ID * (tot_lay * 2) * numnp ) + node_ID + lay * (numnp * 2) + numnp;
            stress_ID =                                      node_ID + lay * (numnp * 2) + numnp - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                 stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
         }
      }/*====================================================== end loop over all layers*/
   }/*======================== end quad4 */

   else if (num_node == 8 || num_node == 9) /*quad8 or quad9*/
   {
      /*--------------------------------- loop over all layers */
      for (lay=0; lay<tot_lay; lay++)
      {
         /*untere Ebene -> Eckknoten*/
         for (k=0; k<4; k++)
         {
            node_ID   = k+1;
            print_ID  = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3) * 9;
            stress_ID =                                  node_ID + (lay*2) * 9 - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                 stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
         }
         /*obere Ebene -> Eckknoten*/
         for (k=0; k<4; k++)
         {
            node_ID   = k+1;
            print_ID  = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+2) * 9;
            stress_ID =                                  node_ID + (lay*2+1) * 9 - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                 stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
         }
         /*untere Ebene -> Seitenmittelknoten*/
         for (k=4; k<8; k++)
         {
            node_ID   = k+1;
            print_ID  = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3) * 9;
            stress_ID =                                  node_ID + (lay*2) * 9 - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                 stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
         }
         /*mittlere Ebene -> Eckknoten   -> mitteln zwischen oben und unten*/
         for (k=0; k<4; k++)
         {
            node_ID     = k+1;
            print_ID    = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+1) * 9;
            stress_ID_u =                                  node_ID + (lay*2)   * 9 - 1;
            stress_ID_o =                                  node_ID + (lay*2+1) * 9 - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,
                        0.5 * (stress[0][stress_ID_u] + stress[0][stress_ID_o]),
                        0.5 * (stress[1][stress_ID_u] + stress[1][stress_ID_o]),
                        0.5 * (stress[2][stress_ID_u] + stress[2][stress_ID_o]),
                        0.5 * (stress[3][stress_ID_u] + stress[3][stress_ID_o]),
                        0.5 * (stress[4][stress_ID_u] + stress[4][stress_ID_o]),
                        0.5 * (stress[5][stress_ID_u] + stress[5][stress_ID_o]));
         }
         /*obere Ebene -> Seitenmittelknoten*/
         for (k=4; k<8; k++)
         {
            node_ID   = k+1;
            print_ID  = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+2) * 9;
            stress_ID =                                  node_ID + (lay*2+1) * 9 - 1;
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                        print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                 stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
         }

         /* write additional stresses for quad9 */
         if (num_node == 9)
         {
            /*untere Ebene -> Mittelknoten*/
            for (k=8; k<9; k++)
            {
               node_ID   = k+1;
               print_ID  = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3) * 9;
               stress_ID =                                  node_ID + (lay*2) * 9 - 1;
               fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                           print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                    stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
            }
            /*mittlere Ebene -> Seitenmittelknoten  -> mitteln zwischen oben und unten */
            for (k=4; k<8; k++)
            {
               node_ID     = k+1;
               print_ID    = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+1) * 9;
               stress_ID_u =                                  node_ID + (lay*2)   * 9 - 1;
               stress_ID_o =                                  node_ID + (lay*2+1) * 9 - 1;
               fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                           print_ID,
                           0.5 * (stress[0][stress_ID_u] + stress[0][stress_ID_o]),
                           0.5 * (stress[1][stress_ID_u] + stress[1][stress_ID_o]),
                           0.5 * (stress[2][stress_ID_u] + stress[2][stress_ID_o]),
                           0.5 * (stress[3][stress_ID_u] + stress[3][stress_ID_o]),
                           0.5 * (stress[4][stress_ID_u] + stress[4][stress_ID_o]),
                           0.5 * (stress[5][stress_ID_u] + stress[5][stress_ID_o]));
            }
            /*obere Ebene -> Mittelknoten*/
            for (k=8; k<9; k++)
            {
               node_ID   = k+1;
               print_ID  = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+2) * 9;
               stress_ID =                                  node_ID + (lay*2+1) * 9 - 1;
               fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                           print_ID,stress[0][stress_ID],stress[1][stress_ID],stress[2][stress_ID],
                                    stress[3][stress_ID],stress[4][stress_ID],stress[5][stress_ID]);
            }
            /*mittlere Ebene -> Mittelknoten  -> mitteln zwischen oben und unten */
            for (k=8; k<9; k++)
            {
               node_ID     = k+1;
               print_ID    = ( ele_ID * (tot_lay * 3) * 9 ) + node_ID + (lay*3+1) * 9;
               stress_ID_u =                                  node_ID + (lay*2)   * 9 - 1;
               stress_ID_o =                                  node_ID + (lay*2+1) * 9 - 1;
               fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                           print_ID,
                           0.5 * (stress[0][stress_ID_u] + stress[0][stress_ID_o]),
                           0.5 * (stress[1][stress_ID_u] + stress[1][stress_ID_o]),
                           0.5 * (stress[2][stress_ID_u] + stress[2][stress_ID_o]),
                           0.5 * (stress[3][stress_ID_u] + stress[3][stress_ID_o]),
                           0.5 * (stress[4][stress_ID_u] + stress[4][stress_ID_o]),
                           0.5 * (stress[5][stress_ID_u] + stress[5][stress_ID_o]));
            }
         }/*================================== end additional values for Quad9 -> Hex27 */
      }/*====================================================== end loop over all layers*/
   } /*==================================================================== end quad8/9 */

   else  dserror("wrong gp type for gid out sol (s9) !!!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_sol_str_smo_unsmo */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

