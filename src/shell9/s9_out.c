/*!----------------------------------------------------------------------
\file
\brief calculates and writes output data for postprocessing with gid
contains the routines:
- s9out_dof:            writes DOFs to the output file -> .out
- s9out_nodal_dis:      writes the displacements of each layer at top,mid,bot to output file -> .out
- s9_out_gid_allcoords: calculates the coordinates for visualization as volume element and
                        writes them to the file -> .flavia.msh
- s9_out_gid_eletop:    writes the element topology of the calculated 3D elements to the
                        file -> .flavia.msh
- s9_out_gid_sol_dis:   calculates the displacements at the nodal points of the 3d element an
                        writes them to the file -> .flavia.res
- s9_out_gid_sol_str:   writes the physical stresses, extrapolated from GP to nodal points
                        to the flavia.res file


<pre>
Maintainer: Stefan Hartmann
            hartmann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/hartmann/
            0771 - 685-6120
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
\brief writes the DOFs in the outputfile

<pre>                                                           sh 012/02
This routine writes the DOFs of a shell9 element to the outputfile
</pre>
\param  *actnode NODE    (i)   actual node for which the DOFs are written
\param   j       INT     (i)   number of actual node

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_general()

*-----------------------------------------------------------------------*/
void s9out_dof(NODE *actnode,INT j)
{
INT        i;
INT        numdf;          /* Number of DOFs to this node */
INT        num_klay;       /* Number of kinematic layers to this node */
FILE      *out = allfiles.out_out;

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9out_dof");
#endif
/*----------------------------------------------------------------------*/
numdf = actnode->numdf;
num_klay = (numdf-3)/3;
/*----------------------------------------------------------------------*/
/* write header line */
if (j==0)
{
  fprintf(out,"                                          X         Y         Z ");
  for (i=1; i<=num_klay; i++) fprintf(out,"     WX(%d)     WY(%d)     WZ(%d)",i,i,i);
  fprintf(out,"\n\n"); /* new line */
}
/* write DOFs */
fprintf(out,"NODE glob_Id %6d loc_Id %6d    %6d    %6d    %6d",
        actnode->Id,actnode->Id_loc,actnode->dof[0],actnode->dof[1],actnode->dof[2]);
for (i=1; i<=num_klay; i++) fprintf(out,"    %6d    %6d    %6d",actnode->dof[0+i*3],actnode->dof[1+i*3],actnode->dof[2+i*3]);
fprintf(out,"\n"); /* new line */

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9out_dof */
/*----------------------------------------------------------------------*/



/*!----------------------------------------------------------------------
\brief writes the nodal displacements to the outputfile

<pre>                                                           sh 012/02
This routine writes the displacements of several layers of a shell9 element
to the outputfile
</pre>
\param  *actnode NODE    (i)   actual node for which the DOFs are written
\param   place   INT     (i)   actual load step

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_sol()

*-----------------------------------------------------------------------*/
void s9out_nodal_dis(NODE *actnode,INT place)
{
INT        i,j,klay,jlay;
DOUBLE     zeta,e3;
INT        numdf;          /* Number of DOFs to this node */
INT        num_klay;       /* Number of kinematic layers to this node */
DOUBLE     dis[3][3];      /* Displacement at top, mid, bot of each layer with x,y,z components */
FILE      *out = allfiles.out_out;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9out_nodal_dis");
#endif
/*----------------------------------------------------------------------*/
numdf = actnode->numdf;
num_klay = (numdf-3)/3;
/*----------------------------------------------------------------------*/
/* write x,y,z displacement of Middle Surface */
fprintf(out,"    NODE         DISPLACEMENT X   DISPLACEMENT Y   DISPLACEMENT Z\n\n");
fprintf(out," %6d           %10.6E     %10.6E     %10.6E\n\n",
        actnode->Id,actnode->sol.a.da[place][0],actnode->sol.a.da[place][1],actnode->sol.a.da[place][2]);

/* write displacements of diff. vector components layerwise */
fprintf(out,"    NODE  LAYER  DIFF.VECTOR WX   DIFF.VECTOR WY   DIFF.VECTOR WZ\n\n");
for (i=1; i<=num_klay; i++) fprintf(out," %6d%6d    %10.6E    %10.6E    %10.6E\n",
                                    actnode->Id,i,actnode->sol.a.da[place][0+3*i],
                                    actnode->sol.a.da[place][1+3*i],actnode->sol.a.da[place][2+3*i]);
fprintf(out,"  DISPLACEMENT ACROSS THE THICKNESS\n");
fprintf(out,"  AT THE TOP, MIDDLE AND BOTTOM OF EACH KINEMATIC LAYER\n\n");
fprintf(out,"    NODE   LAYER   DISPLACEMENT X   DISPLACEMENT Y   DISPLACEMENT Z\n\n");
for (klay=num_klay-1; klay>=0; klay--)
{
   /*initialize the displacements of this layer to zero */
   for (i=0; i<3; i++) for(j=0; j<3; j++) dis[i][j]=0.0;
   /*calculate the relative displacements of this layer */
   for (jlay=0; jlay<num_klay; jlay++)
   {
      for (i=0; i<3; i++) /*bot, mid, top*/
      {
         e3   = i-1.0;
         zeta = s9con(e3,num_klay,klay,jlay,1.0);
         for (j=0; j<3; j++) dis[i][j] = dis[i][j] + zeta*actnode->sol.a.da[place][j+3*(jlay+1)];
      }
   }
   /*add the displacement of reference layer to relative displacement and write this layer*/
   for (i=2; i>=0; i--)
   {
       for (j=0; j<3; j++)
       {
           dis[i][j] = dis[i][j] + actnode->sol.a.da[place][j];
       }
       fprintf(out," %6d%6d     %10.6E     %10.6E     %10.6E\n",
                   actnode->Id,klay+1,dis[i][0],dis[i][1],dis[i][2]);
   }
}
fprintf(out,"\n"); /*Leerzeile*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9out_nodal_dis */

/*!----------------------------------------------------------------------
\brief calcualtes and writes coordinates to visualize a multilayerd shell
       element with Hexahedra elementsto the flavia.msh file

<pre>                                                              sh 012/02
This routine calculates the coordinates of the nodal keypoints in the
different layers and writes them to the flavia.msh file
</pre>
\param  *out     FILE    (i)   file to be written on (flavia.msh)
\param   type    INT     (i)   element type (4/8/90-noded)

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_gid_msh()

*-----------------------------------------------------------------------*/
void s9_out_gid_allcoords(FILE *out,INT type)
{
INT           i,j,k,l,m,mlay,klay,jlay;
DOUBLE        zeta,e3;
DOUBLE        thick;              /*total thicknes of shell element*/
DOUBLE        klayhgt;            /*thickness of aktual layer in per cent of total thickness*/
DOUBLE        hl;                 /*thickness of aktual layer*/
DOUBLE        a3ref_l[3];         /*director a3 of aktual layer in ref. config*/
INT           numnp;              /*number of nodal points*/
INT           num_klay,num_mlay;
DOUBLE       *mlayhgt;
DOUBLE        sum_hgt,sum_hgt_mid;
INT           sum_lay,sum_lay_old;
DOUBLE        x[3],x_u[3],x_o[3]; /*coordinats*/
INT           node_ID;            /*nodal ID*/
INT           is_edge;            /*is_edge=1: actnode is a edge node to the element, else middle node*/
FIELD        *actfield;
NODE         *actnode;
/*-----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_allcoords");
#endif
/*---------------------------------------------------- loop over fields */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   numnp = actfield->dis[0].numnp;

   /*-------------------------------------- loop over all nodal points */
   for (j=0; j<numnp; j++)
   {
      actnode = &(actfield->dis[0].node[j]);

      thick    = actnode->element[0]->e.s9->thick;
      num_klay = actnode->element[0]->e.s9->num_klay;
      node_ID  = actnode->Id+1;

      /*get the director of this node from the first element to this node*/
      for (l=0; l<actnode->element[0]->numnp; l++)
      {
         if (actnode->element[0]->node[l] == actnode)
         {
            a3ref_l[0]= actnode->element[0]->e.s9->a3ref.a.da[0][l];
            a3ref_l[1]= actnode->element[0]->e.s9->a3ref.a.da[1][l];
            a3ref_l[2]= actnode->element[0]->e.s9->a3ref.a.da[2][l];
            break;
         }
      }


      /*switch the elementtype */
      if (type == 4) /* quad4 -> Hex8 */
      {
         /*init sum_lay: sum of all layers, kinematic and material*/
         sum_lay = 0;
         sum_lay_old = 0;

         /*calculate the coordinates of the bottom surface of the first kinematic layer*/
         klay = 0;
         /*initialize the coordinates to values of middle surface*/
         for (k=0; k<3; k++) x[k] = actnode->x[k];
         /*continuity matrix for aktual kin layer at bot*/
         for (jlay=0; jlay<num_klay; jlay++)
         {
            klayhgt = actnode->element[0]->e.s9->klayhgt[jlay];
            hl      = (klayhgt/100)*thick;
            /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
            hl      = A3FAC_SHELL9*hl;
            e3 = -1.0; /*bottom*/
            zeta = s9con(e3,num_klay,klay,jlay,1.0);
            /*calculate coordinates*/
            for (k=0; k<3; k++) x[k] = x[k] + zeta*hl*a3ref_l[k];
         }
         /* write coordinates of first kinematic layer at bottom */
         fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                      node_ID, x[0],x[1],x[2]);
         /*save coordinates for interpolation for the material layers*/
         for (k=0; k<3; k++) x_u[k] = x[k];

         /*calculate the coordinates of the top surface of the each kin layer*/
         for (klay=0; klay<num_klay; klay++)
         {
            /*initialize the coordinates to values of middle surface*/
            for (k=0; k<3; k++) x[k] = actnode->x[k];
            /*continuity matrix for aktual kin layer at top*/
            for (jlay=0; jlay<num_klay; jlay++)
            {
               klayhgt = actnode->element[0]->e.s9->klayhgt[jlay];
               hl      = (klayhgt/100)*thick;
               /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
               hl      = A3FAC_SHELL9*hl;
               e3 = +1.0; /*top*/
               zeta = s9con(e3,num_klay,klay,jlay,1.0);
               /*calculate coordinates*/
               for (k=0; k<3; k++) x[k] = x[k] + zeta*hl*a3ref_l[k];
            }
            /*calculate the coordinates of the material layers of aktual kinematic layer*/
            num_mlay = actnode->element[0]->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actnode->element[0]->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
            sum_lay += num_mlay;
            /*save coordinates for interpolation for the material layers*/
            for (k=0; k<3; k++) x_o[k] = x[k];
            /*interpolate for the material layers within the actual kinematic layer */
            sum_hgt  = 0.0; /*initialize the total hgt*/
            for (mlay=1; mlay<num_mlay; mlay++) /*only if more than one material layer to this kinematic layer*/
            {
               sum_hgt += mlayhgt[mlay-1];
               for (k=0; k<3; k++) x[k] = x_u[k] + (sum_hgt/100.)*(x_o[k]-x_u[k]);
               /* write the bottom coordinates of this material layer  */
               fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                            node_ID+(sum_lay_old+mlay)*numnp, x[0],x[1],x[2]);
            }
            /* write the top coordinates of this kinematic layer */
            fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                         node_ID+sum_lay*numnp, x_o[0],x_o[1],x_o[2]);
            for (k=0; k<3; k++) x_u[k] = x_o[k];
            sum_lay_old = sum_lay;
         }

      }
      else if (type == 8 || type == 9) /* quad8 -> Hex20 / quad9 -> Hex27*/
      {
         if (type == 8)  /*check if actnode is a edge node*/
         {
            is_edge = 0;
            /*check if actnode is a edge node */
            for (m=0; m<actnode->element[0]->numnp; m++)
            {
               if (actnode->element[0]->node[m]->Id == actnode->Id)
               {
                   if (m < 4) is_edge = 1;
               }
            }
         }
         else is_edge = 1; /* for quad9, midnode for all nodes are needed */

         /*init sum_lay: sum of all layers, kinematic and material*/
         sum_lay = 0;
         sum_lay_old = 0;

         /*calculate the coordinates of the bottom surface of the first kinematic layer*/
         klay = 0;
         /*initialize the coordinates to values of middle surface*/
         for (k=0; k<3; k++) x[k] = actnode->x[k];
         /*continuity matrix for aktual kin layer at bot*/
         for (jlay=0; jlay<num_klay; jlay++)
         {
            klayhgt = actnode->element[0]->e.s9->klayhgt[jlay];
            hl      = (klayhgt/100)*thick;
            /*hl      = 0.5*hl; */  /*norm of director is half of kinematic layer hight*/
            hl      = A3FAC_SHELL9*hl;
            e3 = -1.0; /*bottom*/
            zeta = s9con(e3,num_klay,klay,jlay,1.0);
            /*calculate coordinates*/
            for (k=0; k<3; k++) x[k] = x[k] + zeta*hl*a3ref_l[k];
         }
         /* write coordinates of first kinematic layer at bottom */
         fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                      node_ID, x[0],x[1],x[2]);
         /*save coordinates for interpolation for the material layers*/
         for (k=0; k<3; k++) x_u[k] = x[k];

         /*calculate the coordinates of the top surface of the each kin layer*/
         for (klay=0; klay<num_klay; klay++)
         {
            /*initialize the coordinates to values of middle surface*/
            for (k=0; k<3; k++) x[k] = actnode->x[k];
            /*continuity matrix for aktual kin layer at top*/
            for (jlay=0; jlay<num_klay; jlay++)
            {
               klayhgt = actnode->element[0]->e.s9->klayhgt[jlay];
               hl      = (klayhgt/100)*thick;
               /*hl      = 0.5*hl;*/   /*norm of director is half of kinematic layer hight*/
               hl      = A3FAC_SHELL9*hl;
               e3 = +1.0; /*top*/
               zeta = s9con(e3,num_klay,klay,jlay,1.0);
               /*calculate coordinates*/
               for (k=0; k<3; k++) x[k] = x[k] + zeta*hl*a3ref_l[k];
            }
            /*calculate the coordinates of the material layers of aktual kinematic layer*/
            num_mlay = actnode->element[0]->e.s9->kinlay[klay].num_mlay;
            mlayhgt  = actnode->element[0]->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
            sum_lay += num_mlay;
            /*save coordinates for interpolation for the material layers*/
            for (k=0; k<3; k++) x_o[k] = x[k];
            /*interpolate for the material layers within the actual kinematic layer */
            sum_hgt     = 0.0; /*initialize the total hgt*/
            sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/
            for (mlay=0; mlay<num_mlay; mlay++)
            {
               if (mlay > 0) /*if more than one material layer to this kinematic layer*/
               {
                  sum_hgt += mlayhgt[mlay-1];
                  for (k=0; k<3; k++) x[k] = x_u[k] + (sum_hgt/100.)*(x_o[k]-x_u[k]);
                  /* write the bottom coordinates of this material layer  */
                  fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                               node_ID+( (sum_lay_old+mlay)*2 )*numnp, x[0],x[1],x[2]);
               }

               if (is_edge == 1) /*a midnode has to be calculated*/
               {
                  sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
                  for (k=0; k<3; k++) x[k] = x_u[k] + (sum_hgt_mid/100.)*(x_o[k]-x_u[k]);
                  /* write the bottom coordinates of this material layer  */
                  fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                               node_ID+( (sum_lay_old+mlay)*2+1 )*numnp, x[0],x[1],x[2]);
               }
            }
            /* write the top coordinates of this kinematic layer */
            fprintf(out,"%6d %-18.5f %-18.5f %-18.5f\n",
                         node_ID+(sum_lay*2)*numnp, x_o[0],x_o[1],x_o[2]);
            for (k=0; k<3; k++) x_u[k] = x_o[k];
            sum_lay_old = sum_lay;
         }
      }
      else  dserror("wrong element type for gid mesh !!!");
   }
}/* end of (i=0; i<genprob.numfld; i++)*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_allcoords */

/*!----------------------------------------------------------------------
\brief writes the element topology for a multilayerd shell
       element to visualize in gid -> Hexahedra elements

<pre>                                                              sh 012/02
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
void s9_out_gid_eletop(FILE *out,INT type,ELEMENT *actele)
{
INT           i,k,klay,lay;
INT           numnp;            /*number of nodal points*/
INT           numele;           /*number of elements*/
INT           num_klay;
INT           num_mlay;
INT           sum_lay;          /*sum of all layers: kinematic and material*/
INT           node_ID;          /*nodal ID*/
INT           ele_ID;           /*element ID*/
FIELD        *actfield;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_eletop");
#endif
/*----------------------------------------------------------------------*/
ele_ID = actele->Id+1;
num_klay = actele->e.s9->num_klay;

sum_lay = 0;
/*count number of layers to this element*/
for (klay=0; klay<num_klay; klay++)
{
   num_mlay = actele->e.s9->kinlay[klay].num_mlay;
   sum_lay += num_mlay;
}

/*---------------------------------------------------- loop over fields */
for (i=0; i<genprob.numfld; i++)
{
   actfield = &(field[i]);
   numele = actfield->dis[0].numele;
   numnp  = actfield->dis[0].numnp;

   /*switch the elementtype */
   if (type == 4) /*quad4 -> 8-noded Hexahedra element */
   {
      /*--------------------------------- loop over all layers */
      for (lay=0; lay<sum_lay; lay++)
      {
          fprintf(out," %6d ",ele_ID+lay*numele);
          /*untere Ebene*/
          for (k=0; k<actele->numnp; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+lay*numnp);
          }
          /*obere Ebene*/
          for (k=0; k<actele->numnp; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+(lay+1)*numnp);
          }
          fprintf(out,"\n");
      }
   }
   else if (type == 8 || type == 9) /*quad8 -> 20-noded Hexahedra element / quad9 -> Hex27 */
   {
      /*--------------------------------- loop over all layers */
      for (lay=0; lay<sum_lay; lay++)
      {
          fprintf(out," %6d ",ele_ID+lay*numele);
          /*untere Ebene -> Eckknoten*/
          for (k=0; k<4; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+(lay*2)*numnp);
          }
          /*obere Ebene -> Eckknoten*/
          for (k=0; k<4; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+(lay*2+2)*numnp);
          }
          /*untere Ebene -> Seitenmittelknoten*/
          for (k=4; k<8; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+(lay*2)*numnp);
          }
          /*mittlere Ebene -> Eckknoten*/
          for (k=0; k<4; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+(lay*2+1)*numnp);
          }
          /*obere Ebene -> Seitenmittelknoten*/
          for (k=4; k<8; k++)
          {
             node_ID = actele->node[k]->Id+1;
             fprintf(out,"%6d ",node_ID+(lay*2+2)*numnp);
          }

          /* write additional node numbers for quad9 */
          if (type == 9)
          {
             /*untere Ebene -> Mittelknoten*/
             for (k=8; k<9; k++)
             {
                node_ID = actele->node[k]->Id+1;
                fprintf(out,"%6d ",node_ID+(lay*2)*numnp);
             }
             /*mittlere Ebene -> Seitenmittelknoten*/
             for (k=4; k<8; k++)
             {
                node_ID = actele->node[k]->Id+1;
                fprintf(out,"%6d ",node_ID+(lay*2+1)*numnp);
             }
             /*obere Ebene -> Mittelknoten*/
             for (k=8; k<9; k++)
             {
                node_ID = actele->node[k]->Id+1;
                fprintf(out,"%6d ",node_ID+(lay*2+2)*numnp);
             }
             /*mittlere Ebene -> Mittelknoten*/
             for (k=8; k<9; k++)
             {
                node_ID = actele->node[k]->Id+1;
                fprintf(out,"%6d ",node_ID+(lay*2+1)*numnp);
             }
          }

          fprintf(out,"\n");
      }
   }
   else  dserror("wrong element type for gid mesh !!!");
}/* end of (i=0; i<genprob.numfld; i++)*/
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_eletop */

/*!----------------------------------------------------------------------
\brief calculates and writes the nodal displacements of a shell9
       element to visualize in gid -> Hexahedra elements

<pre>                                                              sh 012/02
This routine calculates the nodal displacement from the different vector
component in the different layers and writes them to the flavia.res file
</pre>
\param  *out      FILE    (i)   file to be written on (flavia.res)
\param  *actfield FIELD   (i)   actual field

\warning There is nothing special to this routine
\return void
\sa calling: ---; called by: out_gid_sol()

*-----------------------------------------------------------------------*/
void s9_out_gid_sol_dis(FILE *out,FIELD *actfield,INT place)
{
INT           i,j,m,mlay,klay,jlay;
INT           numnp;                    /*number of nodal points*/
DOUBLE        thick;
INT           num_klay,num_mlay;
DOUBLE       *mlayhgt;
DOUBLE        sum_hgt,sum_hgt_mid;
INT           sum_lay,sum_lay_old;
INT           node_ID;                  /*nodal ID*/
DOUBLE        dis[3],dis_u[3],dis_o[3]; /*displacement x,y,z*/
INT           is_edge;                  /*is_edge=1: actnode is a edge node to the element, else middle node*/
INT           type;
DOUBLE        e3,zeta;
NODE         *actnode;
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_sol_dis");
#endif
/*----------------------------------------------------------------------*/
fprintf(out,"VALUES\n");

numnp = actfield->dis[0].numnp;
for (i=0; i<numnp; i++)
{
   actnode = &(actfield->dis[0].node[i]);

   thick    = actnode->element[0]->e.s9->thick;
   num_klay = actnode->element[0]->e.s9->num_klay;
   node_ID  = actnode->Id+1;
   type     = actnode->element[0]->numnp;

   /*switch the elementtype */
   if (type == 4) /* quad4 -> Hex8 */
   {
      /*init sum_lay: sum of all layers, kinematic and material*/
      sum_lay = 0;
      sum_lay_old = 0;

      /*calculate the displacements of the bottom surface of the first kinematic layer*/
      klay = 0;
      /*initialize the displacements of this layer to zero */
      for(j=0; j<3; j++) dis[j]=0.0;
      /*calculate the relative displacements of this layer */
      for (jlay=0; jlay<num_klay; jlay++)
      {
         e3   = -1.0; /*bottom*/
         zeta = s9con(e3,num_klay,klay,jlay,1.0);
         for (j=0; j<3; j++) dis[j] = dis[j] + zeta*actnode->sol.a.da[place][j+3*(jlay+1)];
      }
      /*add the displacement of reference layer to relative displacement*/
      for (j=0; j<3; j++) dis[j] = dis[j] + actnode->sol.a.da[place][j];

      /* write the displacements of first kinematic layer at bottom */
      fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                   node_ID, dis[0],dis[1],dis[2]);
      /*save coordinates for interpolation for the material layers*/
      for (j=0; j<3; j++) dis_u[j] = dis[j];

      /*calculate the displacements of the top surface of the each kin layer*/
      for (klay=0; klay<num_klay; klay++)
      {
         /*initialize the displacements of this layer to zero */
         for(j=0; j<3; j++) dis[j]=0.0;
         /*calculate the relative displacements of this layer */
         for (jlay=0; jlay<num_klay; jlay++)
         {
            e3   = +1.0; /*top*/
            zeta = s9con(e3,num_klay,klay,jlay,1.0);
            for (j=0; j<3; j++) dis[j] = dis[j] + zeta*actnode->sol.a.da[place][j+3*(jlay+1)];
         }
         /*add the displacement of reference layer to relative displacement*/
         for (j=0; j<3; j++) dis[j] = dis[j] + actnode->sol.a.da[place][j];

         /*calculate the displacements of the material layers of aktual kinematic layer*/
         num_mlay = actnode->element[0]->e.s9->kinlay[klay].num_mlay;
         mlayhgt  = actnode->element[0]->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
         sum_lay += num_mlay;
         /*save displacements for interpolation for the material layers*/
         for (j=0; j<3; j++) dis_o[j] = dis[j];
         /*interpolate for the material layers within the actual kinematic layer */
         sum_hgt  = 0.0; /*initialize the total hgt*/
         for (mlay=1; mlay<num_mlay; mlay++) /*only if more than one material layer to this kinematic layer*/
         {
            sum_hgt += mlayhgt[mlay-1];
            for (j=0; j<3; j++) dis[j] = dis_u[j] + (sum_hgt/100.)*(dis_o[j]-dis_u[j]);
            /* write the bottom displacement of this material layer  */
            fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                         node_ID+(sum_lay_old+mlay)*numnp, dis[0],dis[1],dis[2]);
         }
         /* write the top displacement of this kinematic layer */
         fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                      node_ID+sum_lay*numnp, dis_o[0],dis_o[1],dis_o[2]);
         for (j=0; j<3; j++) dis_u[j] = dis_o[j];
         sum_lay_old = sum_lay;
      }
   }
   else if (type == 8 || type == 9) /* quad8 -> Hex20 / quad9 -> Hex27 */
   {
      if (type == 8)  /*check if actnode is a edge node*/
      {
         is_edge = 0;
         /*check if actnode is a edge node */
         for (m=0; m<actnode->element[0]->numnp; m++)
         {
            if (actnode->element[0]->node[m]->Id == actnode->Id)
            {
                if (m < 4) is_edge = 1;
            }
         }
      }
      else is_edge = 1; /* for quad9, midnode for all nodes are needed */

      /*init sum_lay: sum of all layers, kinematic and material*/
      sum_lay = 0;
      sum_lay_old = 0;

      /*calculate the displacements of the bottom surface of the first kinematic layer*/
      klay = 0;
      /*initialize the displacements of this layer to zero */
      for(j=0; j<3; j++) dis[j]=0.0;
      /*calculate the relative displacements of this layer */
      for (jlay=0; jlay<num_klay; jlay++)
      {
         e3   = -1.0; /*bottom*/
         zeta = s9con(e3,num_klay,klay,jlay,1.0);
         for (j=0; j<3; j++) dis[j] = dis[j] + zeta*actnode->sol.a.da[place][j+3*(jlay+1)];
      }
      /*add the displacement of reference layer to relative displacement*/
      for (j=0; j<3; j++) dis[j] = dis[j] + actnode->sol.a.da[place][j];

      /* write the displacements of first kinematic layer at bottom */
      fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                   node_ID, dis[0],dis[1],dis[2]);
      /*save coordinates for interpolation for the material layers*/
      for (j=0; j<3; j++) dis_u[j] = dis[j];

      /*calculate the displacements of the top surface of the each kin layer*/
      for (klay=0; klay<num_klay; klay++)
      {
         /*initialize the displacements of this layer to zero */
         for(j=0; j<3; j++) dis[j]=0.0;
         /*calculate the relative displacements of this layer */
         for (jlay=0; jlay<num_klay; jlay++)
         {
            e3   = +1.0; /*top*/
            zeta = s9con(e3,num_klay,klay,jlay,1.0);
            for (j=0; j<3; j++) dis[j] = dis[j] + zeta*actnode->sol.a.da[place][j+3*(jlay+1)];
         }
         /*add the displacement of reference layer to relative displacement*/
         for (j=0; j<3; j++) dis[j] = dis[j] + actnode->sol.a.da[place][j];

         /*calculate the displacements of the material layers of aktual kinematic layer*/
         num_mlay = actnode->element[0]->e.s9->kinlay[klay].num_mlay;
         mlayhgt  = actnode->element[0]->e.s9->kinlay[klay].mlayhgt;   /* hgt of material layer in percent of this kin layer */
         sum_lay += num_mlay;
         /*save displacements for interpolation for the material layers*/
         for (j=0; j<3; j++) dis_o[j] = dis[j];
         /*interpolate for the material layers within the actual kinematic layer */
         sum_hgt     = 0.0; /*initialize the total hgt*/
         sum_hgt_mid = 0.0; /*initialize the hgt to a midnode*/
         for (mlay=0; mlay<num_mlay; mlay++)
         {
            if (mlay > 0) /*if more than one material layer to this kinematic layer*/
            {
               sum_hgt += mlayhgt[mlay-1];
               for (j=0; j<3; j++) dis[j] = dis_u[j] + (sum_hgt/100.)*(dis_o[j]-dis_u[j]);
               /* write the bottom displacement of this material layer  */
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            node_ID+( (sum_lay_old+mlay)*2 )*numnp, dis[0],dis[1],dis[2]);
            }

            if (is_edge == 1) /*the displacement at midnode has to be calculated*/
            {
               sum_hgt_mid = sum_hgt + mlayhgt[mlay] * 0.5;
               for (j=0; j<3; j++) dis[j] = dis_u[j] + (sum_hgt_mid/100.)*(dis_o[j]-dis_u[j]);
               /* write the bottom displacement of this material layer  */
               fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                            node_ID+( (sum_lay_old+mlay)*2+1 )*numnp, dis[0],dis[1],dis[2]);
            }
         }
         /* write the top displacement of this kinematic layer */
         fprintf(out," %6d %18.5E %18.5E %18.5E\n",
                      node_ID+(sum_lay*2)*numnp, dis_o[0],dis_o[1],dis_o[2]);
         for (j=0; j<3; j++) dis_u[j] = dis_o[j];
         sum_lay_old = sum_lay;
      }
   }
   else  dserror("wrong element type for gid out sol (s9) !!!");
}

fprintf(out,"END VALUES\n");

/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_sol_dis */


/*!----------------------------------------------------------------------
\brief writes the physical stresses, extrapolated from GP to nodal points
 to the flavia.res file -> Hexahedra elements

<pre>                                                              sh 1/03
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
void s9_out_gid_sol_str(FILE *out,FIELD *actfield,INT place)
{
INT           i,j,k,l,m,klay,lay;
INT           numnp;            /*number of nodal points in actfield */
INT           numele;           /*number of elements to actnode -> in plane*/
INT           num_node;         /*number of nodes to an element*/
DOUBLE      **stress;
DOUBLE        sum_str_u, sum_str_oo, sum_str_uo;
DOUBLE        strK[6],strKo[6],strK_save[6],strK_mid[6];
INT           node_ID;                  /*nodal ID*/
INT           num_klay;
INT           num_mlay;
INT           sum_lay;
NODE         *actnode;
INT           is_edge;                  /*is_edge=1: actnode is a edge node to the element, else middle node*/


/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_enter("s9_out_gid_sol_str");
#endif
/*----------------------------------------------------------------------*/

numnp  = actfield->dis[0].numnp;
for (i=0; i<numnp; i++)
{
   actnode  = &(actfield->dis[0].node[i]);
   numele   = actnode->numele;
   node_ID  = actnode->Id+1;
   num_node = actnode->element[0]->numnp;          /*number of nodes to first element conected to actnode*/
   num_klay = actnode->element[0]->e.s9->num_klay; /*number of kinematic layers to first element conected to actnode*/

   /*------------------------------ count number of layers to this element */
   sum_lay = 0;
   for (klay=0; klay<num_klay; klay++)
   {
      num_mlay = actnode->element[0]->e.s9->kinlay[klay].num_mlay;
      sum_lay += num_mlay;
   }

   /*switch the type */
   if (num_node == 4) /* -> quad4 */
   {
      for (lay=0; lay<sum_lay; lay++)
      {
         for (j=0; j<6; j++)    /*6-stress Komponents*/
         {
            sum_str_u = 0.0;
            sum_str_uo = 0.0;
            sum_str_oo = 0.0;
            for (k=0; k<numele; k++) /*loop over all adjacent elements to actnode*/
            {
               for (l=0; l<num_node; l++)   /*loop over all nodes of the element*/
               {
                  if (actnode->element[k]->node[l]->Id == actnode->Id)
                  {
                      stress  = actnode->element[k]->e.s9->stresses.a.d3[place];
                      sum_str_u  += stress[j][l + (2* lay    + 0) * num_node]; /*lt = 0 von lay  */
                    if (lay == sum_lay-1)
                      sum_str_oo += stress[j][l + (2* lay    + 1) * num_node]; /*lt = 1 von lay   -> oberste Schicht */
                    if (lay != 0)
                      sum_str_uo += stress[j][l + (2*(lay-1) + 1) * num_node]; /*lt = 1 von lay-1 -> mittlere Schichten*/
                  }
               }
            }
            strKo[j] = sum_str_oo / numele;                  /*oberste Schicht*/
            if (lay == 0) strK[j] =  sum_str_u  / numele;                  /*unterste Schicht*/
            else          strK[j] = (sum_str_u + sum_str_uo) / (2*numele); /*mittlere Schichten*/
         }
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                       node_ID + (lay ) *numnp   ,strK[0],strK[1],strK[2],strK[3],strK[4],strK[5]);
         if (lay == sum_lay-1)
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                       node_ID + (lay+1)*numnp  ,strKo[0],strKo[1],strKo[2],strKo[3],strKo[4],strKo[5]);

      }/*====================================================== end loop over all layers*/
   }

   else if (num_node == 8 || num_node == 9) /*quad8 or quad9*/
   {

      if (num_node == 8)  /* -> quad8 check if actnode is a edge node*/
      {
         is_edge = 0;
         /*check if actnode is a edge node */
         for (m=0; m<actnode->element[0]->numnp; m++)
         {
            if (actnode->element[0]->node[m]->Id == actnode->Id)
            {
                if (m < 4) is_edge = 1;
            }
         }
      }
      else is_edge = 1; /* for quad9, midnode for all nodes are needed */

      for (lay=0; lay<sum_lay; lay++)
      {
         for (j=0; j<6; j++)    /*6-stress Komponents*/
         {
            sum_str_u = 0.0;
            sum_str_uo = 0.0;
            sum_str_oo = 0.0;
            for (k=0; k<numele; k++) /*loop over all adjacent elements to actnode*/
            {
               for (l=0; l<num_node; l++)   /*loop over all nodes of the element*/
               {
                  if (actnode->element[k]->node[l]->Id == actnode->Id)
                  {
                      stress  = actnode->element[k]->e.s9->stresses.a.d3[place];
                      sum_str_u  += stress[j][l + (2* lay    + 0) * 9]; /*lt = 0 von lay  */
                    if (lay == sum_lay-1)
                      sum_str_oo += stress[j][l + (2* lay    + 1) * 9]; /*lt = 1 von lay   -> oberste Schicht */
                    if (lay != 0)
                      sum_str_uo += stress[j][l + (2*(lay-1) + 1) * 9]; /*lt = 1 von lay-1 -> mittlere Schichten*/
                  }
               }
            }
            strKo[j] = sum_str_oo / numele;                  /*oberste Schicht*/
            if (lay == 0) strK[j] =  sum_str_u  / numele;                  /*unterste Schicht*/
            else          strK[j] = (sum_str_u + sum_str_uo) / (2*numele); /*mittlere Schichten*/
         }
         /*write middle of actual layer - 1*/
         if (lay != 0)
         {
            if (is_edge == 1) /*average between top and bottom of a layer*/
            {
               for (m=0; m<6; m++) strK_mid[m] = 0.5*(strK[m] + strK_save[m]);
               fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             node_ID + ((lay-1)*2 +1) *numnp   ,strK_mid[0],strK_mid[1],strK_mid[2],strK_mid[3],strK_mid[4],strK_mid[5]);
            }
         }
         /*write bottom of actual layer*/
         fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                       node_ID + (lay*2 ) *numnp   ,strK[0],strK[1],strK[2],strK[3],strK[4],strK[5]);
         for (m=0; m<6; m++) strK_save[m] = strK[m];
         /*write the upper most layer*/
         if (lay == sum_lay-1)
         {
            if (is_edge == 1) /*average between top and bottom of a layer*/
            {
               for (m=0; m<6; m++) strK_mid[m] = 0.5*(strKo[m] + strK_save[m]);
               fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                             node_ID + ((lay)*2 +1) *numnp   ,strK_mid[0],strK_mid[1],strK_mid[2],strK_mid[3],strK_mid[4],strK_mid[5]);
            }
            /*write top of the upper most layer*/
            fprintf(out," %6d %18.5e %18.5e %18.5e %18.5e %18.5e %18.5e \n",
                          node_ID + (lay*2+2)*numnp  ,strKo[0],strKo[1],strKo[2],strKo[3],strKo[4],strKo[5]);
         }

      }/*====================================================== end loop over all layers*/
   }
   else  dserror("wrong gp type for gid out sol (s9) !!!");
}
/*----------------------------------------------------------------------*/
#ifdef DEBUG
dstrc_exit();
#endif
return;
} /* end of s9_out_gid_sol_str */
/*----------------------------------------------------------------------*/
#endif /*D_SHELL9*/
/*! @} (documentation module close)*/

