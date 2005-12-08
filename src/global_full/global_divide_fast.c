/*-----------------------------------------------------------------------*/
/*!
\file
\brief contains the routine 'divide_fast' which sorts the fast elements and
       assignes them into the sets.


<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

 */
/*-----------------------------------------------------------------------*/

/*!
\addtogroup GLOBAL
*//*! @{ (documentation module open)*/

#ifdef D_FLUID3_F


#include "../headers/standardtypes.h"
#include "../fluid3/fluid3.h"

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

/*!----------------------------------------------------------------------
\brief one proc's info about his partition

<pre>                                                         m.gee 8/00
-the partition of one proc (all discretizations)
-the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PARTITION  *partition;


/*-----------------------------------------------------------------------*/
/*!
  \brief assignes the fast elements into the sets

  This routine first loops all fast elements and counts how many of ech typ
  are present in the problem. After that the sets for these elements are
  allocated and the elements are looped again, to assign each fast elements
  into a set.


  \return void

  \author mn
  \date   10/04
 */
/*-----------------------------------------------------------------------*/
void divide_fast(
    )
{
  INT         i,j,k;

  PARTITION      *actpart;            /* pointer to active partition */
  PARTDISCRET    *actpdis;

  ELEMENT    *actele;

  INT         num_f3f_hex8_e = 0;
  INT         vec_f3f_hex8_e = 0;
  INT         cur_vec_f3f_hex8_e = 0;
  INT         cur_ele_f3f_hex8_e = 0;

  INT         num_f3f_hex20_e = 0;
  INT         vec_f3f_hex20_e = 0;
  INT         cur_vec_f3f_hex20_e = 0;
  INT         cur_ele_f3f_hex20_e = 0;

  INT         num_f3f_tet4_e = 0;
  INT         vec_f3f_tet4_e = 0;
  INT         cur_vec_f3f_tet4_e = 0;
  INT         cur_ele_f3f_tet4_e = 0;

  INT         num_f3f_hex8_a = 0;
  INT         vec_f3f_hex8_a = 0;
  INT         cur_vec_f3f_hex8_a = 0;
  INT         cur_ele_f3f_hex8_a = 0;

  INT         num_f3f_hex20_a = 0;
  INT         vec_f3f_hex20_a = 0;
  INT         cur_vec_f3f_hex20_a = 0;
  INT         cur_ele_f3f_hex20_a = 0;

  INT         num_f3f_tet4_a = 0;
  INT         vec_f3f_tet4_a = 0;
  INT         cur_vec_f3f_tet4_a = 0;
  INT         cur_ele_f3f_tet4_a = 0;

  INT         count_vec = 0;

#ifdef DEBUG
  dstrc_enter("divide_fast");
#endif



for(i=0; i<genprob.numfld; i++)  /* loop all fields */
{
  actpart  = &(partition[i]);

  if (actpart->fieldtyp != fluid)
    continue;

  for(j=0; j<actpart->ndis; j++)  /* loop all discretisations */
  {

    count_vec = 0;

    num_f3f_hex8_e = 0;
    vec_f3f_hex8_e = 0;
    cur_vec_f3f_hex8_e = 0;
    cur_ele_f3f_hex8_e = 0;

    num_f3f_hex20_e = 0;
    vec_f3f_hex20_e = 0;
    cur_vec_f3f_hex20_e = 0;
    cur_ele_f3f_hex20_e = 0;

    num_f3f_tet4_e = 0;
    vec_f3f_tet4_e = 0;
    cur_vec_f3f_tet4_e = 0;
    cur_ele_f3f_tet4_e = 0;

    num_f3f_hex8_a = 0;
    vec_f3f_hex8_a = 0;
    cur_vec_f3f_hex8_a = 0;
    cur_ele_f3f_hex8_a = 0;

    num_f3f_hex20_a = 0;
    vec_f3f_hex20_a = 0;
    cur_vec_f3f_hex20_a = 0;
    cur_ele_f3f_hex20_a = 0;

    num_f3f_tet4_a = 0;
    vec_f3f_tet4_a = 0;
    cur_vec_f3f_tet4_a = 0;
    cur_ele_f3f_tet4_a = 0;

    actpdis  =  &(actpart->pdis[j]);


    printf("divide_fast!! field: %d;  dis: %d\n\n",i,j);

    actpdis->num_fele   = 0;
    num_f3f_hex8_e      = 0;
    num_f3f_hex20_e     = 0;
    num_f3f_tet4_e      = 0;
    num_f3f_hex8_e      = 0;
    num_f3f_hex20_a     = 0;
    num_f3f_tet4_a      = 0;

    for (k=0; k<actpdis->numele; k++)  /* loop all elements */
    {

      actele = actpdis->element[k];

      /* count the fluid3_fast elements here */
      if (actele->eltyp == el_fluid3_fast)
      {
        switch (actele->distyp)
        {
          case hex8:
            if (actele->e.f3->is_ale == 1)
              num_f3f_hex8_a += 1;
            else
              num_f3f_hex8_e += 1;
            break;

          case hex20:
            if (actele->e.f3->is_ale == 1)
              num_f3f_hex20_a += 1;
            else
              num_f3f_hex20_e += 1;
            break;

          case tet4:
            if (actele->e.f3->is_ale == 1)
              num_f3f_tet4_a += 1;
            else
              num_f3f_tet4_e += 1;
            break;

          case tet10:
            dserror("tet10 not yet possible for fluid3_fast");
            break;

          default:
            dserror("unknown element typ for fluid3_fast");
            break;
        } /* switch (actele->distyp) */
      } /* if (actele->eltyp == el_fluid3_fast) */

    } /* end of loop all elements */





    /* calculate number of vectors necessary ... */
    printf("LOOPL          = %5i\n\n",LOOPL);

    /* ... for fluid3 hex8 on euler */
    vec_f3f_hex8_e = num_f3f_hex8_e/LOOPL;
    if (num_f3f_hex8_e%LOOPL != 0)
      vec_f3f_hex8_e += 1;
    if (num_f3f_hex8_e != 0)
    {
      printf("num_f3f_hex8_e = %5i\n",num_f3f_hex8_e);
      printf("vec_f3f_hex8_e = %5i\n",vec_f3f_hex8_e);
    }

    /* ... for fluid3 hex8 on ale */
    vec_f3f_hex8_a = num_f3f_hex8_a/LOOPL;
    if (num_f3f_hex8_a%LOOPL != 0)
      vec_f3f_hex8_a += 1;
    if (num_f3f_hex8_a != 0)
    {
      printf("num_f3f_hex8_a = %5i\n",num_f3f_hex8_a);
      printf("vec_f3f_hex8_a = %5i\n",vec_f3f_hex8_a);
    }

    /* ... for fluid3 hex20 on euler */
    vec_f3f_hex20_e = num_f3f_hex20_e/LOOPL;
    if (num_f3f_hex20_e%LOOPL != 0)
      vec_f3f_hex20_e += 1;
    if (num_f3f_hex20_e != 0)
    {
      printf("num_f3f_hex20_e = %5i\n",num_f3f_hex20_e);
      printf("vec_f3f_hex20_e = %5i\n",vec_f3f_hex20_e);
    }

    /* ... for fluid3 hex20 on ale */
    vec_f3f_hex20_a = num_f3f_hex20_a/LOOPL;
    if (num_f3f_hex20_a%LOOPL != 0)
      vec_f3f_hex20_a += 1;
    if (num_f3f_hex20_a != 0)
    {
      printf("num_f3f_hex20_a = %5i\n",num_f3f_hex20_a);
      printf("vec_f3f_hex20_a = %5i\n",vec_f3f_hex20_a);
    }

    /* ... for fluid3 tet4 on euler */
    vec_f3f_tet4_e = num_f3f_tet4_e/LOOPL;
    if (num_f3f_tet4_e%LOOPL != 0)
      vec_f3f_tet4_e += 1;
    if (num_f3f_tet4_e != 0)
    {
      printf("num_f3f_tet4_e = %5i\n",num_f3f_tet4_e);
      printf("vec_f3f_tet4_e = %5i\n",vec_f3f_tet4_e);
    }

    /* ... for fluid3 tet4 on ale */
    vec_f3f_tet4_a = num_f3f_tet4_a/LOOPL;
    if (num_f3f_tet4_a%LOOPL != 0)
      vec_f3f_tet4_a += 1;
    if (num_f3f_tet4_a != 0)
    {
      printf("num_f3f_tet4_a = %5i\n",num_f3f_tet4_a);
      printf("vec_f3f_tet4_a = %5i\n",vec_f3f_tet4_a);
    }



    /* total number of vectors with fast elements */
    actpdis->num_fele = vec_f3f_hex8_e + vec_f3f_hex8_a +
      vec_f3f_hex20_e + vec_f3f_hex20_a +
      vec_f3f_tet4_e + vec_f3f_tet4_a;


    /* allocate the structures */
    actpdis->fast_eles = (FAST_ELES*)CCAMALLOC((actpdis->num_fele*sizeof(FAST_ELES)));



    /* allocate the vectors ... */

    /* ... for fluid3 hex8 on euler */
    if (vec_f3f_hex8_e != 0)
    {
      cur_vec_f3f_hex8_e = count_vec;
      for (k=0;k<vec_f3f_hex8_e;k++)
      {
        actpdis->fast_eles[count_vec].fast_ele_typ = fele_f3f_hex8_e;
        actpdis->fast_eles[count_vec].aloopl  = 0;
        actpdis->fast_eles[count_vec].ele_vec =
          (ELEMENT**) CCAMALLOC(LOOPL*sizeof(ELEMENT*));
        count_vec++;
      }
    }

    /* ... for fluid3 hex8 on ale */
    if (vec_f3f_hex8_a != 0)
    {
      cur_vec_f3f_hex8_a = count_vec;
      for (k=0;k<vec_f3f_hex8_a;k++)
      {
        actpdis->fast_eles[count_vec].fast_ele_typ = fele_f3f_hex8_a;
        actpdis->fast_eles[count_vec].aloopl  = 0;
        actpdis->fast_eles[count_vec].ele_vec =
          (ELEMENT**) CCAMALLOC(LOOPL*sizeof(ELEMENT*));
        count_vec++;
      }
    }

    /* ... for fluid3 hex20 on euler */
    if (vec_f3f_hex20_e != 0)
    {
      cur_vec_f3f_hex20_e = count_vec;
      for (k=0;k<vec_f3f_hex20_e;k++)
      {
        actpdis->fast_eles[count_vec].fast_ele_typ = fele_f3f_hex20_e;
        actpdis->fast_eles[count_vec].aloopl  = 0;
        actpdis->fast_eles[count_vec].ele_vec =
          (ELEMENT**) CCAMALLOC(LOOPL*sizeof(ELEMENT*));
        count_vec++;
      }
    }

    /* ... for fluid3 hex20 on ale */
    if (vec_f3f_hex20_a != 0)
    {
      cur_vec_f3f_hex20_a = count_vec;
      for (k=0;k<vec_f3f_hex20_a;k++)
      {
        actpdis->fast_eles[count_vec].fast_ele_typ = fele_f3f_hex20_a;
        actpdis->fast_eles[count_vec].aloopl  = 0;
        actpdis->fast_eles[count_vec].ele_vec =
          (ELEMENT**) CCAMALLOC(LOOPL*sizeof(ELEMENT*));
        count_vec++;
      }
    }

    /* ... for fluid3 tet4 on euler */
    if (vec_f3f_tet4_e != 0)
    {
      cur_vec_f3f_tet4_e = count_vec;
      for (k=0;k<vec_f3f_tet4_e;k++)
      {
        actpdis->fast_eles[count_vec].fast_ele_typ = fele_f3f_tet4_e;
        actpdis->fast_eles[count_vec].aloopl  = 0;
        actpdis->fast_eles[count_vec].ele_vec =
          (ELEMENT**) CCAMALLOC(LOOPL*sizeof(ELEMENT*));
        count_vec++;
      }
    }

    /* ... for fluid3 tet4 on ale */
    if (vec_f3f_tet4_a != 0)
    {
      cur_vec_f3f_tet4_a = count_vec;
      for (k=0;k<vec_f3f_tet4_a;k++)
      {
        actpdis->fast_eles[count_vec].fast_ele_typ = fele_f3f_tet4_a;
        actpdis->fast_eles[count_vec].aloopl  = 0;
        actpdis->fast_eles[count_vec].ele_vec =
          (ELEMENT**) CCAMALLOC(LOOPL*sizeof(ELEMENT*));
        count_vec++;
      }
    }


    if(count_vec != actpdis->num_fele)
      dserror("Irgendetwas ist beim Allokieren der Vektoren schief gelaufen...");


    /* and fill them ... */
    for (k=0; k<actpdis->numele; k++)  /* loop all elements */
    {

      actele = actpdis->element[k];

      if (actele->eltyp == el_fluid3_fast)
      {
        switch (actele->distyp)
        {
          case hex8:

            if (actele->e.f3->is_ale == 1)
            {
              /* ... with fluid3 hex8 on ale */
              actpdis->fast_eles[cur_vec_f3f_hex8_a].ele_vec[cur_ele_f3f_hex8_a] =
                actele;
              cur_ele_f3f_hex8_a++;
              /* check if vector is full */
              if (cur_ele_f3f_hex8_a == LOOPL)
              {
                actpdis->fast_eles[cur_vec_f3f_hex8_a].aloopl = cur_ele_f3f_hex8_a;
                cur_vec_f3f_hex8_a++;
                cur_ele_f3f_hex8_a = 0;
              }
            }

            else
            {
              /* ... with fluid3 hex8 on euler */
              actpdis->fast_eles[cur_vec_f3f_hex8_e].ele_vec[cur_ele_f3f_hex8_e] =
                actele;
              cur_ele_f3f_hex8_e++;
              /* check if vector is full */
              if (cur_ele_f3f_hex8_e == LOOPL)
              {
                actpdis->fast_eles[cur_vec_f3f_hex8_e].aloopl = cur_ele_f3f_hex8_e;
                cur_vec_f3f_hex8_e++;
                cur_ele_f3f_hex8_e = 0;
              }
            }
            break;

          case hex20:

            if (actele->e.f3->is_ale == 1)
            {
              /* ... with fluid3 hex20 on ale */
              actpdis->fast_eles[cur_vec_f3f_hex20_a].ele_vec[cur_ele_f3f_hex20_a] =
                actele;
              cur_ele_f3f_hex20_a++;
              /* check if vector is full */
              if (cur_ele_f3f_hex20_a == LOOPL)
              {
                actpdis->fast_eles[cur_vec_f3f_hex20_a].aloopl = cur_ele_f3f_hex20_a;
                cur_vec_f3f_hex20_a++;
                cur_ele_f3f_hex20_a = 0;
              }
            }

            else
            {
              /* ... with fluid3 hex20 on euler */
              actpdis->fast_eles[cur_vec_f3f_hex20_e].ele_vec[cur_ele_f3f_hex20_e] =
                actele;
              cur_ele_f3f_hex20_e++;
              /* check if vector is full */
              if (cur_ele_f3f_hex20_e == LOOPL)
              {
                actpdis->fast_eles[cur_vec_f3f_hex20_e].aloopl = cur_ele_f3f_hex20_e;
                cur_vec_f3f_hex20_e++;
                cur_ele_f3f_hex20_e = 0;
              }
            }
            break;

          case tet4:
            if (actele->e.f3->is_ale == 1)
            {
              /* ... with fluid3 tet4 on ale */
              actpdis->fast_eles[cur_vec_f3f_tet4_a].ele_vec[cur_ele_f3f_tet4_a] =
                actele;
              cur_ele_f3f_tet4_a++;
              /* check if vector is full */
              if (cur_ele_f3f_tet4_a == LOOPL)
              {
                actpdis->fast_eles[cur_vec_f3f_tet4_a].aloopl = cur_ele_f3f_tet4_a;
                cur_vec_f3f_tet4_a++;
                cur_ele_f3f_tet4_a = 0;
              }
            }
            else
            {
              /* ... with fluid3 tet4 on euler */
              actpdis->fast_eles[cur_vec_f3f_tet4_e].ele_vec[cur_ele_f3f_tet4_e] =
                actele;
              cur_ele_f3f_tet4_e++;
              /* check if vector is full */
              if (cur_ele_f3f_tet4_e == LOOPL)
              {
                actpdis->fast_eles[cur_vec_f3f_tet4_e].aloopl = cur_ele_f3f_tet4_e;
                cur_vec_f3f_tet4_e++;
                cur_ele_f3f_tet4_e = 0;
              }
            }
            break;

          case tet10:
            dserror("tet10 not yet possible for fluid3_fast");
            break;
          default:
            dserror("unknown element typ for fluid3_fast");
            break;
        } /* switch (actele->distyp) */

      } /* if (actele->eltyp == el_fluid3_fast) */

    } /* end of loop all elements */


    if (vec_f3f_hex8_e != 0 && cur_ele_f3f_hex8_e != 0)
      actpdis->fast_eles[cur_vec_f3f_hex8_e].aloopl = cur_ele_f3f_hex8_e;

    if (vec_f3f_hex8_a != 0 && cur_ele_f3f_hex8_a != 0)
      actpdis->fast_eles[cur_vec_f3f_hex8_a].aloopl = cur_ele_f3f_hex8_a;

    if (vec_f3f_hex20_e != 0 && cur_ele_f3f_hex20_e != 0)
      actpdis->fast_eles[cur_vec_f3f_hex20_e].aloopl = cur_ele_f3f_hex20_e;

    if (vec_f3f_hex20_a != 0 && cur_ele_f3f_hex20_a != 0)
      actpdis->fast_eles[cur_vec_f3f_hex20_a].aloopl = cur_ele_f3f_hex20_a;

    if (vec_f3f_tet4_e != 0 && cur_ele_f3f_tet4_e != 0)
      actpdis->fast_eles[cur_vec_f3f_tet4_e].aloopl = cur_ele_f3f_tet4_e;

    if (vec_f3f_tet4_a != 0 && cur_ele_f3f_tet4_a != 0)
      actpdis->fast_eles[cur_vec_f3f_tet4_a].aloopl = cur_ele_f3f_tet4_a;


  }  /* for(j=0; j<actpart->ndis; j++) */

}  /* for(i=0; i<genprob.numfld; i++) */




#ifdef DEBUG
  dstrc_exit();
#endif

  return;
} /* end of divide_fast */



#endif /* ifdef D_FLUID3_F */

/*! @} (documentation module close)*/

