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


#ifndef MATH_PROTOTYPES_H
#define MATH_PROTOTYPES_H

/*----------------------------------------------------------------------*
 |  math1.c                                               m.gee 11/01   |
 *----------------------------------------------------------------------*/
void math_array_copy(DOUBLE **from, INT n, INT m, DOUBLE **to);
void math_inv3(DOUBLE **a, DOUBLE *det);
void math_tran(DOUBLE **a, INT n);
void math_unvc(DOUBLE *enorm,DOUBLE *vec, INT n);
void math_matvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor);
void math_mattrnvecdense(DOUBLE  *r,
                         DOUBLE **A,
                         DOUBLE  *b,
                         INT      ni,
                         INT      nk,
                         INT      init,
                         DOUBLE   factor);
void math_matmatdense(DOUBLE **R,
                         DOUBLE **A,
                         DOUBLE **B,
                         INT      ni,
                         INT      nk,
                         INT      nj,
                         INT      init,
                         DOUBLE   factor);
void math_mattrnmatdense(DOUBLE **R,
                            DOUBLE **A,
                            DOUBLE **B,
                            INT      ni,
                            INT      nk,
                            INT      nj,
                            INT      init,
                            DOUBLE   factor);
void math_matmattrndense(DOUBLE **R,
                            DOUBLE **A,
                            DOUBLE **B,
                            INT      ni,
                            INT      nk,
                            INT      nj,
                            INT      init,
                            DOUBLE   factor);
void math_sym_inv(DOUBLE **A, INT dim);
void math_unsym_inv(DOUBLE **A, INT dimr, INT dimc);
void math_unsym_inv6x6(DOUBLE **A);
void math_sppr(DOUBLE *spat, DOUBLE *a, DOUBLE *b, DOUBLE *c);
void math_addab(DOUBLE **a, DOUBLE **b, INT dim1, INT dim2, DOUBLE fact);



/*----------------------------------------------------------------------*
 |  geometry.c                                           chfoe 09/04    |
 *----------------------------------------------------------------------*/
DOUBLE area_lin_2d(ELEMENT *ele, DOUBLE xyz[2][MAXNOD]);
 
/*!---------------------------------------------------------------------
\brief extract digits from integer number

<pre>                                                         genk 04/02
</pre>
\param  num	 INT   (i)    integer number
\param *it	 INT   (o)    integer on position "thousand"
\param *ih       INT   (o)    integer on position "hundred"
\param *id       INT   (o)    integer on position "ten"
\param *id       INT   (o)    integer on position "one"
\return void

------------------------------------------------------------------------*/
void math_intextract(
                    INT num,
                    INT *it,
		    INT *ih,
		    INT *id,
		    INT *io
	            );

/*----------------------------------------------------------------------*
 |  sort_find.c                                          m.gee 11/01    |
 *----------------------------------------------------------------------*/
INT cmp_int(const void *a, const void *b );
DOUBLE cmp_double(const void *a, const void *b );
void mg_sort(INT list[], INT N, INT list2[], DOUBLE list3[]);
INT quick_find(INT key, INT list[], INT length, INT shift, INT bins[]);
void init_quick_find(INT list[], INT length, INT *shift, INT *bins);
INT find_index(INT key, INT list[], INT length);


#endif

