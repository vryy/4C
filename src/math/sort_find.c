/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

------------------------------------------------------------------------*/
#ifndef CCADISCRET
#include "../headers/standardtypes.h"

/*----------------------------------------------------------------------*
 |  compare the integers - qsort routine                  a.lipka 5/01  |
 |                                                                      |
 |  the call for the sorter of an INT vector is then                    |
 |                                                                      |
 |  qsort((INT*) vector, lenght, sizeof(INT), cmp_int);                 |
 |                                                                      |
 *----------------------------------------------------------------------*/
INT cmp_int(const void *a, const void *b )
{
    return *(INT *)a - * (INT *)b;
}
/*----------------------------------------------------------------------*
 |  compare the doubles - qsort routine                   a.lipka 5/01  |
 |                                                                      |
 |  the call for the sorter of a DOUBLE vector is then                  |
 |                                                                      |
 |  qsort((DOUBLE*) vector, lenght, sizeof(DOUBLE), cmp_double);        |
 |                                                                      |
 *----------------------------------------------------------------------*/
DOUBLE cmp_double(const void *a, const void *b )
{
    return *(DOUBLE *)a - * (DOUBLE *)b;
}


/******************************************************************************/

void mg_sort(INT list[], INT N, INT list2[], DOUBLE list3[])

/*******************************************************************************

  This routine was taken from Knuth: Sorting and Searching. It puts the input
  data list into a heap and then sorts it.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     none
  ============

  Parameter list:
  ===============

  list:            On input, values to be sorted. On output, sorted values
                   (i.e., list[i] <= list[i+1]).

  N:               length of vector 'vec'.

  list2:           If on input,
                   a) list2 = NULL: it is unchanged on output,
                   b) list2 is a list associated with 'list':
                   on output, if list[k] on input is now element 'j' on output,
                   list2[j] on output is list2[k].

  list3:           If on input,
                   a) list3 = NULL: it is unchanged on output,
                   b) list3 is a list associated with 'list':
                   on output, list3[j] is assigned the input value of list3[k],
                   if list[j] has been assigned the input value of list[k].

*******************************************************************************/

{

  /* local variables */

  INT    l, r, RR, K, j, i, flag;
  INT    RR2;
  DOUBLE RR3;

  /**************************** execution begins ******************************/

  if (N <= 1) return;

  l   = N / 2 + 1;
  r   = N - 1;
  l   = l - 1;
  RR  = list[l - 1];
  K   = list[l - 1];

  if ((list2 != NULL) && (list3 != NULL)) {
    RR2 = list2[l - 1];
    RR3 = list3[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list2[i - 1] = list2[j - 1];
            list3[i - 1] = list3[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list2[i - 1] = RR2;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list [r];
        RR2 = list2[r];
        RR3 = list3[r];

        K = list[r];
        list[r ] = list[0];
        list2[r] = list2[0];
        list3[r] = list3[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR2 = list2[l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list2[0] = RR2;
    list3[0] = RR3;
  }
  else if (list2 != NULL) {
    RR2 = list2[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list2[i - 1] = list2[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list2[i - 1] = RR2;

      if (l == 1) {
        RR  = list [r];
        RR2 = list2[r];

        K = list[r];
        list[r ] = list[0];
        list2[r] = list2[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR2 = list2[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list2[0] = RR2;
  }
  else if (list3 != NULL) {
    RR3 = list3[l - 1];
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
            list3[i - 1] = list3[j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;
      list3[i - 1] = RR3;

      if (l == 1) {
        RR  = list [r];
        RR3 = list3[r];

        K = list[r];
        list[r ] = list[0];
        list3[r] = list3[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        RR3 = list3[l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
    list3[0] = RR3;

  }
  else {
    while (r != 0) {
      j = l;
      flag = 1;

      while (flag == 1) {
        i = j;
        j = j + j;

        if (j > r + 1)
          flag = 0;
        else {
          if (j < r + 1)
            if (list[j] > list[j - 1]) j = j + 1;

          if (list[j - 1] > K) {
            list[ i - 1] = list[ j - 1];
          }
          else {
            flag = 0;
          }
        }
      }

      list[ i - 1] = RR;

      if (l == 1) {
        RR  = list [r];

        K = list[r];
        list[r ] = list[0];
        r = r - 1;
      }
      else {
        l   = l - 1;
        RR  = list[ l - 1];
        K   = list[l - 1];
      }
    }

    list[ 0] = RR;
  }

} /* mg_sort */


/******************************************************************************/

INT quick_find(INT key, INT list[], INT length, INT shift, INT bins[])

/*******************************************************************************

  Find 'key' in 'list' and return the indices number. On exit, find_index()
  returns:
            -1 ==> key not found
             i ==> list[i] = key

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List (assumed to be in ascending order) to be searched.

  length:          Length of list.

  shift:           Integer used to compute an indices into bins[]. Computed by
                   init_quick_find().

  bins:            Used to obtain a sublist where should contain 'key'. In
                   particular, list[bins[k] ... bins[k+1]-1] should contain
                   'key'. Computed by init_quick_find().

*******************************************************************************/

{

  /* local variables */

  INT i, loc, oldkey;

  /**************************** execution begins ******************************/

  if (length == 0)            return -1;
  if (key > list[length - 1]) return -1;

  oldkey = key;
  key   -= list[0];

  if (key < 0) return -1;

  loc = key >> shift;

  i = find_index(oldkey, &(list[bins[loc]]), bins[loc + 1] - bins[loc]);

  if (i == -1) return -1;

  return (i + bins[loc]);

} /* quick_find */
/******************************************************************************/

void init_quick_find(INT list[], INT length, INT *shift, INT *bins)

/*******************************************************************************

  The value 'shift' and the array 'bins' are initialized so that subsequent
  calls to quick_find() will work properly. In particular, the array 'bins'
  and which should be 1/4 the size of list is set so that

    1)  range>>shift     > length/4
              and
        range>>(shift+1) < length/4

    where range = list[length-1] - list[0].

    2)  list[j] = value  ==> bins[k] <= j < bins[k+1]
                             where k = (value-list[0]) >> shift

 Note: list[] is sorted in ascending order. This routine is used in conjunction
 with quick_find(). The idea is to use bins[] to get a good initial guess as
 to the location of 'value' in list[].

  Author:          Ray S. Tuminaro, SNL, 1422
  =======

  Return code:     void
  ============

  Parameter list:
  ===============

  list:            List (assumed to be in ascending order) to be searched.

  length:          Length of list.

  shift:           Integer used to compute an indices into bins[]. Computed by
                   init_quick_find().

  bins:            Used to obtain a sublist where should contain 'key'. In
                   particular, list[bins[k] ... bins[k+1]-1] should contain
                   'key'. Computed by init_quick_find().

*******************************************************************************/

{

  /* local variables */

  register INT i, j = 0;
  INT          range, temp;

  /**************************** execution begins ******************************/

  if (length == 0) return;

  range  = list[length - 1] - list[0];
  *shift = 0;

  while ((range >> (*shift)) > length / 4)
    (*shift)++;

  bins[j++] = 0;

  for (i = 0; i < length; i++) {
    temp = list[i] - list[0];

    while ((temp >> (*shift)) >= j)
      bins[j++] = i;
  }

  bins[j] = length;

} /* init_quick_find */
/******************************************************************************/

INT find_index(INT key, INT list[], INT length)

/*******************************************************************************

  Find 'key' in 'list' and return the index number.

  Author:          Ray Tuminaro, SNL, 1422
  =======

  Return code:     INT, -1 = key not found, i = list[i] = key
  ============

  Parameter list:
  ===============

  key:             Element to be search for in list.

  list:            List to be searched.

  length:          Length of list.

*******************************************************************************/

{

  /* local variables */

  INT start, end;
  INT mid;

  /**************************** execution begins ******************************/

  if (length == 0) return -1;

  start = 0;
  end   = length - 1;

  while (end - start > 1) {
    mid = (start + end) / 2;
    if (list[mid] < key) start = mid;
    else end = mid;
  }

  if (list[start] == key) return start;
  if (list[end] == key)   return end;
  return -1;

} /* end of find_index */
#endif
