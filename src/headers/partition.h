/*!---------------------------------------------------------------------
\file
\brief domain decomposition and metis structures

<pre>
Maintainer: Malte Neumann
            neumann@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/neumann/
            0711 - 685-6121
</pre>

---------------------------------------------------------------------*/

/*!
\addtogroup PARALLEL
*//*! @{ (documentation module open)*/

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
typedef struct _PAR
{
INT               myrank;                /*!< the individual processor number */
INT               nprocs;                /*!< total number of processors */
} PAR;

/*! @} (documentation module close)*/


