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
/*----------------------------------------------------------------------*
 |                                                          al 08/02    |
 | general eigensolution-variables                                      |
 *----------------------------------------------------------------------*/
typedef struct _ALLEIG                 
{
EIG_SOL_TYPE       soltyp;    /* eigenvalue solution process                        */
INT                sturm ;    /* perform sturm sequence check  0=no                 */
INT                subtyp;    /* subspace: 1 = jacobi rotation 2 = qz-algorithm     */
INT                ifsh;      /* flag for shift                                     */
INT                ilmp;      /* use lumped mass matrix (0=consistent)              */
INT                range;     /* search in a range                                  */
INT                numvec;    /* number of iteration vectors                        */
INT                nroot;     /* number of requested eigenvalues to be converged    */
INT                itemax;    /* maximum number of iterations                       */
INT                ifctr;     /* 1 = output to console 2 = print iteration on file  */ 
DOUBLE             toleig;    /* tolerance to be used in convergence check          */
DOUBLE             shift;     /* relative shift (removed after iteration)           */
DOUBLE             boulo;     /* lower boundary                                     */
DOUBLE             bouup;     /* upper boundary                                     */
} ALLEIG;

