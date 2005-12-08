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
 | curve types                                            m.gee 2/02    |
 | there can be different types of load curves                          |
 *----------------------------------------------------------------------*/
typedef enum _CURVTYP
{
                       curve_none,
                       curve_polygonal,  /* polygonal, piecewise linear curve */
                       curve_explicit    /* excplicit function */
} CURVTYP;
/*----------------------------------------------------------------------*
 | general dynamic-curves                                 m.gee 4/01    |
 | the arrays time and value are dependent on type of load curve and    |
 | size of the curve                                                    |
 *----------------------------------------------------------------------*/
typedef struct _CURVE
{
INT                    Id;            /* Id of the load curve */
enum _CURVTYP          curvetyp;      /* type of load curve */
INT                    bystep;        /* flag whether curve operated by number of steps or in absolut time */
INT                    numex;         /* number of explicit function */
ARRAY                  time;          /* array for time steps */
ARRAY                  value;         /* array for values at time steps */
DOUBLE                 T;             /* forgot about it..... */
DOUBLE                 c1;            /* constant for explicit functions */
DOUBLE                 c2;            /* constant for explicit functions */
} CURVE;


/*!----------------------------------------------------------------------
\brief enum of FUNCT types

<pre>                                                              mn 02/04
This is the enumeration of all types for FUNCT
</pre>

*----------------------------------------------------------------------*/
typedef enum _FUNCTTYP
{
  /* functions on a straight line */
  funct_line_lin,      /* linear function on a line */
  funct_line_quad,     /* quadratic parabola on a line */

  /* functions on a surface, axisymmetric */
  funct_radius_lin,      /* linear function on a circle */
  funct_radius_quad,     /* quadratic parabola on a circle */

  /* special functions */
  funct_bel,           /* beltrami flow function */
  funct_kim,           /* Kim-Moin flow function */
  funct_cyl            /* 3d cylinder */
} FUNCTTYP;


/*!----------------------------------------------------------------------
\brief structure FUNCT

<pre>                                                              mn 02/04
This structure contains the data of one spatial loading function.
</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT
{
  INT                    Id;               /* Id of the function */
  enum _FUNCTTYP         functtyp;         /* type of load function */
  union                                    /* union pointer to funct properties */
  {
    struct _FUNCT_LINE_LIN    *funct_line_lin;   /* linear function */
    struct _FUNCT_LINE_QUAD   *funct_line_quad;  /* quadratic parabola */

    struct _FUNCT_RADIUS_LIN    *funct_radius_lin;   /* linear function */
    struct _FUNCT_RADIUS_QUAD   *funct_radius_quad;  /* quadratic parabola */

    struct _FUNCT_BEL         *funct_bel;        /* Beltrami flow function */
    struct _FUNCT_KIM         *funct_kim;        /* Kim-Moin flow function */
    struct _FUNCT_CYL         *funct_cyl;        /* 3D cylinder */
  }                            typ;              /* name of union */

} FUNCT;



/*!----------------------------------------------------------------------
\brief structure FUNCT_LINE_LIN

<pre>                                                              mn 02/04
This structure contains all information about a linear loading function.

f(xi) = m * xi + b  ,  0 <= xi <= 1

</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_LINE_LIN
{
  DOUBLE               x1[3];             /* first point of line  (xi = 0) */
  DOUBLE               x2[3];             /* second point of line (xi = 1) */
  DOUBLE               val1;              /* function value at first point */
  DOUBLE               val2;              /* function value at second point */

  DOUBLE               length;            /* distance between x1 and x2 */
  DOUBLE               m;                 /* slope of this function */
  DOUBLE               b;                 /* offset of this function */

} FUNCT_LINE_LIN;




/*!----------------------------------------------------------------------
\brief structure FUNCT_LINE_QUAD

<pre>                                                              mn 02/04
This structure contains all information about a quadratic loading function.
</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_LINE_QUAD
{
  DOUBLE               x1[3];             /* first point of line  (xi = 0) */
  DOUBLE               x2[3];             /* second point of line (xi = 1) */

  DOUBLE               length;            /* distance between x1 and x2 */

} FUNCT_LINE_QUAD;





/*!----------------------------------------------------------------------
\brief structure FUNCT_RADIUS_LIN

<pre>                                                              mn 02/04
This structure contains all information about a linear loading function.

f(xi) = m * xi + b  ,  0 <= xi <= 1

</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_RADIUS_LIN
{
  DOUBLE               x1[3];             /* center of the circle =
                                             first point of line  (xi = 0) */
  DOUBLE               x2[3];             /* any point on the circle =
                                             second point of line (xi = 1) */
  DOUBLE               val1;              /* function value at center */
  DOUBLE               val2;              /* function value at outside */

  DOUBLE               length;            /* distance between x1 and x2 */
  DOUBLE               m;                 /* slope of this function */
  DOUBLE               b;                 /* offset of this function */

} FUNCT_RADIUS_LIN;






/*!----------------------------------------------------------------------
\brief structure FUNCT_RADIUS_QUAD

<pre>                                                              mn 02/04
This structure contains all information about a quadratic loading function.
</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_RADIUS_QUAD
{
  DOUBLE               x1[3];             /* center of the circle =
                                             first point of line  (xi = 0) */
  DOUBLE               x2[3];             /* second point of line (xi = 1) */

  DOUBLE               length;            /* distance between x1 and x2 */

} FUNCT_RADIUS_QUAD;



/*!----------------------------------------------------------------------
\brief structure FUNCT_KIM

<pre>                                                              mn 02/04
This structure contains all information about a Kim-Moin loading function.
</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_KIM
{
  INT                 dummy;
} FUNCT_KIM;


/*!----------------------------------------------------------------------
\brief structure FUNCT_BEL

<pre>                                                              mn 02/04
This structure contains all information about a Beltrami loading function.
</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_BEL
{
  INT                 dummy;
} FUNCT_BEL;


/*!----------------------------------------------------------------------
\brief structure FUNCT_CYL

<pre>                                                              mn 02/04
This structure contains all information about a loading function.
</pre>

*----------------------------------------------------------------------*/
typedef struct _FUNCT_CYL
{
  DOUBLE                 um;
} FUNCT_CYL;


