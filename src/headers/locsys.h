/*----------------------------------------------------------------------*/
/*!
\file
\brief

<pre>
Maintainer: Steffen Genkinger
            genkinger@statik.uni-stuttgart.de
            http://www.uni-stuttgart.de/ibs/members/genkinger/
            0711 - 685-6127
</pre>
*/

/*----------------------------------------------------------------------*/
/* different definitions of local systems for input
 * bborn 06/07 */
typedef enum _LOCSYSTYP
{
                       locsys_none,        /* no locsys                */
                       locsys_basevec,     /* three basevectors        */
                       locsys_line_plane,  /* line + plane             */
                       locsys_line_line,   /* two lines                */
                       locsys_fmc,         /* fluid mass consistent    */
                       locsys_fsnd         /* free surface normal def. */
} LOCSYSTYP;

/*----------------------------------------------------------------------*/
/* type to differentiate if an element has local systems or none
 * bborn 06/07 */
typedef enum _ELECOSYS
{
                       locsys_no,
                       locsys_yes
} ELESOSYS;

/*----------------------------------------------------------------------*/
/* transformation/rotation matrix 
   Transformation of displacements (or other unbound vectors)
     | (Dx*) |   | cos(Xx*)   cos(Yx*)   cos(Zx*) | | (DX) |
     | (Dy*) | = | cos(Xy*)   cos(Yy*)   cos(Zy*) | | (DY) |
     | (Dz*) |   | cos(Xz*)   cos(Yz*)   cos(Zz*) | | (DZ) |
   in which cos(Xx*) is the direction cosine, i.e. the cosine of
   the angle enclosed by X- and x*-unit-vectors
*/
typedef struct _LOCSYS
{
INT                    Id;        /*!< Id                               */
enum _LOCSYSTYP        locsystyp; /*!< type of local co-sys             */
ARRAY                  xloc;      /*!< x - base vector                  */
ARRAY                  yloc;      /*!< y - base vector                  */
ARRAY                  zloc;      /*!< z - base vector                  */
DOUBLE                 lXx;       /*! cos(Xx*)                          */
DOUBLE                 lXy;       /*! cos(Xy*)                          */
DOUBLE                 lXz;       /*! cos(Xz*)                          */
DOUBLE                 lYx;       /*! cos(Yx*)                          */
DOUBLE                 lYy;       /*! cos(Yy*)                          */
DOUBLE                 lYz;       /*! cos(Yz*)                          */
DOUBLE                 lZx;       /*! cos(Zx*)                          */
DOUBLE                 lZy;       /*! cos(Zy*)                          */
DOUBLE                 lZz;       /*! cos(Zz*)                          */
} LOCSYS;

/*----------------------------------------------------------------------*/
/* transformation direction
 * bborn 06/07 */
typedef enum _LOCSYS_TRF_KIND
{
  locsys_trf_XYZ_to_xyz = 0,      /* from XYZ-system to xyz*-system */
  locsys_trf_xyz_to_XYZ           /* from xyz*-system to XYZ-system */
} LOCSYS_TRF_KIND;

