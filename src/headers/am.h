/*!------------------------------------------------------------------------
\file am.h
\brief C-style definitions of arrays
\level 3
\maintainer Martin Kronbichler
------------------------------------------------------------------------*/
#ifndef AM_H
#define AM_H

#include "ccarat_datatypes.h"

/*!
\addtogroup AMSYSTEM
*//*! @{ (documentation module open)*/

/*!------------------------------------------------------------------------
\brief main structure all kinds of fields are kept

m.gee 6/01

main structure all kinds of fields are kept with

-------------------------------------------------------------------------*/
typedef struct _ARRAY
{
char                name[9];           /*!< name of the field (just for fun) */
INT                 fdim;              /*!< first dimension of field         */
INT                 sdim;              /*!< scnd dimension of field          */
enum
   {
    cca_XX,                            /*!< not defined    */
    cca_DA,                            /*!< DOUBLE array   */
    cca_DV,                            /*!< DOUBLE vector  */
    cca_IA,                            /*!< integer array  */
    cca_IV                             /*!< integer vector */
   }                Typ;               /*!< enum type of field */
union
   {
    INT     *iv;                       /*!< integer vector */
    DOUBLE  *dv;                       /*!< DOUBLE vector  */
    INT    **ia;                       /*!< integer array  */
    DOUBLE **da;                       /*!< DOUBLE array   */
   }                a;                 /*!< ptr used for calculations        */
#ifdef DEBUG
struct _TRACEARRAY  *mytracer;         /*!< bugtracing information */
#endif
} ARRAY;

/*!-------------------------------------------------------------------------
\brief main structure all kinds of 3D & 4D fields are kept with

m.gee 12/01

main structure all kinds of 3D & 4D fields are kept with

-------------------------------------------------------------------------*/
typedef struct _ARRAY4D
{
char                name[9];           /*!< name of the field (just for fun) */
INT                 fdim;              /*!< first dimension of field         */
INT                 sdim;              /*!< scnd dimension of field          */
INT                 tdim;              /*!< third dimension of field         */
INT                 fodim;             /*!< fourth dimension of field        */
enum
   {
    cca_XX4D,                              /*!< not defined    */
    cca_D3,                                /*!< DOUBLE 3D-array   */
    cca_D4,                                /*!< DOUBLE 4D-array  */
    cca_I3,                                /*!< integer 3D-array  */
    cca_I4                                 /*!< integer 4D-array  */
   }                Typ;               /*!< enum type of field */
union
   {
     DOUBLE   ***d3;                   /*!< 3D - DOUBLE array */
     DOUBLE  ****d4;                   /*!< 4D - DOUBLE array */
     INT      ***i3;                   /*!< 3D - integer array */
     INT     ****i4;                   /*!< 4D - integer array */
   }                a;                 /*!< name of union */
#ifdef DEBUG
struct _TRACEARRAY  *mytracer;         /*!< bugtracing information */
#endif
} ARRAY4D;
/*! @} (documentation module close)*/

#endif
