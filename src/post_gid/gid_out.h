/*!
\file
\brief Output functions that write binary GiD files.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is a port of the gidpost library version 1.5. Most parts are just
copied.

\author u.kue
\date 09/04

*/


#ifndef GID_OUT_H
#define GID_OUT_H

#include "post_gid.h"
#include "../headers/standardtypes.h"


/* domain dimension */

typedef enum { GiD_2D = 2, GiD_3D = 3 } GiD_Dimension;

/* ---------------------------------------------------------------------------
 *
 *  Post Mesh Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Begin a new mesh. After that you should write the nodes and the
 *  elements. When writing the elements it is assumed a connectivity
 *  of size given in the paramenter NNode.
 */

int GiD_BeginMesh( char * MeshName,
                   GiD_Dimension Dim, GiD_ElementType EType, int NNode );

/*
 *  End current mesh.
 */

int GiD_EndMesh();

/*
 *  Start a coordinate block in the current mesh. All nodes
 *  coordinates must be written between a to this function and
 *  GiD_EndCoordinates().
 */

int GiD_BeginCoordinates();

/*
 *  Close the current coordinate block
 */

int GiD_EndCoordinates();

/*
 *  Write a coordinate member at the current Coordinates Block
 */

int GiD_WriteCoordinates(int id, double x, double y, double z);
int GiD_WriteCoordinates2D(int id, double x, double y);

/*
 *  Start a elements block in the current mesh
 */

int GiD_BeginElements();

/*
 *  Close the current elements block
 */

int GiD_EndElements();

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to GiD_BeginMesh.
 *
 */

int GiD_WriteElement( int id, int nid[] );

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *
 */

int GiD_WriteElementMat( int id, int nid[] );

/* ---------------------------------------------------------------------------
 *
 *  Post Result Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post result file. All subsequent call to functions
 *  write the information into this file. If there is no mesh file
 *  opened then the output of the mesh is writen into this file also.
 *
 */

int GiD_OpenPostResultFile( char * FileName );

/*
 *  Close the current post result file
 */

int GiD_ClosePostResultFile();

/*
 *  Begin Gauss Points definition. The gauss point definition should
 *  have a name wich may be referenced in futher results blocks. The
 *  gauss points could be internal (InternalCoord=1) or given
 *  (InternalCoord=0). If the gauss points are given then the list of
 *  its natural coordinates should be written using the function
 *  GiD_WriteGaussPoint2D or GiD_WriteGaussPoint3D depending on the
 *  dimension of the element type.
 */

int GiD_BeginGaussPoint( char * name, GiD_ElementType EType, char * MeshName,
			 int GP_number, int NodesIncluded, int InternalCoord );

/*
 *  End current Gauss Points definition
 */

int GiD_EndGaussPoint();

/*
 *  Write internal gauss point coordinate.
 */

int GiD_WriteGaussPoint2D( double x, double y );

int GiD_WriteGaussPoint3D( double x, double y, double z );

/*
 *  Begin a Range Table definition. With a range table you can group
 *  the result values into intervals and label each interval with a
 *  name. Inside GiD this can be visualized with a contour range.
 */

int GiD_BeginRangeTable( char * name );

/*
 *  End a Range Table definition
 */

int GiD_EndRangeTable();

/*
 *  Write Range functions --
 *
 *   WriteMinRange : write a range with an implicit minimum value, the
 *   minimum absolute in the result set.
 *
 *   WriteRange : write an explicit range.
 *
 *   WritemaxRange: write a range with an implicit maximum value, the
 *   maximum absolute in the result set.
 */

int GiD_WriteMinRange( double max, char * name );
int GiD_WriteRange( double min, double max, char * name );
int GiD_WriteMaxRange( double min, char * name );

typedef enum {
  GiD_Scalar = 0,
  GiD_Vector,
  GiD_Matrix,
  GiD_PlainDeformationMatrix,
  GiD_MainMatrix,
  GiD_LocalAxes
} GiD_ResultType;

typedef enum { GiD_OnNodes, GiD_OnGaussPoints } GiD_ResultLocation;

/*
 *  Begin Result Block. This function open a result block. All the result
 *  information is provided the function call. If you define a result using this
 *  function you don't need to call the function GiD_ResultValues because it is
 *  done implicitly.
 */

int GiD_BeginResult(char * Result, char * Analysis, double step,
		    GiD_ResultType Type, GiD_ResultLocation Where,
		    char * GaussPointsName, char * RangeTable,
		    int compc, char * compv[]);

/*
 *  Begin Result Block. This function open a result block. Only the result,
 *  analisys, location and location name are provided in the function call. The
 *  other result attributes as range table or components names are provided in a
 *  separated function calls.
 */

int GiD_BeginResultHeader(char * Result, char * Analysis, double step,
			  GiD_ResultType Type, GiD_ResultLocation Where,
			  char * GaussPointsName);

/*
 *  Define the range table associated to the current result, either a single
 *  result block or the current result defined in a result group.
 */

int GiD_ResultRange(char * RangeTable);

/*
 *  Define the components names associated to the current result, either a single
 *  result block or the current result defined in a result group.
 */

int GiD_ResultComponents(int compc, char* compv[]);

/*
 *  Define the unit string associated to the current result, either a single
 *  result block or the current result defined in a result group.
 */

/* This function is not supported yet inside GiD */
/* int GiD_ResultUnit(char * UnitName); */

/*
 *  Begin a result group. All grouped in the same analysis and step. Also the
 *  result location and location name is provided.
 */

int GiD_BeginResultGroup(char * Analysis, double step, GiD_ResultLocation Where,
			 char * GaussPointsName);

/*
 *  Define a result member of a result group given the name and result type.
 */

int GiD_ResultDescription(char * Result, GiD_ResultType Type);

/*
 *  Mark the starting point for writing the values of the current result either
 *  single or group. This function is not needed if the current result was
 *  defined with a call to GiD_BeginResult.
 */

int GiD_ResultValues();

/*
 *  Close a previously opened result either single or group.
 */

int GiD_EndResult();

/*
 *  Flushes all pending output into the postprocess file. This function should
 *  be called only when strictly necessary when writing in GiD_PostAsciiZipped *
 *  or GiD_PostBinary modes because it can degrade compression.
 */

int GiD_FlushPostFile();

/*
 *  Write result functions
 */

int GiD_WriteScalar( int id, double v );
int GiD_WriteVector( int id, double x, double y, double z );
int GiD_WriteVectorModule( int id, double x, double y, double z, double mod );
int GiD_Write2DMatrix( int id, double Sxx, double Syy, double Sxy );
int GiD_Write3DMatrix( int id, double Sxx, double Syy, double Szz,
		       double Sxy, double Syz, double Sxz);
int GiD_WritePlainDefMatrix( int id, double Sxx, double Syy, double Sxy, double Szz );
int GiD_WriteMainMatrix( int id,
			 double Si, double Sii, double Siii,
			 double Vix, double Viy, double Viz,
			 double Viix, double Viiy, double Viiz,
			 double Viiix, double Viiiy, double Viiiz );
int GiD_WriteLocalAxes( int id, double euler_1, double euler_2, double euler_3 );


#endif
