/*!
\file
\brief Output functions that write text GiD files.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

This is a port of the gidpost library version 1.5.

\author u.kue
\date 09/04

*/

#include <assert.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "gid_out.h"

#define LINE_SIZE 8192

/* although the following function declaration is contained in <stdio.h>
 * it is not available when compiling with compile flag -ansi as we  normally do.
 * To remove the compiler warning, the function declaration was added explicitely.
 * gjb 11/11
 */
int snprintf (char *s, size_t size, const char *template, ...);

/*----------------------------------------------------------------------*/
/*!
  \brief A pretty simple internal structure.

  This has been the CPostAcsii class in gidpost.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/

typedef struct _GID_POST_ASCII {
  FILE* file;
  int fail;
  int LastID;
  int connectivity;
} GID_POST_ASCII;


void init_gid_post_acsii(GID_POST_ASCII* gid)
{
  gid->file = NULL;
  gid->fail = 0;
  gid->LastID = -1;
  gid->connectivity = 0;
}

int gid_close(GID_POST_ASCII* gid);

void destroy_gid_post_acsii(GID_POST_ASCII* gid)
{
  gid_close(gid);
}


int gid_match_connectivity(GID_POST_ASCII* gid, int written)
{
  if (written == gid->connectivity) {
    return 0;
  }
  else if ((written-1) == gid->connectivity) {
    return 1;
  }
  else {
    return 2;
  }
}


int gid_open(GID_POST_ASCII* gid, char* name)
{
  gid->file = fopen(name, "wb");
  return gid->file == NULL;
}

int gid_close(GID_POST_ASCII* gid)
{
  if (gid->file) {
    gid->fail = fclose(gid->file);
    gid->file = NULL;
  }
  else
    gid->fail = 1;
  return gid->fail;
}

int gid_flush(GID_POST_ASCII* gid)
{
  return gid->file ? fflush(gid->file) : 1;
}

int gid_write_string(GID_POST_ASCII* gid, char* str)
{
  fprintf(gid->file, "%s\n", str );
  return 0;
}

int gid_begin_coordinates(GID_POST_ASCII* gid)
{
  return gid_write_string(gid, "Coordinates");
}

int gid_begin_elements(GID_POST_ASCII* gid)
{
  return gid_write_string(gid, "Elements");
}

int gid_begin_values(GID_POST_ASCII* gid)
{
  return gid_write_string(gid, "Values");
}

int gid_end_values(GID_POST_ASCII* gid)
{
  return gid_write_string(gid, "End Values");
}

int gid_write_values(GID_POST_ASCII* gid, int id, int n, ... )
{
  va_list ap;
  int i;
  double value;

  if (gid->LastID != id) {
    fprintf(gid->file, "%d", id);
    gid->LastID = id;
  }
  va_start(ap, n);
  for ( i = 0; i < n; i++ ) {
    value = va_arg(ap, double);
    fprintf(gid->file, " %12.5E", value);
  }
  fprintf(gid->file, "\n");
  va_end(ap);
  return 0;
}

int gid_write_buffer(GID_POST_ASCII* gid, int id, int n, double* buffer)
{
  int i;

  if ( gid->LastID != id ) {
    fprintf(gid->file, "%d", id);
    gid->LastID = id;
  }
  for ( i = 0; i < n; i++ ) {
    fprintf(gid->file, " %12.5E", buffer[i]);
  }
  return 0;
}

int gid_write2d(GID_POST_ASCII* gid, double x, double y)
{
  char line[256];

  sprintf( line, "%12.5E %12.5E", x, y );
  return gid_write_string(gid, line);
}

int gid_write3d(GID_POST_ASCII* gid, double x, double y, double z)
{
  char line[256];

  sprintf( line, "%12.5E %12.5E %12.5E", x, y, z );
  return gid_write_string(gid, line);
}

int gid_write_element(GID_POST_ASCII* gid, int id, int n, int nid[])
{
  int i;

  fprintf(gid->file, "%d", id);
  for ( i = 0; i < n; i++ ) {
    fprintf(gid->file, " %d", nid[i]);
  }
  fprintf(gid->file, "\n");
  return 0;
}

int gid_write_post_header(GID_POST_ASCII* gid)
{
  return gid_write_string(gid, "GiD Post Results File 1.0");
}

/*----------------------------------------------------------------------*/
/*!
  \brief Another support structure.

  This has been the CBufferValues class in gidpost.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/

typedef struct _GID_BUFFER {
  double * buffer;
  int last;
  int values_size;
  int buffer_size;
  GID_POST_ASCII* file;
} GID_BUFFER;

#define buffer_is_empty(buf) (!((buf)->buffer) || ((buf)->last==-1))


void init_gid_buffer(GID_BUFFER* buf)
{
  buf->file = NULL;
  buf->buffer = NULL;
  buf->last = -1;
  buf->values_size = 0;
  buf->buffer_size = 0;
}

void destroy_gid_buffer(GID_BUFFER* buf)
{
  if (buf->buffer) {
    free(buf->buffer);
    buf->buffer = NULL;
  }
  buf->last = -1;
  buf->values_size = 0;
  buf->buffer_size = 0;
}

void gid_buffer_init(GID_BUFFER* buf, int sz)
{
  if (!buf->buffer) {
    buf->buffer = (double*)malloc((buf->buffer_size=sz)*sizeof(double));
  }
  else if (sz > buf->buffer_size) {
    buf->buffer =
      (double*)realloc(buf->buffer,
                       (buf->buffer_size=sz)*sizeof(double));
  }
  buf->last = -1;
  buf->values_size = sz;
}

int gid_buffer_write_values(GID_BUFFER* buf, int id, int n, ...)
{
  va_list ap;
  int i;

  assert(buf->last+n < buf->values_size);

  va_start(ap, n);
  for (i = 0; i < n; i++) {
    buf->buffer[++buf->last] = va_arg(ap, double);
  }
  va_end(ap);
  if (buf->last == buf->values_size-1) {
    assert(buf->file);

    if (gid_write_values(buf->file, id, buf->values_size, buf->buffer)) {
      /* could not write buffer values */
      return 1;
    }
    /* prepare for next group of values */
    buf->last = -1;
  }
  return 0;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Higher level functions that make up the userinterface.

  This has been the CPostAcsii class in gidpost.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/


static GID_POST_ASCII MeshFile;
static GID_POST_ASCII ResultFile;


typedef enum {
  POST_UNDEFINED,
  POST_S0,           /* TOP level */
  POST_MESH_S0,      /* MESH header */
  POST_MESH_COORD0,  /* inside a Coordinate block */
  POST_MESH_COORD1,  /* after a Coordinate block but inside a MESH */
  POST_MESH_ELEM,    /* inside an Element block */
  POST_GAUSS_S0,     /* GAUSS point block: implicit */
  POST_GAUSS_GIVEN,  /* GAUSS point block: explicit */
  POST_RANGE_S0,     /* RANGE table block */
  POST_RESULT_SINGLE,    /* Result block */
  POST_RESULT_GROUP, /* Result group block */
  POST_RESULT_DESC,  /* Result description block */
  POST_RESULT_VALUES /* writing values */
} post_state;

static post_state level_mesh = POST_UNDEFINED;
static post_state level_res  = POST_UNDEFINED;

static const char * level_desc[] = {
  "UNDEFINED level",
  "TOP level",
  "MESH header",
  "inside a Coordinate block",
  "after a Coordinate block but inside a MESH",
  "inside an Element block",
  "GAUSS point block: implicit",
  "GAUSS point block: explicit",
  "RANGE table block",
  "Result block",
  "Result group block",
  "Result description block",
  "writing values",
  "unknown"
};

static const char * GetStateDesc(post_state i)
{
  static int last = sizeof(level_desc)/sizeof(level_desc[0]) - 1;
  return (i >= last) ? level_desc[last] : level_desc[i];
}

static int GP_number_check = 0;   /* number of gauss points to be written */
static int gauss_written = 0;     /* number of gauss points written */
static int values_location = 0;   /* number of values per ID to be written */
static int flag_isgroup = 0;      /* is this a result group */
static int flag_begin_values = 0; /* Is values written ? */

static GID_BUFFER buffer_values;


static const char * strElementType[]= {
  "",
  "Point",
  "Linear",
  "Triangle",
  "Quadrilateral",
  "Tetrahedra",
  "Hexahedra"
};

const char * GetElementTypeName( GiD_ElementType type )
{
  return strElementType[type];
}

typedef struct {
  const char * str;
  int values;
} SResultTypeInfo;

static SResultTypeInfo _ResultTypeInfo[] = {
  {"Scalar", 1},
  {"Vector", 3},
  {"Matrix", 6},
  {"PlainDeformationMatrix", 4},
  {"MainMatrix", 12},
  {"LocalAxes", 3}
};

static const char * GetResultTypeName(GiD_ResultType type)
{
  return _ResultTypeInfo[type].str;
}

static int GetResultTypeValues(GiD_ResultType type)
{
  return _ResultTypeInfo[type].values;
}

static int CheckState(post_state s_req, post_state s_cur)
{
  if (s_req!=s_cur) {
    printf("invalid state '%s' should be '%s'\n", GetStateDesc(s_cur), GetStateDesc(s_req));
    return 0;
  }
  return 1;
}

static int ValidateConnectivity(GiD_ElementType etype , int NNode)
{
  int error;

  switch (etype) {
  case GiD_Point:
    error = (NNode != 1);
    break;
  case GiD_Linear:
    error = (NNode != 2 && NNode != 3);
    break;
  case GiD_Triangle:
    error = (NNode != 3 && NNode != 6);
    break;
  case GiD_Quadrilateral:
    error = (NNode != 4 && NNode != 8 && NNode != 9);
    break;
  case GiD_Tetrahedra:
    error = (NNode != 4 && NNode != 10);
    break;
  case GiD_Hexahedra:
    error = (NNode != 8 && NNode != 20 && NNode != 27);
    break;
   default:
    printf("invalid type of element code %d", etype);
    return 0;
  }
  if (error) {
    printf("invalid number of nodes '%d' for element type %s\n",
	   NNode, GetElementTypeName(etype));
    return 0;
  }
  return 1;
}

/* ---------------------------------------------------------------------------
 *
 *  Post Mesh Interface
 *
 * ---------------------------------------------------------------------------
 */

static int GiD_ClosePostMeshFile();

/*
 *  Open a new post mesh file
 */

static int GiD_OpenPostMeshFile( char * FileName )
{
  level_mesh = POST_UNDEFINED;
  /* now outputMesh points to MeshFile */
  if (gid_open(&MeshFile, FileName)) {
    /* Open failed */
    return 4;
  }
  if (gid_write_post_header(&MeshFile)) {
    /* WritePostHeader failed */
    GiD_ClosePostMeshFile();
    return 5;
  }
  level_mesh = POST_S0;
  return 0;
}

/*
 *  Close the current post mesh file
 */

static int GiD_ClosePostMeshFile()
{
  assert(CheckState(POST_S0, level_mesh));

  destroy_gid_post_acsii(&MeshFile);
  level_mesh = POST_UNDEFINED;
  return 0;
}

/*
 *  Begin a new mesh --
 *
 *    Start a mesh block
 */

int GiD_BeginMesh(const char* MeshName,
                  GiD_Dimension Dim,
                  GiD_ElementType EType,
                  int NNode)
{
  GID_POST_ASCII * mesh;
  int fail = 1;
  char line[LINE_SIZE];

  assert(CheckState(POST_S0,level_mesh));
  /* here we sould validate EType & NNode */
  assert(ValidateConnectivity(EType,NNode));

  mesh = &MeshFile;

  snprintf(line, LINE_SIZE-1,
           "MESH \"%s\" dimension %d ElemType %s Nnode %d",
           MeshName, Dim, GetElementTypeName(EType), NNode);
  if ( !(fail = gid_write_string(mesh, line)) ) {
    mesh->connectivity = NNode;
  }
  level_mesh = POST_MESH_S0;
  return fail;
}

/*
 *  End current mesh
 */

int GiD_EndMesh()
{
  /* check & change state */

  /*
  assert(CheckState(POST_MESH_ELEM,level_mesh));
  level_mesh = POST_S0;
  */

  fflush(MeshFile.file);

  return 0;
}

/*
 *  Start a coordinate block in the current mesh
 */


int GiD_BeginCoordinates()
{
  /* state checking */
  assert(CheckState(POST_MESH_S0, level_mesh));
  level_mesh = POST_MESH_COORD0;
  return gid_begin_coordinates(&MeshFile);
}

/*
 *  Close the current coordinate block
 */

int GiD_EndCoordinates()
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, level_mesh));
  level_mesh = POST_MESH_COORD1;
  return gid_write_string(&MeshFile, "End Coordinates");
}

/*
 *  Write a coordinate member at the current Coordinates Block
 */

int GiD_WriteCoordinates( int id, double x, double y, double z )
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, level_mesh));
  /* keep in the same level */
  return gid_write_values(&MeshFile, id, 3, x, y, z);
}

int GiD_WriteCoordinates2D(int id, double x, double y)
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD0, level_mesh));
  /* keep in the same level */
  return gid_write_values(&MeshFile, id, 3, x, y, 0.0);
}

/*
 *  Start a elements block in the current mesh
 */

int GiD_BeginElements()
{
  /* state checking */
  assert(CheckState(POST_MESH_COORD1, level_mesh));
  level_mesh = POST_MESH_ELEM;
  return gid_begin_elements(&MeshFile);
}

/*
 *  Close the current elements block
 */

int GiD_EndElements()
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));
  level_mesh = POST_S0;
  return gid_write_string(&MeshFile, "End Elements");
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh.
 *
 */

int GiD_WriteElement( int id, int nid[] )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));
  /* keep on the same state */
  return gid_write_element(&MeshFile, id, MeshFile.connectivity, nid);
}

/*
 *  Write an element member at the current Elements Block. The rest of
 *  the arguments define the connectivity of the element and its
 *  number depends on the NNode parameter given in the previous call
 *  to BeginPostMesh. The last id correspond to a material number.
 *
 */

int GiD_WriteElementMat( int id, int nid[] )
{
  /* state checking */
  assert(CheckState(POST_MESH_ELEM, level_mesh));
  /* keep on the same state */
  return gid_write_element(&MeshFile, id, MeshFile.connectivity+1, nid);
}

/* ---------------------------------------------------------------------------
 *
 *  Post Result Interface
 *
 * ---------------------------------------------------------------------------
 */

/*
 *  Open a new post result file
 */

int GiD_OpenPostResultFile( char * FileName )
{
  level_res = POST_UNDEFINED;
  level_mesh = POST_UNDEFINED;
  init_gid_post_acsii(&ResultFile);
  if (gid_open(&ResultFile, FileName)) {
    /* could not open file */
    return 4;
  }
  if (gid_write_post_header(&ResultFile)) {
    /* WritePostHeader failed */
    GiD_ClosePostResultFile();
    return 5;
  }
  level_res = POST_S0;
  level_mesh = POST_S0;
  buffer_values.file = &ResultFile;

  /* this is an ugly hack */
  /* It's needed to have one interface (set of functions) for both
   * text and binary files. */
  {
    int i;
    i = strlen(FileName);
    assert(i>4);
    assert(FileName[i-4] == '.');
    FileName[i-3] = 'm';
    FileName[i-2] = 's';
    FileName[i-1] = 'h';
  }
  return GiD_OpenPostMeshFile(FileName);
}

/*
 *  Close the current post result file
 */

int GiD_ClosePostResultFile()
{
  int fail = 1;

  GiD_ClosePostMeshFile();

  assert(CheckState(POST_S0, level_res));

  destroy_gid_post_acsii(&ResultFile);

  fail = ResultFile.fail;
  level_res = POST_UNDEFINED;
  level_mesh = POST_UNDEFINED;

  return fail;
}

/*
 *  Begin Gauss Points definition
 */

int GiD_BeginGaussPoint( const char * name, GiD_ElementType EType, const char * MeshName,
			 int GP_number, int NodesIncluded, int InternalCoord )
{
  char line[LINE_SIZE];

  /* check state & validation */
  assert(CheckState(POST_S0, level_res));

  snprintf( line, LINE_SIZE-1,
	    "GaussPoints \"%s\" ElemType %s", name, GetElementTypeName(EType));
  if ( MeshName && *MeshName ) {
    strcat(line, " \"");
    strcat(line, MeshName);
    strcat(line, "\"");
  }
  if ( gid_write_string(&ResultFile, line) )
    return 1;
  snprintf(line, LINE_SIZE, "Number Of Gauss Points: %d", GP_number);
  if ( gid_write_string(&ResultFile, line) )
    return 1;
  /* here we could save the number of GP in order to check at
     EndGaussPoint */
  GP_number_check = GP_number;
  if ( EType == GiD_Linear ) {
    if ( NodesIncluded ) {
      if ( gid_write_string(&ResultFile, "  Nodes included") )
	return 1;
    }
    else
      if ( gid_write_string(&ResultFile, "  Nodes not included") )
        return 1;
  }
  if (InternalCoord) {
    if (gid_write_string(&ResultFile, "Natural Coordinates: Internal"))
      return 1;
    level_res = POST_GAUSS_S0;
  }
  else {
    if (gid_write_string(&ResultFile, "Natural Coordinates: Given"))
      return 1;
    level_res = POST_GAUSS_GIVEN;
    /* here we can save the size of the coordinates to check later
       in WriteGaussPointXX*/
  }
  return 0;
}

/*
 *  End current Gauss Points definition
 */

static int CheckGaussPointEnd()
{
  if (level_res != POST_GAUSS_S0 && level_res != POST_GAUSS_GIVEN) {
    printf("Invalid call of GiD_EndGaussPoint. Current state is '%s' and should be '%s' or '%s'\n",
	   GetStateDesc(level_res), GetStateDesc(POST_GAUSS_S0), GetStateDesc(POST_GAUSS_GIVEN));
    return 0;
  }
  return 1;
}

static int CheckGaussPointGiven()
{
  if (gauss_written !=  GP_number_check) {
    printf("missmatch in gauss point given, written %d and %d were requiered",
	   gauss_written, GP_number_check);
    return 0;
  }
  return 1;
}

int GiD_EndGaussPoint()
{
  /* check state */
  assert(CheckGaussPointEnd());
#ifndef NDEBUG
  if (level_res == POST_GAUSS_GIVEN)
    assert(CheckGaussPointGiven());
#endif
  level_res = POST_S0;
  GP_number_check = gauss_written = 0;
  return gid_write_string(&ResultFile, "End GaussPoints");
}

/*
 *  Write internal gauss point coordinate.
 */

int GiD_WriteGaussPoint2D( double x, double y )
{
  /* check state */
  assert(CheckState(POST_GAUSS_GIVEN, level_res));
  if (!gid_write2d(&ResultFile, x, y)) {
    ++gauss_written;
    return 0;
  }
  return 1;
}

int GiD_WriteGaussPoint3D( double x, double y, double z )
{
  /* check state */
  assert(CheckState(POST_GAUSS_GIVEN, level_res));
  if (!gid_write3d(&ResultFile, x, y, z)) {
    ++gauss_written;
    return 0;
  }
  return 1;
}

/*
 *  Begin a Range Table definition
 */

int GiD_BeginRangeTable( char * name )
{
  char line[LINE_SIZE];
  /* check & update state */
  assert(CheckState(POST_S0, level_res));

  snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", name);
  if (!gid_write_string(&ResultFile, line)) {
    level_res = POST_RANGE_S0;
    return 0;
  }
  return 1;
}

/*
 *  End a Range Table definition
 */

int GiD_EndRangeTable()
{
  /* check & update state */
  assert(CheckState(POST_RANGE_S0, level_res));

  if (!gid_write_string(&ResultFile, "End ResultRangesTable")) {
    level_res = POST_S0;
    return 0;
  }
  return 1;
}

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

int GiD_WriteMinRange( double max, char * name )
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckState(POST_RANGE_S0, level_res));

  snprintf(line, LINE_SIZE-1, " - %12.5E : \"%s\"", max, name);
  return gid_write_string(&ResultFile, line);
}

int GiD_WriteRange( double min, double max, char * name )
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckState(POST_RANGE_S0, level_res));

  snprintf(line, LINE_SIZE-1, " %12.5E - %12.5E : \"%s\"", min, max, name);
  return gid_write_string(&ResultFile, line);
}

int GiD_WriteMaxRange( double min, char * name )
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckState(POST_RANGE_S0, level_res));

  snprintf(line, LINE_SIZE-1, "%12.5E - : \"%s\"", min, name);
  return gid_write_string(&ResultFile, line);
}

/*
 *  Begin Result Block
 */

int GiD_BeginResult(const char * Result, const char * Analysis, double step,
		    GiD_ResultType Type, GiD_ResultLocation Where,
		    const char * GaussPointsName, const char * RangeTable,
		    int compc, const char * compv[])
{
  char line[LINE_SIZE];
  const char * loc;
  int i;

  /* check & change state */
  assert(CheckState(POST_S0, level_res));

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  snprintf(line, LINE_SIZE-1, "Result \"%s\" \"%s\" %12.5E %s %s",
	   Result, Analysis, step, GetResultTypeName(Type), loc);
  if (Where == GiD_OnGaussPoints) {
    strcat(line, " \"");
    strcat(line, GaussPointsName);
    strcat(line, "\"");
  }
  if (gid_write_string(&ResultFile, line))
    return 1;
  if (RangeTable) {
    snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", RangeTable);
    if  (gid_write_string(&ResultFile, line))
      return 1;
  }
  if ( compc > 0 ) {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for ( i = 0; i < compc; i++ ) {
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, compv[i]);
      strcat(line, "\"");
    }
    if (gid_write_string(&ResultFile, line))
      return 1;
  }
  flag_isgroup = 0;
  flag_begin_values = 1;
  level_res = POST_RESULT_VALUES;
  return gid_begin_values(&ResultFile);
}

int GiD_BeginResultHeader(const char * Result, const char * Analysis, double step,
			  GiD_ResultType Type, GiD_ResultLocation Where,
			  const char * GaussPointsName)
{
  char line[LINE_SIZE];
  const char * loc;

  /* check & change state */
  assert(CheckState(POST_S0, level_res));
  assert(Result);
  assert(Analysis);

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  snprintf(line, LINE_SIZE-1, "Result \"%s\" \"%s\" %12.5E %s %s",
	   Result, Analysis, step, GetResultTypeName(Type), loc);
  if (Where == GiD_OnGaussPoints) {
    strcat(line, " \"");
    assert(GaussPointsName);
    strcat(line, GaussPointsName);
    strcat(line, "\"");
  }
  if (gid_write_string(&ResultFile, line)) {
    /* could not write result header */
    return 1;
  }
  level_res = POST_RESULT_SINGLE;
  flag_isgroup = 0;
  flag_begin_values = 0;
  return 0;
}

static int CheckResultHeaderState()
{
  if (level_res == POST_RESULT_SINGLE || level_res == POST_RESULT_DESC)
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc(level_res),
	 GetStateDesc(POST_RESULT_SINGLE),
	 GetStateDesc(POST_RESULT_DESC));
  return 0;
}

int GiD_ResultRange(char * RangeTable)
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckResultHeaderState());
  assert(RangeTable);

  if (RangeTable) {
    snprintf(line, LINE_SIZE-1, "ResultRangesTable \"%s\"", RangeTable);
    return gid_write_string(&ResultFile, line);
  }
  return 1;
}

int GiD_ResultComponents(int compc, char* compv[])
{
  char line[LINE_SIZE];
  int i;

  /* check state */
  assert(CheckResultHeaderState());
  assert(compc>0);

  if (compc > 0) {
    snprintf(line, LINE_SIZE-1, "ComponentNames");
    for (i = 0; i < compc; i++) {
      strcat(line, " ");
      strcat(line, "\"");
      strcat(line, compv[i]);
      strcat(line, "\"");
    }
    return gid_write_string(&ResultFile, line);
  }
  return 1;
}

int GiD_ResultUnit(char * UnitName)
{
  char line[LINE_SIZE];

  /* check state */
  assert(CheckResultHeaderState());
  assert(UnitName);

  if (UnitName) {
    snprintf(line, LINE_SIZE-1, "UnitName \"");
    strcat(line, UnitName);
    strcat(line, "\"");
    return gid_write_string(&ResultFile, line);
  }
  return 1;
}

int GiD_BeginResultGroup(const char * Analysis, double step, GiD_ResultLocation Where,
			 const char * GaussPointsName)
{
  char line[LINE_SIZE];
  const char * loc;

  /* at this moment result group could only be writen in ASCII format! */
  assert(0);

  /* check & change state */
  assert(CheckState(POST_S0, level_res));
  assert(Analysis);

  loc = (Where == GiD_OnNodes) ? "OnNodes" : "OnGaussPoints";
  snprintf(line, LINE_SIZE-1, "ResultGroup \"%s\" %12.5E %s", Analysis, step, loc);
  if (Where == GiD_OnGaussPoints) {
    strcat(line, " \"");
    assert(GaussPointsName);
    strcat(line, GaussPointsName);
    strcat(line, "\"");
  }
  if (gid_write_string(&ResultFile, line)) {
    /* could not write result header */
    return 1;
  }
  level_res = POST_RESULT_GROUP;
  /* initialize values counter */
  values_location = 0;
  flag_isgroup = 1;
  return 0;
}

static int CheckStateDesc()
{
  if (level_res == POST_RESULT_GROUP || level_res == POST_RESULT_DESC)
    return 1;
  printf("Invalid result state '%s'. Should be: '%s' or '%s'\n",
	 GetStateDesc(level_res),
	 GetStateDesc(POST_RESULT_GROUP),
	 GetStateDesc(POST_RESULT_DESC));
  return 0;
}

int GiD_ResultDescription(const char * Result, GiD_ResultType Type)
{
  char line[LINE_SIZE];

  /* check & change state */
  assert(CheckStateDesc());
  assert(Result);

  snprintf(line, LINE_SIZE-1, "ResultDescription \"%s\" %s", Result, GetResultTypeName(Type));
  if (gid_write_string(&ResultFile, line)) {
    /* could not write result description */
    return 1;
  }
  level_res = POST_RESULT_DESC;
  /* update number of values to check per location */
  values_location += GetResultTypeValues(Type);
  return 0;
}

int GiD_ResultValues()
{
  /* check & change state */
  assert(CheckResultHeaderState());
  assert(!flag_begin_values);

  if (!gid_begin_values(&ResultFile)) {
    level_res = POST_RESULT_VALUES;
    if (flag_isgroup) {
      assert(values_location>0);
      gid_buffer_init(&buffer_values, values_location);
    }
    flag_begin_values = 1;
    return 0;
  }
  return 1;
}

/*
 *  End Result Block
 */

int GiD_EndResult()
{
  int _fail = 0;

  /* check & change state */
  assert(CheckState(POST_RESULT_VALUES, level_res));
  assert(buffer_is_empty(&buffer_values));

  _fail = gid_end_values(&ResultFile);
  fflush(ResultFile.file);
  ResultFile.LastID = -1;
  level_res = POST_S0;
  return _fail;
}

/*
 * Flushes all pending output into the compressed file.
 */

int GiD_FlushPostFile()
{
  return gid_flush(&ResultFile);
}

/*
 *  Write result functions
 */

int GiD_WriteScalar( int id, double v )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  if (flag_isgroup)
    return gid_buffer_write_values(&buffer_values, id, 1, v);
  else
    return gid_write_values(&ResultFile, id, 1, v);
}

int GiD_WriteVector( int id, double x, double y, double z )
{
  double mod;

  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  /* compute the module */
  mod = sqrt(x*x + y*y + z*z);

  return flag_isgroup ? gid_buffer_write_values(&buffer_values, id, 4, x, y, z, mod) :
    gid_write_values(&ResultFile, id, 4, x, y, z, mod);
}

int GiD_WriteVectorModule( int id, double x, double y, double z, double mod )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ? gid_buffer_write_values(&buffer_values, id, 4, x, y, z, mod) :
    gid_write_values(&ResultFile, id, 4, x, y, z, mod );
}

int GiD_Write3DMatrix( int id,
                       double Sxx, double Syy, double Szz,
                       double Sxy, double Syz, double Sxz )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ?
    gid_buffer_write_values(&buffer_values, id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz) :
    gid_write_values(&ResultFile, id, 6, Sxx, Syy, Szz, Sxy, Syz, Sxz);
}

int GiD_Write2DMatrix( int id,
                       double Sxx, double Syy, double Sxy )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ?
    gid_buffer_write_values(&buffer_values, id, 3, Sxx, Syy, Sxy) :
    gid_write_values(&ResultFile, id, 3, Sxx, Syy, Sxy);
}

int GiD_WritePlainDefMatrix( int id, double Sxx, double Syy, double Sxy, double Szz )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ?
    gid_buffer_write_values(&buffer_values, id, 4, Sxx, Syy, Sxy, Szz) :
    gid_write_values(&ResultFile, id, 4, Sxx, Syy, Sxy, Szz);
}

int GiD_WriteMainMatrix( int id,
                         double Si, double Sii, double Siii,
                         double Vix, double Viy, double Viz,
                         double Viix, double Viiy, double Viiz,
                         double Viiix, double Viiiy, double Viiiz )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ?
    gid_buffer_write_values(&buffer_values, id, 12, Si, Sii, Siii,
                            Vix, Viy, Viz,
                            Viix, Viiy, Viiz,
                            Viiix, Viiiy, Viiiz) :
    gid_write_values(&ResultFile, id, 12, Si, Sii, Siii,
                     Vix, Viy, Viz,
                     Viix, Viiy, Viiz,
                     Viiix, Viiiy, Viiiz);
}

int GiD_WriteLocalAxes( int id, double euler_1, double euler_2, double euler_3 )
{
  /* check state */
  assert(CheckState(POST_RESULT_VALUES, level_res));

  return flag_isgroup ?
    gid_buffer_write_values(&buffer_values, id, 3, euler_1, euler_2, euler_3) :
    gid_write_values(&ResultFile, id, 3, euler_1, euler_2, euler_3);
}
