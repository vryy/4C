/*!
\file
\brief A very simple symbol table implementation.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

*/

#ifndef PSS_TABLE_H
#define PSS_TABLE_H

#include "../headers/standardtypes.h"


typedef enum _SYMBOL_TYPES {
  sym_undefined,
  sym_string,
  sym_int,
  sym_real,
  sym_map
} SYMBOL_TYPES;


typedef struct _SYMBOL {
  SYMBOL_TYPES type;
  union {
    CHAR* string;
    INT integer;
    DOUBLE real;
    struct _MAP* dir;
  } s;
  struct _SYMBOL* next;
} SYMBOL;


typedef struct _MAP_NODE {
  CHAR* key;
  SYMBOL* symbol;
  INT count;
  struct _MAP_NODE* lhs;
  struct _MAP_NODE* rhs;
} MAP_NODE;


typedef struct _MAP {
  MAP_NODE root;
  INT count;
} MAP;


/*----------------------------------------------------------------------*/
/*!
  \brief Bring a map variable up to a clean state.

  That's needed before anything can be done with a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_map(MAP* map);


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_map(MAP* map);


/* Find the first symbol with the given key. Use this if you have to
 * travel all symbols with that key. */
SYMBOL* map_find_symbol(MAP* map, CHAR* key);


/* Find the last symbols value. The value has to be of the given
 * type. Returns false on failture. */
INT map_find_string(MAP* map, CHAR* key, CHAR** string);
INT map_find_int(MAP* map, CHAR* key, INT* integer);
INT map_find_real(MAP* map, CHAR* key, DOUBLE* real);
INT map_find_map(MAP* map, CHAR* key, MAP** dir);


/* Find the last symbols value. The value has to be of the given
 * type. Calls dserror on failture. */
CHAR* map_read_string(MAP* map, CHAR* key);
INT map_read_int(MAP* map, CHAR* key);
DOUBLE map_read_real(MAP* map, CHAR* key);
MAP* map_read_map(MAP* map, CHAR* key);


/* Tell whether there is a symbol with given key and value. Only the
 * last symbol with that key is checked. */
INT map_has_string(MAP* map, CHAR* key, CHAR* value);
INT map_has_int(MAP* map, CHAR* key, INT value);
INT map_has_real(MAP* map, CHAR* key, DOUBLE value);
INT map_has_map(MAP* map, CHAR* key);


/* Insert a new symbol. */
void map_insert_string(MAP* map, CHAR* string, CHAR* key);
void map_insert_int(MAP* map, INT integer, CHAR* key);
void map_insert_real(MAP* map, DOUBLE real, CHAR* key);
void map_insert_map(MAP* map, MAP* dir, CHAR* key);


/* Insert a new symbol. Make a copy of all strings. */
void map_insert_string_cpy(MAP* map, CHAR* string, CHAR* key);
void map_insert_int_cpy(MAP* map, INT integer, CHAR* key);
void map_insert_real_cpy(MAP* map, DOUBLE real, CHAR* key);
void map_insert_map_cpy(MAP* map, MAP* dir, CHAR* key);


/* Tell the number of symbols under this key. */
INT map_symbol_count(MAP* map, CHAR* key);


/* Tell whether this symbol has the given type. */
INT symbol_is_string(SYMBOL* symbol);
INT symbol_is_int(SYMBOL* symbol);
INT symbol_is_real(SYMBOL* symbol);
INT symbol_is_map(SYMBOL* symbol);


/* Extract the value of this symbol. Returns false on failture. */
INT symbol_get_string(SYMBOL* symbol, CHAR** string);
INT symbol_get_int(SYMBOL* symbol, INT* integer);
INT symbol_get_real(SYMBOL* symbol, DOUBLE* real);
INT symbol_get_map(SYMBOL* symbol, MAP** map);


/* Extract the value of this symbol. Call dserror on failture. */
CHAR* symbol_string(SYMBOL* symbol);
INT symbol_int(SYMBOL* symbol);
DOUBLE symbol_real(SYMBOL* symbol);
MAP* symbol_map(SYMBOL* symbol);


/* Write the map to file f in a readable way. */
void map_print(FILE* f, MAP* map, INT indent);


/* Read the control file given by name. Put its contents into the map. */
void parse_control_file(MAP* map, CHAR* filename);


#endif
