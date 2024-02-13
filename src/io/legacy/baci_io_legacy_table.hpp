/*----------------------------------------------------------------------*/
/*! \file
\brief A very simple symbol table implementation.


\level 1

---------------------------------------------------------------------*/
/*!

ccarat writes control files that describe its binary output. Those
files are meant to be human readable and very flexible but still must
be read back. To accomplish this we need a small parser and a way to
represent these files internally. And here it is.

A control file consists of definitions. Simple definitions look like
this:

attr_name = value

where value might be a string, an integer or a double constant. (It's
easily possible to enhance the parser to understand mathematical
expressions. Speak up if you need this.) Additionally you can group
your definitions:

group_name:
    attr_name_1 = value
    attr_name_2 = value

This looks a little like python and indeed indentions do matter
here. But on the other hand you can have many groups using the same
name. That's convenient because restart information is written every
some steps into identical groups.

Files that follow this scheme can be read into a simple symbol
table. This table can be queried for those values quite easily.

*/

#ifndef BACI_IO_LEGACY_TABLE_HPP
#define BACI_IO_LEGACY_TABLE_HPP

#include "baci_config.hpp"

#include "baci_io_legacy_types.hpp"

BACI_NAMESPACE_OPEN

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
SYMBOL* map_find_symbol(MAP* map, const char* key);


/* Find the last symbols value. The value has to be of the given
 * type. Returns false on failure. */
int map_find_string(MAP* map, const char* key, const char** string);
int map_find_int(MAP* map, const char* key, int* integer);
int map_find_real(MAP* map, const char* key, double* real);
int map_find_map(MAP* map, const char* key, MAP** dir);


/* Find the last symbols value. The value has to be of the given
 * type. Calls dserror on failture. */
const char* map_read_string(MAP* map, const char* key);
int map_read_int(MAP* map, const char* key);
double map_read_real(MAP* map, const char* key);
MAP* map_read_map(MAP* map, const char* key);


/* Tell whether there is a symbol with given key and value. Only the
 * last symbol with that key is checked. */
int map_has_string(MAP* map, const char* key, const char* value);
int map_has_int(MAP* map, const char* key, const int value);
int map_has_real(MAP* map, const char* key, const double value);
int map_has_map(MAP* map, const char* key);


/* Insert a new symbol. */
void map_insert_string(MAP* map, char* string, char* key);
void map_insert_int(MAP* map, int integer, char* key);
void map_insert_real(MAP* map, double real, char* key);
void map_insert_map(MAP* map, MAP* dir, char* key);


/* Tell the number of symbols under this key. */
int map_symbol_count(MAP* map, const char* key);


/* Take a symbol chain out of the map. Leave the symbol alive. */
void map_disconnect_symbols(MAP* map, const char* key);


/* Prepend the symbol chain to one under the given key. */
void map_prepend_symbols(MAP* map, const char* key, SYMBOL* symbol, int count);


/* Tell whether this symbol has the given type. */
int symbol_is_map(const SYMBOL* symbol);


/* Extract the value of this symbol. Returns false on failture. */
int symbol_get_string(const SYMBOL* symbol, const char** string);
int symbol_get_int(const SYMBOL* symbol, int* integer);
int symbol_get_real(const SYMBOL* symbol, double* real);
int symbol_get_real_as_float(const SYMBOL* symbol, float* real);
int symbol_get_map(const SYMBOL* symbol, MAP** map);


/* Extract the value of this symbol. Call dserror on failure. */
MAP* symbol_map(const SYMBOL* symbol);


/* Read the control file given by name. Put its contents into the map. */
void parse_control_file(MAP* map, const char* filename, MPI_Comm comm);


/* Read the control file given by name. Put its contents into the map.
 * (serial only!)*/
void parse_control_file_serial(MAP* map, const char* filename);

BACI_NAMESPACE_CLOSE

#endif
