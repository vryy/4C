/*!
\file
\brief A very simple symbol table implementation.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

ccarat writes control files that describe its binary output. Those
files are meant to be human readable and very flexible but still must
be read back. To accomplish this we need a small parser and a way to
represent these files internally. And here it is.

A control file consists of definitions. Simple definitions look like
this:

<pre>
attr_name = value
</pre>

where value might be a string, an integer or a double constant. (It's
easily possible to enhance the parser to understand mathematical
expressions. Speak up if you need this.) Additionally you can group
your definitions:

<pre>
group_name:
    attr_name_1 = value
    attr_name_2 = value
</pre>

This looks a little like python and indeed indentions do matter
here. But on the other hand you can have many groups using the same
name. That's convenient because restart information is written every
some steps into identical groups.

Files that follow this scheme can be read into a simple symbol
table. This table can be queried for those values quite easily.

*/

#include "../headers/standardtypes.h"

#include "pss_table.h"


/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
extern struct _PAR   par;


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void destroy_symbol(SYMBOL* symbol)
{
#ifdef DEBUG
  dstrc_enter("destroy_symbol");
#endif

  while (symbol != NULL) {
    SYMBOL* next;

    if (symbol->type == sym_string) {
      CCAFREE(symbol->s.string);
    }
    if (symbol->type == sym_map) {
      destroy_map(symbol->s.dir);
      CCAFREE(symbol->s.dir);
    }
    next = symbol->next;
    CCAFREE(symbol);
    symbol = next;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void destroy_node(MAP_NODE* node)
{
#ifdef DEBUG
  dstrc_enter("destroy_node");
#endif

  if (node != NULL) {
    destroy_node(node->lhs);
    destroy_node(node->rhs);
    destroy_symbol(node->symbol);
    CCAFREE(node->key);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Bring a map variable up to a clean state.

  That's needed before anything can be done with a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void init_map(MAP* map)
{
#ifdef DEBUG
  dstrc_enter("init_map");
#endif

  /* We have a dummy node at the root to make life easier. The empty
   * key is not legal. */
  map->root.key = "";
  map->root.symbol = NULL; /*CCACALLOC(1, sizeof(SYMBOL));*/
  map->root.lhs = NULL;
  map->root.rhs = NULL;
  map->count = 0;

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void destroy_map(MAP* map)
{
#ifdef DEBUG
  dstrc_enter("destroy_map");
#endif

  destroy_node(map->root.lhs);
  destroy_node(map->root.rhs);
  /*CCAFREE(map->root.symbol);*/

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief See whether a node matches a certain key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static INT map_cmp_nodes(MAP_NODE* lhs, CHAR* rhs_key)
{
  return strcmp(lhs->key, rhs_key);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the node in the map that matches the \a key.

  \return NULL if there's no such node.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
static MAP_NODE* map_find_node(MAP* map, CHAR* key)
{
  MAP_NODE* node;

#ifdef DEBUG
  dstrc_enter("map_find_node");
#endif

  node = &(map->root);

  /* we keep hitting this problem */
  dsassert(node->key != NULL, "input system not initialized");

  for (;;) {
    INT cmp;
    cmp = map_cmp_nodes(node, key);
    if (cmp < 0) {
      if (node->rhs == NULL) {
        node = NULL;
        goto end;
      }
      else {
        node = node->rhs;
      }
    }
    else if (cmp > 0) {
      if (node->lhs == NULL) {
        node = NULL;
        goto end;
      }
      else {
        node = node->lhs;
      }
    }
    else {
      goto end;
    }
  }

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return node;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol in the map that matches the \a key.

  \return NULL if there's no such symbol.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
SYMBOL* map_find_symbol(MAP* map, CHAR* key)
{
  MAP_NODE* node;
  SYMBOL* symbol = NULL;

#ifdef DEBUG
  dstrc_enter("map_find_symbol");
#endif

  node = map_find_node(map, key);
  if (node != NULL) {
    symbol = node->symbol;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return symbol;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a string.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_string(MAP* map, CHAR* key, CHAR** string)
{
  SYMBOL* symbol;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_find_string");
#endif

  symbol = map_find_symbol(map, key);
  ret = symbol_get_string(symbol, string);

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a integer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_int(MAP* map, CHAR* key, INT* integer)
{
  SYMBOL* symbol;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_find_int");
#endif

  symbol = map_find_symbol(map, key);
  ret = symbol_get_int(symbol, integer);

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a real.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_real(MAP* map, CHAR* key, DOUBLE* real)
{
  SYMBOL* symbol;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_find_real");
#endif

  symbol = map_find_symbol(map, key);
  ret = symbol_get_real(symbol, real);

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_map(MAP* map, CHAR* key, MAP** dir)
{
  SYMBOL* symbol;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_find_map");
#endif

  symbol = map_find_symbol(map, key);
  ret = symbol_get_map(symbol, dir);

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a string.

  Stops if no string is found.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
CHAR* map_read_string(MAP* map, CHAR* key)
{
  CHAR* string;
#ifdef DEBUG
  dstrc_enter("map_read_string");
#endif

  if (!map_find_string(map, key, &string)) {
    dserror("no string attribute '%s' in map", key);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return string;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a integer.

  Stops if no integer is found.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_read_int(MAP* map, CHAR* key)
{
  INT integer;
#ifdef DEBUG
  dstrc_enter("map_read_int");
#endif

  if (!map_find_int(map, key, &integer)) {
    dserror("no int attribute '%s' in map", key);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return integer;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a real.

  Stops if no real is found.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
DOUBLE map_read_real(MAP* map, CHAR* key)
{
  DOUBLE real;
#ifdef DEBUG
  dstrc_enter("map_read_real");
#endif

  if (!map_find_real(map, key, &real)) {
    dserror("no real attribute '%s' in map", key);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return real;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a map.

  Stops if no map is found.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
MAP* map_read_map(MAP* map, CHAR* key)
{
  MAP* dir;
#ifdef DEBUG
  dstrc_enter("map_read_map");
#endif

  if (!map_find_map(map, key, &dir)) {
    dserror("no dir attribute '%s' in map", key);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return dir;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key has this \a value.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_string(MAP* map, CHAR* key, CHAR* value)
{
  SYMBOL* symbol;
  CHAR* string;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_has_string");
#endif

  symbol = map_find_symbol(map, key);
  if (symbol != NULL) {
    ret = symbol_get_string(symbol, &string);
    if (ret) {
      ret = strcmp(string, value) == 0;
    }
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key has this \a value.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_int(MAP* map, CHAR* key, INT value)
{
  SYMBOL* symbol;
  INT integer;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_has_int");
#endif

  symbol = map_find_symbol(map, key);
  if (symbol != NULL) {
    ret = symbol_get_int(symbol, &integer);
    if (ret) {
      ret = integer == value;
    }
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key has this \a value.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_real(MAP* map, CHAR* key, DOUBLE value)
{
  SYMBOL* symbol;
  DOUBLE real;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_has_real");
#endif

  symbol = map_find_symbol(map, key);
  if (symbol != NULL) {
    ret = symbol_get_real(symbol, &real);
    if (ret) {
      ret = real == value;
    }
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key is a map.

  No value comparison here.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_map(MAP* map, CHAR* key)
{
  SYMBOL* symbol;
  INT ret;

#ifdef DEBUG
  dstrc_enter("map_has_real");
#endif

  symbol = map_find_symbol(map, key);
  if (symbol != NULL) {
    ret = symbol_is_map(symbol);
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert a symbol.

  Any new symbol becomes the first one with that key.

  Ownership of the symbol and the key is taken. Both have to be
  allocated using CCAMALLOC or the like.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void map_insert_symbol(MAP* map, SYMBOL* symbol, CHAR* key)
{
  MAP_NODE* node;

#ifdef DEBUG
  dstrc_enter("map_insert_symbol");
#endif

  node = &(map->root);
  for (;;) {
    INT cmp;
    cmp = map_cmp_nodes(node, key);
    if (cmp < 0) {
      if (node->rhs == NULL) {
        node->rhs = CCACALLOC(1, sizeof(MAP_NODE));
        node->rhs->key = key;
        node->rhs->symbol = symbol;
        node->rhs->count = 1;
        map->count++;
        goto end;
      }
      else {
        node = node->rhs;
      }
    }
    else if (cmp > 0) {
      if (node->lhs == NULL) {
        node->lhs = CCACALLOC(1, sizeof(MAP_NODE));
        node->lhs->key = key;
        node->lhs->symbol = symbol;
        node->lhs->count = 1;
        map->count++;
        goto end;
      }
      else {
        node = node->lhs;
      }
    }
    else {
      /* This key is already there. Free the duplicated memory. */
      CCAFREE(key);

      /* append symbol */
      symbol->next = node->symbol;
      node->symbol = symbol;
      node->count++;
      map->count++;
      goto end;
    }
  }

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that string with this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_string(MAP* map, CHAR* string, CHAR* key)
{
  SYMBOL* symbol;

#ifdef DEBUG
  dstrc_enter("map_insert_string");
#endif

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_string;
  symbol->s.string = string;

  map_insert_symbol(map, symbol, key);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that integer with this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_int(MAP* map, INT integer, CHAR* key)
{
  SYMBOL* symbol;

#ifdef DEBUG
  dstrc_enter("map_insert_int");
#endif

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_int;
  symbol->s.integer = integer;

  map_insert_symbol(map, symbol, key);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that real with this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_real(MAP* map, DOUBLE real, CHAR* key)
{
  SYMBOL* symbol;

#ifdef DEBUG
  dstrc_enter("map_insert_real");
#endif

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_real;
  symbol->s.real = real;

  map_insert_symbol(map, symbol, key);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that map with this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_map(MAP* map, MAP* dir, CHAR* key)
{
  SYMBOL* symbol;

#ifdef DEBUG
  dstrc_enter("map_insert_map");
#endif

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_map;
  symbol->s.dir = dir;

  map_insert_symbol(map, symbol, key);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert a copy of that string with a copy of this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_string_cpy(MAP* map, CHAR* string, CHAR* key)
{
  CHAR* string_cpy;
  CHAR* key_cpy;

#ifdef DEBUG
  dstrc_enter("map_insert_string_cpy");
#endif

  string_cpy = (CHAR*)CCAMALLOC((strlen(string)+1)*sizeof(CHAR));
  strcpy(string_cpy, string);
  key_cpy = (CHAR*)CCAMALLOC((strlen(key)+1)*sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_string(map, string_cpy, key_cpy);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that integer with a copy of this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_int_cpy(MAP* map, INT integer, CHAR* key)
{
  CHAR* key_cpy;

#ifdef DEBUG
  dstrc_enter("map_insert_int_cpy");
#endif

  key_cpy = (CHAR*)CCAMALLOC((strlen(key)+1)*sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_int(map, integer, key_cpy);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that real with a copy of this key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_real_cpy(MAP* map, DOUBLE real, CHAR* key)
{
  CHAR* key_cpy;

#ifdef DEBUG
  dstrc_enter("map_insert_real_cpy");
#endif

  key_cpy = (CHAR*)CCAMALLOC((strlen(key)+1)*sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_real(map, real, key_cpy);

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Insert that map with a copy of this key.

  \warning The map to be inserted itself is not copied. It must have
  been allocated using CCACALLOC (or the like) and will be owned by
  the map from now on.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_insert_map_cpy(MAP* map, MAP* dir, CHAR* key)
{
  CHAR* key_cpy;

#ifdef DEBUG
  dstrc_enter("map_insert_map_cpy");
#endif

  key_cpy = (CHAR*)CCAMALLOC((strlen(key)+1)*sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_map(map, dir, key_cpy);

#ifdef DEBUG
  dstrc_exit();
#endif
}


INT map_symbol_count(MAP* map, CHAR* key)
{
  INT count = 0;
  MAP_NODE* node;

#ifdef DEBUG
  dstrc_enter("map_symbol_count");
#endif

  node = map_find_node(map, key);
  if (node != NULL) {
    count = node->count;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return count;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is a string.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_string(SYMBOL* symbol)
{
  return (symbol != NULL) && (symbol->type == sym_string);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is an integer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_int(SYMBOL* symbol)
{
  return (symbol != NULL) && (symbol->type == sym_int);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is a real.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_real(SYMBOL* symbol)
{
  return (symbol != NULL) && (symbol->type == sym_real);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_map(SYMBOL* symbol)
{
  return (symbol != NULL) && (symbol->type == sym_map);
}



/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a string.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_string(SYMBOL* symbol, CHAR** string)
{
  INT ret;

#ifdef DEBUG
  dstrc_enter("symbol_get_string");
#endif

  if (symbol && (symbol->type == sym_string)) {
    *string = symbol->s.string;
    ret = 1;
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its an integer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_int(SYMBOL* symbol, INT* integer)
{
  INT ret;

#ifdef DEBUG
  dstrc_enter("symbol_get_int");
#endif

  if (symbol && (symbol->type == sym_int)) {
    *integer = symbol->s.integer;
    ret = 1;
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a real.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_real(SYMBOL* symbol, DOUBLE* real)
{
  INT ret;

#ifdef DEBUG
  dstrc_enter("symbol_get_real");
#endif

  if (symbol && (symbol->type == sym_real)) {
    *real = symbol->s.real;
    ret = 1;
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_map(SYMBOL* symbol, MAP** map)
{
  INT ret;

#ifdef DEBUG
  dstrc_enter("symbol_get_map");
#endif

  if (symbol && (symbol->type == sym_map)) {
    *map = symbol->s.dir;
    ret = 1;
  }
  else {
    ret = 0;
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a string.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
CHAR* symbol_string(SYMBOL* symbol)
{
  CHAR* ret;
#ifdef DEBUG
  dstrc_enter("symbol_string");
#endif

  if (symbol->type == sym_string) {
    ret = symbol->s.string;
  }
  else {
    dserror("Wrong symbol type %d", symbol->type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its an integer.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT symbol_int(SYMBOL* symbol)
{
  INT ret;
#ifdef DEBUG
  dstrc_enter("symbol_int");
#endif

  if (symbol->type == sym_int) {
    ret = symbol->s.integer;
  }
  else {
    dserror("Wrong symbol type %d", symbol->type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a double.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
DOUBLE symbol_real(SYMBOL* symbol)
{
  DOUBLE ret;
#ifdef DEBUG
  dstrc_enter("symbol_real");
#endif

  if (symbol->type == sym_real) {
    ret = symbol->s.real;
  }
  else {
    dserror("Wrong symbol type %d", symbol->type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a map.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
MAP* symbol_map(SYMBOL* symbol)
{
  MAP* ret;
#ifdef DEBUG
  dstrc_enter("symbol_map");
#endif

  if (symbol->type == sym_map) {
    ret = symbol->s.dir;
  }
  else {
    dserror("Wrong symbol type %d", symbol->type);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}


#define PRINT_INDENT(f, count) { INT i; for (i=0; i<(count); ++i) fprintf(f, " "); }


/*----------------------------------------------------------------------*/
/*!
  \brief Write all symbols with the given key.

  Here recursion is used in order to write the last symbol first. This
  way we'll get an identical map when we read it back instead of a map
  with it's symbols in reverse order. Of course there are limits when
  using recursion but we don't expect to hit them here. There maps are
  supposed to be relatively small.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void symbol_print(FILE* f, CHAR* key, SYMBOL* symbol, INT indent)
{
#ifdef DEBUG
  dstrc_enter("symbol_print");
#endif

  dsassert(strstr(key, "\n")==NULL, "malformed key");

  if (symbol != NULL) {
    symbol_print(f, key, symbol->next, indent);

    PRINT_INDENT(f, indent);
    switch (symbol->type) {
    case sym_string:
      fprintf(f, "%s = \"%s\"\n", key, symbol->s.string);
      break;
    case sym_int:
      fprintf(f, "%s = %d\n", key, symbol->s.integer);
      break;
    case sym_real:
      /* we don't care for beauty but we need many digits. */
      fprintf(f, "%s = %.20e\n", key, symbol->s.real);
      break;
    case sym_map:
      fprintf(f, "%s:\n", key);
      map_print(f, symbol->s.dir, indent+4);
      break;
    default:
      dserror("Ups!");
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write this node and its subnodes.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void map_node_print(FILE* f, MAP_NODE* node, INT indent)
{
#ifdef DEBUG
  dstrc_enter("map_node_print");
#endif

  if (node != NULL) {
    map_node_print(f, node->lhs, indent);
    symbol_print(f, node->key, node->symbol, indent);
    map_node_print(f, node->rhs, indent);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write the maps content to file \a f using \a indent as basic
  indention level.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_print(FILE* f, MAP* map, INT indent)
{
#ifdef DEBUG
  dstrc_enter("map_print");
#endif

  map_node_print(f, map->root.rhs, indent);
  fprintf(f, "\n");

#ifdef DEBUG
  dstrc_exit();
#endif
}



/*----------------------------------------------------------------------*/
/*!
  \brief The types of tokens recognized by the lexer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef enum _TOKEN_TYPE {
  tok_none,
  tok_done,
  tok_name,
  tok_string,
  tok_int,
  tok_real,
  tok_colon,
  tok_equal,
  tok_indent,
  tok_dedent
} TOKEN_TYPE;


/*----------------------------------------------------------------------*/
/*!
 \brief The parsers internal state

 Here we have the parsers internal state. It consists of the current
 token (depending on the token type these variables have different
 meanings), the file buffer with the current read position, the
 current line number and indention level. These variables are very
 internal and only used while a control file is read.

  \author u.kue
  \date 08/04
 */
/*----------------------------------------------------------------------*/
typedef struct _PARSER_DATA {
  TOKEN_TYPE tok;
  CHAR* token_string;
  INT token_int;
  DOUBLE token_real;

  CHAR* file_buffer;
  INT file_size;
  CHAR* filename;

  INT pos;
  INT lineno;
  INT indent_level;
  int indent_step;
} PARSER_DATA;


/*----------------------------------------------------------------------*/
/*!
  \brief Init the data structure needed to read a file.

  The file is read on processor 0 and broadcasted to the others.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void init_parser_data(struct _PARSER_DATA* data, CHAR* filename)
{
#ifdef DEBUG
  dstrc_enter("init_parser_data");
#endif

  data->tok = tok_none;
  data->lineno = 1;
  data->pos = 0;
  data->indent_level = 0;
  data->indent_step = -1;

  /* No copy here. Valid only as long as the calling functions
   * filename is valid. */
  data->filename = filename;

  /* We need to have the information on all processes. That's why we
   * read the file on process 0 and broadcast it. The other way would
   * be to use MPI IO, but then we'd have to implement a separat
   * sequential version. */
#ifdef PARALLEL
  if (par.myrank == 0)
#endif
  {
    INT bytes_read;
    FILE* file;
    file = fopen(filename, "rb");

    if (file==NULL) {
      dserror("cannot read file '%s'", filename);
    }

    /* find out the control file size */
    fseek(file, 0, SEEK_END);
    data->file_size = ftell(file);

    /* read file to local buffer */
    data->file_buffer = CCAMALLOC((data->file_size+1)*sizeof(CHAR));
    fseek(file, 0, SEEK_SET);
    bytes_read = fread(data->file_buffer, sizeof(CHAR), data->file_size, file);
    if (bytes_read != data->file_size) {
      dserror("failed to read file %s", filename);
    }
    /* a trailing zero helps a lot */
    data->file_buffer[data->file_size] = '\0';

    fclose(file);
  }

#ifdef PARALLEL
  if (par.nprocs > 1) {
    INT err;
    err = MPI_Bcast(&data->file_size,1,MPI_INT,0,MPI_COMM_WORLD);
    if (err != 0) {
      dserror("MPI_Bcast failed: %d", err);
    }
    if (par.myrank > 0) {
      data->file_buffer = CCAMALLOC((data->file_size+1)*sizeof(CHAR));
    }
    err = MPI_Bcast(data->file_buffer, data->file_size+1, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (err != 0) {
      dserror("MPI_Bcast failed: %d", err);
    }
  }
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void destroy_parser_data(struct _PARSER_DATA* data)
{
#ifdef DEBUG
  dstrc_enter("destroy_parser_data");
#endif

  CCAFREE(data->file_buffer);

#ifdef DEBUG
  dstrc_exit();
#endif
}


#if 0

/*----------------------------------------------------------------------*/
/*!
  \brief Debug output.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void print_token(PARSER_DATA* data)
{
  CHAR buf[256];

  switch (data->tok) {
  case tok_none:
    printf("<none>");
    break;
  case tok_done:
    printf("<eof>");
    break;
  case tok_name:
    strncpy(buf, data->token_string, data->token_int);
    buf[data->token_int] = '\0';
    printf("<name: %s, %d>", buf, data->token_int);
    break;
  case tok_string:
    strncpy(buf, data->token_string, data->token_int);
    buf[data->token_int] = '\0';
    printf("<string: %s, %d>", buf, data->token_int);
    break;
  case tok_int:
    printf("<int: %d>", data->token_int);
    break;
  case tok_real:
    printf("<real: %f>", data->token_real);
    break;
  case tok_colon:
    printf("<colon>");
    break;
  case tok_equal:
    printf("<equal>");
    break;
  case tok_indent:
    printf("<indent: %d, %d>", data->token_int, data->indent_level);
    break;
  case tok_dedent:
    printf("<dedent: %d, %d>", data->token_int, data->indent_level);
    break;
  default:
    printf("<unknown>");
  }
}

#endif


/*----------------------------------------------------------------------*/
/*!
  \brief Get the next char.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static int getnext(PARSER_DATA* data)
{
  if (data->pos < data->file_size) {

    /* ignore dos line endings */
    if (data->file_buffer[data->pos] == '\r') {
      data->pos++;
    }

    /* Increment the counter and return the char at the old position. */
    return data->file_buffer[data->pos++];
  }
  return EOF;
}


#define TABWIDTH 8


/*----------------------------------------------------------------------*/
/*!
  \brief Get the next token.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void lexan(PARSER_DATA* data)
{
  INT line_begin = 0;
  INT t;
  INT indention = data->indent_level;

#ifdef DEBUG
  dstrc_enter("lexan");
#endif

  for (;;) {
    t = getnext(data);
    if (t == ' ') {
      /* ignore whitespaces */
      if (line_begin) {
        indention++;
      }
    }
    else if (t == '\t') {
      /* ignore whitespaces */
      if (line_begin) {
        indention = ((indention + TABWIDTH - 1) / TABWIDTH) * TABWIDTH;
      }
    }
    else if (t == '\n') {
      data->lineno++;
      line_begin = 1;
      indention = 0;
    }
    else if (t == '#') {
      for (;;) {
        t = getnext(data);
        if (t == '\n') {
          break;
        }
      }
      data->lineno++;
      line_begin = 1;
      indention = 0;
    }
    else if (t == EOF) {
      data->tok = tok_done;
      goto end;
    }
    else {

      if (line_begin && (indention != data->indent_level)) {
        if (data->indent_step == -1) {
          if (indention > data->indent_level) {
            dsassert(data->indent_level == 0, "non-zero intention at first line?!");
            data->indent_step = indention;
            data->indent_level = indention;
            data->token_int = 1;
            data->tok = tok_indent;
          }
          else {
            dserror("dedent at toplevel!");
          }
        }
        else {
          if (indention > data->indent_level) {
            data->tok = tok_indent;
            dsassert((indention - data->indent_level) % data->indent_step == 0, "malformed indention");
            data->token_int = (indention - data->indent_level) / data->indent_step;
            data->indent_level = indention;
          }
          else {
            data->tok = tok_dedent;
            dsassert((data->indent_level - indention) % data->indent_step == 0, "malformed dedention");
            data->token_int = (data->indent_level - indention) / data->indent_step;
            data->indent_level = indention;
          }
        }
        data->pos--;
        goto end;
      }
      else {

        line_begin = 0;

        if ((t == '-') || isdigit(t)) {
          data->token_string = &(data->file_buffer[data->pos-1]);
          if (t == '-') {
            t = getnext(data);
          }
          while (isdigit(t)) {
            t = getnext(data);
          }
          if ((t != '.') && (t != 'E') && (t != 'e')) {
            if (t != EOF) {
              data->pos--;
            }
            data->token_int = atoi(data->token_string);
            data->tok = tok_int;
            goto end;
          }
          if (t == '.') {
            t = getnext(data);
            if (isdigit(t)) {
              while (isdigit(t)) {
                t = getnext(data);
              }
            }
            else {
              dserror("no digits after point at line %d", data->lineno);
            }
          }
          if ((t == 'E') || (t == 'e')) {
            t = getnext(data);
            if ((t == '-') || (t == '+')) {
              t = getnext(data);
            }
            if (isdigit(t)) {
              while (isdigit(t)) {
                t = getnext(data);
              }
            }
            else {
              dserror("no digits after exponent at line %d", data->lineno);
            }
          }
          if (t != EOF) {
            data->pos--;
          }
          data->token_real = strtod(data->token_string, NULL);
          data->tok = tok_real;
          goto end;
        }
        else if (isalpha(t) || (t == '_')) {
          data->token_string = &(data->file_buffer[data->pos-1]);
          while (isalnum(t) || (t == '_')) {
            t = getnext(data);
          }
          if (t != EOF) {
            data->pos--;
          }
          data->tok = tok_name;
          data->token_int = &(data->file_buffer[data->pos]) - data->token_string;
          goto end;
        }
        else if (t == '"') {
          data->token_string = &(data->file_buffer[data->pos]);
          t = getnext(data);
          while (t != '"') {
            t = getnext(data);
            if (t==EOF) {
              dserror("expected closing \" on line %d", data->lineno);
            }
          }
          data->tok = tok_string;
          data->token_int = &(data->file_buffer[data->pos-1]) - data->token_string;
          goto end;
        }
        else if (t == ':') {
          data->tok = tok_colon;
          goto end;
        }
        else if (t == '=') {
          data->tok = tok_equal;
          goto end;
        }
        else {
          dserror("unexpected char '%c' at line %d", t, data->lineno);
          data->tok = tok_none;
          goto end;
        }
      }
    }
  }

end:

#if 0
  if (par.myrank == 0) {
    print_token(data);
    printf("\n");
  }
#endif

#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief The top down parser.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void parse_definitions(PARSER_DATA* data, MAP* dir)
{
#ifdef DEBUG
  dstrc_enter("parse_definitions");
#endif

  lexan(data);

  while (data->tok != tok_done) {
    switch (data->tok) {
    case tok_name: {
      CHAR* name;

      /*
       * The string is not null terminated as it's a simple pointer
       * into the file buffer. However, we know its length so we can
       * handle that. */
      name = CCAMALLOC((data->token_int+1)*sizeof(CHAR));
      strncpy(name, data->token_string, data->token_int);
      name[data->token_int] = '\0';

      lexan(data);
      switch (data->tok) {
      case tok_colon: {
        MAP* map;

        lexan(data);
        if ((data->tok != tok_indent) || (data->token_int != 1)) {
          dserror("Syntaxerror at line %d: single indention expected", data->lineno);
        }

        map = CCACALLOC(1, sizeof(MAP));
        init_map(map);
        parse_definitions(data, map);

        map_insert_map(dir, map, name);

        if ((data->tok == tok_dedent) && (data->token_int > 0)) {
          data->token_int--;
          goto end;
        }

        break;
      }
      case tok_equal:
        lexan(data);
        switch (data->tok) {
        case tok_string: {
          CHAR* string;

          /* Again, be carefully with those pointers... */
          string = CCAMALLOC((data->token_int+1)*sizeof(CHAR));
          strncpy(string, data->token_string, data->token_int);
          string[data->token_int] = '\0';

          map_insert_string(dir, string, name);
          break;
        }
        case tok_int:
          map_insert_int(dir, data->token_int, name);
          break;
        case tok_real:
          map_insert_real(dir, data->token_real, name);
          break;
        default:
          dserror("Syntaxerror at line %d: string, int or real expected", data->lineno);
        }
        break;
      default:
        dserror("Syntaxerror at line %d: ':' or '=' expected", data->lineno);
      }
      break;
    }
    case tok_dedent:
      data->token_int--;
      goto end;
    default:
      dserror("Syntaxerror at line %d: name expected", data->lineno);
    }

    lexan(data);
  }

end:
#ifdef DEBUG
  dstrc_exit();
#endif
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Parse the file given by name and fill the map with this
  file's content.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void parse_control_file(MAP* map, CHAR* filename)
{
  PARSER_DATA data;

#ifdef DEBUG
  dstrc_enter("parse_control_file");
#endif

  /*
   * So here we are. Before the symbol table can be filled with values
   * it has to be initialized. That is we expect to get an
   * uninitialized (virgin) map. */
  init_map(map);

  init_parser_data(&data, filename);
  parse_definitions(&data, map);
  destroy_parser_data(&data);

#ifdef DEBUG
  dstrc_exit();
#endif
}
