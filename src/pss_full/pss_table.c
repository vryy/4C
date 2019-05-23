/*!---------------------------------------------------------------------
\brief A very simple symbol table implementation.

\maintainer Martin Kronbichler

\level 1

---------------------------------------------------------------------*/
/*!
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
#include "../drt_lib/drt_dserror.H"

#include "pss_table.h"
#include "pss_prototypes.h"


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void destroy_symbol(SYMBOL* symbol)
{
  while (symbol != NULL)
  {
    SYMBOL* next;

    if (symbol->type == sym_string)
    {
      CCAFREE(symbol->s.string);
    }
    if (symbol->type == sym_map)
    {
      destroy_map(symbol->s.dir);
      CCAFREE(symbol->s.dir);
    }
    next = symbol->next;
    CCAFREE(symbol);
    symbol = next;
  }
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
  if (node != NULL)
  {
    destroy_node(node->lhs);
    destroy_node(node->rhs);
    if (node->symbol != NULL)
    {
      destroy_symbol(node->symbol);
    }
    CCAFREE(node->key);

    node->lhs = 0;
    node->rhs = 0;
    node->symbol = 0;
    node->key = 0;

    CCAFREE(node);
  }
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
  /* We have a dummy node at the root to make life easier. The empty
   * key is not legal. */
  map->root.key = "";
  map->root.symbol = NULL; /*CCACALLOC(1, sizeof(SYMBOL));*/
  map->root.lhs = NULL;
  map->root.rhs = NULL;
  map->count = 0;
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
  destroy_node(map->root.lhs);
  destroy_node(map->root.rhs);
  /*CCAFREE(map->root.symbol);*/
}


/*----------------------------------------------------------------------*/
/*!
  \brief See whether a node matches a certain key.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static INT map_cmp_nodes(const MAP_NODE* lhs, const CHAR* rhs_key)
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
static MAP_NODE* map_find_node(MAP* map, const CHAR* key)
{
  MAP_NODE* node;

  node = &(map->root);

  /* we keep hitting this problem */
  dsassert(node->key != NULL, "input system not initialized");

  for (;;)
  {
    INT cmp;
    cmp = map_cmp_nodes(node, key);
    if (cmp < 0)
    {
      if (node->rhs == NULL)
      {
        node = NULL;
        goto end;
      }
      else
      {
        node = node->rhs;
      }
    }
    else if (cmp > 0)
    {
      if (node->lhs == NULL)
      {
        node = NULL;
        goto end;
      }
      else
      {
        node = node->lhs;
      }
    }
    else
    {
      goto end;
    }
  }

end:
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
SYMBOL* map_find_symbol(MAP* map, const CHAR* key)
{
  MAP_NODE* node;
  SYMBOL* symbol = NULL;

  node = map_find_node(map, key);
  if (node != NULL)
  {
    symbol = node->symbol;
  }

  return symbol;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a string.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_string(MAP* map, const CHAR* key, CHAR** string)
{
  SYMBOL* symbol;
  INT ret;

  symbol = map_find_symbol(map, key);
  ret = symbol_get_string(symbol, string);

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a integer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_int(MAP* map, const CHAR* key, INT* integer)
{
  SYMBOL* symbol;
  INT ret;

  symbol = map_find_symbol(map, key);
  ret = symbol_get_int(symbol, integer);

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a real.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_real(MAP* map, const CHAR* key, DOUBLE* real)
{
  SYMBOL* symbol;
  INT ret;

  symbol = map_find_symbol(map, key);
  ret = symbol_get_real(symbol, real);

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Find the first symbol, return its value if it's a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_find_map(MAP* map, const CHAR* key, MAP** dir)
{
  SYMBOL* symbol;
  INT ret;

  symbol = map_find_symbol(map, key);
  ret = symbol_get_map(symbol, dir);

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
CHAR* map_read_string(MAP* map, const CHAR* key)
{
  CHAR* string;

  if (!map_find_string(map, key, &string))
  {
    dserror("no string attribute '%s' in map", key);
  }

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
INT map_read_int(MAP* map, const CHAR* key)
{
  INT integer;

  if (!map_find_int(map, key, &integer))
  {
    dserror("no int attribute '%s' in map", key);
  }

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
DOUBLE map_read_real(MAP* map, const CHAR* key)
{
  DOUBLE real;

  if (!map_find_real(map, key, &real))
  {
    INT value;
    if (!map_find_int(map, key, &value)) dserror("no real attribute '%s' in map", key);
    real = value;
  }

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
MAP* map_read_map(MAP* map, const CHAR* key)
{
  MAP* dir;

  if (!map_find_map(map, key, &dir))
  {
    dserror("no dir attribute '%s' in map", key);
  }

  return dir;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key has this \a value.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_string(MAP* map, const CHAR* key, const CHAR* value)
{
  SYMBOL* symbol;
  CHAR* string;
  INT ret;

  symbol = map_find_symbol(map, key);
  if (symbol != NULL)
  {
    ret = symbol_get_string(symbol, &string);
    if (ret)
    {
      ret = strcmp(string, value) == 0;
    }
  }
  else
  {
    ret = 0;
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key has this \a value.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_int(MAP* map, const CHAR* key, const INT value)
{
  SYMBOL* symbol;
  INT integer;
  INT ret;

  symbol = map_find_symbol(map, key);
  if (symbol != NULL)
  {
    ret = symbol_get_int(symbol, &integer);
    if (ret)
    {
      ret = integer == value;
    }
  }
  else
  {
    ret = 0;
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell if the first symbol with that \a key has this \a value.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT map_has_real(MAP* map, const CHAR* key, const DOUBLE value)
{
  SYMBOL* symbol;
  DOUBLE real;
  INT ret;

  symbol = map_find_symbol(map, key);
  if (symbol != NULL)
  {
    ret = symbol_get_real(symbol, &real);
    if (ret)
    {
      ret = real == value;
    }
  }
  else
  {
    ret = 0;
  }

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
INT map_has_map(MAP* map, const CHAR* key)
{
  SYMBOL* symbol;
  INT ret;

  symbol = map_find_symbol(map, key);
  if (symbol != NULL)
  {
    ret = symbol_is_map(symbol);
  }
  else
  {
    ret = 0;
  }

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

  node = &(map->root);
  for (;;)
  {
    INT cmp;
    cmp = map_cmp_nodes(node, key);
    if (cmp < 0)
    {
      if (node->rhs == NULL)
      {
        node->rhs = CCACALLOC(1, sizeof(MAP_NODE));
        node->rhs->key = key;
        node->rhs->symbol = symbol;
        node->rhs->count = 1;
        node->rhs->lhs = NULL;
        node->rhs->rhs = NULL;
        map->count++;
        goto end;
      }
      else
      {
        node = node->rhs;
      }
    }
    else if (cmp > 0)
    {
      if (node->lhs == NULL)
      {
        node->lhs = CCACALLOC(1, sizeof(MAP_NODE));
        node->lhs->key = key;
        node->lhs->symbol = symbol;
        node->lhs->count = 1;
        node->lhs->lhs = NULL;
        node->lhs->rhs = NULL;
        map->count++;
        goto end;
      }
      else
      {
        node = node->lhs;
      }
    }
    else
    {
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

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_string;
  symbol->s.string = string;

  map_insert_symbol(map, symbol, key);
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

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_int;
  symbol->s.integer = integer;

  map_insert_symbol(map, symbol, key);
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

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_real;
  symbol->s.real = real;

  map_insert_symbol(map, symbol, key);
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

  symbol = CCACALLOC(1, sizeof(SYMBOL));
  symbol->type = sym_map;
  symbol->s.dir = dir;

  map_insert_symbol(map, symbol, key);
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

  string_cpy = (CHAR*)CCAMALLOC((strlen(string) + 1) * sizeof(CHAR));
  strcpy(string_cpy, string);
  key_cpy = (CHAR*)CCAMALLOC((strlen(key) + 1) * sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_string(map, string_cpy, key_cpy);
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

  key_cpy = (CHAR*)CCAMALLOC((strlen(key) + 1) * sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_int(map, integer, key_cpy);
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

  key_cpy = (CHAR*)CCAMALLOC((strlen(key) + 1) * sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_real(map, real, key_cpy);
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

  key_cpy = (CHAR*)CCAMALLOC((strlen(key) + 1) * sizeof(CHAR));
  strcpy(key_cpy, key);
  map_insert_map(map, dir, key_cpy);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell how many symbols of the given name there are.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT map_symbol_count(MAP* map, const CHAR* key)
{
  INT count = 0;

  const MAP_NODE* node = map_find_node(map, key);
  if (node != NULL)
  {
    count = node->count;
  }

  return count;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Take a symbol chain out of the map. Leave the symbol alive.

  This is for the experienced user only. A symbol chain is removed
  from the map, but the key (the node behind it) stays alive. Also the
  symbols are not deallocated. The caller must already have a pointer
  to the symbol chain and takes responsibility for it.

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void map_disconnect_symbols(MAP* map, const CHAR* key)
{
  MAP_NODE* node = map_find_node(map, key);
  if (node != NULL)
  {
    node->symbol = NULL;
    node->count = 0;
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Prepend the symbol chain to one under the given key.

  \param map    (i/o) map we work with
  \param key      (i) key to those chain we want to prepend
  \param symbol   (i) start of the new symbol chain
  \param count    (i) number of symbol in the chain

  \author u.kue
  \date 12/04
*/
/*----------------------------------------------------------------------*/
void map_prepend_symbols(MAP* map, const CHAR* key, SYMBOL* symbol, INT count)
{
  MAP_NODE* node;

  node = map_find_node(map, key);
  if (node != NULL)
  {
    if (node->symbol != NULL)
    {
      SYMBOL* s;
      s = node->symbol;

      while (s->next != NULL)
      {
        s = s->next;
      }
      s->next = symbol;
      node->count += count;
    }
    else
    {
      node->symbol = symbol;
      node->count = count;
    }

    map->count += count;
  }
  else
  {
    dserror("no node for key '%s'", key);
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is a string.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_string(const SYMBOL* symbol)
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
INT symbol_is_int(const SYMBOL* symbol) { return (symbol != NULL) && (symbol->type == sym_int); }


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is a real.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_real(const SYMBOL* symbol) { return (symbol != NULL) && (symbol->type == sym_real); }


/*----------------------------------------------------------------------*/
/*!
  \brief Tell whether the symbol is a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_is_map(const SYMBOL* symbol) { return (symbol != NULL) && (symbol->type == sym_map); }



/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a string.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_string(const SYMBOL* symbol, CHAR** string)
{
  INT ret;

  if (symbol && (symbol->type == sym_string))
  {
    *string = symbol->s.string;
    ret = 1;
  }
  else
  {
    *string = "";
    ret = 0;
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its an integer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_int(const SYMBOL* symbol, INT* integer)
{
  INT ret;

  if (symbol && (symbol->type == sym_int))
  {
    *integer = symbol->s.integer;
    ret = 1;
  }
  else
  {
    *integer = 0;
    ret = 0;
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a real.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_real(const SYMBOL* symbol, DOUBLE* real)
{
  INT ret;

  if (symbol && (symbol->type == sym_real))
  {
    *real = symbol->s.real;
    ret = 1;
  }
  else
  {
    *real = 0.0;
    ret = 0;
  }

  return ret;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a real.

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
INT symbol_get_real_as_float(const SYMBOL* symbol, float* real)
{
  INT ret;

  if (symbol && (symbol->type == sym_real))
  {
    *real = (float)symbol->s.real;
    ret = 1;
  }
  else
  {
    *real = 0.0;
    ret = 0;
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a map.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
INT symbol_get_map(const SYMBOL* symbol, MAP** map)
{
  INT ret;

  if (symbol && (symbol->type == sym_map))
  {
    *map = symbol->s.dir;
    ret = 1;
  }
  else
  {
    *map = NULL;
    ret = 0;
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a string.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
CHAR* symbol_string(const SYMBOL* symbol)
{
  CHAR* ret = NULL;

  if (symbol->type == sym_string)
  {
    ret = symbol->s.string;
  }
  else
  {
    dserror("Wrong symbol type %d", symbol->type);
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its an integer.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
INT symbol_int(const SYMBOL* symbol)
{
  INT ret = 0;

  if (symbol->type == sym_int)
  {
    ret = symbol->s.integer;
  }
  else
  {
    dserror("Wrong symbol type %d", symbol->type);
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a double.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
DOUBLE symbol_real(const SYMBOL* symbol)
{
  DOUBLE ret = 0.0;

  if (symbol->type == sym_real)
  {
    ret = symbol->s.real;
  }
  else
  {
    dserror("Wrong symbol type %d", symbol->type);
  }

  return ret;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Extract the value if its a map.

  \author u.kue
  \date 09/04
*/
/*----------------------------------------------------------------------*/
MAP* symbol_map(const SYMBOL* symbol)
{
  MAP* ret = NULL;

  if (symbol->type == sym_map)
  {
    ret = symbol->s.dir;
  }
  else
  {
    dserror("Wrong symbol type %d", symbol->type);
  }

  return ret;
}


#define PRINT_INDENT(f, count)                     \
  {                                                \
    INT i;                                         \
    for (i = 0; i < (count); ++i) fprintf(f, " "); \
  }


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
void symbol_print(FILE* f, CHAR* key, SYMBOL* symbol, INT indent)
{
  dsassert(strstr(key, "\n") == NULL, "malformed key");

  if (symbol != NULL)
  {
    symbol_print(f, key, symbol->next, indent);

    PRINT_INDENT(f, indent);
    switch (symbol->type)
    {
      case sym_string:
        fprintf(f, "%s = \"%s\"\n", key, symbol->s.string);
        break;
      case sym_int:
        fprintf(f, "%s = %d\n", key, symbol->s.integer);
        break;
      case sym_real:
      {
        /* We don't care for beauty but we need many digits. */
        fprintf(f, "%s = %.20e\n", key, symbol->s.real);
        break;
      }
      case sym_map:
        fprintf(f, "%s:\n", key);
        map_print(f, symbol->s.dir, indent + 4);
        break;
      default:
        dserror("Ups!");
    }
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Write this node and its subnodes.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void map_node_print(FILE* f, MAP_NODE* node, INT indent)
{
  if (node != NULL)
  {
    map_node_print(f, node->lhs, indent);
    symbol_print(f, node->key, node->symbol, indent);
    map_node_print(f, node->rhs, indent);
  }
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
  map_node_print(f, map->root.rhs, indent);
  fprintf(f, "\n");
}



/*----------------------------------------------------------------------*/
/*!
  \brief The types of tokens recognized by the lexer.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
typedef enum _TOKEN_TYPE
{
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
typedef struct _PARSER_DATA
{
  TOKEN_TYPE tok;
  CHAR* token_string;
  INT token_int;
  DOUBLE token_real;

  CHAR* file_buffer;
  INT file_size;
  /*CHAR* filename;*/

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
static void init_parser_data(struct _PARSER_DATA* data, const CHAR* filename, MPI_Comm comm)
{
  data->tok = tok_none;
  data->lineno = 1;
  data->pos = 0;
  data->indent_level = 0;
  data->indent_step = -1;

  int myrank = 0;
  int nprocs = 1;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nprocs);

  /* No copy here. Valid only as long as the calling functions
   * filename is valid. */
  /*data->filename = filename;*/

  /* We need to have the information on all processes. That's why we
   * read the file on process 0 and broadcast it. The other way would
   * be to use MPI IO, but then we'd have to implement a separat
   * sequential version. */
  if (myrank == 0)
  {
    INT bytes_read;
    FILE* file;
    file = fopen(filename, "rb");

    if (file == NULL)
    {
      dserror("cannot read file '%s'", filename);
    }

    /* find out the control file size */
    fseek(file, 0, SEEK_END);
    data->file_size = ftell(file);

    /* read file to local buffer */
    data->file_buffer = CCAMALLOC((data->file_size + 1) * sizeof(CHAR));
    fseek(file, 0, SEEK_SET);
    /*bytes_read = fread(data->file_buffer, sizeof(CHAR), data->file_size, file);*/
    bytes_read = fread(data->file_buffer, sizeof(CHAR), (size_t)data->file_size, file);
    if (bytes_read != data->file_size)
    {
      dserror("failed to read file %s", filename);
    }
    /* a trailing zero helps a lot */
    data->file_buffer[data->file_size] = '\0';

    fclose(file);
  }

  if (nprocs > 1)
  {
    INT err;
    err = MPI_Bcast(&data->file_size, 1, MPI_INT, 0, comm);
    if (err != 0)
    {
      dserror("MPI_Bcast failed: %d", err);
    }
    if (myrank > 0)
    {
      data->file_buffer = CCAMALLOC((data->file_size + 1) * sizeof(CHAR));
    }
    err = MPI_Bcast(data->file_buffer, data->file_size + 1, MPI_CHAR, 0, comm);
    if (err != 0)
    {
      dserror("MPI_Bcast failed: %d", err);
    }
  }
}


/*----------------------------------------------------------------------*/
/*!
  \brief Clean up.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
static void destroy_parser_data(struct _PARSER_DATA* data) { CCAFREE(data->file_buffer); }


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
  if (data->pos < data->file_size)
  {
    /* ignore dos line endings */
    if (data->file_buffer[data->pos] == '\r')
    {
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

  for (;;)
  {
    t = getnext(data);
    if (t == ' ')
    {
      /* ignore whitespaces */
      if (line_begin)
      {
        indention++;
      }
    }
    else if (t == '\t')
    {
      /* ignore whitespaces */
      if (line_begin)
      {
        indention = ((indention + TABWIDTH - 1) / TABWIDTH) * TABWIDTH;
      }
    }
    else if (t == '\n')
    {
      data->lineno++;
      line_begin = 1;
      indention = 0;
    }
    else if (t == '#')
    {
      for (;;)
      {
        t = getnext(data);
        if (t == '\n')
        {
          break;
        }
      }
      data->lineno++;
      line_begin = 1;
      indention = 0;
    }
    else if (t == EOF)
    {
      data->tok = tok_done;
      goto end;
    }
    else
    {
      if (line_begin && (indention != data->indent_level))
      {
        if (data->indent_step == -1)
        {
          if (indention > data->indent_level)
          {
            dsassert(data->indent_level == 0, "non-zero intention at first line?!");
            data->indent_step = indention;
            data->indent_level = indention;
            data->token_int = 1;
            data->tok = tok_indent;
          }
          else
          {
            dserror("dedent at toplevel!");
          }
        }
        else
        {
          if (indention > data->indent_level)
          {
            data->tok = tok_indent;
            dsassert(
                (indention - data->indent_level) % data->indent_step == 0, "malformed indention");
            data->token_int = (indention - data->indent_level) / data->indent_step;
            data->indent_level = indention;
          }
          else
          {
            data->tok = tok_dedent;
            dsassert(
                (data->indent_level - indention) % data->indent_step == 0, "malformed dedention");
            data->token_int = (data->indent_level - indention) / data->indent_step;
            data->indent_level = indention;
          }
        }
        data->pos--;
        goto end;
      }
      else
      {
        line_begin = 0;

        if ((t == '-') || isdigit(t))
        {
          data->token_string = &(data->file_buffer[data->pos - 1]);
          if (t == '-')
          {
            t = getnext(data);
          }
          while (isdigit(t))
          {
            t = getnext(data);
          }
          if ((t != '.') && (t != 'E') && (t != 'e'))
          {
            if (t != EOF)
            {
              data->pos--;
            }
            data->token_int = atoi(data->token_string);
            data->tok = tok_int;
            goto end;
          }
          if (t == '.')
          {
            t = getnext(data);
            if (isdigit(t))
            {
              while (isdigit(t))
              {
                t = getnext(data);
              }
            }
            else
            {
              dserror("no digits after point at line %d", data->lineno);
            }
          }
          if ((t == 'E') || (t == 'e'))
          {
            t = getnext(data);
            if ((t == '-') || (t == '+'))
            {
              t = getnext(data);
            }
            if (isdigit(t))
            {
              while (isdigit(t))
              {
                t = getnext(data);
              }
            }
            else
            {
              dserror("no digits after exponent at line %d", data->lineno);
            }
          }
          if (t != EOF)
          {
            data->pos--;
          }
          data->token_real = strtod(data->token_string, NULL);
          data->tok = tok_real;
          goto end;
        }
        else if (isalpha(t) || (t == '_'))
        {
          data->token_string = &(data->file_buffer[data->pos - 1]);
          while (isalnum(t) || (t == '_'))
          {
            t = getnext(data);
          }
          if (t != EOF)
          {
            data->pos--;
          }
          data->tok = tok_name;
          data->token_int = &(data->file_buffer[data->pos]) - data->token_string;
          goto end;
        }
        else if (t == '"')
        {
          data->token_string = &(data->file_buffer[data->pos]);
          t = getnext(data);
          while (t != '"')
          {
            t = getnext(data);
            if (t == EOF)
            {
              dserror("expected closing \" on line %d", data->lineno);
            }
          }
          data->tok = tok_string;
          data->token_int = &(data->file_buffer[data->pos - 1]) - data->token_string;
          goto end;
        }
        else if (t == ':')
        {
          data->tok = tok_colon;
          goto end;
        }
        else if (t == '=')
        {
          data->tok = tok_equal;
          goto end;
        }
        else
        {
          if (t >= 32)
            dserror("unexpected char '%c' at line %d", t, data->lineno);
          else
            dserror("unexpected char '%d' at line %d", t, data->lineno);
          data->tok = tok_none;
          goto end;
        }
      }
    }
  }

end:


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
  lexan(data);

  while (data->tok != tok_done)
  {
    switch (data->tok)
    {
      case tok_name:
      {
        CHAR* name;

        /*
         * The string is not null terminated as it's a simple pointer
         * into the file buffer. However, we know its length so we can
         * handle that. */
        name = CCAMALLOC((data->token_int + 1) * sizeof(CHAR));
        /*strncpy(name, data->token_string, data->token_int);*/
        strncpy(name, data->token_string, (size_t)data->token_int);
        name[data->token_int] = '\0';

        lexan(data);
        switch (data->tok)
        {
          case tok_colon:
          {
            MAP* map;

            lexan(data);
            if ((data->tok != tok_indent) || (data->token_int != 1))
            {
              dserror("Syntaxerror at line %d: single indention expected", data->lineno);
            }

            map = CCACALLOC(1, sizeof(MAP));
            init_map(map);
            parse_definitions(data, map);

            map_insert_map(dir, map, name);

            if ((data->tok == tok_dedent) && (data->token_int > 0))
            {
              data->token_int--;
              goto end;
            }

            break;
          }
          case tok_equal:
            lexan(data);
            switch (data->tok)
            {
              case tok_string:
              {
                CHAR* string;

                /* Again, be carefully with those pointers... */
                string = CCAMALLOC((data->token_int + 1) * sizeof(CHAR));
                /*strncpy(string, data->token_string, data->token_int);*/
                strncpy(string, data->token_string, (size_t)data->token_int);
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
  return;
}


/*----------------------------------------------------------------------*/
/*!
  \brief Parse the file given by name and fill the map with this
  file's content (serial only!)

  \author lw
  \date 05/08
*/
/*----------------------------------------------------------------------*/
void parse_control_file_serial(MAP* map, const CHAR* filename)
{
  PARSER_DATA data;

  /*
   * So here we are. Before the symbol table can be filled with values
   * it has to be initialized. That is we expect to get an
   * uninitialized (virgin) map. */
  init_map(map);

  data.tok = tok_none;
  data.lineno = 1;
  data.pos = 0;
  data.indent_level = 0;
  data.indent_step = -1;

  INT bytes_read;
  FILE* file;
  file = fopen(filename, "rb");

  if (file == NULL)
  {
    dserror("cannot read file '%s'", filename);
  }

  /* find out the control file size */
  fseek(file, 0, SEEK_END);
  data.file_size = ftell(file);

  /* read file to local buffer */
  data.file_buffer = CCAMALLOC((data.file_size + 1) * sizeof(CHAR));
  fseek(file, 0, SEEK_SET);
  /*bytes_read = fread(data.file_buffer, sizeof(CHAR), data.file_size, file);*/
  bytes_read = fread(data.file_buffer, sizeof(CHAR), (size_t)data.file_size, file);
  if (bytes_read != data.file_size)
  {
    dserror("failed to read file %s", filename);
  }
  /* a trailing zero helps a lot */
  data.file_buffer[data.file_size] = '\0';

  fclose(file);

  parse_definitions(&data, map);
  destroy_parser_data(&data);
}


/*----------------------------------------------------------------------*/
/*!
  \brief Parse the file given by name and fill the map with this
  file's content.

  \author u.kue
  \date 08/04
*/
/*----------------------------------------------------------------------*/
void parse_control_file(MAP* map, const CHAR* filename, MPI_Comm comm)
{
  PARSER_DATA data;

  /*
   * So here we are. Before the symbol table can be filled with values
   * it has to be initialized. That is we expect to get an
   * uninitialized (virgin) map. */
  init_map(map);

  init_parser_data(&data, filename, comm);
  parse_definitions(&data, map);
  destroy_parser_data(&data);
}
