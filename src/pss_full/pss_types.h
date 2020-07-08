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

#ifndef PSS_TYPES_H
#define PSS_TYPES_H

#include "../headers/standardtypes.h"
#include "../headers/am.h"

/*!
  \brief Types of the symbols we store in a map.

  We can have strings, integer and doubles as well as a subgroup.

  \author u.kue
  \date 08/04
*/
typedef enum _SYMBOL_TYPES
{
  sym_undefined,
  sym_string,
  sym_int,
  sym_real,
  sym_map
} SYMBOL_TYPES;


/*!
  \brief Data struct of a map's symbol.

  Depending on its type a symbol contains different values. You don't
  need to deal with these explicitly, there are access
  functions. However, you need to access the ``next`` pointer.

  Inside a map all symbols with the same key are linked by their next
  pointers. Thus a map node gives the first symbol and from there you
  can follow the next pointer until you reach a NULL.

  \warning String and map values are assumed to be allocated
  dynamically and are owned by the symbol. That is when a symbol is
  destroyed its string or map value will be free'd, too.

  \author u.kue
  \date 08/04
*/
typedef struct _SYMBOL
{
  /* put the union first to have the double value properly aligned */
  union {
    CHAR* string;
    INT integer;
    DOUBLE real;
    struct _MAP* dir;
  } s;
  SYMBOL_TYPES type;
  struct _SYMBOL* next;
} SYMBOL;


/*!
  \brief A map node that contains a key and a symbol.

  The backbone structure of a map is the map node. Map nodes
  constitute a binary tree (unbalanced by now), so there are left and
  right subnodes. The symbol is the start of a list of symbols that
  share one key.

  The node knows how many symbols there are.

  \warning The map node owns its key and symbols as well as any
  subnodes. It'll destroy them when it is destroyed itself.

  \author u.kue
  \date 08/04
*/
typedef struct _MAP_NODE
{
  CHAR* key;
  SYMBOL* symbol;
  INT count;
  struct _MAP_NODE* lhs;
  struct _MAP_NODE* rhs;
} MAP_NODE;


/*!
  \brief The central map structure.

  You need to initialize an object of this structure in order to have
  a map. Internally it's nothing but a map node and a total count.

  \author u.kue
  \date 08/04
*/
typedef struct _MAP
{
  MAP_NODE root;
  INT count;
} MAP;


#endif
