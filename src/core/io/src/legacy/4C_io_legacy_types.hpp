/*----------------------------------------------------------------------*/
/*! \file
\brief A very simple symbol table implementation.


\level 1

---------------------------------------------------------------------*/
/*!

Control files describe binary output. The control files are meant to be human readable and very
flexible but still must be read back. To accomplish this we need a small parser and a way to
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

#ifndef FOUR_C_IO_LEGACY_TYPES_HPP
#define FOUR_C_IO_LEGACY_TYPES_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

/*!
  \brief Types of the symbols we store in a map.

  We can have strings, integer and doubles as well as a subgroup.

  \author u.kue
  \date 08/04
*/
enum SymbolTypes
{
  sym_undefined,
  sym_string,
  sym_int,
  sym_real,
  sym_map
};


/*!
  \brief Data struct of a map's symbol.

  Depending on its type a symbol contains different values. You don't
  need to deal with these explicitly, there are access
  functions. However, you need to access the ``next`` pointer.

  Inside a map all symbols with the same key are linked by their next
  pointers. Thus a map node gives the first symbol and from there you
  can follow the next pointer until you reach a nullptr.

  \warning String and map values are assumed to be allocated
  dynamically and are owned by the symbol. That is when a symbol is
  destroyed its string or map value will be free'd, too.

  \author u.kue
  \date 08/04
*/
struct SYMBOL
{
  /* put the union first to have the double value properly aligned */
  union
  {
    char* string;
    int integer;
    double real;
    struct MAP* dir;
  } s;
  SymbolTypes type;
  struct SYMBOL* next;
};


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
struct MapNode
{
  char* key;
  SYMBOL* symbol;
  int count;
  struct MapNode* lhs;
  struct MapNode* rhs;
};


/*!
  \brief The central map structure.

  You need to initialize an object of this structure in order to have
  a map. Internally it's nothing but a map node and a total count.

  \author u.kue
  \date 08/04
*/
struct MAP
{
  MapNode root;
  int count{};
};

FOUR_C_NAMESPACE_CLOSE

#endif
