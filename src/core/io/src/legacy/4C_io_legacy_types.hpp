// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
