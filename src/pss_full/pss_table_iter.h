/*!---------------------------------------------------------------------
\file pss_table_iter.h
\brief Iterators for the symbol table.

\maintainer Martin Kronbichler

\level 1

---------------------------------------------------------------------*/

#ifndef PSS_TABLE_ITER_H
#define PSS_TABLE_ITER_H

#include "../headers/standardtypes.h"
#include "pss_table.h"


/*----------------------------------------------------------------------*/
/*!
  \brief map node stack element

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
typedef struct _STACK_ELEMENT
{
  struct _STACK_ELEMENT* snext;
  MAP_NODE* map_node;
} STACK_ELEMENT;


/*----------------------------------------------------------------------*/
/*!
  \brief stack of map nodes

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
typedef struct _STACK
{
  INT count;
  STACK_ELEMENT head;
} STACK;


/*----------------------------------------------------------------------*/
/*!
  \brief map iterator

  Visit all maps inside a map. This is a tree iterator. The map is
  implemented as a tree. Hence there is a stack inside this iterator.

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
typedef struct _MAP_ITERATOR
{
  MAP* map;
  STACK stack;
} MAP_ITERATOR;


/*----------------------------------------------------------------------*/
/*!
  \brief map node iterator

  Visit all values of a given symbol.

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
typedef struct _MAP_NODE_ITERATOR
{
  MAP_NODE* map_node;
  SYMBOL* symbol;
} MAP_NODE_ITERATOR;


void init_map_iterator(MAP_ITERATOR* iterator, MAP* map);

INT next_map_node(MAP_ITERATOR* iterator);

void init_map_node_iterator(MAP_NODE_ITERATOR* iterator, MAP_NODE* map_node);

INT next_symbol(MAP_NODE_ITERATOR* iterator);

DOUBLE result_map_find_key(MAP* map, CHAR* key);

DOUBLE result_map_read_key(MAP* map, CHAR* key);

MAP_NODE* iterator_get_node(MAP_ITERATOR* iterator);
INT iterator_find_symbol(MAP_ITERATOR* iterator, char* name);
MAP* iterator_get_map(MAP_ITERATOR* iterator);
INT node_iterator_get_real_as_float(MAP_NODE_ITERATOR* node_iterator, float* real);

#endif
