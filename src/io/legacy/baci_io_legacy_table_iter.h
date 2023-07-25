/*----------------------------------------------------------------------*/
/*! \file
\brief Iterators for the symbol table.


\level 1

---------------------------------------------------------------------*/

#ifndef BACI_IO_LEGACY_TABLE_ITER_H
#define BACI_IO_LEGACY_TABLE_ITER_H

#include "baci_io_legacy_table.h"


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
  int count;
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


void init_map_iterator(MAP_ITERATOR* iterator, MAP* map);

int next_map_node(MAP_ITERATOR* iterator);

MAP_NODE* iterator_get_node(MAP_ITERATOR* iterator);

#endif
