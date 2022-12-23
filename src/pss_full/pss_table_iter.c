/*----------------------------------------------------------------------*/
/*! \file
\brief Iterators for the symbol table.


\level 1

Iterate the symbol table and visit all entries.

---------------------------------------------------------------------*/

#include "pss_table_iter.h"
#include "pss_prototypes.h"
#include "dserror.H"

/*----------------------------------------------------------------------*/
/*!
  \brief map iterator constructor

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
void init_map_iterator(MAP_ITERATOR* iterator, MAP* map)
{
  iterator->stack.count = 0;
  iterator->map = map;
  iterator->stack.head.map_node = NULL;
  iterator->stack.head.snext = NULL;
}

/*----------------------------------------------------------------------*/
/*!
  \brief map iterator push

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
static void push_map_node(MAP_ITERATOR* iterator, MAP_NODE* map_node)
{
  STACK_ELEMENT* new_element;

  new_element = (STACK_ELEMENT*)CCACALLOC(1, sizeof(STACK_ELEMENT));
  new_element->map_node = map_node;
  new_element->snext = iterator->stack.head.snext;
  iterator->stack.head.snext = new_element;
  iterator->stack.count++;
}

/*----------------------------------------------------------------------*/
/*!
  \brief map iterator pop

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
static void pop_map_node(MAP_ITERATOR* iterator)
{
  STACK_ELEMENT* tmp_free;

  if (iterator->stack.count == 0)
  {
    dserror("map iterator stack empty");
  }
  else
  {
    tmp_free = iterator->stack.head.snext;
    iterator->stack.head.snext = iterator->stack.head.snext->snext;
    iterator->stack.count--;
    CCAFREE(tmp_free);
  }
}

/*----------------------------------------------------------------------*/
/*!
  \brief map iterator

  \param iterator (i/o) the map iterator to be advanced
  \return true if a new node was found

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
INT next_map_node(MAP_ITERATOR* iterator)
{
  INT result = 0;

  /* if the map is empty there is nothing to iterate */
  if (iterator->map != NULL)
  {
    /*first call of this iterator*/
    if (iterator->stack.head.map_node == NULL)
    {
      /* we actually dont need the map->root information, we just use it
       * to show that the iterator is finally initalized*/
      iterator->stack.head.map_node = &iterator->map->root;

      if (iterator->map->root.rhs != NULL) push_map_node(iterator, iterator->map->root.rhs);
      if (iterator->map->root.lhs != NULL) push_map_node(iterator, iterator->map->root.lhs);

      /*if iterator is still empty return 0*/
      result = iterator->stack.head.snext != NULL;
    }
    else
    {
      if (iterator->stack.head.snext != NULL)
      {
        MAP_NODE* tmp;
        MAP_NODE* lhs;
        MAP_NODE* rhs;

        /* we remove the first member of the stack and add his rhs and lhs */
        tmp = iterator->stack.head.snext->map_node;
        lhs = tmp->lhs;
        rhs = tmp->rhs;
        tmp = NULL;

        /* caution! tmp is freed at this point! */
        pop_map_node(iterator);

        if (rhs != NULL) push_map_node(iterator, rhs);
        if (lhs != NULL) push_map_node(iterator, lhs);

        /*if iterator is empty now return 0*/
        result = iterator->stack.head.snext != NULL;
      }
    }
  }
  return result;
}


/*----------------------------------------------------------------------*/
/*!
  \brief map iterator current node

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
MAP_NODE* iterator_get_node(MAP_ITERATOR* iterator) { return iterator->stack.head.snext->map_node; }
