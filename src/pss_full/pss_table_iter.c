/*!
\file
\brief Iterators for the symbol table.

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>

Iterate the symbol table and visit all entries.

*/

#include "pss_table_iter.h"


/*----------------------------------------------------------------------*/
/*!
  \brief map iterator constructor

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
void init_map_iterator(MAP_ITERATOR* iterator, MAP* map)
{
  iterator->stack.count=0;
  iterator->map=map;
  iterator->stack.head.map_node=NULL;
  iterator->stack.head.snext=NULL;
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

#ifdef DEBUG
  dstrc_enter("push_map_node");
#endif

  new_element=(STACK_ELEMENT*)CCACALLOC(1, sizeof(STACK_ELEMENT));
  new_element->map_node=map_node;
  new_element->snext=iterator->stack.head.snext;
  iterator->stack.head.snext=new_element;
  iterator->stack.count++;

#ifdef DEBUG
  dstrc_exit();
#endif
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

#ifdef DEBUG
  dstrc_enter("pop_map_node");
#endif

  if (iterator->stack.count==0)
  {
    dserror("map iterator stack empty");
  }
  else
  {
    tmp_free=iterator->stack.head.snext;
    iterator->stack.head.snext=iterator->stack.head.snext->snext;
    iterator->stack.count--;
    CCAFREE(tmp_free);
  }

#ifdef DEBUG
  dstrc_exit();
#endif
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

#ifdef DEBUG
  dstrc_enter("next_map_node");
#endif

  /*first call of this iterator*/
  if (iterator->stack.head.map_node==NULL)
  {
    /* we actually dont need the map->root information, we just use it
     * to show that the iterator is finally initalized*/
    iterator->stack.head.map_node=&iterator->map->root;

    if (iterator->map->root.rhs!=NULL)  push_map_node(iterator, iterator->map->root.rhs);
    if (iterator->map->root.lhs!=NULL)  push_map_node(iterator, iterator->map->root.lhs);

    /*if iterator is still empty return 0*/
    result = iterator->stack.head.snext!=NULL;
  }
  else
  {
    if (iterator->stack.head.snext!=NULL)
    {
      MAP_NODE* tmp;
      MAP_NODE* lhs;
      MAP_NODE* rhs;

      /* we remove the first member of the stack and add his rhs and lhs */
      tmp=iterator->stack.head.snext->map_node;
      lhs = tmp->lhs;
      rhs = tmp->rhs;
      tmp = NULL;

      /* caution! tmp is freed at this point! */
      pop_map_node(iterator);

      if (rhs!=NULL)  push_map_node(iterator, rhs);
      if (lhs!=NULL)  push_map_node(iterator, lhs);

      /*if iterator is empty now return 0*/
      result = iterator->stack.head.snext!=NULL;
    }
  }

#ifdef DEBUG
  dstrc_exit();
#endif
  return result;
}

/*----------------------------------------------------------------------*/
/*!
  \brief map node iterator constructor

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
void init_map_node_iterator(MAP_NODE_ITERATOR* iterator, MAP_NODE* map_node)
{
  iterator->map_node=map_node;
  iterator->symbol=NULL;
}


/*----------------------------------------------------------------------*/
/*!
  \brief map node iterator next

  \param iterator (i/o)
  \return true if a new value has been found

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
INT next_symbol(MAP_NODE_ITERATOR* iterator)
{
  /*first call of this iterator*/
  if (iterator->symbol==NULL)
  {
    iterator->symbol=iterator->map_node->symbol;
    return 1;
  }
  else
  {
    if (iterator->symbol->next==NULL)
      return 0;
    else
    {
      iterator->symbol=iterator->symbol->next;
      return 1;
    }
  }
}

/*----------------------------------------------------------------------*/
/*!
  \brief

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
DOUBLE result_map_read_key(MAP* map, CHAR* key)
{
  MAP_ITERATOR iterator;   /* used for stepping through all map_nodes*/

#ifdef DEBUG
  dstrc_enter("result_map_read_key");
#endif

  init_map_iterator(&iterator, map);
  while (next_map_node(&iterator)) /*loop over all map_nodes*/
  {
    MAP_NODE_ITERATOR node_iterator;
    MAP_NODE* actnode;
    actnode=iterator.stack.head.snext->map_node;
    init_map_node_iterator(&node_iterator, actnode);
    while (next_symbol(&node_iterator)) /*loop over all symbols*/
    {
      if (node_iterator.symbol->type==sym_map)
      {
        if (map_find_symbol(node_iterator.symbol->s.dir, key)!=NULL)
        {
#ifdef DEBUG
          dstrc_exit();
#endif
          return map_read_real(node_iterator.symbol->s.dir, key);
        }
      }
    }
  }
#ifdef DEBUG
  dstrc_exit();
#endif
  return 0;
}

/*----------------------------------------------------------------------*/
/*!
  \brief map iterator current node

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
MAP_NODE* iterator_get_node(MAP_ITERATOR* iterator)
{
  return iterator->stack.head.snext->map_node;
}

/*----------------------------------------------------------------------*/
/*!
  \brief map iterator constructor current map

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
MAP* iterator_get_map(MAP_ITERATOR* iterator)
{
#ifdef DEBUG
  dstrc_enter("iterator_get_node");
#endif

  MAP_NODE* actnode = iterator_get_node(iterator);
  MAP* ret = NULL;

  if (actnode->symbol->type==sym_map)
    ret = actnode->symbol->s.dir;

  else
    dserror("iterator has no map symbol");

#ifdef DEBUG
  dstrc_exit();
#endif

  return ret;

}

/*----------------------------------------------------------------------*/
/*!
  \brief find a symbol in the current map

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
INT iterator_find_symbol(MAP_ITERATOR* iterator, char* name)
{
  int ret = 0;
#ifdef DEBUG
  dstrc_enter("iterator_get_node");
#endif

  if (map_find_symbol(iterator_get_map(iterator), name)!=NULL)
    ret = 1;

#ifdef DEBUG
  dstrc_exit();
#endif
  return ret;
}

/*----------------------------------------------------------------------*/
/*!
  \brief get the float value

  \warning This will always return the first value in the symbol
  chain. The iterator position is ignored.

  \author m.geppert (u.kue)
  \date 08/06
 */
/*----------------------------------------------------------------------*/
INT node_iterator_get_real_as_float(MAP_NODE_ITERATOR* node_iterator, float* real)
{
  INT ret;
#ifdef DEBUG
  dstrc_enter("iterator_get_node");
#endif

  ret = symbol_get_real_as_float(node_iterator->symbol, real);

#ifdef DEBUG
  dstrc_exit();
#endif

  return ret;
}
