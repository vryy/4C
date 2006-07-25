#include "post_map_functions.h"

void init_map_iterator(MAP_ITERATOR* iterator, MAP* map)
{
  iterator->stack.count=0;
  iterator->map=map;
  iterator->stack.head.map_node=NULL;
  iterator->stack.head.snext=NULL;
}

void push_map_node(MAP_ITERATOR* iterator, MAP_NODE* map_node)
{
  STACK_ELEMENT* new_element=(STACK_ELEMENT*)CCACALLOC(1, sizeof(STACK_ELEMENT));
  new_element->map_node=map_node;
  new_element->snext=iterator->stack.head.snext;
  iterator->stack.head.snext=new_element;
  iterator->stack.count++;
}

void pop_map_node(MAP_ITERATOR* iterator)
{
  STACK_ELEMENT* tmp_free;

  if (iterator->stack.count==0)
  {
  }
  else
  {
    tmp_free=iterator->stack.head.snext;
    iterator->stack.head.snext=iterator->stack.head.snext->snext;
    iterator->stack.count--;
    CCAFREE(tmp_free);
  }
}

INT next_map_node(MAP_ITERATOR* iterator)
{
  MAP_NODE* tmp;

  /*first call of this iterator*/
  if (iterator->stack.head.map_node==NULL)
  {
    /* we actually dont need the map->root information, we just use it
     * to show that the iterator is finally initalized*/
    iterator->stack.head.map_node=&iterator->map->root;

    if (iterator->map->root.rhs!=NULL)  push_map_node(iterator, iterator->map->root.rhs);
    if (iterator->map->root.lhs!=NULL)  push_map_node(iterator, iterator->map->root.lhs);

    /*if iterator is still empty return 0*/
    return (iterator->stack.head.snext!=NULL);
  }

  else
  {
    if (iterator->stack.head.snext!=NULL)
    {
      /* we remove the first member of the stack and add his rhs and lhs */
      tmp=iterator->stack.head.snext->map_node;

      pop_map_node(iterator);

      if (tmp->rhs!=NULL)  push_map_node(iterator, tmp->rhs);
      if (tmp->lhs!=NULL)  push_map_node(iterator, tmp->lhs);

      /*if iterator is empty now return 0*/
      return(iterator->stack.head.snext!=NULL);
    }

    else
      return 0;
  }
}

void init_map_node_iterator(MAP_NODE_ITERATOR* iterator, MAP_NODE* map_node)
{
  iterator->map_node=map_node;
  iterator->symbol=NULL;
}

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
    if (iterator->symbol->next==NULL) return 0;
    else
    {
      iterator->symbol=iterator->symbol->next;
      return 1;
    }
  }
}

DOUBLE result_map_read_key(MAP* map, CHAR* key)
{
  MAP_NODE* actnode;
  MAP_ITERATOR iterator;   /* used for stepping through all map_nodes*/
  MAP_NODE_ITERATOR node_iterator;

  init_map_iterator(&iterator, map);
  while (next_map_node(&iterator)) /*loop over all map_nodes*/
  {
    actnode=iterator.stack.head.snext->map_node;
    init_map_node_iterator(&node_iterator, actnode);
    while (next_symbol(&node_iterator)) /*loop over all symbols*/
    {
      if (node_iterator.symbol->type==4)
      {
        if (map_find_symbol(node_iterator.symbol->s.dir, key)!=NULL)
          return map_read_real(node_iterator.symbol->s.dir, key);
      }
    }
  }
  return 0;
}

MAP_NODE* post_map_find_node(MAP* map, CHAR* key)
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
    cmp = post_map_cmp_nodes(node, key);
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

INT post_map_cmp_nodes(MAP_NODE* lhs, CHAR* rhs_key)
{
  return strcmp(lhs->key, rhs_key);
}


INT symbol_get_real_as_float(SYMBOL* symbol, float* real)
{
  INT ret;

#ifdef DEBUG
  dstrc_enter("symbol_get_real_as_float");
#endif

  if (symbol && (symbol->type == sym_real)) {
    *real = (float)symbol->s.real;
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

MAP_NODE* iterator_get_node(MAP_ITERATOR* iterator)
{
  return iterator->stack.head.snext->map_node;
}

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
