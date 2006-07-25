#ifndef POST_MAP_FUNCTIONS_H
#define POST_MAP_FUNCTIONS_H

#include "../post_common/post_common.h"

typedef struct _STACK_ELEMENT
{
  struct _STACK_ELEMENT* snext;
  MAP_NODE* map_node;

}STACK_ELEMENT;

typedef struct _STACK
{
  INT count;
  STACK_ELEMENT head;

}STACK;

typedef struct _MAP_ITERATOR
{
  MAP* map;
  STACK stack;

}MAP_ITERATOR;

typedef struct _MAP_NODE_ITERATOR
{
  MAP_NODE* map_node;
  SYMBOL* symbol;
}MAP_NODE_ITERATOR;

void stack_push_map_node(STACK* stack,
                        MAP_NODE* map_node);


void init_map_iterator(MAP_ITERATOR* iterator, MAP* map);

void stack_pop_map_node(STACK* stack);

INT next_map_node(MAP_ITERATOR* iterator);

void init_map_node_iterator(MAP_NODE_ITERATOR* iterator, MAP_NODE* map_node);

INT next_symbol(MAP_NODE_ITERATOR* iterator);

DOUBLE result_map_find_key(MAP* map, CHAR* key);

DOUBLE result_map_read_key(MAP* map,CHAR* key);

MAP_NODE* post_map_find_node(MAP* map, CHAR* key);

INT post_map_cmp_nodes(MAP_NODE* lhs, CHAR* rhs_key);

MAP_NODE* iterator_get_node(MAP_ITERATOR* iterator);
INT iterator_find_symbol(MAP_ITERATOR* iterator, char* name);
MAP* iterator_get_map(MAP_ITERATOR* iterator);
INT node_iterator_get_real_as_float(MAP_NODE_ITERATOR* node_iterator, float* real);


#endif
