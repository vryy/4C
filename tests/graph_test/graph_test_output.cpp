
#include "graph_test.H"

#ifdef DEBUG
void print_cycle(const cycle_t* c)
{
  std::cout << "    ";
  std::copy(c->begin(), c->end(), std::ostream_iterator<vertex_t>(std::cout, " "));
  std::cout << "\n";
}
#endif
