
#include <fstream>
#include <iostream>
#include <string>

#include "graph_test.H"

void read_coords(const std::string& name, std::map<int, std::vector<double>>& coords)
{
  std::ifstream file(name.c_str());
  if (not file)
  {
    throw std::runtime_error("no point file");
  }

  double x, y, z;
  int n;
  while (file >> x >> y >> z >> n)
  {
    std::vector<double>& c = coords[n];
    c.reserve(3);
    c.push_back(x);
    c.push_back(y);
    c.push_back(z);
  }
}

void read_graph_data(const std::string& name, std::map<int, std::set<int>>& graph_data)
{
  std::ifstream file(name.c_str());
  if (not file)
  {
    throw std::runtime_error("no graph file");
  }

  std::set<int>* current = NULL;

  std::string word;
  while (file >> word)
  {
    if (word[word.length() - 1] == ':')
    {
      word = word.substr(0, word.length() - 1);
      int n = std::atoi(word.c_str());
      current = &graph_data[n];
    }
    else
    {
      int n = std::atoi(word.c_str());
      current->insert(n);
    }
  }
}

void construct_graph(Graph& g, std::map<int, std::set<int>>& graph_data)
{
  name_map_t name_map = boost::get(boost::vertex_name, g);

  std::map<int, vertex_t> vertex_map;

  for (std::map<int, std::set<int>>::iterator i = graph_data.begin(); i != graph_data.end(); ++i)
  {
    int n = i->first;

    vertex_t u = add_vertex(g);
    name_map[u] = n;
    vertex_map[n] = u;
  }

  for (std::map<int, std::set<int>>::iterator i = graph_data.begin(); i != graph_data.end(); ++i)
  {
    int u = i->first;
    std::set<int>& row = i->second;

    for (std::set<int>::iterator i = row.begin(); i != row.end(); ++i)
    {
      int v = *i;
      if (u < v)
      {
        boost::add_edge(vertex_map[u], vertex_map[v], g);
      }
    }
  }
}
