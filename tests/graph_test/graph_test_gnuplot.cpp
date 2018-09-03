
#include "graph_test.H"

void gnuplot_line(std::ostream& stream, std::vector<double>& c)
{
  std::copy(c.begin(), c.end(), std::ostream_iterator<double>(stream, " "));
  stream << "\n";
}

void gnuplot_graph(const std::string& name, Graph& g, std::map<int, std::vector<double>>& coords)
{
  std::ofstream file(name.c_str());
  if (not file)
  {
    throw std::runtime_error("no output file");
  }

  name_map_t name_map = boost::get(boost::vertex_name, g);

  edge_iterator ei, ei_end;
  for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ++ei)
  {
    edge_t e = *ei;
    vertex_t u = boost::source(e, g);
    vertex_t v = boost::target(e, g);

    gnuplot_line(file, coords[name_map[u]]);
    gnuplot_line(file, coords[name_map[v]]);
    file << "\n\n";
  }
}

void gnuplot_cycle(
    const std::string& name, Graph& g, cycle_t* c, std::map<int, std::vector<double>>& coords)
{
  std::ofstream file(name.c_str());
  if (not file)
  {
    throw std::runtime_error("no output file");
  }

  name_map_t name_map = boost::get(boost::vertex_name, g);

  for (unsigned i = 0; i < c->size(); ++i)
  {
    gnuplot_line(file, coords[name_map[c->at(i)]]);
  }
  gnuplot_line(file, coords[name_map[c->at(0)]]);
  file << "\n\n";
}
