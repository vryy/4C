
#include <iostream>
#include <string>

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#include "graph_test.H"

int main(int ac, char** av)
{
  // Declare the supported options.
  boost::program_options::options_description desc("Allowed options");
  desc.add_options()("help", "produce help message")("cmd",
      boost::program_options::value<std::string>()->default_value("find_cycles"),
      "set graph command")("graph",
      boost::program_options::value<std::string>()->default_value("graph.txt"), "set graph file")(
      "points", boost::program_options::value<std::string>()->default_value("all_points.plot"),
      "set point file");

  boost::program_options::variables_map vm;
  boost::program_options::store(boost::program_options::parse_command_line(ac, av, desc), vm);
  boost::program_options::notify(vm);

  if (vm.count("help"))
  {
    std::cout << desc << "\n";
    return 1;
  }

  std::map<int, std::vector<double>> coords;
  std::map<int, std::set<int>> graph_data;

  read_coords(vm["points"].as<std::string>(), coords);
  read_graph_data(vm["graph"].as<std::string>(), graph_data);

  Graph g;
  construct_graph(g, graph_data);

  std::string cmd = vm["cmd"].as<std::string>();

  if (cmd == "find_cycles")
  {
    std::set<cycle_t*> base_cycles;
    find_cycles(g, coords, base_cycles);
    cleanup(g, coords, base_cycles);
  }
  else
  {
    std::cout << "unknown command '" << cmd << "'\n";
  }

  return 0;
}
