#ifndef GRAPH_TEST_H
#define GRAPH_TEST_H

#include <fstream>
#include <iostream>

#include <map>
#include <set>
#include <string>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/planar_face_traversal.hpp>
#include <boost/graph/boyer_myrvold_planar_test.hpp>

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::property<boost::vertex_name_t, int,
        boost::property<boost::vertex_color_t, boost::default_color_type,
            boost::property<boost::vertex_index_t, int>>>,
    boost::property<boost::edge_index_t, int>>
    Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor vertex_t;
typedef boost::graph_traits<Graph>::edge_descriptor edge_t;

typedef boost::graph_traits<Graph>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<Graph>::edge_iterator edge_iterator;
typedef boost::graph_traits<Graph>::adjacency_iterator adjacency_iterator;

typedef boost::property_map<Graph, boost::vertex_name_t>::type name_map_t;
typedef boost::property_map<Graph, boost::vertex_color_t>::type color_map_t;
typedef boost::property_map<Graph, boost::vertex_index_t>::type vertex_index_map_t;
typedef boost::property_map<Graph, boost::edge_index_t>::type edge_index_map_t;

typedef boost::color_traits<boost::property_traits<color_map_t>::value_type> color_t;

typedef std::vector<vertex_t> cycle_t;

void read_coords(const std::string& name, std::map<int, std::vector<double>>& coords);
void read_graph_data(const std::string& name, std::map<int, std::set<int>>& graph_data);
void construct_graph(Graph& g, std::map<int, std::set<int>>& graph_data);

void gnuplot_graph(const std::string& name, Graph& g, std::map<int, std::vector<double>>& coords);
void gnuplot_cycle(
    const std::string& name, Graph& g, cycle_t* c, std::map<int, std::vector<double>>& coords);

#ifdef DEBUG
void print_cycle(const cycle_t* c);
#endif

void find_cycles(
    Graph& g, std::map<int, std::vector<double>>& coords, std::set<cycle_t*>& base_cycles);
void cleanup(Graph& g, std::map<int, std::vector<double>>& coords, std::set<cycle_t*>& base_cycles);

#endif
