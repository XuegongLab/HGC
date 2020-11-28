#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace RcppEigen;

typedef std::vector<int> IVector;
typedef std::vector<double> NVector;

// global functions
NumericMatrix reorder_dendrogram(std::vector<NVector> D);

class Edge {
public:
  // members
  int node1;
  int node2;
  double weight;

  // init Edge
  Edge();
  Edge(int n1, int n2, double w);

  // functions
  // the member operator <
  bool operator <(const Edge& e) const;
  // the member operator ==
  bool operator ==(const Edge& e) const;
};

class Graph {
public:
  // adjacency matrix
  std::vector<std::unordered_map<int, double> > adj;
  // weight of nodes
  NVector wn;
  // total weight
  double wtot;

  //init
  Graph();
  Graph(Eigen::SparseMatrix<double> m);
  // Graph(NumericMatrix m);

  // functions
  // add new node to the back
  void add_node(int node);
  // add Edge
  void add_edge(Edge e);
  // get neighbors of a node
  IVector get_neighbors(int node);
  // get Edge weight
  double get_edge_weight(int node1, int node2);
  // add edge weight
  // void add_edge_weight(int node1, int node2, double weight);
  // del node
  void del_node(int node);
  // merge two nodes and add a new node
  void merge_node(int node1, int node2, int new_node);
  // used to observe the iteration progress
  double average_edge();
};

class Order_pair {
public:
  int order;
  double dist;

  Order_pair();
  Order_pair(int a, double b);

  // the member operator <
  bool operator <(const Order_pair& op) const;
};
