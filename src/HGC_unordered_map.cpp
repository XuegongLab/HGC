#include <iostream>
#include <vector>
#include <set>
#include <unordered_map>
#include <algorithm>
#include <RcppEigen.h>
#include <Rcpp.h>

#include "HGC_unordered_map.h"

using namespace Rcpp;
using namespace RcppEigen;

NumericMatrix reorder_dendrogram(std::vector<NVector> D){
  int n = D.size() + 1;
  std::vector<Order_pair> order;
  for (int i = 0; i < n - 1; i++) {
    order.push_back(Order_pair(i, D[i][2]));
  }
  std::sort(order.begin(), order.end());

  std::vector<IVector> nindex;
  IVector inputiv;
  for (int i = 0; i < 2 * n - 1; i++) {
    inputiv.push_back(i);
  }
  nindex.push_back(inputiv);
  nindex.push_back(inputiv);

  // nindex[0]: key
  // nindex[1]: value
  for (int i = n; i < 2 * n - 1; i++) {
    int neworder = n + order[i - n].order;
    unsigned int neworder_abs = abs(neworder);
    nindex[1][neworder_abs] = i;
    //std::cout << "nindex " << neworder << " " << i << std::endl;
  }

  std::vector<NVector> mlarge;
  for (int i = 0; i < n - 1; i++) {
    int a1 = D[i][0];
    int a2 = D[i][1];
    unsigned int a1_abs = abs(a1);
    unsigned int a2_abs = abs(a2);

    double tempa_ele0 = static_cast<double>(nindex[1][a1_abs]);
    double tempa_ele1 = static_cast<double>(nindex[1][a2_abs]);
    double tempa_ele2 = static_cast<double>(D[i][2]);
    double tempa_ele3 = static_cast<double>(D[i][3]);
    //double tempa[4] = { nindex[1][a1_abs], nindex[1][a2_abs],D[i][2], D[i][3] };
    double tempa[4] = { tempa_ele0, tempa_ele1, tempa_ele2, tempa_ele3 };
    NVector tempnv(tempa, tempa + 4);
    mlarge.push_back(tempnv);
  }

  std::vector<NVector> m;
  for (int i = 0; i < n - 1; i++) {
    int index = order[i].order;
    unsigned int index_abs = abs(index);
    m.push_back(mlarge[index_abs]);
    // std::cout << m[i][0] << " " << m[i][1] << " " << m[i][2] << " " << m[i][3] << " " << std::endl;
  }

  NumericMatrix mout(n - 1, 4);
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < 4; j++) {
      mout(i, j) = m[i][j];
    }
  }

  return(mout);
}

Order_pair::Order_pair() {
  order = -1;
  dist = 0;
}

Order_pair::Order_pair(int a, double b) {
  order = a;
  dist = b;
}

// the member operator <
bool Order_pair::operator <(const Order_pair& op) const {
  if (dist < op.dist)
    return true;
  else if (dist == op.dist && order < op.order)
    return true;
  else
    return false;
}

// init Edge
Edge::Edge() {
  node1 = -1;
  node2 = -1;
  weight = 0;
}

Edge::Edge(int n1, int n2, double w) {
  node1 = n1;
  node2 = n2;
  weight = w;
}

// the member operator <
bool Edge::operator <(const Edge& e) const {
  if (weight > e.weight)
    return true;
  else if (weight == e.weight && node2 < e.node2)
    return true;
  else
    return false;
}

// the member operator ==
bool Edge::operator ==(const Edge& e) const {
  return(node2 == e.node2);
}

Graph::Graph() {
  wtot = 0;
}

Graph::Graph(Eigen::SparseMatrix<double> m) {
  wtot = 0;

  if (m.cols() == m.rows()) {
    int nodes_num = m.cols();

    // initialize adj
    std::unordered_map<int, double> temp_vg(128);
    for (int i = 0; i < nodes_num; i++) {
      adj.push_back(temp_vg);
      wn.push_back(0);
    }

    for (int k = 0; k < m.outerSize(); ++k) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(m, k); it; ++it) {
        if (it.col() >= it.row()) {
          continue;
        }

        add_edge(Edge(it.col(), it.row(), it.value()));
        wn[it.col()] = wn[it.col()] + it.value();
        wn[it.row()] = wn[it.row()] + it.value();
        wtot = wtot + 2 * it.value();
      }
    }

    for (int i = 0; i < nodes_num; i++) {
      if (wn[i] == 0)
        wn[i] = 1;
    } // check the single points

  }
  else {
    stop("The input matrix should be a symmetric matrix.");
  }
}

// add new node to the back
void Graph::add_node(int node) {
  if (node == adj.size()) {
    std::unordered_map<int, double> temp_vg(128);
    adj.push_back(temp_vg);
  }
  else {
    // Rcout << "node:" << node << "\n";
    // stop("add node error");
  }
}

// add Edge
void Graph::add_edge(Edge e) {
  if (e.node1 < adj.size() && e.node2 < adj.size()) {
    if (e.node1 == e.node2) {
      std::pair<int, double> tempin(e.node1, e.weight);
      adj[e.node1].emplace(tempin);
    }
    else {
      std::pair<int, double> tempin(e.node2, e.weight);
      adj[e.node1].emplace(tempin);
      std::pair<int, double> tempin2(e.node1, e.weight);
      adj[e.node2].emplace(tempin2);
    }
  }
  else {
    // Rcout << "node1:" << e.node1 << "\tnode2:" << e.node2 << "\n";
    // stop("edge node out of range");
  }
}

// get neighbors of a node
IVector Graph::get_neighbors(int node) {
  IVector neighbors;
  for (std::unordered_map<int, double>::iterator iter = adj[node].begin(); iter != adj[node].end(); iter++) {
    neighbors.push_back(iter->first);
  }
  return neighbors;
}

// get Edge weight
double Graph::get_edge_weight(int node1, int node2) {

  for (std::unordered_map<int, double>::iterator iter = adj[node1].begin(); iter != adj[node1].end(); iter++) {
    if (iter->first == node2)
      return iter->second;
  }
  // stop("get_edge_weight error: can't find the edge");
  return 0;
}

// del node
void Graph::del_node(int node) {
  for (std::unordered_map<int, double>::iterator iter = adj[node].begin(); iter != adj[node].end(); iter++) {
    int node2 = iter->first;
    std::unordered_map<int, double>::iterator iter2 = adj[node2].find(node);
    if (iter2 != adj[node2].end())
      adj[node2].erase(iter2);
  }
  adj[node].clear();
}

// merge two nodes and add a new node
void Graph::merge_node(int node1, int node2, int new_node) {
  add_node(new_node);

  for (std::unordered_map<int, double>::iterator iter = adj[node1].begin(); iter != adj[node1].end(); iter++) {
    if (iter->first != node2)
      add_edge(Edge(new_node, iter->first, iter->second));
  }

  for (std::unordered_map<int, double>::iterator iter = adj[node2].begin(); iter != adj[node2].end(); iter++) {
    if (iter->first == node1)
      continue;

    std::unordered_map<int, double>::iterator iter2 = adj[new_node].find(iter->first);
    if (iter2 != adj[new_node].end()) {
      // renew the weight in new_node's neighbor set
      iter2->second = iter2->second + iter->second;

      int node3 = iter->first;
      // renew the weight in node3's neighbor set
      std::unordered_map<int, double>::iterator iter3 = adj[node3].find(new_node);
      if (iter3 != adj[node3].end()) {
        iter3->second = iter2->second;
      }
    }
    else {
      add_edge(Edge(new_node, iter->first, iter->second));
    }
  }
  del_node(node1);
  del_node(node2);
}


// used to observe the iteration progress
double Graph::average_edge() {
  int node_num = 0;
  int node_edge = 0;
  for (std::vector<std::unordered_map<int, double> >::iterator iter = adj.begin(); iter != adj.end(); iter++) {
    if (iter->size() > 0) {
      node_num = node_num + 1;
      node_edge = node_edge + iter->size();
    }
  }
  if (node_num == 0) {
    return(0);
  }
  else {
    double avg_edge = node_edge / node_num;
    return(avg_edge);
  }
}
