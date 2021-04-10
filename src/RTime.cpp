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

// define Inf
double const inf = R_PosInf;

// [[Rcpp::export]]
NumericMatrix HierarCluster_paris_time(Eigen::SparseMatrix<double> m) {
  Graph G(m);
  int n = G.adj.size();
  // cluster sizes
  IVector s(n, 1);

  // connected components
  std::vector<IVector> cc;
  // dendrogram as list of merges
  std::vector<NVector> D;
  IVector chain_record;
  NVector depth_record;

  int u = n;
  int record_a = 0;
  while (n > 0) {
    IVector chain;
    bool flag_single = true;
    for (int i = 0; i < G.adj.size(); i++) {
      if (!G.adj[i].empty()) {
        chain.push_back(G.adj[i].begin()->first);
        flag_single = false;
        break;
      }
    }
    if (flag_single) {
      for (int i = 0; i < G.adj.size(); i++) {
        if (G.wn[i] > 0) {
          chain.push_back(i);
          break;
        }
      }
    }
    record_a++;

    while (chain.size() != 0) {
      int a = chain[chain.size() - 1];
      std::vector<int>::iterator itera = chain.begin();
      chain.erase(itera + chain.size() - 1);
      double dmin = inf;
      int b = -1;
      for (std::unordered_map<int, double>::iterator iter = G.adj[a].begin(); iter != G.adj[a].end(); iter++) {
        int v = iter->first;
        double w_va = iter->second;
        if(w_va == 0){
          stop("There is 0 weight!");
        }
        double d = G.wn[v] * G.wn[a] / w_va / G.wtot;
        if (d < dmin) {
          b = v;
          dmin = d;
        }
        else if (d == dmin) {
          b = (b <= v) ? b : v;
        }
      }
      record_a++;

      double d = dmin;
      if (chain.size() != 0) {
        int c = chain[chain.size() - 1];
        std::vector<int>::iterator iterc = chain.begin();
        chain.erase(iterc + chain.size() - 1);
        if (b == c) {
          // merge a, b
          double a_double = static_cast<double>(a);
          double b_double = static_cast<double>(b);
          int s_apb = s[a] + s[b];
          double s_apb_double = static_cast<double>(s_apb);
          // double tempa[4] = { a, b, d, s[a] + s[b] };
          double tempa[4] = { a_double, b_double, d, s_apb_double };
          NVector tempnv(tempa, tempa + 4);
          D.push_back(tempnv);
          // std::cout << "acc: " << tempnv[0] << " b: " << tempnv[1] << " d: " << tempnv[2] << " scc: " << tempnv[3] << "\n";
          // update graph
          G.merge_node(a, b, u);
          n--;
          // update weight and size
          G.wn.push_back(G.wn[a] + G.wn[b]);
          G.wn[a] = 0;
          G.wn[b] = 0;
          s.push_back(s[a] + s[b]);
          s[a] = 0;
          s[b] = 0;
          u++;

          chain_record.push_back(record_a);
          depth_record.push_back(G.average_edge());
          record_a = 0;
        }
        else {
          chain.push_back(c);
          chain.push_back(a);
          chain.push_back(b);
        }
      }
      else if (b >= 0) {
        chain.push_back(a);
        chain.push_back(b);
      }
      else {
        // remove the connected component
        int tempa[2] = { a, s[a] };
        IVector tempiv(tempa, tempa + 2);
        cc.push_back(tempiv);
        G.del_node(a);
        G.wn[a] = 0;
        s[a] = 0;
        n--;
      }
    }
  }

  NumericMatrix output(2,chain_record.size());
  for(int i=0; i<chain_record.size(); i++){
    output(0,i) = chain_record[i];
  }
  for(int i=0; i<chain_record.size(); i++){
    output(1,i) = depth_record[i];
  }

  return(output);
}
