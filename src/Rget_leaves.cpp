// [[Rcpp::depends(Rcpp)]]
#include <vector>
#include <stack>

#include <Rcpp.h>

using namespace Rcpp;

typedef std::vector<int> IVector;
typedef std::vector<double> NVector;

// print the leaves of the binary tree with preorder
// [[Rcpp::export]]
IntegerVector get_leaves(NumericMatrix hglink) {
  int n = hglink.nrow() + 1;
  int i = 2 * n - 2;
  std::vector<int> result;
  std::stack<int> sitestack;

  while (i >= n || !sitestack.empty())
  {
    while (i >= n)
    {
      sitestack.push(i);
      i = hglink(i - n, 0);
    }
    result.push_back(i);
    if (!sitestack.empty())
    {
      i = sitestack.top();
      sitestack.pop();
      i = hglink(i - n, 1);
    }
  }
  result.push_back(i);

  IntegerVector output = wrap(result);
  return(output);
}
