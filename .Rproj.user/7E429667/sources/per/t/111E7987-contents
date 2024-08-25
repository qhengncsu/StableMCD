#pragma once
#include <time.h>
#include <algorithm>
#include <math.h>
#include <float.h>
#include <vector>
#include <set>
#include <stdlib.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random.hpp>
#include <random>
#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;


typedef vector<double> TPoint;
typedef vector<vector<double> > TMatrix;
typedef vector<vector<int> > TIntMatrix;
typedef vector<int> TVariables;

typedef double** TDMatrix;
typedef double*** T3DMatrix;

// by rows
TDMatrix asMatrix(double* arr, int n, int d);
T3DMatrix as3DMatrix(double* arr, int n, int t, int d);

double** newM(int n, int d);
void deleteM(TDMatrix X);
TDMatrix copyM(TDMatrix X, int n, int d);
void printMatrix(TDMatrix mat, int n, int d);

void GetDirections(TDMatrix directions, int k, int d);
void GetProjections(TDMatrix points, int n, int d, TDMatrix directions, int k, TDMatrix projections);

int GetDepthsPrj(TDMatrix points, int n, int d, TDMatrix objects, int m,
                 TVariables cardinalities, int k, bool newDirs,
                 TDMatrix depths, TDMatrix directions, TDMatrix projections);

struct OrderRec {
  int order;
  double value;
  OrderRec(int order = -1, double value = 0) {
    this->order = order;
    this->value = value;
  }
};

#define ran(x) rEngine()%x
#define setseed(x) rEngine.seed(x)