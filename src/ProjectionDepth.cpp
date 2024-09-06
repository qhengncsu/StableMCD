#include "ProjectionDepth.h"
boost::random::rand48 rEngine;
boost::random::normal_distribution<double> normDist;
// 3D-array structures
T3DMatrix as3DMatrix(double* arr, int n, int t, int d){
  T3DMatrix mat = new double**[n];
  for (int i = 0; i < n; i++){
    mat[i] = new double*[t];
    for (int j = 0; j < t; j++)
    {
      mat[i][j] = arr + i*t*d + j*d;
    }
  }
  return mat;
}



// by rows
TDMatrix asMatrix(double* arr, int n, int d){
  TDMatrix mat = new double*[n];
  for (int i = 0; i < n; i++)
    mat[i] = arr + i*d;
  return mat;
}

double** newM(int n, int d){
  double* a = new double[n*d];
  return asMatrix(a, n, d);
}

void deleteM(TDMatrix X){
  delete[] X[0];
  delete[] X;
}

TDMatrix copyM(TDMatrix X, int n, int d){
  double* a = new double[n*d];
  memcpy(a, X[0], n*d*sizeof(double));
  return asMatrix(a, n, d);
}

void printMatrix(TDMatrix mat, int n, int d){
  for (int i = 0; i < n; i++){
    for (int j = 0; j < d; j++)
      Rcout << mat[i][j] << "\t";
    Rcout << endl;
  }
  Rcout << endl;
}

void GetDirections(TDMatrix directions, int k, int d, int dirStyle, 
                   TDMatrix points, int n){
  if (dirStyle == 1){
    for (int i = 0; i < k; i++){
      double* direction = directions[i];
      double sqrSum = 0;
      for (int j = 0; j < d; j++){
        direction[j] = normDist(rEngine);
        sqrSum += direction[j]*direction[j];
      }
      sqrSum = sqrt(sqrSum);
      for (int j = 0; j < d; j++){
        direction[j] = direction[j]/sqrSum;
      }
    }
  } else if(dirStyle == 2){
    std::uniform_int_distribution<int> dis(0, n-1);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = 0; i < k; i++){
      double* direction = directions[i];
      int index1 = dis(gen);
      int index2 = dis(gen);
      double sqrSum = 0;
      while (index1==index2){
        index2 = dis(gen);
      }
      double* row1 = points[index1];
      double* row2 = points[index2];
      for (int j = 0; j < d; j++){
        direction[j] = row2[j] - row1[j];
        sqrSum += direction[j]*direction[j];
      }
      sqrSum = sqrt(sqrSum);
      for (int j = 0; j < d; j++){
        direction[j] = direction[j]/sqrSum;
      }
    }
  } else if(dirStyle == 3){
    for (int i = 0; i < k/2; i++){
      double* direction = directions[i];
      double sqrSum = 0;
      for (int j = 0; j < d; j++){
        direction[j] = normDist(rEngine);
        sqrSum += direction[j]*direction[j];
      }
      sqrSum = sqrt(sqrSum);
      for (int j = 0; j < d; j++){
        direction[j] = direction[j]/sqrSum;
      }
    }
    std::uniform_int_distribution<int> dis(0, n-1);
    std::random_device rd;
    std::mt19937 gen(rd());
    for (int i = k/2; i < k; i++){
      double* direction = directions[i];
      int index1 = dis(gen);
      int index2 = dis(gen);
      double sqrSum = 0;
      while (index1==index2){
        index2 = dis(gen);
      }
      double* row1 = points[index1];
      double* row2 = points[index2];
      for (int j = 0; j < d; j++){
        direction[j] = row2[j] - row1[j];
        sqrSum += direction[j]*direction[j];
      }
      sqrSum = sqrt(sqrSum);
      for (int j = 0; j < d; j++){
        direction[j] = direction[j]/sqrSum;
      }
    }
  }
}

void GetProjections(TDMatrix points, int n, int d, TDMatrix directions, int k, TDMatrix projections){
  for (int i = 0; i < k; i++){
    double* projection = projections[i];
    for (int j = 0; j < n; j++){
      double sum = 0;
      for (int l = 0; l < d; l++){
        sum += points[j][l]*directions[i][l];
      }
      projection[j] = sum;
    }
  }
}

// static int CompareAsc(OrderRec x, OrderRec y)
// {
//   return (x.value < y.value);
// }
// 
// static int CompareDec(OrderRec x, OrderRec y)
// {
//   return (x.value > y.value);
// }

static void GetMedMad(TPoint &points, double &median, double &mad){
  /* First, determine median */
  int n = points.size();
  //sort(points.begin(), points.end());
  //median = (points[(n + 1)/2 - 1] + points[(n + 2)/2 - 1])/2.;
  nth_element(points.begin(), points.begin() + n/2, points.end());
  median = points[n/2];
  /* Obtain median absolute deviation (from median) (MAD) */
  TPoint deviations(n);
  for (int i = 0; i < n; i++){deviations[i] = abs(points[i] - median);}
  //sort(deviations.begin(), deviations.end());
  //median = (deviations[(n + 1)/2 - 1] + deviations[(n + 2)/2 - 1])/2.;
  nth_element(deviations.begin(), deviations.begin() + n/2, deviations.end());
  mad = deviations[n/2];
}

void GetPtsPrjDepths(double* projection, int n, double* objectsProjection, int m,
                     TVariables cardinalities, TMatrix *ptsPrjDepths){
  /* Collect basic statistics */
  int q = cardinalities.size();
  for (int i = 0; i < q; i++){
    /* Prepare data and obtain median and mad*/
    int beginIndex = 0;
    for (int j = 0; j < q; j++){
      if (j >= i){break;}
      beginIndex += cardinalities[j];
    }
    int endIndex = beginIndex + cardinalities[i];
    TPoint curClassProjection(projection + beginIndex, projection + endIndex);
    double median, mad;GetMedMad(curClassProjection, median, mad);
    /* Calculate i-class projectional univariate depths */
    for (int j = 0; j < m; j++){
      (*ptsPrjDepths)[i][j] = (objectsProjection[j] - median)/mad;
    }
  }
}

int GetDepthsPrj(TDMatrix points, int n, int d, TDMatrix objects, int m,
                 TVariables cardinalities, int k, bool newDirs, 
                 TDMatrix depths, TDMatrix directions, TDMatrix projections, int dirStyle){
  /* 1. Collecting basic statistics */
  int q = cardinalities.size();
  TDMatrix objectsProjections = newM(k,m);
  if (newDirs){
    GetDirections(directions, k, d, dirStyle, points, n);
    GetProjections(points, n, d, directions, k, projections);
  }
  GetProjections(objects, m, d, directions, k, objectsProjections);
  /* 2. Calculate projection depths */
  vector<vector<vector<double> > > prjDepths(k, vector<vector<double> >(q, vector<double > (m)));
  for (int i = 0; i < k; i++){
    GetPtsPrjDepths(projections[i], n, objectsProjections[i], m, cardinalities,
                    &prjDepths[i]);
  }
  /* 3. Merge depths */
  for (int i = 0; i < m; i++){
    for (int j = 0; j < q; j++){
      depths[i][j] = DBL_MIN;
    }
  }
  for (int i = 0; i < k; i++){
    for (int j = 0; j < q; j++){
      for (int l = 0; l < m; l++){
        if (prjDepths[i][j][l] > depths[l][j]){
          depths[l][j] = prjDepths[i][j][l];
        }
      }
    }
  }
  for (int i = 0; i < m; i++){
    for (int j = 0; j < q; j++){		
      depths[i][j] = 1/(1 + depths[i][j]);
    }
  }
  deleteM(objectsProjections);
  return 0;
}
// #ifdef __cplusplus
// extern "C" {
// #endif


void Sum(double *a, double *b, double *res){
  res[0] = a[0] + b[0];
}

void setSeed(int random_seed){
  if (random_seed != 0) {
    setseed(random_seed);
    rEngine.seed(random_seed);
  }
  else {
    setseed(time(NULL));
    rEngine.seed(time(NULL));
  }
}

void ProjectionDepth(double *points, double *objects, int *numObjects,
                     int *dimension, int *cardinalities, int *numClasses,
                     double *directions, double *projections, int *k,
                     int *newDirs, int *seed, double *depths, int *dirStyle){
  setSeed(*seed);
  TVariables cars(numClasses[0]);
  int numPoints = 0;
  for (int i = 0; i < numClasses[0]; i++){
    numPoints += cardinalities[i];
    cars[i] = cardinalities[i];
  }
  TDMatrix x = asMatrix(points, numPoints, *dimension);
  TDMatrix z = asMatrix(objects, *numObjects, *dimension);
  
  TDMatrix dirs = asMatrix(directions, k[0], *dimension);
  TDMatrix prjs = asMatrix(projections, k[0], numPoints);
  TDMatrix _depths = asMatrix(depths, *numObjects, *numClasses);
  GetDepthsPrj(x, numPoints, *dimension, z, *numObjects, cars, 
               *k, *newDirs, _depths, dirs, prjs, *dirStyle);
  /*	for (int i = 0; i < numObjects[0]; i++){
   for (int j = 0; j < numClasses[0]; j++){
   depths[i * numClasses[0] + j] = _depths[i][j];
   }
  }
   if (newDirs[0]){
   for (int i = 0; i < k[0]*dimension[0]; i++){
   directions[i] = dirs[i/dimension[0]][i%dimension[0]];
   }
   for (int i = 0; i < k[0]*numPoints; i++){
   projections[i] = prjs[i/numPoints][i%numPoints];
   }
   }
   */
  
  delete[] x;
  delete[] z;
  delete[] dirs;
  delete[] prjs;
  delete[] _depths;
}

// #ifdef __cplusplus
// }
// #endif

// [[Rcpp::export]]
Rcpp::NumericVector proj_depth(const arma::mat& X, const arma::mat& data, int style=3, int multiplier = 10){
  int m = X.n_rows;
  int p = X.n_cols;
  int n = data.n_rows;
  int numClasses = 1;
  int k = max(1000,multiplier*p);
  int newdirs = 1;
  int seed = 0;
  double* points = new double[n*p];
  for (int i = 0; i < n; i++){
    for (int j = 0; j < p; j++)
      points[i*p+j] = data(i,j);
  }
  double* objects = new double[m*p];
  for (int i = 0; i < m; i++){
    for (int j = 0; j < p; j++)
      objects[i*p+j] = X(i,j);
  }
  double* dirs = new double[k*p];
  double* projs = new double[k*n];
  double* depths_ptr = new double[m];
  ProjectionDepth(points, objects, &m, &p, &n, &numClasses, dirs, projs, &k, &newdirs, &seed, depths_ptr, &style);
  arma::colvec depths = arma::zeros<arma::colvec>(m);
  for (int i = 0; i < m; i++){
    depths(i) = depths_ptr[i];
  }
  delete[] points;
  delete[] objects;
  delete[] dirs;
  delete[] projs;
  delete[] depths_ptr;
  return wrap(depths);
}