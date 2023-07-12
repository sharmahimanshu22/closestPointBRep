#include<iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Dense>


struct Surface {

  Surface() {
  }

  void eval_surf(double uv[2], double coord[3], double deriv1[2][3], double deriv2[3][3], double unit_norm[3]) {
    double R = 1.0;

    double u = uv[0];
    double v = uv[1];

    coord[0] = R*sin(u)*cos(v);
    coord[1] = R*sin(u)*sin(v);
    coord[2] = R*cos(u);

    deriv1[0][0] = R*cos(u)*cos(v);
    deriv1[0][1] = R*cos(u)*sin(v);
    deriv1[0][2] = -R*sin(u);

    deriv1[1][0] = -R*sin(u)*sin(v);
    deriv1[1][1] = R*sin(u)*cos(v);
    deriv1[1][2] = 0.0;
  
    deriv2[0][0] = -R*sin(u)*cos(v);
    deriv2[0][1] = -R*sin(u)*sin(v);
    deriv2[0][2] = -R*cos(u);

    deriv2[1][0] = -R*cos(u)*sin(v);
    deriv2[1][1] = R*cos(u)*cos(v);
    deriv2[1][2] = 0.0;

    deriv2[2][0] = -R*sin(u)*cos(v);
    deriv2[2][1] = -R*sin(u)*sin(v);
    deriv2[2][2] = 0.0;

    unit_norm[0] = -sin(u)*cos(v);
    unit_norm[1] = -sin(u)*sin(v);
    unit_norm[2] = -cos(u);

    return;    
  }

  double closestPointOnCircle(double point[3], double surfacePoint[3]) {
    double len = sqrt(point[0]*point[0] + point[1]*point[1] + point[2]*point[2]);
    surfacePoint[0] = point[0]/len;
    surfacePoint[1] = point[1]/len;
    surfacePoint[2] = point[2]/len;
  }
};
  

double norm2(double p[2]) {
  return sqrt(p[0]*p[0] + p[1]*p[1]);
}

double isNormZero(double p[2]) {
  if (norm2(p) < 0.000001) return true;
  return false;
}

void computeInverse(double mat[2][2], double invMat[2][2]) {
  Eigen::Matrix<float, 2, 2, Eigen::RowMajor> eigenMat;
  eigenMat.setZero();
  eigenMat(0,0) = mat[0][0];
  eigenMat(0,1) = mat[0][1];
  eigenMat(1,0) = mat[1][0];
  eigenMat(1,1) = mat[1][1];

  Eigen::Matrix<float, 2, 2, Eigen::RowMajor> invMatEigen = eigenMat.inverse();

  invMat[0][0] = invMatEigen.coeff(0,0);
  invMat[0][1] = invMatEigen.coeff(0,1);
  invMat[1][0] = invMatEigen.coeff(1,0);
  invMat[1][1] = invMatEigen.coeff(1,1);

}

void computeDeltaU(double DPrime[2], double D2PrimeInverse[2][2], double deltaU[2]) {
  deltaU[0] = D2PrimeInverse[0][0] * DPrime[0] + D2PrimeInverse[0][1]*DPrime[1];
  deltaU[1] = D2PrimeInverse[1][0] * DPrime[0] + D2PrimeInverse[1][1]*DPrime[1];
}

void computeDPrime(double point[3], double currentPointOnSurface[3], double deriv1[2][3], double DPrime[2]) {
  DPrime[0] = 0.0;
  DPrime[1] = 0.0;
  
  double fsubtractp[3] {currentPointOnSurface[0] - point[0], currentPointOnSurface[1] - point[1], currentPointOnSurface[2] - point[2]};
  for(int i = 0; i < 3; i++) {
    DPrime[0] += 2 * fsubtractp[i] * deriv1[0][i];
    DPrime[1] += 2 * fsubtractp[i] * deriv1[1][i];
  }
}

void computeD2Prime(double point[3], double currentPointOnSurface[3], double deriv1[2][3], double deriv2[3][3], double D2Prime[2][2]) {

  D2Prime[0][0] = 0.0;
  D2Prime[0][1] = 0.0;
  D2Prime[1][0] = 0.0;
  D2Prime[1][1] = 0.0;
  
  
  for(int i = 0; i < 3; i++) {
    D2Prime[0][0] += 2*deriv1[0][i]*deriv1[0][i];
    D2Prime[0][1] += 2*deriv1[0][i]*deriv1[1][i];
    D2Prime[1][0] += 2*deriv1[1][i]*deriv1[0][i];
    D2Prime[1][1] += 2*deriv1[1][i]*deriv1[1][i];
  }  

  double fsubtractp[3] {currentPointOnSurface[0] - point[0], currentPointOnSurface[1] - point[1], currentPointOnSurface[2] - point[2]};
  double V[3] {0.0, 0.0, 0.0};
  for(int i = 0; i < 3; i++) {
    V[0] += deriv2[0][i]*fsubtractp[i];
    V[1] += deriv2[1][i]*fsubtractp[i];
    V[2] += deriv2[2][i]*fsubtractp[i];
  }


  D2Prime[0][0] += 2*V[0];
  D2Prime[0][1] += 2*V[1];
  D2Prime[1][0] += 2*V[1];
  D2Prime[1][1] += 2*V[2];
  
}

void updateParameters(double point[3], double currentPoint[3], double deriv1[2][3], double deriv2[3][3], double currentParameters[2]) {

  double DPrime[2] = {0.0};
  computeDPrime(point, currentPoint, deriv1, DPrime);
  double D2Prime[2][2] = {0};
  computeD2Prime(point, currentPoint, deriv1, deriv2, D2Prime);
  double D2PrimeInverse[2][2];
  computeInverse(D2Prime, D2PrimeInverse);
  double deltaU[2] {0.0, 0.0};
  computeDeltaU(DPrime, D2PrimeInverse, deltaU);
  
  currentParameters[0] = currentParameters[0] - deltaU[0];
  currentParameters[1] = currentParameters[1] - deltaU[1];

  /*
  if(currentParameters[0] < 0.0) {
    currentParameters[0] = 0.0;
  }
  if(currentParameters[0] > M_PI) {
    currentParameters[0] = M_PI;
  }
  if(currentParameters[1] < 0.0) {
    currentParameters[1] = 0.0;
  }
  if(currentParameters[1] > M_PI) {
    currentParameters[1] = M_PI;
  }
  */
  
  return;
}

void findProjection(Surface* s, double point[3], double guessParameters[2], double finalParameters[2], double projectedPoint[3] ) {

  // check if the initial guess lies within domain.
  // if the guess is null, come up with an initial guess by yourself

  int MAX_ITERATIONS = 100;
  
  double currentParameters[2] {guessParameters[0], guessParameters[1]};
  double currentPoint[3];
  double deriv1[2][3];
  double deriv2[3][3];
  double unit_norm[3];

  int iter = 0;
  while(iter < MAX_ITERATIONS) {
  
    iter = iter+1;
    s->eval_surf(currentParameters, currentPoint, deriv1, deriv2, unit_norm);

    // We compute DPrime even inside updateParameters method. Should check if we can refactor.
    double DPrime[2];
    computeDPrime(point, currentPoint, deriv1, DPrime);

    if(isNormZero(DPrime)) {
      finalParameters[0] = currentParameters[0];
      finalParameters[1] = currentParameters[1];
      break;
    } 
    updateParameters(point,currentPoint, deriv1, deriv2, currentParameters);
  }

  double formulaPoint[3];
  s->closestPointOnCircle(point, formulaPoint);

  std::cout << "Nearest Point computed via Newton's method: " << currentPoint[0] << " " << currentPoint[1] << " " << currentPoint[2] << "\n";
  std::cout << "Nearest Point computed from formula:        " << formulaPoint[0] << " " << formulaPoint[1] << " " << formulaPoint[2] << "\n";
  std::cout << "Final Parameters: " << currentParameters[0] << " " << currentParameters[1] << "\n"; 

  std::cout << "Number of iterations: " << iter << "\n"; 
}


int main(int argc, char **argv) {

  if (argc < 6) {
    return -1;
  }

  double point[3];
  Surface* s = new Surface();
  double finalParameters[2];
  double projectedPoint[3];  
  double guessParameters[2];

  point[0] = atof(argv[1]);
  point[1] = atof(argv[2]);
  point[2] = atof(argv[3]);
  guessParameters[0] = atof(argv[4]);
  guessParameters[1] = atof(argv[5]);

  
  findProjection(s, point, guessParameters, finalParameters, projectedPoint);
  
  return EXIT_SUCCESS;
}
