// output stream
#include <iostream>
#include <fstream>
#include <cantera/numerics/eigen_dense.h>
#include <cantera/numerics/eigen_sparse.h>
#include <string>

//We use the Thomas algorithm for solving the block tridiagonal matrix system
Eigen::MatrixXd ThomasSolver(Eigen::MatrixXd A, Eigen::MatrixXd B, int nsp, int nSize){
//std::cout << "Starting Thomas solver" << std::endl;
//Initialize output vector
Eigen::MatrixXd Xvec = Eigen::VectorXd::Zero(nsp*nSize);
std::vector<Eigen::MatrixXd> C_prime(nSize);
std::vector<Eigen::MatrixXd> D_prime(nSize);
//std::cout << "Initialization done!" << std::endl;

for(int i = 0; i < nSize; i++){


//Forward elimination
  if(i == 0){
  Eigen::MatrixXd B1 = A.block(0, 0, nsp, nsp);
  //std::cout << "B1 done" << std::endl;
  Eigen::MatrixXd C1 = A.block(0, nsp, nsp, nsp);
  //std::cout << "C1 done" << std::endl;
  Eigen::MatrixXd D1 = B.block(0, 0, nsp, 1);
  //std::cout << "D1 done"  << std::endl;
  
  C_prime[0] = B1.inverse()*C1;
  D_prime[0] = B1.inverse()*D1;
  //std::cout << "Forward elimination at point 0 complete" << std::endl;
  }
  
  if(i > 0 && i < nSize - 1){
  Eigen::MatrixXd Bi = A.block(i*nsp, i*nsp, nsp, nsp);
  //std::cout << "Bi done" <<  std::endl;
  Eigen::MatrixXd Ai = A.block((i-1)*nsp, i*nsp, nsp, nsp); 
  //std::cout << "Ai done" << std::endl;
  Eigen::MatrixXd Ci = A.block(i*nsp, (i+1)*nsp, nsp, nsp); 
  //std::cout << "Ci done" << std::endl;
  Eigen::MatrixXd Di = B.block(i*nsp, 0, nsp, 1);
  //std::cout << "Di done" << std::endl;
  
  C_prime[i] = (Bi - Ai*C_prime[i-1]).inverse()*Ci;
  //std::cout << "C prime done"  << std::endl;
  D_prime[i] = (Bi - Ai*C_prime[i-1]).inverse()*(Di - Ai*D_prime[i-1]);
  //std::cout << "D prime done"  << std::endl;

  //std::cout << "Forward elimination at point " << i << " complete" << std::endl;  
  
  }
  
  if(i == nSize - 1){
  Eigen::MatrixXd Bi = A.block(i*nsp, i*nsp, nsp, nsp);
  //std::cout << "Bi done" << std::endl;
  Eigen::MatrixXd Ai = A.block((i-1)*nsp, i*nsp, nsp, nsp); 
  //std::cout << "Ai done" << std::endl;
  Eigen::MatrixXd Di = B.block(i*nsp, 0, nsp, 1);
  //std::cout << "Di done" << std::endl;
  
  D_prime[i] = (Bi - Ai*C_prime[i-1]).inverse()*(Di - Ai*D_prime[i-1]);
  //std::cout << "D prime done"  << std::endl;

  //std::cout << "Forward elimination at point " << i << " complete" << std::endl;  
  
  }
  
  
}  
 //std::cout << "Forward elimination complete!" << std::endl;
//Back substitution
 //std::cout << "Starting backward substitution" << std::endl;

for(int i = nSize-1; i >= 0; i--){
   //std::cout << "D_prime: " << D_prime[i] << std::endl;
   if(i == nSize-1){
   Xvec.block(i*nsp, 0, nsp, 1) = D_prime[i] ;
   //std::cout << D_prime[i] << std::endl;
   //std::cout << "Xvec complete" << std::endl;
   }

   else{
   
   Eigen::MatrixXd X_next = Xvec.block((i+1)*nsp, 0, nsp, 1);
   Xvec.block(i*nsp, 0, nsp, 1) = D_prime[i] - C_prime[i]*X_next;
   //std::cout << Xvec.block(i*nsp, 0, nsp, 1) << std::endl;
   
   }

}



return Xvec;
} 
