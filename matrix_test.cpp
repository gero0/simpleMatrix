#include "Matrix.h"

int main(){
    mx::Matrix<double> m1 { {2,2}, {2,1,3,7}};
    auto m2 = m1.transpose();

    for (auto v : m2.getVector()){
        std::cout << v << " ";
    }
   
    return 0;
}