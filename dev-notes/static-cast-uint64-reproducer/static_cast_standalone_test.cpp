//compile: 
//c++ static_cast_standalone_test.cpp -std=gnu++17 -O0
#include <cstdint>
#include <iostream>
#include <limits>

int main(){
    uint64_t uvalue64_m=std::numeric_limits<uint64_t>::max();
    uint64_t uvalue64_s=uvalue64_m-10;
    //maximal is 18446744073709551615
    //effect holds for value 18446744073709551600
    //std::cout<<"Max uint64_t:\n"<<std::numeric_limits<uint64_t>::max()<<std::endl;
    std::cout<<"uint64_t value 1 (max) : "<<uvalue64_m<<std::endl;
    std::cout<<"uint64_t value 2 (max-10) : "<<uvalue64_s<<std::endl;
    double dvalue_m=static_cast<double>(uvalue64_m);
    double dvalue_s=static_cast<double>(uvalue64_s);
    std::cout<<"value 1 as double : "<<dvalue_m<<std::endl;
    std::cout<<"value 2 as double : "<<dvalue_s<<std::endl;
    std::cout<<"value 1 as uint64_t casted back from double : "<<static_cast<uint64_t>(dvalue_m)<<std::endl;
    std::cout<<"value 2 as uint64_t casted back from double : "<<static_cast<uint64_t>(dvalue_s)<<std::endl;
    return 0;
}
