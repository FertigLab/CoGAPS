//compile: 
//c++ static_cast_standalone_test.cpp -std=gnu++17 -O0
//clang++ -arch x86_64 -std=gnu++17 -O0 static_cast_standalone_test.cpp
#include <cstdint>
#include <iostream>

void priint(uint64_t a){
    std::cout<<a<<std::endl;
}

int main(){
    uint64_t uvalue64=std::numeric_limits<uint64_t>::max();
    //maximal is 18446744073709551615
    //effect holds for value 18446744073709551600
    //std::cout<<"Max uint64_t:\n"<<std::numeric_limits<uint64_t>::max()<<std::endl;
    std::cout<<uvalue64<<std::endl;
    double dvalue=static_cast<double>(uvalue64);
    priint(dvalue);
    std::cout<<dvalue<<std::endl;
    std::cout<<static_cast<uint64_t>(dvalue)<<std::endl;
    uint32_t uvalue32=std::numeric_limits<uint32_t>::max();
    //maximal is 18446744073709551615
    //effect holds for value 18446744073709551600
    //std::cout<<"Max uint64_t:\n"<<std::numeric_limits<uint64_t>::max()<<std::endl;
    std::cout<<uvalue32<<std::endl;
    dvalue=static_cast<double>(uvalue32);
    priint(dvalue);
    return 0;
}
