#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
using namespace std; 

int main() {
    std::ifstream input_file("example.txt");
    string str2, str3;
    if (!input_file.is_open()) {
        std::cerr << "Failed to open file!" << std::endl;
        return 1;
    }

    string str;

    input_file  >> setw(6) >> str;
    input_file  >> setw(6) >> str2;
    cout<<(str2.length())<<endl;
    input_file  >> setw(6) >> str;
    cout<<(str.length())<<endl;

    std::cout << str << std::endl;
    std::cout << str2 << std::endl;
    std::cout << str3 << std::endl;



    input_file.close();
    return 0;
}


