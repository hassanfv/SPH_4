#include <iostream>
#include <chrono>

using namespace std;

int main()
{

auto start = std::chrono::high_resolution_clock::now();

auto end = std::chrono::high_resolution_clock::now();
chrono::duration<double> elapsed = end - start;
cout << "Elapsed time: " << elapsed.count() << " seconds" << endl;


}
