#include <iostream>
#include <cmath>
#include <algorithm>

using namespace std;


int main()
{

    int x = 2;
    int y = 3;

    int c = max(2, 3);

    float h = -1.1;

    c = (h > 0. ? max(x, y) : min(x, y));

    cout << "c = " << c << endl;

    return 0;
}