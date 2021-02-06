#include <iostream>
#include <vector>
#include "normal.h"

using namespace std;

int main()
{
    int m = 0, N = 0;

    cout << "Enter the number of rows: ";
    cin >> m;
    cout << "Enter the number of columns: ";
    cin >> N;

    vector <vector <double>> A(m, vector <double> (N, 0));
    vector <double> b(m, 0), x;

    scan_m(A, b, m, N);

    x = normal_equation(A, b, N);
    
    print_m(x, N);

    return 0;
}