#include <iostream>
#include <vector>
#include "normal.h"

using namespace std;

vector <vector <double>> transpose(vector <vector <double>> & matrix_a, int m, int N)
{
    vector <vector <double>> matrix_a_t(N, vector <double> (m, 0));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < m; ++j)
            matrix_a_t[i][j] = matrix_a[j][i];

    return matrix_a_t;
}

vector <vector <double>> mult(vector <vector <double>> matrix_left, vector< vector<double>> matrix_right)
{
  int l = matrix_left.size(), m = matrix_left[0].size(), n = matrix_right[0].size();

  vector <vector <double>> matrix(l, vector <double> (n, 0));

  for (int i = 0; i < l; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < m; ++k)
        matrix[i][j] += matrix_left[i][k] * matrix_right[k][j];

  return matrix;
}

vector <double> reverse_Gaussian(vector <vector <double>> & matrix_a, vector <double> & phi, int N)
{
    double buffer_1, buffer_2;
    vector <double> matrix_x(N);
    
    for (int i = N - 1; i >= 0; --i)
    {
        buffer_1 = 0.;
        for (int j = i + 1; j < N; ++j)
        {
            buffer_2 = matrix_a[i][j] * matrix_x[j];
            buffer_1 += buffer_2;
        }
        matrix_x[i] = (phi[i] - buffer_1) / matrix_a[i][i];
    }
    return matrix_x;
}

vector <double> normal_equation(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N)
{
    vector <double> x;
    vector <vector <double>> AT = transpose(matrix_a, m, N);
    vector <vector <double>> ATA = mult(AT, matrix_a);
    
    vector <double> ATphi (N);

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < m; ++j)
        {
            ATphi[i] += AT[i][j] * phi[j];
        }
    
    x = reverse_Gaussian(ATA, ATphi, N);
    return x;
}

void scan_m(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N)
{
    cout << "Enter the matrix A: " << endl;
    for (int i = 0; i < m; ++i)
        for (int j = 0; j < N; ++j)
        {
            cin >> matrix_a[i][j];
        }
  
    cout << "Enter the matrix b: " << endl;
    for (int i = 0; i < m; ++i)
    {
        cin >> phi[i];
    }   
}

void print_m(vector <double> & x, int N) 
{
    cout << "Result: " << endl;
    for (int i = 0; i < N; ++i)
        cout << x[i] << endl;
}
