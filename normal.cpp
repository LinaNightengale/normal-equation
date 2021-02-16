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

vector <double> cholesky_decomposition(vector <vector <double>> & matrix_a, vector <double> & phi, int N)
{
    vector <vector <double>> L(N, vector <double> (N, 0)), T(N, vector <double> (N, 0));
    vector <double> x(N, 0), v(N, 0);
    double sum = 0, buffer = 0;

    L[0][0] = sqrt(matrix_a[0][0]);
    for (int i = 1; i < N; i++)
        L[i][0] = matrix_a[i][0] / L[0][0];
 
    for (int i = 1; i < N; ++i)
    {
        sum = 0;
        for (int j = 0; j < i; ++j)
            sum = sum + pow(L[i][j],2);
        buffer = matrix_a[i][i] - sum;
        L[i][i]=sqrt(buffer);
 
        for (int j = i + 1; j < N; ++j)
        {
            sum = 0;
            for (int k = 0; k < i; ++k)
                sum = sum + L[j][k] * L[i][k];
            buffer = matrix_a[j][i] - sum;
            L[j][i] = buffer / L[i][i];
        }
    }
 
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            T[i][j]=L[j][i];
    
    v[0] = phi[0] / L[0][0];
    for (int i = 1; i < N; ++i)
    {
        sum = 0;
        for (int j = 0; j < i; ++j)
            sum = sum + L[i][j] * v[j];
        buffer = phi[i] - sum;
        v[i] = buffer / L[i][i];
    }

    x[N - 1] = v[N - 1] / T[N - 1][N - 1];
    for (int i = N - 2; i >= 0; --i)
    {
        sum = 0;
        for (int j = i + 1; j < N; j++)
            sum = sum + T[i][j] * x[j];
        buffer = v[i] - sum;
        x[i] = buffer / T[i][i];
    }

    return x;
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
    
    x = cholesky_decomposition(ATA, ATphi, N);
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
