#include <vector>

using namespace std;

vector <vector <double>> transpose(vector <vector <double>> & matrix_a, int m, int N);

vector <vector <double>> mult(vector <vector <double>> matrix_left, vector< vector<double>> matrix_right);

vector <double> reverse_Gaussian(vector <vector <double>> & matrix_a, vector <double> & phi, int N);

vector <double> normal_equation(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N);

void scan_m(vector <vector <double>> & matrix_a, vector <double> & phi, int m, int N);

void print_m(vector <double> & x, int N);