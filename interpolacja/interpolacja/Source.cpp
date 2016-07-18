//Sebastian Sarnecki 155189
#include <iostream>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>
#define N 20
using namespace std;
double u[N + 1], y[N + 1], h[N + 1], delta[N + 1];

// x-vector M, tab - wartosci x, s21 - wartosci funkcji s21, j - odpowiada za dany przedzial, X - konkretny X
double funkcja(vector<double> x, double tab[], double s21[], int j, double X)
{
	double a, b, c, d;
	double xt = X - tab[j];
	a = s21[j];
	b = ((s21[j + 1] - s21[j]) / h[j + 1]) - ((2 * x[j] + x[j + 1]) / 6) * h[j + 1];
	c = x[j] / 2;
	d = (x[j + 1] - x[j]) / (6 * h[j + 1]);

	return abs(a + b*(xt)+c*pow(xt, 2) + d*pow(xt, 3));

}
void zapiszWyniki(vector<double> x, double tab[], double s21[])
{
	int dzielnik = 20; //liczba przedzia³ów w jednym z 19 przedzialow czestotliwosci
	double wynik;
	ofstream o;

	o.open("przedzial49_2.txt");
	for (int j = 0; j < x.size() - 1; j++)
	{
		double roznica = tab[j + 1] - tab[j];
		roznica /= (double)dzielnik;
		wynik = tab[j];
		if (j != x.size()-2)
		{
			for (int i = 0; i < dzielnik; i++)
			{
				o << wynik << " " << 20 * log10(funkcja(x, tab, s21, j, wynik)) << endl;;
				wynik += roznica;
			}
		}
		else
			o << wynik << " " << 20 * log10(funkcja(x, tab, s21, j, wynik)) << endl;;
	}
	o.close();
}
void gauss(vector< vector<double> > A, double tab[], double s21[])
{
	int n = A.size();
	

	for (int i = 0; i<n; i++) {
		// Szukam maksymalnego elementu w kolumnie
		double maxEl = abs(A[i][i]), temp;
		int maxRow = i;
		for (int w = i + 1; w<n; w++) {
			temp = abs(A[w][i]);
			if (temp > maxEl) {
				maxEl = temp;
				maxRow = w;
			}
		}

		//zamiana obecnego wiersza z maksymalnym
		for (int k = i; k<n + 1; k++)
			swap(A[maxRow][k], A[i][k]);


		//zerowanie wiersza w danej kolumnie(macierz schodkowa)
		for (int k = i + 1; k<n; k++) {
			double c = A[k][i] / A[i][i];
			for (int j = i; j<n + 1; j++) {
				if (i == j)
					A[k][j] = 0;
				else
					A[k][j] -= c * A[i][j];
			}
		}
	}

	// Ax=b, rozwiazywanie rownan, z postaci macierzy schodkowej
	vector<double> x(n);
	for (int i = n - 1; i >= 0; i--)//petla przechodzaca po vectorze wynikow
	{
		x[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--)
			A[k][n] -= A[k][i] * x[i];
	}

	for (int i = 0; i < x.size();i++)
		printf("%d - %.16lf\n", i, x[i]);

	zapiszWyniki(x, tab, s21);
}

int main() {
	vector<double> line(N+2, 0);
	vector< vector<double> > A(N+1, line);

	double x[N] = { 2.160913, 2.184642, 2.208656, 2.232956, 2.257543, 2.282417, 2.307579, 2.333029,
		2.358767, 2.384794, 2.411110, 2.437714, 2.464608, 2.491789, 2.519259, 2.547017, 2.575062,
		2.603393, 2.632010, 2.660913 };

	double s21[N] = { 0.0385736053, 0.0418602337, 0.0406238069, 0.0222182757, 0.0681165327, 0.5274988276,
		0.9992090072, 0.9877714749, 0.9997958915, 0.9872738013, 0.9872738013, 0.9997958915, 0.9877714749,
		0.9992090072, 0.5274988276, 0.0681165327, 0.0222182757, 0.0406238069, 0.0418602337, 0.0385736053 };

	//warunki brzegowe
	u[0] = 0;
	u[N] = 0;
	y[0] = 0;
	y[N] = 0;

	//wypelniam tablice wartosci h
	for (int j = 0; j < N; j++)
	{
		h[j + 1] = x[j + 1] - x[j];
	}

	//wypelniam tablice wartosci ui oraz yj
	for (int j = 1; j < N; j++)
	{
		u[j] = h[j] / (h[j] + h[j + 1]);
		y[j] = h[j + 1] / (h[j] + h[j + 1]);
	}

	//uzupelnienie delty w macierzy
	for (int j = 1; j < N; j++)
	{
		A[j][N + 1] = (6 / (h[j] + h[j + 1])) * ((s21[j + 1] - s21[j]) / h[j + 1] - (s21[j] - s21[j - 1]) / h[j]);
	}

	//uzupelnienie warunkow brzegowych
	A[0][0] = 2;
	A[0][1] = y[0];
	A[N][N-1] = u[N];
	A[N][N] = 2;
	A[0][N+1] = 0;
	A[N][N + 1] = 0;

	//uzupelnjenje macierzy
	for (int j = 1; j < N; j++)
	{
		//j-1
		A[j][j-1] = u[j];
		//j
		A[j][j] = 2;
		//j+1
		A[j][j+1] = y[j];
	}

	gauss(A, x, s21);
}