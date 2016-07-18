#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <cstdlib>
#include<conio.h>
#include <fstream>
using namespace std;

void print(vector< vector<double> > A) {
	int n = A.size();
	for (int i = 0; i<n; i++) {
		
		for (int j = 0; j<n + 1; j++) {
			cout << A[i][j] << "\t";
			if (j == n - 1) {
				cout << "| ";
			}
		}
		cout << "\n\n\n";
	}
	cout << endl;
}

double gauss(vector< vector<double> > A)
{
	int n = A.size();

	for (int i = 0; i<n; i++) {
		// Szukam maksymalnego elementu w kolumnie
		double maxEl = abs(A[i][i]),temp;
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

	return x[0];
}

int zmienpozycje(int poz, int traps[],int n)
{
	int oczka, pulapka;
	oczka = rand() % 6 + 1;
	poz += oczka;
	if (poz < n)
	{
		pulapka = traps[poz];
		poz += pulapka;
	}
	else
		poz = n; //wygrana
		
	return poz;
}
double MonteCarlo(double ile, int traps[], int n)
{

	int g1=0, g2=0; // startowe pozycje graczy
	double wg1=0, wg2 = 0;
	bool koniec = false;
	for (int i = 0; i < ile; i++)
	{			
		koniec = false;
		g1 = 0;
		g2 = 0;
		while (!koniec)
		{

			g1 = zmienpozycje(g1, traps, n);
			if (g1 < n)
			{
				g2 = zmienpozycje(g2, traps, n);
				if (g2 >= n)
				{
					wg2++;
					koniec = true;
				}
			}
			else
			{
				wg1++;
				koniec = true;
			}
	
		}
	}
	return wg1 / ile;
}
double MetodaGaussaSeidla(vector< vector<double> > A, int m, int ilosc)
{
	vector<double>x(m); //prawdopodobienstwa
	int il = 0,licz=0;
	double mian;
	while (il != ilosc)
	{
		
		for (int i = 0; i < m; i++)
		{
			mian = 0;
			mian += A[i][m]; // ostatnia wartosc z macierzy w danym wierszu
			for (int j = 0; j < m; j++)
			{
				if (i!=j)
					mian += A[i][j] * (-1)*x[j];
			}
			x[i] = mian / A[i][i];
		}
		
		if (licz == 9)
		{
			printf("%d %.16lf\n",il+1, x[0]);
			licz = 0;
		}
		else
		{
			licz++;
		}
		il++;
	}

	return x[0];
}

int main() {
	srand(time(NULL));
	int n = 6;
	int fields=2*n*n,ile=0;
	int cube = 6;
	vector<double> line(fields + 1, 0);
	vector< vector<double> > A(fields, line);
	int traps[] = {
		/*0, 0, -1, 0, -4, -2, -5, -6, 0, 0, -1,
		0, 0, -2, 0, -1, 0, 0, -4, 0, -9, 0, -3,
		0, -1, 0, -26, 0, 0*/
		0,0,0,0,-2,-2
	};

	//tworze macierz do gry
	int y;
	for (int a = 0; a < n; a++)
	{
		for (int b = 0; b < n; b++)
		{
			for (int k = 0; k < 2; k++)
			{
				int x = b*n + a + k*n*n;
				A[x][x] -= cube; //przesuniecie danego x na druga strone rownania z przeciwnym znakiem

					for (int l = 1; l <= cube; l++)
					{
						if (k == 0)
						{
							if (a + l < n)
							{
								y = b*n + l + a + traps[a + l] + (k == 0 ? 1 : 0)*n*n;
								A[x][y]++;
							}
							else{
								y = fields;
								A[x][y]--;
							}
						}
						else{
							if (b + l < n)
							{
								y = a*n + l + b + traps[b + l] + (k == 0 ? 1 : 0)*n*n;
								A[x][y]++;
							}
						}
					}
			}
		}
	}
	
	//dzielenie, na dane prawdopodobienstwo
	for (int a = 0; a < fields; a++)
	{
		for (int b = 0; b <= fields; b++)
		{
			A[a][b] /= cube;
		}
	}
	
	printf("Gauss %.16lf", gauss(A));
	printf("Monte Carlo%d %.16lf\n", MonteCarlo(1000000, traps, n));
	printf("Monte Carlo%d %.16lf\n", MetodaGaussaSeidla(A, fields, 400));
	_getch();
}