#include<NTL/ZZ.h> //pre datovy typ ZZ
#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>
#include <iostream>
#include <tuple>
#include <iostream>
#include <NTL/mat_GF2.h>
#include <NTL/vec_GF2.h>
#include <NTL/vec_vec_GF2.h>
#include "riesenie_sustavy.h"
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>
#include "UOV.h"

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

using namespace std;

NTL_CLIENT




Vec<GF2E> riesenia(Vec<Mat<GF2E>> polynomy_Q, Vec<Vec<GF2E>> polynomy_L, Vec<GF2E> polynomy_A, int n, int m) {
	GF2E res;
	Vec<GF2E> solutions;
	Vec<GF2E> h;
	h = random_vec_GF2E(n * 2);

	for (long k = 0; k < m; k++)
	{

		res = h * polynomy_Q[k] * h + polynomy_L[k] * h + polynomy_A[k];
		solutions.append(res);


	}


	//cout << "riesenie:" << solutions << endl;
	//cout << "hodnoty premennych" << h << endl;

	return solutions;
}




void generate_random_polynomials(int m, int n, Vec<Mat<GF2E>>& polynomy_Q,
	Vec<Vec<GF2E>>& polynomy_L,
	Vec<GF2E>& polynomy_A) {
	polynomy_Q.SetLength(m);
	polynomy_L.SetLength(m);
	polynomy_A.SetLength(m);

	for (int i = 0; i < m; i++) {
		// Generate random quadratic part of the polynomial
		Mat<GF2E> Q;
		Q.SetDims(m + n, m + n);
		for (int j = 0; j < m; j++) {
			for (int k = m; k < m + n; k++) {
				Q[j][k] = random_GF2E();
			}
		}
		for (int j = m; j < m + n; j++) {
			for (int k = j; k < m + n; k++) {
				Q[j][k] = random_GF2E();
			}
		}
		polynomy_Q[i] = Q;

		// Generate random linear part of the polynomial
		Vec<GF2E> L;
		L.SetLength(m + n);
		for (int j = 0; j < m + n; j++) {
			L[j] = random_GF2E();
		}
		polynomy_L[i] = L;

		// Generate random absolute part of the polynomial
		polynomy_A[i] = random_GF2E();
	}
}

void KeyGen(publicKey& pk, privateKey& sk, long m_poly, long n_variables)
{


	Vec<Mat<GF2E>> polynomy_Q; //kvadraticke casti polynomov
	Vec<Vec<GF2E>> polynomy_L; //linearne casti polynomov
	Vec<GF2E> polynomy_A; //absolutne casti polynomov

	long m = m_poly; //pocet polynomov a zaroven pocet olej
	long n = n_variables; //pocet neurcitych ocot
	generate_random_polynomials(m, n, polynomy_Q, polynomy_L, polynomy_A);

	
	/*vypis jednotlivych polynomov
	for (long k = 0; k < m; k++)
	{
		cout << "Polynom cislo: " << k + 1 << endl;
		//kvadraticka cast, linearna cast, absolutna cast
		cout << polynomy_Q[k] << endl;
		cout << polynomy_L[k] << endl;
		cout << polynomy_A[k] << endl;
		cout << "**********" << endl;
	}
	*/
	//TRANSFORMACIA T
	Mat<GF2E> A_T; //musi byt invertovatelne nad GF2
	Vec<GF2E> b_T; //nahodny vektor hodnot GF2
	Mat<GF2E> temp_matrix; //pomocna docasna matica
	while (1)
	{
		random(A_T, n + m, n + m); //vytvori nahodnu maticu
		temp_matrix = A_T; //ulozime jej kopiu do docasnej premennej
		//otestujeme, ci je vygenerovana matica invertovatelna
		
		if (gauss(temp_matrix) == n + m)
			break;
	}
	random(b_T, n + m); //nahodny vektor
	Vec<Mat<GF2E>> polynomy_Q_T; //kvadraticke casti polynomov po aplik T
	Vec<Vec<GF2E>> polynomy_L_T; //linearne casti polynomov po aplik T
	Vec<GF2E> polynomy_A_T; //absolutne casti polynomov po aplik T
	
	for (long k = 0; k < m; k++)
	{
		polynomy_Q_T.append(A_T * polynomy_Q[k] * transpose(A_T));
		polynomy_L_T.append(b_T * polynomy_Q[k] * transpose(A_T) + b_T * transpose(polynomy_Q[k]) * transpose(A_T) + polynomy_L[k] * transpose(A_T));
		polynomy_A_T.append(b_T * polynomy_Q[k] * b_T + polynomy_L[k] * b_T + polynomy_A[k]);
	}
	/*vypis jednotlivych polynomov
	for (long k = 0; k < m; k++)
	{
		cout << "Polynom cislo: po trasformacii T" << k + 1 << endl;
		//kvadraticka cast, linearna cast, absolutna cast
		cout << polynomy_Q_T[k] << endl;
		cout << polynomy_L_T[k] << endl;
		cout << polynomy_A_T[k] << endl;
		cout << "**********" << endl;
	}
	cout << "bt" << b_T << endl;
	cout << "At" << A_T << endl;
	*/
	//TRANSFORMACIA S
	Mat<GF2E> A_S; //musi byt invertovatelne nad GF2
	Vec<GF2E> b_S; //nahodny vektor hodnot GF2

	while (1)
	{
		random(A_S, m, m); //vytvori nahodnu maticu s velkostou m podla poctu polynomov
		temp_matrix = A_S; //ulozime jej kopiu do docasnej premennej
		//otestujeme, ci je vygenerovana matica invertovatelna
		if (gauss(temp_matrix) == m)
			break;
	}
	random(b_S, m); //nahodny vektor

	//cout << "bs" << b_S << endl;
	//cout << "As" << A_S << endl;

	Vec<Mat<GF2E>> polynomy_Q_S; //kvadraticke casti polynomov po aplik S
	Vec<Vec<GF2E>> polynomy_L_S; //linearne casti polynomov po aplik S
	Vec<GF2E> polynomy_A_S; //absolutne casti polynomov po aplik S

	for (int i = 0; i < m; i++) {
		Mat<GF2E> q;
		Vec<GF2E> l;
		GF2E a;
		a = b_S[i];

		for (int j = 0; j < m; j++) {
			if (j == 0) {
				q = A_S[j][i] * polynomy_Q_T[j];
				l = A_S[j][i] * polynomy_L_T[j];
				a = a + A_S[j][i] * polynomy_A_T[j];
			}
			else {
				q = q + (A_S[j][i] * polynomy_Q_T[j]);
				l = l + (A_S[j][i] * polynomy_L_T[j]);
				a = a + (A_S[j][i] * polynomy_A_T[j]);
			}

		}
		polynomy_Q_S.append(q);
		polynomy_L_S.append(l);
		polynomy_A_S.append(a);

	}


	/*vypis jednotlivych polynomov po trasformacii S a T
	for (long k = 0; k < m; k++)
	{
		cout << "Polynom cislo: " << k << endl;
		//kvadraticka cast, linearna cast, absolutna cast
		cout << polynomy_Q_S[k] << endl;
		cout << polynomy_L_S[k] << endl;
		cout << polynomy_A_S[k] << endl;
		cout << "**********" << endl;
	}
	*/
	sk.A_S = A_S; sk.b_S = b_S; sk.A_T = A_T; sk.b_T = b_T;
	sk.Q = polynomy_Q; sk.L = polynomy_L; sk.A = polynomy_A;

	pk.Q = polynomy_Q_S; pk.L = polynomy_L_S; pk.A = polynomy_A_S;
}

int sign(Vec<GF2E>& podpis, privateKey& sk, Vec<GF2E>& dokument, int n_variables, int m_poly) {
	
	Vec<GF2E> dokument_inverzia_S;
	Vec<GF2E> x; //vektor neurcitych
	Vec<Vec<GF2E>> Y;
	Mat<GF2E> LS;
	Vec<GF2E> PS;
	Vec<Vec<GF2E>> Z;
	Vec<Vec<GF2E>> riesenia;
	x.SetLength(m_poly + n_variables);
	
	LS.SetDims(m_poly ,  m_poly);

	//inverzia transformacie S
	dokument_inverzia_S = (dokument - sk.b_S) * inv(sk.A_S);
	//inverzia UOV trapdooru
	int count_vinegars =0;
	while (1)
	{
		clear(x);
		count_vinegars+=1;
		
		for (int j =  m_poly; j < n_variables + m_poly; j++) {
			x[j] = random_GF2E();
		}
		
		Y.kill();
		for (long i = 0; i < m_poly; i++) {
			Y.append(sk.Q[i] * x);
		}

		Z.kill();
		for (int i = 0; i < m_poly; i++) {
			Z.append(Y[i] + sk.L[i]);
		}

		for (long i = 0; i < m_poly; i++) {
			for (long j = 0; j < m_poly; j++) {
				LS[i][j] = Z[i][j];
			}
		}
		PS.kill();
		for (long i = 0; i < m_poly; i++) {
			GF2E temp = dokument_inverzia_S[i] - sk.A[i] - x * Z[i];
			PS.append(temp);
		}

		riesenia.kill();
		if (0 == riesenie_sustavy_GF2E(riesenia, LS, PS))
			break;
	}
	for (long i = 0; i < m_poly; i++)
		x[i] = riesenia[0][i];
	//inverzia transformacie T
	podpis = (x - sk.b_T) * inv(sk.A_T);
	cout << "New vinegra were genrated "<< count_vinegars << " times"<<endl;
	return count_vinegars;
}

int verify(Vec<GF2E>& podpis, Vec<GF2E>& dokument, publicKey& pk, long m_poly)
{
	Vec<GF2E> solutions;
	GF2E res;

	for (long k = 0; k < m_poly; k++)
	{
		res = podpis * pk.Q[k] * podpis + pk.L[k] * podpis + pk.A[k];
		solutions.append(res);
	}
	if (solutions == dokument)
	{
		cout << "Podpis je platny" << endl;
		return 1;
	}
	else
	{
		cout << "Podpis nie je platny" << endl;
		return 0;
	}
}










