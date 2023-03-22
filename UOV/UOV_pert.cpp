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
#include "UOV.h"
#include "UOV_pert.h"
#include <NTL/GF2XFactoring.h>

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

using namespace std;

NTL_CLIENT




void KeyGen_p(publicKey_p& pk, privateKey_p& sk, long m_poly, long n_variables, long t)
{


	Vec<Mat<GF2E>> polynomy_Q; //kvadraticke casti polynomov
	Vec<Vec<GF2E>> polynomy_L; //linearne casti polynomov
	Vec<GF2E> polynomy_A; //absolutne casti polynomov

	long m = m_poly; //pocet polynomov a zaroven pocet olej
	long n = n_variables; //pocet neurcitych ocot
	generate_random_polynomials(m, n, polynomy_Q, polynomy_L, polynomy_A);

	Vec<Mat<GF2E>> polynomy_z;  //kvadraticke casti perturbacnych polynomov z_1, z_2, ..., z_t
	polynomy_z.SetLength(t);
	for (long i = 0; i < t; i++)
	{
		Mat<GF2E> Q;
		Q.SetDims(m_poly + n_variables, m_poly + n_variables);
		for (int j = 0; j < m_poly; j++) {
			for (int k = j; k < m_poly; k++) {
				Q[j][k] = random_GF2E();
			}
		}
		polynomy_z[i] = Q;
	}
	Vec<Vec<GF2E>> lambdas;
	lambdas.SetLength(m_poly);
	for (long i = 0; i < m_poly; i++)
	{
		lambdas[i] = random_vec_GF2E(t);
	}

	Vec<Mat<GF2E>> polynomy_Q_plus_z = polynomy_Q;
	for (long i = 0; i < m_poly; i++)
	{
		for (long j = 0; j < t; j++)
		{
			polynomy_Q_plus_z[i] += polynomy_z[j] * lambdas[i][j];
		}
	}

	//vypis jednotlivych polynomov
	for (long k = 0; k < m; k++)
	{
		cout << "Polynom cislo: " << k + 1 << endl;
		//kvadraticka cast, linearna cast, absolutna cast
		cout << polynomy_Q_plus_z[k] << endl;
		cout << polynomy_L[k] << endl;
		cout << polynomy_A[k] << endl;
		cout << "**********" << endl;
	}

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
		polynomy_Q_T.append(A_T * polynomy_Q_plus_z[k] * transpose(A_T));
		polynomy_L_T.append(b_T * polynomy_Q_plus_z[k] * transpose(A_T) + b_T * transpose(polynomy_Q_plus_z[k]) * transpose(A_T) + polynomy_L[k] * transpose(A_T));
		polynomy_A_T.append(b_T * polynomy_Q_plus_z[k] * b_T + polynomy_L[k] * b_T + polynomy_A[k]);
	}
	//vypis jednotlivych polynomov
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

	cout << "bs" << b_S << endl;
	cout << "As" << A_S << endl;

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


	//vypis jednotlivych polynomov po trasformacii S a T
	for (long k = 0; k < m; k++)
	{
		cout << "Polynom cislo: " << k << endl;
		//kvadraticka cast, linearna cast, absolutna cast
		cout << polynomy_Q_S[k] << endl;
		cout << polynomy_L_S[k] << endl;
		cout << polynomy_A_S[k] << endl;
		cout << "**********" << endl;
	}

	sk.A_S = A_S; sk.b_S = b_S; sk.A_T = A_T; sk.b_T = b_T;
	sk.Q = polynomy_Q_plus_z; sk.L = polynomy_L; sk.A = polynomy_A; sk.Q_wo_z = polynomy_Q; sk.lambdas = lambdas; sk.polynomy_z = polynomy_z;

	pk.Q = polynomy_Q_S; pk.L = polynomy_L_S; pk.A = polynomy_A_S;
}



void generuj_vsetky_moznosti(Vec<Vec<GF2E>>& vsetky_moznosti, long t, long mod)
{
	Vec<GF2> temp_vec;
	for (ZZ i = ZZ::zero(); i < power2_ZZ(mod * t); i++)
	{
		Vec<GF2E> vektor;
		for (long j = 0; j < t; j++)
		{
			temp_vec.kill();
			for (long k = 0; k < mod; k++)
			{
				temp_vec.append(conv<GF2>(bit(i, j * mod + k)));
			}
			vektor.append(conv<GF2E>(conv<GF2X>(temp_vec)));
		}
		vsetky_moznosti.append(vektor);
	}
}

void sign_p(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t, Vec<Vec<GF2E>>& vsetky_moznosti) {

	Vec<GF2E> dokument_inverzia_S;
	Vec<GF2E> x; //vektor neurcitych
	Vec<Vec<GF2E>> Y;
	Mat<GF2E> LS;
	Vec<GF2E> PS;
	Vec<Vec<GF2E>> Z;
	Vec<Vec<GF2E>> riesenia;
	x.SetLength(m_poly + n_variables);

	LS.SetDims(m_poly, m_poly);

	//inverzia transformacie S
	dokument_inverzia_S = (dokument - sk.b_S) * inv(sk.A_S);
	//inverzia UOV trapdooru
	while (1)
	{
		clear(x);

		for (int j = m_poly; j < n_variables + m_poly; j++) {
			x[j] = random_GF2E();
		}

		Y.kill();
		for (long i = 0; i < m_poly; i++) {
			//Y.append(sk.Q[i] * x);
			Y.append(sk.Q_wo_z[i] * x);
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

		//prechadzaj vsetkych q^t volieb pre z_1,z_2,...,z_t
		//a odcitaj od pravych stran lambda_1*z_1+lambda_2*z_2...
		for (auto c : vsetky_moznosti)
		{

			riesenia.kill();
			Vec<GF2E> PS_upravene = PS;
			for (long i = 0; i < m_poly; i++)
			{
				PS_upravene[i] -= c * sk.lambdas[i];
			}

			if (0 == riesenie_sustavy_GF2E(riesenia, LS, PS_upravene))
			{

				for (auto riesenie : riesenia)
				{
					riesenie.SetLength(m_poly + n_variables);
					Vec<GF2E> verifikator;
					for (long i = 0; i < t; i++)
					{

						verifikator.append(riesenie * sk.polynomy_z[i] * riesenie);
					}

					if (verifikator == c)
					{
						for (long i = 0; i < m_poly; i++)
							x[i] = riesenie[i];
						goto riesenie_najdene;
					}
				}
			}
		}


	}
riesenie_najdene:
	//inverzia transformacie T
	podpis = (x - sk.b_T) * inv(sk.A_T);

}

long GEM(Mat<GF2E>& M)
{
	long kroky, i, j, k, l;
	long r = M.NumRows();
	long c = M.NumCols();
	if (r < c)
	{
		kroky = r;
	}
	else
	{
		kroky = c;
	}
	Vec<Vec<GF2E>> M_vec = rep(M);
	GF2E pivot; GF2E pivot_inv; GF2E lead;
	for (i = l = 0; i < kroky; i++)
	{
		/*hladanie pivota v i-tom stlpci*/
		/*staci ist od l-teho riadku dalej*/
		for (j = l; j < r; j++)
		{
			if (!IsZero(M_vec[j][i]))
			{
				/*nasli sme pivota*/
				break;
			}
		}
		if (j == r)	//nulovy stlpec, nenasiel sa pivot
			continue;
		/*vymenime j-ty a l-ty riadok*/
		swap(M_vec[l], M_vec[j]);
		/*nulovanie prvkov v stlpci i*/
		for (j = 0; j < r; j++)
		{
			if (j == i) continue;

			if (IsZero(M_vec[j][i]))
				continue;
			else
			{
				/*nulovanie*/
				pivot = M_vec[l][i];
				inv(pivot_inv, pivot);
				pivot_inv *= M_vec[j][i];
				for (k = i; k < c; k++)
				{
					M_vec[j][k] += M_vec[l][k] * pivot_inv;
				}
			}
		}
		l++;
	}

	for (i = 0; i < r; i++)
	{
		if (IsZero(M_vec[i][i]))
		{
			return -1;
		}
		else
		{
			pivot = M_vec[i][i]; inv(pivot_inv, pivot);
			M_vec[i] *= pivot_inv;
		}
	}
	MakeMatrix(M, M_vec);
	return 0;
}

void sign_p_v2(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t, Vec<Vec<GF2E>>& vsetky_moznosti) {

	Vec<GF2E> dokument_inverzia_S;
	Vec<GF2E> x; //vektor neurcitych
	Vec<Vec<GF2E>> Y;
	Mat<GF2E> LS;
	Vec<GF2E> PS;
	Vec<Vec<GF2E>> Z;
	Vec<Vec<GF2E>> riesenia;
	x.SetLength(m_poly + n_variables);

	LS.SetDims(m_poly, m_poly + t);

	//inverzia transformacie S
	dokument_inverzia_S = (dokument - sk.b_S) * inv(sk.A_S);
	//inverzia UOV trapdooru
	while (1)
	{
		clear(x);

		for (int j = m_poly; j < n_variables + m_poly; j++) {
			x[j] = random_GF2E();
		}

		Y.kill();
		for (long i = 0; i < m_poly; i++) {
			//Y.append(sk.Q[i] * x);
			Y.append(sk.Q_wo_z[i] * x);
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

		//upravime LS pridanim novych premennych z_1, z_2, ..., z_t
		for (long i = 0; i < m_poly; i++) {
			for (long j = m_poly; j < m_poly + t; j++) {
				LS[i][j] = sk.lambdas[i][j - m_poly];
			}
		}

		//kontrola, ci hodnost(LS) == pocet_olejov

		if (GEM(LS) == -1)
			continue;
		else
		{
			cout << LS << endl; break;
		}
		//otestujeme, ci v i-tom riadku je prvok na i-tej pozicii
		//ak nie, musime opakovat



		/*
		riesenia.kill();

		if (0 == riesenie_sustavy_GF2E(riesenia, LS, PS_upravene))
		{
			for (auto riesenie : riesenia)
			{
				riesenie.SetLength(m_poly + n_variables);
				Vec<GF2E> verifikator;
				for (long i = 0; i < t; i++)
				{
						verifikator.append(riesenie * sk.polynomy_z[i] * riesenie);
				}
					if (verifikator == c)
				{
					for (long i = 0; i < m_poly; i++)
						x[i] = riesenie[i];
					goto riesenie_najdene;
				}
			}
		}


*/
	}
	/*
riesenie_najdene:
	//inverzia transformacie T
	podpis = (x - sk.b_T) * inv(sk.A_T);
*/
}


int verify_p(Vec<GF2E>& podpis, Vec<GF2E>& dokument, publicKey_p& pk, long m_poly)
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






/*

int main()
{
	GF2X modulus;
	long mod = 3;
	BuildIrred(modulus, mod);
	GF2E::init(modulus);

	Vec<Mat<GF2E>> Q; //kvadraticke casti polynomov
	Vec<Vec<GF2E>> L; //linearne casti polynomov
	Vec<GF2E> A; //absolutne casti polynomov
	Vec<GF2E> h;
	Vec<GF2E> riesenia;


	generuj_vsetky_moznosti();


	//todo exhastive search
	long v = 3; long o = 4; long t = 2;

	publicKey pk;
	privateKey sk;
	KeyGen(pk, sk, o, v, t);



	Vec<GF2E> dokument;
	Vec<GF2E> podpis;
	dokument = random_vec_GF2E(o);
	sign(podpis, sk, dokument, v, o, t);
	verify(podpis, dokument, pk, o);





	SetSeed(ZZ(time(NULL))); // Initialize NTL random number generator with current time

	int m = 3; // Number of polynomials
	int n = 4; // Number of vinegar terms

	Vec<Mat<GF2>> polynomy_Q;
	Vec<Vec<GF2>> polynomy_L;
	Vec<GF2> polynomy_A;

	generate_random_polynomials(m, n, polynomy_Q, polynomy_L, polynomy_A);

	// Print the generated polynomials
	for (int i = 0; i < m; i++) {
		cout << "Polynomial " << i + 1 << ":" << endl;
		cout << "Quadratic part:" << endl << polynomy_Q[i] << endl;
		cout << "Linear part:" << endl << polynomy_L[i] << endl;
		cout << "Absolute part:" << endl << polynomy_A[i] << endl;
		cout << endl;
	}

	return 0;


	return 0;
}

*/
