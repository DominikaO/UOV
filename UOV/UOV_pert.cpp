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

//nasledujuce 3 pridane kvoli libopenf4
#include<libopenf4.h>
#include<vector>
#include<string>

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
	/*
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
		polynomy_Q_T.append(A_T * polynomy_Q_plus_z[k] * transpose(A_T));
		polynomy_L_T.append(b_T * polynomy_Q_plus_z[k] * transpose(A_T) + b_T * transpose(polynomy_Q_plus_z[k]) * transpose(A_T) + polynomy_L[k] * transpose(A_T));
		polynomy_A_T.append(b_T * polynomy_Q_plus_z[k] * b_T + polynomy_L[k] * b_T + polynomy_A[k]);
	}
	/*
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
	/*
		cout << "bs" << b_S << endl;
		cout << "As" << A_S << endl;
	*/
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
	/*

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
	*/
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

int sign_p(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t, Vec<Vec<GF2E>>& vsetky_moznosti) {

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
	int count_new_vinegar = 0;
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
		count_new_vinegar+=1;

		if (count_new_vinegar >= 9) goto riesenie_najdene;
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
	cout <<"Vinegars were generated" << count_new_vinegar << "times" << endl;
	return count_new_vinegar;
}

void sign_p_random(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t) {

	Vec<GF2E> dokument_inverzia_S;
	Vec<GF2E> x; //vektor neurcitych
	Vec<Vec<GF2E>> Y;
	Mat<GF2E> LS;
	Vec<GF2E> PS;
	Vec<GF2E> c;
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
		for (ZZ pokus_cislo = conv<ZZ>(0); pokus_cislo < power2_ZZ(GF2E::degree() * t); pokus_cislo++)
			//for (long pokus_cislo = 0; pokus_cislo < 1000; pokus_cislo++)
		{
			random(c, t);
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
		//hladanie pivota v i-tom stlpci staci ist od l-teho riadku dalej
		for (j = l; j < r; j++)
		{
			if (!IsZero(M_vec[j][i]))
			{
				//nasli sme pivota
				break;
			}
		}
		if (j == r)	//nulovy stlpec, nenasiel sa pivot
			continue;
		//vymenime j-ty a l-ty riadok
		swap(M_vec[l], M_vec[j]);
		//nulovanie prvkov v stlpci i
		for (j = 0; j < r; j++)
		{
			if (j == l) continue;
			if (IsZero(M_vec[j][i]))
				continue;
			else
			{
				//nulovanie
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

void GB(Vec<GF2E>& riesenie, long num_vars, Vec<Mat<GF2E>>& pols_quad, Vec<Vec<GF2E>>& pols_lin, Vec<GF2E>& pols_abs)
{
	long i, j, k, l;
	vector<string> polynomialArray;
	vector<string> variableName;
	for (i = 0; i < num_vars; i++)
	{
		variableName.push_back('x' + to_string(i));
	}

	//tvorba polynomov
	string polynom;
	string koeficient_str;
	GF2X koeficient;
	for (i = 0; i < num_vars; i++)
	{
		polynom.clear();
		//absolutny koeficient
		koeficient = rep(pols_abs[i]);
		koeficient_str.clear();
		if (!IsZero(koeficient))
		{
			for (l = 0; l <= deg(koeficient); l++)
			{
				if (!IsZero(koeficient[l]))
				{
					if (l == 0)
						koeficient_str += "1+";
					else {
						if (l == 1)
							koeficient_str += "t+";
						else
							koeficient_str += "t^" + to_string(l) + "+";
					}
				}
			}
			koeficient_str.pop_back();
			polynom += "(" + koeficient_str + ")+";
		}
		//linearne koeficienty
		for (j = num_vars - 1; j >= 0; j--)
		{
			koeficient = rep(pols_lin[i][j]);
			koeficient_str.clear();
			if (!IsZero(koeficient))
			{
				for (l = 0; l <= deg(koeficient); l++)
				{
					if (!IsZero(koeficient[l]))
					{
						if (l == 0)
							koeficient_str += "1+";
						else {
							if (l == 1)
								koeficient_str += "t+";
							else
								koeficient_str += "t^" + to_string(l) + "+";
						}
					}
				}
				koeficient_str.pop_back();
				polynom += "(" + koeficient_str + ")*x" + to_string(j) + "+";
			}
		}
		//kvadraticke koeficienty
		for (k = num_vars - 1; k >= 0; k--)
		{
			for (j = k; j >= 0; j--)
			{
				koeficient = rep(pols_quad[i][j][k]);
				koeficient_str.clear();
				if (!IsZero(koeficient))
				{
					for (l = 0; l <= deg(koeficient); l++)
					{
						if (!IsZero(koeficient[l]))
						{
							if (l == 0)
								koeficient_str += "1+";
							else {
								if (l == 1)
									koeficient_str += "t+";
								else
									koeficient_str += "t^" + to_string(l) + "+";
							}
						}
					}
					koeficient_str.pop_back();
					if (j != k)
						polynom += "(" + koeficient_str + ")*x" + to_string(j) + "*x" + to_string(k) + "+";
					else
						polynom += "(" + koeficient_str + ")*x" + to_string(j) + "^2+";
				}
			}
		}
		polynom.pop_back();
		polynomialArray.emplace_back(polynom);
	}

	//na zaver modulus
	koeficient = GF2E::modulus();
	koeficient_str.clear();
	if (!IsZero(koeficient))
	{
		for (l = 0; l <= deg(koeficient); l++)
		{
			if (!IsZero(koeficient[l]))
			{
				if (l == 0)
					koeficient_str += "1+";
				else {
					if (l == 1)
						koeficient_str += "t+";
					else
						koeficient_str += "t^" + to_string(l) + "+";
				}
			}
		}
		koeficient_str.pop_back();
	}

	//Tu musi byt cyklus, v ktorom sa budu hrubou silou dosadzat hodnoty pre povedzme x0, x0 = c, doplni sa vztah (c) + x0
	//a pozrie sa, ci to naslo riesenie, alebo nie
	vector<string> polynomialArray_x0;
	vector<string> basis;

	Vec<Vec<GF2E>> vsetky_moznosti;
	string GF2E_str;
	generuj_vsetky_moznosti(vsetky_moznosti, 1, GF2E::degree());
	for (auto a : vsetky_moznosti)
	{
		polynom.clear();
		polynomialArray_x0 = polynomialArray;
		for (j = 0; j >= 0; j--)
		{
			koeficient = rep(a[0]);
			GF2E_str.clear();
			if (!IsZero(koeficient))
			{
				for (l = 0; l <= deg(koeficient); l++)
				{
					if (!IsZero(koeficient[l]))
					{
						if (l == 0)
							GF2E_str += "1+";
						else {
							if (l == 1)
								GF2E_str += "t+";
							else
								GF2E_str += "t^" + to_string(l) + "+";
						}
					}
				}
				GF2E_str.pop_back();
				polynom += "(" + GF2E_str + ")+x" + to_string(j);
			}
			else
			{
				polynom += "x0";
			}
		}
		polynomialArray_x0.emplace_back(polynom.c_str());
		basis = groebnerBasisGF2ExtensionF4(koeficient_str.c_str(), num_vars, variableName, "t", polynomialArray_x0, 1, 0);


		// Spracovanie redukovanej GB
		if (basis.size() == num_vars)
		{
			riesenie.kill();
			int solution_found = 1;
			//ASI sa naslo riesenie

			for (i = 0; i < num_vars; i++)
			{
				string searched;
				string element;
				Vec<GF2> element_bits;
				searched += "((+1)*x" + to_string(i) + "^1)";
				if (basis[i].find(searched) != 0)
				{
					//baza nie je vhodna
					solution_found = 0;
					break;
				}
				else
				{
					if (basis[i].length() == searched.length())
					{
						riesenie.append(GF2E::zero());
					}
					else
					{
						element = basis[i].substr(searched.length() + 5);
						element_bits.kill(); element_bits.SetLength(GF2E::degree());
						for (j = GF2E::degree() - 1; j >= 2; j--)
						{
							searched = "t^" + to_string(j);
							if (element.find(searched) == 0)
							{
								element_bits[j] = 1;
								element = element.substr(searched.length() + 1);
							}
						}
						//j == 1
						if (element[0] == '+')
							element.erase(0, 1);
						if (element[0] == 't')
						{
							element_bits[1] = 1;
							element.erase(0, 1);
						}
						//j == 0
						if (element[0] == '+')
							element.erase(0, 1);
						if (element[0] == '1')
						{
							element_bits[0] = 1;
						}
						riesenie.append(conv<GF2E>(conv<GF2X>(element_bits)));
					}
				}
			}
			if (solution_found)
				break;
		}
	}
}
int sign_p_v2(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t) {

	Vec<GF2E> dokument_inverzia_S;
	Vec<GF2E> x; //vektor neurcitych
	Vec<Vec<GF2E>> Y;
	Mat<GF2E> LS;
	Mat<GF2E> LS_PS;
	Vec<GF2E> PS;
	Vec<Vec<GF2E>> Z;
	Vec<Vec<GF2E>> riesenia;
	x.SetLength(m_poly + n_variables);

	//LS.SetDims(m_poly, m_poly+t);
	LS_PS.SetDims(m_poly, m_poly + t + 1);
	//inverzia transformacie S
	dokument_inverzia_S = (dokument - sk.b_S) * inv(sk.A_S);
	//inverzia UOV trapdooru
	int count_new_vinegar = 0;
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
				LS_PS[i][j] = Z[i][j];
			}
		}
		for (long i = 0; i < m_poly; i++) {
			GF2E temp = dokument_inverzia_S[i] - sk.A[i] - x * Z[i];
			LS_PS[i][m_poly + t] = temp;
		}

		//upravime LS pridanim novych premennych z_1, z_2, ..., z_t
		for (long i = 0; i < m_poly; i++) {
			for (long j = m_poly; j < m_poly + t; j++) {
				LS_PS[i][j] = sk.lambdas[i][j - m_poly];
			}
		}
		count_new_vinegar +=1;

		//kontrola, ci hodnost(LS) == pocet_olejov
		if (GEM(LS_PS) == -1)
		{
			continue;
		}
		else
		{
			//cout << LS_PS << endl; //break;
		}

		//otestujeme, ci v i-tom riadku je prvok na i-tej pozicii
		//ak nie, musime opakovat

		//na zaklade vyjadrenia x_1 = ()*z1 + ()*z2 + ... toto dosadime do perturbacnych polynomov
		//a dostaneme sustavu t kvadratickych rovnic o t premennych z1,z2,...,zt

		Vec<Mat<GF2E>> poly_z_Q;
		Vec<Vec<GF2E>> poly_z_L;
		Vec<GF2E> poly_z_A;
		Vec<GF2E> row_j, row_k;

		for (long i = 0; i < t; i++)
		{
			Mat<GF2E> mat_temp;
			mat_temp.kill();
			mat_temp.SetDims(t, t);
			poly_z_Q.append(mat_temp);

			Vec<GF2E> vec_temp;
			vec_temp.kill();
			vec_temp.SetLength(t);
			vec_temp[i] = GF2E(1);
			poly_z_L.append(vec_temp);

			poly_z_A.append(GF2E(0));

			//vezmeme i-ty perturbacny polynom a dosadime za x_1, x_2, ..., x_n
			for (int j = 0; j < m_poly; j++) {
				for (int k = 0; k < m_poly; k++) {
					if (!IsZero(sk.polynomy_z[i][j][k]))
					{
						row_j = LS_PS[j];
						row_k = LS_PS[k];

						for (int zj = m_poly; zj < m_poly + t; zj++)
						{
							for (int zk = m_poly; zk < m_poly + t; zk++)
							{
								poly_z_Q[i][zj - m_poly][zk - m_poly] += sk.polynomy_z[i][j][k] * row_j[zj] * row_k[zk];
							}
						}

						for (int zj = m_poly; zj < m_poly + t; zj++)
						{
							poly_z_L[i][zj - m_poly] += sk.polynomy_z[i][j][k] * row_j[zj] * row_k[m_poly + t];
							poly_z_L[i][zj - m_poly] += sk.polynomy_z[i][j][k] * row_j[m_poly + t] * row_k[zj];
						}

						poly_z_A[i] += (sk.polynomy_z[i][j][k] * row_j[m_poly + t] * row_k[m_poly + t]);
					}
				}
			}
		}

		//v premennych poly_z_Q, poly_z_L, poly_z_A mame ulozenu sustavu kvadratickych rovnic, ktoru potrebujeme vyriesit nejakym sofistikovanym algoritmom!
		Vec<GF2E> riesenie_z_i;

		//vypis systemu ktory riesime cez F4

		//"ztrojuholnikujeme matice poly_z_Q
		for (int i = 0; i < t; i++)
		{
			for (int k = 0; k < t; k++)
			{
				for (int j = 0; j < k; j++)
				{
					poly_z_Q[i][j][k] += poly_z_Q[i][k][j];
					poly_z_Q[i][k][j] = GF2E(0);
				}
			}
		}
		
				//cout << poly_z_Q << endl;
				//cout << poly_z_L << endl;
				//cout << poly_z_A << endl;
		
		GB(riesenie_z_i, t, poly_z_Q, poly_z_L, poly_z_A);
		//cout << "riesenie: " << riesenie_z_i << endl;

		if (riesenie_z_i.length() != t)
		{
			//skus vygenerovat nove octy
			continue;
		}

		for (int k = 0; k < m_poly; k++) {
			GF2E x_i = LS_PS[k][m_poly + t];
			for (int j = 0; j < t; j++)
			{
				x_i += LS_PS[k][m_poly + j] * riesenie_z_i[j];
			}
			x[k] = x_i;
		}
		//cout << x << endl;
		break;
	}

riesenie_najdene:
	//inverzia transformacie T
	podpis = (x - sk.b_T) * inv(sk.A_T);
		cout <<"New vinegars generated " <<count_new_vinegar <<" times" <<endl;
	return count_new_vinegar;

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
