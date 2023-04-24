#ifndef source_h
#define source_h

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

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

using namespace std;

NTL_CLIENT


Vec<GF2E> riesenia(Vec<Mat<GF2E>> polynomy_Q, Vec<Vec<GF2E>> polynomy_L, Vec<GF2E> polynomy_A, int n, int m);



struct publicKey {
	Vec<Mat<GF2E>> Q;
	Vec<Vec<GF2E>> L;
	Vec<GF2E> A;
};

struct privateKey {
	Vec<Mat<GF2E>> Q;
	Vec<Vec<GF2E>> L;
	Vec<GF2E> A;
	Mat<GF2E> A_T;
	Vec<GF2E> b_T;
	Mat<GF2E> A_S;
	Vec<GF2E> b_S;
};

void generate_random_polynomials(int m, int n, Vec<Mat<GF2E>>& polynomy_Q,
												Vec<Vec<GF2E>>& polynomy_L,
												Vec<GF2E>& polynomy_A);

void KeyGen(publicKey& pk, privateKey& sk, long m_poly, long n_variables);

int sign(Vec<GF2E>& podpis, privateKey& sk, Vec<GF2E>& dokument, int n_variables, int m_poly);

int verify(Vec<GF2E>& podpis, Vec<GF2E>& dokument, publicKey& pk, long m_poly);

#endif

