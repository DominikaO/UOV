#pragma once
#define source_pert_h
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
#include <NTL/GF2XFactoring.h>

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

using namespace std;

NTL_CLIENT




struct publicKey_p {
	Vec<Mat<GF2E>> Q;
	Vec<Vec<GF2E>> L;
	Vec<GF2E> A;
};

struct privateKey_p {
	Vec<Mat<GF2E>> Q;
	Vec<Mat<GF2E>> Q_wo_z;
	Vec<Mat<GF2E>> polynomy_z;
	Vec<Vec<GF2E>> lambdas;
	Vec<Vec<GF2E>> L;
	Vec<GF2E> A;
	Mat<GF2E> A_T;
	Vec<GF2E> b_T;
	Mat<GF2E> A_S;
	Vec<GF2E> b_S;
};



void KeyGen_p(publicKey_p& pk, privateKey_p& sk, long m_poly, long n_variables, long t);



GF2E convert_int_GF2E(long i, long mod);

void generuj_vsetky_moznosti(Vec<Vec<GF2E>>& vsetky_moznosti, long t, long mod);

int sign_p(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t, Vec<Vec<GF2E>>& vsetky_moznosti);

void sign_p_random(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t);

int sign_p_v2(Vec<GF2E>& podpis, privateKey_p& sk, Vec<GF2E>& dokument, int n_variables, int m_poly, int t);

int verify_p(Vec<GF2E>& podpis, Vec<GF2E>& dokument, publicKey_p& pk, long m_poly);

