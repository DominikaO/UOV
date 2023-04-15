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
#include "UOV_pert.h"
#include "UOV.h"
#include "message_hash.h"
#include "sha256.h"
#include "sha512.h"


using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function

using namespace std;

NTL_CLIENT

int main()
{
	GF2X modulus;
	long mod = 6;
	BuildIrred(modulus, mod);
	GF2E::init(modulus);

	Vec<Mat<GF2E>> Q; //kvadraticke casti polynomov
	Vec<Vec<GF2E>> L; //linearne casti polynomov
	Vec<GF2E> A; //absolutne casti polynomov
	Vec<GF2E> h;
	Vec<GF2E> riesenia;
	Vec<Vec<GF2E>> vsetky_moznosti;
	long v = 14; long o = 12; long t = 3;
	//generuj_vsetky_moznosti(vsetky_moznosti, t, mod);

	publicKey_p pk;
	privateKey_p sk;
	KeyGen_p(pk, sk, o, v, t);
	cerr << "keyGen done!" << endl;

	for (int opakovanie = 0; opakovanie < 1; opakovanie++) {
		cout << opakovanie << endl;
		Vec<GF2E> dokument;
		Vec<GF2E> podpis;
		//dokument = random_vec_GF2E(o);
		//TODO LIST SUBOROV
		hash_file512(dokument, "D:\\FEI\\Diplomovka\\Projekt\\UOV\\tst_files\\EPH.pdf", o, mod);

		cout << "--------------------------------------------------- UOV pert v1---------------------------------------------------" << endl;
		//sign_p_random(podpis, sk, dokument, v, o, t);
		sign_p(podpis, sk, dokument, v, o, t, vsetky_moznosti);
		verify_p(podpis, dokument, pk, o);

		cout << "--------------------------------------------------- UOV pert v2---------------------------------------------------" << endl;
		sign_p_v2(podpis, sk, dokument, v, o, t);
		verify_p(podpis, dokument, pk, o);
	}	

	return 0;
}