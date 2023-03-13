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
	long mod = 3;
	BuildIrred(modulus, mod);
	GF2E::init(modulus);

	Vec<Mat<GF2E>> Q; //kvadraticke casti polynomov
	Vec<Vec<GF2E>> L; //linearne casti polynomov
	Vec<GF2E> A; //absolutne casti polynomov
	Vec<GF2E> h;
	Vec<GF2E> riesenia;
	Vec<Vec<GF2E>> vsetky_moznosti;
	long v = 3; long o = 4; long t = 2;
	generuj_vsetky_moznosti(vsetky_moznosti, t, mod);

	publicKey_p pk;
	privateKey_p sk;
	KeyGen_p(pk, sk, o, v, t);



	Vec<GF2E> dokument;
	Vec<GF2E> podpis;
	dokument = random_vec_GF2E(o);
	cout << "--------------------------------------------------- UOV pert ---------------------------------------------------" << endl;
	sign_p(podpis, sk, dokument, v, o, t, vsetky_moznosti);
	verify_p(podpis, dokument, pk, o);

	publicKey pk2;
	privateKey sk2;
	KeyGen(pk2, sk2, o, v);
	cout << "--------------------------------------------------- UOV ---------------------------------------------------" << endl;
	sign(podpis, sk2, dokument, v, o);
	verify(podpis, dokument, pk2, o);


	return 0;
}