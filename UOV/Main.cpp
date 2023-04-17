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
	long v = 4; long o = 3; long t = 2;

	generuj_vsetky_moznosti(vsetky_moznosti, t, mod);

	publicKey_p pk;
	privateKey_p sk;
	
	KeyGen_p(pk, sk, o, v, t);
	cerr << "keyGen done!" << endl;
	const char* dir = "\tst_files\\";
	char* filename = (char*)malloc(strlen(dir) + 9);

	for (int opakovanie = 0; opakovanie < 100; opakovanie++) {
		cout << opakovanie << endl;
		Vec<GF2E> dokument;
		Vec<GF2E> podpis;
		//dokument = random_vec_GF2E(o);
		/*TODO LIST SUBOROV
		
		sprintf(filename, "%s%04d.dat", dir, opakovanie);
		FILE* file = fopen(filename, "w");
		
		if (file == NULL) {
			printf("Error opening file\n");
			return 1;
		}

		// Seed the random number generator with the current time
		srand(time(NULL));

		char text[1024];
		for (int i = 0; i < 1024 - 1; i++) {
			text[i] = 'a' + rand() % 26; // Generate a random lowercase letter
		}
		text[1024 - 1] = '\0'; // Terminate the string with a null character

		// Write the text to the file
		fprintf(file, "%s", text);

		// Close the file
		fclose(file);
		*/
		cout << "Hash suboru" << filename << endl;
		hash_file512(dokument,filename, o, mod);

		cout << "--------------------------------------------------- UOV pert v1---------------------------------------------------" << endl;
		//sign_p_random(podpis, sk, dokument, v, o, t);
		sign_p(podpis, sk, dokument, v, o, t, vsetky_moznosti);
		verify_p(podpis, dokument, pk, o);

		cout << "--------------------------------------------------- UOV pert v2---------------------------------------------------" << endl;
		//sign_p_v2(podpis, sk, dokument, v, o, t);
		verify_p(podpis, dokument, pk, o);
	}	

	return 0;
}