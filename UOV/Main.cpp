#include<NTL/ZZ.h> //pre datovy typ ZZ
#include<NTL/ZZ_p.h>
#include<NTL/vec_ZZ_p.h>
#include<NTL/mat_ZZ_p.h>
#include<NTL/GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2.h>
#include <iostream>
#include <tuple>
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
#include <vector>

using std::cerr;
using std::endl;
#include <fstream>
using std::ofstream;
#include <cstdlib> // for exit function
#include <time.h> //for time measurement
using namespace std;

NTL_CLIENT

int main()
{
	struct timespec begin, end;
	long seconds, nanoseconds;
	GF2X modulus;
	long mod;
	double UOV,UOVp,UOVp2,keygen,keygenp,ver,verp;
	int vin,vin1,vin2;
	cout <<"Modulo ";
	cin >> mod;
	BuildIrred(modulus, mod);
	GF2E::init(modulus);

	Vec<Mat<GF2E>> Q; //kvadraticke casti polynomov
	Vec<Vec<GF2E>> L; //linearne casti polynomov
	Vec<GF2E> A; //absolutne casti polynomov
	Vec<GF2E> h;
	Vec<GF2E> riesenia;
	Vec<Vec<GF2E>> vsetky_moznosti;
	long v = 4; long o = 3; long t;
	cout << "Insert vinegar,oil,t (example: 56 48 ): ";
	cin >> v >> o;
	//cout << "Generating all possible solutions for UOV pertubation v1 ..." << endl;
	//generuj_vsetky_moznosti(vsetky_moznosti, t, mod);

	publicKey_p pk;
	privateKey_p sk;

	publicKey puk;
        privateKey prk;
	for (t = 1; t <= 8; t++){
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
	KeyGen_p(pk, sk, o, v, t);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
	seconds = end.tv_sec - begin.tv_sec;
	nanoseconds = end.tv_nsec - begin.tv_nsec;
	cout << "Time keygen UOVp: " << seconds + nanoseconds*1e-9 <<endl;
	keygenp=seconds + nanoseconds*1e-9;
	
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
	KeyGen(puk,prk,o,v);
	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
	seconds = end.tv_sec - begin.tv_sec;
	nanoseconds = end.tv_nsec - begin.tv_nsec;
	cout << "Time keygen UOV: " << seconds + nanoseconds*1e-9 <<endl;
	keygen=seconds + nanoseconds*1e-9;

	string f;
	//cout << "Enter filename for output stats: ";
        //cin >> f;
	
	cout << "t= "<< t<<endl;
	f = "new_time" + to_string(t) + ".csv";
	ofstream file(f);

	// Write the column headers to the file
	file << "UOV,UOV_vinegars,UOVp,UOVp_vinegars, UOV2, UOVp2_vinegars,keygen,keygenp\n";

	for (int opakovanie = 0; opakovanie < 100; opakovanie++) {
		cout << "Iteracia č."<< opakovanie << endl;
		//
		Vec<GF2E> dokument;
		Vec<GF2E> podpis;

		//cout << "Sign and verify" <<endl;
		const char* filename = ("file" + to_string(opakovanie) + ".txt").c_str();
		//cout << "File " << filename << endl;
		hash_file512(dokument,filename, o, mod);

		//cout << "--------------------------------------------------- UOV pert v1---------------------------------------------------" << endl;
		/*sign_p_random(podpis, sk, dokument, v, o, t);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		vin1 = sign_p(podpis, sk, dokument, v, o, t, vsetky_moznosti);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
		seconds = end.tv_sec - begin.tv_sec;
		nanoseconds = end.tv_nsec - begin.tv_nsec;
		cout << "Time UOVp1: " << seconds + nanoseconds*1e-9 << endl;
		UOVp=seconds + nanoseconds*1e-9;
		verify_p(podpis, dokument, pk, o);
*/

		//cout << "--------------------------------------------------- UOV pert v2---------------------------------------------------" << endl;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		vin2 = sign_p_v2(podpis, sk, dokument, v, o, t);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
		seconds = end.tv_sec - begin.tv_sec;
		nanoseconds = end.tv_nsec - begin.tv_nsec;
		cout << "Time UOVp2: " << seconds + nanoseconds*1e-9 <<endl;
		UOVp2=seconds + nanoseconds*1e-9;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		verify_p(podpis,dokument,pk,o);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
		verp =  ((end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec))*1e-9;

		//cout << "--------------------------------------------------- UOV ---------------------------------------------------" << endl;
	        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		vin = sign(podpis, prk, dokument, v, o);
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
		seconds = end.tv_sec - begin.tv_sec;
		nanoseconds = end.tv_nsec - begin.tv_nsec;
		cout << "Time UOV: " << seconds + nanoseconds*1e-9 << endl;
		UOV=seconds + nanoseconds*1e-9;
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &begin);
		verify(podpis, dokument, puk, o);	
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &end);
		ver =  ((end.tv_sec - begin.tv_sec) + (end.tv_nsec - begin.tv_nsec))*1e-9;

		// Write the data to the file
	         file << UOV << ","<< vin << "," << UOVp<< ","<< vin1<< "," << UOVp2 << "," << vin2 << ","<< keygen << ","<< keygenp << ","<<ver<<"," << verp<< endl;
		

	
	}
	
	file.close();
	}
	return 0;
		}
