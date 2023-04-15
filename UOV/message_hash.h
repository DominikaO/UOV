#ifndef messagehash_H
#define messagehash_H
#include <stdlib.h>
#include <iostream>
#include <memory.h>
#include "message_hash.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <bitset>
#include <string>
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
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <string>
#include "sha256.h"
#include "sha512.h"

void hash_file256(Vec<GF2E>& digest, const char* filename, long oil, long mod);
void hash_file512(Vec<GF2E>& digest, const char* filename, long oil, long mod);

#endif   // messagehash_H