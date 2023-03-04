#ifndef utoky_h
#define utoky_h

#include<NTL/mat_GF2.h>
#include<NTL/vec_vec_GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2E.h>
#include<NTL/vec_GF2E.h>

NTL_CLIENT;
//funkcie vracaju -1 ak neexistuje riesenie
//funkcie vracaju -2 ak nesedia rozmery A a b
//funkcie vracaju 0 ak existuje riesenie
int riesenie_sustavy_GF2(vec_vec_GF2& riesenia, mat_GF2& A, vec_GF2& b);
int riesenie_sustavy_GF2E(Vec<Vec<GF2E>>& riesenia, mat_GF2E& A, vec_GF2E& b);

#endif
