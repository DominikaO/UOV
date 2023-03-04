#include<stdlib.h>
#include<time.h>
#include<NTL/mat_GF2.h>
#include<NTL/vec_vec_GF2.h>
#include<NTL/vec_GF2.h>
#include<NTL/mat_GF2E.h>
#include<NTL/vec_GF2E.h>
#include<NTL/ZZ.h>

NTL_CLIENT;

int riesenie_sustavy_GF2(vec_vec_GF2& riesenia, mat_GF2& A, vec_GF2& b)
{
	//vektor pravych stran b musi mat tolko prvkov, kolko ma matica A riadkov
	if (b.length() != A.NumRows())
	{
		cout << "Nesedia dimenzie A a b" << endl;
		return -2;
	}

	long pocet_prem;
	bool bool_temp = true;
	vec_vec_GF2 pre_mat_temp;
	vec_vec_GF2 pre_mat;
	int i, j, k, l;
	GF2 det;
	GF2 one = conv<GF2>(1);
	mat_GF2 mat, mat2, ker, img, mat3;
	vec_GF2 vec_temp, right_side, solution;
	pocet_prem = A.NumCols();
	vec_temp.SetLength(pocet_prem + 1);
	right_side.SetLength(pocet_prem);

	pre_mat = A._mat__rep;
	//spojime maticu A a vektor b do vektora vektorov "pre_mat"
	//ktory predstavuje maticu (A|b), t.j. cely linearny system aj s pravou stranou
	for (int i = 0; i < A.NumRows(); i++)
	{
		pre_mat[i].SetLength(pocet_prem + 1);
		pre_mat[i][pocet_prem] = b[i];
	}
	//pre_mat doplnime, aby mal tolko riadkov, kolko "premennych"
	for (int i = A.NumRows(); i < pocet_prem; i++)
	{
		pre_mat.append(vec_temp);
	}

	conv(mat, pre_mat);
	mat2.SetDims(pocet_prem, pocet_prem);

	gauss(mat);

	vec_temp.kill();
	vec_temp.SetLength(pocet_prem + 1);
	vec_temp[pocet_prem] = to_GF2(1);

	//ak sustava obsahuje riadok [0 0 ...  0 1] tak neexistuje riesenie
	for (i = mat._mat__rep.length() - 1; i >= 0; i--)
	{
		if (mat._mat__rep[i][pocet_prem - 1] == to_GF2(1))
		{
			break;
		}
		if (mat._mat__rep[i] == vec_temp)
		{
			append(riesenia, solution);
			//exit(0);
			return -1;
		}

	}
	//rozdelenie matice tvaru n*(n+1) na stvorcovu maticu -mat2- a pravu stranu -right_side-
	for (k = 0; k < pocet_prem; k++)
	{
		for (l = 0; l < pocet_prem; l++)
		{
			mat2[k][l] = mat[k][l];
		}
		right_side[k] = mat[k][l];
	}

	//Existuje jednoznacne riesenie ak je hodnost rovna poctu premennych
	if (gauss(mat2) == pocet_prem)
	{
		solve(det, mat2, solution, right_side);
		append(riesenia, solution);
	}

	//Neexistuje jednoznacne riesenie
	else
	{
		image(img, mat2);	//baza matice mat2
		kernel(ker, transpose(img));		//jadro transponovanej bazy (vlastne zoznam linearne zavislych premennych
		gauss(ker);						//uprava zoznamu premennych na gaussov tvar

		long pocet_volenych = (gauss(ker));
		long iter1, iter2, temp;
		ZZ iterator;

		//cyklus, ktory prejde vsetky moznosti volby premennych, ktore mozno volit
		for (iterator = to_ZZ(0); iterator < power_ZZ(2, pocet_volenych); iterator += 1)
		{
			//pre_mat_temp = pre_mat;
			pre_mat_temp = mat._mat__rep;
			temp = -1;
			//v tomto cykle sa pridavaju rovnice pre volitelne premenne, pricom ich hodnoty su podla "iterator"
			for (iter1 = 0; iter1 < pocet_volenych; iter1++)
			{
				vec_temp.kill();
				vec_temp.SetLength(pocet_prem + 1);
				for (iter2 = temp + 1; iter2 <= ker.NumCols() - 1; iter2++)
				{
					if (ker[iter1][iter2] == one)
					{
						vec_temp[iter2] = one;
						vec_temp[pocet_prem] = to_GF2(bit(iterator, iter1));
						temp = iter2;
						break;
					}
				}
				append(pre_mat_temp, vec_temp);
			}
			//riesenie sustavy rovnic
			conv(mat3, pre_mat_temp);
			gauss(mat3);
			for (k = 0; k < pocet_prem; k++)
			{
				for (l = 0; l < pocet_prem; l++)
				{
					mat2[k][l] = mat3[k][l];
				}
				right_side[k] = mat3[k][l];
			}
			solve(det, mat2, solution, right_side);
			append(riesenia, solution);
		}
	}
	return 0;
}

void iterator_to_GF2E(GF2E& prvok, ZZ& iterator_value, long iter_position, long stupen_rozsirenia)
{
	GF2X temp_polynom;
	for (long i = 0; i < stupen_rozsirenia; i++)
	{
		SetCoeff(temp_polynom, i, bit(iterator_value,iter_position*stupen_rozsirenia + i));
	}
	prvok = conv<GF2E>(temp_polynom);
}

int riesenie_sustavy_GF2E(Vec<Vec<GF2E>>& riesenia, mat_GF2E& A, vec_GF2E& b)
{
	//vektor pravych stran b musi mat tolko prvkov, kolko ma matica A riadkov
	if (b.length() != A.NumRows())
	{
		cout << "Nesedia dimenzie A a b" << endl;
		return -2;
	}
		
	long pocet_prem; 
	long stupen_rozsirenia; 
	bool bool_temp = true;
	Vec<Vec<GF2E>> pre_mat_temp;
	Vec<Vec<GF2E>> pre_mat;
	int i,j,k,l;
	GF2E det;
	mat_GF2E mat, mat2,ker,img, mat3;
	vec_GF2E vec_temp, right_side, solution;
	
	pocet_prem = A.NumCols();
	stupen_rozsirenia = GF2E::degree();
	
	vec_temp.SetLength(pocet_prem+1);  
	right_side.SetLength(pocet_prem);

	pre_mat = A._mat__rep;
	
	//spojime maticu A a vektor b do vektora vektorov "pre_mat"
	//ktory predstavuje maticu (A|b), t.j. cely linearny system aj s pravou stranou
	for (int i = 0; i < A.NumRows(); i++)
	{
		pre_mat[i].SetLength(pocet_prem+1);
		pre_mat[i][pocet_prem] = b[i];
	}
	//pre_mat doplnime, aby mal tolko riadkov, kolko "premennych"
	for (int i = A.NumRows(); i < pocet_prem; i++)
	{
		pre_mat.append(vec_temp);
	}
	
	//vytvorime maticu "mat" ako maticu mat = (A|b),t.j. spojenie A a pravych stran b
	MakeMatrix(mat, pre_mat);	
	mat2.SetDims(pocet_prem, pocet_prem);
	
	//Gaussovska eliminacia (GEM) systemu (A|b)
	gauss(mat);

	long existuje_riesenie = 1;	
	//ak sustava po GEM obsahuje niekde riadok v tvare [0 0 ...  0 c], kde c nie je nula, tak neexistuje riesenie
	for (i=pocet_prem-1; i>=0; i--)
	{
		if (!IsZero(mat[i][pocet_prem]))
		{
			existuje_riesenie = 0;
			for (j = 0; j < pocet_prem; j++)
			{
				if (!IsZero(mat[i][j]))
				{
					existuje_riesenie = 1;
					break;
				}		
			}
			if (existuje_riesenie == 0) //matica obsahuje riadok v tvare [0 0 ... 0 c]
				return -1;
			else
				break;
		}
	}

	//rozdelenie matice mat na stvorcovu maticu -mat2- a pravu stranu -right_side-
	for (k=0; k<pocet_prem; k++)
	{	
		for (l=0; l<pocet_prem; l++)
		{
			mat2[k][l] = mat[k][l];		
		}
		right_side[k] = mat[k][l];
	}


	//Existuje jednoznacne riesenie ak je hodnost rovna poctu premennych
	if (gauss(mat2) == pocet_prem)
	{
		solve(det,mat2,solution,right_side);
		append(riesenia,solution);
	}

	//Neexistuje jednoznacne riesenie
	else
	{
		image(img,mat2);	//baza matice mat2
		kernel(ker,transpose(img));		//jadro transponovanej bazy (vlastne zoznam linearne zavislych premennych
		gauss(ker);						//uprava zoznamu premennych na gaussov tvar
		
		long pocet_volenych = (gauss(ker));
		long iter1,iter2,temp;
		ZZ iterator;
		//cyklus, ktory prejde vsetky moznosti volby premennych, ktore mozno volit
		for (iterator = to_ZZ(0); iterator < power_ZZ(2,pocet_volenych*stupen_rozsirenia); iterator+=1)
		{
			//pre_mat_temp = pre_mat;
			pre_mat_temp = mat._mat__rep;
			temp = -1;
			//v tomto cykle sa pridavaju rovnice pre volitelne premenne, pricom ich hodnoty su podla "iterator"
			for (iter1 = 0; iter1 < pocet_volenych; iter1++)
			{
				vec_temp.kill();
				vec_temp.SetLength(pocet_prem+1);
				for (iter2=temp+1; iter2<=ker.NumCols()-1; iter2++)
				{
					if (!IsZero(ker[iter1][iter2]))
					{
						vec_temp[iter2] = ker[iter1][iter2];
						iterator_to_GF2E(vec_temp[pocet_prem], iterator, iter1, stupen_rozsirenia);
						temp = iter2;
						break;
					}
				}
				append(pre_mat_temp,vec_temp);
			}
			//riesenie sustavy rovnic
			MakeMatrix(mat3,pre_mat_temp);
			gauss(mat3);
			for (k=0; k<pocet_prem; k++)
			{
				for (l=0; l<pocet_prem; l++)
				{
					mat2[k][l] = mat3[k][l];			
				}
				right_side[k] = mat3[k][l];
			}
			solve(det,mat2,solution,right_side); 
			if (!IsZero(det)) append(riesenia,solution);
		}
	}

	return 0;
}
