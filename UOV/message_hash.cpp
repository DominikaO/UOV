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
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>



NTL::GF2E stringToGF2E(const std::string& str) {
	NTL::GF2X poly;
	for (int i = 0; i < str.size(); i++) {
		NTL::SetCoeff(poly, i, str[i] == '1' ? 1 : 0);
	}
	NTL::GF2E elem;
	NTL::conv(elem, poly);
	return elem;
}

Vec<GF2E> create_vectors(const char* hash, long oil, long mod) {
	std::string hex(hash, 64);
	std::string binary;

	for (char c : hex) {
		int value = (c <= '9') ? (c - '0') : (c - 'a' + 10);
		std::bitset<4> bits(value);
		binary += bits.to_string();
	}
	//cout << "hash binary " << binary;
	Vec<GF2E> digest;
	for (int i = 0; i < oil; i++) {
		string substr_mod = binary.substr(i*mod, mod);
		
		digest.append(stringToGF2E(substr_mod));
	}
	//cout << digest;
	return digest;
}

Vec<GF2E> create_vectors2(const std::string& hash, long oil, long mod) {
	// Convert hash to binary string
	std::string binary;
	for (char c : hash) {
		int value = (c <= '9') ? (c - '0') : (c - 'a' + 10);
		std::bitset<4> bits(value);
		binary += bits.to_string();
	}

	// Create vectors from binary string
	Vec<GF2E> digest;
	for (int i = 0; i < oil; i++) {
		std::string substr_mod = binary.substr(i * mod, mod);
		digest.append(stringToGF2E(substr_mod));
	}
	//std::cout << digest << std::endl;
	return digest;
}


void hash_file256(Vec<GF2E>& digest, const char *filename, long oil, long mod) {
	

	// Open the file
	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		printf("Unable to open file: %s\n", filename);
		
	}

	// Read the contents of the file into a buffer
	unsigned char buffer[BUFFER_SIZE];
	size_t bytes_read = 0;
	size_t total_bytes_read = 0;
	while ((bytes_read = fread(buffer + total_bytes_read, 1, BUFFER_SIZE - total_bytes_read, fp))) {
		total_bytes_read += bytes_read;
		if (total_bytes_read == BUFFER_SIZE) {
			break;
		}
	}
	

	// Compute the SHA-256 hash of the buffer

	unsigned char hash[32];
	sha256(buffer, total_bytes_read, hash);

	// Convert the hash to a string representation
	char hash_string[65];
	for (int i = 0; i < 32; i++) {
		sprintf(hash_string + 2 * i, "%02x", hash[i]);
	}
	hash_string[64] = '\0';

	//printf("SHA-256 hash of %s: %s\n", filename, hash_string);
	digest =create_vectors(hash_string,oil, mod);
	

	fclose(fp);
	
}


void hash_file512(Vec<GF2E>& digest, const char* filename, long oil, long mod)
{
	// Open the file
	FILE* fp = fopen(filename, "rb");
	if (!fp) {
		printf("Unable to open file: %s\n", filename);

	}

	unsigned char buffer[BUFFER_SIZE];
	size_t bytes_read = 0;
	size_t total_bytes_read = 0;
	while ((bytes_read = fread(buffer + total_bytes_read, 1, BUFFER_SIZE - total_bytes_read, fp))) {
		total_bytes_read += bytes_read;
		if (total_bytes_read == BUFFER_SIZE) {
			break;
		}
	}



	string buffer_string(reinterpret_cast<char*>(buffer), BUFFER_SIZE);
	SHA512 sha512;
	string hash_str = sha512.hash(buffer_string);

//	std::cout << "Message: " << hash_str << std::endl;


	digest = create_vectors2(hash_str, oil, mod);
	

	fclose(fp);


}

