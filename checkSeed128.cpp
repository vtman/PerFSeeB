#include <stdio.h>
#include <iostream>

#include <xmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>

inline long long areSame(__m128i *m1, __m128i *m2, int nR) {
	__m128i m3, m4;
	for (int i = 0; i < nR; i++) {
		m4 = _mm_xor_si128(m1[i], m2[i]);
		if (i == 0) {
			m3 = m4;
		}
		else {
			m3 = _mm_or_si128(m3, m4);
		}
	}
	m4 = _mm_bsrli_si128(m3, 8);
	m3 = _mm_or_si128(m3, m4);
	return _mm_extract_epi64(m3, 0);
}

int print128(__m128i mu, int n) {
	unsigned int u[4];
	u[0] = _mm_extract_epi32(mu, 0);
	u[1] = _mm_extract_epi32(mu, 1);
	u[2] = _mm_extract_epi32(mu, 2);
	u[3] = _mm_extract_epi32(mu, 3);

	for (int i = 0; i < n; i++) {
		if (((u[i / 32] >> (i % 32)) & 1) == 1) {
			printf("1");
		}
		else {
			printf("0");
		}
	}
	printf("\n");

	return 0;
}

int main(int argc, char* argv[]) {
	char cSeed[1000], outLine[1000];
	int ik, k1, k2, isOut[1000], indOld[1000];
	int nOut, nq, ipnq, nLQ, nR, ind;
	int readLen, seedLen, nMismatches, nRows;
	int* ip;

	unsigned int u[4];
	__m128i* uSeed, *uRes, mOne[128], mOO, mON, mZ, * pmOO;
	bool isOk;

	if (argc != 4) {
		printf("Usage:\n");
		printf("\t1) Read length (integer)\n");
		printf("\t2) Number of mismatches (integer)\n");
		printf("\t3) Seed (string of 1 or 0)\n");
		printf("\n");
		return -1;
	}

	readLen = atoi(argv[1]);
	nMismatches = atoi(argv[2]);
	sprintf(cSeed, argv[3]);
	printf("Seed: %s\n", cSeed);
	seedLen = strlen(cSeed);
	printf("Seed length: %i\n", seedLen);
	if (seedLen < 1) {
		printf("Wrong seed length\n");
		return -3;
	}
	printf("Read length: %i\n", readLen);
	if (readLen < 2) {
		printf("Wrong read length\n");
		return -4;
	}

	nRows = readLen - seedLen + 1;

	if (nRows < 1) {
		printf("Number of rows should be positive (%i)\n", nRows);
		return -2;
	}

	printf("Number of mismatches: %i\n", nMismatches);
	if (nMismatches < 1 || nMismatches > readLen) {
		printf("Wrong number of mismatches\n");
		return -5;
	}
	
	printf("\n");

	nR = nRows / 128;
	if (nRows % 128 > 0) nR++;

	mZ = _mm_set1_epi32(0x00000000);
	mOO = _mm_set1_epi32(0xffffffff);

	for (int i = 0; i < 128; i++) {
		u[0] = 0;
		u[1] = 0;
		u[2] = 0;
		u[3] = 0;

		u[i / 32] = 1 << (i % 32);
		mOne[i] = _mm_set_epi32(u[3], u[2], u[1], u[0]);
	}

	mON = mOO;
	for (int i = 0; i < nRows%128; i++) {
		mON = _mm_xor_si128(mON, mOne[i]);
	}

	pmOO = (__m128i*)malloc(sizeof(__m128i) * nR);
	uSeed = (__m128i*)malloc(sizeof(__m128i) * readLen*nR);
	uRes = (__m128i*)malloc(sizeof(__m128i) * (nMismatches + 1) * nR);
	ip = (int*)malloc(sizeof(int) * nMismatches);

	for (int i = 0; i < nR; i++) {
		pmOO[i] = mOO;
	}

	for (int i = 0; i < readLen; i++) {
		for (int k = 0; k < nR - 1; k++) {
			uSeed[i * nR + k] = mZ;
		}
		uSeed[i*nR + (nR-1)] = mON;
		isOut[i] = 0;
	}

	for (int i = 0; i < seedLen; i++) {
		if (cSeed[i] == '0') continue;
		for (int k = 0; k < nRows; k++) {
			ik = i + k;
			k1 = k / 128;
			k2 = k % 128;
			uSeed[ik*nR + k1] = _mm_or_si128(uSeed[ik*nR + k1], mOne[k2]);
		}
	}

	for (int i = 0; i < readLen; i++) {
		if (isOut[i] == 1) continue;
		for (int j = 0; j < readLen; j++) {
			if (isOut[j] == 1) continue;
			if (i == j)continue;
			
			for (int k = 0; k < nR; k++) {
				uRes[k] = _mm_or_si128(uSeed[i * nR + k], uSeed[j * nR + k]);
			}
			if (areSame(uSeed + j * nR, uRes, nR) == 0) {
				isOut[i] = 1;
				break;
			}
		}
	}
		
	nOut = 0;

	for (int i = 0; i < readLen; i++) {
		if (isOut[i] == 1) continue;
		for (int k = 0; k < nR; k++) {
			uSeed[nOut * nR + k] = uSeed[i * nR + k];
		}
		indOld[nOut] = i;
		nOut++;
	}

	nLQ = nOut - nMismatches + 1;

	printf("Excluded: %i\n\n", readLen - nOut);

	outLine[readLen] = '\0';

	for (int i = 0; i < nRows; i++) {
		for (int j = 0; j < readLen; j++) {
			if (j < i) {
				outLine[j] = '0';
			}
			else {
				if (j - i >= seedLen) {
					outLine[j] = '0';
				}
				else {
					outLine[j] = cSeed[j - i];
				}
			}
		}
		printf("%s\n", outLine);
	}

	isOk = true;

	for (int jj = 0; jj < nMismatches; jj++) {
		ip[jj] = jj;
	}

	for (int i = 0; i < nR - 1; i++) {
		uRes[i] = mZ;
	}
	uRes[nR-1] = mON;
	
	for (nq = 0; nq < nMismatches; ++nq) {
		for (int ir = 0; ir < nR; ir++) {
			uRes[(nq + 1)*nR + ir] = _mm_or_si128(uRes[nq*nR + ir], uSeed[ip[nq]*nR + ir]);
		}
		if (areSame(uRes + (nq + 1)*nR, pmOO, nR) == 0) {
			isOk = false;
			break;
		}
		if (areSame(uRes +(nq + 1)*nR, uRes + nq*nR, nR) == 0) break;
		if (nq == nMismatches - 1)break;
	}

	if (isOk) {
		ipnq = ip[nq];
		do {
			do {
				ipnq++;
				if (ipnq - nq == nLQ) {
					nq--;
					if (nq == -1) break;
					ipnq = ip[nq];
				}
				else {
					break;
				}
			} while (true);
			if (nq == -1) break;

			for (int ir = 0; ir < nR; ir++) {
				uRes[(nq + 1)*nR + ir] = _mm_or_si128(uRes[nq*nR + ir], uSeed[ipnq*nR + ir]);
			}
			ip[nq] = ipnq;
			if (areSame(uRes + (nq + 1)*nR, uRes + nq*nR, nR) == 0) continue;

			while (nq < nMismatches - 1) {
				nq++;
				ipnq++;
				ip[nq] = ipnq;
				for (int ir = 0; ir < nR; ir++) {
					uRes[(nq + 1)*nR + ir] = _mm_or_si128(uRes[nq*nR + ir], uSeed[ipnq*nR + ir]);
				}
				if (areSame(uRes + (nq + 1)*nR, uRes + nq*nR, nR) == 0)break;
			}
			//printf("%i\t%i\n", ip[0], ip[1]);
			if (areSame(uRes +(nq + 1)*nR, pmOO, nR) == 0) {
				isOk = false;
				break;
			}
			if (!isOk) break;

		} while (true);
	}

	if (isOk) {
		printf("The seed can be used.\n");
	}
	else {
		ind = 0;
		for (int i = 0; i < nMismatches; i++) {
			while (indOld[ip[i]] > ind) {
				outLine[ind] = '_';
				//printf("_");
				ind++;
			}
			outLine[ind] = 'X';
			//printf("X");
			ind++;
		}
		while (ind < readLen) {
			outLine[ind] = '_';
			//printf("_");
			ind++;
		}
		printf("%s\n", outLine);
		for (int i = 0; i < nMismatches; i++) {
			printf("%i\t", indOld[ip[i]] + 1);
		}
		printf("\n");
		printf("The seed cannot be used.\n");
	}
	
	free(uRes); uRes = nullptr;
	free(uSeed); uSeed = nullptr;
	free(ip); ip = nullptr;
	return 0;
}
