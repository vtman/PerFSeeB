#include <fstream>

#include <chrono>

#include <smmintrin.h>
#include <emmintrin.h>

#include <omp.h>

class SeedCheck{
public:
	SeedCheck(int nLen_in, int nMismatch_in, int *iv, int *ov);
	~SeedCheck();

	void checkCandidate(int istart);

	__m128i* uSeed, uRes[128], mOne[128], mOO, mON;

	int* inVec, * outVec, *ip;
	int nLen, none, nMismatch, nLQ;
	int nFound;

	unsigned int u[4];

};

inline long long areSame(__m128i m1, __m128i m2) {
	__m128i m3, m4;

	m3 = _mm_xor_si128(m1, m2);
	m4 = _mm_bsrli_si128(m3, 8);
	m3 = _mm_or_si128(m3, m4);
	return _mm_extract_epi64(m3, 0);
}


SeedCheck::SeedCheck(int nLen_in, int nMismatch_in, int *iv, int *ov) {
	nLen = nLen_in;
	nMismatch = nMismatch_in;

	nFound = 0;

	inVec = iv;
	outVec = ov;

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
	for (int i = 0; i < nLen; i++) {
		mON = _mm_xor_si128(mON, mOne[i]);
	}

	uSeed = (__m128i*)malloc(sizeof(__m128i) * nLen);
	ip = (int*)malloc(sizeof(int) * nLen);

	nLQ = nLen - nMismatch + 1;
}

SeedCheck::~SeedCheck() {
	if (uSeed != nullptr) {
		free(uSeed); uSeed = nullptr;
	}
	if (ip != nullptr) {
		free(ip); ip = nullptr;
	}
}

void SeedCheck::checkCandidate(int istart) {
	int *ipos;
	int ik, ipnq, nq;
	bool isOk;

	ipos = inVec + none * istart;

	
	
	for (int k = 0; k < nLen; k++) {
		uSeed[k] = mON;
		for (int i = 0; i < none; i++) {
			ik = (ipos[i] + k) % nLen;
			uSeed[k] = _mm_or_si128(uSeed[k], mOne[ik]);
		}
	}

	isOk = true;

	for (int jj = 0; jj < nMismatch; jj++) {
		ip[jj] = jj;
	}

	uRes[0] = mON;

	for (nq = 0; nq < nMismatch; ++nq) {
		uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ip[nq]]);
		if (areSame(uRes[nq + 1], mOO) == 0) {
			isOk = false;
			break;
		}
		if (areSame(uRes[nq + 1], uRes[nq]) == 0) break;
		if (nq == nMismatch - 1)break;
	}

	if (isOk) {
		ipnq = ip[nq];
		do {
			do {
				ipnq++;
				if (ipnq - nq == nLQ) {
					nq--;
					if (nq == 0) break;
					ipnq = ip[nq];
				}
				else {
					break;
				}
			} while (true);
			if (nq == 0) break;

			uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ipnq]);
			ip[nq] = ipnq;
			if (areSame(uRes[nq + 1], uRes[nq]) == 0) continue;

			while (nq < nMismatch - 1) {
				nq++;
				ipnq++;
				ip[nq] = ipnq;
				uRes[nq + 1] = _mm_or_si128(uRes[nq], uSeed[ipnq]);
				if (areSame(uRes[nq + 1], uRes[nq]) == 0)break;
			}

			if (areSame(uRes[nq + 1], mOO) == 0) {
				isOk = false;
				break;
			}
			if (!isOk) break;

		} while (true);
	}

	if (isOk) {
		memcpy(outVec + nFound * none, ipos, none*sizeof(int));
		nFound++;
	}

}


void printP(FILE *fo, int* ipos, int none) {
	int m;
	m = 0;

	for (int i = 0; i < ipos[0]; i++) {
		fprintf(fo, "0");
	}
	fprintf(fo, "1");

	for (int i = 1; i < none; i++) {
		m = ipos[i] - ipos[i - 1] - 1;
		for (int j = 0; j < m; j++) {
			fprintf(fo, "0");
		}
		fprintf(fo, "1");
	}
	fprintf(fo, "\n");
}


void string2Pos(const char* s, int* ipos) {
	int n, k;
	n = strlen(s);
	k = 0;
	for (int i = 0; i < n; i++) {
		if (s[i] == '1') {
			ipos[k] = i;
			k++;
		}
	}
}


int main2(char *outputFolder, int nMismatch, int nLen, int nstart, int nthreads, int nchunk, long long startCount, char *iniString) {
	char fName[1000];

	int* nSolFound, ** pindSol, ** pcandSol;
	int nL, nTotal, nth;
	int * ipos, * gapSize, * indSame;
	int maxOnes, att, maxGap, mgg;
	int nCut1, mgg2, nSame, nn, ndiff, ki, max_threads;

	long long countTest;

	bool tEqual, tLess, isFound, t, isFileClosed;

	FILE* fo;

#pragma omp parallel
	{
		max_threads = omp_get_num_threads();
	}

	printf("Maximum number of threads: %i\n", max_threads);

	nth = max_threads;
	if (nthreads < nth) {
		nth = nthreads;
	}
	if (nth < 1) nth = 1;

	SeedCheck** sc;

	sc = (SeedCheck**)malloc(sizeof(SeedCheck*) * nth);

	nSolFound = (int*)malloc(sizeof(int) * nth);
	pindSol = (int**)malloc(sizeof(int*) * nth);
	pcandSol = (int**)malloc(sizeof(int*) * nth);

	for (int i = 0; i < nth; i++) {
		pindSol[i] = (int*)malloc(sizeof(int) * nLen * nchunk);
		pcandSol[i] = (int*)malloc(sizeof(int) * nLen * nchunk);
	}

	for (int i = 0; i < nth; i++) {
		sc[i] = new SeedCheck(nLen, nMismatch, pcandSol[i], pindSol[i]);
	}

	printf("Number of threads we use: %i\n", nth);

	maxOnes = nLen / nMismatch;
	if (nLen % nMismatch > 0) maxOnes++;

	att = 0;

	ipos = (int*)malloc(sizeof(int) * nLen);
	gapSize = (int*)malloc(sizeof(int) * nLen);
	indSame = (int*)malloc(sizeof(int) * nLen);

	nTotal = 0;

	isFound = false;

	isFileClosed = true;
	fo = nullptr;

		
	int nCount;

	omp_set_num_threads(nth);

	for (int none = nstart; none >= 1; none--) {
		printf("Mismatches: %i, len: %i, level: %i (nOne: %i, nZero: %i)\n", nMismatch, nLen, att, none, nLen - none);

		countTest = 0;

		nCount = 0;

		tEqual = false;
		tLess = false;

		nL = 0;
		ipos[nL] = nMismatch-1;
		ipos[none - 1] = nLen - 1;
		nCut1 = nLen - none + 1;
		nn = (none - 1) / 2;

		if (iniString[0] != 'x') {
			string2Pos(iniString, ipos);
			nL = none - 2;
		}

		do {
			nCount = 0;

			do {

				do {
					ipos[nL]++;
					if (ipos[nL] - nL == nCut1) {
						nL--;
						if (nL < 0) break;
					}
					else {
						break;
					}
				} while (true);

				if (nL < 0) break;
				if (nL == 0) {
					maxGap = ipos[0] + 1;
					tEqual = false;
				}
				else {
					mgg = ipos[nL] - ipos[nL - 1];
					if (mgg > maxGap) {
						nL--;
						tEqual = false;
						tLess = true;
						continue;
					}
					if (mgg == maxGap) {
						tEqual = true;
					}
				}

				while (nL < none - 2) {
					nL++;
					ipos[nL] = ipos[nL - 1] + 1;
				}

				t = false;
				for (int k = 0; k < nn; k++) {
					ndiff = (ipos[none - 1 - k] - ipos[none - 2 - k]) - (ipos[k + 1] - ipos[k]);
					if (ndiff == 0) continue;
					if (ndiff > 0) t = true;
					break;
				}

				if (t) continue;

				mgg2 = ipos[none - 1] - ipos[none - 2];

				if (mgg2 > maxGap) {
					tEqual = false;
					continue;
				}

				if (mgg2 == maxGap) {
					tEqual = true;
				}

				if (tEqual || tLess) {
					gapSize[0] = ipos[0] + 1;
					maxGap = gapSize[0];
					nSame = 1;
					indSame[0] = 0;
					t = false;

					for (int k = 1; k < none; k++) {
						gapSize[k] = ipos[k] - ipos[k - 1];
						if (gapSize[k] == maxGap) {
							indSame[nSame] = k;
							nSame++;
						}
					}
					if (nSame == 1) {
						tLess = false;
						tEqual = false;
					}
					else {
						for (int k = 1; k < nSame; k++) {
							ki = indSame[k];
							for (int i = 1; i < none; i++) {
								if (gapSize[(ki + i) % none] == gapSize[i])continue;
								if (gapSize[(ki + i) % none] > gapSize[i]) t = true;
								break;
							}
							if (t) break;
						}


						if (!t) {
							for (int k = 0; k < nSame; k++) {
								ki = indSame[k];
								for (int i = 0; i < none; i++) {
									if (gapSize[(ki - i + none) % none] == gapSize[i])continue;
									if (gapSize[(ki - i + none) % none] > gapSize[i]) t = true;
									break;
								}
								if (t) break;
							}
						}
					}

					if (t) continue;

				}

				for (int i = 0; i < none; i++) {
					pcandSol[0][nCount * none + i] = ipos[i];
				}

				nCount++;
				countTest++;

			} while (nCount < nchunk);

			printf("%lld\t%i\n", countTest, nCount);

			if (countTest > startCount) {

				for (int i = 1; i < nth; i++) {
					memcpy(pcandSol[i], pcandSol[0], sizeof(int) * nCount * none);
				}
				for (int i = 0; i < nth; i++) {
					sc[i]->none = none;
					sc[i]->nFound = 0;
				}

#pragma omp parallel
				{

					int tid, i;
					SeedCheck* sc_loc;
					tid = omp_get_thread_num();
					sc_loc = sc[tid];

#pragma omp for
					for (i = 0; i < nCount; i++) {
						sc_loc->checkCandidate(i);
					}
				}

				for (int m = 0; m < nth; m++) {
					if (sc[m]->nFound == 0) continue;

					for (int k = 0; k < sc[m]->nFound; k++) {
						isFound = true;
						if (isFileClosed) {
							sprintf(fName, "%s/block_m_%i_len_%i_level_%i.txt", outputFolder, nMismatch, nLen, att);
							fo = fopen(fName, "w");
							isFileClosed = false;
						}
						printP(fo, sc[m]->outVec + k * none, none);
					}
				}
			}

		}while (nL >= 0);
			
		if (!isFound) att = -1;

		att++;
		if (fo != nullptr) {
			fclose(fo); fo = nullptr;
			isFileClosed = true;
		}
		if (att == 1)break;
		if (nL < 0) continue;

	}

	printf("Num: %lld\n", countTest);
		
	free(ipos); ipos = nullptr;
	
	if (fo != nullptr) {
		fclose(fo); fo = nullptr;
	}

	free(nSolFound); nSolFound = nullptr;
	

	for (int i = 0; i < nth; i++) {
		free(pindSol[i]);
		free(pcandSol[i]);
		delete sc[i];
	}
	free(pindSol); pindSol = nullptr;
	free(pcandSol); pcandSol = nullptr;
	free(sc); sc = nullptr;

	return 0;
}

int main(int argc, char* argv[]) {
	int nMismatch, nLen, nstart, nchunk, nthreads;
	char iniString[1000], outputFolder[1000];
	long long startCount;

	printf("Number of arguments: %i\n", argc);

	sprintf(outputFolder, argv[1]);
	nMismatch = atoi(argv[2]);
	nLen = atoi(argv[3]);
	nstart = atoi(argv[4]);
	nchunk = atoi(argv[5]);
	nthreads = atoi(argv[6]);

	startCount = -1;

	if (argc > 7) {
		startCount = atoll(argv[7]);
	}

	sprintf(iniString, "x");
	if (argc > 8) {
		sprintf(iniString, argv[8]);
	}
	
	printf("Output folder: \"%s\"\n", outputFolder);
	printf("Number of mismatches: %i\n", nMismatch);
	printf("Length of the block: %i\n", nLen);
	printf("Initial number of 0-elements: %i\n", nstart);
	printf("Size of chunks for parallelisation: %i\n", nchunk);
	printf("Number of threads: %i\n", nthreads);
	printf("Checking after N candidates, N: %lld\n", startCount);
	printf("Initial string: %s\n", iniString);
	
	auto start = std::chrono::steady_clock::now();

	main2(outputFolder, nMismatch, nLen, nLen - nstart, nthreads, nchunk, startCount, iniString);
	auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("It took me %lld ms.\n", elapsed.count());

	return 0;
}

