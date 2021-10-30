#include <fstream>
#include <chrono>

#include <smmintrin.h>
#include <emmintrin.h>

#define storLen 1000000

class BComb {
public:
	BComb();
	~BComb();

	int nItems[128];
	int startProcessing(char* outputFolder, int inMismatches, int inReadLen, int iminWeight);
	int initiatePatterns();
	int findPatterns(int iLevel);
	bool checkLeftSubstring(__m128i mCand);
	bool checkM(__m128i m, int nSLen);
	int print128(__m128i x, int n);

	__m128i* mtemp, * mVec, * mVecOrig;
	__m128i* mseeds[128];

	int* countOne;
	int nSolutions, nMismatches, nReadLen, minWeight;
	int curLevel, nFound;
	int* indSortVec, * nSort, * ipos;
	int nCounter;
	bool isOut[storLen];
	FILE* fo;
};


int BComb::print128(__m128i x, int n) {
	unsigned int a, b[4], uu;
	int istart;
	bool isO;

	b[0] = _mm_extract_epi32(x, 0);
	b[1] = _mm_extract_epi32(x, 1);
	b[2] = _mm_extract_epi32(x, 2);
	b[3] = _mm_extract_epi32(x, 3);

	istart = 0;
	if (b[0] == 0) {
		istart++;
		if (b[1] == 0) {
			istart++;
			if (b[2] == 0) {
				istart++;
			}
		}
	}

	uu = 1 << 31;

	isO = false;

	for (int i = istart; i < 4; i++) {
		a = b[i];
		for (int k = 0; k < 32; k++) {
			if ((a & 1) == 1) {
				fprintf(fo, "1");
				printf("1");
				isO = true;
			}
			else {
				if (isO) {
					fprintf(fo, "0");
					printf("0");
				}
			}
			a = (a >> 1);
		}
	}
	fprintf(fo, "\t%i\n", n);
	printf("\t%i\n", n);
	return 0;
}


int print128_2(__m128i x) {
	unsigned int a, b[4], uu;
	int istart;
	bool isO;

	b[0] = _mm_extract_epi32(x, 0);
	b[1] = _mm_extract_epi32(x, 1);
	b[2] = _mm_extract_epi32(x, 2);
	b[3] = _mm_extract_epi32(x, 3);

	istart = 0;
	if (b[0] == 0) {
		istart++;
		if (b[1] == 0) {
			istart++;
			if (b[2] == 0) {
				istart++;
			}
		}
	}

	uu = 1 << 31;

	isO = false;

	for (int i = istart; i < 4; i++) {
		a = b[i];
		for (int k = 0; k < 32; k++) {
			if ((a & 1) == 1) {
				printf("1");
				isO = true;
			}
			else {
				if (isO) {
					printf("0");
				}
			}
			a = (a >> 1);
		}
	}

	printf("\n");
	return 0;
}

int shiftOneBitLeft(__m128i m1, __m128i* m2) {
	__m128i a1, a2, a3;
	a1 = _mm_slli_epi32(m1, 1);
	a2 = _mm_srli_epi32(m1, 31);
	a3 = _mm_bslli_si128(a2, 4);
	*m2 = _mm_or_si128(a1, a3);
	return 0;
}

int shiftOneBitRight(__m128i m1, __m128i* m2) {
	__m128i a1, a2, a3;
	a1 = _mm_srli_epi32(m1, 1);
	a2 = _mm_slli_epi32(m1, 31);
	a3 = _mm_bsrli_si128(a2, 4);
	*m2 = _mm_or_si128(a1, a3);
	return 0;
}


BComb::BComb() {
	countOne = nullptr;
	fo = nullptr;
	nCounter = 0;
	mtemp = nullptr;
	mVec = nullptr;
	mVecOrig = nullptr;
	ipos = nullptr;

	indSortVec = nullptr;
	nSort = nullptr;

	for (int i = 0; i < 128; i++) {
		mseeds[i] = nullptr;
	}
}

BComb::~BComb() {
	if (countOne != nullptr) {
		free(countOne); countOne = nullptr;
	}
	if (fo != nullptr) {
		fclose(fo); fo = nullptr;
	}
	if (ipos != nullptr) {
		free(ipos); ipos = nullptr;
	}
	if (indSortVec != nullptr) {
		free(indSortVec); indSortVec = nullptr;
	}
	if (nSort != nullptr) {
		free(nSort); nSort = nullptr;
	}
	if (mVec != nullptr) {
		free(mVec); mVec = nullptr;
	}
	if (mVecOrig != nullptr) {
		free(mVecOrig); mVecOrig = nullptr;
	}

	if (mtemp != nullptr) {
		free(mtemp); mtemp = nullptr;
	}
	for (int i = 0; i < 128; i++) {
		if (mseeds[i] != nullptr) {
			free(mseeds[i]); mseeds[i] = nullptr;
		}
	}

}


int BComb::initiatePatterns() {
	mseeds[0] = (__m128i*)malloc(sizeof(__m128i) * 1);
	mseeds[0][0] = _mm_set_epi32(0x80000000, 0, 0, 0);
	curLevel = 1;
	nItems[0] = 1;

	nSolutions = 1;

	return 0;
}

inline long long areSame(__m128i m1, __m128i m2) {
	__m128i m3, m4;

	m3 = _mm_xor_si128(m1, m2);
	m4 = _mm_bsrli_si128(m3, 8);
	m3 = _mm_or_si128(m3, m4);
	return _mm_extract_epi64(m3, 0);
}

__m128i shiftRight(__m128i m, int n) {
	__m128i a1, a2, a3, mRes;
	int ishift, ns;
	ns = n / 32;
	ishift = n % 32;
	a1 = _mm_srli_epi32(m, ishift);
	a2 = _mm_slli_epi32(m, 32 - ishift);
	a3 = _mm_bsrli_si128(a2, 4);
	mRes = _mm_or_si128(a1, a3);

	if (ns == 3) {
		mRes = _mm_bsrli_si128(mRes, 12);
	}
	if (ns == 2) {
		mRes = _mm_bsrli_si128(mRes, 8);
	}
	if (ns == 1) {
		mRes = _mm_bsrli_si128(mRes, 4);
	}
	return mRes;
}

bool BComb::checkLeftSubstring(__m128i mCand) {
	int nShift, ns, ishift, n1, n2, n3, n;
	unsigned long long ur1, ur0, * pv;
	unsigned int u[4], ua;
	__m128i* mp;

	u[0] = _mm_extract_epi32(mCand, 0);
	ns = 0;
	if (u[0] == 0) {
		u[1] = _mm_extract_epi32(mCand, 1);
		ns = 1;
		if (u[1] == 0) {
			u[2] = _mm_extract_epi32(mCand, 2);
			ns = 2;
			if (u[2] == 0) {
				u[3] = _mm_extract_epi32(mCand, 3);
				ns = 3;
			}
		}
	}
	ua = u[ns];

	ishift = 0;

	if ((ua & 0x0000ffff) == 0) {
		ishift += 16;
		ua >>= 16;
	}
	if ((ua & 0x000000ff) == 0) {
		ishift += 8;
		ua >>= 8;
	}
	if ((ua & 0x0000000f) == 0) {
		ishift += 4;
		ua >>= 4;
	}
	if ((ua & 0x00000003) == 0) {
		ishift += 2;
		ua >>= 2;
	}
	if ((ua & 0x00000001) == 0) {
		ishift += 1;
		ua >>= 1;
	}

	nShift = 127 - (32 * ns + ishift);

	ur0 = _mm_extract_epi64(mCand, 0);
	ur1 = _mm_extract_epi64(mCand, 1);

	mp = mseeds[nShift];
	n = nItems[nShift];

	pv = (unsigned long long*)mp;

	n1 = 0;
	n2 = n - 1;

	if (ur1 < pv[1]) return false;
	if (ur1 > pv[2 * n2 + 1])return false;
	if (ur1 == pv[1]) {
		if (ur0 < pv[0]) return false;
		if (ur0 == pv[0]) return true;
	}
	if (ur1 == pv[2 * n2 + 1]) {
		if (ur0 > pv[2 * n2]) return false;
		if (ur0 == pv[2 * n2]) return true;
	}

	n3 = (n1 + n2) / 2;
	do {
		if (n2 <= n1 + 1) {
			if (ur1 == pv[2 * n3 + 1] && ur0 == pv[2 * n3]) {
				return true;
			}
			return false;
		}
		n3 = (n1 + n2) / 2;
		if (ur1 < pv[2 * n3 + 1]) {
			n2 = n3;
		}
		else if (ur1 > pv[2 * n3 + 1]) {
			n1 = n3;
		}
		else {
			if (ur0 > pv[2 * n3]) {
				n1 = n3;
			}
			else if (ur0 < pv[2 * n3]) {
				n2 = n3;
			}
			else {
				return true;
			}
		}

	} while (true);

	return false;
}


int BComb::findPatterns(int iLevel) {
	__m128i mOne, mCand, mOneLeft, mLeft;
	unsigned int u[4];
	int i1, i2, nTotal, nOnes;
	bool isFound;

	mtemp = (__m128i*)malloc(sizeof(__m128i) * nSolutions);

	nTotal = 0;

	if (mtemp == nullptr) {
		printf("Error: cannot allocate memory\n");
		return -1;
	}

	u[0] = 0;
	u[1] = 0;
	u[2] = 0;
	u[3] = 0;

	i1 = (128 - iLevel) / 32;
	i2 = (128 - iLevel) % 32;

	u[i1] = 1 << i2;

	mOneLeft = _mm_set_epi32(u[3], u[2], u[1], u[0]);
	mOne = _mm_set_epi32(0x80000000, 0, 0, 0);

	for (int ip = 0; ip < iLevel - 1; ip++) {
		for (int k = 0; k < nItems[ip]; k++) {
			mCand = shiftRight(mseeds[ip][k], iLevel - ip - 1);
			mCand = _mm_or_si128(mOne, mCand);
			if (ip == 0) {
				isFound = true;
			}
			else {
				mLeft = _mm_xor_si128(mCand, mOneLeft);
				isFound = checkLeftSubstring(mLeft);
				if (!isFound) {
					continue;
				}
			}
			nCounter++;
			if (!checkM(mCand, iLevel)) {
				continue;
			}
			nOnes = __popcnt64(_mm_extract_epi64(mCand, 0)) + __popcnt64(_mm_extract_epi64(mCand, 1));
			if (nOnes >= minWeight) {
				print128(mCand, nOnes);
			}
			mtemp[nTotal] = mCand;
			nTotal++;
		}
	}

	mseeds[iLevel - 1] = (__m128i*)malloc(sizeof(__m128i) * nTotal);
	memcpy(mseeds[iLevel - 1], mtemp, nTotal * sizeof(__m128i));
	nItems[iLevel - 1] = nTotal;
	free(mtemp); mtemp = nullptr;

	nSolutions += nTotal;

	return 0;
}

int reduceVector(__m128i* vm, int* nL, bool* isOut) {
	int n, nnew;
	__m128i mq;
	n = *nL;
	for (int i = 0; i < n; i++) {
		isOut[i] = false;
	}

	for (int i = 0; i < n; i++) {
		if (isOut[i]) continue;
		for (int j = i + 1; j < n; j++) {
			if (isOut[j]) continue;

			mq = _mm_or_si128(vm[i], vm[j]);
			if (areSame(mq, vm[j]) == 0) {
				isOut[i] = true;
				break;
			}
			else {
				if (areSame(mq, vm[i]) == 0) {
					isOut[j] = true;
				}
			}
		}
	}

	nnew = 0;

	for (int i = 0; i < n; i++) {
		if (isOut[i]) continue;
		vm[nnew] = vm[i];
		nnew++;
	}
	*nL = nnew;

	return 0;
}

bool BComb::checkM(__m128i m, int nSLen) {
	char u[9], b[16];
	int maskLen, curLen, nSmax, nn, ic, nL;
	__m128i mq, m1, m2, mMask;

	u[0] = 0;
	u[1] = 128;
	u[2] = 192;
	u[3] = 224;
	u[4] = 240;
	u[5] = 248;
	u[6] = 252;
	u[7] = 254;
	u[8] = 255;

	maskLen = nReadLen - nSLen + 1;
	for (int i = 0; i < maskLen / 8; i++) {
		b[15 - i] = u[8];
	}
	if (maskLen % 8 == 0) {
		for (int i = maskLen / 8; i < 16; i++) {
			b[15 - i] = 0;
		}
	}
	else {
		b[15 - maskLen / 8] = u[maskLen % 8];
		for (int i = maskLen / 8 + 1; i < 16; i++) {
			b[15 - i] = 0;
		}
	}

	mMask = _mm_set_epi8(b[15], b[14], b[13], b[12], b[11], b[10], b[9], b[8], b[7], b[6], b[5], b[4], b[3], b[2], b[1], b[0]);

	m1 = m;

	for (int i = 0; i < maskLen - 1; i++) {
		shiftOneBitRight(m1, &m2);
		m1 = m2;
		mVecOrig[maskLen - 2 - i] = _mm_and_si128(m1, mMask);
	}

	m1 = m;

	for (int i = 0; i < nSLen; i++) {
		mVecOrig[maskLen - 1 + i] = _mm_and_si128(m1, mMask);
		shiftOneBitLeft(m1, &m2);
		m1 = m2;
	}

	curLen = nReadLen;

	reduceVector(mVecOrig, &curLen, isOut);

	nSmax = -1;

	for (int i = 0; i <= maskLen; i++) {
		nSort[i] = 0;
	}

	for (int i = 0; i < curLen; i++) {
		nn = __popcnt64(_mm_extract_epi64(mVecOrig[i], 0)) + __popcnt64(_mm_extract_epi64(mVecOrig[i], 1));
		indSortVec[128 * nn + nSort[nn]] = i;
		nSort[nn]++;
		if (nn > nSmax) nSmax = nn;
	}

	if (nSmax == maskLen) return false;

	if (nSmax* nMismatches < maskLen) return true;

	ic = 0;
	for (int i = nSmax; i > 0; i--) {
		for (int j = 0; j < nSort[i]; j++) {
			mVec[ic] = mVecOrig[indSortVec[128 * i + j]];
			countOne[ic] = i;
			ic++;
		}
	}

	for (int i = 0; i < nMismatches; i++) {
		ipos[i] = i;
	}

	do {
		mq = mVec[ipos[0]];
		for (int i = 1; i < nMismatches; i++) {
			mq = _mm_or_si128(mq, mVec[ipos[i]]);
		}
		if (areSame(mq, mMask) == 0) {
			return false;
		}

		nL = nMismatches - 1;
		do {
			ipos[nL]++;
			if (ipos[nL] == ic - nMismatches + 1 + nL) {
				nL--;
				if (nL < 0) break;
			}
			else {
				break;
			}
		} while (true);
		if (nL < 0) break;

		while (nL < nMismatches - 1) {
			nL++;
			ipos[nL] = ipos[nL - 1] + 1;
		}

	} while (true);

	return true;
}

int BComb::startProcessing(char* outputFolder, int inMismatches, int inReadLen, int iminWeight) {
	char ofileName[1000];
	nMismatches = inMismatches;
	nReadLen = inReadLen;
	minWeight = iminWeight;

	sprintf(ofileName, "%s/seed_%i_%i_%i.txt", outputFolder, nMismatches, nReadLen, minWeight);
	fo = fopen(ofileName, "w");

	printf("Length of reads: %i\n", nReadLen);
	printf("Number of mimatches: %i\n", nMismatches);
	printf("Minimum weight: %i\n", minWeight);

	fprintf(fo, "Length of reads: %i\n", nReadLen);
	fprintf(fo, "Number of mismatches: %i\n", nMismatches);
	fprintf(fo, "Minimum weight: %i\n", minWeight);

	ipos = (int*)malloc(sizeof(int) * nMismatches);
	nSort = (int*)malloc(sizeof(int) * 128);
	indSortVec = (int*)malloc(sizeof(int) * 129 * 128);
	mVec = (__m128i*)malloc(sizeof(__m128i) * nReadLen);
	countOne = (int*)malloc(sizeof(int) * nReadLen);
	mVecOrig = (__m128i*)malloc(sizeof(__m128i) * nReadLen);

	initiatePatterns();

	for (int i = 2; i <= nReadLen - nMismatches; i++) {
		printf("Length: %i\n", i);
		findPatterns(i);
		printf("Solutions: %i\n\n", nSolutions);
	}

	printf("Total new: %i\n", nSolutions);

	return 0;
}

int main(int argc, char* argv[]) {
	char outputFolder[1000];
	int ires, nMismatches, nReadLen, minWeight;
	BComb* bc;

	if (argc != 5) {
		printf("Usage:");
		printf("\t 1) Output folder\n");
		printf("\t 2) Number of mismatches\n");
		printf("\t 3) Length of reads\n");
		printf("\t 4) Minimum weight\n\n");
		return -999;
	}

	sprintf(outputFolder, argv[1]);
	nMismatches = atoi(argv[2]);
	nReadLen = atoi(argv[3]);
	minWeight = atoi(argv[4]);

	auto start = std::chrono::steady_clock::now();

	bc = new BComb();
	ires = bc->startProcessing(outputFolder, nMismatches, nReadLen, minWeight);
	delete bc; bc = nullptr;

	auto end = std::chrono::steady_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("It took me %lld ms.\n", elapsed.count());

	return ires;
}