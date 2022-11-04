#include <fstream>
#include <stdlib.h>
#include <string.h>

#include "ipp.h"

//D:\Gen2022\seed1 111111111100011101100100100111010011100010100101000010100110000101111000000011

//D:\BMC\temp 111011010100110101100100100101010110001

//D:\BMC\temp 111110111001011111011100101111101110010111110111001011111011100101111101110010111110111

class blockData {
public:
	blockData(char* ioFolder_in, int nL_in, int nW_in);
	~blockData();

	int startProcessing(int indexFile);
	int setFileNames();
	int freeMemory();

	bool is32;

	int ndigits, nW, nL, nE, iprefix, sLen, n8;
	int nSteps, nRem;
	int indFirst, indLast, indF, bufSize;
	bool t;
	unsigned int upos;

	char inputFile[1000], outputFile[1000], ioFolder[1000], fileName[100];
	FILE* fi, * fo;
	char* vi, * vo, * vtemp;
	Ipp32s* vIndex;
	Ipp8u* vBuffer;

	long long fileSize;
};

blockData::blockData(char* ioFolder_in, int nL_in, int nW_in) {
	fi = nullptr;
	fo = nullptr;
	vi = nullptr;
	vo = nullptr;
	vtemp = nullptr;
	memset(fileName, 0, sizeof(char) * 100);

	sprintf(ioFolder, ioFolder_in);
	nL = nL_in;
	nW = nW_in;
}

int blockData::freeMemory() {
	if (vi != nullptr) { free(vi); vi = nullptr;}
	if (vo != nullptr) { free(vo); vo = nullptr;}
	if (vtemp != nullptr) { free(vtemp); vtemp = nullptr;}
	if (fi != nullptr) { fclose(fi); fi = nullptr;}
	if (fo != nullptr) { fclose(fo); fo = nullptr;}
	return 0;
}

blockData::~blockData() {
	freeMemory();
}

int blockData::setFileNames() {
	int ia;
	memset(fileName, 0, 100 * sizeof(char));
	sprintf(fileName, "orig_");
	ia = indF;
	for (int j = 0; j < ndigits; j++) {
		sprintf(fileName + 5 + j, "%1x", (ia & 15));
		ia >>= 4;
	}
	sprintf(fileName + 5 + ndigits, ".bin");

	sprintf(inputFile, "%s/%s", ioFolder, fileName);

	fileName[0] = 's';
	fileName[1] = 'o';
	fileName[2] = 'r';
	fileName[3] = 't';

	sprintf(outputFile, "%s/%s", ioFolder, fileName);

	return 0;
}

int blockData::startProcessing(int indexFile) {
	indF = indexFile;

	iprefix = (2 * nW) % 8;
	while (iprefix + 8 <= 16) {
		iprefix += 8;
	}
	n8 = (2 * nW - iprefix) / 8;
	
	ndigits = iprefix / 4;
	if (iprefix % 4 > 0)ndigits++;

	sLen = 4 + n8;

	if (sLen >= 8) {
		is32 = false;
		nSteps = sLen / 8;
		nRem = sLen % 8;
	}
	else {
		is32 = true;
		nSteps = sLen / 4;
		nRem = sLen % 4;
	}

	

	setFileNames();

	printf("%s\n", inputFile);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open file \"%s\".\n", inputFile);
		return -1;
	}
	
#ifdef _WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);
#else
	fseeko(fi, 0, SEEK_END);
	fileSize = ftello(fi);
	fseeko(fi, 0, SEEK_SET);
#endif

	vi = (char*)malloc(sizeof(char) * fileSize);
	vtemp = (char*)malloc(sizeof(char) * fileSize);

	fread(vi, sizeof(char), fileSize, fi);
	fclose(fi); fi = nullptr;

	nE = (int)(fileSize / (long long)sLen);

	memcpy(vtemp, vi, sizeof(char) * fileSize);

	vIndex = ippsMalloc_32s(nE);
	if (is32) {
		ippsSortRadixIndexGetBufferSize(nE, IppDataType::ipp32u, &bufSize);
		vBuffer = ippsMalloc_8u(bufSize);

		ippsSortRadixIndexAscend_32u((Ipp32u*)(vi + sLen - 4), sLen, vIndex, nE, vBuffer);
		for (int i = 0; i < nE; i++) {
			memcpy(vi + i * sLen, vtemp + vIndex[i] * sLen, sLen * sizeof(char));
		}

		for (int ns = 1; ns < nSteps; ns++) {
			for (int j = 0; j < nE; j++) {
				indFirst = j;

				for (indLast = indFirst + 1; indLast < nE; indLast++) {
					t = true;
					for (int k = 1; k <= ns; k++) {
						if (*(unsigned long long*)(vi + (indFirst + 1) * sLen - k * 4) != *(unsigned long long*)(vi + (indLast + 1) * sLen - k * 4)) {
							t = false;
							break;
						}
					}
					if (!t) {
						indLast--;
						break;
					}
				}
				if (indLast == nE) indLast = nE - 1;
				if (indLast == indFirst) continue;

				memcpy(vtemp, vi + sLen * indFirst, sizeof(char) * (indLast - indFirst + 1) * sLen);
				ippsSortRadixIndexAscend_32u((Ipp32u*)(vi + sLen * (indFirst + 1) - 4 * (ns + 1)), sLen, vIndex, indLast - indFirst +1, vBuffer);
				for (int i = 0; i < indLast - indFirst +1; i++) {
					memcpy(vi + (indFirst + i) * sLen, vtemp + vIndex[i] * sLen, sLen * sizeof(char));
				}
				j = indLast;
			}
		}

		if (nRem > 0) {
			for (int j = 0; j < nE; j++) {
				indFirst = j;

				for (indLast = indFirst + 1; indLast < nE; indLast++) {
					t = true;
					for (int k = 1; k <= nSteps; k++) {
						if (*(unsigned long long*)(vi + (indFirst + 1) * sLen - k * 4) != *(unsigned long long*)(vi + (indLast + 1) * sLen - k * 4)) {
							t = false;
							break;
						}
					}
					if (!t) {
						indLast--;
						break;
					}
				}
				if (indLast == nE) indLast = nE - 1;
				if (indLast == indFirst) continue;

				memcpy(vtemp, vi + sLen * indFirst, sizeof(char) * (indLast - indFirst + 1) * sLen);
				ippsSortRadixIndexAscend_32u((Ipp32u*)(vi + sLen * indFirst), sLen, vIndex, indLast - indFirst + 1, vBuffer);
				for (int i = 0; i < indLast - indFirst + 1; i++) {
					memcpy(vi + (indFirst + i) * sLen, vtemp + vIndex[i] * sLen, sLen * sizeof(char));
				}
				j = indLast;
			}
		}

	}
	else {
		ippsSortRadixIndexGetBufferSize(nE, IppDataType::ipp64u, &bufSize);
		vBuffer = ippsMalloc_8u(bufSize);
		ippsSortRadixIndexAscend_64u((Ipp64u*)(vi + sLen - 8), sLen, vIndex, nE, vBuffer);

		for (int i = 0; i < nE; i++) {
			memcpy(vi + i * sLen, vtemp + vIndex[i] * sLen, sLen * sizeof(char));
		}

		for (int ns = 1; ns < nSteps; ns++) {
			for (int j = 0; j < nE; j++) {
				indFirst = j;

				for (indLast = indFirst + 1; indLast < nE; indLast++) {
					t = true;
					for (int k = 1; k <= ns; k++) {
						if (*(unsigned long long*)(vi + (indFirst + 1) * sLen - k * 8) != *(unsigned long long*)(vi + (indLast + 1) * sLen - k * 8)) {
							t = false;
							break;
						}
					}
					if (!t) {
						indLast--;
						break;
					}
				}
				if (indLast == nE) indLast = nE - 1;
				if (indLast == indFirst) continue;

				memcpy(vtemp, vi + sLen * indFirst, sizeof(char) * (indLast - indFirst + 1) * sLen);
				ippsSortRadixIndexAscend_64u((Ipp64u*)(vi + sLen * (indFirst + 1) - 8 * (ns + 1)), sLen, vIndex, indLast - indFirst +1, vBuffer);
				for (int i = 0; i < indLast - indFirst + 1; i++) {
					memcpy(vi + (indFirst + i) * sLen, vtemp + vIndex[i] * sLen, sLen * sizeof(char));
				}
				j = indLast;
			}
		}

		if (nRem > 0) {
			for (int j = 0; j < nE; j++) {
				indFirst = j;

				for (indLast = indFirst + 1; indLast < nE; indLast++) {
					t = true;
					for (int k = 1; k <= nSteps; k++) {
						if (*(unsigned long long*)(vi + (indFirst + 1) * sLen - k * 8) != *(unsigned long long*)(vi + (indLast + 1) * sLen - k * 8)) {
							t = false;
							break;
						}
					}
					if (!t) {
						indLast--;
						break;
					}
				}
				if (indLast == nE) indLast = nE - 1;
				if (indLast == indFirst) continue;

				memcpy(vtemp, vi + sLen * indFirst, sizeof(char) * (indLast - indFirst + 1) * sLen);
				ippsSortRadixIndexAscend_64u((Ipp64u*)(vi + sLen * indFirst), sLen, vIndex, indLast - indFirst + 1, vBuffer);
				for (int i = 0; i < indLast - indFirst + 1; i++) {
					memcpy(vi + (indFirst + i) * sLen, vtemp + vIndex[i] * sLen, sLen * sizeof(char));
				}
				j = indLast;
			}
		}
	}
	fo = fopen(outputFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open file \"%s\".\n", outputFile);
		return -2;
	}

	

	fwrite(vi, sizeof(char), fileSize, fo);
	fclose(fo); fo = nullptr;
	
	return 0;
}

int main(int argc, char** argv) {
	char inOutFolder[1000], sSeed[1000];
	int nL, nW, ipref, nFiles;
	
	if (argc != 3) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: input_output_folder seed\n\n");
		return -1;
	}

	sprintf(inOutFolder, argv[1]);
	printf("Input/Output folder: \"%s\".\n", inOutFolder);
	sprintf(sSeed, argv[2]);
	printf("Seed: \"%s\".\n", sSeed);

	nL = (int)strlen(sSeed);
	printf("Length: %i\n", nL);
	if (nL < 1) {
		printf("Error: wrong length.\n");
		return -2;
	}
	if (sSeed[0] != '1') {
		printf("Error: seed should start with 1.\n");
		return -3;
	}
	if (sSeed[nL-1] != '1') {
		printf("Error: seed should end with 1.\n");
		return -4;
	}
	nW = 0;
	for (int i = 0; i < nL; i++) {
		if (sSeed[i] == '1') nW++;
	}

	ipref = (2 * nW) % 8;
	while (ipref + 8 <= 16) {
		ipref += 8;
	}
		
	printf("Weight: %i\n", nW);
	printf("Prefix size: %i\n", ipref);
		
	nFiles = 1 << ipref;

#pragma omp parallel
	{
		blockData* bd;
		bd = new blockData(inOutFolder, nL, nW);

#pragma omp for
		for (int iq = 0; iq < nFiles; iq++) {
			bd->startProcessing(iq);
			bd->freeMemory();
		}
		delete bd; bd = nullptr;
	}

	return 0;
}

