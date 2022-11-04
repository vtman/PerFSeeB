#include <stdio.h>
#include <fstream>

#include <emmintrin.h>
#include <immintrin.h>

//C:\Temp2\Genome\T2T\T2T_ D:\Gen2022\out\err263486_pos.bin D:\Gen2022\data\err263486_m128.bin D:\Gen2022\out\err263486_match.bin

int findBlock(unsigned int* uStart, int n, unsigned int upos) {
	int i1, i2, i3;
	i2 = n - 1;
	i1 = 0;
	if (upos >= uStart[i2]) return i2;

	do {
		i3 = (i1 + i2) / 2;
		if (i1 == i3) return i1;
		if (upos >= uStart[i3]) {
			i1 = i3;
		}
		else {
			i2 = i3;
		}
	} while (true);

}

int countOne(__m128i* vm, int n) {
	__m128i ma, mb;
	int itot;
	unsigned int uval;
	itot = 0;
	for (int i = 0; i < n; i++) {
		ma = _mm_or_si128(vm[i], _mm_bsrli_si128(vm[i], 8));
		mb = _mm_or_si128(ma, _mm_bsrli_si128(ma, 4));
		uval = _mm_extract_epi32(mb, 0);
		itot += _mm_popcnt_u32(uval);
	}
	return itot;
}

int printM128(FILE* fo, __m128i* md, int n) {
	unsigned int ia, ic, ig, it;
	for (int ik = 0; ik < n; ik++) {
		ia = _mm_extract_epi32(md[ik], 0);
		ic = _mm_extract_epi32(md[ik], 1);
		ig = _mm_extract_epi32(md[ik], 2);
		it = _mm_extract_epi32(md[ik], 3);

		for (int i = 0; i < 32; i++) {
			if ((ia & 1) == 1) { fprintf(fo, "A"); }
			else if ((ic & 1) == 1) { fprintf(fo, "C"); }
			else if ((ig & 1) == 1) { fprintf(fo, "G"); }
			else if ((it & 1) == 1) { fprintf(fo, "T"); }
			else { fprintf(fo, " "); }
			ia >>= 1;
			ic >>= 1;
			ig >>= 1;
			it >>= 1;
		}
	}
	fprintf(fo, "\n");
	return 0;
}

int main(int argc, char** argv) {
	char inputFile[1000], inputPosFile[1000], inputReadFile[1000];
	char inputPrefix[1000], outputMatchFile[1000];

	int nLen, n128;
	int nRecords, nChunks;
	int ik, mtot, m1, ind, ib, nOne, mbest;
	int* iCount, * iCount128;

	unsigned int* vmTemp, *vmTemp2, *vStartPos;
	unsigned int upos1, uh, ur;

	FILE* fi, * fo, *fir;
		
	long long fileSize;
		
	__m128i** pvm, *mr, *md, *mq;

	vmTemp = (unsigned int*)malloc(sizeof(int) * 100000000);
	vmTemp2 = (unsigned int*)malloc(sizeof(int) * 100000000);

	pvm = nullptr;

	if (argc != 5) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage:\n");
		printf("1) Input reference folder with prefix\n");
		printf("2) Input binary file (positions)\n");
		printf("3) Input read file (m128)\n");
		printf("4) Output binary file (match)\n");
		return -1;
	}
		

	

	sprintf(inputPrefix, argv[1]);
	printf("Input reference folder with prefix: \"%s\".\n", inputPrefix);
	sprintf(inputPosFile, argv[2]);
	printf("Input binary file (positions): \"%s\".\n", inputPosFile);
	sprintf(inputReadFile, argv[3]);
	printf("Input read file (m128): \"%s\".\n", inputReadFile);
	sprintf(outputMatchFile, argv[4]);
	printf("Output binary file (match): \"%s\".\n", outputMatchFile);

	sprintf(inputFile, "%sinfo.bin", inputPrefix);

	//info

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -2;
	}

	_fseeki64(fi, 0, SEEK_END);
	nChunks = _ftelli64(fi)/4;
	_fseeki64(fi, 0, SEEK_SET);

	printf("Number of chunks: %i\n", nChunks);

	vStartPos = (unsigned int*)malloc(sizeof(int) * (nChunks + 1));

	vStartPos[0] = 0;

	iCount = (int*)malloc(sizeof(int) * nChunks);
	iCount128 = (int*)malloc(sizeof(int) * nChunks);
	fread(iCount, sizeof(int), nChunks, fi);

	fclose(fi); fi = nullptr;

	for (int i = 0; i < nChunks; i++) {
		vStartPos[i + 1] = vStartPos[i] + (unsigned int)iCount[i];
	}

	for (int i = 0; i < nChunks; i++) {
		iCount128[i] = iCount[i] / 32;
		if (iCount[i] % 32 > 0)iCount128[i]++;
	}

	//reference bin
	sprintf(inputFile, "%sm128.bin", inputPrefix);
	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -3;
	}

	pvm = (__m128i**)malloc(sizeof(__m128i*) * nChunks);

	for (int k = 0; k < nChunks; k++) {
		pvm[k] = nullptr;
		pvm[k] = (__m128i*)malloc(sizeof(__m128i) * iCount128[k]);
		if (pvm[k] == nullptr) {
			printf("Error: cannot allocate memory.\n");
			return -4;
		}
		fread(pvm[k], sizeof(__m128i), iCount128[k], fi);
	}

	fclose(fi); fi = nullptr;

	//processing;

	fo = fopen(outputMatchFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputMatchFile);
		return -5;
	}

	fi = fopen(inputPosFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputPosFile);
		return -6;
	}

	fir = fopen(inputReadFile, "rb");
	if (fir == nullptr) {
		printf("Error: cannot open read file \"%s\"\n", inputReadFile);
		return -7;
	}

	fread(&nLen, sizeof(int), 1, fir);
	_fseeki64(fir, 0, SEEK_END);
	fileSize = _ftelli64(fir);
	_fseeki64(fir, 4, SEEK_SET);
	
	n128 = nLen / 32;
	if (nLen % 32 > 0) n128++;

	mr = (__m128i*)malloc(sizeof(__m128i) * 2*n128);
	md = (__m128i*)malloc(sizeof(__m128i) * n128);
	mq = (__m128i*)malloc(sizeof(__m128i) * n128);

	nRecords = (int)((fileSize - 4) /((long long)(n128 * 16)));

	int maxC;

	for(ik = 0; ik < nRecords; ik++){
		fread(mr, sizeof(__m128i), n128, fir);

		maxC = -1;
		mbest = 0;
		fread(&ind, sizeof(int), 1, fi);
		fread(&mtot, sizeof(int), 1, fi);
		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mtot, sizeof(int), 1, fo);
		

		if (mtot == 0) continue; 
			
		fread(&m1, sizeof(int), 1, fi);
		fread(vmTemp, sizeof(int), m1, fi);

		fwrite(&m1, sizeof(int), 1, fo);
		fwrite(vmTemp, sizeof(int), m1, fo);

		for (int i = 0; i < m1; i++) {
			ib = findBlock(vStartPos, nChunks, vmTemp[i]);
			upos1 = vmTemp[i] - vStartPos[ib];
			uh = upos1 >> 5;
			ur = upos1 & 31;
			for (int j = 0; j < n128; j++) {
				md[j] = _mm_or_si128(_mm_srli_epi32(pvm[ib][uh + j], ur), _mm_slli_epi32(pvm[ib][uh + j + 1], 32 - ur));
				mq[j] = _mm_and_si128(md[j], mr[j]);
			}

			nOne = countOne(mq, n128);
			if (nOne < maxC)continue;
			if (nOne > maxC) {
				maxC = nOne;
				mbest = 0;
			}
			vmTemp2[mbest] = vmTemp[i];
			mbest++;
		}

		fwrite(&maxC, sizeof(int), 1, fo);
		fwrite(&mbest, sizeof(int), 1, fo);
		fwrite(vmTemp2, sizeof(int), mbest, fo);
	}
	
	free(mr); mr = nullptr;
	free(md); md = nullptr;
	free(mq); mq = nullptr;
	fclose(fi); fi = nullptr;
	fclose(fir); fir = nullptr;
	fclose(fo); fo = nullptr;

	free(vmTemp); vmTemp = nullptr;
	free(vmTemp2); vmTemp2 = nullptr;

	
	return 0;
}