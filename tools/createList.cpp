#include <fstream>
#include <string.h>
#include <stdlib.h>

//C:\Temp2\Genome\T2T\T2T_ D:\Gen2022\seed1 111111111100011101100100100111010011100010100101000010100110000101111000000011

//D:\BMC\Ref\T2T_ D:\BMC\temp 111011010100110101100100100101010110001



const int itemsPerBlock = 50000;

int main(int argc, char** argv) {
	char inputFile[1000], inputPrefix[1000], outputFolder[1000], outputFile[1000], sSeed[1000];
	long long fileSize;
	FILE* fiInfo, * fiBin, *fo;
	int* viCount, *viCount16, *viRem, nChunks, *indOne;
	int** pvD, *iValue, *iRes;
	int nL, nW, ipref, n8;
	int ndigits, nFiles, ia, sLen;
	int *countItems;
	int ntot, ib, ir, ij, iap, indF;
	unsigned int upos, uu, startPos;
	unsigned char *uWord;
	char* fileNames, ** pvo, * vo;
	
	fiInfo = nullptr;
	fiBin = nullptr;

	if (argc != 4) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: input_folder_with_prefix output_folder seed\n\n");
		return -1;
	}

	sprintf(inputPrefix, argv[1]);
	printf("Input path: \"%s\".\n", inputPrefix);
	sprintf(outputFolder, argv[2]);
	printf("Output folder: \"%s\".\n", outputFolder);
	sprintf(sSeed, argv[3]);
	printf("Seed: \"%s\".\n", sSeed);

	uWord = (unsigned char*)malloc(sizeof(char) * 1000);

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
	n8 = (2 * nW - ipref) / 8;
	
	printf("Weight: %i\n", nW);
	indOne = (int*)malloc(sizeof(int) * nW);
	nW = 0;
	for (int i = 0; i < nL; i++) {
		if (sSeed[i] == '1') {
			indOne[nW] = i;
			nW++;
		}
	}

	printf("Prefix size: %i\n", ipref);
	printf("Number of chars: %i\n", n8);

	ndigits = ipref / 4;
	if (ipref % 4 > 0)ndigits++;

	nFiles = 1 << ipref;

	sLen = 4 + n8;

	countItems = (int*)malloc(sizeof(int) * nFiles);

	pvo = (char**)malloc(sizeof(char*) * nFiles);
	for (int i = 0; i < nFiles; i++) {
		pvo[i] = (char*)malloc(sizeof(char) * itemsPerBlock * sLen);
		if (pvo[i] == nullptr) {
			printf("Error: cannot allocate memory (%i)\n", i);
			return -7;
		}
	}

	fileNames = (char*)malloc(sizeof(char) * 20 * nFiles);

	memset(fileNames, 0, sizeof(char) * 20 * nFiles);

	for (int i = 0; i < nFiles; i++) {
		sprintf(fileNames + 20 * i, "orig_");
		ia = i;
		for (int j = 0; j < ndigits; j++) {
			sprintf(fileNames + 20 * i + 5 + j, "%1x", (ia & 15));
			ia >>= 4;
		}
		sprintf(fileNames + 20 * i + 5+ndigits, ".bin");
	}

	for (int i = 0; i < nFiles; i++) {
		sprintf(outputFile, "%s/%s", outputFolder, fileNames + 20 * i);
		fo = fopen(outputFile, "wb");
		if (fo == nullptr) {
			printf("Error: cannot open file \"%s\"\n", outputFile);
			return -5;
		}
		fclose(fo); fo = nullptr;
	}

	for (int i = 0; i < nFiles; i++) {
		countItems[i] = 0;
	}

	iValue = (int*)malloc(sizeof(int) * nL);
	iRes = (int*)malloc(sizeof(int) * nW);
	
	sprintf(inputFile, "%sinfo.bin", inputPrefix);
	fiInfo = fopen(inputFile, "rb");
	if (fiInfo == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -2;
	}

#ifdef _WIN32
	_fseeki64(fiInfo, 0, SEEK_END);
	fileSize = _ftelli64(fiInfo);
	_fseeki64(fiInfo, 0, SEEK_SET);
#else
	fseeko(fiInfo, 0, SEEK_END);
	fileSize = ftello(fiInfo);
	fseeko(fiInfo, 0, SEEK_SET);
#endif

	nChunks = (int)(fileSize / 4);
	printf("Number of chunks: %i\n", nChunks);
	viCount = (int*)malloc(sizeof(int) * nChunks);
	viCount16 = (int*)malloc(sizeof(int) * nChunks);
	viRem = (int*)malloc(sizeof(int) * nChunks);
	
	fread(viCount, sizeof(int), nChunks, fiInfo);
	fclose(fiInfo); fiInfo = nullptr;

	for (int i = 0; i < nChunks; i++) {
		viCount16[i] = viCount[i] / 16;
		if (viCount[i] % 16 == 0) {
			viRem[i] = 16;
		}
		else {
			viCount16[i]++;
			viRem[i] = viCount[i] % 16;
		}
	}

	pvD = nullptr;
	pvD = (int**)malloc(sizeof(int*) * nChunks);
	for (int i = 0; i < nChunks; i++) {
		pvD[i] = nullptr;
		pvD[i] = (int*)malloc(sizeof(int) * 4 * viCount16[i]);
		if (pvD[i] == nullptr) {
			printf("Error: cannot allocate memory (%i, %i)\n", i + 1, viCount16[i] * 4);
			return -1;
		}
	}

	sprintf(inputFile, "%sdata.bin", inputPrefix);
	fiBin = fopen(inputFile, "rb");
	if (fiBin == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -3;
	}
		
	for (int i = 0; i < nChunks; i++) {
		printf("Chunk: %i\n", i);
		fread(pvD[i], sizeof(int), viCount16[i], fiBin);
	}

	startPos = 0;
	for (int ic = 0; ic < nChunks; ic++) {
		printf("Chunk %i\n", ic);
		upos = startPos;
		ntot = viCount[ic] - nL + 1;
		for (int i = 0; i < nL; i++) {
			ib = i / 16;
			ir = (i % 16) * 2;
			iValue[i] = (pvD[ic][ib] >> ir) & 3;
		}
		
		for (int i = 0; i < ntot; i++) {
			uu = upos + (unsigned int)(i);
			
			if (i % 10000000 == 0)printf("%i\n", i);
			for (int j = 0; j < nW; j++) {
				ij = (i + indOne[j]) % nL;
				iRes[j] = iValue[ij];
			}
			for (int j = 0; j < n8; j++) {
				iap = 0;
				for (int k = 0; k < 4; k++) {
					iap |= (iRes[4 * j + k] << (2 * k));
				}
				uWord[4 + j] = *(unsigned char*)(&iap);
			}
			*(unsigned int*)(uWord) = uu;
			indF = 0;
			for (int j = 4 * n8; j < nW; j++) {
				indF |= (iRes[j] << (2 * (j - 4 * n8)));
			}

			vo = pvo[indF] + sLen * countItems[indF];
			memcpy(vo, uWord, sLen * sizeof(char));
			countItems[indF]++;
			if (countItems[indF] == itemsPerBlock) {
				sprintf(outputFile, "%s/%s", outputFolder, fileNames + 20 * indF);
				fo = fopen(outputFile, "ab");
				fwrite(pvo[indF], sizeof(char), itemsPerBlock * sLen, fo);
				fclose(fo); fo = nullptr;
				countItems[indF] = 0;
			}
			ib = (i + nL) / 16;
			ir = ((i + nL) % 16) * 2;
			iValue[i % nL] = (pvD[ic][ib] >> ir) & 3;
		}
		startPos += (unsigned int)(viCount[ic]);
	}

	for (int i = 0; i < nFiles; i++) {
		if (countItems[i] == 0)continue;
		sprintf(outputFile, "%s/%s", outputFolder, fileNames + 20 * i);
		fo = fopen(outputFile, "ab");
		fwrite(pvo[i], sizeof(char), countItems[i] * sLen, fo);
		fclose(fo); fo = nullptr;
	}


	free(uWord); uWord = nullptr;
	free(iValue); iValue = nullptr;
	free(iRes); iRes = nullptr;

	free(viCount); viCount = nullptr;
	free(viCount16); viCount16 = nullptr;
	free(viRem); viRem = nullptr;

	free(fileNames); fileNames = nullptr;
	for (int i = 0; i < nFiles; i++) {
		free(pvo[i]); pvo[i] = nullptr;
	}
	free(pvo); pvo = nullptr;

	for (int i = 0; i < nChunks; i++) {
		if (pvD[i] != 0) {
			free(pvD[i]), pvD[i] = nullptr;
		}
	}
	free(pvD); pvD = nullptr;
	free(countItems); countItems = nullptr;

	fclose(fiBin); fiBin = nullptr;
	return 0;
}
