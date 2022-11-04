#include <stdio.h>
#include <fstream>

const int groupSize = 1000000;

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], inputPrefix[1000];
	int nLen, n4, m, mrem, n128;
	FILE* fi, * fo;
	char* vi, *vo;
	char* vi_loc;
	unsigned int * vo_loc;
	int nRecords, nbr;
	unsigned int mOne[32];
	long long fileSize;
	int nChunks;
	int* iCount, *iCount4, *iCount128;

	if (argc != 2) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: input_output_folder_with_prefix\n\n");
		return -1;
	}

	for (int i = 0; i < 32; i++) {
		mOne[i] = 1 << i;
	}

	sprintf(inputPrefix, argv[1]);
	printf("Input/output folder with prefix: \"%s\".\n", inputPrefix);

	sprintf(inputFile, "%sinfo.bin", inputPrefix);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -4;
	}

	_fseeki64(fi, 0, SEEK_END);
	nChunks = _ftelli64(fi)/4;
	_fseeki64(fi, 0, SEEK_SET);

	printf("Number of chunks: %i\n", nChunks);

	iCount = (int*)malloc(sizeof(int) * nChunks);
	iCount4 = (int*)malloc(sizeof(int) * nChunks);
	iCount128 = (int*)malloc(sizeof(int) * nChunks);
	fread(iCount, sizeof(int), nChunks, fi);

	fclose(fi);

	for (int i = 0; i < nChunks; i++) {
		iCount4[i] = iCount[i] / 4;
		if (iCount[i] % 4 > 0)iCount4[i]++;

		iCount128[i] = iCount[i] / 32;
		if (iCount[i] % 32 > 0)iCount128[i]++;
	}



	sprintf(inputFile, "%sdata.bin", inputPrefix);
	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -5;
	}


	sprintf(outputFile, "%sm128.bin", inputPrefix);
	fo = fopen(outputFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -6;
	}

	vi = (char*)malloc(sizeof(char) * groupSize * 8);
	vo = (char*)malloc(sizeof(char) * groupSize * 16);

	int ng, nr, nL, nb, ia, nc;

	for (int k = 0; k < nChunks; k++) {
		printf("Chunk: %i\n", k);
		ng = iCount[k] / (32 * groupSize);
		nr = iCount[k] % (32 * groupSize);
		nb = nr / 32;
		nL = iCount[k] % 32;
		nc = nL / 16;
		if (nL % 16 > 0) nc++;

		for (int i = 0; i < ng; i++) {
			fread(vi, sizeof(char), groupSize * 8, fi);
			for (int j = 0; j < groupSize; j++) {
				vi_loc = vi + 8 * j;
				vo_loc = (unsigned int *)(vo + 16 * j);
				memset(vo_loc, 0, sizeof(char) * 16);
				for (int ii = 0; ii < 32; ii++) {
					ia = (vi_loc[ii / 4] >> (2 * (ii % 4))) & 3;
					vo_loc[ia] |= mOne[ii];
				}
			}
			fwrite(vo, sizeof(char), groupSize * 16, fo);
		}

		fread(vi, sizeof(char), nb * 8, fi);
		for (int j = 0; j < nb; j++) {
			vi_loc = vi + 8 * j;
			vo_loc = (unsigned int*)(vo + 16 * j);
			memset(vo_loc, 0, sizeof(char) * 16);
			for (int ii = 0; ii < 32; ii++) {
				ia = (vi_loc[ii / 4] >> (2 * (ii % 4))) & 3;
				vo_loc[ia] |= mOne[ii];
			}
		}
		fwrite(vo, sizeof(char), nb * 16, fo);

		if (nc > 0) {
			fread(vi, sizeof(char), 4*nc, fi);
			vi_loc = vi;
			vo_loc = (unsigned int*)(vo);
			memset(vo_loc, 0, sizeof(char) * 16);
			for (int ii = 0; ii < nL; ii++) {
				ia = (vi_loc[ii / 4] >> (2 * (ii % 4))) & 3;
				vo_loc[ia] |= mOne[ii];
			}
			fwrite(vo, sizeof(char), 16, fo);
		}
		

	}

	free(vi); vi = nullptr;
	free(vo); vo = nullptr;

	fclose(fi); fi = nullptr;
	fclose(fo); fo = nullptr;
	return 0;
}