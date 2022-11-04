#include <stdio.h>
#include <fstream>

const int blockSize = 10000000;

//D:\Gen2022\data\err263486.bin D:\Gen2022\data\err263486_m128.bin

//D:\BMC\Reads\err263486_rev.bin D:\BMC\Reads\err263486_rev_m128.bin

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000];
	int nLen, n4, m, mrem, n128;
	FILE* fi, * fo;
	char* vi;
	int nRecords, nbr;
	unsigned int* vo, mOne[32];
	long long fileSize;

	if (argc != 3) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: input_binary_file output_m128_binary_file\n\n");
		return -1;
	}

	for (int i = 0; i < 32; i++) {
		mOne[i] = 1 << i;
	}

	sprintf(inputFile, argv[1]);
	printf("Input file: \"%s\".\n", inputFile);
	sprintf(outputFile, argv[2]);
	printf("Output file: \"%s\".\n", outputFile);

	vi = (char*)malloc(sizeof(char) * blockSize);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -4;
	}

	fo = fopen(outputFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -5;
	}

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi) - 4;
	_fseeki64(fi, 0, SEEK_SET);
	fread(&nLen, sizeof(int), 1, fi);
	fwrite(&nLen, sizeof(int), 1, fo);

	n128 = nLen / 32;
	if (nLen % 32 > 0) n128++;

	

	n4 = nLen / 4;
	if (nLen % 4 > 0) n4++;
	nRecords = (int)((fileSize) / (long long)n4);
	nbr = (blockSize / n4);

	vo = (unsigned int*)malloc(sizeof(int) * n128 * 4 *nbr);

	m = nRecords / nbr + 1;
	mrem = nRecords % nbr;
	if (mrem == 0) {
		m--;
		mrem = nbr;
	}

	printf("Number of records: %i\n", nRecords);
	printf("Records per data block: %i\n", nbr);


	unsigned int* vo_loc;
	char* vi_loc;
	int mm;
	int j1, j2, jr1, jr2, iq;
	for (int k = 0; k < m; k++) {
		if (k < m - 1) {
			mm = nbr;
		}
		else {
			mm = mrem;
		}
		fread(vi, sizeof(char), mm * n4, fi);
		for (int i = 0; i < mm; i++) {
			vi_loc = vi + n4 * i;
			vo_loc = vo + 4 * n128 * i;
			memset(vo_loc, 0, sizeof(int) * n128 * 4);
			for (int j = 0; j < nLen; j++) {
				j1 = j / 4;
				jr1 = 2 * (j % 4);
				j2 = j / 32;
				jr2 = j % 32;
				iq = (vi_loc[j1] >> jr1) & 3;
				vo_loc[4 * j2 + iq] |= mOne[jr2];
			}
		}

		fwrite(vo, sizeof(int), n128 * 4 * mm, fo);
	}

	free(vi); vi = nullptr;
	free(vo); vo = nullptr;

	fclose(fi); fi = nullptr;
	fclose(fo); fo = nullptr;

	

	return 0;
}