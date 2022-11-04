#include <stdio.h>
#include <fstream>

//D:\Gen2022\data\ERR263486_1.fastq D:\Gen2022\data\err263486.bin 100 10000000

const int chunkSize = 10000000;
const int blockSize = 10000;

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000];
	char* vi;
	FILE* fi, * fo;
	long long upos, fileSize, curPosition, readSize;
	int lineStart[5], lineEnd[4], ik, nextra;
	int lastOk, kksave, numRecords, n, jj, jk;
	int maxRecords, nR, nLen, n4, nrem;
	unsigned char* vo, *vecO;
	bool isDone, isOk;

	if (argc != 5) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: inputFASTQfile outputBINfile numberOfLetters maxRecords\n\n");
		return -1;
	}

	sprintf(inputFile, argv[1]);
	sprintf(outputFile, argv[2]);
	nLen = atoi(argv[3]);
	maxRecords = atoi(argv[4]);

	printf("Input file: \"%s\"\n", inputFile);
	printf("Output file: \"%s\"\n", outputFile);
	printf("Number of letters: %i\n", nLen);
	if (nLen < 1 || nLen > 1000) {
		printf("Error: wrong value.\n");
		return -2;
	}
	printf("Max records: %i\n", maxRecords);
	if (maxRecords < 1) {
		printf("Error: wrong value.\n");
		return -3;
	}

	n4 = nLen / 4;
	nrem = nLen % 4;
	if (nrem == 0) {
		nrem = 4;
	}

	vi = (char*)malloc(sizeof(char) * chunkSize);
	vecO = (unsigned char*)malloc(sizeof(char) * 2*n4 * blockSize);

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

	fwrite(&nLen, sizeof(int), 1, fo);

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);

	printf("File size: %lld\n", fileSize);

	curPosition = _ftelli64(fi);
	readSize = fileSize - curPosition;
	if (chunkSize < readSize) readSize = chunkSize;
	fread(vi, sizeof(char), readSize, fi);

	nextra = 10;
	lastOk = 0;
	isDone = false;

	numRecords = 0;
	nR = 0;

	for (;;) {
		lineStart[0] = lastOk;
		isDone = false;
		kksave = -1;
		for (int k = 0; k < 4; k++) {
			for (int i = lineStart[k]; i < readSize; i++) {
				if (vi[i] == '\n' || vi[i] == '\r') {
					lineEnd[k] = i;
					isDone = true;
					break;
				}
			}
			if (!isDone) break;
			isDone = false;
			ik = lineEnd[k];
			for (int i = ik; i < readSize; i++) {
				if (!(vi[i] == '\n' || vi[i] == '\r')) {
					lineStart[k + 1] = i;
					isDone = true;
					break;
				}
			}
			
			if (!isDone)break;
			if (k < 3) isDone = false;
			kksave++;
		}

		if (!isDone) {
			curPosition += lastOk;

			if (curPosition + nextra > fileSize)break;
			_fseeki64(fi, curPosition, SEEK_SET);
			readSize = fileSize - curPosition;
			if (readSize > chunkSize) readSize = chunkSize;

			fread(vi, sizeof(char), readSize, fi);
			if (readSize < chunkSize) {
				vi[readSize] = ' ';
				readSize++;
			}
			lastOk = 0;
			continue;
		}
		
		n = lineEnd[1] - lineStart[1];
		if (n != nLen) {
			lastOk = lineStart[4];
			continue;
		}

		isOk = true;
		vo = vecO + 2*numRecords * n4;
		for (int j = 0; j < nLen; j++) {
			jj = j / 4;
			jk = 2*(j % 4);
			if (jk == 0) vo[jj] = 0;
			switch (vi[lineStart[1] + j]) {
			case 'A':
			case 'a':
				break;
			case 'C':
			case 'c':
				vo[jj] |= (1 << jk);
				break;
			case 'G':
			case 'g':
				vo[jj] |= (2 << jk);
				break;
			case 'T':
			case 't':
				vo[jj] |= (3 << jk);
				break;
			default:
				isOk = false;
				break;
			}

			if (!isOk)break;
		}

		vo = vecO + (2 * numRecords + 1) * n4;

		for (int j = 0; j < nLen; j++) {
			jj = j / 4;
			jk = 2 * (j % 4);
			if (jk == 0) vo[jj] = 0;
			switch (vi[lineStart[1] + nLen - 1 - j]) {
			case 'A':
			case 'a':
				vo[jj] |= (3 << jk);
				break;
			case 'C':
			case 'c':
				vo[jj] |= (2 << jk);
				break;
			case 'G':
			case 'g':
				vo[jj] |= (1 << jk);
				break;
			case 'T':
			case 't':
				break;
			default:
				isOk = false;
				break;
			}

			if (!isOk)break;
		}

		lastOk = lineStart[4];

		if (!isOk) continue;

		numRecords++;
		nR++;
		if (nR == maxRecords) {
			fwrite(vecO, sizeof(unsigned char), 2*n4 * numRecords, fo);
			numRecords = 0;
			break;
		}
		if (numRecords == blockSize) {
			fwrite(vecO, sizeof(unsigned char), 2*n4 * numRecords, fo);
			numRecords = 0;
			printf("%i\n", nR);
		}
	}

	if (numRecords > 0) {
		fwrite(vecO, sizeof(unsigned char), 2*n4 * numRecords, fo);
		numRecords = 0;
		printf("%i\n", nR);
	}

	fclose(fi); fi = nullptr;
	fclose(fo); fo = nullptr;
	free(vecO); vecO = nullptr;
	free(vi); vi = nullptr;

	return 0;
}