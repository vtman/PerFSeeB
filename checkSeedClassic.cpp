#include <iostream>

int main(int argc, char* argv[]) {
	char cSeed[1000];
	int readLen, seedLen, nMismatches, nRows;
	int ind, *vRows, *iPos, indexI;
	bool isOk;
	int* vR;

	bool isFound;

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

	vRows = (int*)malloc(sizeof(int) * nRows * readLen);

	iPos = (int*)malloc(sizeof(int) * nMismatches);

	for (int i = 0; i < nRows * readLen; i++) {
		vRows[i] = 0;
	}

	for (int k = 0; k < seedLen; k++) {
		if (cSeed[k] != '0') {
			vRows[k] = 1;
		}
	}

	for (int n = 1; n < nRows; n++) {
		for (int i = 0; i < seedLen; i++) {
			vRows[n * readLen + i + n] = vRows[i];
		}
	}

	isOk = true;

	for (int i = 0; i < nMismatches; i++) {
		iPos[i] = i;
	}
	indexI = nMismatches - 1;

	for (int k = 0; k < nRows; k++) {
		vR = vRows + k * readLen;
		for (int i = 0; i < readLen; i++) {
			if (vR[i] == 0) {
				printf("0");
			}
			else {
				printf("1");
			}
		}
		printf("\n");
	}

	do {

		isOk = false;
		for (int k = 0; k < nRows; k++) {
			vR = vRows + k * readLen;
			isFound = true;
			for (int i = 0; i < nMismatches; i++) {
				if (vR[iPos[i]] == 1) {
					isFound = false;
					break;
				}
			}
			if (isFound) {
				isOk = true;
				break;
			}
		}

		if (!isOk) {
			break;
		}
		
		do{
			iPos[indexI]++;
			if (iPos[indexI] <= readLen - nMismatches + indexI)break;
			indexI--;
		} while (indexI >= 0);
		if (indexI < 0) break;
		if (indexI < nMismatches - 1) {
			do {
				iPos[indexI + 1] = iPos[indexI] + 1;
				indexI++;
			} while (indexI < nMismatches - 1);
		}

	} while (true);


	if (isOk) {
		printf("The seed can be used.\n");
	}
	else {
		ind = 0;
		for (int i = 0; i < nMismatches; i++) {
			while (iPos[i] > ind) {
				printf("_");
				ind++;
			}
			printf("X");
			ind++;
		}
		while (ind < readLen) {
			printf("_");
			ind++;
		}
		printf("\n");
		for (int i = 0; i < nMismatches; i++) {
			printf("%i\t", iPos[i] + 1);
		}
		printf("\n");
		printf("The seed cannot be used.\n");
	}
	
	free(vRows); vRows = nullptr;
	free(iPos); iPos = nullptr;
	return 0;
}
