#include <stdio.h>
#include <fstream>

#include <emmintrin.h>
#include <immintrin.h>

//D:\BMC\out2 MB_4_ _match.bin D:\BMC\out2\MB_4_all_match.bin 1 4

const int N = 100000000;

int mergeVectors(unsigned int* v1, unsigned int* v2, unsigned int* v3, int n1, int n2, int* n3) {
	int i1, i2, i3;
	if (n1 == 0) {
		memcpy(v3, v2, sizeof(int) * n2);
		*n3 = n2;
		return 0;
	}

	i1 = 0;
	i2 = 0;
	i3 = 0;
	do {
		if (v1[i1] < v2[i2]) {
			v3[i3] = v1[i1];
			i1++;
			i3++;
		}
		else if (v1[i1] > v2[i2]) {
			v3[i3] = v2[i2];
			i2++;
			i3++;
		}
		else {
			v3[i3] = v1[i1];
			i1++;
			i2++;
			i3++;
		}
		if (i1 == n1) {
			for (int k = i2; k < n2; k++) {
				v3[i3] = v2[k];
				i3++;
			}
			break;
		}

		if (i2 == n2) {
			for (int k = i1; k < n1; k++) {
				v3[i3] = v1[k];
				i3++;
			}
			break;
		}
	} while (true);

	*n3 = i3;
	return 0;
}


int main(int argc, char** argv) {
	char inputFolder[1000], inputPrefix[100], inputSuffix[100];
	int fileFirst, fileLast, nTotalFiles;
	char outputFile[1000], inputFile[1000];
	int indF;
	FILE** pfi, * fo;
	unsigned int *pvecU[3], * pvecB[3];
	int ind, mAllTot, *mAllLoc, mUniTot, mUniLoc, maxRes, maxLoc, countLoc, countTot;
	int iq1, iq2;

	if (argc != 7) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage:\n");
		printf("\t1) input folder\n");
		printf("\t2) input file prefix\n");
		printf("\t3) input file suffix\n");
		printf("\t4) output file\n");
		printf("\t5) index of the first file\n");
		printf("\t6) index of the last file\n");
		return -1;
	}

	sprintf(inputFolder, argv[1]);
	sprintf(inputPrefix, argv[2]);
	sprintf(inputSuffix, argv[3]);
	sprintf(outputFile, argv[4]);
	fileFirst = atoi(argv[5]);
	fileLast = atoi(argv[6]);

	printf("Input folder: \"%s\"\n", inputFolder);
	printf("Input prefix: \"%s\"\n", inputPrefix);
	printf("Input suffix: \"%s\"\n", inputSuffix);
	printf("Output file: \"%s\"\n", outputFile);
	printf("First file: %i\n", fileFirst);
	printf("Last file: %i\n", fileLast);

	pvecU[0] = (unsigned int*)malloc(sizeof(int) * N);
	pvecU[1] = (unsigned int*)malloc(sizeof(int) * N);
	pvecU[2] = (unsigned int*)malloc(sizeof(int) * N);

	pvecB[0] = (unsigned int*)malloc(sizeof(int) * N);
	pvecB[1] = (unsigned int*)malloc(sizeof(int) * N);
	pvecB[2] = (unsigned int*)malloc(sizeof(int) * N);

	if (pvecU[0] == nullptr || pvecU[1] == nullptr || pvecU[2] == nullptr || pvecB[0] == nullptr || pvecB[1] == nullptr || pvecB[2] == nullptr) {
		printf("Error: cannot allocate memory\n");
		return -3;
	}
	
	nTotalFiles = fileLast - fileFirst + 1;

	pfi = (FILE**)malloc(sizeof(FILE*) * nTotalFiles);
	mAllLoc = (int*)malloc(sizeof(int) * nTotalFiles);

	for (int i = 0; i < nTotalFiles; i++) {
		indF = fileFirst + i;
		sprintf(inputFile, "%s/%s%i%s", inputFolder, inputPrefix, indF, inputSuffix);
		pfi[i] = fopen(inputFile, "rb");
		if (pfi[i] == nullptr) {
			printf("Error: cannot open input file \"%s\"\n", inputFile);
			return -1;
		}
	}

	fo = fopen(outputFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -2;
	}

	int m1;

	bool isOut;
	isOut = false;

	do {
		mAllTot = 0;
		iq1 = 0;
		iq2 = 0;
		mUniTot = 0;
		for (int i = 0; i < nTotalFiles; i++) {
			if (fread(&ind, sizeof(int), 1, pfi[i]) != 1) {
				isOut = true;
				break;
			}
			fread(&(mAllLoc[i]), sizeof(int), 1, pfi[i]);
			if (mAllLoc[i] == 0)continue;
			mAllTot += mAllLoc[i];
			fread(&mUniLoc, sizeof(int), 1, pfi[i]);
			fread(pvecU[2], sizeof(int), mUniLoc, pfi[i]);
			m1 = mUniTot;
			mergeVectors(pvecU[iq1], pvecU[2], pvecU[1 - iq1], m1, mUniLoc, &mUniTot);
			iq1 = 1 - iq1;
		}
		if (isOut)break;

		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mAllTot, sizeof(int), 1, fo);
		if (mAllTot > 0) {
			fwrite(&mUniTot, sizeof(int), 1, fo);
			fwrite(pvecU[iq1], sizeof(int), mUniTot, fo);
		}

		if (mAllTot == 0)continue;

		countTot = 0;
		maxRes = -1;
		for (int i = 0; i < nTotalFiles; i++) {
			if (mAllLoc[i] == 0) continue;
			fread(&maxLoc, sizeof(int), 1, pfi[i]);
			fread(&countLoc, sizeof(int), 1, pfi[i]);
			fread(pvecB[2], sizeof(int), countLoc, pfi[i]);

			if (maxLoc < maxRes)continue;

			if (maxLoc > maxRes) {
				maxRes = maxLoc;
				iq2 = 0;
				countTot = 0;
			}
			
			m1 = countTot;
			mergeVectors(pvecB[iq2], pvecB[2], pvecB[1 - iq2], m1, countLoc, &countTot);
			iq2 = 1 - iq2;
		}

		fwrite(&maxRes, sizeof(int), 1, fo);
		fwrite(&countTot, sizeof(int), 1, fo);
		fwrite(pvecB[iq2], sizeof(int), countTot, fo);
	} while (true);


	fclose(fo); fo = nullptr;

	if (pfi != nullptr) {
		for (int i = 0; i < nTotalFiles; i++) {
			if (pfi[i] != nullptr) {
				fclose(pfi[i]); pfi[i] = nullptr;
			}
		}
		free(pfi); pfi = nullptr;
	}
	

	free(pvecU[0]); pvecU[0] = nullptr;
	free(pvecU[1]); pvecU[1] = nullptr;
	free(pvecU[2]); pvecU[2] = nullptr;

	free(pvecB[0]); pvecB[0] = nullptr;
	free(pvecB[1]); pvecB[1] = nullptr;
	free(pvecB[2]); pvecB[2] = nullptr;
	
	return 0;
}