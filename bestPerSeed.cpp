#include <fstream>

int main(int argc, char* argv[]) {
	char myString[1000], fName[1000], inTemplate[1000], outTemplate[1000];
	char inputFolder[1000], outputFolder[1000];
	int nOne, bestOne[500], posOne[1000], iRec[1000], countOne[1000];
	int fileFirst, fileLast, nRec;
	int nMismatches, lenMin, lenMax;
	int nWhole, nRem, nds, countLoc;
	int* nDiff, * indDiff;
	FILE* fi, ** pfo, * fol;

	printf("Number of arguments: %i\n", argc);

	if (argc != 8) {
		printf("Usage:\n");
		printf("1) Input folder\n");
		printf("2) Output folder\n");
		printf("3) Number of mismatches\n");
		printf("4) File index (first)\n");
		printf("5) File index (last)\n");
		printf("6) Read length (min)\n");
		printf("7) Read length (max)\n");
		return -1;
	}

	sprintf(inputFolder, argv[1]);
	sprintf(outputFolder, argv[2]);
	nMismatches = atoi(argv[3]);
	fileFirst = atoi(argv[4]);
	fileLast = atoi(argv[5]);
	lenMin = atoi(argv[6]);
	lenMax = atoi(argv[7]);

	for (int i = lenMin; i <= lenMax; i++) {
		bestOne[i] = -1;
	}

	nds = fileLast + 1;
	nDiff = (int*)malloc(sizeof(int) * (lenMax + 1));
	indDiff = (int*)malloc(sizeof(int) * (lenMax + 1) * nds);

	for (int i = lenMin; i <= lenMax; i++) {
		nDiff[i] = 0;
	}

	sprintf(inTemplate, "%s/block_m_%i_len_%%i_level_0.txt", inputFolder, nMismatches);
	sprintf(outTemplate, "%s/out_m_%i_readLen_%%i.txt", outputFolder, nMismatches);
	sprintf(fName, "%s/res_%i.txt", outputFolder, nMismatches);
	
	fol = fopen(fName, "w");

	pfo = (FILE**)malloc(sizeof(FILE*) * (lenMax + 1));

	for (int i = lenMin; i <= lenMax; i++) {
		sprintf(fName, outTemplate, i);
		pfo[i] = fopen(fName, "w");
	}


	for (int nLen = fileFirst; nLen <= fileLast; nLen++) {
		sprintf(fName, inTemplate, nLen);
		fi = fopen(fName, "r");
		if (fi == nullptr) {
			printf("%i\t no file\n", nLen);
			continue;
		}

		nRec = 0;
		while(fscanf(fi, "%[^\n]\n", myString) > 0) {
			if (strlen(myString) != nLen) break;
			
			nRec++;
			nOne = 0;
			for (int j = 0; j < nLen; j++) {
				if (myString[j] == '1') {
					iRec[j] = 1;
					posOne[nOne] = j;
					nOne++;
				}
				else {
					iRec[j] = 0;
				}
				countOne[j] = nOne;
			}
			for (int j = 0; j < nLen; j++) {
				countOne[j + nLen] = countOne[j] + nOne;
				iRec[j + nLen] = iRec[j];
			}

			for (int iL = lenMin; iL <= lenMax; iL++) {
				nWhole = iL / nLen - 1;
				if (nWhole < 1) continue;
				nRem = iL % nLen;

				for (int ip = 0; ip < nOne; ip++) {
					if (iRec[posOne[ip] + nRem] == 0)continue;
					countLoc = nWhole * nOne + countOne[posOne[ip] + nRem] - countOne[posOne[ip]] + 1;
					if (countLoc > bestOne[iL]) {
						bestOne[iL] = countLoc;
					}
				}
			}

		}
		fclose(fi); fi = nullptr;
		printf("nRec: %i\t%i\n", nLen, nRec);
	}

	for (int nLen = fileFirst; nLen <= fileLast; nLen++) {
		sprintf(fName, inTemplate, nLen);
		fi = fopen(fName, "r");
		if (fi == nullptr) {
			printf("%i\t no file\n", nLen);
			continue;
		}

		nRec = 0;
		while (fscanf(fi, "%[^\n]\n", myString) > 0) {
			if (strlen(myString) != nLen) break;

			nRec++;
			nOne = 0;
			for (int j = 0; j < nLen; j++) {
				if (myString[j] == '1') {
					iRec[j] = 1;
					posOne[nOne] = j;
					nOne++;
				}
				else {
					iRec[j] = 0;
				}
				countOne[j] = nOne;
			}
			for (int j = 0; j < nLen; j++) {
				countOne[j + nLen] = countOne[j] + nOne;
				iRec[j + nLen] = iRec[j];
			}


			for (int iL = lenMin; iL <= lenMax; iL++) {
				nWhole = iL / nLen - 1;
				if (nWhole < 1) continue;
				nRem = iL % nLen;

				for (int ip = 0; ip < nOne; ip++) {
					if (iRec[posOne[ip] + nRem] == 0)continue;
					countLoc = nWhole * nOne + countOne[posOne[ip] + nRem] - countOne[posOne[ip]] + 1;
					if (countLoc < bestOne[iL])continue;

					if (nDiff[iL] == 0) {
						indDiff[iL * nds] = nLen;
						nDiff[iL]++;
					}
					else {
						if (indDiff[iL * nds + nDiff[iL] - 1] != nLen) {
							indDiff[iL * nds + nDiff[iL]] = nLen;
							nDiff[iL]++;
						}
					}
					fprintf(pfo[iL], "%i\t%i\t%i\t%i\t", countLoc, nLen, nWhole, nRem + 1);
					for (int ii = 0; ii < nLen; ii++) {
						fprintf(pfo[iL], "%i", iRec[(posOne[ip] + ii) % nLen]);
					}
					fprintf(pfo[iL], "\t");
					for (int kk = 0; kk < nWhole; kk++) {
						for (int ii = 0; ii < nLen; ii++) {
							fprintf(pfo[iL], "%i", iRec[(posOne[ip] + ii) % nLen]);
						}
					}
					for (int ii = 0; ii <= nRem; ii++) {
						fprintf(pfo[iL], "%i", iRec[(posOne[ip] + ii) % nLen]);
					}
					fprintf(pfo[iL], "\n");
				}
			}

		}
		fclose(fi); fi = nullptr;
		printf("nRec: %i\t%i\n", nLen, nRec);
	}

	for (int iL = lenMin; iL <= lenMax; iL++) {
		if (bestOne[iL] < 1) continue;
		fprintf(fol, "%i\t%i\t%i\t|", iL, bestOne[iL], nDiff[iL]);
		for (int i = 0; i < nDiff[iL]; i++) {
			fprintf(fol, "\t%i", indDiff[iL * nds + i]);
		}
		fprintf(fol, "\n");
	}

	fclose(fol);

	for (int i = lenMin; i <= lenMax; i++) {
		fclose(pfo[i]);
	}

	free(nDiff); nDiff = nullptr;
	free(indDiff); indDiff = nullptr;
	free(pfo); pfo = nullptr;
	
	return 0;
}