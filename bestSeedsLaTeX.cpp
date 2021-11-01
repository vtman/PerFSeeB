#include <fstream>

int main(int argc, char* argv[]) {
	char myString[1000], fName[1000], inTemplate[1000], outTemplate[1000];
	char* outData, mLine[10000], mLine1[10000], mLine2[10000], mLine3[10000], mLine4[10000], mLine5[10000];
	char inputFolder[1000], outputFolder[1000];
	int* indLines;
	int nOne, bestOne[1000], smallestRead[1000], longestSeed[1000], posOne[1000], iRec[1000], countOne[1000];
	int fileFirst, fileLast, nRec, countLoc, qL;
	int nMismatches, lenMin, lenMax;
	int nWhole, nRem, nLines, totLines, nds;
	int* nDiff, * indDiff;
	bool isFound[1000];
	FILE* fi, * fo;

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

	nLines = 10000;
	outData = (char*)malloc(sizeof(char) * 10000 * nLines);
	indLines = (int*)malloc(sizeof(int) * nLines);

	totLines = 0;

	for (int i = lenMin; i <= lenMax; i++) {
		bestOne[i] = -1;
	}

	for (int i = 0; i < 1000; i++) {
		smallestRead[i] = 9999;
		longestSeed[i] = -1;
		isFound[i] = false;
	}

	nds = fileLast + 1;
	nDiff = (int*)malloc(sizeof(int) * (lenMax + 1));
	indDiff = (int*)malloc(sizeof(int) * (lenMax + 1) * nds);

	for (int i = lenMin; i <= lenMax; i++) {
		nDiff[i] = 0;
	}

	sprintf(inTemplate, "%s/block_m_%i_len_%%i_level_0.txt", inputFolder, nMismatches);
	sprintf(outTemplate, "%s/bestAll_%i.tex", outputFolder, nMismatches);

	fo = fopen(outTemplate, "w");
	
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

					qL = nWhole * nLen + nRem + 1;
					if (iL < smallestRead[countLoc]) {
						smallestRead[countLoc] = iL;
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

					qL = nWhole * nLen + nRem + 1;
					if (iL == smallestRead[countLoc]) {
						if (qL > longestSeed[countLoc]) {
							longestSeed[countLoc] = qL;
						}
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

					qL = nWhole * nLen + nRem + 1;
					if (iL == smallestRead[countLoc] && qL == longestSeed[countLoc] && countLoc >= 16 && countLoc<= 160) {
						sprintf(mLine1, "\\item $n_m=%i$, $w=%i$, $L_{min}=%i$, $n_b=%i$, $n_d=%i$,\\newline\n$b={}$\\texttt{", 
							nMismatches, countLoc, iL, nWhole, nRem + 1);
						for (int ii = 0; ii < nLen; ii++) {
							mLine2[ii] = '0' + iRec[(posOne[ip] + ii) % nLen];
						}
						mLine2[nLen] = '\0';
						sprintf(mLine3, "},\\newline\n$s={}$\\texttt{\\seqsplit{");
						for (int k = 0; k < nWhole+1; k++) {
							for (int ii = 0; ii < nLen; ii++) {
								mLine4[k*nLen + ii] = '0' + iRec[(posOne[ip] + ii) % nLen];
							}
						}
						mLine4[nWhole*nLen + nRem+1] = '\0';
						
						sprintf(mLine5, "}}\n\n");

						sprintf(mLine, "%s%s%s%s%s", mLine1, mLine2, mLine3, mLine4, mLine5);
						sprintf(outData + 10000*totLines, mLine);
						indLines[totLines] = countLoc;
						totLines++;
					}
				}
			}

		}
		fclose(fi); fi = nullptr;
		printf("nRec: %i\t%i\n", nLen, nRec);
	}

	for (int k = 0; k < 500; k++) {
		for (int i = 0; i < totLines; i++) {
			if (indLines[i] != k) continue;
			fprintf(fo, outData + 10000 * i);
		}
	}
	
	fclose(fo);

	free(nDiff); nDiff = nullptr;
	free(indDiff); indDiff = nullptr;
	
	return 0;
}