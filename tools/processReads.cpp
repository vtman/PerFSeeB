#include <stdio.h>
#include <fstream>
#include <stdlib.h>

//C:\Temp2\Genome\T2T\T2T_ D:\Gen2022\data\err263486 C:\Temp2\Genome\temp D:\Gen2022\out C:\Temp2\Genome\temp\seeds.txt

//D:\BMC\Ref\T2T_ D:\BMC\Reads\err263486 D:\BMC\temp D:\BMC\out D:\BMC\seeds.txt

int main(int argc, char** argv) {
	char outputFolder[1000], commandLine[5000], tempFolder[1000];
	char inputRef[1000], inputRead[1000], seedFile[1000];
	char* vi;
	int nPairs, *tPos;
	long long fileSize;

	if (argc != 6) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage:\n");
		printf("\t1) input Ref folder + prefix\n");
		printf("\t2) input Read folder + prefix\n");
		printf("\t3) temp folder\n");
		printf("\t4) output folder\n");
		printf("\t5) seed text file\n\n");
		return -1;
	}

	sprintf(inputRef, argv[1]);
	printf("Input reference path: \"%s\".\n", inputRef);
	sprintf(inputRead, argv[2]);
	printf("Input read path: \"%s\".\n", inputRead);

	sprintf(tempFolder, argv[3]);
	printf("Temp folder: \"%s\".\n", tempFolder);

	sprintf(outputFolder, argv[4]);
	printf("Output folder: \"%s\".\n", outputFolder);
	sprintf(seedFile, argv[5]);
	printf("Seed file: \"%s\".\n", seedFile);

	FILE* fi;
	fi = fopen(seedFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open file \"%s\"\n", seedFile);
		return -1;
	}

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);

	vi = (char*)malloc(sizeof(char) * fileSize);
	fread(vi, sizeof(char), fileSize, fi);
	fclose(fi); fi = nullptr;

	for (int i = 0; i < fileSize; i++) {
		if (vi[i] <= 32) vi[i] = 0;
	}

	int u1, u2;
	u1 = 0;
	u2 = 0;
	while (vi[u1] == 0) {
		u1++;
	}

	vi[u2] = vi[u1];
	u1++;
	u2++;

	for (; u1 < fileSize; u1++) {
		if (vi[u2 - 1] == 0) {
			if (vi[u1] == 0) continue;
			vi[u2] = vi[u1];
			u2++;
		}else {
			vi[u2] = vi[u1];
			u2++;
		}
	}
	if (vi[u2] == 0) u2--;

	nPairs = 1;
	for (int i = 0; i < u2; i++) {
		if (vi[i] == 0)nPairs++;
	}
	nPairs /= 2;

	tPos = (int*)malloc(sizeof(int) * nPairs * 2);
	tPos[0] = 0;
	int n;

	n = 1;
	for (int i = 0; i < u2; i++) {
		if (vi[i] == 0) {
			tPos[n] = i + 1;
			n++;
			if (n == 2 * nPairs) break;
		}
	}

	for (int i = 0; i < nPairs; i++) {
		sprintf(commandLine, "createList.exe %s %s %s", inputRef, tempFolder, vi + tPos[2*i+1]);
		system(commandLine);
		//printf("%s\n", commandLine);

		sprintf(commandLine, "sortList.exe %s %s", tempFolder, vi + tPos[2 * i + 1]);
		system(commandLine);
		//printf("%s\n", commandLine);

		//sprintf(commandLine, "searchPositions.exe %s %s %s.bin %s %s/%s_pos.bin %s 0 999999999 4", inputRef, tempFolder, inputRead, tempFolder, outputFolder, vi + tPos[2 * i], vi + tPos[2 * i + 1]);
		sprintf(commandLine, "searchPositions.exe %s %s %s.bin %s %s/%s_pos.bin %s 0 99999 8", inputRef, tempFolder, inputRead, tempFolder, outputFolder, vi + tPos[2 * i], vi + tPos[2 * i + 1]);
		system(commandLine);
		//printf("%s\n", commandLine);

		sprintf(commandLine, "countMatch.exe %s %s/%s_pos.bin %s_m128.bin %s/%s_match.bin", inputRef, outputFolder, vi + tPos[2 * i], inputRead, outputFolder, vi + tPos[2 * i]);
		system(commandLine);
		//printf("%s\n", commandLine);

		sprintf(commandLine, "printMatch.exe %s/%s_match.bin %s/%s_stat.txt", outputFolder, vi + tPos[2 * i], outputFolder, vi + tPos[2 * i]);
		system(commandLine);
		//printf("%s\n", commandLine);

		//sprintf(commandLine, "del /q %s/orig*.bin", tempFolder);
		//system(commandLine);
		//printf("%s\n", commandLine);

		//sprintf(commandLine, "del /q %s/sort*.bin", tempFolder);
		//system(commandLine);
		//printf("%s\n", commandLine);

		printf("\n");

		

		//sprintf(commandLine, "countMatch.exe %s %s %s/%s.bin %s %s/_pos.bin %s", inputRef, outputFolder, vi + tPos[2 * i], tempFolder, outputFolder, vi + tPos[2 * i], vi + tPos[2 * i + 1]);
		//printf("%s\n", commandLine);

	
	}
	




	/*sprintf(commandLine, "fastq2bin.exe %s %s.bin %i %i", inputFile, outputFile, nLen, nRecords);
	printf("%s\n", commandLine);
	sprintf(commandLine, "bin2m128.exe %s.bin %s_m128.bin", outputFile, outputFile);
	printf("%s\n", commandLine);
	*/
	return 0;
}

