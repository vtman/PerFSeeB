#include <fstream>
#include <stdlib.h>

//C:\Temp2\Genome\T2T\GCF_009914755.1_T2T-CHM13v2.0_genomic.fna C:\Temp2\Genome\T2T\T2T_

const int N = 10000000;
const int M = 10000;

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], outputPrefix[1000];
	FILE* fi, * foInfo, * foBin;
	char* vi;
	unsigned int* vo;
	int oC, oG, oT;
	long long fileSize, fpos, chunkSize;
	int ipos, opos, op1, op2;

	oC = 1;
	oG = 2;
	oT = 3;

	fi = nullptr;
	foInfo = nullptr;
	foBin = nullptr;

	if (argc != 3) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: path_to_FNA_file output_folder_with_prefix\n\n");
		return -1;
	}

	sprintf(inputFile, argv[1]);
	sprintf(outputPrefix, argv[2]);
	printf("Input file: \"%s\".\n", inputFile);
	printf("Output path: \"%s\".\n", outputPrefix);

	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -2;
	}

	sprintf(outputFile, "%sdata.bin", outputPrefix);
	foBin = fopen(outputFile, "wb");
	if (foBin == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -3;
	}

	sprintf(outputFile, "%sinfo.bin", outputPrefix);
	foInfo = fopen(outputFile, "wb");
	if (foInfo == nullptr) {
		printf("Error: cannot open output file \"%s\"\n", outputFile);
		return -4;
	}

	

	vi = (char*)malloc(sizeof(int) * N);
	vo = (unsigned int*)malloc(sizeof(int) * M);

	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);

	printf("File size: %lld\n", fileSize);

	chunkSize = N;

	
	op1 = 0;
	op2 = 0;
	opos = 0;
	vo[0] = 0;

	do {
		fpos = _ftelli64(fi);
		//printf("Posistion: %lld\n", fpos);
		if (fpos == fileSize) break;
		if (fpos + N > fileSize) {
			chunkSize = fileSize - fpos;
		}
		ipos = 0;
		fread(vi, sizeof(char), chunkSize, fi);
		if (vi[0] == '>') {
			while (vi[ipos] != '\n' && vi[ipos] != '\r') {
				ipos++;
			}
			while (vi[ipos] == '\n' || vi[ipos] == '\r') {
				ipos++;
			}
		}
		while (ipos < chunkSize) {
			if (vi[ipos] == '>') {
				_fseeki64(fi, ipos - chunkSize, SEEK_CUR);
				if (op2 == 0) {
					fwrite(vo, sizeof(int), op1, foBin);
				}
				else {
					fwrite(vo, sizeof(int), op1 + 1, foBin);
				}
				fwrite(&opos, sizeof(int), 1, foInfo);
				printf("%i\n", opos);
				opos = 0;
				op1 = 0;
				op2 = 0;
				vo[0] = 0;
				break;
			}
			if (vi[ipos] == 'A' || vi[ipos] == 'a') {
				
			}else if (vi[ipos] == 'C' || vi[ipos] == 'c') {
				vo[op1] |= (oC << op2);
			}else if (vi[ipos] == 'G' || vi[ipos] == 'g') {
				vo[op1] |= (oG << op2);
			}else if (vi[ipos] == 'T' || vi[ipos] == 't') {
				vo[op1] |= (oT << op2);
			}else {
				opos--;
				op2 -= 2;
			}
			opos++;
			op2 += 2;
			if (op2 == 32) {
				op1++;
				op2 = 0;
				if (op1 == M) {
					fwrite(vo, sizeof(int), op1, foBin);
					op1 = 0;
				}
				vo[op1] = 0;
			}
			ipos++;
		}
	} while (true);

	if (opos > 0) {
		if (op2 == 0) {
			fwrite(vo, sizeof(int), op1, foBin);
		}
		else {
			fwrite(vo, sizeof(int), op1 + 1, foBin);
		}
		fwrite(&opos, sizeof(int), 1, foInfo);
		printf("%i\n", opos);
	}

	free(vi); vi = nullptr;
	free(vo); vo = nullptr;
	

	fclose(fi); fi = nullptr;
	fclose(foInfo); foInfo = nullptr;
	fclose(foBin); foBin = nullptr;
	return 0;
}