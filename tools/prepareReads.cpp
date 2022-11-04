#include <stdio.h>
#include <fstream>
#include <stdlib.h>

//D:\Gen2022\data\ERR263486_1.fastq D:\Gen2022\data\err263486 100 10000000

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], commandLine[5000];
	int nLen, nRecords;
	if (argc != 5) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage:\n");
		printf("\t1) input_FASTQ_file output_folder_with_prefix\n");
		printf("\t2) output_folder_with_prefix\n");
		printf("\t3) length of reads\n");
		printf("\t4) maximum records\n");
		return -1;
	}

	sprintf(inputFile, argv[1]);
	printf("Input file: \"%s\".\n", inputFile);
	sprintf(outputFile, argv[2]);
	printf("Output file: \"%s\".\n", outputFile);
	nLen = atoi(argv[3]);
	printf("Length of reads: %i.\n", nLen);
	nRecords = atoi(argv[4]);
	printf("Maximum records: %i.\n\n", nRecords);


	sprintf(commandLine, "fastq2bin.exe %s %s.bin %i %i", inputFile, outputFile, nLen, nRecords);
	printf("%s\n", commandLine);
	sprintf(commandLine, "bin2m128.exe %s.bin %s_m128.bin", outputFile, outputFile);
	printf("%s\n", commandLine);

	return 0;
}

