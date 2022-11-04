#include <stdio.h>
#include <fstream>
#include <stdlib.h>

//C:\Temp2\Genome\T2T\GCF_009914755.1_T2T-CHM13v2.0_genomic.fna C:\Temp2\Genome\T2T\T2T_

int main(int argc, char** argv) {
	char inputFile[1000], outputFile[1000], commandLine[5000];
	if (argc != 3) {
		printf("Error: wrong number of input parameters.\n");
		printf("Usage: input_FNA_file output_folder_with_prefix\n\n");
		return -1;
	}

	sprintf(inputFile, argv[1]);
	printf("Input file: \"%s\".\n", inputFile);
	sprintf(outputFile, argv[2]);
	printf("Output file: \"%s\".\n", outputFile);

	sprintf(commandLine, "fna2bin.exe %s %s", inputFile, outputFile);
	printf("%s\n", commandLine);
	sprintf(commandLine, "ref2m128.exe %s", outputFile);
	printf("%s\n", commandLine);

	return 0;
}

