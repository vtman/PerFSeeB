#include <stdio.h>
#include <fstream>

//C:\Temp2\Genome\T2T\T2T_ D:\Gen2022\seed1 D:\Gen2022\data\err263486.bin C:\Temp2\Genome\temp D:\Gen2022\out\err263486_pos.bin 111111111100011101100100100111010011100010100101000010100110000101111000000011 0 9990000 4

//D:\BMC\Ref\T2T_ D:\BMC\temp D:\BMC\Reads\err263486.bin D:\BMC\temp D:\BMC\out\err263486_pos.bin 111011010100110101100100100101010110001 0 9990000 4

const int M = 100000000;

class searchClass {
public:
	searchClass();
	~searchClass();

	FILE* fi, * fo, **pfi;

	char inputPrefix[1000];
	char inputFile[1000], outputFile[1000];
	char inputReadFile[1000], outputResFile[1000], outputTempFolder[1000], strSeed[1000];
	char inputRefFolder[1000];
	char* vecO, * vi, fileName[100], ** pvd;
	unsigned char* vr;
	unsigned int * vStartBlock, * vCount;
	char* uWord;
	unsigned int* iValue, * iRes;

	int nRows, nRecords, nLen, n4, nrem;
	int indStart, indEnd, nChunks, nBlocks, nL, nW, ipref, n8;
	int ndigits, nFiles, sLen, ia, chunkSize;
	int* vnRef, nF, indFileFirst, indFileLast;
	int* indOne;
	int iqStart[100], nQ;
	int nMerge;

	unsigned int* pvm[2], * vmTemp;

	long long fileSize;

	int startProcessing(int n, char **v);
	int readInputParameters(int n, char **v);
	int printCommandLine();
	int checkInputParameters();
	int findParameters();
	int readReads();
	int readRefPref();
	int processReads();
	int processReads4();
	int processReads2();
	int processReads1();
	int compareRecords(char* r1, char* r2);
	int compareRecords4(char* r1, char* r2);
	int compareRecords2(char* r1, char* r2);
	int compareRecords1(char* r1, char* r2);
	int mergeBlocks();
	int printOutputData();
	int readRef();

	int printWord(char *vu, int n);

	int findLocation(char *vd, int n, int *ipos1, int *ipos2);
	int findLocation4(char* vd, int n, int* ipos1, int* ipos2);
	int findLocation2(char* vd, int n, int* ipos1, int* ipos2);
	int findLocation1(char* vd, int n, int* ipos1, int* ipos2);
	int mergeVectors(unsigned int* v1, unsigned int* v2, unsigned int* v3, int n1, int n2, int* n3);
};

searchClass::searchClass() {
	vStartBlock = nullptr;
	vCount = nullptr;

	fi = nullptr;
	fo = nullptr;
	pfi = nullptr;
	nChunks = 0;

	vr = nullptr;
	iValue = nullptr;
	iRes = nullptr;
	vnRef = nullptr;
	pvd = nullptr;

	uWord = (char*)malloc(sizeof(char) * 100);
	vmTemp = nullptr;
	pvm[0] = nullptr;
	pvm[1] = nullptr;

}

int findBlock(unsigned int* uStart, int n, unsigned int upos) {
	int i1, i2, i3;
	i2 = n - 1;
	i1 = 0;
	if (upos >= uStart[i2]) return i2;

	do {
		i3 = (i1 + i2) / 2;
		if (i1 == i3) return i1;
		if (upos >= uStart[i3]) {
			i1 = i3;
		}
		else {
			i2 = i3;
		}
	} while (true);
}

int searchClass::readRefPref() {
	sprintf(inputFile, "%sinfo.bin", inputPrefix);
	fi = fopen(inputFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open input file \"%s\"\n", inputFile);
		return -1;
	}

#ifdef _WIN32
	_fseeki64(fi, 0, SEEK_END);
	fileSize = _ftelli64(fi);
	_fseeki64(fi, 0, SEEK_SET);
#else
	fseeko(fi, 0, SEEK_END);
	fileSize = ftello(fi);
	fseeko(fi, 0, SEEK_SET);
#endif

	nBlocks = (int)(fileSize / 4);
	printf("Number of blocks: %i\n", nBlocks);
	vCount = (unsigned int*)malloc(sizeof(int) * nBlocks);

	fread(vCount, sizeof(int), nBlocks, fi);
	fclose(fi); fi = nullptr;

	vStartBlock = (unsigned int*)malloc(sizeof(int) * (nBlocks + 1));

	vStartBlock[0] = 0;
	for (int i = 0; i < nBlocks; i++) {
		vStartBlock[i + 1] = vStartBlock[i] + vCount[i];
	}
	return 0;
}

int searchClass::mergeVectors(unsigned int* v1, unsigned int* v2, unsigned int* v3, int n1, int n2, int* n3) {
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
		}else if (v1[i1] > v2[i2]) {
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

int searchClass::printWord(char* vu, int n) {
	int iaa, ibb;
	for (int i = 0; i < n; i++) {
		iaa = vu[i] & 15;
		ibb = (vu[i] >> 16) & 15;
		fprintf(fo, "%x%x", iaa, ibb);
	}
	return 0;
}

int searchClass::compareRecords(char* r1, char* r2) {
	for (int i = 0; i < nQ; i++) {
		if ((*(unsigned long long*)(r1 + iqStart[i])) < (*(unsigned long long*)(r2 + iqStart[i]))) return -1;
		if ((*(unsigned long long*)(r1 + iqStart[i])) > (*(unsigned long long*)(r2 + iqStart[i]))) return 1;
	}
	return 0;
}

int searchClass::compareRecords4(char* r1, char* r2) {
	for (int i = 0; i < nQ; i++) {
		if ((*(unsigned int*)(r1 + iqStart[i])) < (*(unsigned int*)(r2 + iqStart[i]))) return -1;
		if ((*(unsigned int*)(r1 + iqStart[i])) > (*(unsigned int*)(r2 + iqStart[i]))) return 1;
	}
	return 0;
}


int searchClass::compareRecords2(char* r1, char* r2) {
	for (int i = 0; i < nQ; i++) {
		if ((*(unsigned short*)(r1 + iqStart[i])) < (*(unsigned short*)(r2 + iqStart[i]))) return -1;
		if ((*(unsigned short*)(r1 + iqStart[i])) > (*(unsigned short*)(r2 + iqStart[i]))) return 1;
	}
	return 0;
}


int searchClass::compareRecords1(char* r1, char* r2) {
	for (int i = 0; i < nQ; i++) {
		if ((*(unsigned char*)(r1 + iqStart[i])) < (*(unsigned char*)(r2 + iqStart[i]))) return -1;
		if ((*(unsigned char*)(r1 + iqStart[i])) > (*(unsigned char*)(r2 + iqStart[i]))) return 1;
	}
	return 0;
}


searchClass::~searchClass() {
	if (vStartBlock != nullptr) { free(vStartBlock); vStartBlock = nullptr; }
	if (vCount != nullptr) { free(vCount); vCount = nullptr; }
	if (pfi != nullptr) {
		for (int i = 0; i < nChunks; i++) {
			if (pfi[i] != nullptr) { fclose(pfi[i]); pfi[i] = nullptr; }
		}
		free(pfi); pfi = nullptr;
	}
	if (pvm[0] != nullptr) { free(pvm[0]); pvm[0] = nullptr; }
	if (pvm[1] != nullptr) { free(pvm[1]); pvm[1] = nullptr; }
	if (vmTemp != nullptr) { free(vmTemp); vmTemp = nullptr; }
	if (uWord != nullptr) { free(uWord); uWord = nullptr; }
	if (fi != nullptr) { fclose(fi); fi = nullptr;}
	if (fo != nullptr) { fclose(fo); fo = nullptr; }
	if (vr != nullptr) { free(vr); vr = nullptr; }
	if (iValue != nullptr) { free(iValue); iValue = nullptr; }
	if (iValue != nullptr) { free(iValue); iValue = nullptr; }
	if (pvd != nullptr) {
		for (int k = 0; k < nF; k++) {
			if (pvd[k] != nullptr) { free(pvd[k]); pvd[k] = nullptr; }
		}
		free(pvd); pvd = nullptr;
	}
	if (vnRef != nullptr) { free(vnRef); vnRef = nullptr; }
}

int searchClass::findLocation(char* vd, int n, int* ipos1, int* ipos2) {
	int i1, i2, i3, i4, i5, ic1, ic2, ic3, ic5, iz1, iz2;;
	bool isZero1, isZero2;
	
	char *vd_loc;

	isZero1 = false;
	isZero2 = false;
	vd_loc = vd + 4;

	i1 = 0;
	ic1 = compareRecords(vd_loc, uWord);

	if (ic1 == 1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic1 == 0) {
		isZero1 = true;
		i3 = i1;
		iz1 = i1;
	}

	i2 = n - 1;
	ic2 = compareRecords(vd_loc + sLen * i2, uWord);

	if (ic2 == -1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic2 == 0) {
		isZero2 = true;
		i3 = i2;
		iz2 = i2;
	}

	i1 = 0;
	
	if (isZero1 == false && isZero2 == false) {
		do {
			i3 = (i1 + i2) / 2;
			if (i1 == i3)break;
			ic3 = compareRecords(vd_loc + sLen * i3, uWord);
			if (ic3 == 0) break;

			if (ic3 == 1) {
				i2 = i3;
			}
			else {
				i1 = i3;
			}
		} while (true);


		if (ic3 != 0) {
			*ipos1 = -1;
			return 0;
		}
	}

	if (isZero1 == false) {
		i4 = i3;
		do {
			i5 = (i1 + i4) / 2;
			if (i1 == i5)break;
			ic5 = compareRecords(vd_loc + sLen * i5, uWord);
			if (ic5 == -1) {
				i1 = i5;
			}else {
				i4 = i5;
			}
		} while (true);
		iz1 = i4;
	}

	if (isZero2 == false) {
		i4 = i3;
		do {
			i5 = (i2 + i4) / 2;
			if (i4 == i5)break;
			ic5 = compareRecords(vd_loc + sLen * i5, uWord);
			if (ic5 == 1) {
				i2 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz2 = i4;
	}

	*ipos1 = iz1;
	*ipos2 = iz2;

	return 0;
}

int searchClass::findLocation4(char* vd, int n, int* ipos1, int* ipos2) {
	int i1, i2, i3, i4, i5, ic1, ic2, ic3, ic5, iz1, iz2;;
	bool isZero1, isZero2;

	char* vd_loc;

	isZero1 = false;
	isZero2 = false;
	vd_loc = vd + 4;

	i1 = 0;
	ic1 = compareRecords4(vd_loc, uWord);

	if (ic1 == 1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic1 == 0) {
		isZero1 = true;
		i3 = i1;
		iz1 = i1;
	}

	i2 = n - 1;
	ic2 = compareRecords4(vd_loc + sLen * i2, uWord);

	if (ic2 == -1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic2 == 0) {
		isZero2 = true;
		i3 = i2;
		iz2 = i2;
	}

	i1 = 0;

	if (isZero1 == false && isZero2 == false) {
		do {
			i3 = (i1 + i2) / 2;
			if (i1 == i3)break;
			ic3 = compareRecords4(vd_loc + sLen * i3, uWord);
			if (ic3 == 0) break;

			if (ic3 == 1) {
				i2 = i3;
			}
			else {
				i1 = i3;
			}
		} while (true);


		if (ic3 != 0) {
			*ipos1 = -1;
			return 0;
		}
	}

	if (isZero1 == false) {
		i4 = i3;
		do {
			i5 = (i1 + i4) / 2;
			if (i1 == i5)break;
			ic5 = compareRecords4(vd_loc + sLen * i5, uWord);
			if (ic5 == -1) {
				i1 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz1 = i4;
	}

	if (isZero2 == false) {
		i4 = i3;
		do {
			i5 = (i2 + i4) / 2;
			if (i4 == i5)break;
			ic5 = compareRecords4(vd_loc + sLen * i5, uWord);
			if (ic5 == 1) {
				i2 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz2 = i4;
	}

	*ipos1 = iz1;
	*ipos2 = iz2;

	return 0;
}

int searchClass::findLocation2(char* vd, int n, int* ipos1, int* ipos2) {
	int i1, i2, i3, i4, i5, ic1, ic2, ic3, ic5, iz1, iz2;;
	bool isZero1, isZero2;

	char* vd_loc;

	isZero1 = false;
	isZero2 = false;
	vd_loc = vd + 4;

	i1 = 0;
	ic1 = compareRecords2(vd_loc, uWord);

	if (ic1 == 1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic1 == 0) {
		isZero1 = true;
		i3 = i1;
		iz1 = i1;
	}

	i2 = n - 1;
	ic2 = compareRecords2(vd_loc + sLen * i2, uWord);

	if (ic2 == -1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic2 == 0) {
		isZero2 = true;
		i3 = i2;
		iz2 = i2;
	}

	i1 = 0;

	if (isZero1 == false && isZero2 == false) {
		do {
			i3 = (i1 + i2) / 2;
			if (i1 == i3)break;
			ic3 = compareRecords2(vd_loc + sLen * i3, uWord);
			if (ic3 == 0) break;

			if (ic3 == 1) {
				i2 = i3;
			}
			else {
				i1 = i3;
			}
		} while (true);


		if (ic3 != 0) {
			*ipos1 = -1;
			return 0;
		}
	}

	if (isZero1 == false) {
		i4 = i3;
		do {
			i5 = (i1 + i4) / 2;
			if (i1 == i5)break;
			ic5 = compareRecords2(vd_loc + sLen * i5, uWord);
			if (ic5 == -1) {
				i1 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz1 = i4;
	}

	if (isZero2 == false) {
		i4 = i3;
		do {
			i5 = (i2 + i4) / 2;
			if (i4 == i5)break;
			ic5 = compareRecords2(vd_loc + sLen * i5, uWord);
			if (ic5 == 1) {
				i2 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz2 = i4;
	}

	*ipos1 = iz1;
	*ipos2 = iz2;

	return 0;
}

int searchClass::findLocation1(char* vd, int n, int* ipos1, int* ipos2) {
	int i1, i2, i3, i4, i5, ic1, ic2, ic3, ic5, iz1, iz2;;
	bool isZero1, isZero2;

	char* vd_loc;

	isZero1 = false;
	isZero2 = false;
	vd_loc = vd + 4;

	i1 = 0;
	ic1 = compareRecords1(vd_loc, uWord);

	if (ic1 == 1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic1 == 0) {
		isZero1 = true;
		i3 = i1;
		iz1 = i1;
	}

	i2 = n - 1;
	ic2 = compareRecords1(vd_loc + sLen * i2, uWord);

	if (ic2 == -1) {
		*ipos1 = -1;
		return 0;
	}
	else if (ic2 == 0) {
		isZero2 = true;
		i3 = i2;
		iz2 = i2;
	}

	i1 = 0;

	if (isZero1 == false && isZero2 == false) {
		do {
			i3 = (i1 + i2) / 2;
			if (i1 == i3)break;
			ic3 = compareRecords1(vd_loc + sLen * i3, uWord);
			if (ic3 == 0) break;

			if (ic3 == 1) {
				i2 = i3;
			}
			else {
				i1 = i3;
			}
		} while (true);


		if (ic3 != 0) {
			*ipos1 = -1;
			return 0;
		}
	}

	if (isZero1 == false) {
		i4 = i3;
		do {
			i5 = (i1 + i4) / 2;
			if (i1 == i5)break;
			ic5 = compareRecords1(vd_loc + sLen * i5, uWord);
			if (ic5 == -1) {
				i1 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz1 = i4;
	}

	if (isZero2 == false) {
		i4 = i3;
		do {
			i5 = (i2 + i4) / 2;
			if (i4 == i5)break;
			ic5 = compareRecords1(vd_loc + sLen * i5, uWord);
			if (ic5 == 1) {
				i2 = i5;
			}
			else {
				i4 = i5;
			}
		} while (true);
		iz2 = i4;
	}

	*ipos1 = iz1;
	*ipos2 = iz2;

	return 0;
}

int searchClass::processReads() {
	unsigned char* vr_loc;
	int j1, j2, indF, iap, indF_loc;
	int ipos1, ipos2;
	int myM, m1, m2, m3, ind, mtot, ib;

	for (int ir = 0; ir < nRows; ir++) {
		ind = ir + indStart;
		vr_loc = vr + n4 * ir;
		mtot = 0;
		for (int ii = 0; ii < nLen; ii++) {
			j1 = ii / 4;
			j2 = 2 * (ii % 4);
			iValue[ii] = (vr_loc[j1] >> j2) & 3;
		}
		m1 = 0;
		myM = 0;

		for (int ik = 0; ik <= nLen - nL; ik++) {
			for (int jj = 0; jj < nW; jj++) {
				iRes[jj] = iValue[ik + indOne[jj]];
			}

			indF = 0;
			for (int j = 4 * n8; j < nW; j++) {
				indF |= (iRes[j] << (2 * (j - 4 * n8)));
			}
			indF_loc = indF - indFileFirst;
			if (indF_loc < 0) continue;
			if (indF_loc >= nF)continue;

			for (int j = 0; j < n8; j++) {
				iap = 0;
				for (int k = 0; k < 4; k++) {
					iap |= (iRes[4 * j + k] << (2 * k));
				}
				uWord[j] = *(unsigned char*)(&iap);
			}
			findLocation(pvd[indF_loc], vnRef[indF_loc], &ipos1, &ipos2);
			if (ipos1 < 0) continue;
			m2 = 0;
			for (int ii = ipos1; ii <= ipos2; ii++) {
				vmTemp[m2] = *(unsigned int*)(pvd[indF_loc] + ii * sLen);
				ib = findBlock(vStartBlock, nBlocks, vmTemp[m2]);
				if (vStartBlock[ib] + ik > vmTemp[m2])continue;
				vmTemp[m2] -= ik;
				m2++;
			}
			if (m2 > 0) {
				mtot += m2;
				mergeVectors(pvm[myM], vmTemp, pvm[1 - myM], m1, m2, &m3);
				m1 = m3;
				myM = 1 - myM;
			}
		}

		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mtot, sizeof(int), 1, fo);
		if (m1 > 0) {
			fwrite(&m1, sizeof(int), 1, fo);
			fwrite(pvm[myM], sizeof(int), m1, fo);
		}
	}

	for (int k = 0; k < nF; k++) {
		free(pvd[k]); pvd[k] = nullptr;
	}
	free(pvd); pvd = nullptr;
	free(vnRef); vnRef = nullptr;
	return 0;
}



int searchClass::processReads4() {
	unsigned char* vr_loc;
	int j1, j2, indF, iap, indF_loc;
	int ipos1, ipos2;
	int myM, m1, m2, m3, ind, mtot, ib;

	for (int ir = 0; ir < nRows; ir++) {
		ind = ir + indStart;
		vr_loc = vr + n4 * ir;
		mtot = 0;
		for (int ii = 0; ii < nLen; ii++) {
			j1 = ii / 4;
			j2 = 2 * (ii % 4);
			iValue[ii] = (vr_loc[j1] >> j2) & 3;
		}
		m1 = 0;
		myM = 0;

		for (int ik = 0; ik <= nLen - nL; ik++) {
			for (int jj = 0; jj < nW; jj++) {
				iRes[jj] = iValue[ik + indOne[jj]];
			}

			indF = 0;
			for (int j = 4 * n8; j < nW; j++) {
				indF |= (iRes[j] << (2 * (j - 4 * n8)));
			}
			indF_loc = indF - indFileFirst;
			if (indF_loc < 0) continue;
			if (indF_loc >= nF)continue;

			for (int j = 0; j < n8; j++) {
				iap = 0;
				for (int k = 0; k < 4; k++) {
					iap |= (iRes[4 * j + k] << (2 * k));
				}
				uWord[j] = *(unsigned char*)(&iap);
			}
			findLocation4(pvd[indF_loc], vnRef[indF_loc], &ipos1, &ipos2);
			if (ipos1 < 0) continue;
			m2 = 0;
			for (int ii = ipos1; ii <= ipos2; ii++) {
				vmTemp[m2] = *(unsigned int*)(pvd[indF_loc] + ii * sLen);
				ib = findBlock(vStartBlock, nBlocks, vmTemp[m2]);
				if (vStartBlock[ib] + ik > vmTemp[m2])continue;
				vmTemp[m2] -= ik;
				m2++;
			}
			if (m2 > 0) {
				mtot += m2;
				mergeVectors(pvm[myM], vmTemp, pvm[1 - myM], m1, m2, &m3);
				m1 = m3;
				myM = 1 - myM;
			}
			
		}

		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mtot, sizeof(int), 1, fo);
		if (m1 > 0) {
			fwrite(&m1, sizeof(int), 1, fo);
			fwrite(pvm[myM], sizeof(int), m1, fo);
		}
	}

	for (int k = 0; k < nF; k++) {
		free(pvd[k]); pvd[k] = nullptr;
	}
	free(pvd); pvd = nullptr;
	free(vnRef); vnRef = nullptr;
	return 0;
}


int searchClass::processReads2() {
	unsigned char* vr_loc;
	int j1, j2, indF, iap, indF_loc;
	int ipos1, ipos2;
	int myM, m1, m2, m3, ind, mtot, ib;

	for (int ir = 0; ir < nRows; ir++) {
		ind = ir + indStart;
		vr_loc = vr + n4 * ir;
		mtot = 0;
		for (int ii = 0; ii < nLen; ii++) {
			j1 = ii / 4;
			j2 = 2 * (ii % 4);
			iValue[ii] = (vr_loc[j1] >> j2) & 3;
		}
		m1 = 0;
		myM = 0;

		for (int ik = 0; ik <= nLen - nL; ik++) {
			for (int jj = 0; jj < nW; jj++) {
				iRes[jj] = iValue[ik + indOne[jj]];
			}

			indF = 0;
			for (int j = 4 * n8; j < nW; j++) {
				indF |= (iRes[j] << (2 * (j - 4 * n8)));
			}
			indF_loc = indF - indFileFirst;
			if (indF_loc < 0) continue;
			if (indF_loc >= nF)continue;

			for (int j = 0; j < n8; j++) {
				iap = 0;
				for (int k = 0; k < 4; k++) {
					iap |= (iRes[4 * j + k] << (2 * k));
				}
				uWord[j] = *(unsigned char*)(&iap);
			}
			findLocation2(pvd[indF_loc], vnRef[indF_loc], &ipos1, &ipos2);
			if (ipos1 < 0) continue;
			m2 = 0;
			for (int ii = ipos1; ii <= ipos2; ii++) {
				vmTemp[m2] = *(unsigned int*)(pvd[indF_loc] + ii * sLen);
				ib = findBlock(vStartBlock, nBlocks, vmTemp[m2]);
				if (vStartBlock[ib] + ik > vmTemp[m2])continue;
				vmTemp[m2] -= ik;
				m2++;
			}
			if (m2 > 0) {
				mtot += m2;
				mergeVectors(pvm[myM], vmTemp, pvm[1 - myM], m1, m2, &m3);
				m1 = m3;
				myM = 1 - myM;
			}

		}

		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mtot, sizeof(int), 1, fo);
		if (m1 > 0) {
			fwrite(&m1, sizeof(int), 1, fo);
			fwrite(pvm[myM], sizeof(int), m1, fo);
		}
	}

	for (int k = 0; k < nF; k++) {
		free(pvd[k]); pvd[k] = nullptr;
	}
	free(pvd); pvd = nullptr;
	free(vnRef); vnRef = nullptr;
	return 0;
}

int searchClass::processReads1() {
	unsigned char* vr_loc;
	int j1, j2, indF, iap, indF_loc;
	int ipos1, ipos2;
	int myM, m1, m2, m3, ind, mtot, ib;

	for (int ir = 0; ir < nRows; ir++) {
		ind = ir + indStart;
		vr_loc = vr + n4 * ir;
		mtot = 0;
		for (int ii = 0; ii < nLen; ii++) {
			j1 = ii / 4;
			j2 = 2 * (ii % 4);
			iValue[ii] = (vr_loc[j1] >> j2) & 3;
		}
		m1 = 0;
		myM = 0;

		for (int ik = 0; ik <= nLen - nL; ik++) {
			for (int jj = 0; jj < nW; jj++) {
				iRes[jj] = iValue[ik + indOne[jj]];
			}

			indF = 0;
			for (int j = 4 * n8; j < nW; j++) {
				indF |= (iRes[j] << (2 * (j - 4 * n8)));
			}
			indF_loc = indF - indFileFirst;
			if (indF_loc < 0) continue;
			if (indF_loc >= nF)continue;

			for (int j = 0; j < n8; j++) {
				iap = 0;
				for (int k = 0; k < 4; k++) {
					iap |= (iRes[4 * j + k] << (2 * k));
				}
				uWord[j] = *(unsigned char*)(&iap);
			}
			findLocation1(pvd[indF_loc], vnRef[indF_loc], &ipos1, &ipos2);
			if (ipos1 < 0) continue;
			m2 = 0;
			for (int ii = ipos1; ii <= ipos2; ii++) {
				vmTemp[m2] = *(unsigned int*)(pvd[indF_loc] + ii * sLen);
				ib = findBlock(vStartBlock, nBlocks, vmTemp[m2]);
				if (vStartBlock[ib] + ik > vmTemp[m2])continue;
				vmTemp[m2] -= ik;
				m2++;
			}
			if (m2 > 0) {
				mtot += m2;
				mergeVectors(pvm[myM], vmTemp, pvm[1 - myM], m1, m2, &m3);
				m1 = m3;
				myM = 1 - myM;
			}

		}

		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mtot, sizeof(int), 1, fo);
		if (m1 > 0) {
			fwrite(&m1, sizeof(int), 1, fo);
			fwrite(pvm[myM], sizeof(int), m1, fo);
		}
	}

	for (int k = 0; k < nF; k++) {
		free(pvd[k]); pvd[k] = nullptr;
	}
	free(pvd); pvd = nullptr;
	free(vnRef); vnRef = nullptr;
	return 0;
}



int searchClass::readRef() {
	vnRef = (int*)malloc(sizeof(int) * nF);
	pvd = (char**)malloc(sizeof(char*) * nF);
	for (int k = 0; k < nF; k++) {
		sprintf(fileName, "sort_");
		ia = k + indFileFirst;
		for (int j = 0; j < ndigits; j++) {
			sprintf(fileName + 5 + j, "%1x", (ia & 15));
			ia >>= 4;
		}
		sprintf(fileName + 5 + ndigits, ".bin");

		sprintf(inputFile, "%s/%s", inputRefFolder, fileName);
		printf("%s\n", inputFile);

		fi = fopen(inputFile, "rb");
		_fseeki64(fi, 0, SEEK_END);
		fileSize = _ftelli64(fi);
		_fseeki64(fi, 0, SEEK_SET);

		vnRef[k] = fileSize / ((long long)sLen);
		pvd[k] = (char*)malloc(sizeof(char) * fileSize);
		if (pvd[k] == nullptr) {
			printf("Error: cannot allocate memory.\n");
			return -1;
		}
		fread(pvd[k], sizeof(char), fileSize, fi);
		fclose(fi); fi = nullptr;
	}

	return 0;
}

int searchClass::startProcessing(int n, char** v) {
	if (readInputParameters(n, v) != 0) return -1;
	if (checkInputParameters() != 0) return -2;
	if (findParameters() != 0) return -3;
	if (readRefPref() != 0) return -4;
	if (readReads() != 0) return -5;

	
	for (int ic = 0; ic < nChunks; ic++) {
		printf("Chunk: %i\n", ic);
		indFileFirst = ic * chunkSize;
		indFileLast = indFileFirst + chunkSize - 1;
		if (indFileLast >= nFiles) indFileLast = nFiles - 1;
		nF = indFileLast - indFileFirst + 1;

		sprintf(outputFile, "%s/temp_%i.bin", outputTempFolder, ic);
		fo = fopen(outputFile, "wb");
		if (fo == nullptr) {
			printf("Error: cannot write to \"%s\"\n", outputFile);
			return -6;
		}

		if (readRef() != 0) return -7;
		if (n8 >= 8) {
			if (processReads() != 0) return -8;
		}else if (n8 >= 4) {
			if (processReads4() != 0) return -8;
		}else if (n8 >= 2) {
			if (processReads2() != 0) return -8;
		}else {
			if (processReads1() != 0) return -8;
		}

		fclose(fo); fo = nullptr;
	}
	
	
	if (mergeBlocks() != 0) return -9;
	//if (printOutputData() != 0) return -9;

	return 0;
}


int searchClass::printOutputData() {
	fi = fopen(outputResFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open file \"%s\".\n", outputResFile);
		return -1;
	}

	sprintf(outputFile, "%s/out.txt", outputTempFolder);
	fo = fopen(outputFile, "w");
	if (fo == nullptr) {
		printf("Error: cannot open file \"%s\".\n", outputFile);
		return -2;
	}

	int m2, mtot, ind, mloc;

	for (int k = indStart; k <= indEnd; k++) {
		mtot = 0;

		fread(&ind, sizeof(int), 1, fi);
		fread(&mloc, sizeof(int), 1, fi);
		if (mloc == 0) {
			fprintf(fo, "%d\t%d\n", ind, mloc);
			continue;
		}
		fread(&m2, sizeof(int), 1, fi);
		fprintf(fo, "%d\t%d\t%d\n", ind, mloc, m2);
		fread(vmTemp, sizeof(int), m2, fi);
		/*for (int i = 0; i < m2; i++) {
			fprintf(fo, "\t%u\n", vmTemp[i]);
		}
		fprintf(fo, "\n");*/
	}

	fclose(fo); fo = nullptr;

	fclose(fi); fi = nullptr;

	return 0;
}

int searchClass::mergeBlocks() {
	int m1, m2, m3, myM, mtot, ind, mloc;

	pfi = (FILE**)malloc(sizeof(FILE*) * nChunks);

	for (int i = 0; i < nChunks; i++) {
		sprintf(inputFile, "%s/temp_%i.bin", outputTempFolder, i);
		pfi[i] = fopen(inputFile, "rb");
		if (pfi[i] == nullptr) {
			printf("Error: cannot open file \"%s\".\n", inputFile);
			return -1;
		}
	}

	fo = fopen(outputResFile, "wb");
	if (fo == nullptr) {
		printf("Error: cannot open file \"%s\".\n", outputResFile);
		return -2;
	}

	for (int k = indStart; k <= indEnd; k++) {
		mtot = 0;
		m1 = 0;
		myM = 0;

		for (int i = 0; i < nChunks; i++) {
			fread(&ind, sizeof(int), 1, pfi[i]);
			fread(&mloc, sizeof(int), 1, pfi[i]);

			if (mloc == 0) continue;
			mtot += mloc;
			fread(&m2, sizeof(int), 1, pfi[i]);
			fread(vmTemp, sizeof(int), m2, pfi[i]);
			mergeVectors(pvm[myM], vmTemp, pvm[1 - myM], m1, m2, &m3);
			m1 = m3;
			myM = 1 - myM;
		}
		fwrite(&ind, sizeof(int), 1, fo);
		fwrite(&mtot, sizeof(int), 1, fo);
		if (mtot == 0)continue;
		fwrite(&m1, sizeof(int), 1, fo);
		fwrite(pvm[myM], sizeof(int), m1, fo);
	}

	fclose(fo); fo = nullptr;
	for (int i = 0; i < nChunks; i++) {
		if (pfi[i] != nullptr) { fclose(pfi[i]); pfi[i] = nullptr; }
	}
	free(pfi); pfi = nullptr;
	return 0;
}

int searchClass::readInputParameters(int n, char** v) {
	if (n != 10) {
		printf("Error: wrong number of input parameters.\n");
		printCommandLine();
		return -1;
	}

	sprintf(inputPrefix, v[1]);
	sprintf(inputRefFolder, v[2]);
	sprintf(inputReadFile, v[3]);
	sprintf(outputTempFolder, v[4]);
	sprintf(outputResFile, v[5]);
	sprintf(strSeed, v[6]);
	indStart = atoi(v[7]);
	indEnd = atoi(v[8]);
	nChunks = atoi(v[9]);

	return 0;
}

int searchClass::printCommandLine() {
	printf("Usage: input_binary_file output_text_file\n\n");
	printf("\t1) Reference folder + prefix (reference)\n");
	printf("\t2) Input folder (reference)\n");
	printf("\t3) Input read binary file\n");
	printf("\t4) Temp folder\n");
	printf("\t5) Output file\n");
	printf("\t6) Seed\n");
	printf("\t7) First row to process (the first row = 0)\n");
	printf("\t8) Last row to process\n");
	printf("\t9) Number of chunks\n\n");
	return 0;
}

int searchClass::checkInputParameters() {
	printf("Reference folder + prefix: \"%s\".\n", inputPrefix);
	printf("Input reference folder: \"%s\".\n", inputRefFolder);
	printf("Input read file: \"%s\".\n", inputReadFile);
	printf("Temp folder: \"%s\".\n", outputTempFolder);
	printf("Output file: \"%s\".\n", outputResFile);
	printf("Seed: \"%s\".\n", strSeed);

	nL = (int)strlen(strSeed);
	printf("Length: %i\n", nL);
	if (nL < 1) {
		printf("Error: wrong length.\n");
		return -2;
	}
	if (strSeed[0] != '1') {
		printf("Error: seed should start with 1.\n");
		return -3;
	}
	if (strSeed[nL - 1] != '1') {
		printf("Error: seed should end with 1.\n");
		return -4;
	}

	printf("Row first: %i\n", indStart);
	if (indStart < 0) {
		indStart = 0;
		printf("Row first (corrected): %i\n", indStart);
	}
	printf("Row last: %i\n", indEnd);
	if (indEnd < indStart) {
		indEnd = indStart;
		printf("Row last (corrected): %i\n", indEnd);
	}
	printf("Number of chunks: %i\n", nChunks);
	if (nChunks < 1) {
		printf("Error: wrong number.\n");
		return -2;
	}

	return 0;
}

int searchClass::findParameters() {
	pvm[0] = (unsigned int*)malloc(sizeof(int) * M);
	pvm[1] = (unsigned int*)malloc(sizeof(int) * M);
	vmTemp = (unsigned int*)malloc(sizeof(int) * M);

	if (vmTemp == nullptr || pvm[0] == nullptr || pvm[1] == nullptr) {
		printf("Error: cannot allocate memory (merge).\n");
		return -1;
	}

	nW = 0;
	for (int i = 0; i < nL; i++) {
		if (strSeed[i] == '1') nW++;
	}

	ipref = (2 * nW) % 8;
	while (ipref + 8 <= 16) {
		ipref += 8;
	}
	n8 = (2 * nW - ipref) / 8;

	printf("Weight: %i\n", nW);
	indOne = (int*)malloc(sizeof(int) * nW);
	nW = 0;
	for (int i = 0; i < nL; i++) {
		if (strSeed[i] == '1') {
			indOne[nW] = i;
			nW++;
		}
	}

	printf("Prefix size: %i\n", ipref);
	printf("Number of chars: %i\n", n8);

	ndigits = ipref / 4;
	if (ipref % 4 > 0)ndigits++;

	nFiles = 1 << ipref;
	sLen = 4 + n8;

	chunkSize = nFiles / nChunks;
	if (nFiles % nChunks > 0)chunkSize++;

	iRes = (unsigned int*)malloc(sizeof(int) * nW);

	if (n8 >= 8) {
		nQ = n8 / 8;
		if (n8 % 8 > 0)nQ++;
		for (int i = 0; i < nQ; i++) {
			iqStart[i] = n8 - 8 * (i + 1);
			if (iqStart[i] < 0) iqStart[i] = 0;
		}
	}
	else if (n8 >= 4) {
		nQ = n8 / 4;
		if (n8 % 4 > 0)nQ++;
		for (int i = 0; i < nQ; i++) {
			iqStart[i] = n8 - 4 * (i + 1);
			if (iqStart[i] < 0) iqStart[i] = 0;
		}
	}else if (n8 >= 2) {
		nQ = n8 / 2;
		if (n8 % 2 > 0)nQ++;
		for (int i = 0; i < nQ; i++) {
			iqStart[i] = n8 - 2 * (i + 1);
			if (iqStart[i] < 0) iqStart[i] = 0;
		}
	}else{
		nQ = 1;
		iqStart[0]  = 0;
	}

	return 0;
}

int searchClass::readReads() {
	fi = fopen(inputReadFile, "rb");
	if (fi == nullptr) {
		printf("Error: cannot open file \"%s\"\n", inputReadFile);
		return -1;
	}
	fread(&nLen, sizeof(int), 1, fi);
	printf("Length of reads: %i\n", nLen);

	n4 = nLen / 4;
	nrem = nLen % 4;
	if (nrem > 0) n4++;
	if (nrem == 0) nrem = 4;

	_fseeki64(fi, 0, SEEK_END);
	nRecords = (int)((_ftelli64(fi) - 4) / ((long long)n4));

	printf("Number of records: %i\n", nRecords);
	if (indStart >= nRecords) {
		printf("Error: not enough records\n");
		return -2;
	}

	if (indEnd >= nRecords) {
		indEnd = nRecords - 1;
		printf("Row last (corrected): %i\n", indEnd);
	}

	nRows = indEnd - indStart + 1;

	_fseeki64(fi, 4 + n4 * indStart, SEEK_SET);

	vr = (unsigned char*)malloc(sizeof(char) * n4 * nRows);
	fread(vr, sizeof(char), n4 * nRows, fi);

	fclose(fi); fi = nullptr;

	iValue = (unsigned int*)malloc(sizeof(int) * nLen);

	return 0;
}

int main(int argc, char** argv) {
	int ires;
	searchClass* sc;

	sc = new searchClass();
	ires = sc->startProcessing(argc, argv);
	delete sc;

	return ires;
}
	

	