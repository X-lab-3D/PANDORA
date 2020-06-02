//Aligns PDB 2 to PDB 1 and renumbers accordingly, provided that all residues of PDB 2 are present in PDB 1. If they are not present, a warning is generated, and the residue is named X 1, 2, ...

#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>
using namespace std;

inline char getresidue(char code[3]) {
    if (!strcmp(code, "ALA")) return 'A';
    if (!strcmp(code, "CYS")) return 'C';
    if (!strcmp(code, "ASP")) return 'D';
    if (!strcmp(code, "GLU")) return 'E';
    if (!strcmp(code, "PHE")) return 'F';
    if (!strcmp(code, "GLY")) return 'G';
    if (!strcmp(code, "HIS")) return 'H';
    if (!strcmp(code, "ILE")) return 'I';
    if (!strcmp(code, "LYS")) return 'K';
    if (!strcmp(code, "LEU")) return 'L';
    if (!strcmp(code, "MET")) return 'M';
    if (!strcmp(code, "ASN")) return 'N';
    if (!strcmp(code, "PRO")) return 'P';
    if (!strcmp(code, "GLN")) return 'Q';
    if (!strcmp(code, "ARG")) return 'R';
    if (!strcmp(code, "SER")) return 'S';
    if (!strcmp(code, "THR")) return 'T';
    if (!strcmp(code, "VAL")) return 'V';
    if (!strcmp(code, "TRP")) return 'W';
    if (!strcmp(code, "TYR")) return 'Y';
    return 'X';
}

int extractpdb(const char *pdb, char chainid, char *&seq, int *&num) {
  bool has_chain = 0;
  FILE *fil = fopen(pdb, "r");
  if (fil == NULL) {fprintf(stderr, "PDB file %s does not exist\n", pdb); exit(1); }
  int num0[10000];
  char seq0[10000];
  int counter = 0;
  int curr_resnr = -1000;

  while (!feof(fil)) {
    char buf[100];
    char code[6];
    char name[3];
    fgets(buf, 100, fil);
    sscanf(buf, "%s %*s %*s %s", code, name);
    if (!strncmp(code,"ATOM", 4)) {
      char currchainid = buf[21];
      if (currchainid != chainid) continue;
      has_chain = 1;
      int resnr = atoi(buf+22);
      if (resnr != curr_resnr) {
        seq0[counter] = getresidue(name);
	num0[counter] = resnr;
	curr_resnr = resnr;
	counter++;
      }
    }
  }
  fclose(fil);
  if (!counter) {fprintf(stderr, "%s is not a PDB file\n", pdb); exit(1); }
  if (!has_chain) {fprintf(stderr, "%s has no chainid %c\n", chainid); exit(1);}
  seq = new char[counter+1];
  memcpy(seq, seq0, counter);
  seq[counter] = 0;
  num = new int[counter];
  memcpy(num, num0, counter *sizeof(int));
  return counter;
}

int match = 40;
int mismtch = -20;
int gap_open = -100;
int gap_cont = -2;


struct AlignmentPos {
  int score;
  int direc; //0 = left, 1 = up, 2 = left-up
};

void initialize(AlignmentPos **aln, int seqlen1, int seqlen2){
  int n;
  for (n = 0; n < seqlen1 + 1; n++) {
    for (int nn = 0; nn < seqlen2 + 1; nn++) aln[nn][n].score = -1000000;
  }
  aln[0][0].score = 0;
  aln[0][0].direc = 0;
}

bool feedforward(AlignmentPos **aln, int seqlen1, int seqlen2, int row) {
  int n;
  bool change = 0;
  for (n = 0; n < seqlen1; n++) {
    int difscore = gap_open;
    if (n > 0 && aln[row][n].direc == 0) difscore = gap_cont;
    int newscore = aln[row][n].score + difscore;
    if (newscore > aln[row][n+1].score) {
      change = 1;
      aln[row][n+1].score = newscore;
      aln[row][n+1].direc = 0;
    } 
  }
  return change;
}

bool feeddown(AlignmentPos **aln, int seqlen1, int seqlen2,int column) {
  int n;
  bool change = 0;
  for (n = 0; n < seqlen2; n++) {
    int difscore = gap_open;
    if (n > 0 && aln[n][column].direc == 1) difscore = gap_cont;
    int newscore = aln[n][column].score + difscore;
    if (newscore > aln[n+1][column].score) {
      change = 1;
      aln[n+1][column].score = newscore;
      aln[n+1][column].direc = 1;
    } 
  }
  return change;
}
bool feedforwarddown(AlignmentPos **aln, char *seq1, char *seq2, int seqlen1, int seqlen2, bool startleft, int diagonal) {
  bool change = 0;
  int start1 = 0, start2 = 0;
  if (startleft) start1 = diagonal; else start2 = diagonal;
  for (int n1 = start1, n2 = start2; n1 < seqlen2 && n2 < seqlen1; n1++, n2++) {
    int difscore = mismtch;
    if (seq2[n1] == seq1[n2]) difscore = match;
    int newscore = aln[n1][n2].score + difscore;
    if (newscore > aln[n1+1][n2+1].score) {
      change = 1;
      aln[n1+1][n2+1].score = newscore;
      aln[n1+1][n2+1].direc = 2;
    }
  }
  return change;
}

map<int, int> getconv(AlignmentPos **aln, char *seq1, char *seq2, int *num1, int *num2, int seqlen1, int seqlen2) {
  
  map<int, int> conv;
  
  char pr[10000][30];
  char er[1000][100];
  int ersize = 0;
  int pos1 = seqlen1, pos2 = seqlen2;
  int counter = 0;
  while (pos1 >= 0 && pos2 >= 0) {
    switch(aln[pos2][pos1].direc) {
      case 0:
        pos1--;
        break;
      case 1:
	sprintf(er[ersize], "Warning: Residue %c%d not present in first PDB\n", seq2[pos2-1], num2[pos2-1]);
	ersize++;
	pos2--;
	break;	
      case 2:
	if (seq1[pos1-1] != seq2[pos2-1]) {
	  sprintf(er[ersize], "Warning: Mismatch %c%d => %c\n", seq2[pos2-1], num2[pos2-1], seq1[pos1-1]);
	  ersize++;
	}
	conv[num2[pos2-1]] = num1[pos1-1];
        pos1--;
	pos2--;	
	break;
    }
  }
  int n;
  for (n = 0; n < ersize; n++) {
    fprintf(stderr, "%s", er[ersize-n-1]);
  }
  return conv;
}

void doconv(map<int,int> &conv, const char *pdb, char chainid1, char chainid2) {
  FILE *fil = fopen(pdb, "r");
  if (fil == NULL) {fprintf(stderr, "PDB file %s does not exist\n", pdb); exit(1); }

  int curr_resnr = -1000;
  int xcount = 0;
  while (!feof(fil)) {
    char buf[100];
    char pbuf[100];
    if (!(fgets(buf, 100, fil))) continue;
    if (buf[21] != chainid2) continue;
    if (strncmp(buf, "ATOM", 4)) {
      printf("%s", buf);
      continue;
    }
    int resnr = atoi(buf+22);
    strcpy(pbuf, buf);
    if (conv.count(resnr)) {
      pbuf[21] = chainid1;
      sprintf(&pbuf[22], "%4d", conv[resnr]);
      memcpy(&pbuf[26], &buf[26], 73);    
    }
    else {
      if (resnr != curr_resnr) xcount++;
      pbuf[21] = 'X';
      sprintf(&pbuf[22],"%4d", xcount);
      memcpy(&pbuf[26], &buf[26], 73);
    }
    curr_resnr = resnr;
    printf("%s", pbuf);
  }
  fclose(fil);
}

int main(int argc, char *argv[]) {
  int n;

  if (argc < 5) {
    cerr << "Usage: pdb-pdbalign <pdb file 1> <chainid 1> <pdb file 2> <chainid2> " << endl;
    return 1;
  }
  char *seq1; int *num1;
  char *seq2; int *num2;
  
  char chainid1 = argv[2][0];
  char chainid2 = argv[4][0];
  int seqlen1 = extractpdb(argv[1], chainid1, seq1, num1);
  int seqlen2 = extractpdb(argv[3], chainid2, seq2, num2);
  AlignmentPos **aln = new AlignmentPos*[seqlen2+1]; //seqlen2+1 rows, ranging from 0 to seqlen
  for (n = 0; n < seqlen2 + 1; n++) {
    aln[n] = new AlignmentPos[seqlen1+1];
  }
  
  initialize(aln, seqlen1, seqlen2);
  bool change = 1;
  int counter = 0;
  while (change) {
    change = 0;
    for (n = 0; n < seqlen2; n++) {
      if (feedforwarddown(aln, seq1, seq2, seqlen1, seqlen2, 1, n)) change = 1;
    }
    for (n = 1; n < seqlen1; n++) {
      if (feedforwarddown(aln, seq1, seq2, seqlen1, seqlen2, 0, n)) change = 1;
    }    
    for (n = 0; n < seqlen2+1; n++) {
      if (feedforward(aln,seqlen1, seqlen2, n)) change = 1;
    }    
    for (n = 0; n < seqlen1+1; n++) {
      if (feeddown(aln,seqlen1, seqlen2, n)) change = 1;
    }        
    counter++;
  }
  
  char *aln1; char *aln2;
  map<int, int> conv = getconv(aln, seq1, seq2, num1, num2, seqlen1, seqlen2);
  doconv(conv, argv[3], chainid1, chainid2);  

  return 0;
}
