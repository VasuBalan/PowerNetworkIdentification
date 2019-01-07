/***********************************************************************************
*   Combining graph.c and pso.c for Reconfiguration of Distribution System    	   *
*   written by Vasudevan.B (vasubdevan@yahoo.com) on  05.04.2016 @ PS Lab  IIT KGP *
************************************************************************************/ 

#include <stdio.h>
#include <stdlib.h>
#define Miter 5

// For PSO 
#define nt 2
#define np 10

// File Details
typedef struct {
  char input[20];
  char output[20];
}FileDetails;

// System Info Details
typedef struct {
  int nb;
  int nl;
  int ntr;
  int nv;
  int ns;
  int nsvs;
  int nphs;
  int nsb;
  double pbase;
  int kmax;
  double epsln;
  int lfdc;
  
  int nfd;
}SysInfo;

// Bus data details
typedef struct {
  int code;
  char name[7];
  int zone;
  double Basekv;
  double RealPg;
  double ReaQg;
  double RealL;
  double ReacL;
  double v;
  double DegAngle;
  int BType;
  int LType;
  int BStat;
  int LStat;
}BusData;

// Feeder data details
typedef struct {
  int sno;
  char sym[2];
  int fb;
  int tb;
  int CktNo;
  double Rpu;
  double Xpu;
  double Bpu;
  double MVA;
  double fbXpu;
  double tbXpu;
  int FStat;
}FeederData;

// Transformer details
typedef struct {
  int sno;
  char sym[2];
  int fb;
  int tb;
  int CktNo;
  double Rpu;
  double Xpu;
  double MVA;
  double IncTap;
  double MinTapPos; // Have to check whether double or integer required.
  double MaxTapPos;
  double CtTapPos;
  double Tapratio;
  double PhShiAng;
  int TranType;
  int TranStat;
}TransfData;

// PV Bus details
typedef struct {
  int sno;
  int busno;
  double MinPG;
  double MaxPG;
  double MinQG;
  double MaxQG;
  double Vpu;
  double MinpuV;
  double MaxpuV;
}PVBus;

// Slack bus details
typedef struct {
  int sno;
  int busno;
  double MinPG;
  double MaxPG;
  double MinQG;
  double MaxQG;
  double Vpu;  
}SlackBus;

// Shunt bus details
typedef struct {
  int sno;
  int busno;
  double G;
  double B;
  int ShStat;
}ShuntBus;

// ----------------------------- Input File Reading Functions -----------------

FileDetails ReadFileDetails(FileDetails fpd){
  // Function is to read the input and output file from offline.txt
  FILE *fp1 = fopen("Offline.txt","r");
  if(fp1 == NULL) {
    printf("\nOffline.txt doesn't exist\n");
    exit(1);
  } else {
    fscanf(fp1, "%s", fpd.input);
    fscanf(fp1, "%s", fpd.output);
  }
  fclose(fp1);
  return fpd;
}

FILE *FilePointerToRead(FILE *fp, char *fname){
  // To read the file pointer name  
  fp = fopen(fname,"r");
  if(fp == NULL){
    printf("\nError in file openeing\n");
    exit(1);
  }
  return fp;
}

SysInfo ReadEntireSystemData(SysInfo sys, FILE *fp) {
  // Reading entire system data  
  fscanf(fp,"\n%d %d %d %d %d %d %d %d %lf %d %lf %d", &sys.nb, &sys.nl, &sys.ntr, &sys.nv, &sys.ns,
	 &sys.nsvs, &sys.nphs, &sys.nsb, &sys.pbase, &sys.kmax, &sys.epsln, &sys.lfdc);
  sys.nfd = sys.nl - sys.ntr;
  return sys;
}

void PrintEntireSystemInfo(SysInfo sys){
  printf("\n\nEntire System Detail is as follows\n\n");
  printf("\n%d %d %d %d %d %d %d %d %lf %d %lf %d", sys.nb, sys.nl, sys.ntr, sys.nv, sys.ns,
	 sys.nsvs, sys.nphs, sys.nsb, sys.pbase, sys.kmax, sys.epsln, sys.lfdc);
  printf("\nNo of feeders = %d\n", sys.nfd);
}

BusData *initializeBusData(BusData *bus, int nb){
  // Initializing the bus pointer
  int i, j;
  for( i = 0; i < nb; i++) {
    bus[i].code= 0;
    for(j = 0; j < 6; j++) {
      bus[i].name[j] = 0;
    }
    bus[i].zone = 0;
    bus[i].Basekv = 0;
    bus[i].RealPg = 0;
    bus[i].ReaQg = 0;
    bus[i].RealL = 0;
    bus[i].ReacL = 0;
    bus[i].v = 0;
    bus[i].DegAngle = 0;
    bus[i].BType = 0;
    bus[i].LType = 0;
    bus[i].BStat = 0;
    bus[i].LStat = 0;
  }
  return bus;
}

void printBusData(BusData *bus, int n) {
  int i;
  printf("\n\nBus data is as follows\n\n");
  for (i = 0; i < n; i++) {
  printf("\n%d %s %d %lf %lf %lf %lf %lf %lf %lf %d %d %d %d",bus[i].code, bus[i].name, bus[i].zone, bus[i].Basekv, bus[i].RealPg, bus[i].ReaQg, 
	 bus[i].RealL, bus[i].ReacL, bus[i].v, bus[i].DegAngle, bus[i].BType, bus[i].LType, bus[i].BStat, bus[i].LStat);
  }
  printf("\n");
}

BusData *ReadBusData(BusData *bus, FILE *fp, int n){
  // Reading bus data
  int i;
  bus = (BusData*)malloc(n * sizeof(BusData));
  if(bus == NULL) {
    printf("\nError in memory allocation of bus pointer\n");
    exit(1);
  }
  
  bus = initializeBusData(bus, n);
   
  for (i = 0; i < n; i++) {
     fscanf(fp, "%d %s %d %lf %lf %lf %lf %lf %lf %lf %d %d %d %d", &bus[i].code, bus[i].name, &bus[i].zone, &bus[i].Basekv, &bus[i].RealPg, &bus[i].ReaQg, 
	 &bus[i].RealL, &bus[i].ReacL, &bus[i].v, &bus[i].DegAngle, &bus[i].BType, &bus[i].LType, &bus[i].BStat, &bus[i].LStat);
  }
  return bus;
}

FeederData *initializeFeederData(FeederData *fd, int n) {
  // Initializing the feeder pointer
  int i, j;
  for(i = 0; i < n; i++) {
    fd[i].sno = 0;
    for(j = 0; j < 2; j++) {
    fd[i].sym[j] = 0;
    }
    fd[i].fb = 0;
    fd[i].tb = 0;
    fd[i].CktNo = 0;
    fd[i].Rpu = 0.0;
    fd[i].Xpu = 0.0;
    fd[i].MVA = 0.0;
    fd[i].fbXpu = 0.0;
    fd[i].tbXpu = 0.0;
    fd[i].FStat = 0;
  }
  return fd;
}

FeederData *ReadFeederData(FeederData *fd, FILE *fp, int n) {
  // Reading feeder data
  int i;
  fd = (FeederData *)malloc(n * sizeof(FeederData));
  if(fd == NULL) {
    printf("\nError in memory allocation of Feeder pointer\n");
    exit(1);
  }
  
  fd = initializeFeederData(fd, n);
  for (i = 0; i < n; i++) {
    fscanf(fp, "%d %s %d %d %d %lf %lf %lf %lf %lf %lf %d",&fd[i].sno, fd[i].sym, &fd[i].fb, &fd[i].tb, &fd[i].CktNo, &fd[i].Rpu, &fd[i].Xpu, &fd[i].Bpu, &fd[i].MVA, 
	   &fd[i].fbXpu, &fd[i].tbXpu, &fd[i].FStat);
  }
  return fd;
}

void printFeederData(FeederData *fd, int n) {
  int i;
  printf("\n\nFeeder Data is as follows\n\n");
  for (i = 0; i < n; i++) {
     printf("\n%d %s %d %d %d %lf %lf %lf %lf %lf %lf %d",fd[i].sno, fd[i].sym, fd[i].fb, fd[i].tb, fd[i].CktNo, fd[i].Rpu, fd[i].Xpu, fd[i].Bpu, fd[i].MVA, 
	      fd[i].fbXpu, fd[i].tbXpu, fd[i].FStat);
  }
  printf("\n");
}

TransfData *initializeTransformerData(TransfData *td, int n) {
  // Initialize the transformer data
  int i, j;
  for(i = 0; i < n; i++) {
    td[i].sno = 0;
    for(j = 0; j < 2; j++) {
    td[i].sym[j] = 0;
    }
    td[i].fb = 0;
    td[i].tb = 0;
    td[i].CktNo = 0;
    td[i].Rpu = 0.0;
    td[i].Xpu = 0.0;
    td[i].MVA = 0.0;
    td[i].IncTap = 0.0;
    td[i].MinTapPos = 0.0;
    td[i].MaxTapPos = 0.0;
    td[i].CtTapPos = 0.0;
    td[i].Tapratio = 0.0;
    td[i].PhShiAng = 0.0;
    td[i].TranType = 0;
    td[i].TranStat = 0;
  }
  return td;
}

TransfData *ReadTransfData(TransfData *td, FILE *fp, int n) {
  // Reading transformer data
  int i;
  td = (TransfData *)malloc(n * sizeof(TransfData));
  if(td == NULL) {
    printf("\nError in memory allocation of Feeder pointer\n");
    exit(1);
  }
  
  td = initializeTransformerData(td, n);
  
  for (i = 0; i < n; i++) {
    fscanf(fp, "%d %s %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", &td[i].sno, td[i].sym, &td[i].fb, &td[i].tb, &td[i].CktNo, &td[i].Rpu, &td[i].Xpu, 
	   &td[i].MVA, &td[i].IncTap, &td[i].MinTapPos, &td[i].MaxTapPos, &td[i].CtTapPos, &td[i].Tapratio, &td[i].PhShiAng, &td[i].TranType, &td[i].TranStat);
  }
  
  return td;
}

void PrintTransformerData(TransfData *td, int n) {
  int i;
  if(n == 0) 
    printf("\nNo transformer in the given system\n");
  else {
    printf("\n\nTransformer Data is as follows\n\n");
    for (i = 0; i < n; i++) {
      printf("\n%d %s %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %d %d", td[i].sno, td[i].sym, td[i].fb, td[i].tb, td[i].CktNo, td[i].Rpu, td[i].Xpu, 
	   td[i].MVA, td[i].IncTap, td[i].MinTapPos, td[i].MaxTapPos, td[i].CtTapPos, td[i].Tapratio, td[i].PhShiAng, td[i].TranType, td[i].TranStat);
    }
    printf("\n");
  }
}

PVBus *initializePVBusData(PVBus *pv, int n) {
  // Initialize pv bus data
  int i;
  for(i = 0; i < n; i++) {
    pv[i].sno = 0;
    pv[i].busno = 0;
    pv[i].MinPG = 0.0;
    pv[i].MaxPG = 0.0;
    pv[i].MinQG = 0.0;
    pv[i].MaxQG = 0.0;
    pv[i].Vpu = 0.0;
    pv[i].MinpuV = 0.0;
    pv[i].MaxpuV = 0.0;
  }
  return pv;
}

PVBus *ReadPVBusData(PVBus *pv, FILE *fp, int n) {
  // Reading PV bus data
  int i;
  pv = (PVBus *)malloc(n * sizeof(PVBus));
  if(pv == NULL) {
    printf("\nError in memory allocation of Feeder pointer\n");
    exit(1);
  }
  
  pv = initializePVBusData(pv, n);
  
  for(i = 0; i < n; i++) {
    fscanf(fp, "%d %d %lf %lf %lf %lf %lf %lf %lf", &pv[i].sno, &pv[i].busno, &pv[i].MinPG, &pv[i].MaxPG, &pv[i].MinQG, 
	   &pv[i].MaxQG, &pv[i].Vpu, &pv[i].MinpuV, &pv[i].MaxpuV);
  }
  return pv;
}

SlackBus initializeSlackBusData(SlackBus slack) {
  // Initialize slack bus variables
  slack.sno = 0;
  slack.busno = 0;
  slack.MinPG = 0;
  slack.MaxPG = 0;
  slack.MinQG = 0;
  slack.MaxQG = 0;
  return slack;
}

void printPVBusData(PVBus *pv, int n) {
  int i;
  if(n == 0)
    printf("\nNo PV bus in the given system\n");
  else {
    printf("\n\nPV Bus data is as follows\n");
    for(i = 0; i < n; i++) {
      printf("\n%d %d %lf %lf %lf %lf %lf %lf %lf", pv[i].sno, pv[i].busno, pv[i].MinPG, pv[i].MaxPG, pv[i].MinQG, 
	   pv[i].MaxQG, pv[i].Vpu, pv[i].MinpuV, pv[i].MaxpuV);
    }
    printf("\n");
  }
}

SlackBus ReadSlackBusData(SlackBus slack, FILE *fp) {
  // Read slack bus data
  slack = initializeSlackBusData(slack);
  
  fscanf(fp, "%d %d %lf %lf %lf %lf %lf", &slack.sno, &slack.busno, &slack.MinPG, &slack.MaxPG, &slack.MinQG, &slack.MaxQG, &slack.Vpu);
  return slack;
}

void printSlackBusData(SlackBus slack) {
  printf("\n\nSlack bus data is as follows\n\n");
  printf("\n%d %d %lf %lf %lf %lf %lf\n", slack.sno, slack.busno, slack.MinPG, slack.MaxPG, slack.MinQG, slack.MaxQG, slack.Vpu);
}

ShuntBus *initializeShuntBusData(ShuntBus *shunt, int n) {
  // Initialize shunt bus data
  int i;
  for(i = 0; i < n; i++) {
    shunt[i].sno = 0;
    shunt[i].busno = 0;
    shunt[i].B = 0;
    shunt[i].G = 0;
    shunt[i].ShStat = 0;
  }
  return shunt;
}

ShuntBus *ReadShuntBusData(ShuntBus *shunt, FILE *fp, int n) {
  // Reading shunt bus data
  int i;
  
  shunt = (ShuntBus *)malloc(n * sizeof(ShuntBus));
  if(shunt == NULL) {
    printf("\nError in memory allocation of shunt bus pointer\n");
    exit(1);
  }
  
  shunt = initializeShuntBusData(shunt, n);
  for(i = 0; i < n; i++) {
    fscanf(fp, "%d %d %lf %lf %d", &shunt[i].sno, &shunt[i].busno, &shunt[i].G, &shunt[i].B, &shunt[i].ShStat);
  }
  
  return shunt;
}

void printShuntBusData(ShuntBus *shunt, int n) {
  int i;
  if(n == 0)
    printf("\nNo Shunt data\n\n");
  else {
  printf("\nShunt data is as follows\n");
    for(i = 0; i < n; i++) {
      printf("\n%d %d %lf %lf %d", shunt[i].sno, shunt[i].busno, shunt[i].G, shunt[i].B, shunt[i].ShStat);
    }
    printf("\n");
  }
}

// ----------------------------- Matrix Operation Starts -----------------
int **initializeAdjMatrixN(int **A, int n) {
    int i,j;
    for(i = 0; i < n; i++) {
	for(j = 0; j < n; j++) {
	    A[i][j] = 0;
	}
    }
    return A;
}

int **createAdjMatrixN(int **A, int n) {
    int i, j;
    A = (int **) malloc(n * sizeof(int *) );
    if(A == NULL) {
      printf("\nOut of memory\n");
      exit(1);
    }

    for(i = 0; i < n; i++) {
      	A[i] = (int *) malloc ( n * sizeof(int) );
	if(A[i] == NULL) {
	    printf("\nMemory error\n");
	    exit(1);
	}
    }
    A = initializeAdjMatrixN(A, n);
    return A;
}

void printAdjMatrixN(int **A, int n) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	int i, j;
	printf("\nPrinting the Square Matrix\n-----------------\n");
	for(i = 0; i < n; i++) {
	    for(j = 0; j < n; j++) {
		printf("%d\t ", A[i][j]);
	    }
	    printf("\n");
	}
	printf("-----------------\n\n");
    }
}

int **populateMatrixN(int **A, SysInfo sys, FeederData *fd) {
    int i, j, l;
    // Creating the adjacency matrix for the given network
    for(i = 0; i < sys.nb; i++){
	for(j = 0; j < sys.nb; j++) {
	      if( i == j)
		  A[i][j] = 0;
	      else {
		for(l = 0; l < sys.nl; l++) {
		    if(i+1 == fd[l].fb && j+1 == fd[l].tb && fd[l].FStat == 1) {
		      A[i][j] = 1; }
		    if(i+1==fd[l].tb && j+1 == fd[l].fb && fd[l].FStat ==1) {
		      A[i][j] = 1; }
		  }
	      }
	  }
	} 
   return A;
}

int *initializeArrayValuesNint(int *A, int n) {
    int i;
    for(i = 0; i < n; i++) {
      A[i] = 0;
    }
    return A;
}

int *createArrayValuesNint(int *A, int n) {
  int *array;
  array = (int *) malloc (n * sizeof(int));
  initializeArrayValuesNint(array, n);
  return array;
}

void printArrayValuesNint(int *A, int n) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	printf("\nArray Element is as follows\n");
	int i;
	for(i = 0; i < n; i++) {
	  printf("%d\n", A[i]);
	}
    }
}

void FreeMemory(int **A, int n) {
    int i;
    for (i = 0; i < n; i++) free(A[i]);
    free(A);
}

void printArrayValuesN(double *A, int n) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	printf("\nArray Element is as follows\n");
	int i;
	for(i = 0; i < n; i++) {
	  printf("%lf\n", A[i]);
	}
    }
}

double *initializeArrayValuesN(double *A, int n) {
    int i;
    for(i = 0; i < n; i++) {
      A[i] = 0.0;
    }
    return A;
}

double *createArrayValuesN(double *A, int n ) {
  double *array;
  array = (double *) malloc (n * sizeof(double));
  initializeArrayValuesN(array, n);
  return array;
}

double *getGenerationData(double *generation, char *fname, int n) {
    int i;
    
    generation = createArrayValuesN(generation, n);
    
    //Reading pdt file 
    FILE *fpg = fopen(fname, "r");
    if (fpg == NULL) {
	printf("The file %s was not opened\n", fname);
	exit(1);
    } else {
	for (i = 0; i < n; i++) {
	    fscanf(fpg, "%lf",&generation[i] );
	}
    } fclose(fpg);
//     printArrayValuesN(generation, n);
    return generation;
}

double *getLoadWeightage(double *weightage, char *fname, int n) {
    int i;
    
    weightage = createArrayValuesN(weightage, n);
    
    //Reading pdt file 
    FILE *fpw = fopen(fname, "r");
    if (fpw == NULL) {
	printf("The file %s was not opened\n", fname);
	exit(1);
    } else {
	for (i = 0; i < n; i++) {
	    fscanf(fpw, "%lf",&weightage[i] );
	}
    } fclose(fpw);
//     printArrayValuesN(generation, n);
    return weightage; 
}

int **initializeAdjMatrixNM(int **A, int n, int m) {
    int i,j;
    for(i = 0; i < n; i++) {
	for(j = 0; j < m; j++) {
	    A[i][j] = 0;
	}
    }
    return A;
}

int **createAdjMatrixNM(int **A, int n, int m) {
    int i;
    A = (int **) malloc(n * sizeof(int *) );

    for(i = 0; i < n; i++) {
	A[i] = (int *) malloc ( m * sizeof(int) );
    }   
    A = initializeAdjMatrixNM(A, n, m);
    
    return A;
}

void printAdjMatrixNM(int **A, int n, int m) {
    if (A == NULL) {
	printf("Matrix is NULL\n");
    } else {
	int i, j;
	printf("\n\nPrinting the matrix\n-----------------\n");
	for(i = 0; i < n; i++) {
	    for(j = 0; j < m; j++) {
		printf("%d\t", A[i][j]);
	    }
	    printf("\n");
	}
	printf("-----------------\n\n");
    }
}

double rand_number() {
    int x, high = 100, low = 0;
    double r;
    x = rand() % (high - low + 1) + 0;		//high = 100, low =0
    r = x * 0.01;
    return r;
}

int random_bus(int nb) {
    int x;
    label1:
    x = rand() % (nb+1);
    if (x <= nb && x > 1)
	return x;
    else {
	goto label1;
    }
    return x;
}

int **populateRandMatrix(int **B, int nb, double *ini_vel, int ntie, int NPART) {
    int i, j, k;
    for (i = 0; i < NPART; i++) {
	for (j = 0; j < (2 * ntie); j++) {
	    B[i][j] = random_bus(nb);
	}
	//Condition to check whether both from and to buses are to be same.
	for (k = 0; k < ((2 * ntie) - 1);) {
	    if (B[i][k] == B[i][k + 1]) {
		label2:
		B[i][k + 1] = random_bus(nb);
		if (B[i][k + 1] == B[i][k]) {
		    goto label2;
		}
	    }
	    k = k + 2;
	}
	ini_vel[i] = rand_number();
    }
    return B;
}

void FreeMemoryMN(int **B, int m, int n) {
    int i;
    for (i = 0; i < m; i++) free(B[i]);
    free(B);
}

int *copyRandMatrixtoDArray(int **B, int a, int *D, int b) {
    int i,j;
    for(i=0; i < a; i++) {
      D[i] = B[b][i];
    }
    return D;
}

int **copyAdjacencyMatrix(int **A, int **B, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
	for (j = 0; j < n; j++) {
	    B[i][j] = A[i][j];
	}
    }
    return B;
}
// ------------------------------ Main Program Starts here --------------------

int main(int argc, char *argv[]) {
  
  int **A;
  
  // For island detection Purpose variables
  int **B, *x, *y, *prod, *madd, *LC, *Is1_Node, *Is2_Node, **C, *D;
  double *generation, *weightage, *sirf, *ini_vel;
  
  int iter, k;
  
  // --------------------------- Reading input datas -----------------------
  
  FileDetails fpd = ReadFileDetails(fpd);
  
  FILE *fp = FilePointerToRead(fp,fpd.input);
  
  SysInfo sys = ReadEntireSystemData(sys, fp);
  // PrintEntireSystemInfo(sys);
  
  BusData *bus = ReadBusData(bus, fp, sys.nb);
  // printBusData(bus, sys.nb);
    
  FeederData *fd = ReadFeederData(fd, fp, sys.nfd);
  // printFeederData(fd, sys.nfd);
  
  TransfData *td = ReadTransfData(td, fp, sys.ntr);
  // PrintTransformerData(td, sys.ntr);
  
  PVBus *pv = ReadPVBusData(pv, fp, sys.nv);
  // printPVBusData(pv, sys.nv);
  
  SlackBus slack = ReadSlackBusData(slack, fp);
  // printSlackBusData(slack);
  
  ShuntBus *shunt = ReadShuntBusData(shunt, fp, sys.ns);
  // printShuntBusData(shunt, sys.ns);
  
  fclose(fp);   // Close the file pointer after read entire data.
  
  // --------------------------- Topology checking process  -----------------------
  
  A = createAdjMatrixN(A, sys.nb);   // A - Adjacency Matrix 
  // printAdjMatrixN(A, sys.nb); 
  
  A = populateMatrixN(A, sys, fd);
  // printAdjMatrixN(A, sys.nb);
  
  B = createAdjMatrixN(B, sys.nb);  // B - Temporary Matrix to copy Adjacency matrix
  // printAdjMatrixN(B, sys.nb);
  
  x = createArrayValuesNint(x, sys.nb);
  // printArrayValuesNint(x, sys.nb);
  
  y = createArrayValuesNint(y, sys.nb);
  // printArrayValuesNint(y, nb);
    
  prod = createArrayValuesNint(prod, sys.nb);
  //printArrayValuesNint(prod, nb);
    
  madd = createArrayValuesNint(madd, sys.nb);
  // printArrayValuesNint(madd, nb);
    
  LC = createArrayValuesNint(LC, sys.nb);
  // printArrayValuesNint(LC, nb);
    
  Is1_Node = createArrayValuesNint(Is1_Node, sys.nb);
  // printArrayValuesNint(Is1_Node, nb);
    
  Is2_Node = createArrayValuesNint(Is2_Node, sys.nb);
  // printArrayValuesNint(Is2_Node, nb);
  
  generation = getGenerationData(generation, "gendata15.pdt", sys.nb);
  // printArrayValuesN(generation, sys.nb);
    
  weightage = getLoadWeightage(weightage, "eweight15.pdt", sys.nb);
  // printArrayValuesN(weightage, sys.nb); 
  
  sirf = createArrayValuesN(sirf, sys.nb);
  // printArrayValuesN(sirf, sys.nb);
  
  ini_vel = createArrayValuesN(ini_vel, sys.nb);
  // printArrayValuesN(ini_vel, nb);
  
  // --------------------------- Random Number Generation -----------------------
  C = createAdjMatrixNM(C, np, (2*nt));    // M - represents the number of particle as rows.
  // printAdjMatrixNM(C, np, (2*nt));      // N -represents the number of tie switch pairs as coluuns.
  
  C = populateRandMatrix(C, sys.nb, ini_vel, nt, np);    // Generating random number for this matrix
  // printAdjMatrixNM(C, np, (2*nt));
  
   D = createArrayValuesNint(D, (2*nt));  // Creating a temporary array in order to copy a single row elements from the above matrix. 
   // printArrayValuesNint(D, (2*nt));
  
  iter = 4;
  do {
    printf("\nIteration = %d", iter+1);
    for(k = 0; k < 1; k++) {  // changed np as1 ---> Vasudevan
      D = initializeArrayValuesNint(D, (2*nt));
      D = copyRandMatrixtoDArray(C, (2*nt), D, k);
      printArrayValuesNint(D, 2*nt);
      
      B = initializeAdjMatrixN(B, sys.nb);
            
      B = copyAdjacencyMatrix(A, B, sys.nb); 
      printAdjMatrixN(B, sys.nb);
      
      
      if(nt == 2) {  // Kept two switches in operation at any time 
	B[D[0]-1][D[1]-1] = 1; 	B[D[1]-1][D[1]-1] = 1;
	B[D[2]-1][D[3]-1] = 1; 	B[D[3]-1][D[2]-1] = 1;
      }
      
      printAdjMatrixN(B, sys.nb);
      
    }    
    iter = iter + 1;    
  }while(iter < Miter); 
  printf("\nReached Final\n");
  // --------------------------- Memory Clearing Starts -----------------------
  
  free(bus);
  free(fd);
  free(td);
  free(pv);
  free(shunt);
  FreeMemory(A, sys.nb);
  FreeMemory(B, sys.nb);
  free(x);
  free(y);
  free(prod);
  free(madd);
  free(LC);
  free(Is1_Node);
  free(Is2_Node);
  free(generation);
  free(weightage);
  free(sirf);
  free(ini_vel);
  FreeMemoryMN(C, np, (2*nt));
  free(D);
  
  
  return 0;
}