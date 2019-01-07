/****************************************************************************************
 *                         Graph Reduction using Graph Theory                           *
 * Mailid: vasubdevan@yahoo.com                                                         *
 * **************************************************************************************/
#include <stdio.h>
#include <stdlib.h>

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

// Node details
typedef struct {
  int NoEnodes;
  int *path;
  int **PathEle;
  int *Enodes;
}NodeDetails;

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
  printf("\nNo of feeders = %d", sys.nfd);
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
  printf("\n%d %d %lf %lf %lf %lf %lf", slack.sno, slack.busno, slack.MinPG, slack.MaxPG, slack.MinQG, slack.MaxQG, slack.Vpu);
  
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
  }
}

// --------------------  Generalized matrix Fucntions -------------------------

int **InitializeIntMatrixmxm(int **A, int n) {
  // Initializing the matrix pointer 
  int i, j; 
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A;
}

int ** CreateIntMatrixmxm(int ** A, int n) {
  // Creating Integer matrix of size n*n
  int i;
  
  A  = (int **)malloc(n * sizeof(int *));
  
  if(A == NULL) {
    printf("\nError in memory allocation\n\n");
    exit(1);
  }  else {
    for(i = 0; i < n; i++) {
      A[i] = (int *)malloc(n * sizeof(int));
      if(A[i] == NULL) {
	printf("\nError in memory allocation\n\n");
	exit(1);
      }
    }
  }
  
  A = InitializeIntMatrixmxm(A, n);
      
  return A;
}

void PrintSquareMatrix(int **pLMat, int n) {
  int i, j;
  printf("\n\nMatrix is as follows\n\n");
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      printf("%d\t", pLMat[i][j]);
    }
    printf("\n");
  }
}

int **InitializeIntMatrixmxn(int **A,int m, int n) {
  // Initializing the matrix pointer 
  int i, j; 
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      A[i][j] = 0;
    }
  }
  return A;
}

int ** CreateIntMatrixmxn(int ** A, int m, int n) {
  // Creating Integer matrix of size n*n
  int i;
  
  A  = (int **)malloc(m * sizeof(int *));
  
  if(A == NULL) {
    printf("\nError in memory allocation\n\n");
    exit(1);
  }  else {
    for(i = 0; i < m; i++) {
      A[i] = (int *)malloc(n * sizeof(int));
      if(A[i] == NULL) {
	printf("\nError in memory allocation\n\n");
	exit(1);
      }
    }
  }
  
  A = InitializeIntMatrixmxn(A, m, n);
      
  return A;
}

void PrintMatrixmxn(int **A, int m, int n) {
  int i, j;
  printf("\n\nMatrix is as follows\n\n");
  for(i = 0; i < m; i++) {
    for(j = 0; j < n; j++) {
      printf("%d\t", A[i][j]);
    }
    printf("\n");
  }
}

int *initializeIntArray(int *A, int n) {
  // Initializing the array pointer
  int i;
  
  for(i = 0; i < n; i++) {
    A[i] = 0;
  }
  return A;
}

int *CreateIntArray( int *A, int n) {
  // Creating an integer array
  int i;
  
  A = (int *)malloc( n * sizeof(int) );
  if(A == NULL)
    printf("\nError in memory allocation\n\n");
 
  A = initializeIntArray(A, n);
  
  return A;
}

void printIntArray(int *A, int n) {
  int i;
  printf("\nInteger Array is as follows\n\n");
  for(i = 0; i < n; i++) {
    printf("%d\t", A[i]);
  }
}

// ---------------------  Determine Nodes, branches & loop---------------------

int **AdjMatFormation(int **adj, FeederData *fd, SysInfo sys){
  // Creating adjacency matrix for any given network
  int i, j, k;
  
  adj = CreateIntMatrixmxm(adj, sys.nb);
  
  for( i = 0; i < sys.nb; i++) {
    for( j = 0 ; j < sys.nb; j++) {
	if( i == j)
	  adj[i][j] = 0; // diagonal elemnets.
	if(j < i) { // Lower diagonal elements
         for( k = 0; k < sys.nfd; k++) {
	   if( (fd[k].fb == i+1) && ( fd[k].tb == j+1 ) && (fd[k].FStat == 1) ) {
	     adj[i][j] = 1;
	     adj[j][i] = 1;
	   }
	   if( (fd[k].tb == i+1) && ( fd[k].fb == j+1 ) && (fd[k].FStat == 1) ) {
	     adj[i][j] = 1;
	     adj[j][i] = 1;
	   }
	 }
	}
    }
  }
  return adj;
}

int **IncedenceMatFormation(int **InMat, FeederData *fd, SysInfo sys) {
  
  int i, j, t1, t2, t3;
  t1 = 0; t2 = 0; t3 = 0;
    
  
  InMat = CreateIntMatrixmxn(InMat, (sys.nb-1), sys.nb);
  // PrintMatrixmxn(InMat, (sys.nb-1), sys.nb);
  
  for(i = 0; i < (sys.nb-1); i++) {
    t1 = fd[i].fb;
    t2 = fd[i].tb;
    for(j = 0; j < sys.nb; j++) {
      t3 = j+1;
      if(t1 == t3) {
	InMat[i][j] = -1;
      }
      if(t2 == t3) {
	InMat[i][j] = 1; 
      }
      t3 = 0;
    }
  }
  
  // PrintMatrixmxn(InMat, (sys.nb-1), sys.nb);
  
  return InMat;
}

int **DeterminPathLowerTriMatrix(int **pLMat, FeederData *fd, SysInfo sys) {
  // Identifying the path lower triangular matrix
  int i, j, k;
  pLMat = CreateIntMatrixmxm(pLMat, sys.nb);
  
  for( i = 0; i < sys.nb; i++) {
    for( j = 0 ; j < sys.nb; j++) {
	if( i == j)
	  pLMat[i][j] = 1; // Making diagonal elemnet to be 1.
	if(j < i) { // Lower diagonal elements
         for( k = 0; k < sys.nfd; k++) {
	   if( (fd[k].fb == i+1) && ( fd[k].tb == j+1 ) && (fd[k].FStat == 1) )
	     pLMat[i][j] = 1;
	   if( (fd[k].tb == i+1) && ( fd[k].fb == j+1 ) && (fd[k].FStat == 1) )
	     pLMat[i][j] = 1;
	 }
	}
    }
  }
  return pLMat;
}

int *SumLowerTraMatrix(int *SumLM, int **pLMat, SysInfo sys) {
  // Determine the sum of lower traingular matrix
  int i, j, count;
  
  count =0;
  
  SumLM = CreateIntArray(SumLM, sys.nb);
  for(i = 0; i < sys.nb; i++) {
    for(j = 0; j < sys.nb; j++) {
      count = pLMat[j][i] + count;      
    }
    SumLM[i] = count;
    count = 0;
  }
  return SumLM;
}

int *CalculateDegree(int *degree, int **A, int n) {
  int i, j;
  int cnt;
    
  degree = CreateIntArray(degree, n);
    
  for(i = 0; i < n; i++) {
    cnt = 0;
    for(j = 0; j < n; j++) {
      cnt += A[j][i];
    }
    degree[i] = cnt;
  }
  return degree;
}

int DetermineSelfLoop(int NoSelfLoop, int **B, int n) {
  int i, j, nloop;
  nloop = 0;
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      if( (i == j) && ( B[i][j] == 1 ) )
	nloop++;
    }
  }
  NoSelfLoop = nloop;
  return NoSelfLoop;
}

int CheckSymmetryNature(int Sym, int n, int **B) {
  int i, j, k;
  Sym = 0;
  
  int **trans = CreateIntMatrixmxm(trans, n);
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      trans[j][i] = B[i][j];
    }
  }
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
      if(trans[i][j] != B[i][j]) {
	Sym++;
      }
    }
  }
  
  // Memory Clearing 
  for(i = 0; i < n; i++) {
    free(trans[i]);
  }
  free(trans);
  
 
  return Sym;
}

void ConditionCheck(int a) {
  if(a == 0)
    printf("\nThe given network is symmetric and undirected\n");
  else
    printf("\nThe given network is asymmetric and directed network\n");
}

int *DetermineSumAdjaMatrix(int *SumAdjMat, int **B, int n) {
  int i, j, k;
  
  k = 0;
  SumAdjMat = CreateIntArray(SumAdjMat, n);
  
  for(i = 0; i < n; i++) {
    for(j = 0; j < n; j++) {
       k += B[i][j]; 
    }
    SumAdjMat[i] = k;
    k = 0;
  }
  
  return SumAdjMat;
}

int DetermineNoofEdges(int NoofEdges, int NoSelfLoop, int Sym, int *SumAdjMat, int n ) {
  int i, j;
  
  NoofEdges = 0;
  
  if( (NoSelfLoop == 0) && (Sym == 0) ) { // Undirected Simple graph - Symmetric graph
    for(i = 0; i < n; i++) {
      NoofEdges += SumAdjMat[i];
    }
    NoofEdges = NoofEdges/2;
  }
  
  if( (NoSelfLoop != 0) && (Sym ==0) ) {  // Symmetric along with self loop
    for( i = 0; i < n; i++)
      NoofEdges += SumAdjMat[i];
    
    NoofEdges =  (NoofEdges - NoSelfLoop)/2 + NoSelfLoop;
  }

  return NoofEdges;
}

int DetermineNoofEnodes(int b, int *degree, int n) {
  int i;
  
  b = 0;
  
  for(i = 1; i < n; i++) {
    if(degree[i] == 1)
      b++;
  }
  
  return b;
}

int *DetermineEndNodes(NodeDetails nd,  int *degree, SysInfo sys) {
 // Determine the End node numbers
  int i, j;
  
  // Allocate memory for the pointer inside that structure
  nd.Enodes = CreateIntArray(nd.Enodes, nd.NoEnodes);
    
 for(i = 0,j = 0; i < sys.nb; i++) {
   if(i!=0) {
    if(degree[i] == 1) {
       nd.Enodes[j] = i+1;
      j++;
    }
  }
 }
 
  // printIntArray(nd.Enodes, nd.NoEnodes);
  
  return nd.Enodes;
}

int *DetermineNoofPath(NodeDetails nd) {
  // Which initialize the path based on the number of Endnodes
  int i;
  nd.path = CreateIntArray(nd.path, nd.NoEnodes);
  // printIntArray(nd.path, nd.NoEnodes);
  
  for(i = 0; i < nd.NoEnodes; i++) {
    nd.path[i] = i+1;
  }
  
  return nd.path;
}

int **DeterminePathElements(NodeDetails nd, int ncol, SysInfo sys, int **A) {
  
  int i, j, k, k1, l;
  int temp = 0;
  
  nd.PathEle = CreateIntMatrixmxn(nd.PathEle, nd.NoEnodes, ncol);
  
  //PrintMatrixmxn(A, sys.nb, sys.nb);
  
  
  for(i = nd.NoEnodes, j = 0; i >= 1, j < nd.NoEnodes; i--, j++) {
    nd.PathEle[j][0] = nd.Enodes[i-1];   
  }
  
  
  for(j=0, i = sys.nb; i >= 1 ; i--) {
    if( A[i-1][nd.Enodes[j]-1] == 1)  {
      temp = nd.Enodes[j]-1;
      l = 1; 
      for(k = i-1; k >= 1 ;k--) {
	// printf("\nk = %d",k);
	if(A[k][temp] == 1) {
	  nd.PathEle[j][l] = k;
	  l++;
	  
	}
      }
    }
  }
  
  
 
  PrintMatrixmxn(nd.PathEle, nd.NoEnodes, ncol);

  return nd.PathEle; 
}

NodeDetails ObtainEndNodeDetails(NodeDetails nd, SysInfo sys, int *degree, int ncol, int **A) {
  
  nd.NoEnodes = DetermineNoofEnodes(nd.NoEnodes, degree, sys.nb);
  // printf("\nNumber of end nodes is %d", nd.NoEnodes);
    
  nd.Enodes = DetermineEndNodes(nd, degree, sys);
  // printIntArray(nd.Enodes, nd.NoEnodes);
  
  nd.path = DetermineNoofPath(nd);
  // printIntArray(nd.path, nd.NoEnodes);
  
  nd.PathEle = DeterminePathElements(nd, ncol, sys, A);
      
   
  return nd;
}

void PrintNodeDetails(NodeDetails nd) {
  int i;
  printf("\nNumber of End nodes in the systems are %d\nEnd nodes are as follows:\n\n", nd.NoEnodes);
  for(i = 0; i < nd.NoEnodes; i++) {
    printf("%d\t",nd.Enodes[i]);
  }
}


// ------------------------------ Main Program Starts here --------------------

int main (int argc, char *argv[]) {
  
  int i;
  int ncol = 20;  // Randomly considering the twenty elements in a path.
  
  FileDetails fpd = ReadFileDetails(fpd);
  
  FILE *fp = FilePointerToRead(fp,fpd.input);
  
  SysInfo sys = ReadEntireSystemData(sys, fp);
  // PrintEntireSystemInfo(sys);
  
  BusData *bus = ReadBusData(bus, fp, sys.nb);
  // printBusData(bus, sys.nb);r
    
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
  
  // ------------------------------ Graph Theory - Radial check ---------------
  int **adj = AdjMatFormation(adj, fd,sys);
  // PrintSquareMatrix(adj, sys.nb);
  
  int **InMat = IncedenceMatFormation(InMat, fd, sys);
  
  
  int **pLMat =  DeterminPathLowerTriMatrix(pLMat, fd, sys);
  // PrintSquareMatrix(pLMat, sys.nb);
  
  int *SumLM = SumLowerTraMatrix(SumLM, pLMat, sys);
  // printIntArray(SumLM, sys.nb);
  
  int *degree = CalculateDegree(degree, adj, sys.nb);
  // printIntArray(degree, sys.nb);

  int NoSelfLoop = DetermineSelfLoop(NoSelfLoop, adj, sys.nb);
  // printf("\nSelf loop in the given network is %d", NoSelfLoop);
  
  int Sym = CheckSymmetryNature(Sym, sys.nb, adj);
  ConditionCheck(Sym);
  
  int *SumAdjMat = DetermineSumAdjaMatrix(SumAdjMat, adj, sys.nb);
  // printIntArray(SumAdjMat, sys.nb);
  
  int NoofEdges = DetermineNoofEdges(NoofEdges, NoSelfLoop, Sym, SumAdjMat, sys.nb);
  // printf("\nNumber of Edges in the network is %d", NoofEdges);

  NodeDetails nd = ObtainEndNodeDetails(nd, sys, degree, ncol, pLMat);
  // PrintNodeDetails(nd);
  
 // --------------------------- Memory Clearing Starts -----------------------
  
  free(bus);
  free(fd);
  free(td);
  free(pv);
  free(shunt);
  // free(PathMat);
  
  for(i = 0; i < sys.nb; i++) {
    free(adj[i]);
    free(pLMat[i]);
  }
  free(pLMat);
  free(adj);
  free(SumLM);
  free(degree);
  free(SumAdjMat);
  
  free(nd.Enodes);
  free(nd.path);
  for(i = 0; i < nd.NoEnodes; i++) {
    free(nd.PathEle[i]);
  }
  free(nd.PathEle); 
  for(i = 0; i < (sys.nb-1); i++) {
    free(InMat[i]);
  }
  free(InMat);
  
  printf("\nDone!!\n");
  return 0;
}