
#include <stdlib.h> 
#include <stdio.h>

#include <libecp.h>
/* dimensions.h is part of the non-public libecp interface */
#include <dimensions.h>

#include <getIntegrals.h>

int
loadBS(const char *filename, const int nrAtoms, 
       int **BSShells, int **BSl, int **BSK, 
       double **BSa, double **BSd) {
  FILE *fp = fopen(filename, "r");
  int i, j, k, l, idx; 
  int nsh = 0;
  int nprim = 0;
  int *NShells = calloc(nrAtoms, sizeof(int));
  int *am, *contraction;
  double *d, *a;
  char buf[512];
  for(i=0; i<nrAtoms; i++) {
    fscanf(fp, "%d %d\n", &idx, &(NShells[i]));
    for(j=0; j<NShells[i]; j++) {
      fscanf(fp, "%d %d\n", &l, &idx);
      for(k=0; k<idx; k++) {
	fgets(buf, sizeof(buf), fp);
      }
      nprim += idx;
    }
    nsh += NShells[i];
  }
  rewind(fp);
  *BSShells = NShells;
  *BSl = am = calloc(nsh, sizeof(int));
  *BSK = contraction = calloc(nsh, sizeof(int));
  *BSa = a = calloc(nprim, sizeof(double));
  *BSd = d = calloc(nprim, sizeof(double));  
  nsh = nprim = 0;
  for(i=0; i<nrAtoms; i++) {
    fscanf(fp, "%d %d\n", &idx, &(NShells[i]));
    printf("BS %d:  %d shells\n", idx, NShells[i]);
    for(j=0; j<NShells[i]; j++) {
      fscanf(fp, "%d %d\n", &(am[nsh]), &(contraction[nsh]));
      printf("  l=%d, K=%d\n", am[nsh], contraction[nsh]);
      for(k=0; k<contraction[nsh]; k++) {
	fscanf(fp, "%d %lf %lf", &idx, &(a[nprim]), &(d[nprim]));
	printf("    %12.6f  % 12.6f\n", a[nprim], d[nprim]);
	nprim++;
      }
      nsh++;
    }
  }
  fclose(fp);
  return nsh;
}

void
loadGeometry(const char *xyzFile, int *nrAtoms, double **geometry) {
  /* load from xyz file */
  FILE *fp = fopen(xyzFile, "r"); 
  int i;
  char buf[512];
  double *geom;
  fscanf(fp, "%d\n", nrAtoms); 
  fscanf(fp, "%[^\n]\n", buf); /* potential buffer overrun! */
  *geometry = geom = calloc((*nrAtoms)*3, sizeof(double));
  for(i=0; i<(*nrAtoms); i++) {
    fscanf(fp, "%s %lf %lf %lf\n", buf, &(geom[i*3]), &(geom[i*3+1]), &(geom[i*3+2]));
  }
  fclose(fp);
}

void
loadECP(const char *filename, const int nrAtoms, 
	int **ECPShells, int **ECPl, int **ECPK, 
	double **ECPa, double **ECPd, double **ECPn) {
  FILE *fp = fopen(filename, "r");
  int i, j, k, l, idx; 
  int nsh = 0;
  int nprim = 0;
  int *NShells = calloc(nrAtoms, sizeof(int));
  int *nCore = calloc(nrAtoms, sizeof(int));
  int *am, *contraction;
  double *d, *a, *n;
  char buf[512];
  for(i=0; i<nrAtoms; i++) {
    fscanf(fp, "%d %d %d\n", &idx, &(nCore[i]), &(NShells[i]));
    for(j=0; j<NShells[i]; j++) {
      fscanf(fp, "%d %d\n", &l, &idx);
      for(k=0; k<idx; k++) {
	fgets(buf, sizeof(buf), fp);
      }
      nprim += idx;
    }
    nsh += NShells[i];
  }
  rewind(fp);
  *ECPShells = NShells;
  *ECPl = am = calloc(nsh, sizeof(int));
  *ECPK = contraction = calloc(nsh, sizeof(int));
  *ECPa = a = calloc(nprim, sizeof(double));
  *ECPd = d = calloc(nprim, sizeof(double));  
  *ECPn = n = calloc(nprim, sizeof(double));
  nsh = nprim = 0;
  for(i=0; i<nrAtoms; i++) {
    fscanf(fp, "%d %d %d\n", &idx, &(nCore[i]), &(NShells[i]));
    printf("ECP %d: nCore=%d  %d shells\n", idx, nCore[i], NShells[i]);
    for(j=0; j<NShells[i]; j++) {
      fscanf(fp, "%d %d\n", &(am[nsh]), &(contraction[nsh]));
      printf("  l=%d, K=%d\n", am[nsh], contraction[nsh]);
      for(k=0; k<contraction[nsh]; k++) {
	fscanf(fp, "%d %lf %lf %lf", &idx, &(a[nprim]), &(d[nprim]), &(n[nprim]));
	printf("    %12.6f  % 12.6f  % 1.6f\n", a[nprim], d[nprim], n[nprim]);
	nprim++;
      }
      nsh++;
    }
  }
  fclose(fp);
  free(nCore);
}

int 
main(int argc, char **argv) {
  double *geometry;
  int nrAtoms;
  /* effective core potentials */
  int *ECPShells, *ECPl, *ECPK;
  double *ECPa, *ECPd, *ECPn;
  /* basis set */
  int *BSShells, *BSl, *BSK;
  double *BSa, *BSd;
  /* total number of BS shells := sum(shells) */
  int nrShells, dim;
  int i, j, idx;
  int s1, s2, shIdx1, shIdx2;
  int a, b;
  /* integral matrix */
  double *I, *tmp;
  /* offsets for ao matrix */
  int *ao, aoDim, lstart, maxShells;

  if (4 != argc) {
    printf("Usage: %d %s structure.xyz ECP BS\n", argc,argv[0]);
    return -1;
  }
  /* initialize */
  loadGeometry(argv[1], &nrAtoms, &geometry);  
  loadECP(argv[2], nrAtoms, &ECPShells, &ECPl, &ECPK, &ECPa, &ECPd, &ECPn);
  nrShells = loadBS(argv[3], nrAtoms, &BSShells, &BSl, &BSK, &BSa, &BSd);
  dim = idx = 0;
  for(i=0; i<nrAtoms; i++) {
    for(j=0; j<BSShells[i]; j++) {
      dim += IJK_DIM(BSl[idx]);
      idx++;
    }
  }
  I = calloc(dim*dim, sizeof(double));
  
  getIntegrals (nrAtoms, geometry, 
		ECPShells, ECPK, ECPl, 
		ECPn, ECPd, ECPa, 
		BSShells, BSl, BSK, BSd, BSa, 
		1024, 1.0E-12, 1.0E-14, 
		dim, I);

  /* FIXME - print integral matrix */
  printf("\n\nECP matrix elements\n\n");
  /* initialize aoDim */
  idx = lstart = maxShells = 0;
  for(i=0; i<nrAtoms; i++) {
    if (BSShells[i] > maxShells)
      maxShells = BSShells[i];
  }
  aoDim = maxShells;
  ao = calloc(nrAtoms*dim, sizeof(int));
  for (i=0; i<nrAtoms; i++) {
    for (j=0; j<BSShells[i]; j++) {
      ao[i*aoDim+j] = lstart;
      lstart += IJK_DIM(BSl[idx]);
      idx++;
    }
  }
  idx = 0;
  shIdx1 = shIdx2 = 0;
  for(i=0; i<nrAtoms; i++) {
    for(s1=0; s1<BSShells[i]; s1++) {
      for(j=i; j<nrAtoms; j++) {     
	shIdx1 = idx;
	shIdx2 = shIdx1;
	for(s2=s1; s2<BSShells[j]; s2++) {
	  printf("shell block (%d,%d: l=%d) x (%d,%d: l=%d)\n", i, s1, BSl[shIdx1], j, s2, BSl[shIdx2]);
	  /* print shell block */
	  for(a=0; a<IJK_DIM(BSl[shIdx1]); a++) {
	    for(b=0; b<IJK_DIM(BSl[shIdx2]); b++) {
	      printf(" % 10.6f  ", I[(a+ao[i*aoDim+s1])*dim+(b+ao[j*aoDim+s2])]);
	    }
	    printf("\n");
	  }
	  printf("\n");
	  shIdx2++;
	}
	shIdx1++;
      }
      idx = shIdx1;
    }
  }

  /* step through integral matrix element wise */
  printf("\n\nstepping through integral matrix element wise\n\n");
  tmp = I;
  idx = shIdx1 = shIdx2 = 0;
  for(i=0; i<nrAtoms; i++) {
    for(s1=0; s1<BSShells[i]; s1++) {
      for(a=0; a<IJK_DIM(BSl[shIdx1]); a++) {

	shIdx2 = 0;
	for(j=0; j<nrAtoms; j++) {    
	  for(s2=0; s2<BSShells[j]; s2++) {
	    printf("shell block (%d,%d: l=%d) x (%d,%d: l=%d)\n", i, s1, BSl[shIdx1], j, s2, BSl[shIdx2]);
	    for(b=0; b<IJK_DIM(BSl[shIdx2]); b++) {
	      printf(" % 10.6f  ", *tmp);
	      tmp++;
	    }
	    printf("\n");
	    shIdx2++;
	  }
	  printf("\n");
	}

      }
      shIdx1++;
    }
  }


  /* cleanup */
  free(ao);
  free(geometry);
  free(ECPShells);
  free(ECPK);
  free(ECPl);
  free(ECPn);
  free(ECPd);
  free(ECPa);
  free(BSShells);
  free(BSl);
  free(BSK);
  free(BSd);
  free(BSa);
  free(I);

  return 0;
}
