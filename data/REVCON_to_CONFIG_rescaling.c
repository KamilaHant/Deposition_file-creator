
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAXATOMS            10000
#define T_expected          $T  /*kazdeho atomu zvlast*/
#define                  k  1.38e-23
#define                 au  1.67e-27

#define mass1               65.38  /*Zn*/
#define mass2               16.0   /*O*/

FILE *fconfig,*fconfig_new;
int i, poradi,atoms, atoms_nonfrozen = 0, cislo, citac =0,start_of_last_snapshot;
double faktor, M;
double bxl,bxh,byl,byh,bzl,bzh;
char radka[201];
short druh[MAXATOMS+1];
double q[MAXATOMS+1],x[MAXATOMS+1],y[MAXATOMS+1],z[MAXATOMS+1],vx[MAXATOMS+1],vy[MAXATOMS+1],vz[MAXATOMS+1],fx[MAXATOMS+1],fy[MAXATOMS+1],fz[MAXATOMS+1];
double v_aktualni, E_total = 0, T_aktualni; 

main() {

fconfig = fopen("REVCON","r");

while (fgets(radka,200,fconfig) != NULL) {
  citac++;
  if (strstr(radka,"ITEM: NUMBER OF ATOMS") != NULL) start_of_last_snapshot = citac;
}

fclose(fconfig);
fconfig = fopen("REVCON","r");
for(citac = 1;citac<=start_of_last_snapshot;citac++) fgets(radka,200,fconfig);
/*ted jsem na zacatku posledniho snapshotu*/

fgets(radka,200,fconfig);sscanf(radka,"%d",&atoms);

do {
  fgets(radka,200,fconfig);
} while (strstr(radka,"ITEM: BOX BOUNDS") == NULL);
fgets(radka,200,fconfig);sscanf(radka,"%lf %lf",&bxl,&bxh);
fgets(radka,200,fconfig);sscanf(radka,"%lf %lf",&byl,&byh);
fgets(radka,200,fconfig);sscanf(radka,"%lf %lf",&bzl,&bzh);

do {
  fgets(radka,200,fconfig);
} while (strstr(radka,"ITEM: ATOMS") == NULL);

for (i=1;i<=atoms;i++) {
  fgets(radka,200,fconfig);
  sscanf(radka,"%d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&cislo,&druh[i],&q[i],&x[i],&y[i],&z[i],&vx[i],&vy[i],&vz[i],&fx[i],&fy[i],&fz[i]);
  if (druh[i] == 1) M = mass1; if (druh[i] == 2) M = mass2;
  v_aktualni = sqrt(vx[i]*vx[i]+vy[i]*vy[i]+vz[i]*vz[i]);
  if (v_aktualni > 0) {
    E_total += 0.5*M*au*v_aktualni*v_aktualni*10000*1000000; /*J*/
    atoms_nonfrozen++;
  }  
}

T_aktualni = 2*E_total/(3*k*atoms_nonfrozen);
faktor = sqrt(T_expected/T_aktualni);

for (i=1;i<=atoms;i++) {
  vx[i] *= faktor; vy[i] *= faktor; vz[i] *= faktor;
}

fclose(fconfig);

fconfig_new = fopen("CONFIG","w");

fprintf(fconfig_new,"from REVCON\n\n%d atoms\n2 atom types\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nMasses\n\n1 %f\n2 %f\n\nAtoms\n\n",atoms,bxl,bxh,byl,byh,bzl,bzh,mass1,mass2);

for (i=1;i<=atoms;i++) {
  fprintf(fconfig_new,"%d %d %.4f %.6f %.6f %.6f\n",i,druh[i],q[i],x[i],y[i],z[i]);
}

fprintf(fconfig_new,"\nVelocities\n\n");

for (i=1;i<=atoms;i++) {
  fprintf(fconfig_new,"%d %.10f %.10f %.10f\n",i,vx[i],vy[i],vz[i]);
}

fclose(fconfig_new);
}
