#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define input   "CONFIG"
#define output  "SNAPSHOT"

#define druh1 "Zn"
#define druh2 "O"

FILE *fconfig, *fconfig_new;
char radka[201];
char druh[5];
double x,y,z,q;
int i,atoms,druh_ciselne,cislo;

main() {

fconfig = fopen("CONFIG","r");
fconfig_new = fopen("SNAPSHOT","w");

fgets(radka,200,fconfig);
fgets(radka,200,fconfig);
fgets(radka,200,fconfig);

sscanf(radka,"%d",&atoms);

do {
  fgets(radka,200,fconfig);
} while (strstr(radka,"Atoms") == NULL);
fgets(radka,200,fconfig);

for (i=1;i<=atoms;i++) {
  fgets(radka,200,fconfig);
  sscanf(radka,"%d %d %lf %lf %lf %lf",&cislo,&druh_ciselne,&q,&x,&y,&z);
  if (druh_ciselne == 1) strcpy(druh,druh1);
  if (druh_ciselne == 2) strcpy(druh,druh2);
  fprintf(fconfig_new,"%s %f %f %f\n",druh,x,y,z);
}

fclose(fconfig);
fclose(fconfig_new);
}

