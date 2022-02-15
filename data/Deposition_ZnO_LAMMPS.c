#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#define  frozen_substrate_1             $frozen1
#define  frozen_substrate_2             $frozen2
#define  unfrozen_substrate_1           $unfrozen1
#define  unfrozen_substrate_2           $unfrozen2

#define  frozen_substrate               (frozen_substrate_1+frozen_substrate_2)
#define  unfrozen_substrate             (unfrozen_substrate_1+unfrozen_substrate_2)

#define atoms                      $defineAtoms   /*novych Zn+ a O-*/

#define ions                       0    
#define Ehigh                      $Ehigh
#define Elow                       $Elow    /* eV */
#define R                          $R  /*procento high, pokud neni pocet ions urcen explicitne*/
#define x                          $x

#define mass1                      $mass1 /*Zn*/
#define mass2                      $mass2  /*O*/

#define charge1                     $charge1
#define charge2                     $charge2


#define druh1                      "$name1" /*zde fakticky nanic*/
#define druh2                      "$name2" /*zde fakticky nanic*/

#define bond12                     $bond12  /*pokud se nedeponuji molekuly tak nanic*/

#define MAXATOMS                   20000

#define cutoff_11                  $cutoff11 /* */
#define cutoff_12                  $cutoff12  /* */
#define cutoff_22                  $cutoff22 /* */

/*ZrN ma vzdalenos MN o 0.168 Angstromu delsi nez TiN*/
/*Foster et al: vazby ZrO v ZrO2 maji 2.07-2.28 A*/

/*cutoff11 a 22 nesmi byt nula i kdyz mne vazby nezajimaji; kvuli vzdalenosti nove priletajicich*/

/*empirical covalent radius 136 Ti, 73 O*/
/*bond length (rutile, anatase) 1.92-1.98*/

#define vyska_priletani            $depositionHigh /*nad maxz*/
#define tloustka                   $thickness /*angstroms nad maxz*/

#define max_kamaradu 10

double bxl,bxh,byl,byh,bzl,bzh;
FILE   *fconfig_old, *fconfig, *fsummary;
double souradnice[4*atoms+1][4], v[4*atoms+1][4];
double souradnice_old[MAXATOMS+1][4];
double velocities_old[MAXATOMS+1][4];
double q_old[MAXATOMS+1],q[atoms+1];
int    atoms_old = 0;
short  druh[4*atoms+1], druh_old[MAXATOMS+1], ion[4*atoms+1];
short  vyradit[MAXATOMS+1], vazan_k_substratu_old[MAXATOMS+1];
int    vyradit1=0, vyradit2=0;
int    new1=0, new2=0;
int    old1=0, old2=0;
int    deposited1, deposited2;
double nahodne;
int    i,j,k,l,m,molecule,nejblizsi,poradi = 1, citac = 0, citac_bez_vyhozenych = 0;
double maxz = -10000.0;
char   radka[201], temp[5];
double cutoff[4][4], maxcutoff;
short  koordinace[MAXATOMS+1],koordinace_noF[MAXATOMS+1];
int    kamaradi[MAXATOMS+1][max_kamaradu+1];

/****************************************************************************/

double vzdalenost_2stare(int stary, int novy) {
double pom1, pom2, pom3, vysledek=0;

pom1 = fabs(souradnice_old[stary][1] - souradnice_old[novy][1]);
while (pom1 > (bxh-bxl)/2) pom1 = pom1-(bxh-bxl);
pom2 = fabs(souradnice_old[stary][2] - souradnice_old[novy][2]);
while (pom2 > (byh-byl)/2) pom2 = pom2-(byh-byl);
pom3 = fabs(souradnice_old[stary][3] - souradnice_old[novy][3]);

vysledek = sqrt(pom1*pom1 + pom2*pom2 + pom3*pom3);

return vysledek;
}

/****************************************************************************/

double horizontalni_vzdalenost_2nove(int nejnovejsi, int novy) {
double pom1, pom2, pom3, vysledek;

pom1 = fabs(souradnice[nejnovejsi][1] - souradnice[novy][1]);
while (pom1 > (bxh-bxl)/2) pom1 = pom1-(bxh-bxl);
pom2 = fabs(souradnice[nejnovejsi][2] - souradnice[novy][2]);
while (pom2 > (byh-byl)/2) pom2 = pom2-(byh-byl);

vysledek = sqrt(pom1*pom1 + pom2*pom2);

return vysledek;
}

/****************************************************************************/

zjisti_stare_informace() {
int cislo,citac=0,start_of_last_snapshot;

fconfig_old = fopen("REVCON","r");
while (fgets(radka,200,fconfig_old) != NULL) {
  citac++;
  if (strstr(radka,"ITEM: NUMBER OF ATOMS") != NULL) start_of_last_snapshot = citac;
}
fclose(fconfig_old);

fconfig_old = fopen("REVCON","r");
do {
  fgets(radka,200,fconfig_old);
} while (strstr(radka,"ITEM: BOX BOUNDS") == NULL);
fgets(radka,200,fconfig_old);sscanf(radka,"%lf %lf",&bxl,&bxh);
fgets(radka,200,fconfig_old);sscanf(radka,"%lf %lf",&byl,&byh);
fgets(radka,200,fconfig_old);sscanf(radka,"%lf %lf",&bzl,&bzh); /*bzl se nacita celkme zbytecne: vim ze to je nula*/

fclose(fconfig_old);

fconfig_old = fopen("REVCON","r");

for(citac = 1;citac<=start_of_last_snapshot;citac++) fgets(radka,200,fconfig_old);
/*ted jsem na zacatku posledniho snapshotu*/

fgets(radka,200,fconfig_old);sscanf(radka,"%d",&atoms_old);

do {
  fgets(radka,200,fconfig_old);
} while (strstr(radka,"ITEM: ATOMS") == NULL);

for (i=1;i<=atoms_old;i++) {
  fgets(radka,200,fconfig_old);
  sscanf(radka,"%d %d %lf %lf %lf %lf %lf %lf %lf",&cislo,&druh_old[i],&q_old[i],&souradnice_old[i][1],&souradnice_old[i][2],&souradnice_old[i][3],&velocities_old[i][1],&velocities_old[i][2],&velocities_old[i][3]);
  if (druh_old[i] == 1) old1++;
  if (druh_old[i] == 2) old2++;
  while (souradnice_old[i][1] < -0.5*(bxh-bxl)) souradnice_old[i][1]+=(bxh-bxl);
  while (souradnice_old[i][1] >  0.5*(bxh-bxl)) souradnice_old[i][1]-=(bxh-bxl); 
  while (souradnice_old[i][2] < -0.5*(byh-byl)) souradnice_old[i][2]+=(byh-byl);
  while (souradnice_old[i][2] >  0.5*(byh-byl)) souradnice_old[i][2]-=(byh-byl);  
  while (souradnice_old[i][3] < 0)              souradnice_old[i][3]+=(bzh-bzl);
  while (souradnice_old[i][3] > (bzh-bzl))      souradnice_old[i][3]-=(bzh-bzl);  
}

fclose(fconfig_old);
}

/****************************************************************************/

vypocti_koordinace() {

  for (i=1;i<=atoms_old;i++) koordinace[i] = 0;
	
  for (i=1;i<=atoms_old;i++) {
    for (j=i+1;j<=atoms_old;j++) {
      if (vzdalenost_2stare(i,j) < cutoff[druh_old[i]][druh_old[j]]) {
        koordinace[i]++;kamaradi[i][koordinace[i]] = j;
        koordinace[j]++;kamaradi[j][koordinace[j]] = i;
      }
    }
  }
}

/****************************************************************************/
eliminuj_odprasene() { /*i ty prilepene zespoda na substrat (dan za periodicitu)*/
short pokrok;

for (i=1;i<=frozen_substrate;i++)    vazan_k_substratu_old[i] = 1;
for (i=frozen_substrate+1;i<=atoms_old;i++) vazan_k_substratu_old[i] = 0;

do {
  pokrok = 0;
  for (i=1;i<=atoms_old;i++) {
    for (j=i+1;j<=atoms_old;j++) {
    /*if (vazan_k_substratu_old[i] == 1 && vazan_k_substratu_old[j] == 0 && vzdalenost_2stare(i,j) < cutoff[druh_old[i]][druh_old[j]] && !(i<=frozen_substrate && souradnice_old[j][3] > (bzh-bzl)/2) ) {*/
      if (vazan_k_substratu_old[i] == 1 && vazan_k_substratu_old[j] == 0 && vzdalenost_2stare(i,j) < cutoff[druh_old[i]][druh_old[j]] && !(i<=0.333*(frozen_substrate+unfrozen_substrate) && j>(frozen_substrate+unfrozen_substrate) ) ) {
        vazan_k_substratu_old[j] = 1;
        pokrok = 1;
      }
    /*if (vazan_k_substratu_old[i] == 0 && vazan_k_substratu_old[j] == 1 && vzdalenost_2stare(i,j) < cutoff[druh_old[i]][druh_old[j]] && !(j<=frozen_substrate && souradnice_old[i][3] > (bzh-bzl)/2)) {*/
      if (vazan_k_substratu_old[i] == 0 && vazan_k_substratu_old[j] == 1 && vzdalenost_2stare(i,j) < cutoff[druh_old[i]][druh_old[j]] && !(j<=0.333*(frozen_substrate+unfrozen_substrate) && i>(frozen_substrate+unfrozen_substrate) ) ) {      	
        vazan_k_substratu_old[i] = 1;
        pokrok = 1;
      }
    }
  }
} while (pokrok == 1);

for (i=1;i<=atoms_old;i++) {vyradit[i] = 0; if (vazan_k_substratu_old[i] == 0) vyradit[i] = 1;}

for (i=1;i<=atoms_old;i++) if (souradnice_old[i][3] > maxz && vyradit[i] == 0) maxz = souradnice_old[i][3];

}
/****************************************************************************/

double blizko(int ktery) {
double vysledek = 0;
for (l=1;l<ktery;l++) {if (horizontalni_vzdalenost_2nove(ktery,l) <= maxcutoff) vysledek++;}
return vysledek;
}

/****************************************************************************/
vygeneruj_nove_atomy() {

/*pripravne prace*/
for (i=1;i<=atoms_old;i++) if (vyradit[i] == 1) {
  if (druh_old[i] == 1) vyradit1++;
  if (druh_old[i] == 2) vyradit2++;
}
deposited1 = old1  - vyradit1 - unfrozen_substrate_1 - frozen_substrate_1;
deposited2 = old2  - vyradit2 - unfrozen_substrate_2 - frozen_substrate_2;
/*konec pripravnych praci*/

for (i=1;i<=atoms;i++) {

  do {
    nahodne = rand();
    souradnice[i][1] = (bxh-bxl) * (nahodne + 0.0000001) / (RAND_MAX + 0.0000002) - (bxh-bxl)/2.0;
    nahodne = rand();
    souradnice[i][2] = (byh-byl) * (nahodne + 0.0000001) / (RAND_MAX + 0.0000002) - (byh-byl)/2.0;
    souradnice[i][3] = maxz + vyska_priletani;
  } while (blizko(i) > 0);
  /*  dopadajici atomy navzajem interagovat nebudou*/

  nahodne = rand();
  ion[i] = 0;
  if ( (nahodne + 0.0000001)/RAND_MAX < R ) ion[i] = 1;
  if ( i<=ions )                            ion[i] = 1;

  v[i][1] = 0;
  v[i][2] = 0;
  
  if (x*(deposited1+new1) < (deposited2+new2)) { 
    druh[i] = 1; 
    new1++;
    v[i][3] = -0.001 * ( sqrt( (Elow + Ehigh*ion[i] - Elow*ion[i]) * 19294.5 / mass1 ) );
    q[i] = charge1;
  }
  
  else  { 
    druh[i] = 2; 
    new2++;
    v[i][3] = -0.001 * ( sqrt( (Elow + Ehigh*ion[i] - Elow*ion[i]) * 19294.5 / mass2 ) );
    q[i] = charge2;
  }
  
}

}

/****************************************************************************/
tisk_Config() {

fconfig = fopen("CONFIG","w");

fprintf(fconfig,"newly generated\n\n%d atoms\n2 atom types\n%f %f xlo xhi\n%f %f ylo yhi\n%f %f zlo zhi\n\nMasses\n\n1 %f\n2 %f\n\nAtoms\n\n",atoms_old+new1+new2-vyradit1-vyradit2,bxl,bxh,byl,byh,0.0,maxz+tloustka,mass1,mass2);

for (i=1;i<=atoms_old;i++) {
  if (vyradit[i] == 0) {
    citac_bez_vyhozenych++;
    fprintf(fconfig,"%d %d %.4f %.6f %.6f %.6f\n",citac_bez_vyhozenych,druh_old[i],q_old[i],souradnice_old[i][1],souradnice_old[i][2],souradnice_old[i][3]);
  }
}

for (i=1;i<=new1+new2;i++) {
  citac_bez_vyhozenych++;
  fprintf(fconfig,"%d %d %.4f %.6f %.6f %.6f\n",citac_bez_vyhozenych,druh[i],q[i],souradnice[i][1],souradnice[i][2],souradnice[i][3]);
}

fprintf(fconfig,"\nVelocities\n\n");
citac_bez_vyhozenych = 0;

for (i=1;i<=atoms_old;i++) {
  if (vyradit[i] == 0) {
    citac_bez_vyhozenych++;
    fprintf(fconfig,"%d %.10f %.10f %.10f\n",citac_bez_vyhozenych,velocities_old[i][1],velocities_old[i][2],velocities_old[i][3]);
  }
}

for (i=1;i<=new1+new2;i++) {
  citac_bez_vyhozenych++;
  fprintf(fconfig,"%d %.10f %.10f %.10f\n",citac_bez_vyhozenych,v[i][1],v[i][2],v[i][3]);
}

}

/****************************************************************************/

tisk_Summary() {
fsummary = fopen("Summary","a");
fprintf(fsummary,"vyrazuji %2d %2d  pridavam %2d %2d  slozeni  %2d %2d\n",vyradit1, vyradit2, new1, new2, deposited1+new1, deposited2+new2);
fclose(fsummary);
}

/****************************************************************************/

cutoffs() {
cutoff[1][1] = cutoff_11; cutoff[1][2] = cutoff_12;
cutoff[2][1] = cutoff_12; cutoff[2][2] = cutoff_22;

maxcutoff = cutoff_11;
}

/****************************************************************************/

main() {
srand ( 10000*time(NULL) );

cutoffs();
zjisti_stare_informace();
vypocti_koordinace();

eliminuj_odprasene();

vygeneruj_nove_atomy();

tisk_Summary();
tisk_Config();

/*system("PAUSE");*/

}
