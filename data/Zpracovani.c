#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define printrings                       0

#define input                           "SNAPSHOT"

#define MAXATOMS                        10000

#define RDF_max                         10.0
#define RDF_krok                        0.05
#define RDF_prihradek                   (int) (RDF_max / RDF_krok)

#define uhel_max                        180.0
#define uhel_krok                       2.0
#define uhel_prihradek                  (int) (uhel_max / uhel_krok)

#define  frozen_substrate_1             $unfrozen1
#define  frozen_substrate_2             $unfrozen2
#define  unfrozen_substrate_1           $unfrozen1
#define  unfrozen_substrate_2           $unfrozen2

#define  frozen_substrate               (frozen_substrate_1+frozen_substrate_2)
#define  unfrozen_substrate             (unfrozen_substrate_1+unfrozen_substrate_2)

#define  substrate_1                    (frozen_substrate_1+unfrozen_substrate_1)
#define  substrate_2                    (frozen_substrate_2+unfrozen_substrate_2)

#define prvni                      "$name1"
#define druhy                      "$name2"

#define cutoff_11                  0.0
#define cutoff_12                  $cutoff12
#define cutoff_22                  0.0

#define max_kamaradu 10

#define cellx                      $cellx
#define celly                      $celly

#define max_ringlength             12
#define max_distance               (max_ringlength / 2)

short  **useky, useky_velikost, pocet_rings[max_ringlength + 1];
FILE   *fsnapshot, *fresults;
double souradnice[MAXATOMS+1][4];
int    atoms, atoms1, atoms2, deposited1, deposited2, bulk1, bulk2;
short  druh[MAXATOMS+1], vazan_k_substratu[MAXATOMS+1], povrchovy_atom[MAXATOMS+1], kamaradu_nahore[MAXATOMS+1], bulk[MAXATOMS+1], substrat[MAXATOMS+1];
double nahodne;
int    i,j,k,l,m,n,nejblizsi;
double maxz, minz;
char   radka[201], temp[5];
double cutoff[4][4], maxcutoff;
short  koordinace[MAXATOMS+1];
int    kamaradi[MAXATOMS+1][max_kamaradu+1];
int    bonds[4][4]; 
double bond_length[4][4],coord[4];
short  stack_cisla[(max_ringlength + 1)], stack_vzdalenosti[(max_ringlength + 1)], uvazovany_kamarad[(max_ringlength + 1)], stack_delka;

/****************************************************************************/

double vzdalenost(int stary, int novy) {
double pom1, pom2, pom3, vysledek;

pom1 = fabs(souradnice[stary][1] - souradnice[novy][1]);
while (pom1 > cellx/2) pom1 = pom1-cellx;
pom2 = fabs(souradnice[stary][2] - souradnice[novy][2]);
while (pom2 > celly/2) pom2 = pom2-celly;
pom3 = fabs(souradnice[stary][3] - souradnice[novy][3]);

vysledek = sqrt(pom1*pom1 + pom2*pom2 + pom3*pom3);

/*printf("[%d %.2f %.2f %.2f] - [%d %.2f %.2f %.2f]: %.2f\n",stary,souradnice[stary][1],souradnice[stary][2],souradnice[stary][3],novy,souradnice[novy][1],souradnice[novy][2],souradnice[novy][3],vysledek);*/

return vysledek;
}

/****************************************************************************/
double horizontalni_vzdalenost(int stary, int novy) {
double pom1, pom2, vysledek;

pom1 = fabs(souradnice[stary][1] - souradnice[novy][1]);
while (pom1 > cellx/2) pom1 = pom1-cellx;
pom2 = fabs(souradnice[stary][2] - souradnice[novy][2]);
while (pom2 > celly/2) pom2 = pom2-celly;

vysledek = sqrt(pom1*pom1 + pom2*pom2);

return vysledek;
}

/****************************************************************************/
zjisti_stare_informace() {
char druh_temp[5];

atoms = 0; atoms1 = 0; atoms2 = 0;
fsnapshot = fopen(input,"r"); 
while (fgets(radka,200,fsnapshot) != NULL) {
  atoms++;      
  sscanf(radka,"%s %lf %lf %lf",druh_temp,&souradnice[atoms][1],&souradnice[atoms][2],&souradnice[atoms][3]);
  if (strcmp(druh_temp,"Zn") == 0) {druh[atoms] = 1; atoms1++;if (atoms1 <= substrate_1) substrat[atoms] = 1;}
  if (strcmp(druh_temp,"O") == 0)  {druh[atoms] = 2; atoms2++;if (atoms2 <= substrate_2) substrat[atoms] = 1;}
  while (souradnice[atoms][1] < -0.5*cellx) souradnice[atoms][1]+=cellx;
  while (souradnice[atoms][1] >  0.5*cellx) souradnice[atoms][1]-=cellx; 
  while (souradnice[atoms][2] < -0.5*celly) souradnice[atoms][2]+=celly;
  while (souradnice[atoms][2] >  0.5*celly) souradnice[atoms][2]-=celly;    
}
fclose(fsnapshot);
}

/****************************************************************************/

vypocti_koordinace() {

useky_velikost = atoms + 1;
useky = (short **) malloc(useky_velikost * sizeof(short *));
for (i = 1; i<=atoms; i++) useky[i] = (short *) malloc(useky_velikost * sizeof(short));

for (i=1;i<=atoms;i++) for (j=1;j<=atoms;j++) useky[i][j] = 1000;
for (i=1;i<=atoms;i++) useky[i][i] = 0;
for (i=1;i<=atoms;i++) koordinace[i] = 0;
for (i=1;i<=atoms;i++) {
  for (j=i+1;j<=atoms;j++) {
    if (vzdalenost(i,j) < cutoff[druh[i]][druh[j]]) {
      koordinace[i]++;kamaradi[i][koordinace[i]] = j;
      koordinace[j]++;kamaradi[j][koordinace[j]] = i;
      useky[i][j] = 1;
      useky[j][i] = 1;
    }
  }
}

}

/****************************************************************************/
bonding_statistics() {
                     
for (i=1;i<=2;i++) for (j=1;j<=2;j++) {bonds[i][j] = 0; bond_length[i][j] = 0;}
for (i=1;i<=2;i++) coord[i] = 0;
  
for (i=1;i<=atoms;i++) if (bulk[i] == 1) coord[druh[i]] += koordinace[i]; 
coord[1] /= bulk1; coord[2] /= bulk2;

for (i=1;i<=atoms;i++) for (j=i+1;j<=atoms;j++) {
  if ((useky[i][j] == 1) && (bulk[i] == 1 || bulk[j] == 1)) {
    bond_length[druh[i]][druh[j]] += vzdalenost(i,j);
    bonds[druh[i]][druh[j]]++;
  }                               
}

bonds[1][2]       += bonds[2][1];       bonds[2][1] = bonds[1][2];
bond_length[1][2] += bond_length[2][1]; bond_length[2][1] = bond_length[1][2];    
for (i=1;i<=2;i++) for (j=1;j<=2;j++) if (bonds[i][j] > 0) bond_length[i][j] /= bonds[i][j];
}

/****************************************************************************/

eliminuj_odprasene() { /*i ty prilepene zespoda na substrat (dan za periodicitu)*/
short pokrok;

minz = 10000; maxz = -10000;
for (i=1;i<=frozen_substrate;i++)    {
  vazan_k_substratu[i] = 1;
  if (souradnice[i][3] < minz) minz = souradnice[i][3];
}        
for (i=1;i<=atoms;i++) souradnice[i][3] -= minz; /*nyni je spodek zamrzleho substratu na 0 - melo by to tak byt i predtim, ale pro jistotu*/
    
for (i=frozen_substrate+1;i<=atoms;i++) vazan_k_substratu[i] = 0;

do {
  pokrok = 0;
  for (i=1;i<=atoms;i++) {
    for (j=i+1;j<=atoms;j++) {
      if (vazan_k_substratu[i] == 1 && vazan_k_substratu[j] == 0 && vzdalenost(i,j) < cutoff[druh[i]][druh[j]] && !(i<=frozen_substrate && souradnice[j][3] < 0) ) {
        vazan_k_substratu[j] = 1;
        pokrok = 1;
      }
      if (vazan_k_substratu[i] == 0 && vazan_k_substratu[j] == 1 && vzdalenost(i,j) < cutoff[druh[i]][druh[j]] && !(j<=frozen_substrate && souradnice[i][3] < 0)) {
        vazan_k_substratu[i] = 1;
        pokrok = 1;
      }
    }
  }
} while (pokrok == 1);

for (i=1;i<=atoms;i++) if (souradnice[i][3] > maxz && vazan_k_substratu[i] == 1) maxz = souradnice[i][3];

}
/****************************************************************************/

eliminuj_povrchove() { /*povrchove budou atomy s zadnym atomem v kuzeli 45deg nad sebou; a pripadne vsechny k nim shora vazane*/               
FILE *fpovrch;
                     
for (i=1;i<=atoms;i++) povrchovy_atom[i] = 0;
for (i=1;i<=atoms;i++) kamaradu_nahore[i] = 0;

for (i=1;i<=atoms;i++) {
  for (j=i+1;j<=atoms;j++) {
    if (vazan_k_substratu[i] && vazan_k_substratu[j] && horizontalni_vzdalenost(i,j) < souradnice[i][3]-souradnice[j][3]) kamaradu_nahore[j]++;       
    if (vazan_k_substratu[i] && vazan_k_substratu[j] && horizontalni_vzdalenost(i,j) < souradnice[j][3]-souradnice[i][3]) kamaradu_nahore[i]++;     
  }
}    

for (i=1;i<=atoms;i++) if (kamaradu_nahore[i] == 0) {
  povrchovy_atom[i] = 1;
  for (j=1;j<=koordinace[i];j++) if (souradnice[i][3] < souradnice[kamaradi[i][j]][3]) povrchovy_atom[kamaradi[i][j]] = 1;
}  
 
for (i=1;i<=atoms;i++) bulk[i] = 0;
deposited1 = 0; deposited2 = 0; bulk1 = 0; bulk2 = 0;

for (i=1;i<=atoms;i++) {
  if (substrat[i] == 0 && vazan_k_substratu[i] == 1) {
    if (druh[i] == 1) deposited1++;                         
    if (druh[i] == 2) deposited2++;                        
    if (povrchovy_atom[i] == 0) {
      bulk[i] = 1; 
      if (druh[i] == 1) bulk1++;                        
      if (druh[i] == 2) bulk2++;       
    }
  }      
}

fpovrch = fopen("Surface.xyz","w");
for (i=1;i<=atoms;i++) {
  if (substrat[i])           fprintf(fpovrch,"substrate ");
  if (bulk[i])               fprintf(fpovrch,"bulk      ");    
  if (povrchovy_atom[i])     fprintf(fpovrch,"surface   ");
  fprintf(fpovrch,"%f %f %f\n",souradnice[i][1],souradnice[i][2],souradnice[i][3]);
}  
fclose(fpovrch); 
  
}
/****************************************************************************/

tisk_rdf() {
FILE *frdf;
double aktualni_vzdalenost;
short RDF[4][4][RDF_prihradek+1];
int umisteni;

for (i=1;i<=RDF_prihradek;i++) for (j=1;j<=2;j++) for (k=1;k<=2;k++) RDF[j][k][i] = 0;

for (i=1;i<=atoms;i++) {
  for (j=i+1;j<=atoms;j++) {
    if (bulk[i] && bulk[j]) {      
      aktualni_vzdalenost = vzdalenost(i,j);
/*      printf("%d-%d %d-%d : %f\n",i,druh[i],j,druh[j],aktualni_vzdalenost);*/
      umisteni = floor(aktualni_vzdalenost/RDF_krok);  
      if (umisteni <= RDF_prihradek) {RDF[druh[i]][druh[j]][umisteni]++;/*printf("rdf %d %d : %d (%f)\n",druh[i],druh[j],umisteni,aktualni_vzdalenost/RDF_krok))*/;}
    }
  }
}

frdf = fopen("RDF.txt","w");        
fprintf(frdf,"dist  Zn-Zn  Zn-O  O-O  (only_film_bulk)\n");       
for (i=1;i<=RDF_prihradek;i++) fprintf(frdf,"%.3f  %d  %d  %d\n",(i-0.5)*RDF_krok,RDF[1][1][i],RDF[1][2][i]+RDF[2][1][i],RDF[2][2][i]);  
fclose(frdf);

frdf = fopen("RDF_temp.txt","w");  
for (i=1;i<=RDF_prihradek;i++) fprintf(frdf,"%d\n",RDF[1][2][i]+RDF[2][1][i]);  
fclose(frdf);
}

/****************************************************************************/
double uhel(int stred, int levy, int pravy) {
double l,p,s,cos,vysledek;       

l = vzdalenost(stred, pravy);
p = vzdalenost(stred, levy);
s = vzdalenost(levy, pravy);

cos = (p*p+l*l-s*s)/(2*l*p);

vysledek = acos(cos);
vysledek *= (180.0/3.1415926);
if (vysledek < 0) vysledek += 360;
if (vysledek > 180) vysledek = 360 - vysledek;
return vysledek;
}

/****************************************************************************/

tisk_bonding_angles() {
FILE *fangles;
double aktualni_uhel;
int umisteni;
short uhel_212[RDF_prihradek+1], uhel_121[RDF_prihradek+1];

for (i=1;i<=uhel_prihradek;i++) {uhel_212[i] = 0; uhel_121[i] = 0;}
                    
for (i=1;i<=atoms;i++) {
    
    printf("*** %d ***\n",i);
    
  for (j=i+1;j<=atoms;j++) {
    for (k=j+1;k<=atoms;k++) {
      if (bulk[i] && bulk[j] && bulk[k]) {                                              
         if (useky[i][j] == 1 && useky[i][k] == 1 && druh[i] == 1 && druh[j] == 2 && druh[k] == 2) {
            aktualni_uhel = uhel(i,j,k);                                      
            umisteni = floor(aktualni_uhel/uhel_krok);
            if (umisteni <= uhel_prihradek) uhel_212[umisteni]++; /*pro maximum 180 je prvni podminka jen formalni*/
         }
         if (useky[i][j] == 1 && useky[i][k] == 1 && druh[i] == 2 && druh[j] == 1 && druh[k] == 1) {
            aktualni_uhel = uhel(i,j,k);                                      
            umisteni = floor(aktualni_uhel/uhel_krok);
            if (umisteni <= uhel_prihradek) uhel_121[umisteni]++; /*pro maximum 180 je prvni podminka jen formalni*/
         }         
         
         if (useky[j][i] == 1 && useky[j][k] == 1 && druh[j] == 1 && druh[i] == 2 && druh[k] == 2) {
            aktualni_uhel = uhel(j,i,k);                                      
            umisteni = floor(aktualni_uhel/uhel_krok);
            if (umisteni <= uhel_prihradek) uhel_212[umisteni]++; /*pro maximum 180 je prvni podminka jen formalni*/
         }
         if (useky[j][i] == 1 && useky[j][k] == 1 && druh[j] == 2 && druh[i] == 1 && druh[k] == 1) {
            aktualni_uhel = uhel(j,i,k);                                      
            umisteni = floor(aktualni_uhel/uhel_krok);
            if (umisteni <= uhel_prihradek) uhel_121[umisteni]++; /*pro maximum 180 je prvni podminka jen formalni*/
         }                                             

         if (useky[k][j] == 1 && useky[k][i] == 1 && druh[k] == 1 && druh[j] == 2 && druh[i] == 2) {
            aktualni_uhel = uhel(k,j,i);                                      
            umisteni = floor(aktualni_uhel/uhel_krok);
            if (umisteni <= uhel_prihradek) uhel_212[umisteni]++; /*pro maximum 180 je prvni podminka jen formalni*/
         }
         if (useky[k][j] == 1 && useky[k][i] == 1 && druh[k] == 2 && druh[j] == 1 && druh[i] == 1) {
            aktualni_uhel = uhel(k,j,i);                                      
            umisteni = floor(aktualni_uhel/uhel_krok);
            if (umisteni <= uhel_prihradek) uhel_121[umisteni]++; /*pro maximum 180 je prvni podminka jen formalni*/
         }   
      }                              
    }
  }
}
  
fangles = fopen("ANGLES.txt","w");        
fprintf(fangles,"angle  O-Zn-O  Zn-O-Zn (only_film_bulk)\n");       
for (i=1;i<=uhel_prihradek;i++) fprintf(fangles,"%.1f  %d  %d\n",(i-0.5)*uhel_krok,uhel_212[i],uhel_121[i]);  
fclose(fangles);  
                      
}

/****************************************************************************/

short otestuj_ring_SP() {
  short m, polomer, m_kamarad, vysledek;

  vysledek = 1; 
  polomer = stack_delka / 2;
  /*pro delku 7 je vysledek 3, coz je ok*/

  for (m = 1; m<=stack_delka; m++) {
    m_kamarad = m + polomer;
    if (m_kamarad > stack_delka) m_kamarad = m_kamarad - stack_delka;
    if (useky[stack_cisla[m]][stack_cisla[m_kamarad]] < polomer) vysledek = 0;
  }
  return vysledek;
}

/****************************************************************************/

short otestuj_ring_vbunce() {
  short m, vysledek;
  double deltax, deltay, deltaz, deltaxsum, deltaysum, deltazsum; /*souradnice z zde nehraje roli*/
  /*predpokladam ze bunka je vetsi nez dvojnasobek nejdelsi vazby*/

  deltaxsum = 0;
  deltaysum = 0;
  deltazsum = 0;

  stack_delka++;
  /*docasna zmena v ramci procedury*/
  stack_cisla[stack_delka] = stack_cisla[1];

  vysledek = 0;
  
  for (m =2; m<=stack_delka; m++) {
    deltax = souradnice[stack_cisla[m]][1] - souradnice[stack_cisla[(m-1)]][1];
    deltay = souradnice[stack_cisla[m]][2] - souradnice[stack_cisla[(m-1)]][2];
    deltaz = souradnice[stack_cisla[m]][3] - souradnice[stack_cisla[(m-1)]][3];
    if (deltax > (cellx/2) )    deltax = deltax - cellx;
    if (deltax < (-1*cellx/2) ) deltax = deltax + cellx;
    if (deltay > (celly/2) )    deltay = deltay - celly;
    if (deltay < (-1*celly/2) ) deltay = deltay + celly;
/*  if (deltaz > (cellz/2) )    deltaz = deltaz - cellz;
    if (deltaz < (-1*cellz/2) ) deltaz = deltaz + cellz;  */
    deltaxsum = deltaxsum + deltax;
    deltaysum = deltaysum + deltay;
    deltazsum = deltazsum + deltaz;
  }
 
  if (deltaxsum < 0) deltaxsum = deltaxsum * -1;
  if (deltaysum < 0) deltaysum = deltaysum * -1;
  if (deltazsum < 0) deltazsum = deltazsum * -1;

  if ( (deltaxsum < 0.1) && (deltaysum < 0.1) /*&& (deltazsum < 0.1)*/ ) vysledek = 1;

  stack_delka--;

  return vysledek;
}

/****************************************************************************/

short unimodal_labelling() {
  short vysledek, predtimvzestup, predtimsestup, predtimrovnost, tedvzestup, tedsestup, tedrovnost, m, bylarovnost;
  
  bylarovnost = 0;  
  predtimvzestup = 0;
  predtimsestup = 0;
  predtimrovnost = 0;
  tedvzestup = 0;
  tedsestup = 0;
  tedrovnost = 0;
  vysledek = 1;  

  for (m = 2; m<stack_delka; m++) {
    if (stack_vzdalenosti[m] == stack_vzdalenosti[(m-1)]) {
      bylarovnost = 1;
    } 
  }
  
  if (stack_delka > 2) {
    if (stack_vzdalenosti[(stack_delka - 1)] - stack_vzdalenosti[(stack_delka - 2)] == 1) predtimvzestup = 1;
    if (stack_vzdalenosti[(stack_delka - 1)] - stack_vzdalenosti[(stack_delka - 2)] == -1) predtimsestup = 1;
    if (stack_vzdalenosti[(stack_delka - 1)] - stack_vzdalenosti[(stack_delka - 2)] == 0) predtimrovnost = 1;

    if (stack_vzdalenosti[(stack_delka)] - stack_vzdalenosti[(stack_delka - 1)] == 1) tedvzestup = 1;
    if (stack_vzdalenosti[(stack_delka)] - stack_vzdalenosti[(stack_delka - 1)] == -1) tedsestup = 1;
    if (stack_vzdalenosti[(stack_delka)] - stack_vzdalenosti[(stack_delka - 1)] == 0) tedrovnost = 1;

    if (predtimsestup && tedvzestup) vysledek = 0;
    if (bylarovnost && tedrovnost) vysledek = 0;
  }
  return vysledek;
} 

/****************************************************************************/

ring_jen_z_bulkovych_atomu() { /*kvuli vzajemne difuzi staci nadpolovicni vetsina*/
short vysledek = 0;
  for (n = 1; n<=stack_delka; n++) if (bulk[stack_cisla[n]] > 0) vysledek++;
  if (vysledek > stack_delka/2) return 1;
  else return 0;                             
}

/****************************************************************************/
short atom_je_v_aktualnim_ringu(int kdo) {
  short vysledek = 0;
  for (n = 1; n<=stack_delka; n++) if (kdo == stack_cisla[n]) vysledek = 1;
  return vysledek;
}  

/****************************************************************************/

vypocet_rings() {
short nalezen, cislo_ringu = 0;
FILE *frings, *faktualni_ring;
char output[50], prikaz[100];;

for (i=1;i<=max_ringlength;i++) pocet_rings[i] = 0;

frings = fopen("RINGS.txt","w");

for (i = 2; i<=max_distance; i++) {     
  for (j = 1; j<=atoms; j++) {
    for (k = (j+1); k<=atoms; k++) {
      if (useky[j][k] > i) {
        nalezen = 0;
        for (l = 1; l<=koordinace[j]; l++) if ( (useky[kamaradi[j][l]][k] + 1) == i ) nalezen = 1;
        if (nalezen == 1) {useky[j][k] = i; useky[k][j] = i; }
      }
    }
  }
}
  
  /*nyni jsou v tabulce "useky" vsechny vzdalenosti az do pozadovaneho maxima*/        
     
for (i = 1; i<=atoms; i++) {
  stack_cisla[1] = i;
  stack_vzdalenosti[1] = 0;    
  stack_delka = 1;
  uvazovany_kamarad[1] = 0;
    
  while (stack_delka >0) {
    uvazovany_kamarad[stack_delka]++;
    if ( uvazovany_kamarad[stack_delka] > koordinace[stack_cisla[stack_delka]] ) stack_delka--;

/*ted - uvnitr nasledujiciho cyklu else - mam otestovany stack s cislem dosud nepouziteho existujiciho kamarada na konci*/
/*tento kamarad bud uzavre ring, nebo - pokud je jeste sance vzhledem k max_ringlength - se prida jako novy konec stacku*/
/*pri scitani ringu nezapomenout na to, ze krome primitivnich delky 2 maji snahu zapocitat se dvakrat - kvuli tomu je ta 3. podminka*/

    else {
      if (kamaradi[stack_cisla[stack_delka]][uvazovany_kamarad[stack_delka]] == i) {
        if ( (otestuj_ring_SP() > 0) && (otestuj_ring_vbunce() > 0) && (stack_cisla[2] <= stack_cisla[stack_delka]) && ring_jen_z_bulkovych_atomu() ) {
	      if (stack_delka > 2) {
            pocet_rings[stack_delka]++;
            cislo_ringu++;
            fprintf(frings,"Ring %d of length %d \nAtoms : ", cislo_ringu, stack_delka);
            for (j = 1; j<=stack_delka; j++) fprintf(frings," %2d", stack_cisla[j]);
            fprintf(frings," %2d\n", stack_cisla[1]);
            
            if (printrings) {
              sprintf(output,"Ring_%d_%d",cislo_ringu,stack_delka);
              faktualni_ring = fopen(output,"w");
              for (j = 1; j <= atoms; j++) {
                if (druh[j] == 1 &&  atom_je_v_aktualnim_ringu(j)) fprintf(faktualni_ring,"Zn_r ");                
                if (druh[j] == 2 &&  atom_je_v_aktualnim_ringu(j)) fprintf(faktualni_ring,"O_r "); 
                if (druh[j] == 1 && !atom_je_v_aktualnim_ringu(j)) fprintf(faktualni_ring,"Zn "); 
                if (druh[j] == 2 && !atom_je_v_aktualnim_ringu(j)) fprintf(faktualni_ring,"O ");                             
                fprintf(faktualni_ring,"%f %f %f\n",souradnice[j][1],souradnice[j][2],souradnice[j][3]);  
              }                                         
              fclose(faktualni_ring);
              sprintf(prikaz,"move %s Rings",output);
              system(prikaz);
            }              
          }/*konec radovani se z nalezeneho ringu*/
        }         
      }
      
      else {
        if ( (stack_delka < max_ringlength) && (kamaradi[stack_cisla[stack_delka]][uvazovany_kamarad[stack_delka]] > i) ) { 
          stack_cisla[(stack_delka + 1)] = kamaradi[stack_cisla[stack_delka]][uvazovany_kamarad[stack_delka]];   
          stack_delka++; if (stack_delka > 12) printf("*** %d ***\n",stack_delka);
          stack_vzdalenosti[stack_delka] = useky[i][stack_cisla[stack_delka]];
          uvazovany_kamarad[stack_delka] = 0;
          if (unimodal_labelling() == 0) {stack_delka--;}
        }
      }
    }
  }
} /*to byl konec klicoveho cyklu for*/

fclose(frings);
}

/****************************************************************************/

tisk_results() {

fresults = fopen("Results.txt","w");
fprintf(fresults,"Zn atoms %d (%d substrate  %d bulk  %d surface %d other)\n",atoms1,frozen_substrate_1+unfrozen_substrate_1,bulk1,deposited1-bulk1,atoms1-(frozen_substrate_1+unfrozen_substrate_1+deposited1));
fprintf(fresults," O atoms %d (%d substrate  %d bulk  %d surface %d other)\n",atoms2,frozen_substrate_2+unfrozen_substrate_2,bulk2,deposited2-bulk2,atoms2-(frozen_substrate_2+unfrozen_substrate_2+deposited2));

fprintf(fresults,"\nBonding statistics (bulk atoms)\n\n");
fprintf(fresults,"Zn-Zn bonds %d  (average length %.4f)\n",bonds[1][1],bond_length[1][1]);
fprintf(fresults,"Zn-O  bonds %d  (average length %.4f)\n",bonds[1][2],bond_length[1][2]);
fprintf(fresults," O-O  bonds %d  (average length %.4f)\n",bonds[2][2],bond_length[2][2]);

fprintf(fresults,"Zn coordination  %.2f\n",coord[1]);
fprintf(fresults," O coordination  %.2f\n",coord[2]);

fprintf(fresults,"\nRing statistics (per 1 bulk Zn)\n\n");
for(i=4;i<=max_ringlength;i+=2) fprintf(fresults,"length %2d : %.3f rings\n",i, pocet_rings[i]/(bulk1+bulk2+0.0));

fprintf(fresults, "\n %d %d %d %d %.4f %.2f %.2f . . ", bulk1,deposited1-bulk1,bulk2,deposited2-bulk2, bond_length[1][2], coord[1],coord[2]);

for(i=4;i<=max_ringlength;i+=2) fprintf(fresults,"%.3f ", pocet_rings[i]/(bulk1+bulk2+0.0));

fclose(fresults);
}

/**********************************************************************/

tisk_details() {
fresults = fopen("Details.txt","w");               
for (i = 1; i<=atoms; i++) {
  if (druh[i] == 1) fprintf(fresults,"%2s ",prvni);
  if (druh[i] == 2) fprintf(fresults,"%2s ",druhy);   
  fprintf(fresults,"%3d %6f %6f %6f ",i,souradnice[i][1],souradnice[i][2],souradnice[i][3]);
  if (vazan_k_substratu[i]) fprintf(fresults,"vazan ");  else fprintf(fresults,"      ");  
  if (bulk[i])              fprintf(fresults,"bulk ");   else fprintf(fresults,"     ");  
  if (povrchovy_atom[i])    fprintf(fresults,"povrch "); else fprintf(fresults,"       ");  
  fprintf(fresults,"%d kam_nahore ", kamaradu_nahore[i]);
  fprintf(fresults,"\n");
}
fclose(fresults);
}  

/****************************************************************************/

cutoffs() {
cutoff[1][1] = cutoff_11; cutoff[1][2] = cutoff_12;
cutoff[2][1] = cutoff_12; cutoff[2][2] = cutoff_22;

maxcutoff = cutoff_11;
}

/**********************************************************************/

main() {
srand ( 10000*time(NULL) );

cutoffs();
zjisti_stare_informace(); 
vypocti_koordinace();
eliminuj_odprasene(); 
eliminuj_povrchove(); 
bonding_statistics(); 
tisk_rdf();           
tisk_bonding_angles(); /*tady byl koment*/
vypocet_rings();
tisk_results();
tisk_details();

}
