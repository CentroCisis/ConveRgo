/*SPDX-License-Identifier: AGPL-3.0-only*/
Copyright (C) 2024 CISIS - Centro Interregionale per I Sistemi Informatici Geografici e Statistici
This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, version 3.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

// ConveRgo_funzioni.cpp
//

#include <stdlib.h>
#include <stdio.h>
#include <io.h>
#include <math.h>

#include "ConveRgo_funzioni.h"

#define VERS_1 2
#define VERS_2 06


// Generali
double g_piGreco,g_tolDen;

double g_dltLaGreen;

BOOL g_snLettiGrid;

char g_pathIgmGrid[MY_MAX_PATH+3];

// Ellissoidi
double A_WGS,F_WGS,EQ_WGS,EP_WGS,E_WGS,B_WGS,K0;
double A_HAY,F_HAY,EQ_HAY,EP_HAY,E_HAY,B_HAY,EPRIQUA_HAY,C_HAY;

// Bessel
double FI_ROMAH,A_BE,C_BE,B_BE,EQUA_BE,EPRIQUA_BE;
double FI_ROMAB,LA_ROMAB,DIFF_FI,DIFF_ALFA;
double COSTA,COSTB,COSTC,COSTD,COSTE,COSTF,COSTG,COSTH;
double COSTI,COSTL,COSTM,COSTN,COSTO,COSTP,COSTQ;

// Griglie
double g_oriFiDeciP,g_oriLaDeciP; // origin of planimetric grids [deg]
double g_pasFiDeciP,g_pasLaDeciP; // pitches (gaps) of planimetric grids

double g_oriFiDeciH,g_oriLaDeciH; // origin of altimetric grid [deg]
double g_pasFiDeciH,g_pasLaDeciH; // pitches (gaps) of altimetric grid

double g_oriFiDeciK,g_oriLaDeciK; // origin of K grids [deg]
double g_pasFiDeciK,g_pasLaDeciK; // pitches (gaps) of K grids

double matFiWgs[NUM_RIG_P+7][NUM_COL_P+7]; // grid with latitude  differences between Roma40 and Wgs84
double matLaWgs[NUM_RIG_P+7][NUM_COL_P+7]; // grid with longitude differences between Roma40 and Wgs84
double matFiEd[NUM_RIG_P+7][NUM_COL_P+7];  // grid with latitude  differences between Roma40 and Ed50
double matLaEd[NUM_RIG_P+7][NUM_COL_P+7];  // grid with longitude differences between Roma40 and Ed50
char   matMaskP[NUM_RIG_P+7][NUM_COL_P+7]; // mask grid (1,2=filled node, 0=empty)

double matHeight[NUM_RIG_H+3][NUM_COL_H+3]; // grid with height differences (ellissoidic Wgs84 and orthometric)
char   matMaskH[NUM_RIG_H+3][NUM_COL_H+3];  // mask grid (1,2=filled node, 0=empty)

double matFiGrK[NUM_RIG_K+7][NUM_COL_K+7]; // grid with latitude  differences between ETRF89 and ETRF2000
double matLaGrK[NUM_RIG_K+7][NUM_COL_K+7]; // grid with longitude differences between ETRF89 and ETRF2000
double matQuGrK[NUM_RIG_K+7][NUM_COL_K+7]; // grid with h ell  differences between ETRF89 and ETRF2000
int    matMaskK[NUM_RIG_K+7][NUM_COL_K+7]; // mask grid (1=filled node, 0=empty)

// Limitazione geografica
double g_fiMinAssDec,g_fiMaxAssDec,g_laMinAssDec,g_laMaxAssDec;
double g_fiMinAssRad,g_fiMaxAssRad,g_laMinAssRad,g_laMaxAssRad;

// Quanti grigliati
int g_quantiGr1,g_quantiGr2,g_quantiGrA,g_quantiGk1,g_quantiGk2;

// File REP (solo per debug)
BOOL g_snRep;
char g_nomeFileRep[260]; //REP
FILE *g_idflRep; //REP
char strPerRep[260]; //REP


// GRA
class CMagliaGra
{
public:
  double m_maxFiDeci;   // Latitudine dell'origine (max perché riga 0 a nord)
  double m_oriLaDeci;   // Longitudine dell'origine
  double m_pasFiDeci;   // Passo latitudine
  double m_pasLaDeci;   // Passo longitudine
  int    m_numRig;      // Numero righe
  int    m_numCol;      // Numero colonne
  double *m_pMatEnne;   // Puntatore per allocazione matrice N
  char   *m_pMatMask;   // Puntatore per allocazione matrice mask
  bool   m_snAllocato;  // Flag per si/no allocazione
  bool   m_snAssegnato; // Flag per si/no assegnazione

  CMagliaGra(void);
  CMagliaGra(int nRig,int nCol);
  ~CMagliaGra(void);

  bool allocaRigCol(int nRig,int nCol);

}; // CMagliaGra

CMagliaGra::CMagliaGra(void)
{
  m_maxFiDeci=0.0;
  m_oriLaDeci=0.0;
  m_pasFiDeci=0.0;
  m_pasLaDeci=0.0;
  m_numRig=0;
  m_numCol=0;
  m_pMatEnne=NULL;
  m_pMatMask=NULL;
  m_snAllocato=false;
  m_snAssegnato=false;
} // CMagliaGra

CMagliaGra::CMagliaGra(int nRig,int nCol)
{
  m_maxFiDeci=0.0;
  m_oriLaDeci=0.0;
  m_pasFiDeci=0.0;
  m_pasLaDeci=0.0;
  allocaRigCol(nRig,nCol);
} // CMagliaGra

CMagliaGra::~CMagliaGra(void)
{
  if (m_snAllocato) {
    delete [] m_pMatEnne;
    delete [] m_pMatMask;
  } // m_snAllocato
} // ~CMagliaGra

bool CMagliaGra::allocaRigCol(int nRig,int nCol)
{
  int quanti;
  quanti=(nRig+1)*(nCol+1);

  m_numRig=nRig;
  m_numCol=nCol;
  m_pMatEnne=new double[quanti];
  m_pMatMask=new char[quanti];
  m_snAllocato=true;
  m_snAssegnato=false;

  if (m_pMatEnne==NULL) return false;
  if (m_pMatMask==NULL) return false;

  return true;
} // allocaRigCol

CMagliaGra g_aryMaglieGra[NMAX_ARYGRA+3];
int g_numAryGra;
bool g_snGra;


// ******************************************************************************

int ConveRgo_INIZIALIZZA(void)
{
  char pathGrid[520];
  int i,j,quantiGrid;

  g_pathIgmGrid[0]='\0';

  g_piGreco=atan(1.0)*4.0;
  g_tolDen=1.0E-06;

  g_dltLaGreen=37.357*g_piGreco/540.0; // 12°27'08.40" in rad

  K0=0.9996;

  g_snRep=leggiSnRepDaReg();

  strcpy_s(g_nomeFileRep,255,"C:\\ConveRgo_report.txt"); //REP

  if (g_snRep) mettiRep("* ConveRgo_INIZIALIZZA\n");

  // Limitazione geografica (appena interna alla griglia H, per EGM96 mondiale)
  g_fiMinAssDec=35.40;  // griglia H 35°20' 35.33333333°
  g_fiMaxAssDec=47.35;  // griglia H 47°24' 47.40000000°
  g_laMinAssDec= 6.05;  // griglia H  6°00'  6.00000000°
  g_laMaxAssDec=19.05;  // griglia H 19°08' 19.13333333°

  g_laMinAssRad=g_laMinAssDec*g_piGreco/180.0;
  g_laMaxAssRad=g_laMaxAssDec*g_piGreco/180.0;
  g_fiMinAssRad=g_fiMinAssDec*g_piGreco/180.0;
  g_fiMaxAssRad=g_fiMaxAssDec*g_piGreco/180.0;

  // Wgs84
  A_WGS=6378137.0;
  F_WGS=1.0/298.257222101;
  B_WGS=A_WGS-A_WGS*F_WGS;
  EQ_WGS=(A_WGS*A_WGS-B_WGS*B_WGS)/(A_WGS*A_WGS);
  EP_WGS=sqrt(EQ_WGS/(1.0-EQ_WGS));
  E_WGS=sqrt(EQ_WGS);

  // Hayford
  A_HAY=6378388.0;
  F_HAY=1.0/297.0;
  B_HAY=A_HAY-A_HAY*F_HAY;
  C_HAY=A_HAY*A_HAY/B_HAY;
  EQ_HAY=(A_HAY*A_HAY-B_HAY*B_HAY)/(A_HAY*A_HAY);
  EPRIQUA_HAY=EQ_HAY/(1.0-EQ_HAY);
  EP_HAY=sqrt(EPRIQUA_HAY);
  E_HAY=sqrt(EQ_HAY);

  // Griglie
  g_oriFiDeciP=35.0;        // 35°00'00.00"
  g_oriLaDeciP=17.8570/3.0; //  5°57'08.40"
  g_pasFiDeciP=5.0/60.0;    //     5'00"
  g_pasLaDeciP=7.5/60.0;    //     7'30"

  for (i=0;i<NUM_RIG_P;i++) {
    for (j=0;j<NUM_COL_P;j++) {
      matMaskP[i][j]=0;
    } // j
  } // i

  g_oriFiDeciH=106.0/3.0;   // 35°20'00.00"
  g_oriLaDeciH=6.0;         //  6°00'00.00"
  g_pasFiDeciH=0.0333333333333;//1.0/30.0;    //     2'00"
  g_pasLaDeciH=0.0333333333333;//1.0/30.0;    //     2'00"

  for (i=0;i<NUM_RIG_H;i++) {
    for (j=0;j<NUM_COL_H;j++) {
      matMaskH[i][j]=0;
    } // j
  } // i

  g_oriFiDeciK=35.0006;
  g_oriLaDeciK=5.9522;
  g_pasFiDeciK=5.0/60.0;  //     5'00"
  g_pasLaDeciK=7.5/60.0;  //     7'30"

  for (i=0;i<NUM_RIG_K;i++) {
    for (j=0;j<NUM_COL_K;j++) {
      matMaskK[i][j]=0;
    } // j
  } // i

  // Altro
  g_snLettiGrid=false;

  g_numAryGra=0;
  g_snGra=false;

  quantiGrid=0;
  pathGrid[0]='\0';

  g_quantiGr1=0;
  g_quantiGr2=0;
  g_quantiGrA=0;
  g_quantiGk1=0;
  g_quantiGk2=0;

  if (leggiPathGridDaReg(pathGrid)) {
    quantiGrid=ConveRgo_SET_PATH_GRI(pathGrid);
    //if (quantiGrid>0) ConveRgo_LOAD_GRI();
  } // leggiPathGridDaReg

  return 1;

} // ConveRgo_INIZIALIZZA

// ******************************************************************************

double ConveRgo_GESI_DECI(double gesiIn)
{
  int segno,gradi,primi;
  double secondi;

  //if (g_snRep) mettiRep("* ConveRgo_GESI_DECI\n");

  segno=1;
  if (gesiIn<0.0) segno=-1;

  gesiIn=fabs(gesiIn)+1.0E-11;
  gradi=(int) gesiIn;
  primi=(int) ((gesiIn-gradi)*100.0);
  secondi=((gesiIn-gradi)*100.0-primi)/36.0;

  return (segno*(gradi+primi/60.0+secondi));
} // ConveRgo_GESI_DECI

// ******************************************************************************

double ConveRgo_DECI_GESI(double deciIn,int numCifreSec)
{
  int segno,gradi,primi;
  double limite,secondi;

  //if (g_snRep) mettiRep("* ConveRgo_DECI_GESI\n");

  if (numCifreSec<0) return 0.0;

  limite=60.0-5.0/pow(10.0,(numCifreSec+1.0));

  segno=1;
  if (deciIn<0.0) segno=-1;

  deciIn=fabs(deciIn);
  gradi=(int) deciIn;
  primi=(int) ((deciIn-gradi)*60.0);
  secondi=(deciIn-gradi)*3600.0-primi*60.0;

  if (secondi>=limite) {
    secondi=0.0;
    primi++;
  } // secondi>=limite

  if (primi==60) {
    primi=0;
    gradi++;
  } // primi==60

  return (segno*(gradi+primi/100.0+secondi/10000.0));
} // ConveRgo_DECI_GESI

// ******************************************************************************

double ConveRgo_GPD_DECI(double gpdIn)
{
  int segno,gradi;
  double primi;

  if (g_snRep) mettiRep("* ConveRgo_GPD_DECI\n");

  segno=1;
  if (gpdIn<0.0) segno=-1;

  gpdIn=fabs(gpdIn)+1.0E-10;
  gradi=(int) gpdIn;
  primi=(gpdIn-gradi)*100.0;

  return segno*(gradi+primi/60.0);

} // ConveRgo_GPD_DECI

// ******************************************************************************

double ConveRgo_DECI_GPD(double deciIn,int numCifrePri)
{
  int segno,gradi;
  double limite,primi;

  if (g_snRep) mettiRep("* ConveRgo_DECI_GPD\n");

  if (numCifrePri<0) return 0.0;

  limite=60.0-5.0/pow(10.0,(numCifrePri+1.0));

  segno=1;
  if (deciIn<0.0) segno=-1;

  deciIn=fabs(deciIn);
  gradi=(int) deciIn;
  primi=(deciIn-gradi)*60.0;

  if (primi>=limite) {
    primi=0;
    gradi++;
  } // primi==60

  return (segno*(gradi+primi/100.0));

} // ConveRgo_DECI_GPD

// ******************************************************************************

int ConveRgo_GEO_PIA(double fiIn,double laIn,char sistema,char fusoOut,double *N_Out,double *E_Out)
{
  double falsOri,aMag,effe,eQua,eta,dLa,biFi;
  double A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2,ARC1,ARC2,ARC3,ARC4;
  double AX,BX,CX,DX,EX,FX,AY,BY,CY,DY,EY,FY;
  double vuDop,enne,ti,ti2,ti4,ti6,ti8,eta2,eta4,eta6,eta8;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_GEO_PIA in (%.8lf,%.8lf sist=%hi)\n",fiIn,laIn,sistema);
    mettiRep(strPerRep);
  } // g_snRep

  *N_Out=0.0;
  *E_Out=0.0;

  if (fiIn<30.0||fiIn>50.0) return -2;
  if (laIn<5.0||laIn>20.0) return -2;

  falsOri=500000.0;

  if (fusoOut==0) {
    fusoOut=32;
    if (laIn>12.0) fusoOut=33;
    if (laIn>18.0) fusoOut=34;
  } // fuso automatico

  switch (fusoOut) {
  case 32:
    dLa=laIn-9.0;
    if (sistema==ROMA) falsOri=1500000.0;
    break;
  case 33:
    dLa=laIn-15.0;
    if (sistema==ROMA) falsOri=2520000.0;
    break;
  case 34:
    dLa=laIn-21.0;
    if (sistema==ROMA) {
      dLa=laIn-15.0;
      falsOri=2520000.0;
    } // ROMA
    break;
  default:
    return -1;
  } // fusoOut;

  fiIn *= (g_piGreco/180.0);
  dLa *= (g_piGreco/180.0);

  if (sistema==ETRF89||sistema==ETRF2000) {
    aMag=A_WGS*K0;
    effe=F_WGS;
    eQua=EQ_WGS;
    eta=EP_WGS*cos(fiIn);
  } // WGS
  else {
    aMag=A_HAY*K0;
    effe=F_HAY;
    eQua=EQ_HAY;
    eta=EP_HAY*cos(fiIn);
  } // ROMA oppure ED

  vuDop=sqrt(1.0-eQua*sin(fiIn)*sin(fiIn));
  enne=aMag/vuDop;
  ti=tan(fiIn);
  ti2=ti*ti;
  ti4=ti2*ti2;
  ti6=ti4*ti2;
  ti8=ti4*ti4;
  eta2=eta*eta;
  eta4=eta2*eta2;
  eta6=eta4*eta2;
  eta8=eta4*eta4;

  //Bonifacino 1
  A1=(1.0/2.0)*enne*sin(fiIn)*cos(fiIn);
  B1=(1.0/24.0)*enne*sin(fiIn)*pow(cos(fiIn),3.0)*(5.0-ti2+9.0*eta2+4*eta4);
  C1=(1.0/720.0)*enne*sin(fiIn)*pow(cos(fiIn),5.0)*(61.0-58.0*ti2+ti4+270.0*eta2-330.0*eta2*ti2+
    445.0*eta4-680.0*eta4*ti2+324.0*eta6-600.0*eta6*ti2+88.0*eta8-192.0*eta8*ti2);
  D1=(1.0/40320.0)*enne*sin(fiIn)*pow(cos(fiIn),7.0)*(1385.0-3111.0*ti2+543.0*ti4-ti6+
    10899.0*eta2-32802.0*eta2*ti2+9219.0*eta2*ti4+34419.0*eta4-129087.0*eta4*ti2+
    49644.0*eta4*ti4+56385.0*eta6-252084.0*eta6*ti2+121800.0*eta6*ti4+53740.0*eta8-
    310408.0*eta8*ti2-200256.0*eta8*ti4);
  E1=(1.0/3628800.0)*enne*sin(fiIn)*pow(cos(fiIn),9.0)*(50521.0-206276.0*ti2+101166.0*ti4-
    4916.0*ti6+ti8+612540.0*eta2-3277980.0*eta2*ti2+2402100.0*eta2*ti4-239220.0*eta2*ti6+
    3043190.0*eta4-19954380.0*eta4*ti2+19210830.0*eta4*ti4-2879440.0*eta4*ti6);
  F1=(1.0/479001600.0)*enne*sin(fiIn)*pow(cos(fiIn),11.0)*(2702765.0-17460701.0*ti2+
    16535834.0*ti4-2996242.0*ti6+221257.0*ti8-ti8*ti2);

  //Bonifacino 2
  A2=enne*cos(fiIn);
  B2=(1.0/6.0)*enne*pow(cos(fiIn),3.0)*(1.0-ti2+eta2);
  C2=(1.0/120.0)*enne*pow(cos(fiIn),5.0)*(5.0+14.0*eta2-18.0*ti2-58.0*eta2*ti2+ti4+
    13.0*eta4-64.0*eta4*ti2+4.0*eta6-24.0*eta6*ti2);
  D2=(1.0/5040.0)*enne*pow(cos(fiIn),7.0)*(61.0-479.0*ti2+179.0*ti4-ti6+331.0*eta2-
    3298.0*eta2*ti2+1771.0*eta2*ti4+715.0*eta4-8655.0*eta4*ti2+6080.0*eta4*ti4+769.0*eta6-
    10964.0*eta6*ti2+9480.0*eta6*ti4+412.0*eta8-6760.0*eta8*ti2+6912.0*eta8*ti4+88.0*eta8*eta2-
    1632.0*eta8*eta2*ti2+1920.0*eta8*eta2*ti4);
  E2=(1.0/362880.0)*enne*pow(cos(fiIn),9.0)*(1385.0-19028.0*ti2+18270.0*ti4-1636.0*ti6+ti8+
    12284.0*eta2-214140.0*eta2*ti2+290868.0*eta2*ti4-47188.0*eta2*ti6+45318.0*eta4-
    951468.0*eta4*ti2+1652910.0*eta4*ti4-384384.0*eta4*ti6+90804.0*eta6-2220708.0*eta6*ti2+
    4662840.0*eta6*ti4-1394064.0*eta6*ti6);
  F2=(1.0/39916800.0)*enne*pow(cos(fiIn),11.0)*(50521.0-1073517.0*ti2+1949762.0*ti4-
    481250.0*ti6+73749.0*ti8-ti8*ti2+663061.0*eta2-17594876.0*eta2*ti2+43255806.0*eta2*ti4-
    23673524.0*eta2*ti6+1264933.0*eta2*ti8);

  ARC1=1.0-eQua/4.0-3.0*eQua*eQua/64.0-5.0*eQua*eQua*eQua/256.0;
  ARC2=3.0*eQua/8.0+3.0*eQua*eQua/32.0+45.0*eQua*eQua*eQua/1024.0;
  ARC3=15.0*eQua*eQua/256.0+45.0*eQua*eQua*eQua/1024.0;
  ARC4=35.0*eQua*eQua*eQua/3072.0;
  biFi=aMag*(ARC1*fiIn-ARC2*sin(2.0*fiIn)+ARC3*sin(4.0*fiIn)-ARC4*sin(6.0*fiIn));

  //Nord
  AX=A1*pow(dLa,2.0);
  BX=B1*pow(dLa,4.0);
  CX=C1*pow(dLa,6.0);
  DX=D1*pow(dLa,8.0);
  EX=E1*pow(dLa,10.0);
  FX=F1*pow(dLa,12.0);
  *N_Out=biFi+AX+BX+CX+DX+EX+FX;

  //Est
  AY=A2*dLa;
  BY=B2*pow(dLa,3.0);
  CY=C2*pow(dLa,5.0);
  DY=D2*pow(dLa,7.0);
  EY=E2*pow(dLa,9.0);
  FY=F2*pow(dLa,11.0);
  *E_Out=falsOri+(AY+BY+CY+DY+EY+FY);

  if (g_snRep) mettiRep("* ConveRgo_GEO_PIA ok\n");

  return 1;
} // ConveRgo_GEO_PIA

// ******************************************************************************

int ConveRgo_GEO_PIA_U(double fiIn,double laIn,double *N_Out,double *E_Out)
{
  double laCt,fatSca,falsaEst,falsaNord;
  double aMag,effe,eQua,eta,dLa,biFi;
  double A1,B1,C1,D1,E1,F1,A2,B2,C2,D2,E2,F2,ARC1,ARC2,ARC3,ARC4;
  double AX,BX,CX,DX,EX,FX,AY,BY,CY,DY,EY,FY;
  double vuDop,enne,ti,ti2,ti4,ti6,ti8,eta2,eta4,eta6,eta8;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_GEO_PIA_U in (%.8lf,%.8lf)\n",fiIn,laIn);
    mettiRep(strPerRep);
  } // g_snRep

  *N_Out=0.0;
  *E_Out=0.0;

  laCt=12.0;
  fatSca=0.9985;
  falsaEst=7000000.0;
  falsaNord=0.0;

  dLa=laIn-laCt;

  fiIn *= (g_piGreco/180.0);
  dLa *= (g_piGreco/180.0);

  aMag=A_WGS*fatSca;
  effe=F_WGS;
  eQua=EQ_WGS;
  eta=EP_WGS*cos(fiIn);

  vuDop=sqrt(1.0-eQua*sin(fiIn)*sin(fiIn));
  enne=aMag/vuDop;
  ti=tan(fiIn);
  ti2=ti*ti;
  ti4=ti2*ti2;
  ti6=ti4*ti2;
  ti8=ti4*ti4;
  eta2=eta*eta;
  eta4=eta2*eta2;
  eta6=eta4*eta2;
  eta8=eta4*eta4;

  //Bonifacino 1
  A1=(1.0/2.0)*enne*sin(fiIn)*cos(fiIn);
  B1=(1.0/24.0)*enne*sin(fiIn)*pow(cos(fiIn),3.0)*(5.0-ti2+9.0*eta2+4*eta4);
  C1=(1.0/720.0)*enne*sin(fiIn)*pow(cos(fiIn),5.0)*(61.0-58.0*ti2+ti4+270.0*eta2-330.0*eta2*ti2+
    445.0*eta4-680.0*eta4*ti2+324.0*eta6-600.0*eta6*ti2+88.0*eta8-192.0*eta8*ti2);
  D1=(1.0/40320.0)*enne*sin(fiIn)*pow(cos(fiIn),7.0)*(1385.0-3111.0*ti2+543.0*ti4-ti6+
    10899.0*eta2-32802.0*eta2*ti2+9219.0*eta2*ti4+34419.0*eta4-129087.0*eta4*ti2+
    49644.0*eta4*ti4+56385.0*eta6-252084.0*eta6*ti2+121800.0*eta6*ti4+53740.0*eta8-
    310408.0*eta8*ti2-200256.0*eta8*ti4);
  E1=(1.0/3628800.0)*enne*sin(fiIn)*pow(cos(fiIn),9.0)*(50521.0-206276.0*ti2+101166.0*ti4-
    4916.0*ti6+ti8+612540.0*eta2-3277980.0*eta2*ti2+2402100.0*eta2*ti4-239220.0*eta2*ti6+
    3043190.0*eta4-19954380.0*eta4*ti2+19210830.0*eta4*ti4-2879440.0*eta4*ti6);
  F1=(1.0/479001600.0)*enne*sin(fiIn)*pow(cos(fiIn),11.0)*(2702765.0-17460701.0*ti2+
    16535834.0*ti4-2996242.0*ti6+221257.0*ti8-ti8*ti2);

  //Bonifacino 2
  A2=enne*cos(fiIn);
  B2=(1.0/6.0)*enne*pow(cos(fiIn),3.0)*(1.0-ti2+eta2);
  C2=(1.0/120.0)*enne*pow(cos(fiIn),5.0)*(5.0+14.0*eta2-18.0*ti2-58.0*eta2*ti2+ti4+
    13.0*eta4-64.0*eta4*ti2+4.0*eta6-24.0*eta6*ti2);
  D2=(1.0/5040.0)*enne*pow(cos(fiIn),7.0)*(61.0-479.0*ti2+179.0*ti4-ti6+331.0*eta2-
    3298.0*eta2*ti2+1771.0*eta2*ti4+715.0*eta4-8655.0*eta4*ti2+6080.0*eta4*ti4+769.0*eta6-
    10964.0*eta6*ti2+9480.0*eta6*ti4+412.0*eta8-6760.0*eta8*ti2+6912.0*eta8*ti4+88.0*eta8*eta2-
    1632.0*eta8*eta2*ti2+1920.0*eta8*eta2*ti4);
  E2=(1.0/362880.0)*enne*pow(cos(fiIn),9.0)*(1385.0-19028.0*ti2+18270.0*ti4-1636.0*ti6+ti8+
    12284.0*eta2-214140.0*eta2*ti2+290868.0*eta2*ti4-47188.0*eta2*ti6+45318.0*eta4-
    951468.0*eta4*ti2+1652910.0*eta4*ti4-384384.0*eta4*ti6+90804.0*eta6-2220708.0*eta6*ti2+
    4662840.0*eta6*ti4-1394064.0*eta6*ti6);
  F2=(1.0/39916800.0)*enne*pow(cos(fiIn),11.0)*(50521.0-1073517.0*ti2+1949762.0*ti4-
    481250.0*ti6+73749.0*ti8-ti8*ti2+663061.0*eta2-17594876.0*eta2*ti2+43255806.0*eta2*ti4-
    23673524.0*eta2*ti6+1264933.0*eta2*ti8);

  ARC1=1.0-eQua/4.0-3.0*eQua*eQua/64.0-5.0*eQua*eQua*eQua/256.0;
  ARC2=3.0*eQua/8.0+3.0*eQua*eQua/32.0+45.0*eQua*eQua*eQua/1024.0;
  ARC3=15.0*eQua*eQua/256.0+45.0*eQua*eQua*eQua/1024.0;
  ARC4=35.0*eQua*eQua*eQua/3072.0;
  biFi=aMag*(ARC1*fiIn-ARC2*sin(2.0*fiIn)+ARC3*sin(4.0*fiIn)-ARC4*sin(6.0*fiIn));

  //Nord
  AX=A1*pow(dLa,2.0);
  BX=B1*pow(dLa,4.0);
  CX=C1*pow(dLa,6.0);
  DX=D1*pow(dLa,8.0);
  EX=E1*pow(dLa,10.0);
  FX=F1*pow(dLa,12.0);
  *N_Out=falsaNord+biFi+AX+BX+CX+DX+EX+FX;

  //Est
  AY=A2*dLa;
  BY=B2*pow(dLa,3.0);
  CY=C2*pow(dLa,5.0);
  DY=D2*pow(dLa,7.0);
  EY=E2*pow(dLa,9.0);
  FY=F2*pow(dLa,11.0);
  *E_Out=falsaEst+AY+BY+CY+DY+EY+FY;

  if (g_snRep) mettiRep("* ConveRgo_GEO_PIA_U ok\n");

  return 1;

} // ConveRgo_GEO_PIA_U

// ******************************************************************************

int ConveRgo_PIA_GEO(double N_In,double E_In,char fusoIn,char sistema,double *fiOut,double *laOut)
{
  double aMag,effe,eQua,epQua,eta;
  double vuDop,enne,ti,ti2,ti4,ti6,ti8,eta2,eta4,eta6,eta8;
  double mu,emmeP,fiP,laCt,eSuEnneQua,eUno,cosFiP;
  double csi1,csi2,csi3,csi4,csi5,zeta1,zeta2,zeta3,zeta4,zeta5;
  double A1P,B1P,C1P,D1P,E1P,F1P,A2P,B2P,C2P,D2P,E2P,F2P;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_PIA_GEO in (%.3lf,%.3lf fuso=%hi sist=%hi)\n",N_In,E_In,fusoIn,sistema);
    mettiRep(strPerRep);
  } // g_snRep

  *fiOut=0.0;
  *laOut=0.0;

  switch (fusoIn) {
  case 0:
    if (sistema!=ROMA) return -1;
    laCt=9.0; fusoIn=32;
    if (E_In>2.0E+06) { laCt=15.0; fusoIn=33; }
    break;
  case 32:
    laCt=9.0;
    break;
  case 33:
    laCt=15.0;
    break;
  case 34:
    laCt=21.0;
    if (sistema==ROMA) laCt=15.0;
    break;
  default:
    return -1;
  } // fusoIn;

  E_In -= 500000.0;
  if (sistema==ROMA) {
    E_In -= 1000000.0;
    if (fusoIn!=32) E_In -= 1020000.0;
  } // ROMA

  if (sistema==ETRF89||sistema==ETRF2000) {
    aMag=A_WGS*K0;
    effe=F_WGS;
    eQua=EQ_WGS;
    epQua=EP_WGS;
  } // WGS
  else {
    aMag=A_HAY*K0;
    effe=F_HAY;
    eQua=EQ_HAY;
    epQua=EP_HAY;
  } // ROMA oppure ED

  eUno=(1.0-sqrt(1.0-eQua))/(1.0+sqrt(1.0-eQua));
  emmeP=aMag*(1.0-eQua/4.0-3.0*eQua*eQua/64.0-5.0*eQua*eQua*eQua/256.0)*(g_piGreco/2.0);
  mu=(g_piGreco*N_In)/(2.0*emmeP);
  fiP=mu+(3.0*eUno/2.0-27.0*eUno*eUno*eUno/32.0)*sin(2.0*mu)+
    (21.0*eUno*eUno/16.0-55.0*eUno*eUno*eUno*eUno/32.0)*sin(4.0*mu)+
    (151.0*eUno*eUno*eUno/96.0)*sin(6.0*mu)+(1097.0*eUno*eUno*eUno*eUno/512.0)*sin(8.0*mu);

  vuDop=sqrt(1.0-eQua*sin(fiP)*sin(fiP));
  enne=aMag/vuDop;
  eta=epQua*cos(fiP);
  ti=tan(fiP);

  ti2=ti*ti;
  ti4=ti2*ti2;
  ti6=ti4*ti2;
  ti8=ti4*ti4;
  eta2=eta*eta;
  eta4=eta2*eta2;
  eta6=eta4*eta2;
  eta8=eta4*eta4;

  cosFiP=cos(fiP);

  // SVILUPPI BONIFACINO A1P, B1P, C1P, D1P, E1P, F1P
  A1P=-ti*(1.0+eta2)/2.0;
  B1P= (1.0/24.0)*ti*(5.0+3.0*ti2+6.0*eta2-6.0*eta2*ti2-3.0*eta4-9.0*eta4*ti2-4.0*eta6);
  C1P=-(1.0/720.0)*ti*(61.0+90.0*ti2+45.0*ti4+107.0*eta2-162.0*eta2*ti2-45.0*eta2*ti4+
    43.0*eta4-318.0*eta4*ti2+135.0*eta4*ti4+97.0*eta6+18.0*eta6*ti2+225.0*eta6*ti4+
    188.0*eta8-108.0*eta8*ti2+88.0*eta8*eta2-192.0*eta8*eta2*ti2);
  D1P= (1.0/40320.0)*ti*(1385.0+3633.0*ti2+4095.0*ti4+1575.0*ti6+3116.0*eta2-
    5748.0*eta2*ti2-3276.0*eta2*ti4-1260.0*eta2*ti6+1158.0*eta4-17826.0*eta4*ti2+
    12954.0*eta4*ti4+1890.0*eta4*ti6-3500.0*eta6+1164.0*eta6*ti2+30876.0*eta6*ti4-
    6300.0*eta6*ti6-11735.0*eta8+29001.0*eta8*ti2+3519.0*eta8*ti4-11025.0*eta8*ti6);
  E1P=-(1.0/3628800.0)*ti*(50521.0+204180.0*ti2+383670.0*ti4+321300.0*ti6+99225.0*ti8+
    138933.0*eta2+264684.0*eta2*ti2-211410.0*eta2*ti4-192780.0*eta2*ti6-70875.0*eta2*ti8+
    97362.0*eta4-1265040.0*eta4*ti2+777024.0*eta4*ti4+354240.0*eta4*ti6+85050.0*eta4*ti8);
  F1P= (1.0/479001600.0)*ti*(2702765.0+15487263.0*ti2+42660090.0*ti4+
    57900150.0*ti6+38201625.0*ti8+9823275.0*ti8*ti2);

  // SVILUPPI BONIFACINO A2P, B2P, C2P, D2P, E2P, F2P
  A2P= 1.0/cosFiP;
  B2P=-1.0/(6.0*cosFiP)*(1.0+2.0*ti2+eta2);
  C2P= 1.0/(120.0*cosFiP)*(5.0+28.0*ti2+24.0*ti4+6.0*eta2+8.0*eta2*ti2-3.0*eta4+
    4.0*eta4*ti2-4.0*eta6+24.0*eta6*ti2);
  D2P=-1.0/(5040.0*cosFiP)*(61.0+662.0*ti2+1320.0*ti4+720.0*ti6+107.0*eta2+
    440.0*eta2*ti2+336.0*eta2*ti4+43.0*eta4-234.0*eta4*ti2-192.0*eta4*ti4+
    97.0*eta6-772.0*eta6*ti2+408.0*eta6*ti4+188.0*eta8-2392.0*eta8*ti2+
    1536.0*eta8*ti4+88.0*eta8*eta2-1632.0*eta8*eta2*ti2+1920.0*eta8*eta2*ti4);
  E2P= 1.0/(362880.0*cosFiP)*(1385.0+24568.0*ti2+83664.0*ti4+100800.0*ti6+
    40320.0*ti8+3116.0*eta2+26736.0*eta2*ti2+47808.0*eta2*ti4+24192.0*eta2*ti6+
    1158.0*eta4-4944.0*eta4*ti2-20736.0*eta4*ti4-13824.0*eta4*ti6-3500.0*eta6+
    27104.0*eta6*ti2+576.0*eta6*ti4+9216.0*eta6*ti6);
  F2P=-1.0/(39916800.0*cosFiP)*(50521.0+1326122.0*ti2+6749040.0*ti4+13335840.0*ti6+
    11491200.0*ti8+3628800.0*ti8*ti2+138933.0*eta2+2036560.0*eta2*ti2+
    6269472.0*eta2*ti4+7032960.0*eta2*ti6+2661120.0*eta2*ti8);

  eSuEnneQua=(E_In/enne)*(E_In/enne);

  //Calcolo della fi
  csi1=E1P+F1P*eSuEnneQua;
  csi2=D1P+csi1*eSuEnneQua;
  csi3=C1P+csi2*eSuEnneQua;
  csi4=B1P+csi3*eSuEnneQua;
  csi5=A1P+csi4*eSuEnneQua;
  *fiOut=(fiP+csi5*eSuEnneQua)*180.0/g_piGreco;

  //Calcolo della lambda
  zeta1=E2P+F2P*eSuEnneQua;
  zeta2=D2P+zeta1*eSuEnneQua;
  zeta3=C2P+zeta2*eSuEnneQua;
  zeta4=B2P+zeta3*eSuEnneQua;
  zeta5=A2P+zeta4*eSuEnneQua;
  *laOut=laCt+(zeta5*(E_In/enne))*180.0/g_piGreco;

  if (g_snRep) mettiRep("* ConveRgo_PIA_GEO ok\n");

  return 1;
} // ConveRgo_PIA_GEO

// ******************************************************************************

int ConveRgo_PIA_GEO_U(double N_In,double E_In,double *fiOut,double *laOut)
{
  double laCt,fatSca,falsaEst,falsaNord;
  double aMag,effe,eQua,epQua,eta;
  double vuDop,enne,ti,ti2,ti4,ti6,ti8,eta2,eta4,eta6,eta8;
  double mu,emmeP,fiP,eSuEnneQua,eUno,cosFiP;
  double csi1,csi2,csi3,csi4,csi5,zeta1,zeta2,zeta3,zeta4,zeta5;
  double A1P,B1P,C1P,D1P,E1P,F1P,A2P,B2P,C2P,D2P,E2P,F2P;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_PIA_GEO_U in (%.3lf,%.3lf)\n",N_In,E_In);
    mettiRep(strPerRep);
  } // g_snRep

  *fiOut=0.0;
  *laOut=0.0;

  laCt=12.0;
  fatSca=0.9985;
  falsaEst=7000000.0;
  falsaNord=0.0;

  E_In -= falsaEst;
  N_In -= falsaNord;

  aMag=A_WGS*fatSca;
  effe=F_WGS;
  eQua=EQ_WGS;
  epQua=EP_WGS;

  eUno=(1.0-sqrt(1.0-eQua))/(1.0+sqrt(1.0-eQua));
  emmeP=aMag*(1.0-eQua/4.0-3.0*eQua*eQua/64.0-5.0*eQua*eQua*eQua/256.0)*(g_piGreco/2.0);
  mu=(g_piGreco*N_In)/(2.0*emmeP);
  fiP=mu+(3.0*eUno/2.0-27.0*eUno*eUno*eUno/32.0)*sin(2.0*mu)+
    (21.0*eUno*eUno/16.0-55.0*eUno*eUno*eUno*eUno/32.0)*sin(4.0*mu)+
    (151.0*eUno*eUno*eUno/96.0)*sin(6.0*mu)+(1097.0*eUno*eUno*eUno*eUno/512.0)*sin(8.0*mu);

  vuDop=sqrt(1.0-eQua*sin(fiP)*sin(fiP));
  enne=aMag/vuDop;
  eta=epQua*cos(fiP);
  ti=tan(fiP);

  ti2=ti*ti;
  ti4=ti2*ti2;
  ti6=ti4*ti2;
  ti8=ti4*ti4;
  eta2=eta*eta;
  eta4=eta2*eta2;
  eta6=eta4*eta2;
  eta8=eta4*eta4;

  cosFiP=cos(fiP);

  // SVILUPPI BONIFACINO A1P, B1P, C1P, D1P, E1P, F1P
  A1P=-ti*(1.0+eta2)/2.0;
  B1P= (1.0/24.0)*ti*(5.0+3.0*ti2+6.0*eta2-6.0*eta2*ti2-3.0*eta4-9.0*eta4*ti2-4.0*eta6);
  C1P=-(1.0/720.0)*ti*(61.0+90.0*ti2+45.0*ti4+107.0*eta2-162.0*eta2*ti2-45.0*eta2*ti4+
    43.0*eta4-318.0*eta4*ti2+135.0*eta4*ti4+97.0*eta6+18.0*eta6*ti2+225.0*eta6*ti4+
    188.0*eta8-108.0*eta8*ti2+88.0*eta8*eta2-192.0*eta8*eta2*ti2);
  D1P= (1.0/40320.0)*ti*(1385.0+3633.0*ti2+4095.0*ti4+1575.0*ti6+3116.0*eta2-
    5748.0*eta2*ti2-3276.0*eta2*ti4-1260.0*eta2*ti6+1158.0*eta4-17826.0*eta4*ti2+
    12954.0*eta4*ti4+1890.0*eta4*ti6-3500.0*eta6+1164.0*eta6*ti2+30876.0*eta6*ti4-
    6300.0*eta6*ti6-11735.0*eta8+29001.0*eta8*ti2+3519.0*eta8*ti4-11025.0*eta8*ti6);
  E1P=-(1.0/3628800.0)*ti*(50521.0+204180.0*ti2+383670.0*ti4+321300.0*ti6+99225.0*ti8+
    138933.0*eta2+264684.0*eta2*ti2-211410.0*eta2*ti4-192780.0*eta2*ti6-70875.0*eta2*ti8+
    97362.0*eta4-1265040.0*eta4*ti2+777024.0*eta4*ti4+354240.0*eta4*ti6+85050.0*eta4*ti8);
  F1P= (1.0/479001600.0)*ti*(2702765.0+15487263.0*ti2+42660090.0*ti4+
    57900150.0*ti6+38201625.0*ti8+9823275.0*ti8*ti2);

  // SVILUPPI BONIFACINO A2P, B2P, C2P, D2P, E2P, F2P
  A2P= 1.0/cosFiP;
  B2P=-1.0/(6.0*cosFiP)*(1.0+2.0*ti2+eta2);
  C2P= 1.0/(120.0*cosFiP)*(5.0+28.0*ti2+24.0*ti4+6.0*eta2+8.0*eta2*ti2-3.0*eta4+
    4.0*eta4*ti2-4.0*eta6+24.0*eta6*ti2);
  D2P=-1.0/(5040.0*cosFiP)*(61.0+662.0*ti2+1320.0*ti4+720.0*ti6+107.0*eta2+
    440.0*eta2*ti2+336.0*eta2*ti4+43.0*eta4-234.0*eta4*ti2-192.0*eta4*ti4+
    97.0*eta6-772.0*eta6*ti2+408.0*eta6*ti4+188.0*eta8-2392.0*eta8*ti2+
    1536.0*eta8*ti4+88.0*eta8*eta2-1632.0*eta8*eta2*ti2+1920.0*eta8*eta2*ti4);
  E2P= 1.0/(362880.0*cosFiP)*(1385.0+24568.0*ti2+83664.0*ti4+100800.0*ti6+
    40320.0*ti8+3116.0*eta2+26736.0*eta2*ti2+47808.0*eta2*ti4+24192.0*eta2*ti6+
    1158.0*eta4-4944.0*eta4*ti2-20736.0*eta4*ti4-13824.0*eta4*ti6-3500.0*eta6+
    27104.0*eta6*ti2+576.0*eta6*ti4+9216.0*eta6*ti6);
  F2P=-1.0/(39916800.0*cosFiP)*(50521.0+1326122.0*ti2+6749040.0*ti4+13335840.0*ti6+
    11491200.0*ti8+3628800.0*ti8*ti2+138933.0*eta2+2036560.0*eta2*ti2+
    6269472.0*eta2*ti4+7032960.0*eta2*ti6+2661120.0*eta2*ti8);

  eSuEnneQua=(E_In/enne)*(E_In/enne);

  //Calcolo della fi
  csi1=E1P+F1P*eSuEnneQua;
  csi2=D1P+csi1*eSuEnneQua;
  csi3=C1P+csi2*eSuEnneQua;
  csi4=B1P+csi3*eSuEnneQua;
  csi5=A1P+csi4*eSuEnneQua;
  *fiOut=(fiP+csi5*eSuEnneQua)*180.0/g_piGreco;

  //Calcolo della lambda
  zeta1=E2P+F2P*eSuEnneQua;
  zeta2=D2P+zeta1*eSuEnneQua;
  zeta3=C2P+zeta2*eSuEnneQua;
  zeta4=B2P+zeta3*eSuEnneQua;
  zeta5=A2P+zeta4*eSuEnneQua;
  *laOut=laCt+(zeta5*(E_In/enne))*180.0/g_piGreco;

  if (g_snRep) mettiRep("* ConveRgo_PIA_GEO_U ok\n");

  return 1;

} // ConveRgo_PIA_GEO_U

// ******************************************************************************

int ConveRgo_HELL_QSLM(double fiWgs,double laWgs,double hWgs,double *quoOrto)
{
  char cheGeoide;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_HELL_QSLM (%.8lf,%.8lf,%.3lf)\n",fiWgs,laWgs,hWgs);
    mettiRep(strPerRep);
  } // g_snRep

  cheGeoide=0;

  return ConveRgo_HELL_QSLM_G(fiWgs,laWgs,hWgs,cheGeoide,quoOrto);

} // ConveRgo_HELL_QSLM

// ******************************************************************************

int ConveRgo_QSLM_HELL(double fiIn,double laIn,double quoOrto,char sistIn,double *hWgs)
{
  char cheGeoide;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_QSLM_HELL (%.8lf,%.8lf,%.3lf sist=%hi)\n",fiIn,laIn,quoOrto,sistIn);
    mettiRep(strPerRep);
  } // g_snRep

  cheGeoide=0;

  return ConveRgo_QSLM_HELL_G(fiIn,laIn,quoOrto,sistIn,cheGeoide,hWgs);

} // ConveRgo_QSLM_HELL

// ******************************************************************************

int ConveRgo_HELL_QSLM_G(double fiWgs,double laWgs,double hWgs,char cheGeoide,double *quoOrto)
{
  int idxRig,idxCol,idxGra,idxAS,idxAD,idxBS,idxBD,tipoGri;
  double prcFi,prcLa,prcFL,valSep;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_HELL_QSLM (%.8lf,%.8lf,%.3lf geoid=%hi)\n",fiWgs,laWgs,hWgs,cheGeoide);
    mettiRep(strPerRep);
  } // g_snRep

  *quoOrto=0.0;

  tipoGri=0;

  if (cheGeoide<0) return 0;

  if (fiWgs<g_fiMinAssDec||fiWgs>g_fiMaxAssDec||laWgs<g_laMinAssDec||laWgs>g_laMaxAssDec) {
    return -2;
  } // fuori limiti

  if (!g_snLettiGrid) {
    return -1;
  } // !g_snLettiGrid

  if (g_snGra&&(cheGeoide==0||cheGeoide==3)) {
    if (posToIndexGra(fiWgs,laWgs,&idxGra,&idxAS,&idxAD,&idxBS,&idxBD,&prcFi,&prcLa)) {
      // AS=rig,col AD=rig,col+1 BS=rig+1,col BD=rig+1,col+1
      prcFL=prcFi*prcLa;

      *quoOrto=hWgs-(g_aryMaglieGra[idxGra].m_pMatEnne[idxAS]*(1.0-prcFi-prcLa+prcFL)+g_aryMaglieGra[idxGra].m_pMatEnne[idxAD]*(prcLa-prcFL)+
        g_aryMaglieGra[idxGra].m_pMatEnne[idxBS]*(prcFi-prcFL)+g_aryMaglieGra[idxGra].m_pMatEnne[idxBD]*prcFL);

      if (g_snRep) mettiRep("* ConveRgo_HELL_QSLM ok 3\n");

      return 3;
    } // gra
    else if (cheGeoide==3) return -1;
  } // g_snGra

  if (!posToIndexH(fiWgs,laWgs,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

  if (matMaskH[idxRig][idxCol]!=0&&matMaskH[idxRig+1][idxCol]!=0&&
    matMaskH[idxRig][idxCol+1]!=0&&matMaskH[idxRig+1][idxCol+1]!=0) {

    prcFL=prcFi*prcLa;

    *quoOrto=hWgs-(matHeight[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matHeight[idxRig][idxCol+1]*(prcLa-prcFL)+
      matHeight[idxRig+1][idxCol]*(prcFi-prcFL)+matHeight[idxRig+1][idxCol+1]*prcFL);

    if (g_snRep) mettiRep("* ConveRgo_HELL_QSLM ok 2\n");

    return 2;
  } // gri
  else return -1;

} // ConveRgo_HELL_QSLM_G

// ******************************************************************************

int ConveRgo_QSLM_HELL_G(double fiIn,double laIn,double quoOrto,char sistIn,char cheGeoide,double *hWgs)
{
  int idxRig,idxCol,idxGra,idxAS,idxAD,idxBS,idxBD,tipoGri;
  double fiWgs,laWgs,prcFi,prcLa,prcFL,valSep;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_QSLM_HELL (%.8lf,%.8lf,%.3lf sist=%hi geoid=%hi)\n",fiIn,laIn,quoOrto,sistIn,cheGeoide);
    mettiRep(strPerRep);
  } // g_snRep

  *hWgs=0.0;

  tipoGri=0;

  if (cheGeoide<0) return 0;

  if (!g_snLettiGrid) {
    return -1;
  } // !g_snLettiGrid

  if (sistIn==ETRF89||sistIn==ETRF2000) {
    fiWgs=fiIn;
    laWgs=laIn;
  } // WGS
  else if (sistIn==ROMA||sistIn==ED) {
    if (ConveRgo_SISTEMI_GRI(fiIn,laIn,sistIn,ETRF89,&fiWgs,&laWgs)<=0) return 0;
  } // ROMA o ED50
  else return 0;

  if (fiWgs<g_fiMinAssDec||fiWgs>g_fiMaxAssDec||laWgs<g_laMinAssDec||laWgs>g_laMaxAssDec) {
    return -2;
  } // fuori limiti

  if (g_snGra&&(cheGeoide==0||cheGeoide==3)) {
    if (posToIndexGra(fiWgs,laWgs,&idxGra,&idxAS,&idxAD,&idxBS,&idxBD,&prcFi,&prcLa)) {
      // AS=rig,col AD=rig,col+1 BS=rig+1,col BD=rig+1,col+1
      prcFL=prcFi*prcLa;

      *hWgs=quoOrto+(g_aryMaglieGra[idxGra].m_pMatEnne[idxAS]*(1.0-prcFi-prcLa+prcFL)+g_aryMaglieGra[idxGra].m_pMatEnne[idxAD]*(prcLa-prcFL)+
        g_aryMaglieGra[idxGra].m_pMatEnne[idxBS]*(prcFi-prcFL)+g_aryMaglieGra[idxGra].m_pMatEnne[idxBD]*prcFL);

      if (g_snRep) mettiRep("* ConveRgo_QSLM_HELL ok 3\n");

      return 3;
    } // gra
    else if (cheGeoide==3) return -1;
  } // g_snGra

  if (!posToIndexH(fiWgs,laWgs,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

  if (matMaskH[idxRig][idxCol]!=0&&matMaskH[idxRig+1][idxCol]!=0&&
    matMaskH[idxRig][idxCol+1]!=0&&matMaskH[idxRig+1][idxCol+1]!=0) {

    prcFL=prcFi*prcLa;

    *hWgs=quoOrto+(matHeight[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matHeight[idxRig][idxCol+1]*(prcLa-prcFL)+
      matHeight[idxRig+1][idxCol]*(prcFi-prcFL)+matHeight[idxRig+1][idxCol+1]*prcFL);

    if (g_snRep) mettiRep("* ConveRgo_QSLM_HELL ok 2\n");

    return 2;
  } // gri
  else return -1;

} // ConveRgo_QSLM_HELL_G

// ******************************************************************************

BOOL posToIndexP(double fi,double la,int *rig,int *col)
{
  // arrotonda al nodo più vicino

  *rig=0;
  *col=0;

  *rig=(int) ((fi-g_oriFiDeciP)/g_pasFiDeciP+0.5);
  *col=(int) ((la-g_oriLaDeciP)/g_pasLaDeciP+0.5);

  if (*rig<0||*rig>NUM_RIG_P) return false;
  if (*col<0||*col>NUM_COL_P) return false;

  return true;
} // posToIndexP

// ***************************************************************************

BOOL posToIndexP(double fi,double la,int *rig,int *col,double *prcFi,double *prcLa)
{
  double quantiFi,quantiLa;

  *rig=0;
  *col=0;
  *prcFi=0.0;
  *prcLa=0.0;

  quantiFi=(fi-g_oriFiDeciP)/g_pasFiDeciP;
  quantiLa=(la-g_oriLaDeciP)/g_pasLaDeciP;

  *rig=(int) quantiFi;
  *col=(int) quantiLa;

  if (*rig<0||*rig>NUM_RIG_P) return false;
  if (*col<0||*col>NUM_COL_P) return false;

  *prcFi = quantiFi - *rig;
  *prcLa = quantiLa - *col;

  return true;
} // posToIndexP

// ***************************************************************************

BOOL posToIndexH(double fi,double la,int *rig,int *col)
{
  // arrotonda al nodo più vicino

  *rig=0;
  *col=0;

  *rig=(int) ((fi-g_oriFiDeciH)/g_pasFiDeciH+0.5);
  *col=(int) ((la-g_oriLaDeciH)/g_pasLaDeciH+0.5);

  if (*rig<0||*rig>NUM_RIG_H) return false;
  if (*col<0||*col>NUM_COL_H) return false;

  return true;
} // posToIndexH

// ***************************************************************************

BOOL posToIndexH(double fi,double la,int *rig,int *col,double *prcFi,double *prcLa)
{
  double quantiFi,quantiLa;

  *rig=0;
  *col=0;
  *prcFi=0.0;
  *prcLa=0.0;

  quantiFi=(fi-g_oriFiDeciH)/g_pasFiDeciH;
  quantiLa=(la-g_oriLaDeciH)/g_pasLaDeciH;

  *rig=(int) quantiFi;
  *col=(int) quantiLa;

  if (*rig<0||*rig>NUM_RIG_H) return false;
  if (*col<0||*col>NUM_COL_H) return false;

  *prcFi = quantiFi - *rig;
  *prcLa = quantiLa - *col;

  return true;
} // posToIndexH

// ***************************************************************************

BOOL posToIndexK(double fi,double la,int *rig,int *col)
{
  // Parametri:
  // - fi,la = coordinate geografiche sessadecimali del punto
  // - rig,col = indici della matrice

  // Calcola nelle matrici planimetriche gli indici (riga, colonna)
  // del nodo più vicino a dove cade la posizione del punto

  *rig=0;
  *col=0;

  *rig=(int) ((fi-g_oriFiDeciK)/g_pasFiDeciK+0.5);
  *col=(int) ((la-g_oriLaDeciK)/g_pasLaDeciK+0.5);

  if (*rig<0||*rig>NUM_RIG_K) return false;
  if (*col<0||*col>NUM_COL_K) return false;

  return true;
} // posToIndexK

// ***************************************************************************

BOOL posToIndexK(double fi,double la,int *rig,int *col,double *prcFi,double *prcLa)
{
  // Parametri:
  // - fi,la = coordinate geografiche sessadecimali del punto
  // - rig,col = indici della matrice
  // - prcFi,prcLa = distanza in latitudine e longitudine fra la
  //   posizione del punto e il nodo della matrice

  // Calcola nelle matrici planimetriche gli indici (riga, colonna)
  // del nodo più vicino a dove cade la posizione del punto e la quantità
  // residua, cioè quanto il punto sia spostato rispetto al nodo.

  double quantiFi,quantiLa; // quante righe e quante colonne per arrivare al punto (valori non interi)

  *rig=0;
  *col=0;
  *prcFi=0.0;
  *prcLa=0.0;

  quantiFi=(fi-g_oriFiDeciK)/g_pasFiDeciK;
  quantiLa=(la-g_oriLaDeciK)/g_pasLaDeciK;

  *rig=(int) quantiFi;
  *col=(int) quantiLa;

  if (*rig<0||*rig>NUM_RIG_K) return false;
  if (*col<0||*col>NUM_COL_K) return false;

  *prcFi = quantiFi - *rig;
  *prcLa = quantiLa - *col;

  return true;
} // posToIndexK

// ***************************************************************************

int ConveRgo_SET_PATH_GRI(char grPath[520])
{
  char pathPerFind[MY_MAX_PATH+3];
  int k,quanti;
  long int findId;
  _finddata_t f_t;

  if (g_snRep) mettiRep("* ConveRgo_SET_PATH_GRI in\n");
  if (g_snRep) mettiRep(grPath); // occhio che è 260
  if (g_snRep) mettiRep("\n");

  k=(int) strlen(grPath);
  if (k<1) return 0;
  if (grPath[k-1]!='\\') strcat_s(grPath,512,"\\");

  quanti=0;

  // GR1
  strcpy_s(pathPerFind,520,grPath);
  strcat_s(pathPerFind,520,"*.GR1");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    quanti++;
    while (_findnext(findId,&f_t)==0) {
      quanti++;
    } // while
  } // findId

  _findclose(findId);

  // GR2
  strcpy_s(pathPerFind,520,grPath);
  strcat_s(pathPerFind,520,"*.GR2");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    quanti++;
    while (_findnext(findId,&f_t)==0) {
      quanti++;
    } // while
  } // findId

  _findclose(findId);

  // GRA
  strcpy_s(pathPerFind,520,grPath);
  strcat_s(pathPerFind,520,"*.GRA");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    quanti++;
    while (_findnext(findId,&f_t)==0) {
      quanti++;
    } // while
  } // findId

  _findclose(findId);

  // GK1
  strcpy_s(pathPerFind,520,grPath);
  strcat_s(pathPerFind,520,"*.GK1");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    quanti++;
    while (_findnext(findId,&f_t)==0) {
      quanti++;
    } // while
  } // findId

  _findclose(findId);

  // GK2
  strcpy_s(pathPerFind,520,grPath);
  strcat_s(pathPerFind,520,"*.GK2");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    quanti++;
    while (_findnext(findId,&f_t)==0) {
      quanti++;
    } // while
  } // findId

  _findclose(findId);

  //if (quanti<1) return 0;

  strcpy_s(g_pathIgmGrid,520,grPath);

  k=(int) strlen(g_pathIgmGrid);
  if (g_pathIgmGrid[k-1]!='\\') strcat_s(g_pathIgmGrid,520,"\\");

  if (g_snRep) mettiRep("%i - trovati %i file \n",1,quanti);

  if (g_snRep) mettiRep("* ConveRgo_SET_PATH_GRI ok\n");

  return quanti;

} // ConveRgo_SET_PATH_GRI

// ***************************************************************************

int ConveRgo_LOAD_GRI(void)
{
  char strNomeFile[MY_MAX_PATH+3],strSoloNome[MY_MAX_PATH+3],pathPerFind[MY_MAX_PATH+3];
  _finddata_t f_t;
  long int findId;
  BOOL snTro;

  // -1 = percorso non impostato
  // -2 = non trovati file *.gr? nel folder

  if (g_snRep) mettiRep("* ConveRgo_LOAD_GRI\n");

  if (strlen(g_pathIgmGrid)<1) return -1;

  g_quantiGr1=0;
  g_quantiGr2=0;
  g_quantiGrA=0;
  g_quantiGk1=0;
  g_quantiGk2=0;

  // GR1
  strcpy_s(pathPerFind,520,g_pathIgmGrid);
  strcat_s(pathPerFind,520,"*.GR1");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    snTro=true;
    while (snTro) {

      strncpy_s(strSoloNome,520,f_t.name,40);

      strcpy_s(strNomeFile,520,g_pathIgmGrid);
      strNomeFile[MY_MAX_PATH-40]='\0';
      strcat_s(strNomeFile,520,strSoloNome);

      if (g_snRep) mettiRep("Trovato file:  ");
      if (g_snRep) mettiRep(strNomeFile);
      if (g_snRep) mettiRep("\n");

      if (leggiFileIgmGr12(strNomeFile,strSoloNome,1)) g_quantiGr1++;

      snTro=false;
      if (_findnext(findId,&f_t)==0) snTro=true;
    } // while
  } // findId

  _findclose(findId);

  // GR2
  strcpy_s(pathPerFind,520,g_pathIgmGrid);
  strcat_s(pathPerFind,520,"*.GR2");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    snTro=true;
    while (snTro) {

      strncpy_s(strSoloNome,520,f_t.name,40);

      strcpy_s(strNomeFile,520,g_pathIgmGrid);
      strNomeFile[MY_MAX_PATH-40]='\0';
      strcat_s(strNomeFile,520,strSoloNome);

      if (g_snRep) mettiRep("Trovato file:  ");
      if (g_snRep) mettiRep(strNomeFile);
      if (g_snRep) mettiRep("\n");

      if (leggiFileIgmGr12(strNomeFile,strSoloNome,2)) g_quantiGr2++;

      snTro=false;
      if (_findnext(findId,&f_t)==0) snTro=true;
    } // while
  } // findId

  _findclose(findId);

  // GRA
  strcpy_s(pathPerFind,520,g_pathIgmGrid);
  strcat_s(pathPerFind,520,"*.GRA");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    snTro=true;
    while (snTro) {

      strncpy_s(strSoloNome,520,f_t.name,40);

      strcpy_s(strNomeFile,520,g_pathIgmGrid);
      strNomeFile[MY_MAX_PATH-40]='\0';
      strcat_s(strNomeFile,520,strSoloNome);

      if (g_snRep) mettiRep("Trovato file:  ");
      if (g_snRep) mettiRep(strNomeFile);
      if (g_snRep) mettiRep("\n");

      if (leggiFileGra(strNomeFile)) g_quantiGrA++;

      snTro=false;
      if (_findnext(findId,&f_t)==0) snTro=true;
    } // while
  } // findId

  _findclose(findId);

  // GK1
  strcpy_s(pathPerFind,520,g_pathIgmGrid);
  strcat_s(pathPerFind,520,"*.GK1");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    snTro=true;
    while (snTro) {

      strncpy_s(strSoloNome,520,f_t.name,40);

      strcpy_s(strNomeFile,520,g_pathIgmGrid);
      strNomeFile[MY_MAX_PATH-40]='\0';
      strcat_s(strNomeFile,520,strSoloNome);

      if (g_snRep) mettiRep("Trovato file:  ");
      if (g_snRep) mettiRep(strNomeFile);
      if (g_snRep) mettiRep("\n");

      if (leggiFileIgmGrdK(strNomeFile,strSoloNome,1)) g_quantiGk1++;

      snTro=false;
      if (_findnext(findId,&f_t)==0) snTro=true;
    } // while
  } // findId

  _findclose(findId);

  // GK2
  strcpy_s(pathPerFind,520,g_pathIgmGrid);
  strcat_s(pathPerFind,520,"*.GK2");

  findId=_findfirst(pathPerFind,&f_t);

  if (findId>=0) {
    snTro=true;
    while (snTro) {

      strncpy_s(strSoloNome,520,f_t.name,40);

      strcpy_s(strNomeFile,520,g_pathIgmGrid);
      strNomeFile[MY_MAX_PATH-40]='\0';
      strcat_s(strNomeFile,520,strSoloNome);

      if (g_snRep) mettiRep("Trovato file:  ");
      if (g_snRep) mettiRep(strNomeFile);
      if (g_snRep) mettiRep("\n");

      if (leggiFileIgmGrdK(strNomeFile,strSoloNome,2)) g_quantiGk2++;

      snTro=false;
      if (_findnext(findId,&f_t)==0) snTro=true;
    } // while
  } // findId

  _findclose(findId);

  if (g_quantiGr1==0&&g_quantiGr2==0&&g_quantiGrA==0&&g_quantiGk1==0&&g_quantiGk2==0) return -2;

  g_snLettiGrid=true;

  if (g_numAryGra>0) g_snGra=true;

  if (g_snRep) mettiRep("%i - caricati %i file \n",2,g_quantiGr1+g_quantiGr2+g_quantiGrA+g_quantiGk1+g_quantiGk2);

  return 1;

} // ConveRgo_LOAD_GRI

// ******************************************************************************

int ConveRgo_QUANTI_GRI(int *quantiGr1,int *quantiGr2,int *quantiGrA)
{
  if (g_snRep) mettiRep("* ConveRgo_QUANTI_GRI\n");

  *quantiGr1=g_quantiGr1;
  *quantiGr2=g_quantiGr2;
  *quantiGrA=g_quantiGrA;

  return g_quantiGr1+g_quantiGr2+g_quantiGrA;

} // ConveRgo_QUANTI_GRI

// ******************************************************************************

int ConveRgo_QUANTI_GK(int *quantiGk1,int *quantiGk2)
{
  if (g_snRep) mettiRep("* ConveRgo_QUANTI_GK\n");

  *quantiGk1=g_quantiGk1;
  *quantiGk2=g_quantiGk2;

  return g_quantiGk1+g_quantiGk2;

} // ConveRgo_QUANTI_GRI

// ******************************************************************************

int ConveRgo_TEST_GRI(char sistema,double laMin,double fiMin,double laMax,double fiMax)
{
  FILE *idflDxf;
  char recs[260];
  int i,j,idxRigMin,idxColMin,idxRigMax,idxColMax;
  long int quantiPieni;
  unsigned long int valHdl;
  double prcFi,prcLa,fiNod,laNod,valZ;
  BOOL snVuoti;

  if (g_snRep) mettiRep("* ConveRgo_TEST_GRI\n");

  quantiPieni=0l;
  snVuoti=false;

  if (!posToIndexP(fiMin,laMin,&idxRigMin,&idxColMin,&prcFi,&prcLa)) return 0;

  if (!posToIndexP(fiMax,laMax,&idxRigMax,&idxColMax,&prcFi,&prcLa)) return 0;
  if (prcLa>g_tolDen) idxColMax++;
  if (prcFi>g_tolDen) idxRigMax++;

  for (i=idxRigMin;i<=idxRigMax;i++) {
    for (j=idxColMin;j<=idxColMax;j++) {
      if (matMaskP[i][j]>0) quantiPieni++;
      else snVuoti=true;
    } // j
  } // i

  if (quantiPieni==0l) return -2;

  if (snVuoti) return -1;

  return 1;

} // ConveRgo_TEST_GRI

// ******************************************************************************

BOOL leggiFileIgmGr12(char strFileGrid[MY_MAX_PATH+3],char strSoloNome[43],int cheTipoR)
{
  // Parametri:
  // - strFileGrid = nome completo del file
  // - strSoloNome = nome del file senza percorso

  // Legge il contenuto del file ricevuto come parametro
  // assegnando i valori contenuti nel file agli appositi elementi delle matrici
  // (parcheggia temporaneamente i valori dentro ad array locali, poi li copia
  // nelle matrici vere).

  // Usa le funzioni:
  // - mettiZero(stringa) che inserice un carattere zero se manca (es. .123 divente 0.123)
  // - GESI_DECI(valore) che passa gli angoli da sessagesimali a sessadecimali
  // - posToIndexP e posToIndexH che calcolano gli indici delle matrici che corrispondono
  //   alle coordinate del punto (P matrici planimetriche, H matrice altimetrica)

  // Fa uso di vecchie funzioni C (fgets, sscanf, atof) per motivi di compatibilità

  FILE *idflGr;               // vecchia struttura FILE (per motivi di riutilizzo)
  char recs[260];             // stringa per i record del file
  char strErr[260];           // stringa per i messaggi di errore
  int i,j;                    // di servizio
  int nriga;                  // numero di riga del file
  int rigZeroPla,colZeroPla;  // indici nelle matrici planimetriche dell'origine dell'area del file
  int rigZeroQuo,colZeroQuo;  // indici nella matrie altimetrica dell'origine dell'area del file
  int codErr;                 // codice di errore
  double doppia;
  double fiZeroPla,laZeroPla; // coordinate dell'origine dell'area planimetrica del file
  double fiZeroQuo,laZeroQuo; // coordinate dell'origine dell'area altimetrica del file
  double aryFiEd[40],aryLaEd[40];   // array di servizio per le coordinate ED50
  double aryFiWgs[40],aryLaWgs[40]; // array di servizio per le coordinate WGS84
  double aryQuo[150];               // array di servizio per le quote
  int numRigPla,numColPla,numRigQuo,numColQuo,quantiNodiPla,quantiNodiQuo,idxNodoZeroQuo;
  BOOL snTipoPunto;

  codErr=111;

  snTipoPunto=false;
  numRigPla=6;
  numColPla=6;
  numRigQuo=10;
  numColQuo=14;
  quantiNodiPla=36;
  quantiNodiQuo=140;
  idxNodoZeroQuo=126;

  // apertura file
  idflGr=NULL;
  idflGr=fopen(strFileGrid,"r");
  if (idflGr==NULL) {
    if (g_snRep) mettiRep("*** Errore nell'apertura del file\n");
    return false;
  } // !Open

  // lettura e controllo dell'intestazione
  nriga=1;
  if (fgets(recs,250,idflGr)==NULL) {
    if (idflGr!=NULL) fclose(idflGr);
    if (g_snRep) mettiRep("*** Errore nell'accesso al file\n");
    return false;
  } // !ReadString
  _strupr_s(recs,255);

  if (strncmp(recs,"FOGLIO",6)==0||strncmp(recs+1,"FOGLIO",6)==0) {
    snTipoPunto=false;
    numRigPla=6;
    numColPla=6;
    numRigQuo=10;
    numColQuo=14;
    quantiNodiPla=36;   // numRigPla*numColPla
    quantiNodiQuo=140;  // numRigQuo*numColQuo
    idxNodoZeroQuo=126; // (numRigQuo-1)*numColQuo
  } // foglio
  else if (strncmp(recs,"PUNTO",5)==0||strncmp(recs+1,"PUNTO",5)==0) {
    snTipoPunto=true;
    numRigPla=4;
    numColPla=4;
    numRigQuo=8;
    numColQuo=10;
    quantiNodiPla=16;  // numRigPla*numColPla
    quantiNodiQuo=80;  // numRigQuo*numColQuo
    idxNodoZeroQuo=70; // (numRigQuo-1)*numColQuo
  } // punto
  else {
    if (idflGr!=NULL) fclose(idflGr);
    if (g_snRep) mettiRep("*** Errore nel formato del file\n");
    return false;
  } // err

  // data ED50
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=112; goto LEGGIGRID_ERR; }

  // lati ED50
  for (i=0;i<quantiNodiPla;i++) { // era 36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=113; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=114; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryFiEd[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // longi ED50
  for (i=0;i<quantiNodiPla;i++) { // era 36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=115; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=116; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryLaEd[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // data WGS84
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=117; goto LEGGIGRID_ERR; }

  // lati WGS84
  for (i=0;i<quantiNodiPla;i++) { // era 36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=118; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=119; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryFiWgs[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // longi WGS84
  for (i=0;i<quantiNodiPla;i++) { // era 36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=120; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=121; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryLaWgs[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // lati ROMA40 nodo 1
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=122; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=123; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  fiZeroPla=doppia;

  // longi ROMA40 nodo 1
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=124; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=125; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  laZeroPla=doppia;

  // data geoide
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=126; goto LEGGIGRID_ERR; }

  // geoide
  for (i=0;i<quantiNodiQuo;i++) { // era 140
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=127; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=128; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryQuo[i]=doppia;
  } // i

  // lati ROMA40 nodo 127 (idx 126) quote - se tipo punto => nodo 73 (idx 72) - è idxNodoZeroQuo
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=129; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=130; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  fiZeroQuo=doppia;

  // longi ROMA40 nodo 127 (idx 126) quote - se tipo punto => nodo 73 (idx 72) - è idxNodoZeroQuo
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=131; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=132; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  laZeroQuo=doppia;

  // chiusura
  if (idflGr!=NULL) { fclose(idflGr); idflGr=NULL; }

  // assegnazione matrice
  if (!posToIndexP(fiZeroPla,laZeroPla,&rigZeroPla,&colZeroPla)) {
    codErr=133;
    goto LEGGIGRID_ERR;
  } // err

  for (i=0;i<numRigPla;i++) { // era i<6
    for (j=0;j<numColPla;j++) { // era j<6
      matFiEd[rigZeroPla+i][colZeroPla+j]=aryFiEd[i*numColPla+j];   // era i*6
      matLaEd[rigZeroPla+i][colZeroPla+j]=aryLaEd[i*numColPla+j];   // era i*6
      matFiWgs[rigZeroPla+i][colZeroPla+j]=aryFiWgs[i*numColPla+j]; // era i*6
      matLaWgs[rigZeroPla+i][colZeroPla+j]=aryLaWgs[i*numColPla+j]; // era i*6
      matMaskP[rigZeroPla+i][colZeroPla+j]=cheTipoR;
    } // j
  } // i

  // assegnazione matrice quote
  if (!posToIndexH(fiZeroQuo,laZeroQuo,&rigZeroQuo,&colZeroQuo)) {
    codErr=134;
    goto LEGGIGRID_ERR;
  } // err

  for (i=0;i<numRigQuo;i++) { // era i<10
    for (j=0;j<numColQuo;j++) { // era j<14
      matHeight[rigZeroQuo+i][colZeroQuo+j]=(float) aryQuo[idxNodoZeroQuo-i*numColQuo+j]; // era aryQuo[126-i*14+j]
      matMaskH[rigZeroQuo+i][colZeroQuo+j]=cheTipoR;
    } // j
  } // i

  if (g_snRep) mettiRep("Caricato file: ");
  if (g_snRep) mettiRep(strFileGrid); // occhio che è 260
  if (g_snRep) mettiRep("\n");

  return true;

LEGGIGRID_ERR:

  if (idflGr!=NULL) fclose(idflGr);

  sprintf_s(strErr,255,"*** Errore %i nella lettura del file\n",codErr);
  if (g_snRep) mettiRep(strErr);

  recs[255]='\0';
  if (g_snRep) mettiRep(">");
  if (g_snRep) mettiRep(recs);

  return false;

} // leggiFileIgmGr12

// ***************************************************************************

BOOL leggiFileIgmGrdK(char strFileGrid[MY_MAX_PATH+3],char strSoloNome[43],int cheTipoK)
{
  // Parametri:
  // - strFileGrid = nome completo del file
  // - strSoloNome = nome del file senza percorso

  // Legge il contenuto del file ricevuto come parametro
  // assegnando i valori contenuti nel file agli appositi elementi delle matrici
  // (parcheggia temporaneamente i valori dentro ad array locali, poi li copia
  // nelle matrici vere).

  // Usa le funzioni:
  // - mettiZero(stringa) che inserice un carattere zero se manca (es. .123 divente 0.123)
  // - GESI_DECI(valore) che passa gli angoli da sessagesimali a sessadecimali
  // - posToIndexP e posToIndexH che calcolano gli indici delle matrici che corrispondono
  //   alle coordinate del punto (P matrici planimetriche, H matrice altimetrica)

  // Fa uso di vecchie funzioni C (fgets, sscanf, atof) per motivi di compatibilità

  FILE *idflGr;               // vecchia struttura FILE (per motivi di riutilizzo)
  char recs[260];             // stringa per i record del file
  char strErr[260];           // stringa per i messaggi di errore
  int i,j;                    // di servizio
  int nriga;                  // numero di riga del file
  int rigZeroPla,colZeroPla;  // indici nelle matrici planimetriche dell'origine dell'area del file
  int rigZeroQuo,colZeroQuo;  // indici nella matrie altimetrica dell'origine dell'area del file
  int rigZeroGrK,colZeroGrK;  // indici nelle matrici plano-altim. dell'origine dell'area del file per le griglie K
  int codErr;                 // codice di errore
  double doppia;
  double fiZeroPla,laZeroPla; // coordinate dell'origine dell'area planimetrica del file
  double fiZeroQuo,laZeroQuo; // coordinate dell'origine dell'area altimetrica del file
  double fiZeroGrK,laZeroGrK; // coordinate dell'origine dell'area plano-altim. del file per le griglie K
  double aryFiEd[40],aryLaEd[40];   // array di servizio per le coordinate ED50
  double aryFiWgs[40],aryLaWgs[40]; // array di servizio per le coordinate ETRF89
  double aryFiGrK[40],aryLaGrK[40],aryQuGrK[40]; // array di servizio per le coordinate ETRF2000
  double aryQuo[150];               // array di servizio per le quote
  int numRigPla,numColPla,numRigQuo,numColQuo,quantiNodiPla,quantiNodiQuo,idxNodoZeroQuo;
  BOOL snTipoPunto;

  codErr=211;

  snTipoPunto=false;
  numRigPla=6;
  numColPla=6;
  numRigQuo=10;
  numColQuo=14;
  quantiNodiPla=36;
  quantiNodiQuo=140;
  idxNodoZeroQuo=126;

  // apertura file
  idflGr=NULL;
  idflGr=fopen(strFileGrid,"r");
  if (idflGr==NULL) {
    if (g_snRep) mettiRep("*** Errore nell'apertura del file\n");
    return false;
  } // !Open

  // lettura e controllo dell'intestazione
  nriga=1;
  if (fgets(recs,250,idflGr)==NULL) {
    if (g_snRep) mettiRep("*** Errore nell'accesso al file\n");
    if (idflGr!=NULL) fclose(idflGr);
    return false;
  } // !ReadString
  _strupr_s(recs,255);

  if (strncmp(recs,"FOGLIO",6)==0||strncmp(recs+1,"FOGLIO",6)==0) {
    snTipoPunto=false;
    numRigPla=6;
    numColPla=6;
    numRigQuo=10;
    numColQuo=14;
    quantiNodiPla=36;   // numRigPla*numColPla
    quantiNodiQuo=140;  // numRigQuo*numColQuo
    idxNodoZeroQuo=126; // (numRigQuo-1)*numColQuo
  } // foglio
  else if (strncmp(recs,"PUNTO",5)==0||strncmp(recs+1,"PUNTO",5)==0) {
    snTipoPunto=true;
    numRigPla=4;
    numColPla=4;
    numRigQuo=8;
    numColQuo=10;
    quantiNodiPla=16;  // numRigPla*numColPla
    quantiNodiQuo=80;  // numRigQuo*numColQuo
    idxNodoZeroQuo=70; // (numRigQuo-1)*numColQuo
  } // punto
  else {
    if (idflGr!=NULL) fclose(idflGr);
    if (g_snRep) mettiRep("*** Errore nel formato del file\n");
    return false;
  } // err

  // K
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=212; goto LEGGIGRID_ERR; }
  _strupr_s(recs,255);
  if (strncmp(recs,"K",1)!=0&&strncmp(recs+1,"K",1)!=0)  { codErr=213; goto LEGGIGRID_ERR; }

  // ED50-ROMA40
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=214; goto LEGGIGRID_ERR; }
  _strupr_s(recs,255);
  if (strncmp(recs,"ED50-ROMA40",11)!=0&&strncmp(recs+1,"ED50-ROMA40",11)!=0)  { codErr=215; goto LEGGIGRID_ERR; }

  // data ED50
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=216; goto LEGGIGRID_ERR; }

  // lati ED50
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=217; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=218; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryFiEd[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // longi ED50
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=219; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=220; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryLaEd[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // ETRF89-ROMA40
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=221; goto LEGGIGRID_ERR; }
  _strupr_s(recs,255);
  if (strncmp(recs,"ETRF89-ROMA40",13)!=0&&strncmp(recs+1,"ETRF89-ROMA40",13)!=0)  { codErr=222; goto LEGGIGRID_ERR; }

  // data WGS84
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=223; goto LEGGIGRID_ERR; }

  // lati WGS84
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=224; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=225; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryFiWgs[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // longi WGS84
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=226; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=227; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryLaWgs[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // lati ROMA40 nodo 1
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=228; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=229; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  fiZeroPla=doppia;

  // longi ROMA40 nodo 1
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=230; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=231; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  laZeroPla=doppia;

  // ETRF89-ETRF2000
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=232; goto LEGGIGRID_ERR; }
  _strupr_s(recs,255);
  if (strncmp(recs,"ETRF89-ETRF2000",15)!=0&&strncmp(recs+1,"ETRF89-ETRF2000",15)!=0)  { codErr=233; goto LEGGIGRID_ERR; }

  // data etrf2000
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=234; goto LEGGIGRID_ERR; }

  // lati ETRF2000
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=235; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=236; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryFiGrK[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // longi ETRF2000
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=237; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=238; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryLaGrK[i]=ConveRgo_GESI_DECI(doppia/10000.0);
  } // i

  // quota ETRF2000
  for (i=0;i<quantiNodiPla;i++) { // era i<36
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=239; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=240; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryQuGrK[i]=doppia;
  } // i

  // lati ETRF89 nodo 1
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=241; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=242; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  fiZeroGrK=doppia;

  // longi ETRF89 nodo 1
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=243; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=244; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  laZeroGrK=doppia;

  // GEOIDE
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=245; goto LEGGIGRID_ERR; }
  _strupr_s(recs,255);
  if (strncmp(recs,"GEOIDE",6)!=0&&strncmp(recs+1,"GEOIDE",6)!=0)  { codErr=246; goto LEGGIGRID_ERR; }

  // data geoide
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=247; goto LEGGIGRID_ERR; }

  // geoide
  for (i=0;i<quantiNodiQuo;i++) { // era i<140
    nriga++;
    if (fgets(recs,250,idflGr)==NULL) { codErr=248; goto LEGGIGRID_ERR; }
    mettiZero(recs);
    if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=249; goto LEGGIGRID_ERR; }
    doppia=atof(recs); // aggiunto per pb. sep. decimale
    aryQuo[i]=doppia;
  } // i

  // lati ROMA40 nodo 127 (idx 126) quote - se tipo punto => nodo 73 (idx 72) - è idxNodoZeroQuo
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=250; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=251; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  fiZeroQuo=doppia;

  // longi ROMA40 nodo 127 (idx 126) quote - se tipo punto => nodo 73 (idx 72) - è idxNodoZeroQuo
  nriga++;
  if (fgets(recs,250,idflGr)==NULL) { codErr=252; goto LEGGIGRID_ERR; }
  mettiZero(recs);
  if (sscanf_s(recs,"%lf",&doppia)!=1) { codErr=253; goto LEGGIGRID_ERR; }
  doppia=atof(recs); // aggiunto per pb. sep. decimale
  laZeroQuo=doppia;

  // chiusura
  if (idflGr!=NULL) { fclose(idflGr); idflGr=NULL; }

  // assegnazione matrice
  if (!posToIndexP(fiZeroPla,laZeroPla,&rigZeroPla,&colZeroPla)) {
    codErr=254;
    goto LEGGIGRID_ERR;
  } // err

  for (i=0;i<numRigPla;i++) { // era i<6
    for (j=0;j<numColPla;j++) { // era j<6
      if (cheTipoK==1) {
        if (matMaskP[rigZeroPla+i][colZeroPla+j]>0) continue; // i gk1 non devono sovrascrivere i gr2
      } // 1

      matFiEd[rigZeroPla+i][colZeroPla+j]=aryFiEd[i*numColPla+j];   // era i*6
      matLaEd[rigZeroPla+i][colZeroPla+j]=aryLaEd[i*numColPla+j];   // era i*6
      matFiWgs[rigZeroPla+i][colZeroPla+j]=aryFiWgs[i*numColPla+j]; // era i*6
      matLaWgs[rigZeroPla+i][colZeroPla+j]=aryLaWgs[i*numColPla+j]; // era i*6
      matMaskP[rigZeroPla+i][colZeroPla+j]=cheTipoK;
    } // j
  } // i

  // assegnazione matrice quote
  if (!posToIndexH(fiZeroQuo,laZeroQuo,&rigZeroQuo,&colZeroQuo)) {
    codErr=255;
    goto LEGGIGRID_ERR;
  } // err

  for (i=0;i<numRigQuo;i++) { // era i<10
    for (j=0;j<numColQuo;j++) { // era j<14
      if (cheTipoK==1) {
        if (matMaskH[rigZeroQuo+i][colZeroQuo+j]>0) continue; // i gk1 non devono sovrascrivere i gr2
      } // 1

      matHeight[rigZeroQuo+i][colZeroQuo+j]=(float) aryQuo[idxNodoZeroQuo-i*numColQuo+j]; // era aryQuo[126-i*14+j]
      matMaskH[rigZeroQuo+i][colZeroQuo+j]=cheTipoK;
    } // j
  } // i

  // assegnazione matrici griglie K
  if (!posToIndexK(fiZeroGrK,laZeroGrK,&rigZeroGrK,&colZeroGrK)) {
    codErr=256;
    goto LEGGIGRID_ERR;
  } // err

  for (i=0;i<numRigPla;i++) { // era i<6
    for (j=0;j<numColPla;j++) { // era i<6
      matFiGrK[rigZeroGrK+i][colZeroGrK+j]=aryFiGrK[i*numColPla+j]; // era i*6
      matLaGrK[rigZeroGrK+i][colZeroGrK+j]=aryLaGrK[i*numColPla+j]; // era i*6
      matQuGrK[rigZeroGrK+i][colZeroGrK+j]=aryQuGrK[i*numColPla+j]; // era i*6
      matMaskK[rigZeroGrK+i][colZeroGrK+j]=cheTipoK;
    } // j
  } // i

  if (g_snRep) mettiRep("Caricato file: ");
  if (g_snRep) mettiRep(strFileGrid); // occhio che è 260
  if (g_snRep) mettiRep("\n");

  return true;

LEGGIGRID_ERR:

  if (idflGr!=NULL) fclose(idflGr);

  sprintf_s(strErr,255,"*** Errore %i nella lettura del file\n",codErr);
  if (g_snRep) mettiRep(strErr);

  recs[255]='\0';
  if (g_snRep) mettiRep(">");
  if (g_snRep) mettiRep(recs);

  return false;

} // leggiFileIgmGrdK

// ***************************************************************************

int ConveRgo_SISTEMI_GRI(double fiIn,double laIn,char sistIn,char sistOut,double *fiOut,double *laOut)
{
  // errore => return 0;
  // sistIn|sistOut: ROMA ED WGS84

  int idxRig,idxCol;
  double prcFi,prcLa,prcFL,dltFi,dltLa;
  double fiWgs,laWgs,fiRoma,laRoma,fiEd50,laEd50;

  if (g_snRep) {
    sprintf_s(strPerRep,255,"* ConveRgo_SISTEMI_GRI in (%.8lf,%.8lf sistIn=%hi sistOut=%hi)\n",fiIn,laIn,sistIn,sistOut);
    mettiRep(strPerRep);
  } // g_snRep

  if (!g_snLettiGrid) return -1;

  *fiOut=0.0;
  *laOut=0.0;

  if (!g_snAuth) {
    if (!ctrlAuth(NULL,true)) return 0;
  } // !g_snAuth

  if (fiIn<g_fiMinAssDec||fiIn>g_fiMaxAssDec) return -2;
  if (laIn<g_laMinAssDec||laIn>g_laMaxAssDec) return -2;

  fiWgs=0.0;
  laWgs=0.0;
  fiEd50=0.0;
  laEd50=0.0;

  // ctrl etrf2000
  if (sistIn==ETRF2000) {
    fiWgs=fiIn;
    laWgs=laIn;

    if (!posToIndexK(fiWgs,laWgs,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    if (matMaskK[idxRig][idxCol]>0&&matMaskK[idxRig+1][idxCol]>0&&
        matMaskK[idxRig][idxCol+1]>0&&matMaskK[idxRig+1][idxCol+1]>0) {

      prcFL=prcFi*prcLa;

      dltFi=(matFiGrK[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiGrK[idxRig][idxCol+1]*(prcLa-prcFL)+
        matFiGrK[idxRig+1][idxCol]*(prcFi-prcFL)+matFiGrK[idxRig+1][idxCol+1]*prcFL);

      dltLa=(matLaGrK[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaGrK[idxRig][idxCol+1]*(prcLa-prcFL)+
        matLaGrK[idxRig+1][idxCol]*(prcFi-prcFL)+matLaGrK[idxRig+1][idxCol+1]*prcFL);

      fiIn=fiWgs-dltFi;
      laIn=laWgs-dltLa;
    } // matMaskK
    else {
      return -3;
    } // !matMaskK

    sistIn=ETRF89;
  } // ETRF2000

  if (sistIn==ROMA) {
    fiRoma=fiIn;
    laRoma=laIn;

    //mettiRep("ROMA40 %.8lf %.8lf\n",fiRoma,laRoma);

    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    //mettiRep("Indici %i %i\n",idxRig,idxCol);

    prcFL=prcFi*prcLa;

    dltFi=(matFiWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matFiWgs[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matLaWgs[idxRig+1][idxCol+1]*prcFL);

    fiWgs=fiRoma+dltFi;
    laWgs=laRoma+dltLa;

    // ED50
    fiEd50=fiRoma+(matFiEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiEd[idxRig+1][idxCol]*(prcFi-prcFL)+matFiEd[idxRig+1][idxCol+1]*prcFL);

    laEd50=laRoma+(matLaEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaEd[idxRig+1][idxCol]*(prcFi-prcFL)+matLaEd[idxRig+1][idxCol+1]*prcFL);
  } // ROMA
  else if (sistIn==ETRF89||sistIn==ETRF2000) {
    fiWgs=fiIn;
    laWgs=laIn;

    //mettiRep("WGS84 %.8lf %.8lf\n",fiWgs,laWgs);

    fiRoma=fiIn; // per innnescare,
    laRoma=laIn; // poi migliora iterando

    // iter 1
    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    //mettiRep("Indici_1 %i %i\n",idxRig,idxCol);

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    prcFL=prcFi*prcLa;

    dltFi=(matFiWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matFiWgs[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matLaWgs[idxRig+1][idxCol+1]*prcFL);

    fiRoma=fiWgs-dltFi;
    laRoma=laWgs-dltLa;

    // iter 2
    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    //mettiRep("Indici_2 %i %i\n",idxRig,idxCol);

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    prcFL=prcFi*prcLa;

    dltFi=(matFiWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matFiWgs[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matLaWgs[idxRig+1][idxCol+1]*prcFL);

    fiRoma=fiWgs-dltFi;
    laRoma=laWgs-dltLa;

    // iter 3
    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    //mettiRep("Indici_3 %i %i\n",idxRig,idxCol);

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    prcFL=prcFi*prcLa;

    dltFi=(matFiWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matFiWgs[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matLaWgs[idxRig+1][idxCol+1]*prcFL);

    fiRoma=fiWgs-dltFi;
    laRoma=laWgs-dltLa;

    // ED50
    fiEd50=fiRoma+(matFiEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiEd[idxRig+1][idxCol]*(prcFi-prcFL)+matFiEd[idxRig+1][idxCol+1]*prcFL);

    laEd50=laRoma+(matLaEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaEd[idxRig+1][idxCol]*(prcFi-prcFL)+matLaEd[idxRig+1][idxCol+1]*prcFL);
  } // WGS
  else if (sistIn==ED) {
    fiEd50=fiIn;
    laEd50=laIn;

    //mettiRep("ED50 %.8lf %.8lf\n",fiEd50,laEd50);

    fiRoma=fiIn; // per innnescare,
    laRoma=laIn; // poi migliora iterando

    // iter 1
    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    //mettiRep("Indici_1 %i %i\n",idxRig,idxCol);

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    prcFL=prcFi*prcLa;

    dltFi=(matFiEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiEd[idxRig+1][idxCol]*(prcFi-prcFL)+matFiEd[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaEd[idxRig+1][idxCol]*(prcFi-prcFL)+matLaEd[idxRig+1][idxCol+1]*prcFL);

    fiRoma=fiEd50-dltFi;
    laRoma=laEd50-dltLa;

    // iter 2
    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    //mettiRep("Indici_2 %i %i\n",idxRig,idxCol);

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    prcFL=prcFi*prcLa;

    dltFi=(matFiEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiEd[idxRig+1][idxCol]*(prcFi-prcFL)+matFiEd[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaEd[idxRig+1][idxCol]*(prcFi-prcFL)+matLaEd[idxRig+1][idxCol+1]*prcFL);

    fiRoma=fiEd50-dltFi;
    laRoma=laEd50-dltLa;

    // iter 3
    if (!posToIndexP(fiRoma,laRoma,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    //mettiRep("Indici_3 %i %i\n",idxRig,idxCol);

    if (matMaskP[idxRig][idxCol]==0||matMaskP[idxRig+1][idxCol]==0||
      matMaskP[idxRig][idxCol+1]==0||matMaskP[idxRig+1][idxCol+1]==0) return -3;

    prcFL=prcFi*prcLa;

    dltFi=(matFiEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiEd[idxRig+1][idxCol]*(prcFi-prcFL)+matFiEd[idxRig+1][idxCol+1]*prcFL);

    dltLa=(matLaEd[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaEd[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaEd[idxRig+1][idxCol]*(prcFi-prcFL)+matLaEd[idxRig+1][idxCol+1]*prcFL);

    fiRoma=fiEd50-dltFi;
    laRoma=laEd50-dltLa;

    // ED50
    fiWgs=fiRoma+(matFiWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matFiWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matFiWgs[idxRig+1][idxCol+1]*prcFL);

    laWgs=laRoma+(matLaWgs[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaWgs[idxRig][idxCol+1]*(prcLa-prcFL)+
      matLaWgs[idxRig+1][idxCol]*(prcFi-prcFL)+matLaWgs[idxRig+1][idxCol+1]*prcFL);
  } // ED
  else return 0;

  if (sistOut==ROMA) {
    *fiOut=fiRoma;
    *laOut=laRoma;
  } // ROMA
  else if (sistOut==ETRF89) {
    *fiOut=fiWgs;
    *laOut=laWgs;
  } // WGS
  else if (sistOut==ED) {
    *fiOut=fiEd50;
    *laOut=laEd50;
  } // ED
  else if (sistOut==ETRF2000) {
    if (!posToIndexK(fiWgs,laWgs,&idxRig,&idxCol,&prcFi,&prcLa)) return -2;

    if (matMaskK[idxRig][idxCol]>0&&matMaskK[idxRig+1][idxCol]>0&&
        matMaskK[idxRig][idxCol+1]>0&&matMaskK[idxRig+1][idxCol+1]>0) {

      prcFL=prcFi*prcLa;

      dltFi=(matFiGrK[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matFiGrK[idxRig][idxCol+1]*(prcLa-prcFL)+
        matFiGrK[idxRig+1][idxCol]*(prcFi-prcFL)+matFiGrK[idxRig+1][idxCol+1]*prcFL);

      dltLa=(matLaGrK[idxRig][idxCol]*(1.0-prcFi-prcLa+prcFL)+matLaGrK[idxRig][idxCol+1]*(prcLa-prcFL)+
        matLaGrK[idxRig+1][idxCol]*(prcFi-prcFL)+matLaGrK[idxRig+1][idxCol+1]*prcFL);

      *fiOut=fiWgs+dltFi;
      *laOut=laWgs+dltLa;
    } // matMaskK
    else {
      *fiOut=fiWgs;
      *laOut=laWgs;
      return -3;
    } // !matMaskK

  } // ETRF2000
  else return 0;

  if (g_snRep) mettiRep("* ConveRgo_SISTEMI_GRI ok\n");

  return 1;
} // ConveRgo_SISTEMI_GRI

// ******************************************************************************

int ConveRgo_GET_VERS(char *mainVer,char *subVer)
{
  if (g_snRep) mettiRep("* ConveRgo_GET_VERS\n");

  *mainVer=VERS_1;
  *subVer=VERS_2;

  return 1;
} // ConveRgo_GET_VERS

// ******************************************************************************

BOOL mettiRep(char strRep[260])  //REP
{
  // return;

  g_idflRep=NULL;
  g_idflRep=fopen(g_nomeFileRep,"a");
  if (g_idflRep==NULL) return false;

  fprintf(g_idflRep,"%s",strRep);

  fclose(g_idflRep);

  return true;
} // mettiRep

// ******************************************************************************

void mettiRep(char strFrm[260],double valUno,double valDue)
{
  char strRep[260];

  sprintf_s(strRep,255,strFrm,valUno,valDue);

  mettiRep(strRep);

} // mettiRep

// ******************************************************************************

void mettiRep(char strFrm[260],int valUno,int valDue)
{
  char strRep[260];

  sprintf_s(strRep,255,strFrm,valUno,valDue);

  mettiRep(strRep);

} // mettiRep

// ******************************************************************************

void mettiZero(char recs[260])
{
  int i,k;

  if (recs[0]=='-'&&recs[1]=='.'&&isdigit(recs[2])) {
    k=(int) strlen(recs);

    for (i=k;i>1;i--) {
      recs[i]=recs[i-1];
    } // i

    recs[1]='0';

    return;
  } // -.123

  if (recs[0]=='.'&&isdigit(recs[1])) {
    k=(int) strlen(recs);

    for (i=k;i>0;i--) {
      recs[i]=recs[i-1];
    } // i

    recs[0]='0';

    return;
  } // .123

  return;
} // mettiZero

// ******************************************************************************

int ConveRgo_REP_ON(void)
{
  if (g_snRep) mettiRep("* ConveRgo_REP_ON\n");

  if (g_snRep) return 0;

  if (!mettiRep("\nReport attivo\n")) return false;

  g_snRep=true;

  return 1;
} // ConveRgo_REP_ON

// ******************************************************************************

int ConveRgo_REP_OFF(void)
{
  if (g_snRep) mettiRep("* ConveRgo_REP_OFF\n");

  g_snRep=false;

  return 1;
} // ConveRgo_REP_OFF

// ******************************************************************************

BOOL leggiFileGra(char strFileGra[MY_MAX_PATH+3])
{
  FILE *idflGr;
  char recs[260],strErr[260];
  int i,j,k,numRig,numCol,numPunti,codErr;
  long int nriga;
  double xMin,yMin,xMax,yMax,passoX,passoY,xWgs,yWgs,x1Wgs,y1Wgs,dX,dY,valEnne;
  BOOL snPrimo;

  codErr=11;

  idflGr=NULL;
  idflGr=fopen(strFileGra,"r");
  if (idflGr==NULL) {
    if (g_snRep) mettiRep("*** Errore nell'apertura del file\n");
    return false;
  } // !Open

  xMin=0.0; yMin=0.0; xMax=0.0; yMax=0.0; passoX=0.0; passoY=0.0;

  snPrimo=true;
  numPunti=0;

  nriga=0l;
  while (fgets(recs,250,idflGr)!=NULL) {
    nriga++;

    k=(int) strlen(recs);

    if (k<3) {
      continue;
    } // vuota

    for (i=0;i<k;i++) {
      if (recs[i]<' ') recs[i]=' ';
    } // i

    if (sscanf_s(recs,"%lf %lf %lf",&xWgs,&yWgs,&valEnne)!=3) { codErr=12; goto LEGGIGRA_ERR; }

    xWgs=ConveRgo_GESI_DECI(xWgs);
    yWgs=ConveRgo_GESI_DECI(yWgs);

    if (snPrimo) {
      xMin=xWgs;
      yMin=yWgs;
      xMax=xWgs;
      yMax=yWgs;
      passoX=1.0E+06;
      passoY=1.0E+06;

      x1Wgs=xWgs;
      y1Wgs=yWgs;

      snPrimo=false;
      numPunti=1;
    } // snPrimo
    else {
      if (xWgs<xMin) xMin=xWgs;
      if (yWgs<yMin) yMin=yWgs;
      if (xWgs>xMax) xMax=xWgs;
      if (yWgs>yMax) yMax=yWgs;
      dX=fabs(xWgs-x1Wgs);
      dY=fabs(yWgs-y1Wgs);
      if (dX>0.0001&&dX<passoX) passoX=dX;
      if (dY>0.0001&&dY<passoY) passoY=dY;

      numPunti++;
    } // !snPrimo

  } // while

  if (numPunti<4) {
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore: numero dei punti xyz insufficiente (%i)\n",numPunti);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  if (passoX>10.0||passoY>10.0) {
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore nel calcolo del passo della maglia (%.4lfx%.4lf)\n",passoX,passoY);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  rewind(idflGr);

  numCol=(long int) ((xMax-xMin)/passoX+0.5)+1;
  numRig=(long int) ((yMax-yMin)/passoY+0.5)+1;

  if (numCol<1||numCol>NMAX_ELEGRA) {
    //g_sMsgErr.Format("Errore nel numero delle colonne (%i) per il file \n\n",numCol);
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore nel numero di colonne della matrice (%i)\n",numCol);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  if (numRig<1||numRig>NMAX_ELEGRA) {
    //g_sMsgErr.Format("Errore nel numero delle righe (%i) per il file \n\n",numRig);
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore nel numero di righe della matrice (%i)\n",numRig);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  if (!g_aryMaglieGra[g_numAryGra].allocaRigCol(numRig,numCol)) {
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore nell'allocazione di memoria (%ix%i celle)\n",numRig,numCol);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  if (!g_aryMaglieGra[g_numAryGra].m_snAllocato) {
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore: memoria non allocata correttamente (gra %i)\n",g_numAryGra);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  g_aryMaglieGra[g_numAryGra].m_maxFiDeci=yMax;
  g_aryMaglieGra[g_numAryGra].m_oriLaDeci=xMin;
  g_aryMaglieGra[g_numAryGra].m_pasFiDeci=passoY;
  g_aryMaglieGra[g_numAryGra].m_pasLaDeci=passoX;

  for (i=0;i<numRig;i++) {
    for (j=0;j<numCol;j++) {
      g_aryMaglieGra[g_numAryGra].m_pMatEnne[i*numCol+j]=0.0;
      g_aryMaglieGra[g_numAryGra].m_pMatMask[i*numCol+j]=0;
    } // j
  } // i

  nriga=0l;
  while (fgets(recs,250,idflGr)!=NULL) {
    nriga++;

    k=(int) strlen(recs);

    if (k<3) {
      continue;
    } // vuota

    for (i=0;i<k;i++) {
      if (recs[i]<' ') recs[i]=' ';
    } // i

    if (sscanf_s(recs,"%lf %lf %lf",&xWgs,&yWgs,&valEnne)!=3) { codErr=13; goto LEGGIGRA_ERR; }

    xWgs=ConveRgo_GESI_DECI(xWgs);
    yWgs=ConveRgo_GESI_DECI(yWgs);

    i=(int) ((yMax-yWgs)/passoY+0.5);  // riga 0 a nord
    j=(int) ((xWgs-xMin)/passoX+0.5);

    if (i<0||i>=numRig||j<0||j>=numCol) { codErr=14; goto LEGGIGRA_ERR; }

    g_aryMaglieGra[g_numAryGra].m_pMatEnne[i*numCol+j]=valEnne;
    g_aryMaglieGra[g_numAryGra].m_pMatMask[i*numCol+j]=1;

  } // while

  fclose(idflGr);
  idflGr=NULL;

  g_aryMaglieGra[g_numAryGra].m_snAssegnato=true;

  g_numAryGra++;
  if (g_numAryGra>NMAX_ARYGRA) {
    fclose(idflGr);
    sprintf_s(strErr,255,"*** Errore: troppi grigliati adattati (max )\n",NMAX_ARYGRA);
    if (g_snRep) mettiRep(strErr);
    return false;
  } // err

  if (g_snRep) mettiRep("Caricato file: ");
  if (g_snRep) mettiRep(strFileGra); // occhio che è 260
  if (g_snRep) mettiRep("\n");

  return true;

LEGGIGRA_ERR:

  if (idflGr!=NULL) fclose(idflGr);

  sprintf_s(strErr,255,"*** Errore %i alla riga %li del file\n",codErr,nriga);
  if (g_snRep) mettiRep(strErr);

  recs[255]='\0';
  if (g_snRep) mettiRep(">");
  if (g_snRep) mettiRep(recs);

  return false;

} // leggiFileGra

// ******************************************************************************

BOOL posToIndexGra(double fi,double la,int *idxGra,int *idxAS,int *idxAD,int *idxBS,int *idxBD,double *prcFi,double *prcLa)
{
  int i,cheGra,rig,col;
  double quantiFi,quantiLa;

  *idxGra=0;
  *idxAS=0;
  *idxAD=0;
  *idxBS=0;
  *idxBD=0;
  *prcFi=0.0;
  *prcLa=0.0;

  cheGra=-1;

  for (i=0;i<g_numAryGra;i++) {
    quantiFi=(g_aryMaglieGra[i].m_maxFiDeci-fi)/g_aryMaglieGra[i].m_pasFiDeci;
    quantiLa=(la-g_aryMaglieGra[i].m_oriLaDeci)/g_aryMaglieGra[i].m_pasLaDeci;

    rig=(int) quantiFi;
    col=(int) quantiLa;

    if (rig<0||rig>g_aryMaglieGra[i].m_numRig-2) continue;  // dovrebbe essere
    if (col<0||col>g_aryMaglieGra[i].m_numCol-2) continue;  // pieno anche -1, ma...

    *idxAS = rig*g_aryMaglieGra[i].m_numCol+col;
    *idxAD = *idxAS+1;
    *idxBS = (rig+1)*g_aryMaglieGra[i].m_numCol+col;
    *idxBD = *idxBS+1;

    if (g_aryMaglieGra[i].m_pMatMask[*idxAS]==0) continue;
    if (g_aryMaglieGra[i].m_pMatMask[*idxAD]==0) continue;
    if (g_aryMaglieGra[i].m_pMatMask[*idxBS]==0) continue;
    if (g_aryMaglieGra[i].m_pMatMask[*idxBD]==0) continue;

    *idxGra=i;

    *prcFi = quantiFi-rig;
    *prcLa = quantiLa-col;

    return true;
  } // i

  return false;

} // posToIndexGra

// ******************************************************************************

BOOL leggiPathGridDaReg(char pathGrid[520])
{
  char *perReg;
  char strRep[260];
  unsigned char strValReg[520];
  unsigned long int bho,il,lunValReg;
  long int risp;
  BOOL snTro;
  HKEY chiaveReg;

  pathGrid[0]='\0';

  perReg=new char [520];

  strcpy_s(perReg,512,"Software\\ConveRgo\\Settings");

  bho=1;

  snTro=false;

  risp=RegOpenKeyEx(HKEY_CURRENT_USER,perReg,0,KEY_QUERY_VALUE,&chiaveReg); //KEY_ALL_ACCESS
  if (risp==ERROR_SUCCESS) {
    lunValReg=512;
    risp=RegQueryValueEx(chiaveReg,"GRID_FOLDER",NULL,&bho,strValReg,&lunValReg);
    RegCloseKey(chiaveReg);

    if (risp==ERROR_SUCCESS) {
      for (il=0l;il<lunValReg;il++) {
        pathGrid[il]=strValReg[il];
      } // il
      pathGrid[lunValReg]='\0';
    } // !err
    else {
      sprintf_s(strRep,255,"Errore %i nella funzione RegQueryValueEx\n",risp);
      if (g_snRep) mettiRep(strRep);
    } // err
  } // !err
  else {
    sprintf_s(strRep,255,"Errore %i nella funzione RegOpenKeyEx\n",risp);
    if (g_snRep) mettiRep(strRep);
  } // err

  delete [] perReg;

  if (strlen(pathGrid)>0) return true;
  else return false;

} // leggiPathGridDaReg

// ******************************************************************************

BOOL scriviPathGridSuReg(char pathGrid[520])
{
  CString sPerReg,sNomeChiave,sValChiave;
  CRegKey regKey;
  LONG lErr;
  DWORD dwDisposition;

  sNomeChiave="GRID_FOLDER";
  sValChiave=pathGrid;

  dwDisposition=NULL;

  if (sNomeChiave.IsEmpty()) return false;

  sPerReg="Software\\ConveRgo\\Settings";

  lErr=regKey.Create(HKEY_CURRENT_USER,sPerReg,REG_NONE,REG_OPTION_NON_VOLATILE,KEY_ALL_ACCESS,NULL,&dwDisposition);
  if (lErr!=ERROR_SUCCESS) {
    return false;
  } // err

  lErr=regKey.SetStringValue(sNomeChiave,sValChiave,REG_SZ);

  regKey.Close();

  if (lErr!=ERROR_SUCCESS) {
    return false;
  } // err

  return true;

} // scriviPathGridSuReg

// ******************************************************************************

BOOL leggiSnRepDaReg(void)
{
  char *perReg;
  unsigned char strValReg[520];
  unsigned long int bho,lunValReg;
  long int risp;
  BOOL snRep;
  HKEY chiaveReg;

  perReg=new char [520];

  strcpy_s(perReg,512,"Software\\ConveRgo\\ReportOn");
  bho=1;

  snRep=false;

  risp=RegOpenKeyEx(HKEY_LOCAL_MACHINE,perReg,0,KEY_QUERY_VALUE,&chiaveReg); //KEY_ALL_ACCESS
  if (risp==ERROR_SUCCESS) {
    lunValReg=512;
    risp=RegQueryValueEx(chiaveReg,"",NULL,&bho,strValReg,&lunValReg);
    RegCloseKey(chiaveReg);

    if (risp==ERROR_SUCCESS) {
      if (strValReg[0]=='T'||strValReg[0]=='Y'||strValReg[0]=='S') snRep=true;
      if (strValReg[0]=='t'||strValReg[0]=='y'||strValReg[0]=='s') snRep=true;
    } // !err
  } // !err

  delete [] perReg;

  return snRep;

} // leggiSnRepDaReg

// ******************************************************************************

int ConveRgo_GET_LIM_GEO(double *fiMin,double *laMin,double *fiMax,double *laMax)
{
  if (g_snRep) mettiRep("* ConveRgo_GET_LIM_GEO\n");

  *fiMin=g_fiMinAssDec;
  *fiMax=g_fiMaxAssDec;
  *laMin=g_laMinAssDec;
  *laMax=g_laMaxAssDec;

  return 1;
} // ConveRgo_GET_LIM_GEO

// ******************************************************************************
// ******************************************************************************
