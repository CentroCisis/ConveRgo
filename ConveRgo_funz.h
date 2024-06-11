/*SPDX-License-Identifier: AGPL-3.0-only*/
Copyright (C) 2024 CISIS - Centro Interregionale per I Sistemi Informatici Geografici e Statistici
This program is free software: you can redistribute it and/or modify it under the terms of the GNU Affero General Public License as published by the Free Software Foundation, version 3.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public License for more details.
You should have received a copy of the GNU Affero General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.

// ConveRgo_funzioni.h

#include "ConveRgo_define.h"

int ConveRgo_INIZIALIZZA(void);

int ConveRgo_SET_PATH_GRI(char grPath[520]);
int ConveRgo_LOAD_GRI(void);
int ConveRgo_TEST_GRI(char sistema,double laMin,double fiMin,double laMax,double fiMax);
int ConveRgo_QUANTI_GRI(int *quantiGr1,int *quantiGr2,int *quantiGrA);
int ConveRgo_QUANTI_GK(int *quantiGk1,int *quantiGk2);

int ConveRgo_SISTEMI_GRI(double fiIn,double laIn,char sistIn,char sistOut,double *fiOut,double *laOut);

int ConveRgo_HELL_QSLM(double fiWgs,double laWgs,double hWgs,double *quoOrto);
int ConveRgo_QSLM_HELL(double fiWgs,double laWgs,double quoOrto,char sistIn,double *hWgs);
int ConveRgo_HELL_QSLM_G(double fiWgs,double laWgs,double hWgs,char cheGeoide,double *quoOrto);
int ConveRgo_QSLM_HELL_G(double fiWgs,double laWgs,double quoOrto,char sistIn,char cheGeoide,double *hWgs);

int ConveRgo_PIA_GEO(double N_In,double E_In,char fusoIn,char sistema,double *fiOut,double *laOut);
int ConveRgo_GEO_PIA(double fiIn,double laIn,char sistema,char fusoOut,double *N_Out,double *E_Out);

int ConveRgo_PIA_GEO_U(double N_In,double E_In,double *fiOut,double *laOut);
int ConveRgo_GEO_PIA_U(double fiIn,double laIn,double *N_Out,double *E_Out);

double ConveRgo_GESI_DECI(double gesiIn);
double ConveRgo_DECI_GESI(double deciIn,int numCifreSec);

double ConveRgo_GPD_DECI(double gpdIn);
double ConveRgo_DECI_GPD(double deciIn,int numCifrePri);

int ConveRgo_GET_VERS(char *mainVer,char *subVer);

int ConveRgo_REP_ON(void);
int ConveRgo_REP_OFF(void);

// -----------------------------------------

BOOL leggiPathGridDaReg(char pathGrid[520]);
BOOL scriviPathGridSuReg(char pathGrid[520]);
BOOL leggiSnRepDaReg(void);

BOOL posToIndexP(double fi,double la,int *rig,int *col);
BOOL posToIndexP(double fi,double la,int *rig,int *col,double *prcFi,double *prcLa);
BOOL posToIndexH(double fi,double la,int *rig,int *col);
BOOL posToIndexH(double fi,double la,int *rig,int *col,double *prcFi,double *prcLa);
BOOL posToIndexK(double fi,double la,int *rig,int *col);
BOOL posToIndexK(double fi,double la,int *rig,int *col,double *prcFi,double *prcLa);

BOOL leggiFileIgmGr12(char strFileGrid[MY_MAX_PATH+3],char strSoloNome[43],int cheTipoR); // CString npFileGrid
BOOL leggiFileIgmGrdK(char strFileGrid[MY_MAX_PATH+3],char strSoloNome[43],int cheTipoK);

BOOL leggiFileGra(char strFileGra[MY_MAX_PATH+3]); // CString npFileGra
BOOL posToIndexGra(double fi,double la,int *idxGra,int *idxAS,int *idxAD,int *idxBS,int *idxBD,double *prcFi,double *prcLa);

BOOL mettiRep(char strRep[260]); //REP
void mettiRep(char strFrm[260],double valUno,double valDue); //REP
void mettiRep(char strFrm[260],int valUno,int valDue); //REP

void mettiZero(char recs[260]);
