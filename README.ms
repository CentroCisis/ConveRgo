# ConveRgo - Conversioni di coordinate per le Regioni

## Descrizione
Il programma ConveRgo serve ad eseguire trasformazioni di coordinate fra i vari sistemi di 
riferimento in cui sono espressi i dati geografici delle Amministrazioni regionali (ROMA40, 
ED50, ETRS89 nelle due realizzazioni ETRF89 e ETRF2000), considerando anche i rispettivi 
sistemi cartografici (Gauss-Boaga, UTM-ED50, UTM-ETRF89 e UTM-ETRF2000). Viene 
considerata anche la componente altimetrica, per le conversioni fra quote ellissoidiche e 
geoidiche.

## Conversioni di coordinate
Il programma permette il passaggio di coordinate fra i sistemi ETRS89 (nelle realizzazioni 
ETRF89 e ETRF2000), ED50 e ROMA40, con le relative rappresentazioni UTM e Gauss-Boaga, 
elaborando singoli punti oppure file con liste di coordinate, file di georeferenziazione (es. tfw), 
shapefile, dxf (di  solo contenuto cartografico) e altri formati.
I calcoli possono essere eseguiti sulla base dei "grigliati" nazionali: il programma richiede 
l'indicazione della cartella che contiene i file *.GRn o *.GKn (con n = 1 oppure 2), quindi carica 
automaticamente tutti quelli presenti. Viene considerata anche la componente altimetrica, 
con le opportune trasformazioni fra quote ellissoidiche e geoidiche.
Il programma non contiene al proprio interno alcun valore dei grigliati: il reperimento dei 
relativi file compete all’utente. Nel caso in cui non siano presenti i grigliati nell’area relativa ai 
file da trasformare, ConveRgo permette di eseguire il calcolo col modello approssimato 
(modalità “CartLab2”).
E’ possibile trasformare le coordinate Cassini-Soldner dei sistemi catastali (Bessel su 
Genova) per le zone di cui si conosca il centro di emanazione. Il programma contiene i centri 
delle zone più estese, che rappresentano circa il 50% del territorio nazionale.
L’impostazione generale del programma è tale da consentire tutte le possibili trasformazioni 
fra i sistemi considerati. La finestra di dialogo principale è quindi stata progettata in modo 
simmetrico: la parte di sinistra è riservata all’input e quella di destra all’output. Le due parti 
appaiono uguali poiché il programma consente, tranne qualche eccezione, tutte le possibili 
combinazioni di scelta.

## Codici sorgenti pubblicati, per le funzioni di calcolo e di gestione dei "grigliati"
Le funzioni presenti nei codici sorgenti costituiscono la parte principale del software, consentendo di eseguire 
le conversioni di coordinate fra i vari sistemi di riferimento nazionali secondo le modalità
ufficiali prescritte dall'allegato tecnico del D.M. del novembre 2011 "Adozione del Sistema
di riferimento geodetico nazionale".
I codici sorgenti sono espressi nel linguaggio C++ e corredati da commenti e note
esplicative di supporto all'eventuale riutilizzo da parte di altri soggetti.
Tramite i sorgenti, gli sviluppatori possono implementare le conversioni di coordinate 
semplicemente integrando le funzioni di calcolo e di gestione dei "grigliati" all'interno del proprio software.
Si precisa che i sorgenti contengono le funzioni di lettura dei file dei grigliati ma non contengono i valori numerici dei modelli. Tali valori,
strutturati per un utilizzo razionale all'interno di un software orientato alle conversioni di
coordinate a livello nazionale (quindi riuniti in un'unica struttura anziché suddivisi in
molteplici file) sono a disposizione del CISIS, ma non sono pubblicabili: uno specifico
accordo fra CISIS e IGM prescrive infatti che possano circolare solo all'interno delle
Regioni o degli Enti subordinati.


## Informazioni generali su ConveRgo
Le funzioni permettono la trasformazione di coordinate fra i vari sistemi (ETRF2000, ETRF89, ED50 e ROMA40) mediante l’utilizzo dei “grigliati” IGM (sia quelli “nazionali” *.gr1-2 e *.gk1-2 sia quelli frutto di adattamento locale del geoide, denominati *.gra).
Per l’utilizzo dei grigliati, i file *.gr1, *.gr2, *.gk1, *.gk2 e *.gra devono essere preventivamente caricati mediante chiamata dell’apposita funzione  ConveRgo_LOAD_GRI.
I file devono risiedere all’interno di una stessa cartella, il cui percorso viene letto nella chiave di registro:
HKEY_LOCAL_MACHINE\Software\ConveRgo\GridFilePath
oppure può essere assegnato direttamente mediante la funzione   ConveRgo_SET_PATH_GRI.
I file *.gra sono derivati dagli adattamenti locali del geoide realizzati da IGM su richiesta specifica. Essi devono essere espressi in formato ASCII, contenenti una riga di intestazione e un numero variabile di record di uguale formato. Ciascun record è formato da 5 valori separati da spazio, contenenti nell’ordine un codice di punto, la latitudine e la longitudine in gradi sessagesimali nella forma gg.ppssdddd, i valori di scostamento fra geoide ed ellissoide (N) nazionale e locale in metri. Tali punti formano una griglia regolare a maglia quadrata, il cui passo ed estensione si ricavano dalla posizione dei punti stessi.
E’ presente una funzione di “test” della presenza dei grigliati (ConveRgo_TEST_GRI), che verifica un’area rettangolare a partire dalle coordinate di due spigoli.
Le funzioni permettono anche le trasformazioni di coordinate fra geografiche e piane all’interno dello stesso sistema, considerando le formule di Gauss, più altre funzioni di utilità, come ad esempio le conversioni fra differenti unità di misura degli angoli.
Negli argomenti delle funzioni di trasformazione, le coordinate geografiche sono attese in gradi sessadecimali, con longitudine da Greenwich anche per il sistema ROMA40.
Al solo scopo di debug sono presenti le funzioni ConveRgo_REP_ON e ConveRgo_REP_OFF, che attivano/disattivano la registrazione delle principali attività delle funzioni in un report di testo accodato al file C:\ConveRgo_debug.txt
La variabile “sistema” è così codificata:
1 = ROMA40 2 = ED50
3  =  ETRF89 6 = ETRF2000
La variabile “fuso” è così codificata:
32 = fuso 32 oppure Ovest
33 = fuso 33 oppure Est
34 = fuso 34
0  = scelta automatica (in base alla longitudine)
Le funzioni di tipo “int” restituiscono un valore che codifica l’esito dell’azione richiesta, come descritto nel seguito. In generale, le funzioni restituiscono 1 in caso di esecuzione regolare e 0 oppure valori negativi in caso di errore.
Segue la descrizione delle funzioni disponibili.
 
## Percorso file con grigliati IGM:
int ConveRgo_SET_PATH_GRI(char grPath[520])
-	restituisce 0 in caso di errore, altrimenti il numero di file dei grigliati presenti in quella cartella

Caricamento grigliati IGM:
int ConveRgo_LOAD_GRI()
-	restituisce:   1 = ok,   0 = errore generico,   -1 = percorso non impostato, -2 = non trovati file dei grigliati

Test grigliati IGM su box:
int ConveRgo_TEST_GRI(char sistema,double laMin,double fiMin,double laMax,double fiMax)
-	restituisce: 0 = errore generico, 1 = presenza dei grigliati su tutta l’area, -1 = assenza parziale,  -2 = assenza  totale

Passaggio fra i sistemi ROMA40, ED50, ETRF89 ed ETRF2000:
int ConveRgo_SISTEMI_GRI(double fiIn,double laIn,char sistIn,char sistOut,double *fiOut,double *laOut)
-	restituisce: 1 = ok, 0 = errore generico, -1 = grigliati non caricati, -2 = punto fuori limiti, -3 = assenza dei grigliati nel punto

Trasformazioni altimetriche:
int ConveRgo_HELL_QSLM_G(double fiIn,double laIn,double H_In,char cheGeoide,double *quoOut)
int ConveRgo_QSLM_HELL_G(double fiIn,double laIn,double quoIn,char sistIn,char cheGeoide,double *H_Out)
-	cheGeoide = scelta del modello di geoide da utilizzare: 0 = il migliore disponibile nel punto, 2 = grigliato nazionale *.gr1-2 o *.gk1-2,  3 = grigliato “adattato” *.gra
-	restituiscono:  3 = ok con grigliato “adattato” *.gra,   2 = ok con grigliato nazionale, 0 = errore generico,  -1 = assenza di grigliati;  -2 = punto fuori  limiti
Sono disponibili anche le funzioni prive del parametro per la scelta del geoide:
int ConveRgo_HELL_QSLM(double fiIn,double laIn,double H_In,double *quoOut)
int ConveRgo_QSLM_HELL(double fiIn,double laIn,double quoIn,char sistIn,double *H_Out)

Da geografiche a piane:
int ConveRgo_GEO_PIA(double fiIn,double laIn,char sistema,char fusoOut,double *N_Out,double *E_Out) int ConveRgo_GEO_PIA_U(double fiIn,double laIn,double *N_Out,double *E_Out)
-	restituiscono: 1 = ok, 0 = errore generico, -1 = errore negli argomenti, -2 = punto esterno ai limiti geografici

Da piane a geografiche:
int ConveRgo_PIA_GEO(double N_In,double E_In,char fusoIn,char sistema,double *fiOut,double *laOut) int ConveRgo_PIA_GEO_U(double N_In,double E_In,double *fiOut,double *laOut)
-	restituiscono: 1 = ok, 0 = errore generico, -1 = errore negli argomenti

Da sessagesimali a sessadecimali:
double ConveRgo_GESI_DECI(double gesiIn)

Da sessadecimali a sessagesimali:
double ConveRgo_DECI_GESI(double deciIn,int numCifreSec)
-	numCifreSec = cifre decimali per i secondi (per evitare che un successivo arrotondamento provochi 60 secondi)

Da gradi primi e decimali di primo a sessadecimali:
double ConveRgo_GPD_DECI(double gpdIn)

Da sessadecimali a gradi primi e decimali di primo:
double ConveRgo_DECI_GPD(double deciIn,int numCifrePri)
-	numCifrePri = cifre decimali per i primi (per evitare che un successivo arrotondamento provochi 60 primi)
Versione delle funzioni
int ConveRgo_GET_VERS(char *mainVer,char *subVer)
-	mainVer,subVer = numero principale e sottonumero della versione (es. versione 3.2 mainVer=3 subVer=2)
 
## Note sulle modalità di conversione di coordinate
Le “griglie” contenute nei file *.GR1 e *.GR2 utilizzati dalle funzioni per i passaggi fra sistemi contengono le differenze di latitudine e longitudine da sommare alle coordinate espresse in un sistema per ottenere l’altro (cfr. la pubblicazione IGM “La trasformazione tra i sistemi di riferimento utilizzati in Italia” di Donatelli-Maseroli-Pierozzi, sul Bollettino di Geodesia e Scienze Affini n. 4 del 2002 e la tesi di dottorato “Dallo statico al network RTK: l’evoluzione del rilievo” di Ronci, Università di Bologna 2007).
Le griglie per i passaggi fra ROMA40 e ED50 e fra ROMA40 e ETRF89 sono georeferenziate nel sistema ROMA40, con longitudini riferite a Greenwich, e sono definite in coordinate geografiche espresse in gradi sessadecimali.
I valori sono assegnati sui nodi di una griglia regolare, di 153 righe per 109 colonne. Il passo è di 0.08333333° (5’) in latitudine e di 0.125° (7’30”) in longitudine, in modo che la maglia risulti quasi quadrata sul terreno. Le righe sono contate a partire dal basso, cioè da sud a nord. Le coordinate del primo nodo in basso a sinistra sono le seguenti.
Origine delle griglie planimetriche: latitudine = 35°00'00.00", longitudine = 5°57'08.40" (Greenwich) Il passaggio fra ED50 e ETRF89 si realizza mediante un doppio calcolo, passando da ROMA40.
Per la componente altimetrica, cioè per il passaggio da altezza ellissoidica a quota geoidica (s.l.m.m.) e viceversa, esiste un’apposita griglia che contiene i valori di separazione, cioè i valori del modello di geoide (ITALGEO99 nei file
*.gr1 e *.gk1, ITALGEO2005 nei file *.gr2 e *.gk2).
La griglia altimetrica è più fitta di quelle planimetriche: il passo è di 0.03333333° (2’), uguale in latitudine e in longitudine. La griglia è georeferenziata nel sistema ETRF89. Le coordinate del primo nodo in basso a sinistra sono le seguenti:
Origine della griglia altimetrica: latitudine = 35°20'00.00", longitudine = 6°00'0.00"
Con la versione “K” dei grigliati si sono aggiunte le griglie per la gestione del sistema ETRF2000. I file *.gk1 e *.gk2 contengono, oltre ai valori già presenti nei corrispondenti file *.gr1 o *.gr2, anche la differenze fra le due realizzazioni ETRF89 e ETRF2000. Cioè, anche in questo caso, i valori da aggiungere alle coordinate espresse in un sistema per ottenere l’altro, sia per la parte planimetrica (latitudine e longitudine) sia per la parte altimetrica (differenze fra le altezze ellissoidiche ETRF2000 e quelle ETRF89).
Le griglie “K” sono anch’esse composte da 153 righe per 109 colonne e il passo è lo stesso delle griglie planimetriche descritte sopra. Le griglie “K” sono georeferenziate nel sistema ETRF89. La griglia altimetrica ha le stesse caratteristiche di georeferenziazione e passo di quelle planimetriche. Le coordinate del primo nodo in basso a sinistra sono le seguenti:
Origine delle griglie “K”: latitudine = 35°00'02.16" (35.0006°), longitudine = 5°57'07.92" (5.9522°)
I passaggi fra le coordinate ETRF2000 e gli altri sistemi si realizzano quindi sempre attraverso le coordinate ETRF89, applicando in sequenza la differenze fornite dalle diverse griglie: ad esempio per passare da ROMA40 a ETRF2000 si sommano i valori delle griglie fra ROMA40 e ETRF89 e poi quelli delle griglie fra ETRF89  e ETRF2000.
Per l’utilizzo nelle funzioni di calcolo, i valori delle griglie vengono memorizzati al momento della lettura dei file *.gr1- 2 e *.gk1-2, cioè quando viene chiamata la funzione ConveRgo_LOAD_GRI (che carica il contenuto di tutti i file dei grigliati contenuti nella cartella precedentemente indicata mediante la funzione ConveRgo_SET_PATH_GRI).
Per comodità, è possibile utilizzare la funzione scriviPathGridSuReg per copiare il percorso della cartella che contiene i grigliati in una chiave del registro di Windows, in modo da recuperarla nelle successive sessioni di utilizzo delle funzioni (chiave “GRID_FOLDER” in “HKEY_CURRENT_USER\Software\ConveRgo\Settings”). In questo modo, quando viene nuovamente eseguito il software che contiene le funzioni di calcolo, la funzione di inizializzazione ConveRgo_INIZIALIZZA può chiamare direttamente la funzione ConveRgo_SET_PATH_GRI e assegnare così automaticamente il percorso della cartella che contiene i grigliati.
I valori vengono archiviati in array bidimensionali di tipo double, uno per ciascuna componente delle differenze fra ogni coppia di sistemi, in gradi sessadecimali. Gli array bidimensionali corrispondono all’organizzazione matriciale (righe e colonne) delle griglie:
double matFiWgs[][] = differenze in latitudine fra Roma40 e ETRF89 double matLaWgs[][] = differenze in longitudine fra Roma40 e ETRF89 double matFiEd[][] = differenze in latitudine fra Roma40 e ED50 double matLaEd[][] = differenze in longitudine fra Roma40 e ED50
 
Un ulteriore array, delle stesse dimensioni dei precedenti ma di tipo char, contiene l’informazione dell’avvenuta assegnazione o meno del valore per ogni nodo:
char matMaskP[][] = maschera (1,2 = nodo assegnato, 0 = vuoto)
Normalmente, l’insieme dei grigliati a disposizione dell’utente copre solo una parte del territorio nazionale; gli array  di valori vengono quindi riempiti solo parzialmente. I nodi che rimangono vuoti mantengono il valore 0 nell’array di “maschera”, mentre i nodi assegnati assumono i valori 1 oppure 2, a seconda del tipo di grigliato da cui hanno preso l’informazione (gr1 oppure gr2).
Per i passaggi altimetrici da altezza ellissoidica ETRF89 a quota geoidica l’organizzazione dei dati è analoga: un array bidimensionale di tipo double per i valori di separazione e una maschera di tipo char per l’informazione sul riempimento. Le dimensioni degli array sono diverse dalle precedenti perché la griglia altimetrica è diversamente georeferenziata e ha una diversa densità:
double matHeight[][] = differenze in metri fra altezza ellissoidica e quota geoidica char matMaskH[][] = maschera (1,2 = nodo assegnato, 0 = vuoto)
Anche le griglie di tipo “K”, che contengono i valori per il passaggio fra ETRF89 e ETRF2000, sono organizzate allo stesso modo:
double matFiGrK[][] = differenze in latitudine fra ETRF89 e ETRF2000 double matLaGrK[][] = differenze in longitudine fra ETRF89 e ETRF2000
double matQuGrK[][] = differenze in metri fra altezza ellissoidica ETRF89 e ETRF2000 char matMaskK[][] = maschera (1 = nodo assegnato, 0 = vuoto)
In alcuni casi particolari, l’utente dispone di un modello di geoide “adattato”, cioè modificato ad hoc su una certa porzione di territorio. Normalmente tale informazione è disponibile in forma di file di testo contenente una lista di punti, già disposti secondo una maglia regolare, ciascuno dei quali esprime il nuovo valore di separazione fra ellissoide e geoide con cui sostituire in quel punto il modello nazionale.
L’estensione geografica dell’adattamento è ovviamente scelta di volta in volta, e va quindi trattata come variabile. Le funzioni permettono di memorizzare e utilizzare vari modelli adattati, archiviati in modo dinamico mediante la classe CMagliaGra e raccolti in un array di istanze di tale classe.
La funzione ConveRgo_LOAD_GRI carica direttamente i modelli di geoide adattati eventualmente presenti nella cartella dei grigliati, riconosciuti dall’estensione del file che deve essere “.GRA”.
Le funzioni per il passaggio della componente altimetrica da altezza ellissoidica a quota geoidica e viceversa richiedono l’indicazione del modello di geoide da utilizzare: adattato, nazionale o scelta automatica in base alla presenza o meno di adattamenti. In caso di eventuale presenza di più adattamenti che insistano sulla stessa area geografica (circostanza che comunque non si verifica in pratica), viene utilizzata l’informazione con indice di array più basso.
Nel caso che venga richiesto l’utilizzo degli adattamenti, le funzioni scorrono l’array dei modelli adattati e per ogni elemento verificano se le coordinate del punto da convertire ricadono o meno all’interno di quel modello. Il valore di ritorno delle funzioni indica il tipo di modello trovato in quel punto, cioè quello utilizzato per il calcolo.

In generale, tutte le griglie vengono interrogate con lo stesso criterio, che consiste nell’interpolazione bilineare all’interno della maglia. Le funzioni posToIndex forniscono, a partire dalle coordinate geografiche di un punto, gli indici di riga e colonna della maglia all’interno della quale ricade il punto e la percentuale di maglia compresa fra i nodi e la posizione del punto, che vengono poi utilizzati per l’interpolazione bilineare.

Alcune conversioni di coordinate richiedono la chiamata di più funzioni in successione. Ad esempio, per convertire una coppia di coordinate Gauss-Boaga in UTM-ETRF2000 si tratterà di chiamare in sequenza tre funzioni:
-	da piane Gauss-Boaga a geografiche ROMA40;
-	da geografiche ROMA40 a geografiche ETRF2000;
-	da geografiche ETRF2000 a piane UTM-ETRF2000;
fornendo in input alle chiamate successive alla prima il risultato della chiamata precedente.
 
Esempio:
double N_Roma,E_Roma,N_Etrf,E_Etrf; // coordinate piane nei due sistemi
double fiRoma,laRoma,fiEtrf,laEtrf; // coordinate geografiche (fi e lambda, sessadecimali) nei due sistemi
Funzioni da chiamare una sola volta all’inizio:
if (ConveRgo_SET_PATH_GRI(percorso_cartella_grigliati)<1) { gestire l’errore } if (1!=ConveRgo_LOAD_GRI()) { gestire l’errore }
Calcoli successivi, ad es. conversione da Gauss-Boaga fuso Ovest a UTM-ETRF2000 fuso 32:
if (1!=ConveRgo_PIA_GEO(N_Roma,E_Roma,32,ROMA,&fiRoma,&laRoma)) { gestire l’errore }
if (1!=ConveRgo_SISTEMI_GRI(fiRoma,laRoma,ROMA,ETRF2000,&fiEtrf,&laEtrf)) { gestire l’errore } if (1!=ConveRgo_GEO_PIA(fiEtrf,laEtrf,ETRF2000,32,&N_Etrf,&E_Etrf)) { gestire l’errore }


## License
[GNU Affero General Public License v3.0](https://choosealicense.com/licenses/agpl-3.0/)