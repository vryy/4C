for (i=0; i<n; i++)
{
   for (j=0; j<m; j++)
   {

   }/* end of for (j=0; j<m; j++) */
}/* end of for (i=0; i<n; i++) */

/* oder */
for (i=0; i<n; i++){
   for (j=0; j<m; j++){

   }/* end of for (j=0; j<m; j++) */
}/* end of for (i=0; i<n; i++) */

/* oder */
for (i=0; i<n; i++)
   for (j=0; j<m; j++){

   }/* end of for (j=0; j<m; j++) and end of for (i=0; i<n; i++) */


/*---------------------------------------------------------------Verzweigungen*/
if (a <= b)
{
    a += 6;
} /* end of if (a <= b) */
else
{
    a += 7;
} /* end of else to if (a <= b) */

/* oder */
if (a <= b)  a += 6;
else         a += 7;

/*-------------------------------------------------------------------Switches */
switch(b)
{
case 1:
   b += 5;
break;
case 2:
   b += 4;
break;
default:
   dserror("Unknown case for b");
break;
}/* end of switch(b) */

/*-------------------------------------- Strukturen und Zeiger auf Strukturen */
/* aus Performancegruenden keine Megastapel von Strukturen 
  beim Rechnen verwenden: 
*/
/* nicht: */
field[i].dis[0].element[k].node[0]->dof[3] = 5;
/* sondern */
NODE *actnode;
actnode = &(field[i].dis[0].element[k].node[0]);
actnode->dof[3] = 5;

/*----------------------------------------------------------------Datentypen */
/*
Als Datentypen sind nur
*/
INT
DOUBLE
CHAR
enum
struct
union  
/*
zulaessig, nachstehende bitte nicht benutzen:
*/
int
double
float
long int
short char
/* 
und andere in ihrer Groesse manipulierte Typen. 
Diese fuehren bei der Portierung auf andere Plattformen zu u.U. grossen
Schwierigkeiten.
*/
/* 
eigene Datentypen (ausser struct, enums und unions) bitte nur in
begruendeten Faellen definieren, und zwar im header definitions.h
z.B. :
*/
typedef int PTRSIZE;
/*
Niemals pointer umdefinieren, das macht den Code vollkommen unverstaendlich
und fuehrt zwangsweise dazu, dass Fehler in der Allokierung und Freigabe von
Speicher gemacht werden.Negativbeispiel:
*/
typedef ELEMENT* EPTR;


EPTR a;/* pointereigenschaft hier unsichtbar */
a->numnp = 5;


/* Positivbeispiel */
ELEMENT *a; /* sichtbar ein Pointer */
a->numnp=5;

/*----------------------------------------------calloc, malloc, realloc, free */
/*
Bitte die Basisroutinen
*/
malloc,
calloc,
realloc,
free
/* nicht benuetzen, sondern */
MALLOC,
CALLOC,
REALLOC,
FREE
/* 
Dies sind Huellroutinen, die eine Ueberwachung des Speichers im DEBUG-exe
unterstuetzen, und ansonsten die identische Funktionalitaet haben 
*/

/*-----------------------------------------dynamische mehrdimensionale Felder */
/*
Moeglichst davon abstand nehmen mehrdimensionale Felder selbst zu allokieren,
da in den AM-Routinen hierfuer bereits optimierte Algorithmen zur Ver-
fuegung stehen, und mittels AM auch sichergestellt wird, dass aller Speicher
auch wieder korrekt freigegeben wird. Auch zur Initialisierung und anderes
bitte AM-Routinen benutzen. Die Vorteile sind
- Fehlerfreiheit
- Geschwindigkeit aufgrund guter Algorithmen in AM
- gute Lesbarkeit fuer alle Mitarbeiter
Der Funktionsumfang von AM kann natuerlich noch erweiter werden
*/
/* Negativbeispiel: */
DOUBLE **a;
a = (DOUBLE**)CALLOC(10,sizeof(DOUBLE*));
for (i=0; i<10; i++) a[i] = (DOUBLE*)CALLOC(10,sizeof(DOUBLE));

/* Positivbeispiel: */
ARRAY a_array;
DOUBLE **a;
a = amdef("name",&a_array,10,10,"DA");
amdel(&a_array);

/*---------------------------------------------Fehlerkontrolle, sicherer Code */
/*
Mit dem DS-System sind einige maechtige Tools zur Fehlervermeidung gegeben, diese
bitte intensiv benutzen:
*/
/*
Fehlermeldung
*/
if (a==0) dserror("a is zero");
/*
Es gibt noch die Moeglichkeit Werte 'sicherzustellen'. Dies kann intensivst benutzt
werden, da die Routine dsassert bei einem nicht-DEBUG-exe verschwindet
*/
/* stelle sicher, dass a nicht 0 ist, tue dies nur im DEBUG-exe, 
   nicht im fast-exe */
dsassert(a!=0,"a is zero");
}
