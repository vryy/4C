/*!---------------------------------------------------------------------
\file
\brief guidelines to programming in ccarat

---------------------------------------------------------------------*/
/*! 
\addtogroup GUIDELINES 
*//*! @{ (documentation module open)*/

/*!---------------------------------------------------------------------
\brief Einige Hinweise zur Performance                                              

<pre>
Man hat beim Programmieren einen sehr grossen Einfluss darauf, ob
ein eine bestimmte Sache langsam oder schnell gerechnet wird, darum hier
einige Hinweise:

----*----
Generelles

Der Computer ist dumm und macht daher gerne dumme Sachen, d.h.

- Schnell abgearbeitet wird alles, was in moeglichst grossen Mengen 
in moeglichst EINFACHER Form abgespeichert ist.
(Das beste sind lange Vektoren)
  
- Der Compiler holt mit seiner Optimierung einen betraechtlichen 
Geschwindigkeitszuwachs aus dem Code, wenn er in der Lage ist
die Strukturen eines Programmes LOKAL zu sehen und die Operationen vorherzusagen
if/switch                                                           :  schlecht
Strukturen                                                          :  schlecht
Pointer                                                             :  sehr problematisch
Pointer auf Srukturen                                               :  schlecht
verschachtelte Daten                                                :  schlecht
lange Vektoren (mit fuer den Compiler sichtbarer Groesse)           :  sehr gut
starr dimensionierte Matrizen                                       :  sehr gut
dynamische Matrizen gut wenn Groesse fuer den Compiler sichtbar     :  gut
einfache, sehr grosse Funktionen                                    :  sehr gut
viele kleine Funktionen                                             :  schlecht
viele Flags, ueber die viele Entscheidungen gefaellt werden         :  schlecht
  

----*----
Man sollte sich in Erinnerung rufen, an welcher Stelle im Code man
gerade arbeitet, da dies bestimmt, wie haeufig etwas ausgefuehrt wird. 
Grundsaetzlich sind

in Steuerungsroutinen:

skalare Operationen            praktisch umsonst
vektorwertige Operationen      teuer
Aufrufe von Elementen          teuer
Aufrufe von solserv-Routinen   teuer
sparse matrix Operationen      teuer - sehr teuer
(Natuerlich kommt es noch darauf an, ob man etwas einmalig macht
 oder sich z.B. in einer Iterationsschleife befindet)
 
in Elementroutinen:            

skalare Operationen               preiswert
vektorwertige Operationen         o.k. wenn gut programmiert
dense matrix Operationen          o.k. wenn gut programmiert
Aufrufe von Eigenwert oder 
Loesern im Element                sehr teuer
Verzweigungen und switches        teuer und behindern den Compiler
kurzlaufende Schleifen            teurer als ausschreiben
Funktionsaufrufe                  teuer aber noetig
aufrufe von am-Funktionen      	  teuer
Operationen mit strings als flags teuer
Zugriffe in Strukturen direkt 
oder ueber Pointer		  teuer und behindern den Compiler


Am Integrationspunkt

skalare Operationen            preiswert in kleinen Mengen
vektorwertige Operationen      o.k. wenn gut programmiert
dense matrix Operationen       o.k. wenn gut programmiert
Aufrufe von Eigenwert oder 
Loesern im Element             sehr teuer
Verzweigungen und switches        teuer und behindern den Compiler
kurzlaufende Schleifen            teurer als ausschreiben
Funktionsaufrufe                  teuer aber noetig
aufrufe von am-Funktionen      	  teuer
Operationen mit strings als flags teuer
Zugriffe in Strukturen direkt 
oder ueber Pointer		  teuer und behindern den Compiler

Mehr Sorgfalt sollte also dort angewandt werden, wo Operationen haeufiger
ausgefuehrt werden. Als Belohnung fuer das Ausprobieren verschiedener
moeglicher Arten der Programmierung koennen u.U. massive Geschwindigkeitszuwaechse
erzielt werden.

----*----
skalare Operationen

Eine Division ist deutlich teuer als jede andere skalare Verknuepfung, daher
wenn moeglich Divisionen vorwegnehmen:
\code
anstatt
for (i=0; i<n; i++)
{
   b[i] = a[i]/2.0;
}
lieber
half = 1.0/2.0;
for (i=0; i<n; i++)
{
   b[i] = a[i]*half;
}
\endcode

----*----
vektor- und matrix-wertige Operationen

Koennen entweder mit Schleifen gemacht, ausgeschrieben oder unter Verwendung
von speziellen Serviceroutinen gemacht werden.
Grundsaetzlich gilt:

- Kurze Schleifen lieber ausschreiben: 
\code
anstatt
for (i=0; i<3; i++)
{
   b[i] = a[i]*c[i];
}
lieber
b[0] = a[0]*c[0];
b[1] = a[1]*c[1];
b[2] = a[2]*c[2];
\endcode
- Wenn schon Schleifen, dann moeglichst buendeln:
\code
anstatt
for (i=0; i<3; i++) b[i] = a[i]*c[i];
for (i=0; i<3; i++) d[i] = e[i]*f[i];
for (i=0; i<3; i++) g[i] = h[i]*j[i];
for (i=0; i<3; i++) n[i] = l[i]*k[i];
lieber
for (i=0; i<3; i++)
{
   b[i] = a[i]*c[i];
   d[i] = e[i]*f[i];
   g[i] = h[i]*j[i];
   n[i] = l[i]*k[i];
}
\endcode

- tensorielle Ausdruecke lieber in Matrixform programmieren
(Tensor 4.Stufe als 2D Matrix speichern)

- Bei komplexeren und/oder groesseren Vektoroperationen
auf eine Serviceroutine zurueckgreifen. Dies kann
eine von uns in der math Sammlung selbst geschriebene, jedoch bevorzugt
eine 'geklaute' Routine sein. Bei groesseren Matrizen (so ab Dimension 30)
lohnt sich der Aufruf von BLAS Routinen, die als lib sowieso dazugelinkt sind.
(BLAS Aufrufe findet man z.B. unter netlib.org erlaeutert)

- Serviceroutinen fuer matrix-vektor Operationen sollten bevorzugt in Fortran sein.

- Gelegentlich ist es auch moeglich 2D Matrizen als Vektoren zu interpretieren
und die Operationen mit einem Vektor auszufuehren anstatt mit einer Matrix
\code
anstatt
for (i=0; i<10; i++)
for (j=0; j<10; j++) 
{
   A[i][j] = A[i][j] * xi;
}
lieber
for (i=0; i<100; i++)
{
  A[0][i] = A[0][i] * xi;
} 
\endcode
Dies geht insbesondere mit 2D-4D Matrizen, die mit dem am-System erzeugt wurden auch,
da diese im letzten Laufindex physikalisch kontinuierlich sind.

- logische Abbrueche und Verzweigungen mit if oder switch 
in Schleifen vermeiden da diese den Compiler an der Optimierung der 
Schleife hindern:
\code
anstatt
for (i=0; i<10; i++)
{
   A[i] = b[i] * c[i]; 

   if (FABS(A[i])<EPS12) break;
}
lieber
for (i=0; i<10; i++)
{
   A[i] = b[i] * c[i]; 
}
for (i=0; i<10; i++)
if (FABS(A[i])<EPS12) break;
\endcode

----*----
logische Steuerungen des Programmablaufs

Grundsaetzlich ist schlichter Code schneller als komplizierter, soll heissen:

- Auf Integrationspunkt- und Elementebene moeglichst wenig if und switch

- Auf Integrationspunkt- und Elementebene haben Aufrufe von 
hoeheren c-Bibliothekroutinen wie printf, strcmp, o.ae. nichts verloren
da diese richtig teuer sind/sein koennen / die Optimierung des Compilers stoeren

- Auf Integrationspunkt- und Elementebene sollten Felder statische Groessen haben und
entweder Lokalvariablen fester Groesse sein, oder mittels am-System EINMAL 
als static deklariert werden.
Das Anlegen und Loeschen von dynamischen Matrizen ist hier sauteuer (und sehr fehlertraechtig). 

- Funktionen sollten spezialisiert sein, also eine extra Funktion fuer jede spezifische
Aufgabe schreiben, anstatt einer grossen, die ueber n flags n*n verschiedene Dinge tun kann.
  
- Andererseits Funktionen, die nur 'dumm' etwas ausrechnen, ohne Entscheidungen zu faellen
nicht in viele kleine Subroutinen splitten. Fuer den Compiler ist es besser, wenn moeglichst viel
in einer Routine steht. 

----*----
selbstverstaendliches, was praktisch doch immer wieder auftaucht

- Nicht die gleichen Dinge mehrfach ausrechnen

- Nicht Sachen ausrechnen, die man nicht haben moechte, und die fuer den
Fortgang einer Simulation auch keine Bedeutung haben

- Nur Daten auf Platte schreiben, die man gedenkt hinterher auch anzukucken und auszuwerten
(Kostet Platz UND ist langsam)

- Kollegen fragen, wenn man sich nicht sicher ist, wie man etwas umsetzen soll

- etc...

</pre>
------------------------------------------------------------------------*/
void performance()
{
#ifdef DEBUG 
dstrc_enter("performance");
#endif
/*----------------------------------------------------------------------*/
/* Code */
/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frint */



/*!---------------------------------------------------------------------
\brief general appearance of a routine                                              

<pre>                                                    author month/year 
There is a very detailed description of the routine here
</pre>
\param firstvar  typ  (i/o) description of firstvar
\param scndvar   typ  (i)   description of scndvar
\param thrdvar   typ  (o)   description of thrdvar

\warning There is something special to this routine
\return void                                               
\sa routine1() , routine2() , routine3() , structure1 , variable1                                    

<pre>
Der Quellcode, so dass er von doxygen richtig erfasst wird sieht so aus:
\include example3.c
wobei nach dem Beginn-Zeichen des Kommentars noch ein '!' folgen muss, um
den Kommenar wirksam werden zu lassen (Hier nicht dargestellt, um Interpretation
durch doxygen zu verhindern).
siehe auch www.doxygen.org fuer weiteres 
</pre>
------------------------------------------------------------------------*/
void general_appearance(INT *var,INT num, INT *ierr)
{
#ifdef DEBUG 
dstrc_enter("general_appearance");
#endif
/*----------------------------------------------------------------------*/


/* Code */


/*----------------------------------------------------------------------*/
#ifdef DEBUG 
dstrc_exit();
#endif
return;
} /* end of frint */





/*!---------------------------------------------------------------------
\brief commenting the code                                              

<pre>                                                   
Kommentare sind fuer alle ausgesprochen wichtig. Auch man selbst weiss haeufig
nach 2 Wochen nicht mehr, was man an einer Stelle gemacht hat und ist daher
froh, wenn man an der Stelle einen Kommentar geschrieben hat. 
Ich denke jeder, der von einem Vorgaenger mal ein grosses, schlecht kommentiertes
aber kompliziertes Stueck Code geerbt hat versteht die Bedeutung der Kommentare
speziell in einem Forschungscode wie unserem. (Man denke sich nur den Sockel
dieses Programms als unkommentiert - er waere wertlos, da kein Mitarbeiter
auf diesen Sockel aufsetzen koennte)

Daher:
-Kommentare koennen gar nicht reichlich, praezise und aufuehrlich genug sein
-Nur eine Datei, die zu mindestens 50% aus Kommentar besteht ist eine gute Datei
-Querverweise zur Literatur und Theorie sind ausgesprochen erwuenscht
-Formeln und Belegungplaene fuer Felder durfen ruhig grafisch aufwendig in
 Kommentaren verewigt werden
-Doxygen Kommandos und html Kommandos in Kommentaren benutzen um eine automatische
 Dokumentationserstellung zu ermoeglichen
  
Beispiele aus dem Sockel (ja, dort findet man auch schlecht kommentierten Code):

<b>Positivbeispiel:</b>
  
   rhs[3]    original load vector
   rhs[2]             load vector at time t-dt
   rhs[1]             load vector at time t
   rhs[0]    interpolated load vector and working array

   fie[2]    internal forces at step t
   fie[1]    internal forces at step t-dt
   fie[0]    interpolated internal forces and working array

   dispi[0]  displacement increment from t-dt to t

   sol[0]    total displacements at time t-dt
   sol[1]    total displacements at time t
   
   vel[0]    velocities    at t-dt
   acc[0]    accelerations at t-dt

   work[2]   working vector for sums and matrix-vector products 
   work[1]   working vector for sums and matrix-vector products 
   work[0]   working vector for sums and matrix-vector products 
   work[0]   is used to hold residual displacements in corrector 
             iteration
             
   in the nodes, displacements are kept in node[].sol[0][0..numdf-1]
                 velocities    are kept in node[].sol[1][0..numdf-1]
                 accelerations are kept in node[].sol[2][0..numdf-1]

Values of the different vectors from above in one loop:
  /......no change in this step
  =,+=...evaluation in this step
  {}....takes values from another vector
vector       Predictor - Start     Precictor - End     Corrector - Start     Corrector - End                 Update - End

rhs[3]       /{=orig. load vect.}    /                       /               /                               /
rhs[2]       /{=rhs(t-dt)}           /                       /               /                               =rhs[1]{=rhs(t)}
rhs[1]       =rhs(t)                 /                       /               /                               /
rhs[0]       /{=rhs(t-2dt)}          =feff_p                 /               =feff_c                         =rhs[2]{=rhs(t-dt)}

fie[2]       /                       /                       =fint(t)        /                               /
fie[1]       =fint(t-dt)             /                       /               /                               /
fie[0]       /                       /                       /               =(1-alpha_f)*fie[2]+alpha_f*fie[1]      /

dispi[0]     =0                      =Keff^-1*feff-p         /               +=work[0]                       /

sol[0]       /{=d(t-dt)}             /                       /               /                               =sol[1]{=d(t)}
sol[1]       {=d(t-dt)}              =sol[0]+dispi[0]{=d(t)} /               =sol[0]+dispi[0]                /

vel[0]       /{=v(t-dt)}             /                       /               /                               =v(t)
acc[0]       /{=a(t-dt)}             /                       /               /                               =a(t)   

work[2]      /{=v(t-2dt)}            /                       /               /                               =v(t-dt)
work[1]      /{=a(t-2dt)}            /                       /               /                               =a(t-dt)
work[0]      /                       /                       /               =Keff^-1*feff-c                 =M*vel[0]


<b>Negativbeispiel:</b>
\code
switch (array_from->Typ)
{
case D3:
   dim = array_from->fdim * array_from->sdim * array_from->tdim;
   am4def(array_from->name,
          array_to,
          array_from->fdim,
          array_from->sdim,
          array_from->tdim,
          0,
          "D3");
   dptr_from = array_from->a.d3[0][0];
   dptr_to   = array_to->a.d3[0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break; 
case D4:
   dim = array_from->fdim * array_from->sdim * array_from->tdim * 
         array_from->fodim;
   am4def(array_from->name,
          array_to,
          array_from->fdim,
          array_from->sdim,
          array_from->tdim,
          array_from->fodim,
          "D4");
   dptr_from = array_from->a.d4[0][0][0];
   dptr_to   = array_to->a.d4[0][0][0];
   for (i=0; i<dim; i++) *(dptr_to++) = *(dptr_from++);
break;
}
\endcode
</pre>

------------------------------------------------------------------------*/
void commenting_code()
{
};  



/*!---------------------------------------------------------------------

\brief style of c-Syntax                                              


<pre>                                                   
<b>Man kann C so scheiben:</b>
\code
if (iptr != NULL)
{
  for (i=0; i<n; i++) iptr[i]++;
}
\endcode
<b>oder so:</b>
\code
iptr?(for (i=0;i<n;(*(iptr++))++)):;
\endcode
da sollte man einen Mittelweg finden , damit es der naechste,
der auch nicht als Informatiker auf die Welt gekommen ist,
auch noch versteht. Deswegen einige Vorgaben zur Gestalt von 
Basistrukturen:
\include example1.c
</pre>

------------------------------------------------------------------------*/
void csyntaxstyle()
{
}

/*!---------------------------------------------------------------------

\brief in general....                                              


<pre>                                                   
- keine Funktionspointer verwenden, das ruiniert die Lesbarkeit und das Debuggen
- keine Bitmuster verwenden

bitte kein 'Schweinskram' mit Datentypen:

<b>Negativbeispiel:</b>
\code
double a;
int *iptr;
iptr = &((int)a);
*iptr=5;
iptr++;
*iptr=6;
\endcode

Rueckgabewert sollte immer void sein, ausser es gibt einen vernuenftigen
Grund hiervon abzuweichen.Rueckgabewert sollte im speziellen nie INT sein 
(=Defaultwert, wird vom Compiler nicht geprueft)

von JEDER Routine muss es einen Prototyp geben, ausgenommen sind nur Aufrufe
von Fortran Routinen
 Prototyp:
\code
void frint_n(CHAR string[],INT *var,INT num, INT *ierr);
\endcode
eigentliche Routine:
\code
void frint_n(CHAR string[],INT *var,INT num, INT *ierr)
\endcode
Wo der Prototyp eingefuegt wird ist zweckabhaengig:

Routinen, die nur von einer anderen Routine der gleichen Datei 
benoetigt werden haben ihren Prototyp im Kopf der Datei und
sind vom Typ static (damit von aussen unsichtbar):
Prototyp der gekapselten Routine:
\code
static void routine2(void);
\endcode
routine1 ruft als einziges routine2, beide sind in der gleichen Datei 
\code
static void routine2()
{
return;
}
void routine1()
{
routine2();
return;
}
\endcode

Routinen, die klar abgegrenzt zu einem Modul gehoeren (z.B. zu einem Finiten Element) 
haben ihren Prototyp in dem modulspezifischen header, der nur in Dateien dieses
Moduls eingefuegt wird:
in shell8.h:
\code
void s8_routine1(void);
void s8_routine2(void);
\endcode
in Datei s8_rout1.c:
\code
#include "shell8.h"
void s8_routine1()
{
return;
}
\endcode
in Datei s8_rout2.c:
\code
#include "shell8.h"
void s8_routine2()
{
return;
}
\endcode
Routinen, die Strukturen aus dem Header solution.h verwenden
haben in irgendeiner Form Kontakt zu Loeserbibliotheken und
haben ihren Prototyp in prototypes_sol.h :
in prototypes_sol.h:
\code
void solver_control(struct _SOLVAR       *actsolv,
                    struct _INTRA        *actintra,
                    enum   _SPARSE_TYP   *sysarray_typ,
                    union  _SPARSE_ARRAY *sysarray,
                    struct _DIST_VECTOR  *sol,
                    struct _DIST_VECTOR  *rhs,
                    INT                   option);
\endcode
in solver_control.c
\code
#include "../headers/standardtypes.h"
#include "../headers/solution.h"
void solver_control(struct _SOLVAR       *actsolv,
                    struct _INTRA        *actintra,
                    enum   _SPARSE_TYP   *sysarray_typ,
                    union  _SPARSE_ARRAY *sysarray,
                    struct _DIST_VECTOR  *sol,
                    struct _DIST_VECTOR  *rhs,
                    INT                   option)
{
return;
}               
\endcode
Routinen, die keine der vorgenannten Bedingungen erfuellen
haben ihren Prototyp in prototypes.h und sind damit von ueberallaus
sichtbar
Diese Moeglichkeit sollte auf die notwendigsten Routinen beschraenkt
bleiben! Der Sockel ist diesbezueglich noch kein alzu grosses Vorbild,
einige Routinen wie z.B. alle input-Routinen kann man noch kapseln

</pre>

------------------------------------------------------------------------*/
void in_general(){};








/*!---------------------------------------------------------------------

\brief global variables                                             

<pre>                                                    
Man sollte je nach Verwendungszweck einer globalen Variablen pruefen, wo diese
deifniert werden muss. Grundsaetzlich gilt, dass eine globale Variable nur in
jenen Dateien sichtbar werden darf, in denen sie auch verwendet wird und nicht
im gesamten Programm definiert ist.

globale Variablen die nur in einer Datei gebraucht werden, werden 
im Kopf dieser datei als static deklariert:
Datei test.c:
\code
static INT nurintest;
void routine1()
{
nurintest=5;
return;
}
\endcode

globale Variablen von 'lokaler Bedeutung', die in mehreren Dateien
gebraucht werden, werden in einer datei definiert und in den anderen dateien
als 'extern' deklariert:
Datei test1.c (hier ist intestmodul definiert):
\code
INT intestmodul;
void routine1()
{
intestmodul=5;
return;
}
\endcode
Datei test2.c (hier ist intestmodul nur 'extern', d.h. als vereis auf test1.c definiert):
\code
extern INT intestmodul;
void routine2()
{
intestmodul+=3;
return;
}
\endcode

globale Variablen von 'groesserer Bedeutung' werden im Kopf derjenigen
Datei definiert, in der diese Variable hauptsaechlich eine Rolle spielt,
in allen anderen Dateien, die diese Variable brauchen wird sie als 'extern'
definiert und es wird im standardtypes.h ein Kommentar zu dieser Variablen
angelegt:
\include example2.c
Wobei nur das 'Original' in global_control.c einen Doxygen-wirksamen Kommentar
enthalten sollte, also ein ! nach Kommentaranfang (hier nicht dargestellt)
</pre>

------------------------------------------------------------------------*/
void global_variables(){}



/*! @} (documentation module close)*/





