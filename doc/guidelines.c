/*!---------------------------------------------------------------------
\file
\brief guidelines to programming in ccarat

---------------------------------------------------------------------*/
/*! 
\addtogroup GUIDELINES 
*//*! @{ (documentation module open)*/

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





