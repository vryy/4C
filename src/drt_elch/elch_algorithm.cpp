
#ifdef CCADISCRET

#include "elch_algorithm.H"

//#include "../drt_lib/drt_globalproblem.H"
//#include "../drt_lib/drt_validparameters.H"
//#include <Teuchos_StandardParameterEntryValidators.hpp>


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::Algorithm(Epetra_Comm& comm)
  :  FluidBaseAlgorithm(),
     comm_(comm)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::Algorithm::~Algorithm()
{
}


#endif // CCADISCRET
