/*!----------------------------------------------------------------------
\file
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "Epetra_SerialDenseMatrix.h"
#include "global_inp_control2.H"



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

#ifdef DEBUG
/*!----------------------------------------------------------------------
  \brief the tracing variable

  <pre>                                                         m.gee 8/00
  defined in pss_ds.c, declared in tracing.h
  </pre>
 *----------------------------------------------------------------------*/
extern struct _CCA_TRACE         trace;
#endif

/*----------------------------------------------------------------------*
  |                                                       m.gee 06/01    |
  | general problem data                                                 |
  | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;

/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate design if needed                                 |
 | defined in global_control.c                                          |
 *----------------------------------------------------------------------*/
extern struct _DESIGN *design;

/*!----------------------------------------------------------------------
\brief ranks and communicators

<pre>                                                         m.gee 8/00
This structure struct _PAR par; is defined in main_ccarat.c
and the type is in partition.h
</pre>

*----------------------------------------------------------------------*/
 extern struct _PAR   par;

#if 0

// Obsolete now that we always read in jumbo mode.

/*----------------------------------------------------------------------*
  |                                                        m.gee 10/06  |
  | replace mpi comms in discretizations                                |
  | make design redundant                                               |
  | distribute grids                                                    |
 *----------------------------------------------------------------------*/
void distribute_drt_grids()
{
  //the comm in the design can stay MPI_COMM_WORLD
  // get discretizations and replace their comm
  for (int i=0; i<genprob.numfld; ++i)
  {
    vector<RefCountPtr<DRT::Discretization> >* discretization =
      (vector<RefCountPtr<DRT::Discretization> >*)field[i].ccadis;

    for (int j=0;j<field[i].ndis;j++)
    {
      RefCountPtr<DRT::Discretization> actdis = (*discretization)[j];
      int ierr=0;
#ifdef PARALLEL
      INTRA* actintra = &(par.intra[i]);
      RefCountPtr<Epetra_MpiComm> comm =
                           rcp(new Epetra_MpiComm(actintra->MPI_INTRA_COMM));
      actdis->SetComm(comm);
      ierr = actdis->FillComplete();
      if (ierr) dserror("FillComplete returned %d",ierr);
#endif
      if (!actdis->Filled()) ierr = actdis->FillComplete();
      if (ierr) dserror("FillComplete returned %d",ierr);

      // partition the discretization using metis
      // create nodal graph of problem
      RefCountPtr<Epetra_CrsGraph> nodegraph = actdis->BuildNodeGraph();

      // repartition the nodal graph using metis
      Epetra_Vector weights(nodegraph->RowMap(),false);
      weights.PutScalar(1.0);
      RefCountPtr<Epetra_CrsGraph> newnodegraph =
                          DRT::Utils::PartGraphUsingMetis(*nodegraph,weights);
      nodegraph = null;

      // the rowmap will become the new distribution of nodes
      const Epetra_BlockMap rntmp = newnodegraph->RowMap();
      Epetra_Map newnoderowmap(-1,rntmp.NumMyElements(),rntmp.MyGlobalElements(),0,actdis->Comm());

      // the column map will become the new ghosted distribution of nodes
      const Epetra_BlockMap cntmp = newnodegraph->ColMap();
      Epetra_Map newnodecolmap(-1,cntmp.NumMyElements(),cntmp.MyGlobalElements(),0,actdis->Comm());
      newnodegraph = null;

      // do the redistribution
      actdis->Redistribute(newnoderowmap,newnodecolmap);

    } // for (int j=0;j<field[i].ndis;j++)
  } // for (int i=0; i<genprob.numfld; ++i)

  return;
} // void distribute_drt_grids

#endif

#endif  // #ifdef CCADISCRET
