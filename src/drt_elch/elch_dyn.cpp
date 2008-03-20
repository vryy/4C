
#ifdef CCADISCRET

#include <string>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "elch_dyn.H"
#include "elch_algorithm.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
// entry point for ELCH in DRT
/*----------------------------------------------------------------------*/
void elch_dyn()
{
    // create a communicator
#ifdef PARALLEL
    Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm comm;
#endif

    // print out the ELCH module logo
    if (comm.MyPID() == 0)
    {
        // more at http://www.ascii-art.de under entry "moose" (or "elk")
        cout<<"     ___            ___    "<<endl;
        cout<<"    /   \\          /   \\ "<<endl;
        cout<<"    \\_   \\        /  __/ "<<endl;
        cout<<"     _\\   \\      /  /__  "<<"     _____ _     _____  _   _   "<<endl;
        cout<<"     \\___  \\____/   __/  "<<"    |  ___| |   /  __ \\| | | |  "<<endl;
        cout<<"         \\_       _/      "<<"   | |__ | |   | /  \\/| |_| |  "<<endl;
        cout<<"           | @ @  \\_      "<<"   |  __|| |   | |    |  _  |   "<<endl;
        cout<<"           |               "<<"  | |___| |___| \\__/\\| | | | "<<endl;
        cout<<"         _/     /\\        "<<"   \\____/\\_____/\\____/\\_| |_/ "<<endl;
        cout<<"        /o)  (o/\\ \\_     "<<endl;
        cout<<"        \\_____/ /         "<<endl;
        cout<<"          \\____/          "<<endl;
        cout<<"                           "<<endl;
    }
    // create an ELCH::Algorithm instance
    Teuchos::RCP<ELCH::Algorithm> elch = Teuchos::rcp(new ELCH::Algorithm(comm));
    // solve the electrochemistry problem
    elch->TimeLoop();

    // summarize the performance measurements
    Teuchos::TimeMonitor::summarize();

} // elch_dyn()


#endif  // CCADISCRET
