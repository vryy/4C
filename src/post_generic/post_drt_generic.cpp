/*!
 \file post_drt_generic.cpp

 \brief main routine of a generic filter

 <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

 */



#include "post_drt_generic_writer.H"
#include "post_drt_generic_single_field_writers.H"

/*!
 \brief filter main routine

 Write binary ensight format.

 The ens_checker that is part of ensight can be used to verify the
 files generated here.

 \author u.kue
 \date 03/07
 */
 
int main(
        int argc,
        char** argv)
{
    Teuchos::CommandLineProcessor My_CLP;
    My_CLP.setDocString("Post DRT ensight Filter\n");

    PostProblem problem(My_CLP, argc, argv);

#if 0
    for (int i = 0; i<problem.num_discr(); ++i)
    {
        PostField* field = problem.get_discretization(i);
        StructureEnsightWriter writer(field, problem.outname());
        writer.WriteFiles();
    }
#endif

    try
    {

    // each problem type is different and writes different results
    switch (problem.Problemtype())
    {
    case prb_structure:
    {
        PostField* field = problem.get_discretization(0);
        StructureEnsightWriter writer(field, problem.outname(), problem.stresstype(), problem.straintype());
        writer.WriteFiles();
        break;
    }
    default:
        dserror("problem type %d not yet supported", problem.Problemtype());
    } // switch

    }
    catch ( std::runtime_error & err )
    {
      char line[] = "=========================================================================\n";
      std::cout << "\n\n"
                << line
                << err.what()
                << "\n"
                << line
                << "\n" << std::endl;
#ifdef DSERROR_DUMP
      abort();
#endif

#ifdef PARALLEL
      MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
#else
      exit(1);
#endif
    }

    return 0;
}

