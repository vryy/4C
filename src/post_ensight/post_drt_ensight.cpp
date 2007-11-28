/*!
 \file post_drt_ensight.cpp

 \brief main routine of the Ensight filter

 <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

 */

#ifdef CCADISCRET

#include "post_drt_ensight_writer.H"
#include "post_drt_ensight_single_field_writers.H"

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

    // each problem type is different and writes different results
    switch (problem.Problemtype())
    {
    case prb_fsi:
    {
        string basename = problem.outname();
        PostField* structfield = problem.get_discretization(0);
        StructureEnsightWriter structwriter(structfield, basename);
        structwriter.WriteFiles();

        PostField* fluidfield = problem.get_discretization(1);
        FluidEnsightWriter fluidwriter(fluidfield, basename);
        fluidwriter.WriteFiles();

        PostField* alefield = problem.get_discretization(2);
        AleEnsightWriter alewriter(alefield, basename);
        alewriter.WriteFiles();
        break;
    }
    case prb_structure:
    {
        PostField* field = problem.get_discretization(0);
        StructureEnsightWriter writer(field, problem.outname());
        writer.WriteFiles();
        break;
    }
    case prb_fluid:
    {
        PostField* field = problem.get_discretization(0);
        FluidEnsightWriter writer(field, problem.outname());
        writer.WriteFiles();
        break;
    }
    case prb_ale:
    {
        PostField* field = problem.get_discretization(0);
        AleEnsightWriter writer(field, problem.outname());
        writer.WriteFiles();
        break;
    }
    case prb_fluid_xfem:
    {
        cout << "Output XFEM Problem" << endl;
        
        PostField* structfield = problem.get_discretization(0);
        StructureEnsightWriter structwriter(structfield, problem.outname());
        structwriter.WriteFiles();

        PostField* fluidfield = problem.get_discretization(1);
        XFluidEnsightWriter fluidwriter(fluidfield, structfield, problem.outname());
        fluidwriter.WriteFiles();
        break;
    }
    default:
        dserror("problem type %d not yet supported", problem.Problemtype());
    }

    return 0;
}

#endif
