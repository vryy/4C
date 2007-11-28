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

/*
 \brief Writer for structural problems
 */
class StructureEnsightWriter : public EnsightWriter
{
public:
    StructureEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual string WriteAllResults(
            PostField* field);
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string StructureEnsightWriter::WriteAllResults(
        PostField* field)
{
    stringstream str;
    str << EnsightWriter::WriteResult("displacement", "displacement", field->problem()->num_dim());
    str << EnsightWriter::WriteResult("velocity", "velocity", field->problem()->num_dim());
    str << EnsightWriter::WriteResult("acceleration", "acceleration", field->problem()->num_dim());
    return str.str();
}

/*----------------------------------------------------------------------*/
/*
 \brief Writer for fluid problems
 */
/*----------------------------------------------------------------------*/
class FluidEnsightWriter : public EnsightWriter
{
public:
    FluidEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual string WriteAllResults(
            PostField* field);
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string FluidEnsightWriter::WriteAllResults(
        PostField* field)
{
    stringstream str;
    str << EnsightWriter::WriteResult("velnp", "velocity", field->problem()->num_dim());
    str << EnsightWriter::WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
    str << EnsightWriter::WriteResult("residual", "residual", field->problem()->num_dim());
    str << EnsightWriter::WriteResult("dispnp", "displacement", field->problem()->num_dim());
    str << EnsightWriter::WriteResult("traction", "traction", field->problem()->num_dim());
    return str.str();
}

/*----------------------------------------------------------------------*/
/*
 \brief Writer for ale problems
 */
/*----------------------------------------------------------------------*/
class AleEnsightWriter : public EnsightWriter
{
public:
    AleEnsightWriter(
            PostField* field,
            string filename) :
        EnsightWriter(field, filename)
    {
    }

protected:

    virtual string WriteAllResults(
            PostField* field);
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string AleEnsightWriter::WriteAllResults(
        PostField* field)
{
    stringstream str;
    str << EnsightWriter::WriteResult("dispnp", "displacement", field->problem()->num_dim());
    return str.str();
}

/*----------------------------------------------------------------------*/
/*
 \brief Writer for fluid problems including xfem interfaces
 
 this might be handled more general in the near future, when integration cells
 are saved along the data for each timestep
 
 At the moment, we just intersect in the post filter again and conlude on the DOF 
 distribution which then should fit to the unknown vector. 
 hence, the xfem post filter is doing quite some work at this moment...
 
 \author a.ger
 \date 11/07
 
 */
/*----------------------------------------------------------------------*/
class XFluidEnsightWriter : public EnsightWriter
{
public:
    XFluidEnsightWriter(
            PostField* field,           ///< Field to be plotted
            PostField* cutterfield,     ///< intersecting field (for fsi)
            string filename) :
        EnsightWriter(field, filename),
        cutterfield_(cutterfield)
    {
    }

protected:

    virtual string WriteAllResults(
            PostField* field);
    
//    /*!
//     \brief write all time steps of a result for xfem problems with varying number of unknowns per node
//
//     Write nodal results. The results are taken from a reconstructed
//     Epetra_Vector. We always assume a dofmanager that knows about the unknowns distribution
//     So we have to specify which XFEM::Physics::Field.
//
//     Finally, after writing to the result file, a string is returned that
//     describes the result for the case file VARIABLE section
//
//     \return string with entry for VARIABLE section in case file
//     \author a.ger
//     \date 11/07
//     */
//    string WriteResult(
//            const string groupname,                   ///< name of the result group in the control file
//            const string name,                        ///< name of the result to be written
//            const vector<XFEM::Physics::Field> fields ///< Fieldnames to print, could be just one (pressure) or a list (Dispx,Dispy,Dispz)
//            );
    
    
    
    /*!
     \brief write all time steps of a result

     Write nodal results. The results are taken from a reconstructed
     Epetra_Vector. In many cases this vector will contain just one
     variable (displacements) and thus is easy to write as a whole. At
     other times, however, there is more than one result (velocity,
     pressure) and we want to write just one part of it. So we have to
     specify which part.

     Finally, after writing to the result file, a string is returned that
     describes the result for the case file VARIABLE section

     \return string with entry for VARIABLE section in case file
     \author u.kue
     \date 03/07
     */
    string WriteResult(
            const string groupname, ///< name of the result group in the control file
            const string name, ///< name of the result to be written
            const int numdf, ///< number of dofs per node to this result
            const int from=0 ///< start position of values in nodes
            );
    
    
    
    
    
    
    PostField* cutterfield_;
};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string XFluidEnsightWriter::WriteAllResults(
        PostField* field)
{
    stringstream str;
    str << WriteResult("velnp", "velocity", field->problem()->num_dim());
    str << WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
    str << WriteResult("residual", "residual", field->problem()->num_dim());
    return str.str();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string XFluidEnsightWriter::WriteResult(
        const string groupname,
        const string name,
        const int numdf,
        const int from)
{
    // maximum file size
    const unsigned FILE_SIZE_LIMIT = 0x7fffffff; // 2GB

    PostResult result = PostResult(field_);
    result.next_result();
    if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
        return "";

    // new for file continuation
    unsigned int fileset = 1;
    const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
    const int numnp = nodemap->NumMyElements();
    const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);

    string filename = filename_ + "_"+ field_->name() + "."+ name;
    ofstream file(filename.c_str());

    WriteResultStep(file, result, resultfilepos_, groupname, name, numdf, from);
    while (result.next_result())
    {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT)
        {
            FileSwitcher(file, fileset, filesets_, resultfilepos_, stepsize, name, filename);
        }
        WriteResultStep(file, result, resultfilepos_, groupname, name, numdf, from);
    }

    // append index table
    WriteIndexTable(file, resultfilepos_[name]);
    resultfilepos_[name].clear();

    if (fileset != 1)
    {
        filesets_[fileset-3].push_back(file.tellp()/stepsize);
        filename += "***";
    }

    return GetVariableEntryForCaseFile(numdf, fileset, name, filename);
}


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
