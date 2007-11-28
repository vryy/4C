/*!
 \file post_drt_ensight_single_field_writers.cpp

 \brief main routine of the Ensight filter

 <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

 */

#ifdef CCADISCRET

#include "post_drt_ensight_single_field_writers.H"

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
/*----------------------------------------------------------------------*/
string AleEnsightWriter::WriteAllResults(
        PostField* field)
{
    stringstream str;
    str << EnsightWriter::WriteResult("dispnp", "displacement", field->problem()->num_dim());
    return str.str();
}

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

#endif
