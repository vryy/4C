/*!
 \file post_drt_ensight_single_field_writers.cpp

 \brief main routine of the Ensight filter

 <pre>
  Maintainer: Axel Gerstenberger
 gerstenberger@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/gerstenberger
 089 - 289-15236
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
    str << EnsightWriter::WriteResult("velnp", "velocity", field->problem()->num_dim());
    str << EnsightWriter::WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
    str << EnsightWriter::WriteResult("residual", "residual", field->problem()->num_dim());
    
    str << XFluidEnsightWriter::WriteResult("velnp", "velocity(physical)", field->problem()->num_dim());
    str << XFluidEnsightWriter::WriteResult("velnp", "pressure(physical)", 1, field->problem()->num_dim());
    str << XFluidEnsightWriter::WriteResult("residual", "residual(physical)", field->problem()->num_dim());
    return str.str();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteFiles()
{
    // loop all results to get number of result timesteps
    PostResult result = PostResult(field_);
    vector<double> soltime; // timesteps when the solution is written
    if (result.next_result())
        soltime.push_back(result.time());
    else
        dserror("no solution found in field '%s'", field_->name().c_str());

    while (result.next_result())
        soltime.push_back(result.time());

    //
    // now do the case file
    //
    const string casefilename = filename_ + "_"+ field_->name() + ".case";
    ofstream casefile;
    casefile.open(casefilename.c_str());
    casefile << "# created using post_drt_ensight\n"<< "FORMAT\n\n"<< "type:\tensight gold\n";

    //
    // GEOMETRY section
    //
    const string geofilename = filename_ + "_"+ field_->name() + ".geo";
    casefile << "\nGEOMETRY\n\n"<< "model:\t2\t2\t"<< geofilename<< "\n";
    WriteGeoFile(geofilename);

    //
    // VARIABLE section
    //
    // whatever result we need, inside the variable entries in the case file are generated
    casefile << "\nVARIABLE\n\n";
    casefile << WriteAllResults(field_);

    //
    // TIME section
    //
    casefile << "\nTIME\n";

    // write time steps for result file (time set 1)
    casefile << "time set:\t\t1\n"<< "number of steps:\t"<< soltime.size() << "\ntime values: ";
    for (unsigned i=0; i<soltime.size(); ++i)
    {
        casefile << soltime[i]<< " ";
        if (i%8==0&& i!=0)
            casefile << "\n";
    }
    casefile << "\n\n";

    // write time steps for geometry file (time set 2)
    // at the moment, we can only print out the first step -> to be changed
    vector<double> geotime; // timesteps when the geometry is written
    geotime.push_back(soltime[0]);
    casefile << "time set:\t\t2\n"<< "number of steps:\t"<< geotime.size() << "\ntime values: ";
    for (unsigned i=0; i<geotime.size(); ++i)
    {
        casefile << geotime[i]<< " ";
        if (i%8==0&& i!=0)
            casefile << "\n";
    }
    casefile << "\n\n";

    //
    // FILE section
    //
    casefile << "FILE\n";
    casefile << "file set:\t\t1\n"<< "number of steps:\t"<< soltime.size()
            << "\n\nfile set:\t\t2\n"<< "number of steps:\t"<< geotime.size() << "\n\n";
    for (unsigned int i = 0; i < filesets_.size(); ++i)
    {
        casefile << "\nfile set:\t\t"<< 3+i << "\n";
        for (unsigned int j = 0; j < filesets_[i].size(); ++j)
        {
            casefile << "filename index:\t"<< 1+j << "\n"<< "number of steps:\t"<< filesets_[i][j]
                    << "\n";
        }
    }

    casefile.close();
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

    map<string, vector<ofstream::pos_type> > resultfilepos;
    WriteResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    while (result.next_result())
    {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT)
        {
            FileSwitcher(file, fileset, filesets_, resultfilepos, stepsize, name, filename);
        }
        WriteResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    }

    // append index table
    WriteIndexTable(file, resultfilepos[name]);
    resultfilepos[name].clear();

    if (fileset != 1)
    {
        filesets_[fileset-3].push_back(file.tellp()/stepsize);
        filename += "***";
    }

    return GetVariableEntryForCaseFile(numdf, fileset, name, filename);
}

#endif
