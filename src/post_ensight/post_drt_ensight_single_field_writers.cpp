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

#include "post_drt_ensight_writer.H"
#include "post_drt_ensight_single_field_writers.H"

using namespace std;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteAllResults(
        PostField* field)
{
    EnsightWriter::WriteResult("displacement", "displacement", field->problem()->num_dim());
    EnsightWriter::WriteResult("velocity", "velocity", field->problem()->num_dim());
    EnsightWriter::WriteResult("acceleration", "acceleration", field->problem()->num_dim());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidEnsightWriter::WriteAllResults(
        PostField* field)
{
    EnsightWriter::WriteResult("velnp", "velocity", field->problem()->num_dim());
    EnsightWriter::WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
    EnsightWriter::WriteResult("residual", "residual", field->problem()->num_dim());
    EnsightWriter::WriteResult("dispnp", "displacement", field->problem()->num_dim());
    EnsightWriter::WriteResult("traction", "traction", field->problem()->num_dim());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleEnsightWriter::WriteAllResults(
        PostField* field)
{
    EnsightWriter::WriteResult("dispnp", "displacement", field->problem()->num_dim());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteAllResults(
        PostField* field)
{
    stringstream str;
    EnsightWriter::WriteResult("velnp", "velocity", field->problem()->num_dim());
    EnsightWriter::WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
    EnsightWriter::WriteResult("residual", "residual", field->problem()->num_dim());
    
    str << XFluidEnsightWriter::WriteResult("velnp", "velocity(physical)", field->problem()->num_dim());
    str << XFluidEnsightWriter::WriteResult("velnp", "pressure(physical)", 1, field->problem()->num_dim());
    str << XFluidEnsightWriter::WriteResult("residual", "residual(physical)", field->problem()->num_dim());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteFiles()
{
    
    PostResult result = PostResult(field_);
    
    // timesteps when the solution is written
    const vector<double> soltime = result.get_result_times(field_->name());

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
    const int geotimeset = 1;
    const int geofileset = 1;
    const string geofilename = filename_ + "_"+ field_->name() + ".geo";
    casefile << "\nGEOMETRY\n\n"<< "model:\t"<<geotimeset<<"\t"<<geofileset<<"\t"<< geofilename<< "\n";
    WriteGeoFile(geofilename);

    //
    // VARIABLE section
    //
    // whatever result we need, inside the variable entries in the case file are generated
    const int soltimeset = 2;
    const int solfileset = 2; /// default, if not split to multiple files
    casefile << "\nVARIABLE\n\n";
    WriteAllResults(field_);

    //
    // TIME section
    //
    casefile << "\nTIME\n";

    // write time steps for geometry file
    // at the moment, we can only print out the first step -> to be changed
    vector<double> geotime; // timesteps when the geometry is written
    geotime.push_back(soltime[0]);
    casefile << GetTimeSectionString(geotimeset, geotime);
    
    // write time steps for result file
    casefile << GetTimeSectionString(soltimeset, soltime);

    //
    // FILE section
    //
    casefile << "FILE\n";
    casefile << "file set:\t\t"<<geofileset<<"\n"<< "number of steps:\t"<< geotime.size() << "\n\n";
    casefile << "file set:\t\t"<<solfileset<<"\n"<< "number of steps:\t"<< soltime.size() << "\n\n";
    //casefile << GetFileSectionStringFromFilesets(solfilesets_);

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
//    PostResult result = PostResult(field_);
//    result.next_result();
//    if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
//        return "";
//
//    // new for file continuation
//    bool multiple_files = false;
//    const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
//    const int numnp = nodemap->NumMyElements();
//    const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);
//
//    const string filename = filename_ + "_"+ field_->name() + "."+ name;
//    ofstream file(filename.c_str());
//
//    map<string, vector<ofstream::pos_type> > resultfilepos;
//    WriteResultStep(file, result, resultfilepos, groupname, name, numdf, from);
//    while (result.next_result())
//    {
//        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
//        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
//        {
//            //FileSwitcher(file, multiple_files, solfilesets_, resultfilepos, stepsize, name, filename);
//        }
//        WriteResultStep(file, result, resultfilepos, groupname, name, numdf, from);
//    }
//
//    // append index table
//    WriteIndexTable(file, resultfilepos[name]);
//    resultfilepos[name].clear();
//
//    string filename_for_casefile;
//    if (multiple_files)
//    {
//        const int last_fileset = solfilesets_.size()-1;
//        solfilesets_[last_fileset].push_back(file.tellp()/stepsize);
//        filename_for_casefile = filename + "***";
//    }
//    else
//    {
//        filename_for_casefile = filename;
//    }
//    const int file_set = solfilesets_.size();
//    return GetVariableEntryForCaseFile(numdf, file_set, name, filename_for_casefile);
}

#endif
