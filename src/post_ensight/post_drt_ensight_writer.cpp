/*!
 \file post_drt_ensight_writer.cpp

 \brief Ensight filter basis class

 <pre>
 Maintainer: Ulrich Kuettler
 kuettler@lnm.mw.tum.de
 http://www.lnm.mw.tum.de/Members/kuettler
 089 - 289-15238
 </pre>

 */

#ifdef CCADISCRET

#include "post_drt_ensight_writer.H"
#include <string>

using namespace std;

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EnsightWriter::EnsightWriter(
        PostField* field,
        const string filename) :
    field_(field),
    filename_(filename),
    myrank_(((field->problem())->comm())->MyPID())
{
    // map between distype in BACI and Ensight string
    //  note: these would be the direct mappings
    //  look in the actual writer, whether we output hex instead of hex27, e.g.
    distype2ensightstring_.clear();
    distype2ensightstring_[DRT::Element::hex8] = "hexa8";
    distype2ensightstring_[DRT::Element::hex20] = "hexa20";
    distype2ensightstring_[DRT::Element::hex27] = "hexa20"; //on purpose, ensight does not know hex27
    distype2ensightstring_[DRT::Element::tet4] = "tetra4";
    distype2ensightstring_[DRT::Element::tet10] = "tetra10";
    distype2ensightstring_[DRT::Element::quad4] = "quad4";
    distype2ensightstring_[DRT::Element::quad8] = "quad8";
    distype2ensightstring_[DRT::Element::quad9] = "quad9";
    distype2ensightstring_[DRT::Element::tri3] = "tria3";
    distype2ensightstring_[DRT::Element::tri6] = "tria6";
    distype2ensightstring_[DRT::Element::wedge6] = "penta6";
    distype2ensightstring_[DRT::Element::wedge15]= "penta15";
    distype2ensightstring_[DRT::Element::pyramid5]= "pyramid5";
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteFiles()
{
    PostResult result = PostResult(field_);
    
#ifndef PARALLEL
    if (myrank_>0) dserror("have serial filter version, but myrank > 0");
#endif
    
    // timesteps when the solution is written
    const vector<double> soltime = result.get_result_times(field_->name());

    //
    // write geometry file
    //
    const string geofilename = filename_ + "_"+ field_->name() + ".geo";
    WriteGeoFile(geofilename);
    vector<int> filesteps;
    filesteps.push_back(1);
    filesetmap_["geo"] = filesteps;
    vector<int> timesteps;
    timesteps.push_back(1);
    timesetmap_["geo"] = timesteps;
    const int geotimeset = 1;
    const int geofileset = 1;
    // at the moment, we can only print out the first step -> to be changed
    vector<double> geotime; // timesteps when the geometry is written
    geotime.push_back(soltime[0]);
    
    
    //
    // solution fields files
    //
    const int soltimeset = 2;
    WriteAllResults(field_);
    
    
    
    
    
    //
    // now write the case file
    //
    const string casefilename = filename_ + "_"+ field_->name() + ".case";
    ofstream casefile;
    casefile.open(casefilename.c_str());
    casefile << "# created using post_drt_ensight\n"<< "FORMAT\n\n"<< "type:\tensight gold\n";
    
    casefile << "\nGEOMETRY\n\n";
    casefile << "model:\t"<<geotimeset<<"\t"<<geofileset<<"\t"<< geofilename<< "\n";
        
    casefile << "\nVARIABLE\n\n";
    casefile << GetVariableSection(filesetmap_, variablenumdfmap_, variablefilenamemap_);
    
    casefile << "\nTIME\n\n";
    casefile << GetTimeSectionString(geotimeset, geotime);
    casefile << GetTimeSectionString(soltimeset, soltime);

    casefile << "\nFILE\n\n";
    casefile << GetFileSectionStringFromFilesets(filesetmap_);

    casefile.close();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteGeoFile(
        const string& geofilename) const
{
    // open file
    ofstream geofile;
    geofile.open(geofilename.c_str());
    if (!geofile)
        dserror("failed to open file: %s", geofilename.c_str());

    // header
    Write(geofile, "C Binary");

    // print out one timestep
    // if more are needed, this has to go into a loop
    vector<ofstream::pos_type> fileposition;
    {
        Write(geofile, "BEGIN TIME STEP");
        fileposition.push_back(geofile.tellp());
        Write(geofile, field_->name() + " geometry");
        Write(geofile, "Comment");
        Write(geofile,"node id given");
        //Write(geofile, "node id assign");
        Write(geofile, "element id off");

        // part + partnumber + comment
        // careful! field_->field_pos() returns the position of the ccarat
        // field, ignoring the discretizations. So if there are many
        // discretizations in one field, we have to do something different...
        Write(geofile, "part");
        Write(geofile, field_->field_pos()+1);
        Write(geofile, field_->name() + " field");

        Write(geofile, "coordinates");
        Write(geofile, field_->num_nodes());

        // write the grid information
        WriteCoordinates(geofile, field_->discretization());
        WriteCells(geofile, field_->discretization());

        Write(geofile, "END TIME STEP");
    }

    // append index table
    // TODO: ens_checker complains if this is turned!!!! but I can't see, whats wrong here a.ger 11/07
    // it is also correct to ommit WriteIndexTable, however the EnsightGold Format manual says, 
    // it would improve performance to have it on...
    // WriteIndexTable(geofile, fileposition); 

    geofile.close();
}


void EnsightWriter::WriteGeoFileOneTimeStep(
        ofstream& file,
        PostResult& result,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string groupname,
        const string name,
        const int from) const
{
    vector<ofstream::pos_type>& filepos = resultfilepos[name];
    Write(file, "BEGIN TIME STEP");
    filepos.push_back(file.tellp());
    
    Write(file, field_->name() + " geometry");
    Write(file, "Comment");
    Write(file,"node id given");
    //Write(geofile, "node id assign");
    Write(file, "element id off");

    // part + partnumber + comment
    // careful! field_->field_pos() returns the position of the ccarat
    // field, ignoring the discretizations. So if there are many
    // discretizations in one field, we have to do something different...
    Write(file, "part");
    Write(file, field_->field_pos()+1);
    Write(file, field_->name() + " field");

    Write(file, "coordinates");
    Write(file, field_->num_nodes());

    // write the grid information
    WriteCoordinates(file, field_->discretization());
    WriteCells(file, field_->discretization());

    Write(file, "END TIME STEP");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCoordinates(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis) const
{
    const Epetra_Map* nodemap = dis->NodeRowMap();
    dsassert(nodemap->NumMyElements() == nodemap->NumGlobalElements(),
            "filter cannot be run in parallel");

    const int numnp = nodemap->NumMyElements();
    // write node ids
    for (int inode=0; inode<numnp; ++inode)
    {
        Write(geofile,nodemap->GID(inode)+1);
    }

    const int NSD = 3; ///< number of space dimensions
    // ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n
    for (int isd=0; isd<NSD; ++isd)
    {
        for (int inode=0; inode<numnp; inode++)
        {
            const int gid = nodemap->GID(inode);
            const DRT::Node* actnode = dis->gNode(gid);
            Write(geofile, static_cast<float>(actnode->X()[isd]));
        }
    }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteCells(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis) const
{
    const Epetra_Map* nodemap = dis->NodeRowMap();
    dsassert(nodemap->NumMyElements() == nodemap->NumGlobalElements(),
            "filter cannot be run in parallel");

    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "filter cannot be run in parallel");

    // get the number of elements for each distype
    const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis);

    // for each found distype write block of the same typed elements
    NumElePerDisType::const_iterator iter;
    for (iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
    {
        const DRT::Element::DiscretizationType distypeiter = iter->first;
        const int ne = GetNumEleOutput(distypeiter, iter->second);
        const string ensightString = GetEnsightString(distypeiter);

        map<DRT::Element::DiscretizationType, string>::const_iterator entry = distype2ensightstring_.find(distypeiter);
        if (entry == distype2ensightstring_.end())
            dserror("no entry in distype2ensightstring_ found");
        const string realcellshape = entry->second;
        
        cout << "writing "<< iter->second<< " "<< realcellshape << " elements"
                << " ("<< ne << " cells) per distype."<< " ensight output celltype: "
                << ensightString<< endl;

        Write(geofile, ensightString);
        Write(geofile, ne);

        // loop all available elements
        for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
        {
            DRT::Element* const actele = dis->gElement(elementmap->GID(iele));
            if (actele->Shape() == distypeiter)
            {
                DRT::Node** const nodes = actele->Nodes();
                switch (actele->Shape())
                {
                case DRT::Element::hex20:
                case DRT::Element::hex8:
                case DRT::Element::quad4:
                case DRT::Element::quad8:
                case DRT::Element::tet4:
                case DRT::Element::tri3:
                case DRT::Element::wedge6:
                case DRT::Element::wedge15:
                case DRT::Element::pyramid5:
                {
                    // standard case with direct support
                    const int numnp = actele->NumNode();
                    for (int inode=0; inode<numnp; ++inode)
                        Write(geofile, nodemap->LID(nodes[inode]->Id())+1);
                    break;
                }
                case DRT::Element::hex27:
                {
                    // write subelements
                    for (int isubele=0; isubele<8; ++isubele)
                        for (int isubnode=0; isubnode<8; ++isubnode)
                            Write(geofile, nodemap->LID(nodes[subhexmap[isubele][isubnode]]->Id())
                                    +1);
                    break;
                }
                case DRT::Element::quad9:
                {
                    // write subelements
                    for (int isubele=0; isubele<4; ++isubele)
                        for (int isubnode=0; isubnode<4; ++isubnode)
                            Write(geofile, nodemap->LID(nodes[subquadmap[isubele][isubnode]]->Id())
                                    +1);
                    break;
                }
                case DRT::Element::tet10:
                {
                    // write just the corner nodes of the tet10
                  for (int isubnode=0; isubnode<4; ++isubnode)
                      Write(geofile, nodemap->LID(nodes[subtet10map[0][isubnode]]->Id())
                              +1);
                  break;
               }
                default:
                    dserror("don't know, how to write this element type as a Cell");
                };
            };
        };
    };
}

/*!
 * \brief parse all elements and get the number of elements for each distype
 */
NumElePerDisType EnsightWriter::GetNumElePerDisType(
        const RefCountPtr<DRT::Discretization> dis) const
{
    const Epetra_Map* elementmap = dis->ElementRowMap();

    NumElePerDisType numElePerDisType;
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
        DRT::Element* actele = dis->gElement(elementmap->GID(iele));
        const DRT::Element::DiscretizationType distype = actele->Shape();
        // update counter for current distype
        numElePerDisType[distype]++;
    }
    return numElePerDisType;
}

int EnsightWriter::GetNumEleOutput(
        const DRT::Element::DiscretizationType distype,
        const int numele) const
{

    int numeleout = 0;
    switch (distype)
    {
    case DRT::Element::hex27:
        numeleout = 8*numele;
        break;
    case DRT::Element::quad9:
        numeleout = 4*numele;
        break;
    default:
        numeleout = numele;
    }
    return numeleout;
}

string EnsightWriter::GetEnsightString(
        const DRT::Element::DiscretizationType distype) const
{
    map<DRT::Element::DiscretizationType, string>::const_iterator entry;
    switch (distype)
    {
    case DRT::Element::hex27:
        entry = distype2ensightstring_.find(DRT::Element::hex8);
        break;
    case DRT::Element::quad9:
        entry = distype2ensightstring_.find(DRT::Element::quad4);
        break;
    case DRT::Element::tet10:
        entry = distype2ensightstring_.find(DRT::Element::tet4);
        break;
    default:
        entry = distype2ensightstring_.find(distype);
    }
    if (entry == distype2ensightstring_.end())
        dserror("no entry in distype2ensightstring_ found");
    return entry->second;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteResult(
        const string groupname,
        const string name,
        const int numdf,
        const int from)
{
    PostResult result = PostResult(field_);
    result.next_result();
    if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
        return;

    // new for file continuation
    bool multiple_files = false;
    const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
    const int numnp = nodemap->NumMyElements();
    const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);

    const string filename = filename_ + "_"+ field_->name() + "."+ name;
    ofstream file(filename.c_str());

    map<string, vector<ofstream::pos_type> > resultfilepos;
    WriteResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    while (result.next_result())
    {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
        {
            FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
        }
        WriteResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    }

    // append index table
    WriteIndexTable(file, resultfilepos[name]);
    resultfilepos[name].clear();
    
    // store information for later case file creation
    filesetmap_[name].push_back(file.tellp()/stepsize);
    variablenumdfmap_[name] = numdf;
    variablefilenamemap_[name] = filename;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::FileSwitcher(
        ofstream& file,
        bool& multiple_files,
        map<string,vector<int> >& filesetmap,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const int stepsize,
        const string name,
        const string filename) const
{
    ostringstream newfilename;
    
    if (multiple_files == false)
    {
        multiple_files = true;

        vector<int> numsteps;
        numsteps.push_back(file.tellp()/stepsize);
        filesetmap[name] = numsteps;
        
        // append index table
        WriteIndexTable(file, resultfilepos[name]);
        resultfilepos[name].clear();
        file.close();
        rename(filename.c_str(), (filename+"001").c_str());
        
        newfilename << filename << "002";
    }
    else
    {
        filesetmap[name].push_back(file.tellp()/stepsize);
        
        // append index table
        WriteIndexTable(file, resultfilepos[name]);
        resultfilepos[name].clear();
        file.close();

        newfilename << filename;
        newfilename.width(3);
        newfilename.fill('0');
        newfilename << filesetmap[name].size()+1;
    }
    file.open(newfilename.str().c_str());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetVariableEntryForCaseFile(
        const int numdf,
        const unsigned int fileset,
        const string name,
        const string filename) const
{
    stringstream str;
    
    const int timeset = 2;
    
    // create variable entry
    switch (numdf)
    {
    case 6:
        str << "tensor symm per node:\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
        break;
    case 3:
    case 2:
        str << "vector per node:\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
        break;
    case 1:
        str << "scalar per node:\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
        break;
    default:
        dserror("unknown number of dof per node");
    };
    return str.str();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetVariableSection(
        map<string,vector<int> > filesetmap,
        map<string,int> variablenumdfmap,
        map<string,string> variablefilenamemap
        ) const
{
    stringstream str;
    
    map<string,int>::const_iterator variable;
    
    int filesetcounter = 2; // we begin with 2, since 1 is for the geometry section
    for (variable = variablenumdfmap.begin(); variable != variablenumdfmap.end(); ++variable)
    {
        const string name = variable->first;
        const int numdf = variable->second;
        const string filename = variablefilenamemap[name];
        
        map<string,vector<int> >::const_iterator entry = filesetmap.find(name);
        if (entry == filesetmap.end()) 
            dserror("filesetmap not defined for '%s'", name.c_str());
        
        const int numsubfilesteps = entry->second.size();
        string filename_for_casefile;
        if (numsubfilesteps > 1)
        {
            filename_for_casefile = filename + "***";
        }
        else
        {
            filename_for_casefile = filename;
        }
            
        str << GetVariableEntryForCaseFile(numdf, filesetcounter, name, filename_for_casefile);
        filesetcounter++;
    }

    return str.str();
}


/*!
 \brief Write nodal values for one timestep

 Each node has to have the same number of dofs.
 */
void EnsightWriter::WriteResultStep(
        ofstream& file,
        PostResult& result,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string groupname,
        const string name,
        const int numdf,
        const int from) const
{
    vector<ofstream::pos_type>& filepos = resultfilepos[name];
    Write(file, "BEGIN TIME STEP");
    filepos.push_back(file.tellp());

    Write(file, "description");
    Write(file, "part");
    Write(file, field_->field_pos()+1);
    Write(file, "coordinates");

    const RefCountPtr<DRT::Discretization> dis = field_->discretization();
    const Epetra_Map* nodemap = dis->NodeRowMap();
    const RefCountPtr<Epetra_Vector> data = result.read_result(groupname);
    const Epetra_BlockMap& datamap = data->Map();

    const int numnp = nodemap->NumMyElements();
    for (int idf=0; idf<numdf; ++idf)
    {
        for (int inode=0; inode<numnp; inode++)
        {
            DRT::Node* n = dis->lRowNode(inode);
            const int lid = datamap.LID(dis->Dof(n, from+idf));
            if (lid > -1)
            {
                Write(file, static_cast<float>((*data)[lid]));
            }
            else
            {
                // Assume we have to write a value here.
                Write<float>(file, 0.);
            }
        }
    }

    // 2 component vectors in a 3d problem require a row of zeros.
    // do we really need this?
    if (numdf==2)
    {
        for (int inode=0; inode<numnp; inode++)
        {
            Write<float>(file, 0.);
        }
    }

    Write(file, "END TIME STEP");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteIndexTable(
        ofstream& file,
        const vector<ofstream::pos_type>& filepos) const
{
    ofstream::pos_type lastpos = file.tellp();
    const unsigned steps = filepos.size();
    Write(file, steps);
    for (unsigned i=0; i<steps; ++i)
    {
        Write<long>(file, filepos[i]);
    }
    Write(file, 0);
    Write<long>(file, lastpos);
    Write(file, "FILE_INDEX");
}

/*----------------------------------------------------------------------*/
/*!
 \brief write strings of exactly 80 chars
 */
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteString(
        ofstream& stream,
        const string str) const
{
    // we need to write 80 bytes per string
    vector<char> s(str.begin(), str.end());
    while (s.size()<80)
    {
        s.push_back('\0');
    }
    stream.write(&s[0], 80);
}

/*----------------------------------------------------------------------*/
/*!
 \brief create string for one TIME section in the case file
 */
/*----------------------------------------------------------------------*/
string EnsightWriter::GetTimeSectionString(
        const int timeset,
        const vector<double>& times) const
{
    stringstream s;
    s << "\ntime set:\t\t" << timeset << "\n"<< "number of steps:\t"<< times.size() << "\ntime values: ";
    for (unsigned i=0; i<times.size(); ++i)
    {
        s << times[i]<< " ";
        if (i%8==0&& i!=0)
        {
            s << "\n";
        }
    }
    s << "\n";
    return s.str();
}



/*----------------------------------------------------------------------*/
/*!
 \brief create string for the FILE section in the case file
 */
/*----------------------------------------------------------------------*/
string EnsightWriter::GetFileSectionStringFromFilesets(
        const map<string,vector<int> >& filesetmap) const
{
    stringstream s;
        
    map<string,vector<int> >::const_iterator fileset;
    int setnumber = 1;

    for (fileset = filesetmap.begin(); fileset != filesetmap.end(); ++fileset)
    {
        vector<int> subsets = fileset->second;
        s << "file set:\t\t"<< setnumber << "\n";
        if (subsets.size() == 1)
        {
            s << "number of steps:\t"<< subsets[0] << "\n\n";
        }
        else
        {
            for (unsigned int j = 0; j < subsets.size(); ++j)
            {
                s << "filename index:\t"<< 1+j << "\n";
                s << "number of steps:\t"<< subsets[j] << "\n";
            }
            s << "\n";
        }
        setnumber++;
    }
    return s.str();
}



#endif
