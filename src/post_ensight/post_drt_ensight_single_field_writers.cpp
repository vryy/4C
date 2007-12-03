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
    
    XFluidEnsightWriter::WriteResult("velnp", "velocity(physical)", field->problem()->num_dim());
    XFluidEnsightWriter::WriteResult("velnp", "pressure(physical)", 1, field->problem()->num_dim());
    XFluidEnsightWriter::WriteResult("residual", "residual(physical)", field->problem()->num_dim());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteFiles()
{
#ifndef PARALLEL
    if (myrank_ > 0) dserror("have serial filter version, but myrank_ > 0");
#endif

    PostResult result = PostResult(field_);

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
    // write solution fields files
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
void XFluidEnsightWriter::WriteGeoFileOneTimeStep(
        ofstream& file,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string name) const
{
    vector<ofstream::pos_type>& filepos = resultfilepos[name];
    Write(file, "BEGIN TIME STEP");
    filepos.push_back(file.tellp());
    
    Write(file, field_->name() + " geometry");
    Write(file, "Comment");
    //Write(file,"node id given");
    Write(file, "node id assign");
    Write(file, "element id off");

    WriteGeoFilePart(file, resultfilepos, name);
    
    Write(file, "END TIME STEP");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFilePart(
        ofstream& file,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string name) const
{
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteCoordinates(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis) const
{
    const Epetra_Map* nodemap = dis->NodeRowMap();
    dsassert(nodemap->NumMyElements() == nodemap->NumGlobalElements(),
            "filter cannot be run in parallel");

    const int numnp = nodemap->NumMyElements();

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
void XFluidEnsightWriter::WriteCoordinatesIntCell(
        ofstream& geofile,                          ///< filestream for the geometry
        const RefCountPtr<DRT::Discretization> dis, ///< discretization where the nodal positions are take from
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "xfem filter cannot be run in parallel");

    const set<int> elegidset = ih->getSetOfIntersectedElementIds();
    
    for (set<int>::const_iterator elegid = elegidset.begin(); elegid != elegidset.end(); ++elegid)
    {
        DRT::Element* const actele = dis->gElement(*elegid);
        //DRT::Element* actele = dis->Element(elegid); 
        XFEM::DomainIntCells domainintcells = ih->domainIntCells(*elegid, actele->Shape());
    }

//    const int NSD = 3; ///< number of space dimensions
//    // ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n
//    for (int isd=0; isd<NSD; ++isd)
//    {
//        for (int inode=0; inode<numnp; inode++)
//        {
//            const int gid = nodemap->GID(inode);
//            const DRT::Node* actnode = dis->gNode(gid);
//            Write(geofile, static_cast<float>(actnode->X()[isd]));
//        }
//    }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteCells(
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteResult(
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

#endif
