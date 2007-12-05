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
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management.H"

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
//    EnsightWriter::WriteResult("velnp", "velocity", field->problem()->num_dim());
//    EnsightWriter::WriteResult("velnp", "pressure", 1, field->problem()->num_dim());
//    EnsightWriter::WriteResult("residual", "residual", field->problem()->num_dim());

    cout << "now to xfem solutions" << endl;
    set<XFEM::PHYSICS::Field> velocity_fieldset;
    velocity_fieldset.insert(XFEM::PHYSICS::Velx);
    velocity_fieldset.insert(XFEM::PHYSICS::Vely);
    velocity_fieldset.insert(XFEM::PHYSICS::Velz);
    set<XFEM::PHYSICS::Field> pressure_fieldset;
    pressure_fieldset.insert(XFEM::PHYSICS::Pres);
    
    XFluidEnsightWriter::WriteResult("velnp", "velocity_physical", velocity_fieldset);
    XFluidEnsightWriter::WriteResult("residual", "residual_physical", velocity_fieldset);
    XFluidEnsightWriter::WriteResult("velnp", "pressure_physical", pressure_fieldset);
    
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
void XFluidEnsightWriter::WriteGeoFile(
        const string& geofilename) const
{
    // open file
    ofstream geofile;
    if (myrank_ == 0)
    {
      geofile.open(geofilename.c_str());
      if (!geofile)
          dserror("failed to open file: %s", geofilename.c_str());
    }

    // header
    Write(geofile, "C Binary");

    // print out one timestep
    // if more are needed, this has to go into a loop
    map<string, vector<ofstream::pos_type> > resultfilepos;
    
    
    vector<ofstream::pos_type> fileposition;
    {
        WriteGeoFileOneTimeStep(geofile, resultfilepos, "geo");
    }

    // append index table
    // TODO: ens_checker complains if this is turned!!!! but I can't see, whats wrong here a.ger 11/07
    // it is also correct to ommit WriteIndexTable, however the EnsightGold Format manual says, 
    // it would improve performance to have it on...
    // WriteIndexTable(geofile, fileposition); 
    if (geofile.is_open())
      geofile.close();
    return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFileOneTimeStep(
        ofstream& file,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string name) const
{
    // initial Intersection
    RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(field_->discretization(),
                                                                  cutterfield_->discretization()));
    // apply enrichments
    RCP<XFEM::DofManager> initialdofmanager = rcp(new XFEM::DofManager(ih));
    
    vector<ofstream::pos_type>& filepos = resultfilepos[name];
    Write(file, "BEGIN TIME STEP");
    filepos.push_back(file.tellp());
    
    Write(file, field_->name() + " geometry");
    Write(file, "Comment");
    //Write(file,"node id given");
    Write(file, "node id assign");
    Write(file, "element id off");

    WriteGeoFilePart(file, resultfilepos, name, ih);
    
    Write(file, "END TIME STEP");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFilePart(
        ofstream& file,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string name,
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    // part + partnumber + comment
    // careful! field_->field_pos() returns the position of the ccarat
    // field, ignoring the discretizations. So if there are many
    // discretizations in one field, we have to do something different...
    Write(file, "part");
    Write(file, field_->field_pos()+1);
    Write(file, field_->name() + " field");
    
    Write(file, "coordinates");
    Write(file, NumNodesPerField(ih));
    
    cout << "writing " << NumNodesPerField(ih) << " nodes" << endl;
    
    // write the grid information
    WriteCoordinates(file, field_->discretization(), ih);
    WriteCells(file, field_->discretization(), ih);
}

int XFluidEnsightWriter::NumNodesPerField(
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    RCP<DRT::Discretization> dis = field_->discretization();
    
    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "xfem filter cannot be run in parallel");
    
    // loop all available elements
    int counter = 0;
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
        const int elegid = elementmap->GID(iele);
        const XFEM::DomainIntCells domainintcells = ih->domainIntCells(elegid, dis->gElement(elegid)->Shape());
        for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
            switch (cell->Shape())
            {
            case DRT::Element::hex20:
            case DRT::Element::hex8:
            case DRT::Element::quad4:
            case DRT::Element::quad8:
            case DRT::Element::tet4:
            case DRT::Element::tet10:
            case DRT::Element::tri3:
            case DRT::Element::wedge6:
            case DRT::Element::wedge15:
            case DRT::Element::pyramid5:
            {
                // standard case with direct support
                const int numnp = cell->NumNode();
                for (int inode=0; inode<numnp; ++inode)
                {
                    counter++;
                }
                break;
            }
            case DRT::Element::hex27:
            {
                // write subelements
                for (int isubele=0; isubele<8; ++isubele)
                    for (int isubnode=0; isubnode<8; ++isubnode)
                    {
                        counter++;
                    }
                break;
            }
            default:
                dserror("don't know, how to write this element type as a Cell");
            };
        };
    };
    return counter;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteCoordinates(
        ofstream& geofile,                          ///< filestream for the geometry
        const RefCountPtr<DRT::Discretization> dis, ///< discretization where the nodal positions are take from
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    const int nsd = 3;
    
    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "xfem filter cannot be run in parallel");

    // get the number of elements for each distype
    const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis, ih);

    for (int isd = 0; isd < nsd; ++isd)
    {
        // for each found distype write block of the same typed elements
        for (NumElePerDisType::const_iterator iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
        {
            const DRT::Element::DiscretizationType distypeiter = iter->first;

            for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
            {
                const int elegid = elementmap->GID(iele);
                DRT::Element* const actele = dis->gElement(elegid);
                const XFEM::DomainIntCells domainintcells = ih->domainIntCells(elegid, actele->Shape());
                for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
                {
                    if (cell->Shape() == distypeiter)
                    {
                        const vector<vector<double> > xarray = cell->GetPhysicalCoord(*actele);
                        for (int inen = 0; inen < cell->NumNode(); ++inen)
                        {
                            Write(geofile, static_cast<float>(xarray[inen][isd]));                    
                        }
                    }
                }
            }
        }
    }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteCells(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis,
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "filter cannot be run in parallel");

    // get the number of elements for each distype
    const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis, ih);

    // for each found distype write block of the same typed elements
    int counter = 0;
    for (NumElePerDisType::const_iterator iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
    {
        const DRT::Element::DiscretizationType distypeiter = iter->first;
        const int ne = GetNumEleOutput(distypeiter, iter->second);
        const string ensightCellType = GetEnsightString(distypeiter);

        map<DRT::Element::DiscretizationType, string>::const_iterator entry = distype2ensightstring_.find(distypeiter);
        if (entry == distype2ensightstring_.end())
            dserror("no entry in distype2ensightstring_ found");
        const string realcellshape = entry->second;
        
        cout << "writing "<< iter->second<< " "<< realcellshape << " elements"
                << " ("<< ne << " cells) per distype."<< " ensight output celltype: "
                << ensightCellType<< endl;

        Write(geofile, ensightCellType);
        Write(geofile, ne);

        // loop all available elements
        for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
        {
            const int elegid = elementmap->GID(iele);
            const XFEM::DomainIntCells domainintcells = ih->domainIntCells(elegid, dis->gElement(elegid)->Shape());
            for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
            {
                if (cell->Shape() == distypeiter)
                {
                    switch (cell->Shape())
                    {
                    case DRT::Element::hex20:
                    case DRT::Element::hex8:
                    case DRT::Element::quad4:
                    case DRT::Element::quad8:
                    case DRT::Element::tet4:
                    case DRT::Element::tet10:
                    case DRT::Element::tri3:
                    case DRT::Element::wedge6:
                    case DRT::Element::wedge15:
                    case DRT::Element::pyramid5:
                    {
                        // standard case with direct support
                        const int numnp = cell->NumNode();
                        for (int inode=0; inode<numnp; ++inode)
                        {
                            Write(geofile, counter+1);
                            counter++;
                        }
                        break;
                    }
                    case DRT::Element::hex27:
                    {
                        // write subelements
                        for (int isubele=0; isubele<8; ++isubele)
                            for (int isubnode=0; isubnode<8; ++isubnode)
                            {
                                Write(geofile, counter+1);
                                counter++;
                            }
                        break;
                    }
                    default:
                        dserror("don't know, how to write this element type as a Cell");
                    };
                };
            };
        };
    };
}


/*!
 * \brief parse all elements and get the number of elements for each distype
 */
NumElePerDisType XFluidEnsightWriter::GetNumElePerDisType(
        const RefCountPtr<DRT::Discretization> dis,
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    
    NumElePerDisType numElePerDisType;
    const Epetra_Map* elementmap = dis->ElementRowMap();
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
        const int elegid = elementmap->GID(iele);
        const XFEM::DomainIntCells domainintcells = ih->domainIntCells(elegid, dis->gElement(elegid)->Shape());
        for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
        {
            // update counter for current distype
            numElePerDisType[cell->Shape()]++;
        }
    }
    return numElePerDisType;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteResult(
        const string groupname,
        const string name,
        const set<XFEM::PHYSICS::Field> fieldset
        )
{
    PostResult result = PostResult(field_);
    result.next_result();
    if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
        return;

    
    // initial Intersection
    RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(field_->discretization(),
                                                                  cutterfield_->discretization()));
    // apply enrichments
    RCP<XFEM::DofManager> dofman = rcp(new XFEM::DofManager(ih));
    
    
    // for file continuation
    bool multiple_files = false;
    const int numdf = fieldset.size();
    const int numnp = NumNodesPerField(ih);
    const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);

    const string filename = filename_ + "_"+ field_->name() + "."+ name;
    ofstream file(filename.c_str());

    map<string, vector<ofstream::pos_type> > resultfilepos;
    WriteResultStep(file, result, resultfilepos, groupname, name, fieldset, ih, dofman);
    while (result.next_result())
    {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
        {
            dserror("nope");
            FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
        }
        WriteResultStep(file, result, resultfilepos, groupname, name, fieldset, ih, dofman);
    }

    // append index table
    WriteIndexTable(file, resultfilepos[name]);
    resultfilepos[name].clear();
    
    // store information for later case file creation
    filesetmap_[name].push_back(file.tellp()/stepsize);
    variablenumdfmap_[name] = numdf;
    variablefilenamemap_[name] = filename;
}

vector<double> computeScalarCellNodeValues(
        DRT::Element&  ele,
        const XFEM::ElementDofManager& dofman,
        const XFEM::DomainIntCell& cell,
        const XFEM::PHYSICS::Field field,
        const vector<double> elementvalues
        )
{
    vector<double> cellvalues;
    
    const int nen_cell = DRT::Utils::getNumberOfElementNodes(cell.Shape());
    const int numparam  = dofman.NumDofPerField(field);
    
    const int maxnod = 27;
    const int nsd = 3;
    
    blitz::Range _  = blitz::Range::all();
    
    DRT::Node** const nodes = ele.Nodes();
    blitz::Array<double,2> xyze(nsd,maxnod,blitz::ColumnMajorArray<2>());
    for (int inode=0; inode<ele.NumNode(); inode++)
    {
      const double* x = nodes[inode]->X();
      xyze(0,inode) = x[0];
      xyze(1,inode) = x[1];
      xyze(2,inode) = x[2];
    }
    
    // cell corner nodes
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        // shape functions
        blitz::Array<double,1> funct(maxnod);
        DRT::Utils::shape_function_3D(
                funct,
                cell.GetDomainCoord()[inen][0],
                cell.GetDomainCoord()[inen][1],
                cell.GetDomainCoord()[inen][2],
                ele.Shape());
    
        vector<double> cellnodeposvector = cell.GetPhysicalCoord(ele)[inen];
        blitz::Array<double,1> cellnodepos(3);
        for (int isd = 0; isd < nsd; ++isd) {
            cellnodepos(isd) = cellnodeposvector[isd];
        }
        
        // simplest case: jump or void enrichment
        // -> no chain rule, since enrichment function derivative is zero
        // TODO: if enrichment function has derivatives not equal to zero, we need the chain rule here
        
        blitz::Array<double,1> enr_funct(numparam);
        int dofcounter = 0;
        for (int inode=0; inode<ele.NumNode(); inode++)
        {
          const int gid = nodes[inode]->Id();
          const blitz::Array<double,1> nodalpos(xyze(_,inode));

          const std::set<XFEM::FieldEnr>  enrfieldset = dofman.FieldEnrSetPerNode(gid);
          for (std::set<XFEM::FieldEnr>::const_iterator enrfield = enrfieldset.begin(); enrfield != enrfieldset.end(); ++enrfield)
          {
              if (enrfield->getField() == XFEM::PHYSICS::Velx)
              {
                  const XFEM::Enrichment enr = enrfield->getEnrichment();
                  const double enrval = dofman.enrValue(enr,cellnodepos,nodalpos);
                  enr_funct(dofcounter) = funct(inode) * enrval;
                  dofcounter += 1;
              }
          }
        }
        
        // interpolate value
        double x = 0.0;
        for (int iparam = 0; iparam < numparam; ++iparam)
        {
            x += elementvalues[iparam] * enr_funct(iparam);
        }
        // store position
        cellvalues.push_back(x);
    }
    return cellvalues;
}




/*!
 \brief Write nodal values for one timestep
 */
void XFluidEnsightWriter::WriteResultStep(
        ofstream& file,
        PostResult& result,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string groupname,
        const string name,
        const set<XFEM::PHYSICS::Field> fieldset,
        const RCP<XFEM::InterfaceHandle> ih,
        const RCP<XFEM::DofManager> dofman
        ) const
{
    vector<ofstream::pos_type>& filepos = resultfilepos[name];
    Write(file, "BEGIN TIME STEP");
    filepos.push_back(file.tellp());

    Write(file, "description");
    Write(file, "part");
    Write(file, field_->field_pos()+1);
    Write(file, "coordinates");

    const RefCountPtr<DRT::Discretization> dis = field_->discretization();
    //const Epetra_Map* nodemap = dis->NodeRowMap();
    const RefCountPtr<Epetra_Vector> data = result.read_result(groupname);
    const Epetra_BlockMap& datamap = data->Map();
    
    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "xfem filter cannot be run in parallel");
    
    const int numdf = fieldset.size();
    
    bool propernumdf;
    if (numdf == 1 or numdf == 3 or numdf == 6) {
        propernumdf = true;
    } else {
        propernumdf = false;
    }
    if (not propernumdf)
    {
        dserror("number of output dofs is not 1, 3 or 6");
    }
    
    
    
    // get the number of elements for each distype
    const NumElePerDisType numElePerDisType = GetNumElePerDisType(dis, ih);
    
    //
    const XFEM::Enrichment enr_std(0, XFEM::Enrichment::typeStandard);
    
    for (set<XFEM::PHYSICS::Field>::const_iterator fielditer=fieldset.begin(); fielditer!=fieldset.end(); ++fielditer)
    {
        const XFEM::PHYSICS::Field field = *fielditer;
        
        // for each found distype write block of the same typed elements
        int counter = 0;
        for (NumElePerDisType::const_iterator iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
        {
            const DRT::Element::DiscretizationType distypeiter = iter->first;
            for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
            {
                const int elegid = elementmap->GID(iele);
                DRT::Element* const actele = dis->gElement(elegid);
                
                // create local copy of information about dofs
                XFEM::ElementDofManager eledofman = dofman->constructElementDofManager(*actele);
                
                vector<int> lm;
                vector<int> lmowner;
                actele->LocationVector(*(field_->discretization()),lm,lmowner);
                
                // extract local values from the global vector
                vector<double> myvelnp(lm.size());
                DRT::Utils::ExtractMyValues(*data,myvelnp,lm);
                
                const int numparam = eledofman.NumDofPerField(field);
                const vector<int> dof = eledofman.LocalDofPosPerField(field);
                
                vector<double> elementvalues(numparam);
                for (int iparam=0; iparam<numparam; ++iparam)   elementvalues[iparam] = myvelnp[dof[iparam]];
                
                const XFEM::DomainIntCells domainintcells = ih->domainIntCells(elegid, actele->Shape());
                for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
                {
                    if (cell->Shape() == distypeiter)
                    {
                        const vector<double> cellvalues = computeScalarCellNodeValues(*actele, eledofman, *cell, field, elementvalues);
                        for (int inode = 0; inode < cell->NumNode(); ++inode)
                        {
                            Write(file, static_cast<float>(cellvalues[inode]));
                            counter++;
                        }
                    }
                }
            }
        }
        cout << "number of entries per field " << counter << endl;
    }

    Write(file, "END TIME STEP");
    
    return;
}




#endif
