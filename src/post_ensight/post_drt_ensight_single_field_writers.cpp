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
#include <string>
#include "post_drt_ensight_single_field_writers.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/integrationcell.H"
#include "../drt_xfem/physics.H"
#include "../drt_xfem/dof_management.H"
#include "../drt_xfem/xfem.H"

using namespace std;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void StructureEnsightWriter::WriteAllResults(
        PostField* field)
{
    EnsightWriter::WriteResult("displacement", "displacement", nodebased, field->problem()->num_dim());
    EnsightWriter::WriteResult("velocity", "velocity", nodebased, field->problem()->num_dim());
    EnsightWriter::WriteResult("acceleration", "acceleration", nodebased, field->problem()->num_dim());
    WriteElementResults(field);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FluidEnsightWriter::WriteAllResults(
        PostField* field)
{
    EnsightWriter::WriteResult("velnp", "velocity", nodebased, field->problem()->num_dim());
    EnsightWriter::WriteResult("velnp", "pressure", nodebased, 1, field->problem()->num_dim());
    EnsightWriter::WriteResult("residual", "residual", nodebased, field->problem()->num_dim());
    EnsightWriter::WriteResult("dispnp", "displacement", nodebased, field->problem()->num_dim());
    EnsightWriter::WriteResult("traction", "traction", nodebased, field->problem()->num_dim());
    WriteElementResults(field);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AleEnsightWriter::WriteAllResults(
        PostField* field)
{
    EnsightWriter::WriteResult("dispnp", "displacement", nodebased, field->problem()->num_dim());
    WriteElementResults(field);
}


/*----------------------------------------------------------------------*
 |                                                             gjb 12/07  |
\*----------------------------------------------------------------------*/
void ConDifEnsightWriter::WriteAllResults(
        PostField* field)
{
    //phinp is a scalar result field with ONE dof per node.
    //Therefore it is NOT possible to hand over field->problem()->num_dim()
    // (equals 2 or 3) as a number of dofs
    EnsightWriter::WriteResult("phinp", "phi", nodebased, 1);
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
    //XFluidEnsightWriter::WriteResult("residual", "residual_physical", velocity_fieldset);
    XFluidEnsightWriter::WriteResult("velnp", "pressure_physical", pressure_fieldset);

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteFiles()
{
#ifndef PARALLEL
    if (myrank_ > 0) dserror("have serial filter version, but myrank_ = %d",myrank_);
#endif

    PostResult result = PostResult(field_);

    // timesteps when the solution is written
    const vector<double> soltime = result.get_result_times(field_->name());


    ///////////////////////////////////
    //  write geometry file          //
    ///////////////////////////////////
    const string geofilename = filename_ + "_"+ field_->name() + ".geo";
    WriteGeoFile(geofilename);
    vector<int> filesteps;
    filesteps.push_back(1);
    filesetmap_["geo"] = filesteps;
    vector<int> timesteps;
    timesteps.push_back(1);
    timesetmap_["geo"] = timesteps;
    const int geotimeset = 1;
    // at the moment, we can only print out the first step -> to be changed
    vector<double> geotime; // timesteps when the geometry is written
    geotime.push_back(soltime[0]);


    ///////////////////////////////////
    //  write solution fields files  //
    ///////////////////////////////////
    const int soltimeset = 2;
    WriteAllResults(field_);

    int counttime = 0;
    for (map<string,vector<int> >::const_iterator entry = timesetmap_.begin(); entry != timesetmap_.end(); ++entry) {
        counttime++;
        string key = entry->first;
        timesetnumbermap_[key] = counttime;
    }
    int countfile = 0;
    for (map<string,vector<int> >::const_iterator entry = filesetmap_.begin(); entry != filesetmap_.end(); ++entry) {
        countfile++;
        string key = entry->first;
        filesetnumbermap_[key] = countfile;
    }



    ///////////////////////////////////
    //  now write the case file      //
    ///////////////////////////////////
    if (myrank_ == 0)
    {
        const string casefilename = filename_ + "_"+ field_->name() + ".case";
        ofstream casefile;
        casefile.open(casefilename.c_str());
        casefile << "# created using post_drt_ensight\n"<< "FORMAT\n\n"<< "type:\tensight gold\n";

        casefile << "\nGEOMETRY\n\n";
        casefile << "model:\t"<<timesetnumbermap_["geo"]<<"\t"<<filesetnumbermap_["geo"]<<"\t"<< geofilename<< "\n";

        casefile << "\nVARIABLE\n\n";
        casefile << GetVariableSection(filesetmap_, variablenumdfmap_, variablefilenamemap_);

        casefile << "\nTIME\n\n";
        casefile << GetTimeSectionString(geotimeset, geotime);
        casefile << GetTimeSectionString(soltimeset, soltime);

        casefile << "\nFILE\n\n";
        casefile << GetFileSectionStringFromFilesets(filesetmap_);

        casefile.close();
    }

    return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFluidEnsightWriter::WriteGeoFile(
        const string& geofilename)
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
        const XFEM::DomainIntCells domainintcells = ih->GetDomainIntCells(elegid, dis->gElement(elegid)->Shape());
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
                counter += cell->NumNode();
                break;
            }
            case DRT::Element::hex27:
            {
                counter += 20;
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
    const Epetra_Map* elementmap = dis->ElementRowMap();
    dsassert(elementmap->NumMyElements() == elementmap->NumGlobalElements(),
            "xfem filter cannot be run in parallel");
    const int nsd = 3; ///< number of space dimensions

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
                const XFEM::DomainIntCells domainintcells = ih->GetDomainIntCells(elegid, actele->Shape());
                for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
                {
                    if (cell->Shape() == distypeiter)
                    {
                        const vector<vector<double> > xarray = cell->GetPhysicalCoord(*actele);
                        int numnode = cell->NumNode();
                        if (distypeiter == DRT::Element::hex27)
                        {
                            numnode = 20;
                        }
                        for (int inen = 0; inen < numnode; ++inen)
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
string XFluidEnsightWriter::GetEnsightString(
        const DRT::Element::DiscretizationType distype) const
{
    map<DRT::Element::DiscretizationType, string>::const_iterator entry;
    switch (distype)
    {
    case DRT::Element::hex27:
        entry = distype2ensightstring_.find(DRT::Element::hex20);
        break;
    case DRT::Element::quad9:
        entry = distype2ensightstring_.find(DRT::Element::quad8);
        break;
    default:
        entry = distype2ensightstring_.find(distype);
    }
    if (entry == distype2ensightstring_.end())
        dserror("no entry in distype2ensightstring_ found");
    return entry->second;
}

/*----------------------------------------------------------------------*
 | write node connectivity for every element                  gjb 12/07 |
 *----------------------------------------------------------------------*/
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
        const int ne = iter->second;
        const string ensightCellType = GetEnsightString(distypeiter);

        if (myrank_ == 0)
        {
            cout << "writing "<< iter->second<< " "<< DistypeToString(distypeiter) << " element(s) as "
            << ne << " " << ensightCellType << " ensight cell(s)..." << endl;
            Write(geofile, ensightCellType);
            Write(geofile, ne);
        }


        // loop all available elements
        for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
        {
            const int elegid = elementmap->GID(iele);
            const XFEM::DomainIntCells domainintcells = ih->GetDomainIntCells(elegid, dis->gElement(elegid)->Shape());
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
                        // standard case with direct support
                        const int numnp = 20;
                        for (int inode=0; inode<numnp; ++inode)
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
 * \brief parse all elements and get the global(!) number of elements for each distype
 */
NumElePerDisType XFluidEnsightWriter::GetNumElePerDisType(
        const RefCountPtr<DRT::Discretization> dis,
        const RefCountPtr<XFEM::InterfaceHandle> ih ///< interfacehandle
        ) const
{
    const Epetra_Map* elementmap = dis->ElementRowMap();

    NumElePerDisType numElePerDisType;
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
        const int elegid = elementmap->GID(iele);
        const XFEM::DomainIntCells domainintcells = ih->GetDomainIntCells(elegid, dis->gElement(elegid)->Shape());
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
    // Intersection
    RCP<XFEM::InterfaceHandle> ih = rcp(new XFEM::InterfaceHandle(field_->discretization(),
                                                                  cutterfield_->discretization()));
    // apply enrichments
    RCP<XFEM::DofManager> dofman = rcp(new XFEM::DofManager(ih));

    // tell elements about the dofs and the integration
    {
        ParameterList eleparams;
        eleparams.set("action","store_xfem_info");
        eleparams.set("dofmanager",dofman);
        eleparams.set("assemble matrix 1",false);
        eleparams.set("assemble matrix 2",false);
        eleparams.set("assemble vector 1",false);
        eleparams.set("assemble vector 2",false);
        eleparams.set("assemble vector 3",false);
        field_->discretization()->Evaluate(eleparams,null,null,null,null,null);
    }

    // ensure that degrees of freedom in the discretization have been set
    field_->discretization()->FillComplete();

    PostResult result = PostResult(field_);
    result.next_result();
    if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
        return;

    // new for file continuation
    bool multiple_files = false;
    const int numdf = fieldset.size();

    // open file
    const string filename = filename_ + "_"+ field_->name() + "."+ name;
    ofstream file;
    if (myrank_==0)
    {
        file.open(filename.c_str());
    }

    map<string, vector<ofstream::pos_type> > resultfilepos;
    const int startfilepos = file.tellp(); // file position should be zero, but we stay flexible

    cout<<"writing node-based field "<<name<<endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

    WriteNodalResultStep(file, result, resultfilepos, groupname, name, fieldset, ih, dofman);
    // how many bits are necessary per time step (we assume a fixed size)?
    const int stepsize = ((int) file.tellp())-startfilepos;
    if (stepsize <= 0) dserror("found invalid step size for result file");

    while (result.next_result())
    {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
        {
            dserror("nope");
            FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
        }
        WriteNodalResultStep(file, result, resultfilepos, groupname, name, fieldset, ih, dofman);
    }

    // store information for later case file creation
    filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
    variablenumdfmap_[name] = numdf;
    variablefilenamemap_[name] = filename;

    // append index table
    WriteIndexTable(file, resultfilepos[name]);
    resultfilepos[name].clear();

    // close result file
    if (file.is_open())
        file.close();

    return;
}


vector<double> computeScalarCellNodeValues(
        DRT::Element&  ele,
        const RCP<XFEM::InterfaceHandle>&  ih,
        const XFEM::ElementDofManager& dofman,
        const XFEM::DomainIntCell& cell,
        const XFEM::PHYSICS::Field field,
        const blitz::Array<double,1> elementvalues
        )
{
    // return value
    vector<double> cellvalues;

    const int nen_cell = DRT::Utils::getNumberOfElementNodes(cell.Shape());
    const int numparam  = dofman.NumDofPerField(field);

    const int maxnod = 27;
    const int nsd = 3;

    blitz::Range _  = blitz::Range::all();

    // if cell node is on the interface, the value is not defined for a jump.
    // however, we approach the interface from one particular side and therefore,
    // -> we use the center of the cell to determine, where we come from
    const blitz::Array<double,1> cellcenterpos(cell.GetPhysicalCenterPosition(ele));

    // cell corner nodes
    for (int inen = 0; inen < nen_cell; ++inen)
    {
        const vector<double> cellnodeposvector = cell.GetPhysicalCoord(ele)[inen];
        blitz::Array<double,1> cellnodepos(3);
        for (int isd = 0; isd < nsd; ++isd) {
            cellnodepos(isd) = cellnodeposvector[isd];
        }

        // shape functions
        blitz::Array<double,1> funct(maxnod);
        DRT::Utils::shape_function_3D(
                funct,
                cell.GetDomainCoord()[inen][0],
                cell.GetDomainCoord()[inen][1],
                cell.GetDomainCoord()[inen][2],
                ele.Shape());

        blitz::Array<double,1> enr_funct(numparam);
        XFEM::ComputeEnrichedShapefunction(ele, ih, dofman, field, cellcenterpos, funct, enr_funct);
        // interpolate value
        const double x = blitz::sum(elementvalues * enr_funct);

        // store position
        cellvalues.push_back(x);
    }
    return cellvalues;
}




/*!
 \brief Write nodal values for one timestep

 No node has to have the same number of dofs.
 */
void XFluidEnsightWriter::WriteNodalResultStep(
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
    //-------------------------------------------
    // write some key words and read result data
    //-------------------------------------------

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

        // for each found distype, write block of the same typed elements
        int counter = 0;
        for (NumElePerDisType::const_iterator iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
        {
            const DRT::Element::DiscretizationType distypeiter = iter->first;
            for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
            {
                const int elegid = elementmap->GID(iele);
                DRT::Element* const actele = dis->gElement(elegid);

                // create local copy of information about dofs
                const XFEM::ElementDofManager eledofman = dofman->constructElementDofManager(*actele);

                vector<int> lm;
                vector<int> lmowner;
                actele->LocationVector(*(field_->discretization()),lm,lmowner);

                // extract local values from the global vector
                vector<double> myvelnp(lm.size());

                DRT::Utils::ExtractMyValues(*data,myvelnp,lm);

                const int numparam = eledofman.NumDofPerField(field);
                const vector<int> dofpos = eledofman.LocalDofPosPerField(field);
                //cout << XFEM::PHYSICS::physVarToString(field) << ": numparam = " << numparam << ": lm.size() = " << lm.size() << endl;

                blitz::Array<double,1> elementvalues(numparam);
                for (int iparam=0; iparam<numparam; ++iparam)   elementvalues(iparam) = myvelnp[dofpos[iparam]];

                const XFEM::DomainIntCells domainintcells = ih->GetDomainIntCells(elegid, actele->Shape());
                for (XFEM::DomainIntCells::const_iterator cell = domainintcells.begin(); cell != domainintcells.end(); ++cell)
                {
                    if (cell->Shape() == distypeiter)
                    {
                        const vector<double> cellvalues = computeScalarCellNodeValues(*actele, ih, eledofman, *cell, field, elementvalues);
                        int numnode = cell->NumNode();
                        if (distypeiter == DRT::Element::hex27)
                        {
                            numnode = 20;
                        }
                        for (int inode = 0; inode < numnode; ++inode)
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
