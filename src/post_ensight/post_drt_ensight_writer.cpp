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
    myrank_(((field->problem())->comm())->MyPID()),
    nodeidgiven_(true)
{
    // initialize proc0map_ correctly    	
    const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
    proc0map_ = DRT::Utils::AllreduceEMap(*nodemap,0);

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
void EnsightWriter::WriteGeoFile(
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
void EnsightWriter::WriteGeoFileOneTimeStep(
        ofstream& file,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string name) 
{
    vector<ofstream::pos_type>& filepos = resultfilepos[name];
    Write(file, "BEGIN TIME STEP");
    filepos.push_back(file.tellp());
    
    Write(file, field_->name() + " geometry");
    Write(file, "Comment");
    
    //nodeidgiven_ is set to true inside the class constructor
    if (nodeidgiven_)
      Write(file,"node id given");
    else
      Write(file, "node id assign");

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
    RefCountPtr<Epetra_Map> proc0map = WriteCoordinates(file, field_->discretization());
	proc0map_=proc0map; // update the internal map 
    WriteCells(file, field_->discretization(), proc0map);

    Write(file, "END TIME STEP");
    return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<Epetra_Map> EnsightWriter::WriteCoordinates(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis) const
{
	const Epetra_Map* nodemap = dis->NodeRowMap();
	const int numnp = nodemap->NumMyElements();
	const int numnpglobal = nodemap->NumGlobalElements();
	RefCountPtr<Epetra_MultiVector> nodecoords = rcp(new Epetra_MultiVector(*nodemap,3));

	const int NSD = 3; ///< number of space dimensions

	// loop over the nodes on this proc and store the coordinate information
	for (int inode=0; inode<numnp; inode++)
	{
		int gid = nodemap->GID(inode);
		const DRT::Node* actnode = dis->gNode(gid);
		for (int isd=0; isd<NSD; ++isd)
		{
			double val = ((actnode->X())[isd]);      
			nodecoords->ReplaceMyValue(inode, isd, val);
		}
	}

	// put all coordinate information on proc 0
	RefCountPtr<Epetra_Map> proc0map;
	proc0map = DRT::Utils::AllreduceEMap(*nodemap,0);

	// import my new values (proc0 gets everything, other procs empty)
	Epetra_Import proc0importer(*proc0map,*nodemap);   
	RefCountPtr<Epetra_MultiVector> allnodecoords = rcp(new Epetra_MultiVector(*proc0map,3));
	int err = allnodecoords->Import(*nodecoords,proc0importer,Insert);
	if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);   

	// write the node coordinates (only proc 0)
	// ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n  
	// this is fulfilled automatically due to Epetra_MultiVector usage (columnwise writing data)
	if (myrank_==0)
	{
		double* coords = allnodecoords->Values();
		int numentries = (3*(allnodecoords->GlobalLength()));
		dsassert(numentries == (3*numnpglobal),"proc 0 has not all of the node coordinates");
		if (nodeidgiven_)
		{
			// first write node global ids (default)
			for (int inode=0; inode<proc0map->NumGlobalElements(); ++inode)
			{
				Write(geofile,static_cast<float>(proc0map->GID(inode))+1); 
				// gid+1 delivers the node numbering of the *.dat file starting with 1
			}
		}
		// now write the coordinate information      
		for (int i=0; i<numentries; ++i) 
		{
			Write(geofile, static_cast<float>(coords[i])); 
		}
	}
	
	return proc0map;
}


/*----------------------------------------------------------------------*
 | write node connectivity for every element                  gjb 12/07 |
 *----------------------------------------------------------------------*/
void EnsightWriter::WriteCells(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis,
        const RefCountPtr<Epetra_Map>& proc0map
        ) const
{
	const Epetra_Map* elementmap = dis->ElementRowMap();

	vector<int> nodevector;
	if (myrank_>0)
	{
		//reserve sufficient memory for storing the node connectivity
		//(ghosted nodes included)
		nodevector.reserve(dis->NumMyColNodes()); 
	}

	// get the number of elements for each distype (global numbers)
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

		if (myrank_ == 0)
		{
			cout << "writing "<< iter->second<< " "<< realcellshape << " elements"
			<< " ("<< ne << " cells) per distype."<< " ensight output celltype: "
			<< ensightString<< endl;       
			Write(geofile, ensightString);
			Write(geofile, ne);
		}

		nodevector.clear();

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
				case DRT::Element::tet10:
				case DRT::Element::tri3:
				case DRT::Element::tri6:
				case DRT::Element::wedge6:
				case DRT::Element::wedge15:
				case DRT::Element::pyramid5:
				{
					// standard case with direct support
					const int numnp = actele->NumNode();
					for (int inode=0; inode<numnp; ++inode)
					{
						if (myrank_==0) // proc0 can write its elements immidiately
						    Write(geofile, proc0map->LID(nodes[inode]->Id())+1);
						else // elements on other procs have to store their global node ids
							nodevector.push_back(nodes[inode]->Id());
					}
					break;
				}
				case DRT::Element::hex27:
				{
					// write subelements
					for (int isubele=0; isubele<8; ++isubele)
						for (int isubnode=0; isubnode<8; ++isubnode)
							if (myrank_==0) // proc0 can write its elements immidiately
								Write(geofile, proc0map->LID(nodes[subhexmap[isubele][isubnode]]->Id())
										+1);
							else // elements on other procs have to store their global node ids
								nodevector.push_back(nodes[subhexmap[isubele][isubnode]]->Id());
					break;
				}
				case DRT::Element::quad9:
				{
					// write subelements
					for (int isubele=0; isubele<4; ++isubele)
						for (int isubnode=0; isubnode<4; ++isubnode)
							if (myrank_==0) // proc0 can write its elements immidiately
							Write(geofile, proc0map->LID(nodes[subquadmap[isubele][isubnode]]->Id())
									+1);
							else // elements on other procs have to store their global node ids
								nodevector.push_back(nodes[subquadmap[isubele][isubnode]]->Id());
					break;
				}

				default:
					dserror("don't know, how to write this element type as a cell");
				}
			}
		}

#ifdef PARALLEL
		// now do some communicative work for the parallel case:
		// proc 1 to proc n have to send their stored node connectivity to proc0 
		// which does the writing
		
		WriteNodeConnectivityPar(geofile, dis, nodevector, proc0map);
#endif	

	}
	return;
}
	
	
/*!
* \brief communicate and write node connectivity in parallel case

  \author gjb
  \date 12/07
*/	
void EnsightWriter::WriteNodeConnectivityPar(       
		ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis,
	    const vector<int>& nodevector,
        const RefCountPtr<Epetra_Map> proc0map) const
{
#ifdef PARALLEL
	// no we have communicate the connectivity infos from proc 1...proc n to proc 0

	vector<char> sblock; // sending block
	vector<char> rblock; // recieving block

	// create an exporter for communication
	DRT::Exporter exporter(dis->Comm());

	// pack my node ids into sendbuffer
	sblock.clear();
	DRT::ParObject::AddtoPack(sblock,nodevector);

	// now we start the communication
	for (int pid=0;pid<(dis->Comm()).NumProc();++pid)
	{
		MPI_Request request;
		int         tag    =0;   
		int         frompid=pid;
		int         topid  =0;   
		int         length=sblock.size();

		//--------------------------------------------------
		// proc pid sends its values to proc 0
		if (myrank_==pid)
		{
			exporter.ISend(frompid,topid,
					&(sblock[0]),sblock.size(),
					tag,request);
		}

		//--------------------------------------------------
		// proc 0 receives from proc pid
		rblock.clear();  
		if (myrank_ == 0)
		{
			exporter.ReceiveAny(frompid,tag,rblock,length);
			if(tag!=0)
			{
				dserror("Proc 0 received wrong message (ReceiveAny)");
			}
			exporter.Wait(request);
		}

		// for safety
		exporter.Comm().Barrier();

		//--------------------------------------------------
		// Unpack received block and write the data
		if (myrank_==0)
		{
			int index = 0;
			vector<int> nodeids;
			// extract data from recieved package
			while (index < (int)rblock.size())
			{

				DRT::ParObject::ExtractfromPack(index,rblock,nodeids);
			}
			// compute node lid based on proc0map and write it to file
			for(int i=0;i<(int) nodeids.size();++i)
			{
				// using the same map as for the writing the node coordinates
				int id = (proc0map->LID(nodeids[i]))+1;
				Write(geofile, id);
			}
			nodeids.clear();
		} // end unpack

		// for safety
		exporter.Comm().Barrier();

	}// for pid

#endif

	return;
}


/*!
 * \brief parse all elements and get the global(!) number of elements for each distype
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
    
#ifndef PARALLEL
     return numElePerDisType; // these are already the global numbers

#else
    // in parallel case we have to sum up the local element distype numbers
     
    // determine maximum number of possible element discretization types
    DRT::Element::DiscretizationType numeledistypes = DRT::Element::max_distype;

    // write the final local numbers into a vector
    vector<int> myNumElePerDisType(numeledistypes);
	NumElePerDisType::const_iterator iter;
	for (iter=numElePerDisType.begin(); iter != numElePerDisType.end(); ++iter)
	{
		const DRT::Element::DiscretizationType distypeiter = iter->first;
		const int ne = iter->second;
        myNumElePerDisType[distypeiter]+=ne;
	}
	
    // wait for all procs before communication is started
    (dis->Comm()).Barrier();
    
    // form the global sum 
	vector<int> globalnumeleperdistype(numeledistypes); 
	(dis->Comm()).SumAll(&(myNumElePerDisType[0]),&(globalnumeleperdistype[0]),numeledistypes);
  
	// create return argument containing the global element numbers per distype
	NumElePerDisType globalNumElePerDisType;
	for(int i =0; i<numeledistypes ;++i)
	{
		if (globalnumeleperdistype[i]>0) // no entry when we have no element of this type
        globalNumElePerDisType[DRT::Element::DiscretizationType(i)]=globalnumeleperdistype[i];
	}

    return globalNumElePerDisType;
    
#endif
    
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
        const int restype,
        const int numdf,
        const int from)
{
    PostResult result = PostResult(field_);
    result.next_result();
    if (!map_has_map(result.group(), const_cast<char*>(groupname.c_str())))
        return;

    // new for file continuation
    bool multiple_files = false;
    
    // open file
    const string filename = filename_ + "_"+ field_->name() + "."+ name;
    ofstream file;
    if (myrank_==0)
    {
    	file.open(filename.c_str());
    }
    
    // distinguish between node- and element-based results
    switch(restype)
    {
    case nodebased:
    {
        if (myrank_==0)     		
        	cout<<"writing node-based field "<<name<<endl;
    	
    	const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
        const int numnp = nodemap->NumGlobalElements();
        const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);
    	
        map<string, vector<ofstream::pos_type> > resultfilepos;
        WriteNodalResultStep(file, result, resultfilepos, groupname, name, numdf, from);
        while (result.next_result())
        {
            const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
            if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
            {
                FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
            }
            WriteNodalResultStep(file, result, resultfilepos, groupname, name, numdf, from);
        }
        
        // append index table
        WriteIndexTable(file, resultfilepos[name]);
        resultfilepos[name].clear();
        
        // store information for later case file creation
        filesetmap_[name].push_back(file.tellp()/stepsize);
        variablenumdfmap_[name] = numdf;
        variablefilenamemap_[name] = filename;
        variableresulttypemap_[name] = "node";
    }
        break;
        
    case elementbased:
    {
    	if (myrank_==0)     		
    		cout<<"writing element-based field "<<name<<endl;

    	const Epetra_Map* elementmap = field_->discretization()->ElementRowMap();
    	const int numele = elementmap->NumGlobalElements();
    	const int stepsize = 5*80+sizeof(int)+numdf*numele*sizeof(float);
    	
        map<string, vector<ofstream::pos_type> > resultfilepos;
        WriteElementResultStep(file, result, resultfilepos, groupname, name, numdf, from);
        while (result.next_result())
        {
        	const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        	if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
        	{
        		FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
        	}
        	WriteElementResultStep(file, result, resultfilepos, groupname, name, numdf, from);
        }
        
        // append index table
        WriteIndexTable(file, resultfilepos[name]);
        resultfilepos[name].clear();
        
        // store information for later case file creation
        filesetmap_[name].push_back(file.tellp()/stepsize);
        variablenumdfmap_[name] = numdf;
        variablefilenamemap_[name] = filename;
        variableresulttypemap_[name] = "element";
    }
    	break;
    	
    case no_restype:
    case max_restype:
    	dserror("found invalid result type");
    } // end of switch(restype)

    // close result file   
    if (file.is_open())
    	file.close();

    return;
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
	if (myrank_==0)
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
	} // if (myrank_==0)
    return;
}


/*!
 \brief Write nodal values for one timestep

 Each node has to have the same number of dofs.
 */
void EnsightWriter::WriteNodalResultStep(
        ofstream& file,
        PostResult& result,
        map<string, vector<ofstream::pos_type> >& resultfilepos,
        const string groupname,
        const string name,
        const int numdf,
        const int from) const
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
	const Epetra_Map* nodemap = dis->NodeRowMap(); //local node row map
	const int numnp = nodemap->NumGlobalElements();

	const RefCountPtr<Epetra_Vector> data = result.read_result(groupname);
	const Epetra_BlockMap& datamap = data->Map();

	// do stupid conversion into Epetra map
	RefCountPtr<Epetra_Map> epetradatamap;
	epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
			datamap.NumMyElements(),
			datamap.MyGlobalElements(),
			0,
			datamap.Comm()));
#if 0
	if (epetradatamap->PointSameAs(*proc0map_))
		cout<<"INFO: proc0map and epetradatamap are identical."<<endl;
	// check if the data is distributed over several processors
	bool isdistributed = (data->DistributedGlobal());
#endif

	//------------------------------------------------------
	// each processor provides its result values for proc 0
	//------------------------------------------------------

	RefCountPtr<Epetra_Map> proc0datamap;
	proc0datamap = DRT::Utils::AllreduceEMap(*epetradatamap,0);

	// contract result values on proc0 (proc0 gets everything, other procs empty)
	Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);   
	RefCountPtr<Epetra_Vector> proc0data = rcp(new Epetra_Vector(*proc0datamap));
	int err = proc0data->Import(*data,proc0dataimporter,Insert);
	if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

	const Epetra_BlockMap& finaldatamap = proc0data->Map();

	//------------------------------------------------------------------
	// each processor provides its dof global id information for proc 0
	//------------------------------------------------------------------

	RefCountPtr<Epetra_MultiVector> dofgidpernodelid = rcp(new Epetra_MultiVector(*nodemap,numdf));
	dofgidpernodelid->PutScalar(-1.0);

	const int mynumnp = nodemap->NumMyElements();
	for (int idf=0; idf<numdf; ++idf)
	{
		for (int inode=0; inode<mynumnp; inode++)
		{
			DRT::Node* n = dis->lRowNode(inode);
			const double dofgid = (double) dis->Dof(n, from + idf);
			if (dofgid > -1.0)
			{    
				dofgidpernodelid->ReplaceMyValue(inode, idf, dofgid);
			}
			else
			{
				dserror("Error while creating Epetra_MultiVector dofgidperlocalnodeid");
			}
		}
	}

	// contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
	RefCountPtr<Epetra_MultiVector> dofgidpernodelid_proc0 = rcp(new Epetra_MultiVector(*proc0map_,numdf));
	Epetra_Import proc0dofimporter(*proc0map_,*nodemap);   
	err = dofgidpernodelid_proc0->Import(*dofgidpernodelid,proc0dofimporter,Insert);
	if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

	//---------------
	// write results
	//---------------

	const int finalnumnode = proc0map_->NumGlobalElements();
	if (myrank_==0) // ensures pointer dofgids is valid
	{
		double* dofgids = (dofgidpernodelid_proc0->Values()); // columnwise data storage
		for (int idf=0; idf<numdf; ++idf)
		{
			for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0map_
			{
				// local storage position of desired dof gid
				const int doflid = inode + (idf*numnp);  
                // get the dof global id 
				const int actdofgid = (int) (dofgids[doflid]);
				dsassert(actdofgid>= 0, "error while getting dof global id");
                // get the dof local id w.r.t. the finaldatamap
				int lid = finaldatamap.LID(actdofgid);
				if (lid > -1)
				{
					Write(file, static_cast<float>((*proc0data)[lid]));
				}
				else
					dserror("recieved illegal dof local id: %d", lid);
			}
		}// for idf
	} // if (myrank_==0)

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
	return;
}

/*----------------------------------------------------------------------*/
/*!
 \brief Write element values for one timestep

 Each element has to have the same number of dofs.
 */
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteElementResultStep(
	       ofstream& file,
	        PostResult& result,
	        map<string, vector<ofstream::pos_type> >& resultfilepos,
	        const string groupname,
	        const string name,
	        const int numdf,
	        const int from) const
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
	//specify the element type
	Write(file, "quad4");
	
	const RefCountPtr<DRT::Discretization> dis = field_->discretization();
	const Epetra_Map* elementmap = dis->ElementRowMap(); //local element row map
	const int numele = elementmap->NumGlobalElements();

	// create dummy output
	for (int iele=0; iele<numele; iele++)
	{
		Write<float>(file, iele);
	}
	
	// finish writing the current time step
	Write(file, "END TIME STEP");
	return;

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
    return;
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
    
    return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetVariableSection(
        map<string,vector<int> >  filesetmap,
        map<string,int>           variablenumdfmap,
        map<string,string>        variablefilenamemap
        ) const
{
    stringstream str;
    
    map<string,int>::const_iterator variable;
    
    for (variable = variablenumdfmap.begin(); variable != variablenumdfmap.end(); ++variable)
    {
        const string key = variable->first;
        const int numdf = variable->second;
        const string filename = variablefilenamemap[key];
        
    	map<string,int>::const_iterator entry1 = filesetnumbermap_.find(key);
    	if (entry1 == filesetnumbermap_.end())
    		dserror("key not found!");
    	const int setnumber = entry1->second;
        
        map<string,vector<int> >::const_iterator entry2 = filesetmap.find(key);
        if (entry2 == filesetmap.end()) 
            dserror("filesetmap not defined for '%s'", key.c_str());
        
        const int numsubfilesteps = entry2->second.size();
        string filename_for_casefile;
        if (numsubfilesteps > 1)
        {
            filename_for_casefile = filename + "***";
        }
        else
        {
            filename_for_casefile = filename;
        }
            
        str << GetVariableEntryForCaseFile(numdf, setnumber, key, filename_for_casefile);
    }

    return str.str();
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
    
    // determine the type of this result variable (node-/element-based)
	map<string,string>::const_iterator entry = variableresulttypemap_.find(name);
	if (entry == variableresulttypemap_.end())
		dserror("key not found!");
	const string restypestring = entry->second;
       
    // create variable entry in the case-file
    switch (numdf)
    {
    case 6:
        str << "tensor symm per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
        break;
    case 3:
    case 2:
        str << "vector per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
        break;
    case 1:
        str << "scalar per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
        break;
    default:
        dserror("unknown number of dof per node");
    };
    return str.str();
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
    s << "time set:\t\t" << timeset << "\n"<< "number of steps:\t"<< times.size() << "\ntime values: ";
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

    for (fileset = filesetmap.begin(); fileset != filesetmap.end(); ++fileset)
    {
    	string key = fileset->first;
    	map<string,int>::const_iterator entry = filesetnumbermap_.find(key);
    	if (entry == filesetnumbermap_.end())
    		dserror("key not found!");
    	const int setnumber = entry->second;
        vector<int> stepsperfile = fileset->second;
        s << "file set:\t\t"<< setnumber << "\n";
        if (stepsperfile.size() == 1)
        {
            s << "number of steps:\t"<< stepsperfile[0] << "\n\n";
        }
        else
        {
            for (unsigned int j = 0; j < stepsperfile.size(); ++j)
            {
                s << "filename index:\t"<< 1+j << "\n";
                s << "number of steps:\t"<< stepsperfile[j] << "\n";
            }
            s << "\n";
        }
    }
    return s.str();
}



#endif
