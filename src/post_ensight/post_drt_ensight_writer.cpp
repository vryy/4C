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
#ifdef PARALLEL
    dserror("parallel configured filter not finished yet. build a serial one!");
    RefCountPtr<Epetra_Map> proc0map = WriteCoordinatesPar(file, field_->discretization());
	// update the internal map 
	proc0map_=proc0map;
    WriteCellsPar(file, field_->discretization(), proc0map);
#else
    WriteCoordinates(file, field_->discretization());
    WriteCells(file, field_->discretization());    
#endif

    Write(file, "END TIME STEP");
    return;
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
    return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
RefCountPtr<Epetra_Map> EnsightWriter::WriteCoordinatesPar(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis) const
{
	const Epetra_Map* nodemap = dis->NodeRowMap();
	const int numnp = nodemap->NumMyElements();
	const int numnpglobal = nodemap->NumGlobalElements();
	RefCountPtr<Epetra_MultiVector> nodecoords = rcp(new Epetra_MultiVector(*nodemap,3));

	const int NSD = 3; ///< number of space dimensions

	// loop over nodes on this proc and store the coordinate information
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

	// import node coordinates from ALL nodes of current discretization
	int err = allnodecoords->Import(*nodecoords,proc0importer,Insert);
	if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);   

	// write the node coordinates
	// ensight format requires x_1 .. x_n, y_1 .. y_n, z_1 ... z_n  
	// this is fulfilled automatically due to Epetra_MultiVector usage (columnwise writing data)
	if (myrank_==0)
	{
		double* coords = allnodecoords->Values();
		int numentries = (3*(allnodecoords->GlobalLength()));
		dsassert(numentries == (3*numnpglobal),"proc 0 has not all of the node coordinates");
		if (nodeidgiven_)
		{
			// write gids first (case: "node id given" is true)
			for (int inode=0; inode<proc0map->NumGlobalElements(); ++inode) // this loop is empty on other procs
			{
				// write node ids
				Write(geofile,static_cast<float>(proc0map->GID(inode))+1); // gid +1 delivers GID numbering
#if 0
				cout<<"LID:"<<inode<<" -- map GID: "<<(proc0map->GID(inode))<<endl; 
#endif
			}
		}
		// go on with the coordinate information      
		for (int i=0; i<numentries; ++i) // this loop is empty on other procs
		{
#if 0
			cout<<"Proc "<<myrank_<<" writing geofile entry "<<i<<" = "<< coords[i]<<endl;
#endif
			Write(geofile, static_cast<float>(coords[i])); // ensures writing on myrank==0
		}
	}
	
	return proc0map;
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
                case DRT::Element::tet10:
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
                default:
                    dserror("don't know, how to write this element type as a cell");
                };
            };
        };
    };
    return;
}


/*----------------------------------------------------------------------*
 | write node connectivity for every element                  gjb 12/07 |
 *----------------------------------------------------------------------*/
void EnsightWriter::WriteCellsPar(
        ofstream& geofile,
        const RefCountPtr<DRT::Discretization> dis,
        const RefCountPtr<Epetra_Map> proc0map
        ) const
{
	//const Epetra_Map* localnodemap = dis->NodeRowMap();
#if 0    
	if (!(nodemap->DistributedGlobal())) dserror("only one proc");
#endif

	const Epetra_Map* elementmap = dis->ElementRowMap();


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

		vector<int> nodevector;
		const int numnodes = 4;
		//reserve memory
		nodevector.reserve(ne*numnodes);

		// loop all available elements
		for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
		{
			DRT::Element* const actele = dis->gElement(elementmap->GID(iele));
			if (actele->Shape() == distypeiter)
			{
				DRT::Node** const nodes = actele->Nodes();
                //vector<int> nodalidsNodeIds ();
				switch (actele->Shape())
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

		// now do some communicative work for the parallel case:
		// proc 1 to proc n have to send their stored node connectivity to proc0 for writing
#ifdef PARALLEL
		WriteNodeConnectivityPar(geofile, dis, nodevector, proc0map);
#endif	

	}
	return;
}
	
	
/*!
* \brief write node connectivity information in parallel case
author gjb
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

	// pack node ids into sendbuffer
	sblock.clear();
	DRT::ParObject::AddtoPack(sblock,nodevector);

	// now we start the the communication starting with proc 1
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
		// proc 0 receives from proc with myrank_==pid
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
			while (index < (int)rblock.size())
			{

				DRT::ParObject::ExtractfromPack(index,rblock,nodeids);
			}
			// compute node id based on proc0map and write it to file
			for(int i=0;i<(int) nodeids.size();++i)
			{
				int id = (proc0map->LID(nodeids[i]))+1;
				Write(geofile, id);
			}
			nodeids.clear();
		} // end unpacking

		// for safety
		exporter.Comm().Barrier();

	}// for pid
	
#endif
	
	return;
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
    
#ifndef PARALLEL
     return numElePerDisType;
#else
    // in parallel case we have to make the ele numbers known globally
     
    // how many element discretization types exist?
    DRT::Element::DiscretizationType numeledistypes = DRT::Element::max_distype;
    vector<int> myNumElePerDisType(numeledistypes);
    // write final local numbers into a vector
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
    const int numnp = nodemap->NumGlobalElements();
    const int stepsize = 5*80+sizeof(int)+numdf*numnp*sizeof(float);


    const string filename = filename_ + "_"+ field_->name() + "."+ name;
    ofstream file;
    if (myrank_==0)
    {
    	file.open(filename.c_str());
    }
    
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

    // close result file   
    if (file.is_open())
    	file.close();

    // some output
    if (myrank_==0)
    	cout<<"writing field "<<name<<endl;

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
    const Epetra_Map* nodemap = dis->NodeRowMap(); //local node map
    const RefCountPtr<Epetra_Vector> data = result.read_result(groupname);

#ifdef PARALLEL
	int numnpglobal = nodemap->NumGlobalElements();
	
	// check if the data is distributed over several processors
	if (data->DistributedGlobal())
		;
	// if not we have to have everything on proc0

	const Epetra_BlockMap& datamapser = data->Map();

	// each processor provides its information
	RefCountPtr<Epetra_MultiVector> dofgidperlocalnodeid = rcp(new Epetra_MultiVector(*nodemap,numdf));
	//initialize
	dofgidperlocalnodeid->PutScalar(-1.0);

	const int numnp = nodemap->NumMyElements();
	for (int idf=0; idf<numdf; ++idf)
	{
		for (int inode=0; inode<numnp; inode++)
		{
			DRT::Node* n = dis->lRowNode(inode);
			//const int nodegid = proc0map_->GID(inode); // old information??? unused
			const double dofgid = (double) dis->Dof(n, from + idf);
			if (dofgid > -1.0)
			{    
				dofgidperlocalnodeid->ReplaceMyValue(inode, idf, dofgid);
			}
			else
			{
				dserror("Error while creating Epetra_MultiVector dofgidperlocalnodeid");
			}
		}
	}

	// put all coordinate information on proc 0

	RefCountPtr<Epetra_Map> epetradatamap;
	epetradatamap = rcp(new Epetra_Map(datamapser.NumGlobalElements(),
			datamapser.NumMyElements(),
			datamapser.MyGlobalElements(),
			0,
			dis->Comm()));

	RefCountPtr<Epetra_Map> proc0datamap;
	proc0datamap = DRT::Utils::AllreduceEMap(*epetradatamap,0);
	
	// import my new values (proc0 gets everything, other procs empty)
	Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);   
	RefCountPtr<Epetra_Vector> proc0data = rcp(new Epetra_Vector(*proc0datamap));

	// import node coordinates from ALL nodes of current discretization
	int err = proc0data->Import(*data,proc0dataimporter,Insert);
	if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);   
	// redirect
	const Epetra_BlockMap& datamap = proc0data->Map();

	// contract on proc0
	RefCountPtr<Epetra_MultiVector> dofgidperlocalnodeid_proc0 = rcp(new Epetra_MultiVector(*proc0map_,numdf));
	// import my new values (proc0 gets everything, other procs empty)
	Epetra_Import proc0dofimporter(*proc0map_,*nodemap);   

	// import node coordinates from ALL nodes of current discretization
	err = dofgidperlocalnodeid_proc0->Import(*dofgidperlocalnodeid,proc0dofimporter,Insert);
	if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

	// now write
	if (myrank_==0)
	{
		numnpglobal = proc0map_->NumGlobalElements();

		dsassert(proc0map_->NumMyElements()==numnpglobal,"have not all data");

		double* dofgids = (dofgidperlocalnodeid_proc0->Values()); // columnwise data storage
		for (int idf=0; idf<numdf; ++idf)
		{
			for (int inode=0; inode<numnpglobal; inode++) // inode == lid of node since use of proc0map_
			{
				// get node gid
				//const int nodegid = proc0map_->GID(inode);
				//dofgidperlocalnodeid_proc0->Map().LID(nodegid);
				// get row lid
				const int doflid =  inode;

				// sorry for converting double to int
				const int actdofgid = (int) (dofgids[doflid]);
				dsassert(actdofgid>= 0, "error");
				// now we do the usual things
				int lid = datamap.LID(actdofgid);
				if (lid > -1)
				{
					Write(file, static_cast<float>((*proc0data)[lid]));
				}
				else
					dserror("error while writing result file");
			}
		}// for idf
	} // myrank==0

      
#else
    const Epetra_BlockMap& datamap = data->Map();

    
    const int numnp = nodemap->NumGlobalElements();
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
#endif
  
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
