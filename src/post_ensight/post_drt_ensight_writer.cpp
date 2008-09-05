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


//! 6 Surfaces of a Hex27 element with 9 nodes per surface
const int Hex20_BaciToEnsightGold[20] = 
        { 0,  1,  2,  3,
          4,  5,  6,  7,
          8,  9, 10, 11,
         16, 17, 18, 19,
         12, 13, 14, 15};

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
EnsightWriter::EnsightWriter(PostField* field,
                             const string filename)
  : field_(field),
    filename_(filename),
    myrank_(((field->problem())->comm())->MyPID()),
    nodeidgiven_(true)
{
  // initialize proc0map_ correctly
  const RCP<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* noderowmap = dis->NodeRowMap();
  proc0map_ = LINALG::AllreduceEMap(*noderowmap,0);

  // sort proc0map_ so that we can loop it and get nodes in ascending order.
  std::vector<int> sortmap;
  sortmap.reserve(proc0map_->NumMyElements());
  sortmap.assign(proc0map_->MyGlobalElements(), proc0map_->MyGlobalElements()+proc0map_->NumMyElements());
  std::sort(sortmap.begin(), sortmap.end());
  proc0map_ = Teuchos::rcp(new Epetra_Map(-1, sortmap.size(), &sortmap[0], 0, proc0map_->Comm()));

  // get the number of elements for each distype (global numbers)
  numElePerDisType_ = GetNumElePerDisType(dis);

  // get the global ids of elements for each distype (global numbers)
  eleGidPerDisType_ = GetEleGidPerDisType(dis, numElePerDisType_);

  // map between distypes in BACI and existing Ensight strings
  // it includes only strings for cell types known in ensight
  // you need to manually switch to other types distypes before querying this map
  distype2ensightstring_.clear();
  distype2ensightstring_[DRT::Element::line2] = "bar2";
  distype2ensightstring_[DRT::Element::hex8] = "hexa8";
  distype2ensightstring_[DRT::Element::hex20] = "hexa20";
  distype2ensightstring_[DRT::Element::tet4] = "tetra4";
  distype2ensightstring_[DRT::Element::tet10] = "tetra10";
  distype2ensightstring_[DRT::Element::nurbs8] = "hexa8";
  distype2ensightstring_[DRT::Element::nurbs27] = "hexa8";
  distype2ensightstring_[DRT::Element::nurbs4] = "quad4";
  distype2ensightstring_[DRT::Element::nurbs9] = "quad4";
  distype2ensightstring_[DRT::Element::quad4] = "quad4";
  distype2ensightstring_[DRT::Element::quad8] = "quad8";
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
  unsigned int numsoltimes = soltime.size();

  ///////////////////////////////////
  //  write geometry file          //
  ///////////////////////////////////
  const string geofilename = filename_ + "_"+ field_->name() + ".geo";
  const size_t found_path = geofilename.find_last_of("/\\");
  const string geofilename_nopath = geofilename.substr(found_path+1);
  WriteGeoFile(geofilename);
  vector<int> filesteps;
  filesteps.push_back(1);
  filesetmap_["geo"] = filesteps;
  vector<double> timesteps;
  if (soltime.size()>0)
    timesteps.push_back(soltime[0]);
  else
    timesteps.push_back(0.0);
  timesetmap_["geo"] = timesteps;
  // at the moment, we can only print out the first step -> to be changed

  ///////////////////////////////////
  //  write solution fields files  //
  ///////////////////////////////////
  WriteAllResults(field_);

  // prepare the time sets and file sets for case file creation
  int setcounter = 0;
  int allresulttimeset = 0;
  for (map<string,vector<double> >::const_iterator entry = timesetmap_.begin(); entry != timesetmap_.end(); ++entry)
  {
    string key = entry->first;
    if ((entry->second).size()== numsoltimes)
    {
      if (allresulttimeset == 0)
      {
        setcounter++;
        allresulttimeset = setcounter;
      }
      timesetnumbermap_[key] = allresulttimeset; // reuse the default result time set, when possible
    }
    else
    {
      setcounter++;
      timesetnumbermap_[key] = setcounter; // a new time set number is needed
    }
  }

  setcounter = 0;
  for (map<string,vector<int> >::const_iterator entry = filesetmap_.begin(); entry != filesetmap_.end(); ++entry)
  {
    setcounter++;
    string key = entry->first;
    filesetnumbermap_[key] = setcounter;
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
    casefile << "model:\t"<<timesetnumbermap_["geo"]<<"\t"<<filesetnumbermap_["geo"]<<"\t"<< geofilename_nopath<< "\n";

    casefile << "\nVARIABLE\n\n";
    casefile << GetVariableSection(filesetmap_, variablenumdfmap_, variablefilenamemap_);

    casefile << "\nTIME\n\n";
    casefile << GetTimeSectionStringFromTimesets(timesetmap_);

    casefile << "\nFILE\n\n";
    casefile << GetFileSectionStringFromFilesets(filesetmap_);

    casefile.close();
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteGeoFile(const string& geofilename)
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
  Write(file, "part");
  Write(file, field_->field_pos()+1);
  Write(file, field_->name() + " field");


  //switch between nurbs an others
  if(field_->problem()->SpatialApproximation()=="Nurbs")
  {

    // cast dis to NurbsDiscretisation
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*(field_->discretization())));

    if(nurbsdis==NULL)
    {
      dserror("This probably isn't a NurbsDiscretization\n");
    }

    // get number of patches
    int npatches = (nurbsdis->GetKnotVector())->ReturnNP();

    int totalnumvisp=0;

    // loop all patches
    for(int np=0;np<npatches;++np)
    {
      // get nurbs dis' knotvector sizes
      vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np));

      int numvisp=1;

      for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
      {
	numvisp*=2*(nele_x_mele_x_lele[rr])+1;
      }
      totalnumvisp+=numvisp;
    }
    if(myrank_==0)
    {
      cout << "Writing coordinates for " << totalnumvisp << " visualisation points\n";
    }
    Write(file, "coordinates");
    Write(file, totalnumvisp);
  }
  else
  {
    Write(file, "coordinates");
    Write(file, field_->num_nodes());
  }

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
  const RefCountPtr<DRT::Discretization> dis
  )
{

  if (myrank_==0)
  {
    cout << "(computing) coordinates for  a ";
    cout << field_->problem()->SpatialApproximation();
    cout << " approximation\n";
  }

  // map for all visualisation points after they have been
  // communicated to proc 0
  RefCountPtr<Epetra_Map> proc0map;

  // refcountpointer to vector of all coordinates
  // distributed among all procs
  RefCountPtr<Epetra_MultiVector> nodecoords;

  const int NSD = 3; // number of space dimensions

  if(field_->problem()->SpatialApproximation()=="Polynomial")
  {
    const Epetra_Map* nodemap = dis->NodeRowMap();
    const int numnp = nodemap->NumMyElements();
    const int numnpglobal = nodemap->NumGlobalElements();
    nodecoords = rcp(new Epetra_MultiVector(*nodemap,3));

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
    proc0map = LINALG::AllreduceEMap(*nodemap,0);

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
  }
  else if(field_->problem()->SpatialApproximation()=="Nurbs")
  {
    // the ids of the visualisation points on this proc
    vector<int> local_vis_point_ids;
    local_vis_point_ids.clear();

    // the coordinates of the visualisation points on this proc
    // used to construct the multivector nodecoords
    vector<vector<double> > local_vis_point_x;
    local_vis_point_x.clear();

    // cast dis to NurbsDiscretisation
    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

    if(nurbsdis==NULL)
    {
      dserror("This probably isn't a NurbsDiscretization\n");
    }

    // get dimension
    int dim = (nurbsdis->Return_nele_x_mele_x_lele(0)).size();

    // get the knotvector itself
    RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    // detrmine number of patches from knotvector
    int npatches=knots->ReturnNP();

    // get vispoint offsets among patches
    vector<int> vpoff(npatches);

    vpoff[0]=0;

    // loop all patches
    for(int np=1;np<npatches;++np)
    {
      // get nurbs dis' knotvector sizes
      vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np-1));

      int numvisp=1;

      for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
      {
	numvisp*=2*nele_x_mele_x_lele[rr]+1;
      }

      vpoff[np]=vpoff[np-1]+numvisp;
    }

    // get element map
    const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));
      DRT::Node**   nodes = actele->Nodes();

      // get gid, location in the patch
      int gid = actele->Id();

      vector<int> ele_cart_id(dim);
      int np=-1;
      knots->ConvertEleGidToKnotIds(gid,np,ele_cart_id);

      // get nurbs dis' element numbers
      vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np));

      // want to loop all control points of the element,
      // so get the number of points
      const int numnp = actele->NumNode();

      // access elements knot span
      std::vector<blitz::Array<double,1> > knots(dim);
      (*((*nurbsdis).GetKnotVector())).GetEleKnots(knots,actele->Id());

      // aquire weights from nodes
      blitz::Array<double,1> weights(numnp);

      for (int inode=0; inode<numnp; ++inode)
      {
	DRT::NURBS::ControlPoint* cp
	  =
	  dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

	weights(inode) = cp->W();
      }

      // get shapefunctions, compute all visualisation point positions
      blitz::Array<double,1> nurbs_shape_funct(numnp);

      switch (actele->Shape())
      {
      case DRT::Element::nurbs4:
      {
	// element local point position
	blitz::Array<double,1> uv(2);

	// standard

	// 3           4
	//  X---------X
	//  |         |
	//  |         |
	//  |         |
	//  |         |
	//  |         |
	//  X---------X
	// 1           2
	// append 4 points
	local_vis_point_ids.push_back((2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]  );
	local_vis_point_ids.push_back((2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+1);
	local_vis_point_ids.push_back((2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]  );
	local_vis_point_ids.push_back((2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+1);

	// temporary x vector
	std::vector<double> x(3);
	x[2]=0;

	// point 1
	uv(0)= -1.0;
	uv(1)= -1.0;
	DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
					      uv               ,
					      knots            ,
					      weights          ,
					      actele->Shape()  );
	for (int isd=0; isd<2; ++isd)
	{
	  double val = 0;
	  for (int inode=0; inode<numnp; ++inode)
	  {
	    val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	  }
	  x[isd]=val;
	}
	local_vis_point_x.push_back(x);

	// point 2
	uv(0)=  1.0;
	uv(1)= -1.0;
	DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
					      uv               ,
					      knots            ,
					      weights          ,
					      actele->Shape()  );
	for (int isd=0; isd<2; ++isd)
	{
	  double val = 0;
	  for (int inode=0; inode<numnp; ++inode)
	  {
	    val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	  }
	  x[isd]=val;
	}
	local_vis_point_x.push_back(x);

	// point 3
	uv(0)= -1.0;
	uv(1)=  1.0;
	DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
					      uv               ,
					      knots            ,
					      weights          ,
					      actele->Shape()  );
	for (int isd=0; isd<2; ++isd)
	{
	  double val = 0;
	  for (int inode=0; inode<numnp; ++inode)
	  {
	    val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	  }
	  x[isd]=val;
	}
	local_vis_point_x.push_back(x);

	// point 4
	uv(0)= 1.0;
	uv(1)= 1.0;
	DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
					      uv               ,
					      knots            ,
					      weights          ,
					      actele->Shape()  );
	for (int isd=0; isd<2; ++isd)
	{
	  double val = 0;
	  for (int inode=0; inode<numnp; ++inode)
	  {
	    val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	  }
	  x[isd]=val;
	}
	local_vis_point_x.push_back(x);

	break;
      }
      case DRT::Element::nurbs9:
      {
	// element local point position
	blitz::Array<double,1> uv(2);

	{
	  // standard

	  //
	  //  +---------+
	  //  |         |
	  //  |         |
	  //  X    X    |
	  // 3|   4     |
	  //  |         |
	  //  X----X----+
	  // 1    2
	  // append 4 points
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]  );
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+1);
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]  );
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+1);

	  // temporary x vector
	  std::vector<double> x(3);
	  x[2]=0;

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 2
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 4
	  uv(0)= 0.0;
	  uv(1)= 0.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	if(ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{
	  // top line

	  //
	  //  X----X----+
	  // 5|   6     |
	  //  |         |
	  //  X    X    |
	  // 3|   4     |
	  //  |         |
	  //  X----X----+
	  // 1    2
	  //

	  // append points 5 and 6
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]  );
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+1);

	  // temporary x vector
	  std::vector<double> x(3);
	  x[2]=0;

	  // point 5
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );

	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 6
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0])
	{
	  // right line

	  //
	  //  +---------+
	  //  |         |
	  //  |         |
	  //  X    X    X
	  // 4|   5    6|
	  //  |         |
	  //  X----X----X
	  // 1    2    3

	  // append points 3 and 6
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+2);
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+2);

	  // temporary x vector
	  std::vector<double> x(3);
	  x[2]=0;

	  // point 3
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 6
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}
	if(ele_cart_id[1]+1==nele_x_mele_x_lele[1]
	   &&
	   ele_cart_id[0]+1==nele_x_mele_x_lele[0])
	{
	  // top right corner

	  //
	  //  X----X----X
	  // 7|   8    9|
	  //  |         |
	  //  X    X    X
	  // 4|   5    6|
	  //  |         |
	  //  X----X----X
	  // 1    2    3

	  // append point 9
	  local_vis_point_ids.push_back(vpoff[np]+(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1)+2*ele_cart_id[0]+2);

	  // temporary x vector
	  std::vector<double> x(3);
	  x[2]=0;

	  // point 9
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<2; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}
	break;
      }
      case DRT::Element::nurbs27:
      {
	// element local point position
	blitz::Array<double,1> uv(3);

	int idu;
	int idv;
	int idw;

	// number of visualisation points in u direction
	int nvpu=2*(nurbsdis->Return_nele_x_mele_x_lele(np))[0]+1;

	// number of visualisation points in v direction
	int nvpv=2*(nurbsdis->Return_nele_x_mele_x_lele(np))[1]+1;

	{
	  // standard

	  //               v
	  //              /
	  //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /         /  |
	  //  +---------+   |
	  //  | A----A  |   |
	  //  |/|   /|  |   +
	  //  A----A |  |  /
	  //  | A--|-A  | /
	  //  |/   |/   |/
	  //  A----A----+ ----->u
	  //
	  // append 8 points

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);


	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);


	  // point 2
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 5
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 6
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 7
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 8
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0])
	{

	  //               v
	  //              /
	  //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /         /  |
	  //  +---------+   |
	  //  | X----X--|-A |
	  //  |/|   /|  |/| +
	  //  X----X----A |/
	  //  | X--|-X--|-A
	  //  |/   |/   |/
	  //  X----X----A ----->u
	  //
	  // append 4 additional points

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);


	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);


	  // point 2
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 3
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 4
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	if(ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{

	  //               v
	  //              /
	  //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /  A----A /  |
	  //  +---------+   |
	  //  | X----X ||   |
	  //  |/| A-/|-A|   +
	  //  X----X |/ |  /
	  //  | X--|-X  | /
	  //  |/   |/   |/
	  //  X----X----+ ----->u
	  //
	  // append 4 additional points

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);


	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);


	  // point 2
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0]
	   &&
	   ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{

	  //               v
	  //              /
	  //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /  X----X-/--A
	  //  +---------+  /|
	  //  | X----X--|-X |
	  //  |/| X-/|-X|/|-A
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  // append 2 additional points

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]  )*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]+1)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 2
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}


	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2])
	{
	  //               v
	  //              /
	  //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | A----A    / |
	  //   /    /|   /  |
	  //  A----A----+   |
	  //  | X--|-X  |   |
	  //  |/|  |/|  |   +
	  //  X----X |  |  /
	  //  | X--|-X  | /
	  //  |/   |/   |/
	  //  X----X----+ ----->u
	  //
	  //
	  // append 4 additional points

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);


	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);


	  // point 2
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	}


	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2]
	   &&
	   ele_cart_id[1]+1==nele_x_mele_x_lele[1]
	  )
	{
	  //               v
	  //              /
	  //  w          /
	  //  ^   A----A----+
	  //  |  /|   /|   /|
	  //  | X----X |  / |
	  //   /| X /| X /  |
	  //  X----X----+   |
	  //  | X--|-X ||   |
	  //  |/|  |/| X|   +
	  //  X----X |/ |  /
	  //  | X--|-X  | /
	  //  |/   |/   |/
	  //  X----X----+ ----->u
	  //
	  //
	  // append 2 additional points

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);
	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 2
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	}

	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2]
	   &&
	   ele_cart_id[0]+1==nele_x_mele_x_lele[0]
	  )
	{
	  //               v
	  //              /
	  //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | X----X----A |
	  //   /    /|   /| |
	  //  X----X----A | |
	  //  | X--|-X--|-X |
	  //  |/|  |/|  |/| +
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  //
	  // append 2 additional points

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]  )*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+1)*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);

	  // point 2
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2]
	   &&
	   ele_cart_id[1]+1==nele_x_mele_x_lele[1]
	   &&
	   ele_cart_id[0]+1==nele_x_mele_x_lele[0]
	  )
	{

	  //               v
	  //              /
	  //  w          /
	  //  ^   X----X----A
	  //  |  /|   /    /|
	  //  | X----X----X |
	  //   /| X-/|-X-/|-X
	  //  X----X----X |/|
	  //  | X--|-X--|-X |
	  //  |/| X|/|-X|/|-X
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  // append 1 additional point


	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+2)*nvpu;
	  idw=(2*ele_cart_id[2]+2)*nvpu*nvpv;
	  local_vis_point_ids.push_back(vpoff[np]+idu+idv+idw);

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						knots            ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<3; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=(((nodes[inode])->X())[isd])*nurbs_shape_funct(inode);
	    }
	    x[isd]=val;
	  }
	  local_vis_point_x.push_back(x);
	}

	break;
      }
      default:
	dserror("Unknown distype for nurbs element output\n");
      }
    }

    // construct map for visualisation points. Store it in
    // class variable for access in data interpolation
    int numvispoints = 0;

    // loop all patches
    for(int np=0;np<npatches;++np)
    {
      // get nurbs dis' knotvector sizes
      vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np));

      int numvisp=1;

      for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
      {
	numvisp*=2*(nele_x_mele_x_lele[rr])+1;
      }
      numvispoints+=numvisp;
    }

    vispointmap_ = Teuchos::rcp(new Epetra_Map(numvispoints,
					       local_vis_point_ids.size(),
					       &local_vis_point_ids[0],
					       0,
					       nurbsdis->Comm()));

    // allocate the coordinates of the vizualisation points
    nodecoords = rcp(new Epetra_MultiVector(*vispointmap_,3));

    // loop over the nodes on this proc and store the coordinate information
    for (int inode=0; inode<(int)local_vis_point_x.size(); inode++)
    {
      for (int isd=0; isd<3; ++isd)
      {
	double val = (local_vis_point_x[inode])[isd];
	nodecoords->ReplaceMyValue(inode, isd, val);
      }
    }

    //new procmap
    proc0map = LINALG::AllreduceEMap(*vispointmap_,0);

    // import my new values (proc0 gets everything, other procs empty)
    Epetra_Import proc0importer(*proc0map,*vispointmap_);
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
  }
  else
  {
    dserror("spatial approximation neither Nurbs nor Polynomial\n");
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

  // for each found distype write block of the same typed elements
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType_.begin(); iter != numElePerDisType_.end(); ++iter)
  {
    const DRT::Element::DiscretizationType distypeiter = iter->first;
    const int ne = GetNumEleOutput(distypeiter, iter->second);
    const string ensightCellType = GetEnsightString(distypeiter);

    if (myrank_ == 0)
    {
      cout << "writing "<< iter->second<< " "<< DRT::DistypeToString(distypeiter) << " element(s) as "
           << ne << " " << ensightCellType << " ensight cell(s)..." << endl;
      Write(geofile, ensightCellType);
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
        case DRT::Element::line2:
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
            if (myrank_==0) // proc0 can write its elements immediately
              Write(geofile, proc0map->LID(nodes[inode]->Id())+1);
            else // elements on other procs have to store their global node ids
              nodevector.push_back(nodes[inode]->Id());
          }
          break;
        }
        case DRT::Element::hex20:
        {
          // standard case with direct support
          const int numnp = actele->NumNode();
          for (int inode=0; inode<numnp; ++inode)
          {
            if (myrank_==0) // proc0 can write its elements immediately
              Write(geofile, proc0map->LID(nodes[Hex20_BaciToEnsightGold[inode]]->Id())+1);
            else // elements on other procs have to store their global node ids
              nodevector.push_back(nodes[Hex20_BaciToEnsightGold[inode]]->Id());
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
	case DRT::Element::nurbs4:
	{
	  // cast dis to NurbsDiscretisation
	  DRT::NURBS::NurbsDiscretization* nurbsdis
	    =
	    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

	  if(nurbsdis==NULL)
	  {
	    dserror("This probably isn't a NurbsDiscretization\n");
	  }

	  // get dimension
	  int dim = (nurbsdis->Return_nele_x_mele_x_lele(0)).size();

	  // get the knotvector itself
	  RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

	  // determine type of element

	  // get gid, location in the patch and the number of the patch
	  int gid = actele->Id();

	  int npatch  =-1;
	  vector<int> ele_cart_id(dim);
	  knots->ConvertEleGidToKnotIds(gid,npatch,ele_cart_id);

	  // number of visualisation points in u direction
	  int nvpu=(nurbsdis->Return_nele_x_mele_x_lele(npatch))[0]+1;
	  {
	    // 3           4
	    //  X---------X
	    //  |         |
	    //  |         |
	    //  |         |
	    //  |         |
	    //  |         |
	    //  X---------X
	    // 1           2

	    // append 4 elements
            if (myrank_==0) // proc0 can write its elements immediately
	    {
              Write(geofile, proc0map->LID(((ele_cart_id[1]  )*(nvpu)+ele_cart_id[0]  ))+1);
	      Write(geofile, proc0map->LID(((ele_cart_id[1]  )*(nvpu)+ele_cart_id[0]+1))+1);
              Write(geofile, proc0map->LID(((ele_cart_id[1]+1)*(nvpu)+ele_cart_id[0]+1))+1);
              Write(geofile, proc0map->LID(((ele_cart_id[1]+1)*(nvpu)+ele_cart_id[0]  ))+1);
	    }
            else // elements on other procs have to store their global node ids
	    {
	      nodevector.push_back((ele_cart_id[1]  )*(nvpu)+ele_cart_id[0]  );
	      nodevector.push_back((ele_cart_id[1]  )*(nvpu)+ele_cart_id[0]+1);
	      nodevector.push_back((ele_cart_id[1]+1)*(nvpu)+ele_cart_id[0]+1);
	      nodevector.push_back((ele_cart_id[1]+1)*(nvpu)+ele_cart_id[0]  );
	    }
	  }

	  break;
	}
	case DRT::Element::nurbs9:
	{
	  // cast dis to NurbsDiscretisation
	  DRT::NURBS::NurbsDiscretization* nurbsdis
	    =
	    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

	  if(nurbsdis==NULL)
	  {
	    dserror("This probably isn't a NurbsDiscretization\n");
	  }

          // get the knotvector itself
	  RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

	  // detrmine number of patches from knotvector
	  int npatches=knots->ReturnNP();

	  // get vispoint offsets among patches
	  vector<int> vpoff(npatches);

	  vpoff[0]=0;

	  // loop all patches
	  for(int np=1;np<npatches;++np)
	  {
	    // get nurbs dis' knotvector sizes
	    vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np-1));

	    int numvisp=1;

	    for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
	    {
	      numvisp*=2*nele_x_mele_x_lele[rr]+1;
	    }

	    vpoff[np]=vpoff[np-1]+numvisp;
	  }

	  // get dimension
	  int dim = (nurbsdis->Return_nele_x_mele_x_lele(0)).size();

	  // get gid, location in the patch
	  int gid = actele->Id();

	  int npatch  =-1;
	  vector<int> ele_cart_id(dim);
	  knots->ConvertEleGidToKnotIds(gid,npatch,ele_cart_id);

	  // number of visualisation points in u direction
	  int nvpu=2*(nurbsdis->Return_nele_x_mele_x_lele(npatch))[0]+1;

	  {
	    //
	    //  X----X----X
	    // 7|   8    9|
	    //  |         |
	    //  X    X    X
	    // 4|   5    6|
	    //  |         |
	    //  X----X----X
	    // 1    2    3

	    // append 4 elements
            if (myrank_==0) // proc0 can write its elements immediately
	    {
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]  ))+1);
	      Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+1))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]  ))+1);

              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]  ))+1);
	      Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+1))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]  ))+1);

              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+1))+1);
	      Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+2))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+2))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1))+1);

              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1))+1);
	      Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+2))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+2))+1);
              Write(geofile, proc0map->LID(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+1))+1);
	    }
            else // elements on other procs have to store their global node ids
	    {
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]  ));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+1));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]  ));

	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]  ));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+1));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]  ));

	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+1));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+2));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+2));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1));

	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+1));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+2));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+2));
	      nodevector.push_back(vpoff[npatch]+((2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+1));
	    }
	  }
	  break;
	}
	case DRT::Element::nurbs27:
	{
	  //               v
	  //              /
          //  w          /
	  //  ^   X----X----A
	  //  |  /|   /    /|
	  //  | X----X----X |
	  //   /| X-/|-X-/|-X
	  //  X----X----X |/|
	  //  | X--|-X--|-X |
	  //  |/| X|/|-X|/|-X
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  // cast dis to NurbsDiscretisation
	  DRT::NURBS::NurbsDiscretization* nurbsdis
	    =
	    dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*dis));

	  if(nurbsdis==NULL)
	  {
	    dserror("This probably isn't a NurbsDiscretization\n");
	  }

          // get the knotvector itself
	  RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

	  // detrmine number of patches from knotvector
	  int npatches=knots->ReturnNP();

	  // get vispoint offsets among patches
	  vector<int> vpoff(npatches);

	  vpoff[0]=0;

	  // loop all patches
	  for(int np=1;np<npatches;++np)
	  {
	    // get nurbs dis' knotvector sizes
	    vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np-1));

	    int numvisp=1;

	    for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
	    {
	      numvisp*=2*nele_x_mele_x_lele[rr]+1;
	    }

	    vpoff[np]=vpoff[np-1]+numvisp;
	  }

	  // get dimension
	  int dim = (nurbsdis->Return_nele_x_mele_x_lele(0)).size();

	  // get gid, location in the patch
	  int gid = actele->Id();

	  int npatch  =-1;
	  vector<int> ele_cart_id(dim);
	  knots->ConvertEleGidToKnotIds(gid,npatch,ele_cart_id);

	  // number of visualisation points in u direction
	  int nvpu=2*(nurbsdis->Return_nele_x_mele_x_lele(npatch))[0]+1;

	  // number of visualisation points in v direction
	  int nvpv=2*(nurbsdis->Return_nele_x_mele_x_lele(npatch))[1]+1;

	  if (myrank_==0) // proc0 can write its elements immediately
	  {
	    int i;

	    // bottom, left front
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    // bottom, right front
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    // bottom, left rear
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    // bottom, right rear
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    //-------------------------------------------------------------------------------

	    // top, left front
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    // top, right front
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    // top, left rear
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    // top, right rear
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;
	    Write(geofile,proc0map->LID(vpoff[npatch]+i)+1);

	  }
	  else // elements on other procs have to store their global node ids
	  {

	    int i;

	    // bottom, left front
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    // bottom, right front
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    // bottom, left rear
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    // bottom, right rear
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]  )*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    //-------------------------------------------------------------------------------

	    // top, left front
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    // top, right front
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]  )*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    // top, left rear
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]  );
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    // top, right rear
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+1)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);

	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+1)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+2);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	    i  = (2*ele_cart_id[0]+1);
	    i += (2*ele_cart_id[1]+2)*nvpu;
	    i += (2*ele_cart_id[2]+2)*nvpu*nvpv;

	    nodevector.push_back(vpoff[npatch]+i);
	  }
	}
	break;
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
 * \author gjb
 * \date 01/08
 */
NumElePerDisType EnsightWriter::GetNumElePerDisType(
  const RefCountPtr<DRT::Discretization> dis
  ) const
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
  case DRT::Element::nurbs27:
    numeleout = 8*numele;
    break;
  case DRT::Element::quad9:
    numeleout = 4*numele;
    break;
  case DRT::Element::nurbs9:
    numeleout = 4*numele;
    break;
  default:
    numeleout = numele;
  }
  return numeleout;
}


/*!
 * \brief parse all elements and get the global ids of the elements for each distype
 */
EleGidPerDisType EnsightWriter::GetEleGidPerDisType(
  const RefCountPtr<DRT::Discretization> dis,
  NumElePerDisType numeleperdistype
  ) const
{
  const Epetra_Map* elementmap = dis->ElementRowMap();

  EleGidPerDisType eleGidPerDisType;

  //allocate memory
  NumElePerDisType::const_iterator iter;
  for (iter=numElePerDisType_.begin(); iter != numElePerDisType_.end(); ++iter)
  {
    eleGidPerDisType[iter->first].reserve(iter->second);
  }

  for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
  {
    const int gid = elementmap->GID(iele);
    DRT::Element* actele = dis->gElement(gid);
    const DRT::Element::DiscretizationType distype = actele->Shape();
    // update counter for current distype
    eleGidPerDisType[distype].push_back(gid);
  }

#ifndef PARALLEL
  return eleGidPerDisType; // these are already the global numbers

#else

  // in parallel case we have to provide the ele gids located on other procs as well
  EleGidPerDisType globaleleGidPerDisType;
  NumElePerDisType::const_iterator iterator;

  for (iterator=numElePerDisType_.begin(); iterator != numElePerDisType_.end(); ++iterator)
  {
    // wait for all procs before communication is started
    (dis->Comm()).Barrier();

    // no we have to communicate everything from proc 1...proc n to proc 0
    vector<char> sblock; // sending block
    vector<char> rblock; // recieving block

    // create an exporter for communication
    DRT::Exporter exporter(dis->Comm());

    // pack my element gids of this discretization type into sendbuffer
    sblock.clear();
    DRT::ParObject::AddtoPack(sblock,eleGidPerDisType[iterator->first]);

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
        vector<int> elegids;
        // extract data from recieved package
        while (index < (int)rblock.size())
        {

          DRT::ParObject::ExtractfromPack(index,rblock,elegids);
        }
        for(int i=0;i<(int) elegids.size();++i)
        {
          globaleleGidPerDisType[iterator->first].push_back(elegids[i]);
        }
        elegids.clear();
      } // end unpack

    }// for pid
  } // for iter over type

    // note: this map is only filled on proc 0 !!!!
  return globaleleGidPerDisType;

#endif
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
  case DRT::Element::tet10:
    entry = distype2ensightstring_.find(DRT::Element::tet10);
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
void EnsightWriter::WriteAnyResults(PostField* field, const char* type, const ResultType restype)
{
  PostResult result = PostResult(field_);
  result.next_result();

  MAP_ITERATOR iter;
  init_map_iterator(&iter,result.group());

  while (next_map_node(&iter))
  {
    // We do not support multiple definitions of the same name here. We just
    // use the map node to obtain the key string. Afterward we can use normal
    // map functions to find out if this key names an element vector group.
    MAP_NODE* node = iterator_get_node(&iter);
    char* key = node->key;
    if (map_has_map(result.group(),key))
    {
      MAP* entry = map_read_map(result.group(),key);
      if (map_has_string(entry, "type", type))
      {
        int dim;
        // This is bad. We should have a generic way to find how many dofs
        // there are. Until then this is remains a special purpose routine
        // that cannot serve everybody.
        if (restype==elementbased)
          // for elements we have the number of columns
          dim = map_read_int(entry, "columns");
        else if (restype==nodebased)
          // for node the number of columns might be a same bet as well
          dim = map_read_int(entry, "columns");
        else
          // Normal dof vectors have ndim dofs per node. (But then there are
          // velocity / pressure vectors and such...)
          dim = field->problem()->num_dim();
        WriteResult(key,key,restype,dim);
      }
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteResult(const string groupname,
                                const string name,
                                const ResultType restype,
                                const int numdf,
                                const int from)
{
  PostResult result = PostResult(field_);
  bool foundit = false;
  while (result.next_result(groupname))
  {
    if (map_has_map(result.group(), groupname.c_str()))
    {
      foundit = true;
      break;
    }
  }
  if (!foundit) return;

  // new for file continuation
  bool multiple_files = false;

  // open file
  const string filename = filename_ + "_"+ field_->name() + "."+ name;
  ofstream file;
  int startfilepos = 0;
  if (myrank_==0)
  {
    file.open(filename.c_str());
    startfilepos = file.tellp(); // file position should be zero, but we stay flexible
  }

  map<string, vector<ofstream::pos_type> > resultfilepos;
  int stepsize = 0;

  // distinguish between node- and element-based results
  switch(restype)
  {
  case dofbased:
  {
    if (myrank_==0)
      cout<<"writing node-based field "<<name<<endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

    //const Epetra_Map* nodemap = field_->discretization()->NodeRowMap();
    //const int numnp = nodemap->NumGlobalElements();
    //int effnumdf = numdf;
    //if (numdf==2) effnumdf=3; // in 2D we still have to write a 3D vector with zero z-components!!!
    // get the number of bits to be written each time step
    //const int stepsize = 5*80+sizeof(int)+effnumdf*numnp*sizeof(float);

    WriteDofResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      if (map_has_map(result.group(), groupname.c_str()))
      {
        const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
        if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
        {
          FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
        }
        WriteDofResultStep(file, result, resultfilepos, groupname, name, numdf, from);
      }
    }
  }
  break;

  case nodebased:
  {
    if (myrank_==0)
      cout<<"writing node-based field "<<name<<endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "node";

    WriteNodalResultStep(file, result, resultfilepos, groupname, name, numdf);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
      }
      WriteNodalResultStep(file, result, resultfilepos, groupname, name, numdf);
    }
  }
  break;

  case elementdof:
  {
    if (myrank_==0)
      cout<<"writing element based field "<<name<<endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "element";

    WriteElementDOFResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
      }
      WriteElementDOFResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    }
  }
  break;

  case elementbased:
  {
    if (myrank_==0)
      cout<<"writing element-based field "<<name<<endl;
    // store information for later case file creation
    variableresulttypemap_[name] = "element";

    WriteElementResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    // how many bits are necessary per time step (we assume a fixed size)?
    if (myrank_==0)
    {
      stepsize = ((int) file.tellp())-startfilepos;
      if (stepsize <= 0) dserror("found invalid step size for result file");
    }
    else
      stepsize = 1; //use dummy value on other procs

    while (result.next_result(groupname))
    {
      const int indexsize = 80+2*sizeof(int)+(file.tellp()/stepsize+2)*sizeof(long);
      if (static_cast<long unsigned int>(file.tellp())+stepsize+indexsize>= FILE_SIZE_LIMIT_)
      {
        FileSwitcher(file, multiple_files, filesetmap_, resultfilepos, stepsize, name, filename);
      }
      WriteElementResultStep(file, result, resultfilepos, groupname, name, numdf, from);
    }
  }
  break;

  case no_restype:
  case max_restype:
    dserror("found invalid result type");
  default:
    dserror("Invalid output type in WriteResult");
  } // end of switch(restype)

  // store information for later case file creation
  filesetmap_[name].push_back(file.tellp()/stepsize);// has to be done BEFORE writing the index table
  variablenumdfmap_[name] = numdf;
  variablefilenamemap_[name] = filename;
  // store solution times vector for later case file creation
  {
    PostResult res = PostResult(field_); // this is needed!
    vector<double> restimes = res.get_result_times(field_->name(),groupname);
    timesetmap_[name] = restimes;
  }

  // append index table
  WriteIndexTable(file, resultfilepos[name]);
  resultfilepos[name].clear();

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
  const string filename
  ) const
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
  \brief Write nodal values for one timestep for dof-based vectors
  Each node has to have the same number of dofs.
*/
void EnsightWriter::WriteDofResultStep(ofstream& file,
                                       PostResult& result,
                                       map<string, vector<ofstream::pos_type> >& resultfilepos,
                                       const string groupname,
                                       const string name,
                                       const int numdf,
                                       const int frompid) const
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

  // determine offset of dofs in case of multiple discretizations in
  // separate files (e.g. multi-scale problems). during calculation,
  // dofs are numbered consecutively for all discretizations. in the
  // post-processing phase, when only one discretization is called,
  // numbering always starts with 0, so a potential offset needs to be
  // taken into account.
  // NOTE 1: for the pressure result vector of FLUID calculations,
  // this offset is 2 or 3, depending on the number of space dimensions.
  // NOTE 2: this command is only valid, if you use NOT MORE processors
  //         for filtering than for computation. Otherwise we have empty procs
  //         owning empty maps, and therefore epetradatamap->MinAllGID()
  //         will always return zero, resulting in a wrong offset value.
  //         This is the only reason, why not to use more (and empty) procs.
  //         All other code parts of post_drt_ensight can handle that.
  if (epetradatamap->NumMyElements()<1)
    dserror("Proc %d is empty. Do not use more procs for postprocessing than for calculation.",myrank_);
  int offset = epetradatamap->MinAllGID() - dis->DofRowMap()->MinAllGID();


  //switch between nurbs an others
  if(field_->problem()->SpatialApproximation()=="Nurbs")
  {
    // a multivector for the interpolated data
    Teuchos::RefCountPtr<Epetra_MultiVector> idata;
    idata = Teuchos::rcp(new Epetra_MultiVector(*vispointmap_,numdf));

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*field_->discretization()));

    if(nurbsdis==NULL)
    {
      dserror("This probably isn't a NurbsDiscretization\n");
    }

    // get number of patches
    int npatches = (nurbsdis->GetKnotVector())->ReturnNP();

    // assuming that dimension of the manifold is
    // equal to spatial dimension
    int dim = (int)(nurbsdis->Return_nele_x_mele_x_lele(0)).size();

    // the number of vizualisation points
    int numvispoints = 0;

    for(int np=0;np<npatches;++np)
    {
      int numvisp=1;

      // get nurbs dis' knotvector sizes
      vector<int> n_x_m_x_l(nurbsdis->Return_n_x_m_x_l(np));

      // get nurbs dis' knotvector sizes
      vector<int> degree(nurbsdis->Return_degree(np));

      for(unsigned rr=0;rr<n_x_m_x_l.size();++rr)
      {
	numvisp*=2*(n_x_m_x_l[rr]-2*degree[rr])-1;
      }
      numvispoints+=numvisp;
    } // end loop over patches

    // get the knotvector itself
    RefCountPtr<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    // get vispoint offsets among patches
    vector<int> vpoff(npatches);

    vpoff[0]=0;

    // loop all patches
    for(int np=1;np<npatches;++np)
    {
      // get nurbs dis' knotvector sizes
      vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(np-1));

      int numvisp=1;

      for(unsigned rr=0;rr<nele_x_mele_x_lele.size();++rr)
      {
	numvisp*=2*nele_x_mele_x_lele[rr]+1;
      }

      vpoff[np]=vpoff[np-1]+numvisp;
    }


    // get element map
    const Epetra_Map* elementmap = nurbsdis->ElementRowMap();

    // construct a colmap for data to have it available at
    // all elements (the critical ones are the ones at the
    // processor boundary)
    // loop all available elements
    std::set<int> coldofset;
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));

      vector<int> lm;
      vector<int> lmowner;

      // extract local values from the global vectors
      actele->LocationVector(*nurbsdis,lm,lmowner);

      for (int inode=0; inode<actele->NumNode(); ++inode)
      {

	if(name == "velocity")
	{
	  for(int rr=0;rr<dim;++rr)
	  {
	    coldofset.insert(lm[inode*(dim+1)+rr]);
	  }
	}
	else if(name == "pressure")
	{
	  coldofset.insert(lm[inode*(dim+1)+dim]);
	}
	else
	{
	  dserror("up to now, I'm only able to write velocity and pressure\n");
	}
      }
    }

    std::vector<int> coldofmapvec;
    coldofmapvec.reserve(coldofset.size());
    coldofmapvec.assign(coldofset.begin(), coldofset.end());
    coldofset.clear();
    Teuchos::RCP<Epetra_Map> coldofmap =
    Teuchos::rcp(new Epetra_Map(-1,
                                coldofmapvec.size(),
                                &coldofmapvec[0],
                                0,
                                nurbsdis->Comm()));
    coldofmapvec.clear();

    const Epetra_Map* fulldofmap = &(*coldofmap);
    const RefCountPtr<Epetra_Vector> coldata
      = Teuchos::rcp(new Epetra_Vector(*fulldofmap,true));

    // create an importer and import the data
    Epetra_Import importer((*coldata).Map(),(*data).Map());
    int imerr = (*coldata).Import((*data),importer,Insert);
    if(imerr)
    {
      dserror("import falied\n");
    }

    // loop all available elements
    for (int iele=0; iele<elementmap->NumMyElements(); ++iele)
    {
      DRT::Element* const actele = nurbsdis->gElement(elementmap->GID(iele));
      DRT::Node**   nodes = actele->Nodes();

      // get gid, location in the patch and the number of the patch
      int gid = actele->Id();

      int npatch  =-1;
      vector<int> ele_cart_id(dim);
      knots->ConvertEleGidToKnotIds(gid,npatch,ele_cart_id);

      // get nele_x_mele_x_lele array
      vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(npatch));

      // number of all control points of the element
      const int numnp = actele->NumNode();

      // access elements knot span
      std::vector<blitz::Array<double,1> > eleknots(dim);
      knots->GetEleKnots(eleknots,actele->Id());

      // aquire weights from nodes
      blitz::Array<double,1> weights(numnp);

      for (int inode=0; inode<numnp; ++inode)
      {
	DRT::NURBS::ControlPoint* cp = dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);
	weights(inode) = cp->W();
      }

      // get shapefunctions, compute all visualisation point positions
      blitz::Array<double,1> nurbs_shape_funct(numnp);

      // element local visualisation point position
      blitz::Array<double,1> uv(dim);

      // extract local values from the global vectors
      vector<int> lm;
      vector<int> lmowner;

      actele->LocationVector(*nurbsdis,lm,lmowner);

      vector<double> my_data(lm.size());
      if(name == "velocity")
      {
	my_data.resize(dim*numnp);

	for (int inode=0; inode<numnp; ++inode)
	{
	  for(int rr=0;rr<dim;++rr)
	  {
	    my_data[dim*inode+rr]=(*coldata)[(*coldata).Map().LID(lm[inode*(dim+1)+rr])];
	  }
	}
      }
      else if(name == "pressure")
      {
	my_data.resize(numnp);

	for (int inode=0; inode<numnp; ++inode)
	{
	  my_data[inode]=(*coldata)[(*coldata).Map().LID(lm[inode*(dim+1)+dim])];
	}
      }
      else
      {
	dserror("up to now, I'm only able to write velocity and pressure\n");
      }

      switch (actele->Shape())
      {
      case DRT::Element::nurbs4:
      {

	// number of visualisation points in u direction
	int nvpu=(nurbsdis->Return_nele_x_mele_x_lele(npatch))[0]+1;

	{
	  // standard

	  // 3           4
	  //  X---------X
	  //  |         |
	  //  |         |
	  //  |         |
	  //  |         |
	  //  |         |
	  //  X---------X
	  // 1           2

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID((ele_cart_id[1]  )*(nvpu)+ele_cart_id[0]  );
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID((ele_cart_id[1]  )*(nvpu)+ele_cart_id[0]+1);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }


	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID((ele_cart_id[1]+1)*(nvpu)+ele_cart_id[0]  );
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 4
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID((ele_cart_id[1]+1)*(nvpu)+ele_cart_id[0]+1);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	}
	break;
      }
      case DRT::Element::nurbs9:
      {

	int idu;
	int idv;

	// number of visualisation points in u direction
	int nvpu=2*(nurbsdis->Return_nele_x_mele_x_lele(npatch))[0]+1;

	{
	  // standard

	  //
	  //  +---------+
	  //  |         |
	  //  |         |
	  //  X    X    |
	  // 3|   4     |
	  //  |         |
	  //  X----X----+
	  // 1    2

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  idu=2*ele_cart_id[0];
	  idv=2*ele_cart_id[1]*(nvpu);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  idu=2*ele_cart_id[0]+1;
	  idv=2*ele_cart_id[1]*(nvpu);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idv+idu);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }


	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  idu=2*ele_cart_id[0];
	  idv=(2*ele_cart_id[1]+1)*(nvpu);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  idu=2*ele_cart_id[0]+1;
	  idv=(2*ele_cart_id[1]+1)*(nvpu);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	}
	// top line

	//
	//  X----X----+
	// 5|   6     |
	//  |         |
	//  X    X    |
	// 3|   4     |
	//  |         |
	//  X----X----+
	// 1    2
	//
	// two additional points

	if(ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{
	  // point 5
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+(2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]  );
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 6
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+(2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+1);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}

	// right line
	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0])
	{

	  //
	  //  +---------+
	  //  |         |
	  //  |         |
	  //  x    x    X
	  // 4|   5    6|
	  //  |         |
	  //  x----x----X
	  // 1    2    3

	  // two additional points
	  // point 5
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+(2*ele_cart_id[1]  )*(nvpu)+2*ele_cart_id[0]+2);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 6
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+(2*ele_cart_id[1]+1)*(nvpu)+2*ele_cart_id[0]+2);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}

	// top right corner
	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0]&&ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{
	  //
	  //  x----x----X
	  // 7|   8    9|
	  //  |         |
	  //  x    x    x
	  // 4|   5    6|
	  //  |         |
	  //  x----x----x
	  // 1    2    3

	  // point 9
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_2D_funct(nurbs_shape_funct,
						uv               ,
						eleknots         ,
						weights          ,
						actele->Shape()  );
	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+(2*ele_cart_id[1]+2)*(nvpu)+2*ele_cart_id[0]+2);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}
	break;
      }
      case DRT::Element::nurbs27:
      {
	// element local point position
	blitz::Array<double,1> uv(3);

	int idu;
	int idv;
	int idw;

	{
	  // standard

	  //               v
	  //              /
          //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /         /  |
	  //  +---------+   |
	  //  | A----A  |   |
	  //  |/|   /|  |   +
	  //  A----A |  |  /
	  //  | A--|-A  | /
	  //  |/   |/   |/
	  //  A----A----+ ----->u
	  //
	  // append 8 points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 5
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 6
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 7
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 8
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}

	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0])
	{

	  //               v
	  //              /
          //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /         /  |
	  //  +---------+   |
	  //  | X----X--|-A |
	  //  |/|   /|  |/| +
	  //  X----X----A |/
	  //  | X--|-X--|-A
	  //  |/   |/   |/
	  //  X----X----A ----->u
	  //
	  // append 4 additional points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );
	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 3
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 4
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}

	if(ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{

	  //               v
	  //              /
          //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /  A----A /  |
	  //  +---------+   |
	  //  | X----X ||   |
	  //  |/| A-/|-A|   +
	  //  X----X |/ |  /
	  //  | X--|-X  | /
	  //  |/   |/   |/
	  //  X----X----+ ----->u
	  //
	  // append 4 additional points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }


	  // point 2
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	}

	if(ele_cart_id[0]+1==nele_x_mele_x_lele[0]
	   &&
	   ele_cart_id[1]+1==nele_x_mele_x_lele[1])
	{

	  //               v
	  //              /
          //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | /         / |
	  //   /  X----X-/--A
	  //  +---------+  /|
	  //  | X----X--|-X |
	  //  |/| X-/|-X|/|-A
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  // append 2 additional points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  uv(2)= -1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );
	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]  )*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  uv(2)=  0.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+1)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}


	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2])
	{
	  //               v
	  //              /
          //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | A----A    / |
	  //   /    /|   /  |
	  //  A----A----+   |
	  //  | X--|-X  |   |
	  //  |/|  |/|  |   +
	  //  X----X |  |  /
	  //  | X--|-X  | /
	  //  |/   |/   |/
	  //  X----X----+ ----->u
	  //
	  //
	  // append 4 additional points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)= -1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  0.0;
	  uv(1)= -1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 3
	  uv(0)= -1.0;
	  uv(1)=  0.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 4
	  uv(0)=  0.0;
	  uv(1)=  0.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}


	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2]
	   &&
	   ele_cart_id[1]+1==nele_x_mele_x_lele[1]
	  )
	{
	  //               v
	  //              /
          //  w          /
	  //  ^   A----A----+
	  //  |  /|   /|   /|
	  //  | X----X |  / |
	  //   /| X /| X /  |
	  //  X----X----+   |
	  //  | X--|-X ||   |
	  //  |/|  |/| X|   +
	  //  X----X |/ |  /
	  //  | X--|-X  | /
	  //  |/   |/   |/
	  //  X----X----+ ----->u
	  //
	  //
	  // append 2 additional points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)= -1.0;
	  uv(1)=  1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]  );
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  0.0;
	  uv(1)=  1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+1);
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	}

	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2]
	   &&
	   ele_cart_id[0]+1==nele_x_mele_x_lele[0]
	  )
	{
	  //               v
	  //              /
          //  w          /
	  //  ^   +---------+
	  //  |  /         /|
	  //  | X----X----A |
	  //   /    /|   /| |
	  //  X----X----A | |
	  //  | X--|-X--|-X |
	  //  |/|  |/|  |/| +
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  //
	  // append 2 additional points

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)= -1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]  )*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	  // point 2
	  uv(0)=  1.0;
	  uv(1)=  0.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+1)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }
	}

	if(ele_cart_id[2]+1==nele_x_mele_x_lele[2]
	   &&
	   ele_cart_id[1]+1==nele_x_mele_x_lele[1]
	   &&
	   ele_cart_id[0]+1==nele_x_mele_x_lele[0]
	  )
	{

	  //               v
	  //              /
          //  w          /
	  //  ^   X----X----A
	  //  |  /|   /    /|
	  //  | X----X----X |
	  //   /| X-/|-X-/|-X
	  //  X----X----X |/|
	  //  | X--|-X--|-X |
	  //  |/| X|/|-X|/|-X
	  //  X----X----X |/
	  //  | X--|-X--|-X
	  //  |/   |/   |/
	  //  X----X----X ----->u
	  //
	  // append 1 additional point

	  // temporary x vector
	  std::vector<double> x(3);

	  // point 1
	  uv(0)=  1.0;
	  uv(1)=  1.0;
	  uv(2)=  1.0;
	  DRT::NURBS::UTILS::nurbs_get_3D_funct(nurbs_shape_funct,
						uv               ,
						eleknots            ,
						weights          ,
						actele->Shape()  );

	  idu=(2*ele_cart_id[0]+2);
	  idv=(2*ele_cart_id[1]+2)*(2*nele_x_mele_x_lele[0]+1);
	  idw=(2*ele_cart_id[2]+2)*(2*nele_x_mele_x_lele[1]+1)*(2*nele_x_mele_x_lele[0]+1);

	  for (int isd=0; isd<numdf; ++isd)
	  {
	    double val = 0;
	    for (int inode=0; inode<numnp; ++inode)
	    {
	      val+=my_data[numdf*inode+isd]*nurbs_shape_funct(inode);
	    }
	    int lid = (*vispointmap_).LID(vpoff[npatch]+idu+idv+idw);
	    (idata)->ReplaceMyValue(lid,isd,val);
	  }

	}

	break;
      }
      default:
       dserror("unable to visualise this as a nurbs discretisation\n");
      }
    }

    // import my new values (proc0 gets everything, other procs empty)
    Epetra_Import proc0importer(*proc0map_,*vispointmap_);
    RefCountPtr<Epetra_MultiVector> allsols = rcp(new Epetra_MultiVector(*proc0map_,numdf));
    int err = allsols->Import(*idata,proc0importer,Insert);
    if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

    // write the node results (only proc 0)
    // ensight format requires u_1 .. u_n, v_1 .. v_n, w_1 ... w_n, as for nodes
    // this is fulfilled automatically due to Epetra_MultiVector usage (columnwise writing data)
    if (myrank_==0)
    {
      double* solvals = allsols->Values();
      int numentries = (numdf*(allsols->GlobalLength()));

      // now write the solution
      for (int i=0; i<numentries; ++i)
      {
	Write(file, static_cast<float>(solvals[i]));
      }

      // 2 component vectors in a 3d problem require a row of zeros.
      // do we really need this?
      if (numdf==2)
      {
        for (int inode=0; inode<numvispoints; inode++)
        {
          Write<float>(file, 0.);
        }
      }
    }
  }
  else if(field_->problem()->SpatialApproximation()=="Polynomial")
  {
    //------------------------------------------------------
    // each processor provides its result values for proc 0
    //------------------------------------------------------

    RefCountPtr<Epetra_Map> proc0datamap;
    proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

    // contract result values on proc0 (proc0 gets everything, other procs empty)
    Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);
    RefCountPtr<Epetra_Vector> proc0data = rcp(new Epetra_Vector(*proc0datamap));
    int err = proc0data->Import(*data,proc0dataimporter,Insert);
    if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

    const Epetra_BlockMap& finaldatamap = proc0data->Map();


    //------------------------------------------------------------------
    // each processor provides its dof global id information for proc 0
    //------------------------------------------------------------------

    // would be nice to have an Epetra_IntMultiVector, instead of casting to doubles
    RefCountPtr<Epetra_MultiVector> dofgidpernodelid = rcp(new Epetra_MultiVector(*nodemap,numdf));
    dofgidpernodelid->PutScalar(-1.0);

    const int mynumnp = nodemap->NumMyElements();
    for (int idf=0; idf<numdf; ++idf)
    {
      for (int inode=0; inode<mynumnp; inode++)
      {
        DRT::Node* n = dis->lRowNode(inode);
        const double dofgid = (double) dis->Dof(n, frompid + idf) + offset;
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
            dserror("received illegal dof local id: %d", lid);
        }
      }// for idf

      // 2 component vectors in a 3d problem require a row of zeros.
      // do we really need this?
      if (numdf==2)
      {
	for (int inode=0; inode<numnp; inode++)
	{
	  Write<float>(file, 0.);
	}
      }
    } // if (myrank_==0)
  }
  else
  {
    dserror("spatial approximation neither Nurbs nor Polynomial\n");
  }


  Write(file, "END TIME STEP");
  return;
}



/*!
  \brief Write nodal values for one timestep for node-based vectors
  Each node has to have the same number of dofs.
*/
void EnsightWriter::WriteNodalResultStep(ofstream& file,
                                         PostResult& result,
                                         map<string, vector<ofstream::pos_type> >& resultfilepos,
                                         const string groupname,
                                         const string name,
                                         const int numdf) const
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

  const RefCountPtr<Epetra_MultiVector> data = result.read_multi_result(groupname);
  const Epetra_BlockMap& datamap = data->Map();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));



  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> data_proc0 = rcp(new Epetra_MultiVector(*proc0map_,numdf));
  Epetra_Import proc0dofimporter(*proc0map_,datamap);
  int err = data_proc0->Import(*data,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  //---------------
  // write results
  //---------------

  const int finalnumnode = proc0map_->NumGlobalElements();

  if (myrank_==0) // ensures pointer dofgids is valid
  {
    for (int idf=0; idf<numdf; ++idf)
    {
      for (int inode=0; inode<finalnumnode; inode++) // inode == lid of node because we use proc0map_
      {
        Write(file, static_cast<float>((*((*data_proc0)(idf)))[inode]));
      }
    }
  } // if (myrank_==0)

  Write(file, "END TIME STEP");
  return;
}


/*!
  \brief Write element dof values for one timestep

  Each element has to have the same number of dofs.
*/
void EnsightWriter::WriteElementDOFResultStep(
  ofstream& file,
  PostResult& result,
  map<string, vector<ofstream::pos_type> >& resultfilepos,
  const string groupname,
  const string name,
  const int numdof,
  const int from
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

  const RefCountPtr<DRT::Discretization> dis = field_->discretization();
  const Epetra_Map* elementmap = dis->ElementRowMap(); //local node row map

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
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

  // contract result values on proc0 (proc0 gets everything, other procs empty)
  Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);
  RefCountPtr<Epetra_Vector> proc0data = rcp(new Epetra_Vector(*proc0datamap));
  int err = proc0data->Import(*data,proc0dataimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  const Epetra_BlockMap& finaldatamap = proc0data->Map();

  //------------------------------------------------------------------
  // each processor provides its dof global id information for proc 0
  //------------------------------------------------------------------

  RefCountPtr<Epetra_MultiVector> dofgidperelementlid = rcp(new Epetra_MultiVector(*elementmap,numdof));
  dofgidperelementlid->PutScalar(-1.0);

  const int nummyelem = elementmap->NumMyElements();
  for (int idof=0; idof<numdof; ++idof)
  {
    for (int ielem=0; ielem<nummyelem; ielem++)
    {
      DRT::Element* n = dis->lRowElement(ielem);
      const double dofgid = (double) dis->Dof(n, from + idof);
      if (dofgid > -1.0)
      {
        dofgidperelementlid->ReplaceMyValue(ielem, idof, dofgid);
      }
      else
      {
        dserror("Error while creating Epetra_MultiVector dofgidperlocalnodeid");
      }
    }
  }

  // contract Epetra_MultiVector on proc0 (proc0 gets everything, other procs empty)
  RefCountPtr<Epetra_MultiVector> dofgidperelementlid_proc0 = rcp(new Epetra_MultiVector(*proc0map_,numdof));
  Epetra_Import proc0dofimporter(*proc0map_,*elementmap);
  err = dofgidperelementlid_proc0->Import(*dofgidperelementlid,proc0dofimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  const int numglobelem = elementmap->NumGlobalElements();

  //-------------------------
  // specify the element type
  //-------------------------
  // loop over the different element types present
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //---------------
    // write results
    //---------------
    if (myrank_==0)
    {
      if (eleGidPerDisType_.empty()==true) dserror("no element types available");
    }

    if (myrank_==0) // ensures pointer dofgids is valid
    {
      double* dofgids = (dofgidperelementlid_proc0->Values()); // columnwise data storage
      for (int idof=0; idof<numdof; ++idof)
      {
        for (int ielem=0; ielem<numelepertype; ielem++) // inode == lid of node because we use proc0map_
        {
          // local storage position of desired dof gid
          const int doflid = ielem + (idof*numglobelem);
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
            dserror("received illegal dof local id: %d", lid);
        }
      }
    }// for idf

    // 2 component vectors in a 3d problem require a row of zeros.
    // do we really need this?
    if (numdof==2)
    {
      for (int ielem=0; ielem<numelepertype; ielem++)
      {
        Write<float>(file, 0.);
      }
    }

  } // eledistype


  Write(file, "END TIME STEP");
  return;
}

/*----------------------------------------------------------------------*/
/*!
  \brief Write element values for one timestep

  Each element has to have the same number of dofs.
  \author gjb
  \date 01/08
*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteElementResultStep(
  ofstream& file,
  PostResult& result,
  map<string, vector<ofstream::pos_type> >& resultfilepos,
  const string groupname,
  const string name,
  const int numdf,
  const int from
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

  // read the results
  const RefCountPtr<Epetra_MultiVector> data = result.read_multi_result(groupname);
  const Epetra_BlockMap& datamap = data->Map();
  const int numcol = data->NumVectors();

  // do stupid conversion into Epetra map
  RefCountPtr<Epetra_Map> epetradatamap;
  epetradatamap = rcp(new Epetra_Map(datamap.NumGlobalElements(),
                                     datamap.NumMyElements(),
                                     datamap.MyGlobalElements(),
                                     0,
                                     datamap.Comm()));

  //------------------------------------------------------
  // each processor provides its result values for proc 0
  //------------------------------------------------------

  RefCountPtr<Epetra_Map> proc0datamap;
  proc0datamap = LINALG::AllreduceEMap(*epetradatamap,0);

  // contract result values on proc0 (proc0 gets everything, other procs empty)
  Epetra_Import proc0dataimporter(*proc0datamap,*epetradatamap);
  RefCountPtr<Epetra_MultiVector> proc0data = rcp(new Epetra_MultiVector(*proc0datamap,numcol));
  int err = proc0data->Import(*data,proc0dataimporter,Insert);
  if (err>0) dserror("Importing everything to proc 0 went wrong. Import returns %d",err);

  const Epetra_BlockMap& finaldatamap = proc0data->Map();

  //-------------------------
  // specify the element type
  //-------------------------
  if (myrank_==0)
  {
    if (eleGidPerDisType_.empty()==true) dserror("no element types available");
  }
  // loop over the different element types present
  EleGidPerDisType::const_iterator iter;
  for (iter=eleGidPerDisType_.begin(); iter != eleGidPerDisType_.end(); ++iter)
  {
    const string ensighteleString = GetEnsightString(iter->first);
    const int numelepertype = (iter->second).size();
    vector<int> actelegids(numelepertype);
    actelegids = iter->second;
    // write element type
    Write(file, ensighteleString);

    //---------------
    // write results
    //---------------
    if (myrank_==0)
    {
      if (numdf+from > numcol) dserror("violated column range of Epetra_MultiVector: %d",numcol);
      for (int col=0; col<numdf; ++col)
      {
        //extract actual column
        Epetra_Vector* datacolumn = (*proc0data)(col+from);

        for (int iele=0; iele<numelepertype; iele++)
        {
          // extract element global id
          const int gid = actelegids[iele];
          // get the dof local id w.r.t. the finaldatamap
          //int lid = datamap.LID(gid);
          int lid = finaldatamap.LID(gid);
          if (lid > -1)
          {
            Write(file, static_cast<float>((*datacolumn)[lid]));
          }
          else
            dserror("received illegal dof local id: %d", lid);
        }
      }
    } // if (myrank_==0)

    // 2 component vectors in a 3d problem require a row of zeros.
    if (numdf==2)
    {
      for (int iele=0; iele<numelepertype; iele++)
      {
        Write<float>(file, 0.);
      }
    }

  } // end iteration over eleGidPerDisType_;

    // finish writing the current time step
  Write(file, "END TIME STEP");
  return;

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void EnsightWriter::WriteIndexTable(
  ofstream& file,
  const vector<ofstream::pos_type>& filepos
  ) const
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

    // Get rid of path
    const size_t found_path = filename.find_last_of("/\\");
    const string filename_nopath = filename.substr(found_path+1);

    map<string,int>::const_iterator timeentry = timesetnumbermap_.find(key);
    if (timeentry == timesetnumbermap_.end())
      dserror("key not found!");
    const int timesetnumber = timeentry->second;

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
      filename_for_casefile = filename_nopath + "***";
    }
    else
    {
      filename_for_casefile = filename_nopath;
    }
    str << GetVariableEntryForCaseFile(numdf, setnumber, key, filename_for_casefile,timesetnumber);
  }

  return str.str();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetVariableEntryForCaseFile(
  const int numdf,
  const unsigned int fileset,
  const string name,
  const string filename,
  const int timeset
  ) const
{
  stringstream str;

  // determine the type of this result variable (node-/element-based)
  map<string,string>::const_iterator entry = variableresulttypemap_.find(name);
  if (entry == variableresulttypemap_.end())
    dserror("key not found!");
  const string restypestring = entry->second;

  // create variable entry in the case-file
  switch (numdf)
  {
  case 9:
    str << "tensor asymm per "<<restypestring<<":\t"<< timeset <<"\t"<< fileset << "\t"<< name << "\t"<< filename << "\n";
    break;
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
  \brief create string for one TIME set in the case file
*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetTimeSectionString(
  const int timeset,
  const vector<double>& times
  ) const
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
  \brief create string for the TIME section in the case file
*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetTimeSectionStringFromTimesets(
  const map<string,vector<double> >& timesetmap
  ) const
{
  stringstream s;
  map<string,vector<double> >::const_iterator timeset;
  set<int> donetimesets;

  for (timeset = timesetmap.begin(); timeset != timesetmap.end(); ++timeset)
  {
    string key = timeset->first;
    map<string,int>::const_iterator entry = timesetnumbermap_.find(key);
    if (entry == timesetnumbermap_.end())
      dserror("key not found!");
    const int timesetnumber = entry->second;
    const vector<double> soltimes = timeset->second;
    if (donetimesets.find(timesetnumber)==donetimesets.end()) // do not write redundant time sets
    {
      donetimesets.insert(timesetnumber);
      string outstring = GetTimeSectionString(timesetnumber, soltimes);
      s<< outstring<<endl;
    }
  }
  return s.str();
}


/*----------------------------------------------------------------------*/
/*!
  \brief create string for the FILE section in the case file
*/
/*----------------------------------------------------------------------*/
string EnsightWriter::GetFileSectionStringFromFilesets(
  const map<string,vector<int> >& filesetmap
  ) const
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
