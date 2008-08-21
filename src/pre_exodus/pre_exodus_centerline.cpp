//"pre_exodus_centerline.cpp"

#include "pre_exodus_centerline.H"
#include <iostream>
#include "pre_exodus_soshextrusion.H"// for Gmsh plot


using namespace EXODUS;

map<int,map<int,vector<vector<double> > > > EXODUS::EleCenterlineInfo(string& cline,EXODUS::Mesh& mymesh, const vector<double> coordcorr)
{

  if (cline=="mesh"){ // create Centerline object from NodeSet
    int centerlineid = -1;
    
    map<int,EXODUS::NodeSet> nss = mymesh.GetNodeSets();
    map<int,EXODUS::NodeSet>::const_iterator i_ns;
    // check for Centerline Nodeset
    for(i_ns=nss.begin();i_ns!=nss.end();++i_ns){
      const string myname = i_ns->second.GetName();
      if (myname.find("centerline") != string::npos) centerlineid = i_ns->first;
    }
    if (centerlineid == -1) dserror("Have not found centerline NodeSet");
    
    EXODUS::Centerline myCLine(nss.find(centerlineid)->second,mymesh.GetNodes());
    myCLine.PlotCL_Gmsh();             //generation of accordant Gmsh-file

    // get rid of helper eb where the centerline ns was based on
    map<int,RCP<EXODUS::ElementBlock> > ebs = mymesh.GetElementBlocks();
    map<int,RCP<EXODUS::ElementBlock> >::const_iterator i_eb;
    vector<int> eb_ids;
    // check for Centerline ElementBlock
    for(i_eb=ebs.begin();i_eb!=ebs.end();++i_eb){
      const string myname = i_eb->second->GetName();
      if (myname.find("centerline") != string::npos) mymesh.EraseElementBlock(i_eb->first);
      else eb_ids.push_back(i_eb->first);
    }
    
    //generation of coordinate systems
    map<int,map<int,vector<vector<double> > > > centlineinfo = EXODUS::element_cosys(myCLine,mymesh,eb_ids);

    EXODUS::PlotCosys(myCLine,mymesh,eb_ids);       //generation of accordant Gmsh-file  
    // plot mesh to gmsh
    string meshname = "centerlinemesh.gmsh";
    mymesh.PlotElementBlocksGmsh(meshname,mymesh);

    
    return centlineinfo;
    
  } else { //creation of a Centerline object from file
    cout << "Reading centerline..." << endl;
    EXODUS::Centerline myCLine(cline,coordcorr);
    cout << "...done" << endl;
    
    //myCLine.PrintPoints();
    
    myCLine.PlotCL_Gmsh();             //generation of accordant Gmsh-file
    
    // get ids of the eblocks you want to calculate the locsys's
    string identifier = "ext";
    map<int,RCP<EXODUS::ElementBlock> > ebs = mymesh.GetElementBlocks();
    map<int,RCP<EXODUS::ElementBlock> >::const_iterator i_eb;
    vector<int> eb_ids;
    for (i_eb=ebs.begin(); i_eb!=ebs.end(); ++i_eb) {
      string actname = i_eb->second->GetName();
      size_t found;
      found = actname.find(identifier);
      if (found!=string::npos) eb_ids.push_back(i_eb->first);
    }

    cout << "Generating local coosys..." << endl;
    //generation of coordinate systems
    map<int,map<int,vector<vector<double> > > > centlineinfo = EXODUS::element_cosys(myCLine,mymesh,eb_ids);
    cout << "...done" << endl;
    
    EXODUS::PlotCosys(myCLine,mymesh,eb_ids);       //generation of accordant Gmsh-file  

    // plot mesh to gmsh
    string meshname = "centerlinemesh.gmsh";
    mymesh.PlotElementBlocksGmsh(meshname,mymesh,eb_ids);

    
    return centlineinfo;
  }
  
  // weirdo impossible case
  map<int,map<int,vector<vector<double> > > > mymap;
  return mymap;
  }


/*------------------------------------------------------------------------*
 |Ctor                                                            SP 06/08|
 *------------------------------------------------------------------------*/
Centerline::Centerline(string filename,vector<double> coordcorr)
{
	//initialization of points_
	points_ = rcp(new map<int,vector<double> >);
	
	//routine to open file
	ifstream infile;
	infile.open(filename.c_str(), ifstream::in);
	
	// check
	if(!infile){
	  cout << "Could not open Centerline file: " << filename << endl;
	  dserror("Could not open Centerline file!");
	}
	
	// read in the whole file into a "table"
	// for large file this might be memory intensive!
	typedef vector<float> Row;
	vector<Row> table;
	
	while(infile){
	  string line;
	  getline(infile, line);
	  istringstream is(line);
	  Row row;
	  while (is){
	    float data;
	    is >> data;
	    row.push_back(data);
	  }
	  table.push_back(row);
	}
  infile.close();
	
  // sort the table
  int clp_id = 0;
  for (unsigned i=0; i<table.size(); i++)  
  {  
    Row row = table[i];
    if(row.size() > 1){ // if true we have a "number-row"
      vector<double> clp(3);
      for (int j = 0; j < 3; ++j){
        // correct possible offset of coords
        clp[j] = table[i][j] + coordcorr[j];  // first 3 numbers are CL coords
      }
      points_->insert(pair<int,vector<double> >(clp_id,clp)); // fill map
      ++ clp_id;
    }
  }
  
  //PrintMap(cout,*points_);
	
  /* Stefans old code to read matlab file and shift coords
	//auxiliary variables
	double d;
	int i=0,j=0;
	vector<double> CLPoint(3,0);
	
	while(infile.read((char*) &d, sizeof(d)))  //reads coordinates of points in d
	{
		//and stores all three coordinates of one point in CLPoint
	  
    // displacement of coordinate systems of the centerline and the mesh must be considered
    // therefore we transform the centerline to fit into the mesh
    double delta = 205;
    double scale = -1.0;
    // mind that coordinates also have switched! x->y y->x
		switch(i%3)
		{
		case 0: CLPoint[1] = scale*d + delta; break; 
		case 1: CLPoint[0] = d; break;
		case 2: 
			{
				CLPoint[2] = d;
				//fill points_ with current point of centerline
				points_->insert(std::pair<int,vector<double> >(j,CLPoint));
				++j;
				break;
			}
		default: break;
		}
		++i;
	}
	infile.close();
	*/
}

Centerline::Centerline(const EXODUS::NodeSet& ns, const RCP<map<int,vector<double> > > nodes)
{
  //initialization of points_
  points_ = rcp(new map<int,vector<double> >);
  set<int> nodeset = ns.GetNodeSet();
  set<int>::const_iterator it;
  int id = 0;
  //ns.Print(cout,true);
  for(it=nodeset.begin(); it!=nodeset.end(); ++it){
    vector<double> node = nodes->find(*it)->second;
    points_->insert(std::pair<int,vector<double> >(id,node));
    ++id;
  }
  
}

/*------------------------------------------------------------------------*
 |Dtor                                                            SP 06/08|
 *------------------------------------------------------------------------*/
Centerline::~Centerline()
{
}

/*------------------------------------------------------------------------*
 |displays points_ on console                                      SP 06/08|
 *------------------------------------------------------------------------*/
void Centerline::PrintPoints()
{
	for(map<int,vector<double> >::const_iterator it = points_->begin(); it !=points_->end(); ++it)
	{
		cout << it->first << ": " << it->second[0] << " " << it->second[1] << " " << it->second[2] << endl;
	}
}

/*------------------------------------------------------------------------*
 |creates gmsh-file to visualize points                            SP 06/08|
 *------------------------------------------------------------------------*/
void Centerline::PlotCL_Gmsh()
{
	ofstream gmshFile("centerline.gmsh");
	gmshFile << "View \" Centerline \" {" << endl;
	
	for(map<int,vector<double> >::const_iterator it = points_->begin(); it != points_->end(); ++it)
	{
		gmshFile << "SP(" << it->second[0] << "," << it->second[1] << "," << it->second[2] << "){" << it->first << "};" << endl;	
	}
	
	gmshFile << "};";
	gmshFile.close();
}

/*------------------------------------------------------------------------*
 |calculates distance of two 3-dim vectors                         SP 06/08|
 *------------------------------------------------------------------------*/
double EXODUS::distance3d(vector<double> v1,vector<double> v2)
{
	double distance;
	distance = sqrt((v1[0]-v2[0])*(v1[0]-v2[0])+(v1[1]-v2[1])*(v1[1]-v2[1])+(v1[2]-v2[2])*(v1[2]-v2[2]));
	return distance;
}

/*------------------------------------------------------------------------*
 |calculates difference of two 3-dim vectors                       SP 06/08|
 *------------------------------------------------------------------------*/
vector<double> EXODUS::substract3d(vector<double> v1,vector<double> v2)
{
	vector<double> result(3,0);
	
	result[0]=v1[0]-v2[0];
	result[1]=v1[1]-v2[1];
	result[2]=v1[2]-v2[2];
	
	return result;
}

/*------------------------------------------------------------------------*
 |calculates sum of two 3-dim vectors                              SP 06/08|
 *------------------------------------------------------------------------*/
vector<double> EXODUS::add3d(vector<double> v1,vector<double> v2)
{
	vector<double> result(3,0);
	
	result[0]=v1[0]+v2[0];
	result[1]=v1[1]+v2[1];
	result[2]=v1[2]+v2[2];
	
	return result;
}


/*------------------------------------------------------------------------*
 |calculates cross product of two 3-dim vectors                    SP 06/08|
 *------------------------------------------------------------------------*/
vector<double> EXODUS::cross_product3d(vector<double> v1,vector<double> v2)
{
	vector<double> result(3,0);
	
	result[0] = v1[1]*v2[2] - v1[2]*v2[1];
	result[1] = v1[2]*v2[0] - v1[0]*v2[2];
	result[2] = v1[0]*v2[1] - v1[1]*v2[0];
	
	return result;
}

/*------------------------------------------------------------------------*
 |normalizes a 3-dim vector                                       SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::normalize3d(vector<double>& v)
{
	double d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	v[0] = v[0]/d;
	v[1] = v[1]/d;
	v[2] = v[2]/d;
}

/*------------------------------------------------------------------------*
 |- creates local coordinate systems for each element             SP 06/08|
 |- returns a map containing calc. directions referred to each element    |
 *------------------------------------------------------------------------*/
map<int,map<int,vector<vector<double> > > > EXODUS::element_cosys(EXODUS::Centerline& mycline,
    const EXODUS::Mesh& mymesh, const vector<int>& eb_ids)
{
  
  map<int,vector<double> > midpoints;  // here midpoints are stored
  //mp_eb_el contains (midpoint-ID, eblock-ID, element-ID)
  map<int,pair<int,int> > mp_eb_el = mymesh.createMidpoints(midpoints,eb_ids);
	//conn_mp_cp will contain (midpoint-ID, centerpoint-ID_1, centerpoint-ID_2)
	map<int,vector<int> > conn_mp_cp;
	//auxiliary variables
	int clID, clID_2;
	vector<int> ids(2,0);
	double min_distance,temp;
	
	map<int,vector<double> > clpoints = *(mycline.GetPoints());
	
	// this search should later be replaced by a nice search-tree!
	//in this section for each element the nearest point on the centerline is searched
	//and the ids of each element midpoint and the accordant centerline points are stored
	//
	//loop over all midpoints of all elements 
	for(map<int,vector<double> >::const_iterator el_iter = midpoints.begin(); el_iter != midpoints.end(); ++el_iter)
	{
		min_distance = -1;
		
		//loop over all points of the centerline to find nearest
		for(map<int,vector<double> >::const_iterator cl_iter = clpoints.begin(); cl_iter != clpoints.end(); ++cl_iter)
		{
			temp = EXODUS::distance3d(el_iter->second,cl_iter->second);
								
			if(min_distance == -1) //just for the first step
			{
				min_distance = temp;
				clID = cl_iter->first;
			}
			else
			{
				if(min_distance > temp)
				{
				min_distance = temp;
				clID = cl_iter->first;
				}
			}
		}
		//storage of IDs in conn_mp_cp
		if (clID == ((int) mycline.GetPoints()->size())-1 )
			clID_2 = clID - 1;
		else
			clID_2 = clID + 1;
		ids[0] = clID;
		ids[1] = clID_2;
		conn_mp_cp[el_iter->first]=ids;		
	}
		
	//in this section the three directions of all local coordinate systems are calculated
	//with the aid of conn_mp_cp
	//
	//map that will be returned containing (eblock-ID, element-ID, directions of local coordinate systems)
	map<int,map<int,vector<vector<double> > > > ebID_elID_local_cosy;
	
	vector<double> r_0(3,0); 
	vector<double> r_1(3,0); 
	vector<double> r_2(3,0); 
	vector<double> r_3(3,0); 
	vector<vector<double> > directions;

	//loop over conn_mp_cp
	for(map<int,vector<int> >::const_iterator it = conn_mp_cp.begin(); it != conn_mp_cp.end(); ++it)
	{
		directions.clear();
				
		//position vector from centerline point 1 to midpoint of element
		r_0 = EXODUS::substract3d(midpoints.find(it->first)->second,mycline.GetPoints()->find(it->second[0])->second);
		normalize3d(r_0);
		//position vector from centerline point 1 to centerline point 2 (axial direction)
		r_1 = EXODUS::substract3d(mycline.GetPoints()->find(it->second[1])->second,mycline.GetPoints()->find(it->second[0])->second);
		normalize3d(r_1);
				
		//if last CLPoint has been reached
		if (it->second[0] == ((int) mycline.GetPoints()->size())-1 )
		{
			r_1[0]=-r_1[0];
			r_1[1]=-r_1[1];
			r_1[2]=-r_1[2];
		}	
		
		//r_2 = r_0 x r_1 (circumferential direction)
		r_2 = EXODUS::cross_product3d(r_0,r_1);
		normalize3d(r_2);
		
		//r_3 = r_1 x r_2 (radial direction)
		r_3 = EXODUS::cross_product3d(r_1,r_2);
		normalize3d(r_3);
				
		//directions = {r_3,r_1,r_2}
		directions.push_back(r_3);
		directions.push_back(r_1);
		directions.push_back(r_2);
		
		//ebID_elID_local_cosy(ebID,elID,directions)
	  pair<int,int> eb_el = mp_eb_el.find(it->first)->second;
    ebID_elID_local_cosy[eb_el.first][eb_el.second] = directions;
				
	}
	return ebID_elID_local_cosy;
}

/*------------------------------------------------------------------------*
 |- creates local coordinate systems for each element like element_cosys   |
 |- generates gmsh-file to visualize coordinate systems            SP 06/08|
 *------------------------------------------------------------------------*/
void EXODUS::PlotCosys(EXODUS::Centerline& mycline,const EXODUS::Mesh& mymesh, const vector<int>& eb_ids)
{
  map<int,vector<double> > midpoints; // here midpoints are stored
  //mp_eb_el contains (midpoint-ID, eblock-ID, element-ID)
  map<int,pair<int,int> > mp_eb_el = mymesh.createMidpoints(midpoints,eb_ids);
	//conn_mp_cp will contain (midpoint-ID, centerpoint-ID_1, centerpoint-ID_2)
	map<int,vector<int> > conn_mp_cp;
	//auxiliary variables
	int clID, clID_2;
	vector<int> ids(2,0);
	double min_distance,temp;
	
	map<int,vector<double> > clpoints = *(mycline.GetPoints());
		
	//in this section for each element the nearest point on the centerline is searched
	//and the ids of each element midpoint and the accordant centerline points are stored
	//
	//loop over all midpoints of all elements 
	for(map<int,vector<double> >::const_iterator el_iter = midpoints.begin(); el_iter != midpoints.end(); ++el_iter)
	{
		min_distance = -1;
		
		//loop over all points of the centerline
		for(map<int,vector<double> >::const_iterator cl_iter = clpoints.begin(); cl_iter != clpoints.end(); ++cl_iter)
		{
			temp = EXODUS::distance3d(el_iter->second,cl_iter->second);
								
			if(min_distance == -1) //just for the first step
			{
				min_distance = temp;
				clID = cl_iter->first;
			}
			else
			{
				if(min_distance > temp)
				{
				min_distance = temp;
				clID = cl_iter->first;
				}
			}
		}
		//storage of IDs in conn_mp_cp
		if (clID == ((int) mycline.GetPoints()->size())-1 )
			clID_2 = clID - 1;
		else
			clID_2 = clID + 1;
		ids[0] = clID;
		ids[1] = clID_2;
		conn_mp_cp[el_iter->first]=ids;		
	}
	
	//in this section the three directions of all local coordinate systems are calculated
	//with the aid of conn_mp_cp
	//auxiliary variables
	vector<double> r_0(3,0); 
	vector<double> r_1(3,0); 
	vector<double> r_2(3,0); 
	vector<double> r_3(3,0); 
	vector<vector<double> > directions;
	map<int,vector<vector<double> > > mpID_directions;

	//loop over conn_mp_cp
	for(map<int,vector<int> >::const_iterator it = conn_mp_cp.begin(); it != conn_mp_cp.end(); ++it)
	{
		directions.clear();
						
		//position vector from centerline point 1 to midpoint of element
		r_0 = EXODUS::substract3d(midpoints.find(it->first)->second,mycline.GetPoints()->find(it->second[0])->second);
		normalize3d(r_0);
		//position vector from centerline point 1 to centerline point 2 (axial direction)
		r_1 = EXODUS::substract3d(mycline.GetPoints()->find(it->second[1])->second,mycline.GetPoints()->find(it->second[0])->second);
		normalize3d(r_1);
				
		//if last CLPoint has been reached
		if (it->second[0] == ((int) mycline.GetPoints()->size())-1 )
		{
			r_1[0]=-r_1[0];
			r_1[1]=-r_1[1];
			r_1[2]=-r_1[2];
		}	
		
		//r_2 = r_0 x r_1 (circumferential direction)
		r_2 = EXODUS::cross_product3d(r_0,r_1);
		normalize3d(r_2);
				
		//r_3 = r_1 x r_2 (radial direction)
		r_3 = EXODUS::cross_product3d(r_1,r_2);
		normalize3d(r_3);
						
		//directions = {r_3,r_1,r_2}
		directions.push_back(r_3);
		directions.push_back(r_1);
		directions.push_back(r_2);
		
		mpID_directions[it->first] = directions;
	}
	
	//
	//generation of gmsh-file to visualize local coordinate systems
	//
	ofstream gmshFile("local_coordinate_systems.gmsh");
					
	gmshFile << "View \" local coordinate systems \" {" << endl;
					
	for(map<int,vector<vector<double> > >::iterator iti = mpID_directions.begin(); iti != mpID_directions.end(); ++iti)
	{
	  vector<double> mp = midpoints.find(iti->first)->second;
	  vector<double> r1 = iti->second[0];
    vector<double> r2 = iti->second[1];
    vector<double> r3 = iti->second[2];
	  gmshFile << "VP(" << mp[0] << "," << mp[1] << "," << mp[2] << "){" << r1[0] << "," << r1[1] << "," << r1[2] << "};" << endl; 
    gmshFile << "VP(" << mp[0] << "," << mp[1] << "," << mp[2] << "){" << 2.* r2[0] << "," << 2.* r2[1] << "," << 2.* r2[2] << "};" << endl; 
    gmshFile << "VP(" << mp[0] << "," << mp[1] << "," << mp[2] << "){" << 3.* r3[0] << "," << 3.* r3[1] << "," << 3.* r3[2] << "};" << endl; 
//		//VL(mp,mp,mp,mp+r_3,mp+r_3,mp+r_3){1,1,1,1,1,1};
//		gmshFile << "SL(" << midpoints.find(iti->first)->second[0] << "," << midpoints.find(iti->first)->second[1] << "," << midpoints.find(iti->first)->second[2] << "," << midpoints.find(iti->first)->second[0] + iti->second[0][0] << "," << midpoints.find(iti->first)->second[1] + iti->second[0][1] << "," << midpoints.find(iti->first)->second[2] + iti->second[0][2] << "){1,1,1,1,1,1};" << endl;
//			
//		//VL(mp,mp,mp,mp+r_1,mp+r_1,mp+r_1){2,2,2,2,2,2};
//		gmshFile << "SL(" << midpoints.find(iti->first)->second[0] << "," << midpoints.find(iti->first)->second[1] << "," << midpoints.find(iti->first)->second[2] << "," << midpoints.find(iti->first)->second[0] + iti->second[1][0] << "," << midpoints.find(iti->first)->second[1] + iti->second[1][1] << "," << midpoints.find(iti->first)->second[2] + iti->second[1][2] << "){2,2,2,2,2,2};" << endl;
//
//		//VL(mp,mp,mp,mp+r_2,mp+r_2,mp+r_2){3,3,3,3,3,3};
//		gmshFile << "SL(" << midpoints.find(iti->first)->second[0] << "," << midpoints.find(iti->first)->second[1] << "," << midpoints.find(iti->first)->second[2] << "," << midpoints.find(iti->first)->second[0] + iti->second[2][0] << "," << midpoints.find(iti->first)->second[1] + iti->second[2][1] << "," << midpoints.find(iti->first)->second[2] + iti->second[2][2] << "){3,3,3,3,3,3};" << endl;					
						
	}
	gmshFile << "};";
	gmshFile.close();
}
