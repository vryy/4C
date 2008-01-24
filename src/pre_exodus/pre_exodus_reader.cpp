/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_reader.cpp

\brief preprocessor reader for exodusII format 

<pre>
Maintainer: Moritz & Georg
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/frenzel
            089 - 289-15240
</pre>

Here everything related with the exodus format and the accessible data
is handed to a c++ object mesh.
*/
/*----------------------------------------------------------------------*/
#ifdef D_EXODUS
#include "pre_exodus_reader.H"
#include <time.h>

#ifdef PARALLEL
#include <mpi.h>
#endif

using namespace std;
using namespace Teuchos;

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
Mesh::Mesh(string exofilename)
{
  int error;
  int CPU_word_size,IO_word_size;
  float exoversion;                   /* version of exodus */
  CPU_word_size = sizeof(float);      /* float or double */
  IO_word_size = 0;                   /* use what is stored in file */

  //cout << "meshfilename: " << exofilename.c_str() << endl;

  const char *exofilenamechar;
  exofilenamechar=exofilename.c_str();

  // open EXODUS II file
  exoid_ = ex_open(exofilenamechar,EX_READ,&CPU_word_size,&IO_word_size,&exoversion);
  if (exoid_<=0)
	  dserror("Error while opening EXODUS II file %s",exofilenamechar);

  // print version
  cout<<exofilename<<" was created with EXODUS II library version "<<exoversion<<endl;

  // read database parameters
  error = ex_get_init (exoid_, title_, &num_dim_, &num_nodes_,&num_elem_, &num_elem_blk_, &num_node_sets_, &num_side_sets_);

  // get nodal coordinates
  float x[num_nodes_];
  float y[num_nodes_];
  float z[num_nodes_];
  error = ex_get_coord(exoid_,x,y,z);
  if (error != 0) dserror("exo error returned");
  
  // create node pointer
  myNodes_.reserve(num_nodes_);
  for (int i = 0; i < num_nodes_; ++i) {
    double coords[3];
    coords[0] = x[i];
    coords[1] = y[i];
    coords[2] = z[i];
    myNodes_[i] = rcp(new PreNode(i,coords));
  }

  // num_entities_ are all ElementBlocks, NodeSets, and SideSets together
  num_entities_ = num_elem_blk_ + num_node_sets_ + num_side_sets_;

  // allocate sufficient memory for the entities
  myEntities_.reserve(num_entities_);
  
  // entitycounter counts all ElementBlocks, NodeSets, and SideSets together
  int entitycounter = 0;
  
  int epropID[num_elem_blk_];
  error = ex_get_prop_array(exoid_, EX_ELEM_BLOCK, "ID", epropID);
  for (int i = 0; i < num_elem_blk_; ++i) 
  {
//	  //RCP<Entity> myentity = rcp(new Entity(exoid_,entitycounter,epropID[i],Entity::elem_blk));
//	  myEntities_[entitycounter] = rcp(new Entity(exoid_,entitycounter,epropID[i],Entity::elem_blk));
//	  //myEntities_[entitycounter]->Print(cout);
//	  entitycounter++;
	  
	  // Read Element Blocks into Map
    char mychar[MAX_STR_LENGTH+1];
    int num_el_in_blk,num_nod_per_elem,num_attr;
    error = ex_get_elem_block (exoid_, epropID[i], mychar, &num_el_in_blk, &num_nod_per_elem, &num_attr);
    // prefer string to store element type
    string ele_type(mychar,int (MAX_STR_LENGTH));
    
    // get ElementBlock name
    error = ex_get_name (exoid_, EX_ELEM_BLOCK, epropID[i], mychar);
    // prefer string to store name
    string blockname(mychar,int (MAX_STR_LENGTH));
    
    // get element connectivity
    int connect[num_nod_per_elem*num_el_in_blk];
    error = ex_get_elem_conn(exoid_,epropID[i],connect);
    if (error != 0) dserror("exo error returned");
	  map<int,vector<int> > eleconn;
	  for (int j = 0; j < num_el_in_blk; ++j) {
	    vector<int> actconn;
	    for (int k = 0; k < num_nod_per_elem; ++k) actconn.push_back(connect[j*k]);
	    eleconn.insert(std::pair<int,vector<int> >(j,actconn));
    }
	  ElementBlock actEleBlock(StringToShape(ele_type), eleconn, blockname);
	  
	  // Add this ElementBlock into Mesh map
	  elementBlocks_.insert(std::pair<int,ElementBlock>(i,actEleBlock));
  }
  
  // get all NodeSets
  map<int,NodeSet> prelimNodeSets;   // prelim due to possible prop names
  int npropID[num_node_sets_];
  error = ex_get_prop_array(exoid_, EX_NODE_SET, "ID", npropID);
  for (int i = 0; i < num_node_sets_; ++i) {
//    myEntities_[entitycounter] = rcp(new Entity(exoid_,entitycounter,npropID[i],Entity::node_set));
//    //myEntities_[entitycounter]->Print(cout);
//    //Entity myentity(exoid_,entitycounter,npropID[i],Entity::node_set);
//    entitycounter++;
    
    
    // Read NodeSet params
    int num_nodes_in_set,num_df_in_set;
    error = ex_get_node_set_param (exoid_, npropID[i], &num_nodes_in_set,&num_df_in_set);
    
    // get NodeSet name
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_name (exoid_, EX_NODE_SET, npropID[i], mychar);
    // prefer string to store name
    string nodesetname(mychar, int(MAX_STR_LENGTH));

    // get nodes in node set
    int node_set_node_list[num_nodes_in_set];
    error = ex_get_node_set (exoid_, npropID[i], node_set_node_list);
    if (error != 0) dserror("error reading node set");
    set<int> nodes_in_set;
    for (int j = 0; j < num_nodes_in_set; ++j) nodes_in_set.insert(node_set_node_list[j]);
    NodeSet actNodeSet(nodes_in_set,nodesetname,"none");
    
    // Add this NodeSet into Mesh map (here prelim due to pro names)
    prelimNodeSets.insert(std::pair<int,NodeSet>(i,actNodeSet));
  }
  
  /* Read NodeSet property names *********************************************
   * They are assigned by ICEM and provide recognition */
  int num_props;
  float fdum;
  char *cdum;
  error = ex_inquire (exoid_, EX_INQ_NS_PROP, &num_props, &fdum, cdum);
  // allocate memory for NodeSet property names
  char** prop_names = new char*[num_props];
  for (int i=0; i<num_props; ++i)
  {
     prop_names[i] = new char[MAX_STR_LENGTH+1];
  }
  // get prop names of node sets  
  error = ex_get_prop_names(exoid_,EX_NODE_SET,prop_names);
  
  // Add prop names to final Mesh NodeSet if available
  if ((num_props-1) == num_node_sets_){
    for (int i = 1; i < num_props; ++i) {
      string propname(prop_names[i], int(MAX_STR_LENGTH));
      
      map<int,NodeSet>::const_iterator blub = prelimNodeSets.find(i-1);
      if (blub == prelimNodeSets.end())
        dserror("impossible");
      const NodeSet actNodeSet = blub->second;
      //NodeSet actNodeSet = nodeSets_[i-1];
      const NodeSet newNodeSet(actNodeSet.GetNodeSet(),actNodeSet.GetName(),propname);
      //nodeSets_[i-1] = newNodeSet;
      nodeSets_.insert(std::pair<int,NodeSet>(blub->first,newNodeSet));
    }
  } else {
    // this is the standard case without prop names
    nodeSets_ = prelimNodeSets;
  }
  
  // get all SideSets
  int spropID[num_side_sets_];
  error = ex_get_prop_array(exoid_, EX_SIDE_SET, "ID", spropID);
  for (int i = 0; i < num_side_sets_; ++i) {
//    myEntities_[entitycounter] = rcp(new Entity(exoid_,entitycounter,spropID[i],Entity::side_set));
//    //myEntities_[entitycounter]->Print(cout);
//    //Entity myentity(exoid_,entitycounter,spropID[i],Entity::side_set);
//    entitycounter++;
    
    
    // get SideSet name
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_name (exoid_, EX_SIDE_SET, spropID[i], mychar);
    // prefer string to store name
    string sidesetname(mychar, int(MAX_STR_LENGTH));
    SideSet actSideSet(sidesetname);
    
    // Add this SideSet into Mesh map
    sideSets_.insert(std::pair<int,SideSet>(i,actSideSet));
  }


  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
Mesh::~Mesh()
{
  CloseExo();
  return;
}

/*----------------------------------------------------------------------*
 |  Close corresponding Exofile(public)                        maf 12/07|
 *----------------------------------------------------------------------*/
void Mesh::CloseExo()
{
  // close exodus II file
  int error = ex_close(exoid_);
  if (error < 0)
    dserror("error while closing exodus II file");
}


/*----------------------------------------------------------------------*
 |  Print method (public)                                      maf 12/07|
 *----------------------------------------------------------------------*/
void Mesh::Print(ostream & os, bool verbose) const
{
  os << "Mesh consists of ";
  os << num_nodes_ << " Nodes, ";
  os << num_elem_ << " Elements, organized in " << endl;
  os << num_elem_blk_ << " ElementBlocks, ";
  os << num_node_sets_ << " NodeSets, ";
  os << num_side_sets_ << " SideSets ";
  os << endl << endl;
  if (verbose){
    os << "ElementBlocks" << endl;
    map<int,ElementBlock>::const_iterator it;
    map<int,ElementBlock> eleBlocks = GetElementBlocks();
    for (it=eleBlocks.begin(); it != eleBlocks.end(); it++){
      os << it->first << ": ";
      it->second.Print(os);
    }
    os << endl << "NodeSets" << endl;
    map<int,NodeSet>::const_iterator it2;
    map<int,NodeSet> nodeSets = GetNodeSets();
    for (it2=nodeSets.begin(); it2 != nodeSets.end(); it2++){
      os << "NodeSet " << it2->first << ": ";
      it2->second.Print(os);
    }
    os << endl << "SideSets" << endl;
    map<int,SideSet>::const_iterator it3;
    map<int,SideSet> sideSets = GetSideSets();
    for (it3=sideSets.begin(); it3 != sideSets.end(); it3++){
      os << "SideSet " << it3->first << ": ";
      it3->second.Print(os);
    }
  }
}

/*----------------------------------------------------------------------*
 |  Write Mesh into exodus file                                maf 01/08|
 *----------------------------------------------------------------------*/
void Mesh::WriteMesh(string newexofilename)
{
  //string newexofile(newexofilename);
  //newexofile += ".exo";
  const char *newexofilechar;
  newexofilechar = newexofilename.c_str();
  
  int CPU_word_size, IO_word_size, exoid, error;
  CPU_word_size = sizeof(float); /* use float or double*/
  IO_word_size = 8; /* store variables as doubles */
  
  /* create EXODUS II file */
  exoid = ex_create (newexofilechar, /* filename path */
                     EX_CLOBBER,     /* create mode */
                     &CPU_word_size, /* CPU float word size in bytes */
                     &IO_word_size); /* I/O float word size in bytes */
  
//  // prefer strings
//  string title = "New Exodus Mesh";
//  const char *titlechar;
//  titlechar = title.c_str();
  
  /* initialize file with parameters */
  error = ex_put_init (exoid, title_, num_dim_, num_nodes_, num_elem_,num_elem_blk_, num_node_sets_, num_side_sets_);
  if (error!=0) dserror("error in exfile init");
  
  /* Write QA record based on original exofile */
  int num_qa_rec;
  char *qa_record[MAX_STR_LENGTH+1][4];
  char *cdum;
  float fdum;
  /* read QA records */
  ex_inquire (exoid_, EX_INQ_QA, &num_qa_rec, &fdum, cdum);/* write QA records */
  for (int i=0; i<num_qa_rec; i++)
    for (int j=0; j<4; j++) qa_record[i][j] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
  error = ex_get_qa (exoid_, qa_record);
  error = ex_put_qa (exoid, num_qa_rec, qa_record);
  free(qa_record);
  
  // Write coord names based on original exofile
  char *coord_names[3];
  for (int i=0; i<num_dim_; i++) coord_names[i] = (char *) calloc ((MAX_STR_LENGTH+1), sizeof(char));
  error = ex_get_coord_names (exoid_, coord_names);
  error = ex_put_coord_names (exoid, coord_names);
  free(coord_names);
  
  // Write nodal coordinates
  float x[num_nodes_];
  float y[num_nodes_];
  float z[num_nodes_];
  for (int i = 0; i < num_nodes_; ++i) {
    RCP<PreNode> actnode = GetNode(i);
    x[i] = actnode->X()[0];
    y[i] = actnode->X()[1];
    z[i] = actnode->X()[2];
  }
  error = ex_put_coord (exoid, x, y, z);
  
  // Write mesh entities
  int eb_counter = 0;
  //int ns_counter = 0;
  int ss_counter = 0;
  for (int i = 0; i < num_entities_; ++i) {
    RCP<Entity> actEntity = GetEntity(i);
    switch (actEntity->GetEntityType()){
    case Entity::elem_blk:{
      eb_counter += 1;
      error = ex_put_elem_block(exoid,
                                eb_counter,
                                actEntity->GetElementType().c_str(),
                                actEntity->GetNumEle(),
                                actEntity->GetNumNodpElem(),
                                actEntity->GetNumAttr());
      if (error!=0) dserror("error writing element block");
      
      int* eleconn[num_nodes_];
      //eleconn.resize(num_nodes_);
      for (int i = 0; i < num_nodes_; ++i) {
        (*eleconn)[i] = (*actEntity->Cont())[i];
      }
      //error = ex_put_elem_conn (exoid, eb_counter, actEntity->Cont());
      error = ex_put_elem_conn (exoid, eb_counter, eleconn[0]);
      if (error!=0) dserror("error writing element conn");
      
      // write block name
      string bname = actEntity->GetEntityName();
      const char* blockname = bname.c_str();
      error = ex_put_name (exoid, EX_ELEM_BLOCK, eb_counter, blockname);

      break;
    }
    case Entity::node_set:{
//      ns_counter += 1;
//      //actEntity->Print(cout);
//      cout << "NS" << endl;
//      error = ex_put_node_set_param (exoid,
//                                     ns_counter,
//                                     actEntity->GetNumNodes(),
//                                     0);
//      if (error!=0) dserror("error writing node set");
//      error = ex_put_node_set (exoid, ns_counter, actEntity->Cont());
      break;
    }
    case Entity::side_set:{
      ss_counter += 1;
      cout << "SS" << endl;
      break;
    }
    default: dserror("unknown entity type");
    }
    
    

  }
  
  
  // close file
  error = ex_close (exoid);
  if (error!=0) dserror("error closing exodus file");
  
  
  
  
}



/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
Entity::Entity(int exoid, int entityID, int typeID, EntityType entitytype)
{
  entityID_ = entityID;
  int error;
  switch (entitytype){
  case elem_blk:
  {
    entitytype_ = Entity::elem_blk;
    entity_type_ = "ElementBlock";
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_elem_block (exoid, typeID, mychar, &num_el_in_blk_, &num_nod_per_elem_, &num_attr_);
    // prefer string to store element type
    string ele_type(mychar,int (MAX_STR_LENGTH));
    elem_type_ = ele_type;
    
    // get ElementBlock name
    error = ex_get_name (exoid, EX_ELEM_BLOCK, typeID, mychar);
    // prefer string to store name
    string blockname(mychar,int (MAX_STR_LENGTH));
    entity_name_ = blockname;
    
    // number of nodes for blocks is:
    num_nodes_ = num_nod_per_elem_ * num_el_in_blk_;
    // property name (not yet supported for ElementBlocks)
    entity_prop_name_ = "None";
    
    // get element connectivity
    int connect[num_nodes_];
    (*entity_cont_).resize(num_nodes_);
    error = ex_get_elem_conn(exoid,typeID,connect);
    if (error != 0) dserror("exo error returned");
    for (int i = 0; i < num_nodes_; ++i) {
      //!TODO: Sort start id from nodes: 0 or 1
      (*entity_cont_)[i] = connect[i];
    }
    break;
  }
  case node_set:
  {
    entitytype_ = Entity::node_set;
    entity_type_ = "NodeSet";
    int num_nodes_in_set=0, num_df_in_set=0;
    error = ex_get_node_set_param (exoid, typeID, &num_nodes_in_set,&num_df_in_set);
    num_nodes_ = num_nodes_in_set;
    
    // get NodeSet name
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_name (exoid, EX_NODE_SET, typeID, mychar);
    // prefer string to store name
    string nodesetname(mychar, int(MAX_STR_LENGTH));
    entity_name_ = nodesetname;

    // get nodes in node set
    //int *node_set_node_list;
    //(*entity_cont_).resize(0);
//    entity_cont_.resize(num_nodes_);
//    error = ex_get_node_set (exoid, typeID, node_set_node_list);
//    if (error != 0) dserror("error reading node set");
//    for (int i = 0; i < num_nodes_; ++i) {
//      //!TODO: Sort start id from nodes: 0 or 1
//      entity_cont_[i] = node_set_node_list[i];
//    }

    // set other variables to default NodeSet
    elem_type_ = "No Element Type";
    num_nod_per_elem_ = 0;
    num_el_in_blk_ = 0;
    entity_prop_name_ = "None";
    break;
  }
  case side_set:
  {
    entitytype_ = Entity::side_set;
    entity_type_ = "SideSet";
    // get SideSet name
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_name (exoid, EX_SIDE_SET, typeID, mychar);
    // prefer string to store name
    string sidesetname(mychar, int(MAX_STR_LENGTH));
    entity_name_ = sidesetname;
    
    // set other variables to default
    elem_type_ = "Side Sets not yet supported";
    entity_prop_name_ = "None";
    (*entity_cont_).resize(0);
    break;
  }
  default: cout << "EntityType not valid" << endl;
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
Entity::~Entity()
{
  return;
}

void Entity::SetPropertyName(string propname)
{
  entity_prop_name_ = propname;
  return;
}

void Entity::Print(ostream& os) const
{
  // do not remove the .c_str() since they are needed for printing into a file stream
  // entityID_ is raised by 1 to match matr numbering of existing exofilter
  os << "Entity " << entityID_+1<< " is of type " << entity_type_;
  os << " is named " << entity_name_.c_str() << endl;
  os << "with " << num_nodes_ << " Nodes in cloud" << endl;
  os << "Additional Info: " << endl;
  os << "Property Name: " << entity_prop_name_.c_str() << endl;
  os << "Element Type: " << elem_type_.c_str() << endl;
  os << "Num Attr: " << num_attr_;
  os << ", ele per block: " << num_el_in_blk_ << ", num per ele: " << num_nod_per_elem_ << endl << endl;
  return;
}

ElementBlock::ElementBlock(ElementBlock::Shape Distype, map<int,vector<int> > &eleconn, string name)
{
  distype_ = Distype;
  eleconn_ = eleconn;
  name_ = name;
  return;
}
ElementBlock::~ElementBlock()
{
  return;
}
void ElementBlock::Print(ostream& os, bool verbose) const
{
  os << "Element Block, named: " << name_ << endl
  << "of Shape: " << ShapeToString(distype_) << endl;
  if (verbose){
    map<int,vector<int> >::const_iterator it;
    for (it=eleconn_.begin(); it != eleconn_.end(); it++){
      os << "Ele " << it->first << ": ";
//      for (vector<int>::const_iterator var = it->second.begin(); var < it->second.end(); ++var) {
//        os << it->second.at(var) << ","; 
//      }
      for (int i=0; i < signed(it->second.size()); i++ ){
        os << it->second.at(i) << ",";
      }
      os << endl;
    }
  }
}


NodeSet::NodeSet(set<int> nodeids, string name, string propname)
{
  nodeids_ = nodeids;
  name_ = name;
  propname_ = propname;
}

NodeSet::~NodeSet()
{
  return;
}

void NodeSet::Print(ostream& os, bool verbose) const
{
  os << "Node Set, named: " << name_ << endl
  << "Property Name: " << propname_ << endl;
  if (verbose){
    //!TODO:
    os << "NodeSet verbose" << endl;
  }
}


SideSet::SideSet(string name)
{
  name_ = name;
}

SideSet::~SideSet()
{
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
PreElement::PreElement(int id, ShapeType shape)
{
  id_ =id;
  shape_ = shape;
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
PreElement::~PreElement()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/08|
 *----------------------------------------------------------------------*/
PreElement::PreElement(const PreElement& old)
{
  id_ = old.id_;
  shape_ = old.shape_;
  nodeid_ = old.nodeid_;
  node_ = old.node_;
  return;
}

/*----------------------------------------------------------------------*
 |  print element (public)                                     maf 01/08|
 *----------------------------------------------------------------------*/
void PreElement::Print(ostream& os) const
{
  os << "Id: " << Id() << " Shape: " ;
  switch (Shape()) {
    case quad4:
      os << "quad4 ";
      break;
    case tri3:
      os << "tri3 ";
      break;
    default:
      dserror("unknown PreElement shape");
      break;
  } 
  const int nnode = NumNode();
  const int* nodes = NodeIds();
  if (nnode)
  {
    os << " Nodes ";
    for (int i=0; i<nnode; ++i) os << setw(10) << nodes[i] << " ";
  }
  os << endl;
}

/*----------------------------------------------------------------------*
 |  set node numbers to element (public)                       maf 01/08|
 *----------------------------------------------------------------------*/
void PreElement::SetNodeIds(const int nnode, const int* nodes)
{
  nodeid_.resize(nnode);
  for (int i=0; i<nnode; ++i) nodeid_[i] = nodes[i];
  node_.resize(0);
  return;
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
PreNode::PreNode(int id, const double* coords)
{
  id_ =id;
  for (int i=0; i<3; ++i) x_[i] = coords[i];
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 01/08|
 *----------------------------------------------------------------------*/
PreNode::~PreNode()
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                         maf 01/08|
 *----------------------------------------------------------------------*/
PreNode::PreNode(const PreNode& old)
{
  id_ = old.id_;
  for (int i=0; i<3; ++i) x_[i] = old.x_[i];
  return;
}

/*----------------------------------------------------------------------*
 |  print node    (public)                                     maf 01/08|
 *----------------------------------------------------------------------*/
void PreNode::Print(ostream& os) const
{
  os << "Id: " << Id() << " Coords: "
     << X()[0] << "," << X()[1] << "," << X()[2] << endl;
}


#endif
