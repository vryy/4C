/*----------------------------------------------------------------------*/
/*!
\file pre_exodus_reader.H

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
#ifdef EXODUS
#include "pre_exodus_reader.H"

using namespace std;

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

  cout << "meshfilename: " << exofilename.c_str() << endl;
  
  const char *exofilenamechar;
  exofilenamechar=exofilename.c_str();
  
  /* open EXODUS II files */
  exoid_ = ex_open(exofilenamechar,EX_READ,&CPU_word_size,&IO_word_size,&exoversion);
  if (exoid_<0){ cout <<"Exo-file does not exist"<< endl; exit(1);}
  // print version
  cout<<"Input file uses EXODUS II library version "<<exoversion<<endl;
  
  /* read database parameters */
  error = ex_get_init (exoid_, title_, &num_dim_, &num_nodes_,&num_elem_, &num_elem_blk_, &num_node_sets_, &num_side_sets_);
  
  // num_entities_ are all ElementBlocks, NodeSets, and SideSets together
  num_entities_ = num_elem_blk_ + num_node_sets_ + num_side_sets_;
  
  
  // entityconuter counts all ElementBlocks, NodeSets, and SideSets together
  int entitycounter = 0;
  
  // get all ElementBlocks
  for (int i = 0; i < num_elem_blk_; ++i) {
    Entity myentity(exoid_,entitycounter,i+1,Entity::elem_blk);
    entitycounter++;
    
    myentity.Print(cout);
  }
  
  // get all NodeSets
  for (int i = 0; i < num_node_sets_; ++i) {
    Entity myentity(exoid_,entitycounter,i+1,Entity::node_set);
    entitycounter++;

    myentity.Print(cout);
  }
  
  // get all SideSets
  for (int i = 0; i < num_side_sets_; ++i) {
    Entity myentity(exoid_,entitycounter,i+1,Entity::side_set);
    entitycounter++;
    
    myentity.Print(cout);
  }

  /* Read NodeSet property names
   * They are assigned by ICEM and provide recognition */
  // read number of node set properties
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

  

  

  // close exofile
  error = ex_close(exoid_);
  if (error);
  
  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
Mesh::~Mesh()
{
  return;
}

void Mesh::Print(ostream & os) const
{
  os << "Mesh consists of ";
  os << num_nodes_ << " Nodes, ";
  os << num_elem_ << " Elements, organized in ";
  os << num_entities_ << " Entities, divided into " << endl;
  os << num_elem_blk_ << " ElementBlocks, ";
  os << num_node_sets_ << " NodeSets, ";
  os << num_side_sets_ << " SideSets, ";
  os << endl << endl;
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
    break;
  }
  case node_set:
  {
    entity_type_ = "NodeSet";
    int num_nodes_in_set, num_df_in_set;
    error = ex_get_node_set_param (exoid, typeID, &num_nodes_in_set,&num_df_in_set);
    num_nodes_ = num_nodes_in_set;
    
    // get NodeSet name
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_name (exoid, EX_NODE_SET, typeID, mychar);
    // prefer string to store name
    string nodesetname(mychar, int(MAX_STR_LENGTH));
    entity_name_ = nodesetname;

    // set other variables to default NodeSet
    elem_type_ = "No Element Type";
    num_nod_per_elem_ = 0;
    num_el_in_blk_ = 0;
    break;
  }
  case side_set:
  {
    entity_type_ = "SideSet";
    // get SideSet name
    char mychar[MAX_STR_LENGTH+1];
    error = ex_get_name (exoid, EX_SIDE_SET, typeID, mychar);
    // prefer string to store name
    string sidesetname(mychar, int(MAX_STR_LENGTH));
    entity_name_ = sidesetname;
    
    elem_type_ = "Side Sets not yet supported";
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

void Entity::Print(ostream& os) const
{
  os << "Entity " << entityID_ << " is of type " << entity_type_;
  os << " is named " << entity_name_ << endl;
  os << "with " << num_nodes_ << " Nodes in cloud" << endl;
  os << "Additional Info: " << endl;
  os << "Element Type: " << elem_type_;
  os << ", Num Attr: " << num_attr_;
  os << ", ele per block: " << num_el_in_blk_ << ", num per ele: " << num_nod_per_elem_ << endl << endl;
  return;
}


#endif
