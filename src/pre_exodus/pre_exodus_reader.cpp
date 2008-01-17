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
#ifdef EXODUS
#include "pre_exodus_reader.H"

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

  // num_entities_ are all ElementBlocks, NodeSets, and SideSets together
  num_entities_ = num_elem_blk_ + num_node_sets_ + num_side_sets_;

  // allocate sufficient memory for the entities
  myEntities_.reserve(num_entities_);
  
  // entitycounter counts all ElementBlocks, NodeSets, and SideSets together
  int entitycounter = 0;
  
  // get all ElementBlocks
  int epropID[num_elem_blk_];
  error = ex_get_prop_array(exoid_, EX_ELEM_BLOCK, "ID", epropID);
  for (int i = 0; i < num_elem_blk_; ++i) 
  {
	  //RCP<Entity> myentity = rcp(new Entity(exoid_,entitycounter,epropID[i],Entity::elem_blk));
	  myEntities_[entitycounter] = rcp(new Entity(exoid_,entitycounter,epropID[i],Entity::elem_blk));
	  //myEntities_[entitycounter]->Print(cout);
	  entitycounter++;
  }
  
  // get all NodeSets
  int npropID[num_node_sets_];
  error = ex_get_prop_array(exoid_, EX_NODE_SET, "ID", npropID);
  for (int i = 0; i < num_node_sets_; ++i) {
    myEntities_[entitycounter] = rcp(new Entity(exoid_,entitycounter,npropID[i],Entity::node_set));
    //myEntities_[entitycounter]->Print(cout);
    //Entity myentity(exoid_,entitycounter,npropID[i],Entity::node_set);
    entitycounter++;
  }
  
  // get all SideSets
  int spropID[num_side_sets_];
  error = ex_get_prop_array(exoid_, EX_SIDE_SET, "ID", spropID);
  for (int i = 0; i < num_side_sets_; ++i) {
    myEntities_[entitycounter] = rcp(new Entity(exoid_,entitycounter,spropID[i],Entity::side_set));
    //myEntities_[entitycounter]->Print(cout);
    //Entity myentity(exoid_,entitycounter,spropID[i],Entity::side_set);
    entitycounter++;
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
//  for (int i=0; i<num_props; ++i)   
//  {
//    printf("%3d   %s\n",i,prop_names[i]);
//  }
  if ((num_props-1) == num_node_sets_){
    for (int i = 1; i < num_props; ++i) {
      string propname(prop_names[i], int(MAX_STR_LENGTH));
      myEntities_[num_elem_blk_-1+i]->SetPropertyName(propname);
      //myEntities_[num_elem_blk_-1+i]->Print(cout);
    }
  }

  // close exodus II file
  error = ex_close(exoid_);
  if (error < 0)
	  dserror("error while closing exodus II file");

  return;
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                              maf 12/07|
 *----------------------------------------------------------------------*/
Mesh::~Mesh()
{
  return;
}

/*----------------------------------------------------------------------*
 |  Print method (public)                                      maf 12/07|
 *----------------------------------------------------------------------*/
void Mesh::Print(ostream & os) const
{
  os << "Mesh consists of ";
  os << num_nodes_ << " Nodes, ";
  os << num_elem_ << " Elements, organized in ";
  os << num_entities_ << " Entities, divided into " << endl;
  os << num_elem_blk_ << " ElementBlocks, ";
  os << num_node_sets_ << " NodeSets, ";
  os << num_side_sets_ << " SideSets ";
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
    // property name (not yet supported for ElementBlocks)
    entity_prop_name_ = "None";
    break;
  }
  case node_set:
  {
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

    // set other variables to default NodeSet
    elem_type_ = "No Element Type";
    num_nod_per_elem_ = 0;
    num_el_in_blk_ = 0;
    entity_prop_name_ = "None";
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
    entity_prop_name_ = "None";
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


#endif
