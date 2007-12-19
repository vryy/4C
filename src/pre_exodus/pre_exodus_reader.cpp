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

#include "pre_exodus_reader.H"

using namespace std;

mesh::mesh(string exofilename)
{
  int error;
  int CPU_word_size,IO_word_size, exoid;
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
  
  cout << "NumNodes:" << num_nodes_ << endl;
  cout << "NumEle:" << num_elem_ << endl;
  cout << "NumEleBlk:" << num_elem_blk_ << endl;
  
  char** eb_names = new char*[num_elem_blk_];
  for (int i=0; i<num_elem_blk_; ++i)
  {
     eb_names[i] = new char[MAX_STR_LENGTH+1];
  }
  
  error = ex_get_names (exoid_,EX_ELEM_BLOCK,eb_names);
  
  for (int i=0; i<num_elem_blk_; ++i)   
  {
    printf("%3d   %s\n",i,eb_names[i]);
  }
  
  for (int i = 1; i <= num_elem_blk_; ++i) {
    char* elem_type = new char[MAX_STR_LENGTH+1];
    int num_el_in_blk, num_nod_per_elem, num_attr;
    error = ex_get_elem_block (exoid_, i, elem_type, &num_el_in_blk, &num_nod_per_elem, &num_attr);
    printf("%3s\n",elem_type);
    cout << num_el_in_blk << endl;
  }
  
  // close exofile
  error = ex_close(exoid_);
  
  return;
}

mesh::~mesh()
{
  return;
}

