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
#ifndef PRE_EXODUS_READER_H_
#define PRE_EXODUS_READER_H_

#include <string>
#include <iostream>

extern "C" /* stuff which is c and is accessed from c++ */
{
#include <exodusII.h>
}

using namespace std;

class mesh
{
public:
  //! constructor
  mesh(string exofilename);
  
  //! destructor
  virtual ~mesh();
  
private:
  //! number of dimensions
  int num_dim_;
  //! number of nodes
  int num_nodes_;
  //! number of elements
  int num_elem_;
  //! number of element blocks
  int num_elem_blk_;
  //! number of node sets
  int num_node_sets_;
  //! number of side sets
  int num_side_sets_;
  //! exoid
  int exoid_;
  //! title
  char title_[MAX_LINE_LENGTH+1];
  
};



#endif /*PRE_EXODUS_READER_H_*/
