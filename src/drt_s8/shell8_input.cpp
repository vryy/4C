/*!----------------------------------------------------------------------
\file shell8_input.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SHELL8
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

extern "C" 
{
#include "../headers/standardtypes.h"
/*!----------------------------------------------------------------------
  \brief file pointers

  <pre>                                                         m.gee 8/00
  This structure struct _FILES allfiles is defined in input_control_global.c
  and the type is in standardtypes.h
  It holds all file pointers and some variables needed for the FRSYSTEM
  </pre>
 *----------------------------------------------------------------------*/
extern struct _FILES  allfiles;
}
#include "shell8.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                              mwgee 11/06|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Shell8::ReadElement()
{
  // read element's nodes
  int ierr=0;
  int nnode=0;
  int nodes[9]; 
  frchk("QUAD4",&ierr);
  if (ierr==1)
  {
    nnode = 4;
    frint_n("QUAD4",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("QUAD8",&ierr);
  if (ierr==1)
  {
    nnode = 8;
    frint_n("QUAD8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("QUAD9",&ierr);
  if (ierr==1)
  {
    nnode = 9;
    frint_n("QUAD9",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI3",&ierr);
  if (ierr==1)
  {
    nnode = 3;
    frint_n("TRI3",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  frchk("TRI6",&ierr);
  if (ierr==1)
  {
    nnode = 6;
    frint_n("TRI6",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;
  
  SetNodeIds(nnode,nodes);
  
  // read number of material model
  material_ = 0;
  frint("MAT",&material_,&ierr);
  if (ierr!=1) dserror("Reading of SHELL8 element failed");
  
  // read shell thickness
  thickness_ = 1.0;
  frdouble("THICK",&thickness_,&ierr);
  if (ierr!=1) dserror("Reading of SHELL8 element failed");

  // read gaussian points
  frint_n("GP",ngp_,3,&ierr);
  if (ierr!=1) dserror("Reading of SHELL8 element failed");
    
  // read gaussian points for triangle element
  frint("GP_TRI",&ngptri_,&ierr);
  if (ierr!=1) dserror("Reading of SHELL8 element failed");
  
  // read local or global forces
  char buffer[50];
  frchar("FORCES",buffer,&ierr);
  if (ierr)
  {
   if      (strncmp(buffer,"XYZ",3)==0)       forcetype_ = s8_xyz;
   else if (strncmp(buffer,"RST",3)==0)       forcetype_ = s8_rst;
   else if (strncmp(buffer,"RST_ortho",9)==0) forcetype_ = s8_rst_ortho;
   else dserror("Reading of SHELL8 element failed");
  }
  
  // read EAS parameters
  for (int i=0; i<5; ++i) eas_[i] = 0;
  char* colpointer = strstr(allfiles.actplace,"EAS");
  colpointer+=3;
  colpointer = strpbrk(colpointer,"Nn");
  ierr = sscanf(colpointer," %s ",buffer);  
  if (ierr!=1) dserror("Reading of shell8 eas failed");
  if (strncmp(buffer,"none",4)==0)  eas_[0]=0;
  if (strncmp(buffer,"N4_1",4)==0)  eas_[0]=1;
  if (strncmp(buffer,"N4_2",4)==0)  eas_[0]=2;
  if (strncmp(buffer,"N4_3",4)==0)  eas_[0]=3;
  if (strncmp(buffer,"N4_4",4)==0)  eas_[0]=4;
  if (strncmp(buffer,"N4_5",4)==0)  eas_[0]=5;
  if (strncmp(buffer,"N4_7",4)==0)  eas_[0]=7;
  if (strncmp(buffer,"N9_7",4)==0)  eas_[0]=7;
  if (strncmp(buffer,"N9_9",4)==0)  eas_[0]=9;
  if (strncmp(buffer,"N9_11",4)==0) eas_[0]=11;
  colpointer += strlen(buffer);
  
  colpointer = strpbrk(colpointer,"Nn");
  ierr = sscanf(colpointer," %s ",buffer);
  if (ierr!=1) dserror("Reading of shell8 eas failed");
  if (strncmp(buffer,"none",4)==0)  eas_[1]=0;
  if (strncmp(buffer,"N4_4",4)==0)  eas_[1]=4;
  if (strncmp(buffer,"N4_5",4)==0)  eas_[1]=5;
  if (strncmp(buffer,"N4_6",4)==0)  eas_[1]=6;
  if (strncmp(buffer,"N4_7",4)==0)  eas_[1]=7;
  if (strncmp(buffer,"N9_9",4)==0)  eas_[1]=9;
  if (strncmp(buffer,"N9_11",4)==0) eas_[1]=11;
  colpointer += strlen(buffer);
    
  colpointer = strpbrk(colpointer,"Nn");
  ierr = sscanf(colpointer," %s ",buffer);
  if (ierr!=1) dserror("Reading of shell8 eas failed");
  if (strncmp(buffer,"none",4)==0)  eas_[2]=0;
  if (strncmp(buffer,"N_1",4)==0)   eas_[2]=1;
  if (strncmp(buffer,"N_3",4)==0)   eas_[2]=3;
  if (strncmp(buffer,"N_4",4)==0)   eas_[2]=4;
  if (strncmp(buffer,"N_6",4)==0)   eas_[2]=6;
  if (strncmp(buffer,"N_8",4)==0)   eas_[2]=8;
  if (strncmp(buffer,"N_9",4)==0)   eas_[2]=9;
  colpointer += strlen(buffer);
  
  colpointer = strpbrk(colpointer,"Nn");
  ierr = sscanf(colpointer," %s ",buffer);
  if (ierr!=1) dserror("Reading of shell8 eas failed");
  if (strncmp(buffer,"none",4)==0)  eas_[3]=0;
  if (strncmp(buffer,"N4_2",4)==0)  eas_[3]=2;
  if (strncmp(buffer,"N4_4",4)==0)  eas_[3]=4;
  if (strncmp(buffer,"N9_2",4)==0)  eas_[3]=2;
  if (strncmp(buffer,"N9_4",4)==0)  eas_[3]=4;
  if (strncmp(buffer,"N9_6",4)==0)  eas_[3]=6;
  colpointer += strlen(buffer);                   
  
  colpointer = strpbrk(colpointer,"Nn");
  ierr = sscanf(colpointer," %s ",buffer);
  if (ierr!=1) dserror("Reading of shell8 eas failed");
  if (strncmp(buffer,"none",4)==0)  eas_[4]=0;
  if (strncmp(buffer,"N4_2",4)==0)  eas_[4]=2;
  if (strncmp(buffer,"N4_4",4)==0)  eas_[4]=4;
  if (strncmp(buffer,"N9_2",4)==0)  eas_[4]=2;
  if (strncmp(buffer,"N9_4",4)==0)  eas_[4]=4;
  if (strncmp(buffer,"N9_6",4)==0)  eas_[4]=6;

  // count no. eas parameters
  nhyb_ = 0;
  for (int i=0; i<5; ++i) nhyb_ += eas_[i];
  // create arrays alfa, Dtildinv, Lt, Rtild in data_
  vector<double> alfa(nhyb_);
  vector<double> Rtild(nhyb_);
  for (int i=0; i<nhyb_; ++i) 
  {
    alfa[i] = 0.0;
    Rtild[i] = 0.0;
  }
  Epetra_SerialDenseMatrix Dtildinv;
  Epetra_SerialDenseMatrix Lt;
  Dtildinv.Shape(nhyb_,nhyb_);
  Lt.Shape(nhyb_,nnode*6);
  data_.Add("alfa",alfa);
  data_.Add("Rtild",Rtild);
  data_.Add("Dtildinv",Dtildinv);                                                  
  data_.Add("Lt",Lt);                                                  
    
  // read ANS
  ans_ = 0;
  frchar("ANS",buffer,&ierr);
  if (ierr!=1) dserror("reading of shell8 ans failed");
  if (strncmp(buffer,"none",4)==0) ans_=0;
  if (strncmp(buffer,"Q",4)==0)    ans_=1;
  if (strncmp(buffer,"T",4)==0)    ans_=2;
  if (strncmp(buffer,"QT",4)==0)   ans_=3;
  if (strncmp(buffer,"TQ",4)==0)   ans_=3;
      
  // read SDC
  sdc_ = 1.0;  
  frdouble("SDC",&sdc_,&ierr);
  if (ierr!=1) dserror("Reading of shell8 sdc failed");
    
  return true;
} // Shell8::ReadElement()
















#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SHELL8
