/*!----------------------------------------------------------------------
\file beam2_input.cpp
\brief

<pre>
Maintainer: Christian Cyron
            cyron@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef D_BEAM2R
#ifdef CCADISCRET

#include "beam2r.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_lib/drt_utils.H"



/*----------------------------------------------------------------------*
 |  read element input (public)                              cyron 01/08|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::Beam2r::ReadElement()
{
	/*
	\The element is capable of using higher order functions from linear to quartic. Please make sure you put the nodes in the right order
	\in the input file.
	\LIN2    1---2
	\LIN3	 1---3---2
	\LIN4 	 1---4---2---3
	\LIN5	 1---5---2---3---4
	*/
	
  //function receives only one line of the input .dat file
  int ierr=0;  //error-variable used throughout the function
 
  //typedef for later conversion of string to distype
  //note: GID gives LINX in the .dat file while pre_exodus gives LINEX
  typedef map<string, DiscretizationType> Beam2rDisType;
  Beam2rDisType beam2rdistype;
  beam2rdistype["LIN2"]     = line2;
  beam2rdistype["LINE2"]    = line2;
  beam2rdistype["LIN3"]     = line3;
  beam2rdistype["LINE3"]    = line3;
  beam2rdistype["LIN4"]     = line4;
  beam2rdistype["LINE4"]    = line4;
  beam2rdistype["LIN5"]     = line5;
  beam2rdistype["LINE5"]    = line5;
  
  
  DiscretizationType distype = dis_none;
  
  //Iterator goes through all possibilities
  Beam2rDisType::const_iterator iter;
      for( iter = beam2rdistype.begin(); iter != beam2rdistype.end(); iter++ )
      {
          const string eletext = iter->first;
          frchk(eletext.c_str(), &ierr);
          if (ierr == 1)
          {
              //Get DiscretizationType
        	  distype = beam2rdistype[eletext];
              //Get Number of Nodes of DiscretizationType
              int nnode = DRT::UTILS::getNumberOfElementNodes(distype);
              
              //Set the applied Gaussrule ( It can be proven that we need 1 GP less than nodes to integrate exact )
              //note: we use a static cast for the enumeration here cf. Practical C++ Programming p.185
              gaussrule_ = static_cast<enum DRT::UTILS::GaussRule1D>(nnode-1);
              //Get an array for the global node numbers
              int nodes[nnode];
              //Read global node numbers
              frint_n(eletext.c_str(), nodes, nnode, &ierr);
              
              dsassert(ierr==1, "Reading of ELEMENT Topology failed\n");
              
              // reduce global node numbers by one because BACI nodes begin with 0 and inputfile nodes begin with 1
              for (int i=0; i<nnode; ++i) nodes[i]--;

              SetNodeIds(nnode,nodes); // has to be executed in here because of the local scope of nodes
              break;
          }
      }

  // read material parameters using structure _MATERIAL which is defined by inclusion of
  // "../drt_lib/drt_timecurve.H"
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of Beam2r element failed");
  SetMaterial(material);

  // read beam cross section
  crosssec_ = 1.0;
  frdouble("CROSS",&crosssec_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2r element failed");

   // read beam cross section including shear correction factor
  double shear_correction = 0.0;
  frdouble("SHEARCORR",&shear_correction,&ierr);
  if (ierr!=1) dserror("Reading of Beam2r element failed");
  crosssecshear_ = crosssec_ * shear_correction;

  // read beam moment of inertia of area
  mominer_ = 1.0;
  frdouble("INERMOM",&mominer_,&ierr);
  if (ierr!=1) dserror("Reading of Beam2r element failed");
  
   
  return true;
} // Beam2r::ReadElement()




#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_BEAM2R
