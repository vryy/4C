/*----------------------------------------------------------------------*/
/*!
\file drt_condition_utils_baci.cpp

\brief

<pre>
Maintainer: Christian Roth
            roth@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15255
</pre>
*/
/*----------------------------------------------------------------------*/

#include <map>
#include <string>

#include "drt_condition_utils_baci.H"
#include "drt_globalproblem.H"

#include "../drt_io/io_control.H"


/*-----------------------------------------------------------------------*
 * Writes boundary surfaces of a volumetrically coupled problem to file  *
 * 'boundarysurfaces.log' storing the condition-Id as surface-Id. For    *
 * visualisation in gmsh and checking for tetrahedra whose four surfaces *
 * are wrongly contained in the boundary surface of the volumetric       *
 * coupling this file can be used.                         (croth 01/15) *
 *-----------------------------------------------------------------------*/
void DRT::UTILS::WriteBoundarySurfacesVolumeCoupling(
    std::map< std::vector<int>, Teuchos::RCP<DRT::Element> > surfmap,
    int condID,
    int numproc,
    int mypid)
{
  if (numproc==1)
  {
    //Get output prefix
    std::string outputprefix = DRT::Problem::Instance()->OutputControlFile()->NewOutputFileName();
    //Create boundary surface file
    std::ostringstream sf;
    sf << outputprefix << "_boundarysurfaces.log";
    std::string boundarysurffilename;
    boundarysurffilename = sf.str();

    std::ofstream myfile;
    myfile.open (boundarysurffilename.c_str(), std::ios_base::app);
    myfile << "Surfaces in Surfmap for Coupling Condition No. " << condID+1 << ":\n";
    myfile << "Format: [Node1, Node2, Node3, CondID] \n";
    for(std::map< std::vector<int>, Teuchos::RCP<DRT::Element> >::const_iterator iterat=surfmap.begin(); iterat!=surfmap.end(); ++iterat)
    {
      myfile << iterat->first[0] << " " << iterat->first[1] << " " << iterat->first[2] << " " << condID+1 <<"\n";
    }
    myfile << "End \n";
    myfile.close();
  }
  else if (mypid==0)
  {
    std::cout << " No 'boundarysurfaces.log' written as number of procs = "<< numproc <<" is bigger than 1." << std::endl;
  }
}

