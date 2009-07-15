/*!----------------------------------------------------------------------
\file so_shw6_input.cpp
\brief

<pre>
Maintainer: Moritz Frenzel
            frenzel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15240
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "so_shw6.H" //**
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_shw6::ReadElement()
{
  // read element's nodes
  int ierr=0;
  const int nnode=6;
  int nodes[6];
  frchk("SOLIDSHW6",&ierr);
  if (ierr==1)
  {
    frint_n("WEDGE6",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDSHW6 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i){
    nodes[i]--;
  }

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SOLIDSHW6 element material failed");
  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod){
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(Material().get());
    remo->Setup(NUMGPT_WEG6, this->Id());
  } else if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic){
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_WEG6);
  }

  // we expect kintype to be total lagrangian
  kintype_ = sow6_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically linear
   if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = sow6_geolin;
   // geometrically non-linear with Total Lagrangean approach
   else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = sow6_totlag;
   // geometrically non-linear with Updated Lagrangean approach
   else if (strncmp(buffer,"Updlag",6)==0)
   {
       kintype_ = sow6_updlag;
       dserror("Updated Lagrange for SOLIDSHW6 is not implemented!");
   }
   else dserror("Reading of SOLIDSHW6 element failed");
  }

  // read EAS technology flag
  frchar("EAS",buffer,&ierr);
  if (ierr){
    // full sohw6 EAS technology
    if      (strncmp(buffer,"soshw6",6)==0){
      eastype_ = soshw6_easpoisthick; // EAS to allow linear thickness strain
      neas_ = soshw6_easpoisthick;    // number of eas parameters
      soshw6_easinit();
    }
    // no EAS technology
    else if (strncmp(buffer,"none",4)==0){
      eastype_ = soshw6_easnone;
      cout << "Warning: Solid-Shell Wegde6 without EAS" << endl;
    }
    else dserror("Reading of SO_SH8 EAS technology failed");
  }
  else
  {
    eastype_ = soshw6_easpoisthick;     // default: EAS to allow linear thickness strain
    neas_ = 1;                          // number of eas parameters
    soshw6_easinit();
  }

  // check for automatically align material space optimally with parameter space
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;
  frchar("OPTORDER",buffer,&ierr);
  if (ierr)
    if (strncmp(buffer,"true",4)==0)
      optimal_parameterspace_map_ = true;



  return true;
} // So_weg6::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
