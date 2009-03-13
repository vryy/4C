/*----------------------------------------------------------------------*/
/*!
\file so_sh8p8_input.cpp
\brief

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
/* macros */
#ifdef D_SOLID3
#ifdef CCADISCRET

/*----------------------------------------------------------------------*/
/* headers */
// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif

#include "so_sh8p8.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                         |
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_sh8p8::ReadElement()
{
  // read element's nodes
  int ierr = 0;
  const int nnode = NUMNOD_SOSH8P8;
  int nodes[NUMNOD_SOSH8P8];
  frchk("SOLIDSH8P8",&ierr);
  if (ierr==1)
  {
    frint_n("HEX8",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDSH8P8 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i){
    nodes[i]--;
  }

  SetNodeIds(nnode,nodes);


  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_SH8P8 element material failed");
  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod){
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(Material().get());
    remo->Setup(NUMGPT_SOH8, this->Id());
  } else if (Material()->MaterialType() == INPAR::MAT::m_anisotropic_balzani){
    MAT::AnisotropicBalzani* balz = static_cast <MAT::AnisotropicBalzani*>(Material().get());
    balz->Setup();
  } else if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic){
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_SOH8);
  } else if (Material()->MaterialType() == INPAR::MAT::m_visconeohooke){
    MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(Material().get());
    visco->Setup(NUMGPT_SOH8);
  }

  // read possible gaussian points, obsolete for computation
  {
    int ngp[3];
    frint_n("GP",ngp,3,&ierr);
    if (ierr==1) for (int i=0; i<3; ++i) if (ngp[i]!=2) dserror("Only 2 GP for So_SH8P8");
  }

  // we expect kintype to be total lagrangian
  kintype_ = soh8_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
    // geometrically linear
    if      (strncmp(buffer,"Geolin",6)==0)    kintype_ = soh8_geolin;
    // geometrically non-linear with Total Lagrangean approach
    else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = soh8_totlag;
    // geometrically non-linear with Updated Lagrangean approach
    else if (strncmp(buffer,"Updlag",6)==0)
    {
      kintype_ = soh8_updlag;
      dserror("Updated Lagrange for SO_SH8P8 is not implemented!");
    }
    else dserror("Reading of SO_SH8P8 element failed");
  }

  // set EAS technology flag
  eastype_ = soh8_easnone;
  neas_ = 0;

  // read global coordinate of shell-thickness direction
  thickdir_ = autoj;           // default: auto by Jacobian
  nodes_rearranged_ = false;
  frchar("THICKDIR",buffer,&ierr);
  if (ierr)
  {
    // global X
    if      (strncmp(buffer,"xdir",4)==0)    thickdir_ = globx;
    // global Y
    else if (strncmp(buffer,"ydir",4)==0)    thickdir_ = globy;
    // global Z
    else if (strncmp(buffer,"zdir",4)==0)    thickdir_ = globz;
    // find automatically through Jacobian of Xrefe
    else if (strncmp(buffer,"auto",4)==0)    thickdir_ = autoj;
    // no noderearrangement
    else if (strncmp(buffer,"none",4)==0){
      thickdir_ = none;
      nodes_rearranged_ = true;
    }
    else dserror("Reading of SO_SH8P8 thickness direction failed");
  }

  // stabilisation
  stab_ = stab_affine;
  frchar("STAB",buffer,&ierr);
  if (ierr)
  {
    if (strncmp(buffer,"Aff",3)==0)
      stab_ = stab_affine;
    else if (strncmp(buffer,"NonAff",6)==0)
      stab_ = stab_nonaffine;
    else if (strncmp(buffer,"PureDisp",8)==0)
      stab_ = stab_puredisp;
    else
      dserror("Reading of SO_SH8P8 stabilisation failed");
  }

  // ANS
  ans_ = ans_lateral;
  frchar("ANS",buffer,&ierr);
  if (ierr)
  {
    if (strncmp(buffer,"Lateral",10)==0)
      ans_ = ans_lateral;
    else if (strncmp(buffer,"None",4)==0)
      ans_ = ans_none;
    else
      dserror("Reading of SO_SH8P8 ANS type failed");
  }

  return true;
} // So_sh8p8::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
