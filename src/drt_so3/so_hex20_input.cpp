/*!----------------------------------------------------------------------
\file so_hex20_input.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/
#ifdef D_SOLID3
#ifdef CCADISCRET

#include "so_hex20.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"

/*----------------------------------------------------------------------*
 |  read element input (public)                                maf 04/07|
 *----------------------------------------------------------------------*/
bool DRT::ELEMENTS::So_hex20::ReadElement()
{
  // read element's nodes
  int ierr=0;
  const int nnode=20;
  int nodes[20];
  frchk("SOLIDH20",&ierr);
  if (ierr==1)
  {
    frint_n("HEX20",nodes,nnode,&ierr);
    if (ierr != 1) dserror("Reading of ELEMENT Topology failed");
  }
  else
  {
    dserror ("Reading of SOLIDH20 failed");
  }
  // reduce node numbers by one
  for (int i=0; i<nnode; ++i) nodes[i]--;

  SetNodeIds(nnode,nodes);

  // read number of material model
  int material = 0;
  frint("MAT",&material,&ierr);
  if (ierr!=1) dserror("Reading of SO_HEX20 element material failed");
  SetMaterial(material);

  // special element-dependent input of material parameters
  if (Material()->MaterialType() == INPAR::MAT::m_artwallremod){
    MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(Material().get());
    remo->Setup(NUMGPT_SOH20, this->Id());
  } else if (Material()->MaterialType() == INPAR::MAT::m_viscoanisotropic){
    MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(Material().get());
    visco->Setup(NUMGPT_SOH20);
  } else if (Material()->MaterialType() == INPAR::MAT::m_visconeohooke){
    MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(Material().get());
    visco->Setup(NUMGPT_SOH20);
  } else if (Material()->MaterialType() == INPAR::MAT::m_charmm){
    MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(Material().get());
    charmm->Setup(data_);
  } else if (Material()->MaterialType() == INPAR::MAT::m_aaaraghavanvorp_damage){
    double strength = 0.0; // section for extracting the element strength 
    frdouble("STRENGTH",&strength,&ierr);
    if (ierr!=1) dserror("Reading of SO_SH8 element strength failed");
    MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(Material().get());
    aaadamage->Setup(NUMGPT_SOH8,strength);
    //aaadamage->Setup(NUMGPT_SOH8);
  }

  // read possible gaussian points, obsolete for computation
  int ngp[3];
  frint_n("GP",ngp,3,&ierr);
  if (ierr==1) for (int i=0; i<3; ++i) if (ngp[i]!=3) dserror("Only version with 3 GP for So_H20 implemented");

  // we expect kintype to be total lagrangian
  kintype_ = soh20_totlag;

  // read kinematic type
  char buffer[50];
  frchar("KINEM",buffer,&ierr);
  if (ierr)
  {
   // geometrically linear
   if      (strncmp(buffer,"Geolin",6)==0)
   {
     kintype_ = soh20_geolin;
     dserror("Only Total Lagrange for SO_HEX20 implemented!");
   }
   // geometrically non-linear with Total Lagrangean approach
   else if (strncmp(buffer,"Totlag",6)==0)    kintype_ = soh20_totlag;
   // geometrically non-linear with Updated Lagrangean approach
   else if (strncmp(buffer,"Updlag",6)==0)
   {
       kintype_ = soh20_updlag;
       dserror("Only Total Lagrange for SO_HEX20 implemented!");
   }
   else dserror("Reading of SO_HEX20 element failed");
  }

  return true;
} // So_hex20::ReadElement()


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
