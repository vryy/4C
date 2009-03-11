/*!----------------------------------------------------------------------
\file so_material.cpp
\brief Everything concerning structural material selection for solid elements

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
#include "so_hex27.H"
#include "so_hex20.H"
#include "so_hex8.H"
#include "so_tet4.H"
#include "so_tet10.H"
#include "so_weg6.H"
#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_serialdensematrix.H"
#include "../drt_lib/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"
#include "../drt_mat/yeoh.H"
#include "../drt_mat/lung_penalty.H"
#include "../drt_mat/lung_ogden.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/contchainnetw.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/biocell.H"
#include "../drt_mat/material.H"
#include "../drt_mat/charmm.H"


using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp,
                    ParameterList&  params)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(glstrain,cmat,stress);
      *density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen = static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(glstrain,cmat,stress);

      *density = lungpen->Density();
      return;
      break;
    }
    case INPAR::MAT::m_visconeohooke: /*----------------- Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      /* Initialization moved to element input. So we can be sure, that material is initialized. */
      //if (!visco->Initialized())
      //  visco->Setup(NUMGPT_SOH8);
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      return;
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    case INPAR::MAT::m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,gp,Id(),time,cmat,stress);

      *density = anba->Density();
      return;
      break;
    }
    case INPAR::MAT::m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();
      return;
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
//      // Check if we use EAS
//      if (eastype_ != soh8_easnone)
//      {
//        // In this case, we have to calculate the "enhanced" deformation gradient
//        // from the enhanced GL strains with the help of two polar decompositions
//
//        // First step: determine enhanced material stretch tensor U_enh from C_enh=U_enh^T*U_enh
//        // -> get C_enh from enhanced GL strains
//        LINALG::Matrix<3,3> C_enh;
//        for (int i = 0; i < 3; ++i) C_enh(i,i) = 2.0 * (*glstrain)(i) + 1.0;
//        // off-diagonal terms are already twice in the Voigt-GLstrain-vector
//        C_enh(0,1) =  (*glstrain)(3);  C_enh(1,0) =  (*glstrain)(3);
//        C_enh(1,2) =  (*glstrain)(4);  C_enh(2,1) =  (*glstrain)(4);
//        C_enh(0,2) =  (*glstrain)(5);  C_enh(2,0) =  (*glstrain)(5);
//
//        // -> polar decomposition of (U^mod)^2
//        LINALG::Matrix<3,3> Q;
//        LINALG::Matrix<3,3> S;
//        LINALG::Matrix<3,3> VT;
//        LINALG::SVD<3,3>(C_enh,Q,S,VT); // Singular Value Decomposition
//        LINALG::Matrix<3,3> U_enh;
//        LINALG::Matrix<3,3> temp;
//        for (int i = 0; i < 3; ++i) S(i,i) = sqrt(S(i,i));
//        temp.MultiplyNN(Q,S);
//        U_enh.MultiplyNN(temp,VT);
//
//        // Second step: determine rotation tensor R from F (F=R*U)
//        // -> polar decomposition of displacement based F
//        LINALG::SVD<3,3>(*defgrd,Q,S,VT); // Singular Value Decomposition
//        LINALG::Matrix<3,3> R;
//        R.MultiplyNN(Q,VT);
//
//        // Third step: determine "enhanced" deformation gradient (F_enh=R*U_enh)
//        defgrd->MultiplyNN(R,U_enh);
//      }

      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_struct_multiscale: /*------------------- multiscale approach */
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

      // Check if we use EAS on this (macro-)scale
      if (eastype_ != soh8_easnone)
      {
        // In this case, we have to calculate the "enhanced" deformation gradient
        // from the enhanced GL strains with the help of two polar decompositions

        // First step: determine enhanced material stretch tensor U_enh from C_enh=U_enh^T*U_enh
        // -> get C_enh from enhanced GL strains
        LINALG::Matrix<3,3> C_enh;
        for (int i = 0; i < 3; ++i) C_enh(i,i) = 2.0 * (*glstrain)(i) + 1.0;
        // off-diagonal terms are already twice in the Voigt-GLstrain-vector
        C_enh(0,1) =  (*glstrain)(3);  C_enh(1,0) =  (*glstrain)(3);
        C_enh(1,2) =  (*glstrain)(4);  C_enh(2,1) =  (*glstrain)(4);
        C_enh(0,2) =  (*glstrain)(5);  C_enh(2,0) =  (*glstrain)(5);

        // -> polar decomposition of (U^mod)^2
        LINALG::Matrix<3,3> Q;
        LINALG::Matrix<3,3> S;
        LINALG::Matrix<3,3> VT;
        LINALG::SVD<3,3>(C_enh,Q,S,VT); // Singular Value Decomposition
        LINALG::Matrix<3,3> U_enh;
        LINALG::Matrix<3,3> temp;
        for (int i = 0; i < 3; ++i) S(i,i) = sqrt(S(i,i));
        temp.MultiplyNN(Q,S);
        U_enh.MultiplyNN(temp,VT);

        // Second step: determine rotation tensor R from F (F=R*U)
        // -> polar decomposition of displacement based F
        LINALG::SVD<3,3>(*defgrd,Q,S,VT); // Singular Value Decomposition
        LINALG::Matrix<3,3> R;
        R.MultiplyNN(Q,VT);

        // Third step: determine "enhanced" deformation gradient (F_enh=R*U_enh)
        defgrd->MultiplyNN(R,U_enh);
      }

      const double time = params.get<double>("total time",-1.0);
      const double dt = params.get<double>("delta time",-1.0);

      micro->Evaluate(defgrd, cmat, stress, density, gp, Id(), time, dt);
      return;
      break;
    }
    case INPAR::MAT::m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();
      return;
      break;
    }
    case INPAR::MAT::m_charmm: /*------------------------------------ CHARmm */
    {
      MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(mat.get());

      LINALG::SerialDenseMatrix XREFE(3,3);
      LINALG::SerialDenseMatrix XCURR(3,3);
      for (int i=0;i<3;i++)
      for (int j=0;j<3;j++) {
          //XREFE(i,j) = (*xrefe)(i,j);
          //XCURR(i,j) = (*xcurr)(i,j);
          XREFE(i,j) = 0.0; // Quick hack, that needs to be resoved
          XCURR(i,j) = 0.0;
      }
      const double time = params.get("total time",-1.0);
      charmm->Evaluate(glstrain,cmat,stress,Id(),gp,data_,time,XREFE,XCURR);
      *density = charmm->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of material");
    break;
  } // switch (mat->MaterialType())

  return;
} // of soh8_mat_sel


/*----------------------------------------------------------------------*
 | material laws for So_weg6                                   gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp,
                    ParameterList&  params)         // algorithmic parameters e.g. time
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore....
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(glstrain,cmat,stress);
      *density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen = static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(glstrain,cmat,stress);

      *density = lungpen->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of material");
    break;
  } // switch (mat->MaterialType())

}  // of sow6_mat_sel

/*----------------------------------------------------------------------*
 | material laws for So_hex27                                    tk 02/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex27::soh27_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp,
                    ParameterList&  params)
{

  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(glstrain,cmat,stress);
      *density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen = static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(glstrain,cmat,stress);

      *density = lungpen->Density();
      return;
      break;
    }
    case INPAR::MAT::m_visconeohooke: /*----------------- Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      return;
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    case INPAR::MAT::m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,gp,Id(),time,cmat,stress);

      *density = anba->Density();
      return;
      break;
    }
    case INPAR::MAT::m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();
      return;
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());

      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_struct_multiscale: /*------------------- multiscale approach */
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

      const double time = params.get<double>("total time",-1.0);
      const double dt = params.get<double>("delta time",-1.0);

      micro->Evaluate(defgrd, cmat, stress, density, gp, Id(), time, dt);
      return;
      break;
    }
    case INPAR::MAT::m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();
      return;
      break;
    }
    case INPAR::MAT::m_charmm: /*------------------------------------ CHARmm */
    {
      MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(mat.get());

      LINALG::SerialDenseMatrix XREFE(3,3);
      LINALG::SerialDenseMatrix XCURR(3,3);
      for (int i=0;i<3;i++)
      for (int j=0;j<3;j++) {
          //XREFE(i,j) = (*xrefe)(i,j);
          //XCURR(i,j) = (*xcurr)(i,j);
          XREFE(i,j) = 0.0; // Quick hack, that needs to be resoved
          XCURR(i,j) = 0.0;
      }
      const double time = params.get("total time",-1.0);
      charmm->Evaluate(glstrain,cmat,stress,Id(),gp,data_,time,XREFE,XCURR);
      *density = charmm->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of material");
    break;
  } // switch (mat->MaterialType())

  return;
} // of soh27_mat_sel

/*----------------------------------------------------------------------*
 | material laws for So_hex20                                   tk 02/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp,
                    ParameterList&  params)
{

  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(glstrain,cmat,stress);
      *density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen = static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(glstrain,cmat,stress);

      *density = lungpen->Density();
      return;
      break;
    }
    case INPAR::MAT::m_visconeohooke: /*----------------- Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      return;
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    case INPAR::MAT::m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,gp,Id(),time,cmat,stress);

      *density = anba->Density();
      return;
      break;
    }
    case INPAR::MAT::m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();
      return;
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());

      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_struct_multiscale: /*------------------- multiscale approach */
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

      const double time = params.get<double>("total time",-1.0);
      const double dt = params.get<double>("delta time",-1.0);

      micro->Evaluate(defgrd, cmat, stress, density, gp, Id(), time, dt);
      return;
      break;
    }
    case INPAR::MAT::m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();
      return;
      break;
    }
    case INPAR::MAT::m_charmm: /*------------------------------------ CHARmm */
    {
      MAT::CHARMM* charmm = static_cast <MAT::CHARMM*>(mat.get());

      LINALG::SerialDenseMatrix XREFE(3,3);
      LINALG::SerialDenseMatrix XCURR(3,3);
      for (int i=0;i<3;i++)
      for (int j=0;j<3;j++) {
          //XREFE(i,j) = (*xrefe)(i,j);
          //XCURR(i,j) = (*xcurr)(i,j);
          XREFE(i,j) = 0.0; // Quick hack, that needs to be resoved
          XCURR(i,j) = 0.0;
      }
      const double time = params.get("total time",-1.0);
      charmm->Evaluate(glstrain,cmat,stress,Id(),gp,data_,time,XREFE,XCURR);
      *density = charmm->Density();
      return;
      break;
    }
    default:
      dserror("Unknown type of material");
    break;
  } // switch (mat->MaterialType())

  return;
} // of soh20_mat_sel

/*----------------------------------------------------------------------*
 | material laws for SoDisp                                   maf 08/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::sodisp_mat_sel(
        LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
        LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
        double* density,
        LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
        ParameterList&            params)         // algorithmic parameters e.g. time
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here
  // the old interface does not exist anymore
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_ogden: /* lung tissue material with Ogden for volumetric part */
    {
      MAT::LungOgden* lungog = static_cast <MAT::LungOgden*>(mat.get());
      lungog->Evaluate(glstrain,cmat,stress);
      *density = lungog->Density();
      return;
      break;
    }
    case INPAR::MAT::m_lung_penalty: /* lung tissue material with penalty function for incompressibility constraint */
    {
      MAT::LungPenalty* lungpen= static_cast <MAT::LungPenalty*>(mat.get());

      lungpen->Evaluate(glstrain,cmat,stress);

      *density = lungpen->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    default:
    break;
  } // switch (mat->MaterialType())

  return;
} // of so_disp_mat_sel



/*----------------------------------------------------------------------*
 | material laws for So_tet4                                  gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here
  // the old interface does not exist anymore
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    default:
      dserror("Unknown material to tet4 element");
    break;
  } // switch (mat->MaterialType())

}

/*----------------------------------------------------------------------*
 | material laws for So_tet10                                   maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here, all others go through to
  // the old interface
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      return;
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      return;
      break;
    }
    default:
      dserror("Unknown material to tet10 element");
    break;
  } // switch (mat->MaterialType())

  return;
} // of So_tet10_mat_sel


#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
