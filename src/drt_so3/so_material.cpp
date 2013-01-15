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
#include "so_nurbs27.H"
#include "so_hex27.H"
#include "so_hex20.H"
#include "so_hex8.H"
#include "so_tet4.H"
#include "so_tet10.H"
#include "so_weg6.H"
#include "so_disp.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"

#include "../drt_mat/micromaterial.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/thermoplasticlinelast.H"
#include "../drt_mat/plasticneohooke.H"
#include "../drt_mat/plastichyperelast.H"
#include "../drt_mat/plasticlinelast.H"
#include "../drt_mat/robinson.H"
#include "../drt_mat/damage.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/aaaneohooke_stopro.H"
#include "../drt_mat/aaagasser.H"
#include "../drt_mat/aaaraghavanvorp_damage.H"
#include "../drt_mat/aaa_mixedeffects.H"
#include "../drt_mat/logneohooke.H"
#include "../drt_mat/mooneyrivlin.H"
#include "../drt_mat/yeoh.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/contchainnetw.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/biocell.H"
#include "../drt_mat/material.H"
#include "../drt_mat/charmm.H"
#include "../drt_mat/itskov.H"
#include "../drt_mat/protein.H"
#include "../drt_mat/holzapfelcardiovascular.H"
#include "../drt_mat/humphreycardiovascular.H"
#include "../drt_mat/growth_ip.H"
#include "../drt_mat/constraintmixture.H"
#include "../drt_mat/structporo.H"
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* plglstrain,
                    LINALG::Matrix<3,3>* defgrd,
                    const int gp,
                    Teuchos::ParameterList&  params)
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

  if(mat->MaterialType() == INPAR::MAT::m_structporo)
  {
	  const MAT::StructPoro* actmat = static_cast<const MAT::StructPoro*>(mat.get());
	  mat = actmat->GetMaterial();
  }

  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      break;
    }
    // added this material only for hex8
    case INPAR::MAT::m_thermostvenant: /*------------------ st.venant-kirchhoff-material with temperature */
    {
      MAT::ThermoStVenantKirchhoff* thrstvk = static_cast <MAT::ThermoStVenantKirchhoff*>(mat.get());
      thrstvk->Evaluate(*glstrain,*cmat,*stress,params);
      *density = thrstvk->Density();
      break;
    }
    // added this material only for hex8
    case INPAR::MAT::m_thermopllinelast: /*-- linear thermo-elastio-plastic Material */
    {
      MAT::ThermoPlasticLinElast* thrpllinelast = static_cast <MAT::ThermoPlasticLinElast*>(mat.get());
      thrpllinelast->Evaluate(*glstrain,*plglstrain,gp,params,*cmat,*stress);
      *density = thrpllinelast->Density();
      break;
    }
    case INPAR::MAT::m_plneohooke: /*-- Plastic NeoHookean Material */
    {
      MAT::PlasticNeoHooke* plastic = static_cast <MAT::PlasticNeoHooke*>(mat.get());
      plastic->Evaluate(defgrd,gp,params,cmat,stress);
      *density = plastic->Density();
      break;
    }
    case INPAR::MAT::m_plhyperelast: /*-- plastic hyperelastic Material */
    {
      MAT::PlasticHyperElast* plastic = static_cast <MAT::PlasticHyperElast*>(mat.get());
      plastic->Evaluate(*glstrain,*cmat,*stress);
      *density = plastic->Density();
      break;
    }
    case INPAR::MAT::m_pllinelast: /*-- plastic linear elastic Material */
    {
      MAT::PlasticLinElast* pllinelast = static_cast <MAT::PlasticLinElast*>(mat.get());
      pllinelast->Evaluate(*glstrain,*plglstrain,gp,params,*cmat,*stress);
      *density = pllinelast->Density();
      break;
    }
    case INPAR::MAT::m_elpldamage: /*-- plastic linear elastic Material */
    {
      MAT::Damage* damage = static_cast <MAT::Damage*>(mat.get());
      damage->Evaluate(*glstrain,*plglstrain,gp,params,*cmat,*stress);
      *density = damage->Density();
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke_stopro: /*-- special case of generalised NeoHookean material see Raghavan, Vorp with stochastic mat parameters*/
    {
      MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(mat.get());
      aaa_stopro->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa_stopro->Density();
      break;
    }
    case INPAR::MAT::m_aaagasser: /*-- AAA thrombus material acc. to GASSER [2008] */
    {
      MAT::AAAgasser* gasser = static_cast<MAT::AAAgasser*>(mat.get());
      double normdist = params.get("iltthick meanvalue",-999.0);
      if (normdist==-999.0) dserror("Aneurysm mean ilt distance not found");
      gasser->Evaluate(*glstrain,*cmat,*stress, normdist);
      *density = gasser->Density();
      break;
    }
    case INPAR::MAT::m_aaa_mixedeffects: /*-- AAA mixed effect model */
    {
      MAT::AAA_mixedeffects* aaamixedeffects = static_cast<MAT::AAA_mixedeffects*>(mat.get());
      double localrad = params.get("localrad meanvalue",-999.0);
      if (localrad==-999.0) dserror("Aneurysm local radii not found");
      aaamixedeffects->Evaluate(*glstrain,*cmat,*stress, localrad);
      *density = aaamixedeffects->Density();
      break;
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage: /*-- special case of generalised NeoHookean material see Raghavan, Vorp, with damage */
    {
      MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
      /* Initialization moved to element input. So we can be sure, that material is initialized. */
      //if (!aaadamage->Initialized())
      //  aaadamage->Setup(NUMGPT_SOH8);
      aaadamage->Evaluate(glstrain,gp,params,cmat,stress);
      *density = aaadamage->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
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
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      //visco->UpdateFiberDirs(gp,defgrd);
      *density = visco->Density();
      break;
    }
    case INPAR::MAT::m_viscogenmax: /*------- Viscous Generalized Maxwell model compatible with hyperelastic toolbox  */
    {
      MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
      viscogenmax->Evaluate(*glstrain,*cmat,*stress,gp,params);
      *density = viscogenmax->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,gp,Id(),time,cmat,stress);

      *density = anba->Density();
      break;
    }
      case INPAR::MAT::m_itskov: //----------------- Itskov Material
    {
      MAT::Itskov* its = static_cast <MAT::Itskov*>(mat.get());
      const double time = params.get("total time",-1.0);
      its->Evaluate(*glstrain,gp,Id(),data_,time,*cmat,*stress);
      *density = its->Density();
      break;
    }
    case INPAR::MAT::m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
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

      micro->Evaluate(defgrd, cmat, stress, density, gp, Id());
      break;
    }
    case INPAR::MAT::m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();
      break;
    }
    case INPAR::MAT::m_charmm: /*------------------------------------ CHARMm */
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
      break;
    }
    case INPAR::MAT::m_protein: /*--------------------------- CHARMm Protein */
    {
      MAT::PROTEIN* protein = static_cast <MAT::PROTEIN*>(mat.get());

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
      protein->Evaluate(glstrain,cmat,stress,Id(),gp,data_,time,XREFE,XCURR);
      *density = protein->Density();
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
      break;
    }
    case INPAR::MAT::m_holzapfelcardiovascular: /*------- Anisotropic Fiber Material for arteries */
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Evaluate(glstrain,gp,cmat,stress);
      //holzcard->UpdateFiberDirs(gp,defgrd);
      *density = holzcard->Density();
      break;
    }
    case INPAR::MAT::m_humphreycardiovascular: /*------- Anisotropic Material for arteries cf Humphrey */
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Evaluate(glstrain,gp,cmat,stress);
      //humcard->UpdateFiberDirs(gp,defgrd);
      *density = humcard->Density();
      break;
    }
    case INPAR::MAT::m_growth: /*------- integration point based growth */
    {
      MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
      grow->Evaluate(glstrain,gp,cmat,stress,params);
      *density = grow->Density();
      break;
    }
    case INPAR::MAT::m_constraintmixture: /*------- growth and remodeling */
    {
      MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
      comix->Evaluate(glstrain,gp,cmat,stress,params);
      *density = comix->Density();
      break;
    }
    default:
      dserror("Unknown type of material");
    break;
  } // switch (mat->MaterialType())

  return;
} // soh8_mat_sel


/*----------------------------------------------------------------------*
 | material laws for So_hex8 with GEMM                        popp 02/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel_gemm(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat,
                    double* density,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_m,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_new,
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain_old,
                    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* rcg_new,
                    LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* rcg_old)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain_m) dserror("No GL strains supplied");
  if (!glstrain_new) dserror("No GL strains supplied");
  if (!glstrain_old) dserror("No GL strains supplied");
#endif

  // strain energy function
  double psi = 0.0;
  double psio = 0.0;

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  RCP<MAT::Material> mat = Material();

  if(mat->MaterialType() == INPAR::MAT::m_structporo)
  {
    const MAT::StructPoro* actmat = static_cast<const MAT::StructPoro*>(mat.get());
    mat = actmat->GetMaterial();
  }

  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------ St.Venant-Kirchhoff material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain_m,*cmat,*stress);
      stvk->StrainEnergy(*glstrain_new,psi);
      stvk->StrainEnergy(*glstrain_old,psio);
      *density = stvk->Density();
      break;
    }
    case INPAR::MAT::m_neohooke: /*--------------------- NeoHookean material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain_m,*cmat,*stress);
      neo->StrainEnergy(*glstrain_new,psi);
      neo->StrainEnergy(*glstrain_old,psio);
      *density = neo->Density();
      break;
    }
    default:
      dserror("Unknown type of material for GEMM");
    break;
  } // switch (mat->MaterialType())

  //**********************************************************************
  // ALGORITHMIC STRESSES AND CMAT FOR GEMM
  //**********************************************************************
  // tensor M = increment of Cauchy-Green tensor
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> M;
  M.Update(1.0,*rcg_new,-1.0,*rcg_old);
  double Mb = M.Dot(M);

  // second term in algorithmic stress only if Mb > 0
  // see: O. Gonzalez, Exact energy and momentum conserving algorithms for
  // general models in nonlinear elasticity, CMAME, 190(2000), pp. 1763-1783
  if (Mb < 1.0e-12) return;

  // derivative of strain energy function dpsi = 0.5*stressm
  // double contraction dpsi : M
  double dpsiM = 0.5*(*stress)(0)*M(0,0) + 0.5*(*stress)(1)*M(1,1) + 0.5*(*stress)(2)*M(2,2)
               +     (*stress)(3)*M(0,1) +     (*stress)(4)*M(1,2) +     (*stress)(5)*M(0,2);

  // extend stressm to algorithmic stress
  (*stress)(0) += 2 * ((psi - psio - dpsiM) / Mb) * M(0,0);
  (*stress)(1) += 2 * ((psi - psio - dpsiM) / Mb) * M(1,1);
  (*stress)(2) += 2 * ((psi - psio - dpsiM) / Mb) * M(2,2);
  (*stress)(3) += 2 * ((psi - psio - dpsiM) / Mb) * M(0,1);
  (*stress)(4) += 2 * ((psi - psio - dpsiM) / Mb) * M(1,2);
  (*stress)(5) += 2 * ((psi - psio - dpsiM) / Mb) * M(0,2);

  // TODO: extend cmat to algorithmic material tensor
  // -> not yet completely implemented!!!
  // -> using only cmat so far, which is ok but not optimal!!!

  //**********************************************************************
  //**********************************************************************
  //**********************************************************************

  return;
} // soh8_mat_sel_gemm


/*----------------------------------------------------------------------*
 | material laws with hyperelastic SEF for So_hex8            popp 02/12|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel_strainenergy(
                    LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
                    double* psi)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!glstrain) dserror("No GL strains supplied");
#endif

  // All materials that have a pure LINALG::Matrix
  // interface go to the material law here.
  // the old interface does not exist anymore...
  RCP<MAT::Material> mat = Material();

  if(mat->MaterialType() == INPAR::MAT::m_structporo)
  {
    const MAT::StructPoro* actmat = static_cast<const MAT::StructPoro*>(mat.get());
    mat = actmat->GetMaterial();
  }

  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------ St.Venant-Kirchhoff material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->StrainEnergy(*glstrain,*psi);
      break;
    }
    case INPAR::MAT::m_neohooke: /*--------------------- NeoHookean material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->StrainEnergy(*glstrain,*psi);
      break;
    }
    default:
      dserror("Hyperelastic strain energy not (yet) implemented for this material");
    break;
  } // switch (mat->MaterialType())

  return;
} // of soh8_mat_sel_strainenergy


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
                    Teuchos::ParameterList&  params)         // algorithmic parameters e.g. time
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
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke_stopro: /*-- special case of generalised NeoHookean material see Raghavan, Vorp with stochastic mat parameters*/
        {
          MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(mat.get());
          aaa_stopro->Evaluate(*glstrain,*cmat,*stress);
          *density = aaa_stopro->Density();
          break;
    }
    case INPAR::MAT::m_aaagasser: /*-- AAA thrombus material acc. to GASSER [2008] */
    {
      MAT::AAAgasser* gasser = static_cast<MAT::AAAgasser*>(mat.get());
      double normdist = params.get("iltthick meanvalue",-999.0);
      if (normdist==-999.0) dserror("Aneurysm mean ilt distance not found");
      gasser->Evaluate(*glstrain,*cmat,*stress, normdist);
      *density = gasser->Density();
      break;
    }
    case INPAR::MAT::m_aaa_mixedeffects: /*-- AAA mixed effect model */
    {
      MAT::AAA_mixedeffects* aaamixedeffects = static_cast<MAT::AAA_mixedeffects*>(mat.get());
      double localrad = params.get("localrad meanvalue",-999.0);
      if (localrad==-999.0) dserror("Aneurysm local radii not found");
      aaamixedeffects->Evaluate(*glstrain,*cmat,*stress, localrad);
      *density = aaamixedeffects->Density();
      break;
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage: /*-- special case of generalised NeoHookean material see Raghavan, Vorp, with damage */
    {
      MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
      aaadamage->Evaluate(glstrain,gp,params,cmat,stress);
      *density = aaadamage->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
      break;
    }
    case INPAR::MAT::m_holzapfelcardiovascular: /*------- Anisotropic Fiber Material for arteries */
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Evaluate(glstrain,gp,cmat,stress);
      *density = holzcard->Density();
      break;
    }
    case INPAR::MAT::m_humphreycardiovascular: /*------- Anisotropic Material for arteries cf Humphrey */
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Evaluate(glstrain,gp,cmat,stress);
      *density = humcard->Density();
      break;
    }
    case INPAR::MAT::m_growth: /*------- integration point based growth */
    {
      MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
      grow->Evaluate(glstrain,gp,cmat,stress,params);
      *density = grow->Density();
      break;
    }
    case INPAR::MAT::m_constraintmixture: /*------- growth and remodeling */
    {
      MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
      comix->Evaluate(glstrain,gp,cmat,stress,params);
      *density = comix->Density();
      break;
    }
    default:
      dserror("Unknown type of material");
    break;
  } // switch (mat->MaterialType())

}  // of sow6_mat_sel


/*----------------------------------------------------------------------*
 | material laws for So_nurbs27                              gammi 04/09|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::NURBS::So_nurbs27::sonurbs27_mat_sel(
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>*                  stress  ,
  LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D>* cmat    ,
  double*                                                density ,
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>*                  glstrain,
  LINALG::Matrix<3,3>*                                   defgrd  ,
  const int                                              gp      ,
  Teuchos::ParameterList&                                         params
  )
{

  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    default:
      dserror("Unknown type of material for nurbs27 implementation");
    break;
  } // switch (mat->MaterialType())

  return;
} // of sonurbs27_mat_sel


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
                    Teuchos::ParameterList&  params)
{

  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage: /*-- special case of generalised NeoHookean material see Raghavan, Vorp, with damage */
    {
      MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
      aaadamage->Evaluate(glstrain,gp,params,cmat,stress);
      *density = aaadamage->Density();
      break;
    }
    case INPAR::MAT::m_aaagasser: /*-- AAA thrombus material acc. to GASSER [2008] */
    {
      MAT::AAAgasser* gasser = static_cast<MAT::AAAgasser*>(mat.get());
      double normdist = params.get("iltthick meanvalue",-999.0);
      if (normdist==-999.0) dserror("Aneurysm mean ilt distance not found");
      gasser->Evaluate(*glstrain,*cmat,*stress, normdist);
      *density = gasser->Density();
      break;
    }
    case INPAR::MAT::m_aaa_mixedeffects: /*-- AAA mixed effect model */
    {
      MAT::AAA_mixedeffects* aaamixedeffects = static_cast<MAT::AAA_mixedeffects*>(mat.get());
      double localrad = params.get("localrad meanvalue",-999.0);
      if (localrad==-999.0) dserror("Aneurysm local radii not found");
      aaamixedeffects->Evaluate(*glstrain,*cmat,*stress, localrad);
      *density = aaamixedeffects->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
      break;
    }
    case INPAR::MAT::m_visconeohooke: /*----------------- Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      break;
    }
    case INPAR::MAT::m_viscogenmax: /*------- Viscous Generalized Maxwell model compatible with hyperelastic toolbox  */
    {
      MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
      viscogenmax->Evaluate(*glstrain,*cmat,*stress,gp,params);
      *density = viscogenmax->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,gp,Id(),time,cmat,stress);

      *density = anba->Density();
      break;
    }
    case INPAR::MAT::m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());

      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      break;
    }
    case INPAR::MAT::m_struct_multiscale: /*------------------- multiscale approach */
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->Evaluate(defgrd, cmat, stress, density, gp, Id());
      break;
    }
    case INPAR::MAT::m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();
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
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
      break;
    }
    case INPAR::MAT::m_holzapfelcardiovascular: /*------- Anisotropic Fiber Material for arteries */
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Evaluate(glstrain,gp,cmat,stress);
      *density = holzcard->Density();
      break;
    }
    case INPAR::MAT::m_humphreycardiovascular: /*------- Anisotropic Material for arteries cf Humphrey */
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Evaluate(glstrain,gp,cmat,stress);
      *density = humcard->Density();
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
                    Teuchos::ParameterList&  params)
{

  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaraghavanvorp_damage: /*-- special case of generalised NeoHookean material see Raghavan, Vorp, with damage */
    {
      MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
      aaadamage->Evaluate(glstrain,gp,params,cmat,stress);
      *density = aaadamage->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_aaagasser: /*-- AAA thrombus material acc. to GASSER [2008] */
    {
      MAT::AAAgasser* gasser = static_cast<MAT::AAAgasser*>(mat.get());
      double normdist = params.get("iltthick meanvalue",-999.0);
      if (normdist==-999.0) dserror("Aneurysm mean ilt distance not found");
      gasser->Evaluate(*glstrain,*cmat,*stress, normdist);
      *density = gasser->Density();
      break;
    }
    case INPAR::MAT::m_aaa_mixedeffects: /*-- AAA mixed effect model */
    {
      MAT::AAA_mixedeffects* aaamixedeffects = static_cast<MAT::AAA_mixedeffects*>(mat.get());
      double localrad = params.get("localrad meanvalue",-999.0);
      if (localrad==-999.0) dserror("Aneurysm local radii not found");
      aaamixedeffects->Evaluate(*glstrain,*cmat,*stress, localrad);
      *density = aaamixedeffects->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
      break;
    }
    case INPAR::MAT::m_visconeohooke: /*----------------- Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      break;
    }
    case INPAR::MAT::m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();
      break;
    }
    case INPAR::MAT::m_viscogenmax: /*------- Viscous Generalized Maxwell model compatible with hyperelastic toolbox  */
    {
      MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
      viscogenmax->Evaluate(*glstrain,*cmat,*stress,gp,params);
      *density = viscogenmax->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,gp,Id(),time,cmat,stress);

      *density = anba->Density();
      break;
    }
    case INPAR::MAT::m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();
      break;
    }
    case INPAR::MAT::m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());

      remo->Evaluate(glstrain,gp,params,cmat,stress,*defgrd);
      *density = remo->Density();
      break;
    }
    case INPAR::MAT::m_struct_multiscale: /*------------------- multiscale approach */
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->Evaluate(defgrd, cmat, stress, density, gp, Id());
      break;
    }
    case INPAR::MAT::m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();
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
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
      break;
    }
    case INPAR::MAT::m_holzapfelcardiovascular: /*------- Anisotropic Fiber Material for arteries */
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Evaluate(glstrain,gp,cmat,stress);
      *density = holzcard->Density();
      break;
    }
    case INPAR::MAT::m_humphreycardiovascular: /*------- Anisotropic Material for arteries cf Humphrey */
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Evaluate(glstrain,gp,cmat,stress);
      *density = humcard->Density();
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
        Teuchos::ParameterList&   params)         // algorithmic parameters e.g. time
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
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_aaa_mixedeffects: /*-- AAA mixed effect model */
    {
      MAT::AAA_mixedeffects* aaamixedeffects = static_cast<MAT::AAA_mixedeffects*>(mat.get());
      double localrad = params.get("localrad meanvalue",-999.0);
      if (localrad==-999.0) dserror("Aneurysm local radii not found");
      aaamixedeffects->Evaluate(*glstrain,*cmat,*stress, localrad);
      *density = aaamixedeffects->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
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
                    const int gp,
                    Teuchos::ParameterList& params
                    )
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


  if(mat->MaterialType() == INPAR::MAT::m_structporo)
  {
	  const MAT::StructPoro* actmat = static_cast<const MAT::StructPoro*>(mat.get());
	  mat = actmat->GetMaterial();
  }

  switch (mat->MaterialType())
  {
    case INPAR::MAT::m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      break;
    }
    // added this material only for hex8
    case INPAR::MAT::m_thermostvenant: /*------------------ st.venant-kirchhoff-material with temperature */
    {
      MAT::ThermoStVenantKirchhoff* thrstvk = static_cast <MAT::ThermoStVenantKirchhoff*>(mat.get());
      thrstvk->Evaluate(*glstrain,*cmat,*stress,params);
      *density = thrstvk->Density();
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke_stopro: /*-- special case of generalised NeoHookean material see Raghavan, Vorp with stochastic mat parameters*/
        {
          MAT::AAAneohooke_stopro* aaa_stopro = static_cast <MAT::AAAneohooke_stopro*>(mat.get());
          aaa_stopro->Evaluate(*glstrain,*cmat,*stress);
          *density = aaa_stopro->Density();
          break;
     }
    case INPAR::MAT::m_aaaraghavanvorp_damage: //-- special case of generalised NeoHookean material see Raghavan, Vorp, with damage
    {
      MAT::AAAraghavanvorp_damage* aaadamage = static_cast <MAT::AAAraghavanvorp_damage*>(mat.get());
      aaadamage->Evaluate(glstrain,gp,params,cmat,stress);
      *density = aaadamage->Density();
      break;
    }
    case INPAR::MAT::m_aaagasser: /*-- AAA thrombus material acc. to GASSER [2008] */
    {
      MAT::AAAgasser* gasser = static_cast<MAT::AAAgasser*>(mat.get());
      double normdist = params.get("iltthick meanvalue",-999.0);
      if (normdist==-999.0) dserror("Aneurysm mean ilt distance not found");
      gasser->Evaluate(*glstrain,*cmat,*stress, normdist);
      *density = gasser->Density();
      break;
    }
    case INPAR::MAT::m_aaa_mixedeffects: /*-- AAA mixed effect model */
    {
      MAT::AAA_mixedeffects* aaamixedeffects = static_cast<MAT::AAA_mixedeffects*>(mat.get());
      double localrad = params.get("localrad meanvalue",-999.0);
      if (localrad==-999.0) dserror("Aneurysm local radii not found");
      aaamixedeffects->Evaluate(*glstrain,*cmat,*stress, localrad);
      *density = aaamixedeffects->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
      break;
    }
    case INPAR::MAT::m_holzapfelcardiovascular: /*------- Anisotropic Fiber Material for arteries */
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Evaluate(glstrain,gp,cmat,stress);
      *density = holzcard->Density();
      break;
    }
    case INPAR::MAT::m_humphreycardiovascular: /*------- Anisotropic Material for arteries cf Humphrey */
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Evaluate(glstrain,gp,cmat,stress);
      *density = humcard->Density();
      break;
    }
    case INPAR::MAT::m_growth: /*------- integration point based growth */
    {
      MAT::Growth* grow = static_cast <MAT::Growth*>(mat.get());
      grow->Evaluate(glstrain,gp,cmat,stress,params);
      *density = grow->Density();
      break;
    }
    case INPAR::MAT::m_constraintmixture: /*------- growth and remodeling */
    {
      MAT::ConstraintMixture* comix = static_cast <MAT::ConstraintMixture*>(mat.get());
      comix->Evaluate(glstrain,gp,cmat,stress,params);
      *density = comix->Density();
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
                    const int gp
                    // ParameterList& params
                    )
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
      break;
    }
    case INPAR::MAT::m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      break;
    }
    case INPAR::MAT::m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      break;
    }
    case INPAR::MAT::m_logneohooke: /*-- logarithmic neo-Hookean material */
    {
      MAT::LogNeoHooke* logneo = static_cast <MAT::LogNeoHooke*>(mat.get());
      logneo->Evaluate(*glstrain,*cmat,*stress);
      *density = logneo->Density();
      break;
    }
    case INPAR::MAT::m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();
      break;
    }
    case INPAR::MAT::m_yeoh: /*----------------- Mooney-Rivlin Material */
    {
      MAT::Yeoh* yeoh = static_cast <MAT::Yeoh*>(mat.get());
      yeoh->Evaluate(glstrain,cmat,stress);
      *density = yeoh->Density();
      break;
    }
    case INPAR::MAT::m_elasthyper: /*----------- general hyperelastic matrial */
    {
      MAT::ElastHyper* hyper = static_cast <MAT::ElastHyper*>(mat.get());
      Teuchos::ParameterList params;
      hyper->Evaluate(*glstrain,*cmat,*stress,params);
      *density = hyper->Density();
      break;
    }
    case INPAR::MAT::m_holzapfelcardiovascular: /*------- Anisotropic Fiber Material for arteries */
    {
      MAT::HolzapfelCardio* holzcard = static_cast <MAT::HolzapfelCardio*>(mat.get());
      holzcard->Evaluate(glstrain,gp,cmat,stress);
      *density = holzcard->Density();
      break;
    }
    case INPAR::MAT::m_humphreycardiovascular: /*------- Anisotropic Material for arteries cf Humphrey */
    {
      MAT::HumphreyCardio* humcard = static_cast <MAT::HumphreyCardio*>(mat.get());
      humcard->Evaluate(glstrain,gp,cmat,stress);
      *density = humcard->Density();
      break;
    }
    default:
      dserror("Unknown material to tet10 element");
    break;
  } // switch (mat->MaterialType())

  return;
} // of So_tet10_mat_sel


