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
#include "../drt_mat/hyperpolyconvex.H"
#include "../drt_mat/neohooke.H"
#include "../drt_mat/anisotropic_balzani.H"
#include "../drt_mat/aaaneohooke.H"
#include "../drt_mat/mooneyrivlin.H"
#include "../drt_mat/hyperpolyconvex_ogden.H"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/contchainnetw.H"
#include "../drt_mat/artwallremod.H"
#include "../drt_mat/biocell.H"

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel(
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* stress,
                    LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                    double* density,
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                    LINALG::FixedSizeSerialDenseMatrix<3,3>* defgrd,
                    const int gp,
                    ParameterList&  params,
                    const string action)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::FixedSizeSerialDenseMatrix
  // interface go to the material law here, all others go through to
  // the old interface
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    default:
    break;
  } // switch (mat->MaterialType())
  
  
  // This is a wrapper for the Epetra style material interface
  Epetra_SerialDenseVector stress_e(View,stress->A(),stress->Rows());
  Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());
  const Epetra_SerialDenseVector glstrain_e(View,glstrain->A(),glstrain->Rows());
  Epetra_SerialDenseMatrix defgrd_e(View,defgrd->A(),defgrd->Rows(),defgrd->Rows(),defgrd->Columns());
  soh8_mat_sel(&stress_e,&cmat_e,density,&glstrain_e,&defgrd_e,gp,params,action);

  return;
}                    

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_mat_sel(
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* stress,
                    LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                    double* density,
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                    LINALG::FixedSizeSerialDenseMatrix<3,3>* defgrd,
                    const int gp)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  Epetra_SerialDenseVector stress_e(View,stress->A(),stress->Rows());
  Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());
  const Epetra_SerialDenseVector glstrain_e(View,glstrain->A(),glstrain->Rows());
  Epetra_SerialDenseMatrix defgrd_e(View,defgrd->A(),defgrd->Rows(),defgrd->Rows(),defgrd->Columns());

  so_tet10_mat_sel(&stress_e,&cmat_e,density,&glstrain_e,&defgrd_e,gp);

  return;
}     

/*----------------------------------------------------------------------*
 | material laws for SoDisp                                   maf 08/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::sodisp_mat_sel(
        LINALG::FixedSizeSerialDenseMatrix<6,1>* stress,
        LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
        double* density,
        LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
        ParameterList&            params)         // algorithmic parameters e.g. time
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
#endif

  Epetra_SerialDenseVector stress_e(View,stress->A(),stress->Rows());
  Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());
  const Epetra_SerialDenseVector glstrain_e(View,glstrain->A(),glstrain->Rows());

  sodisp_mat_sel(&stress_e,&cmat_e,density,&glstrain_e,params);

  return;
}     



/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      Epetra_SerialDenseMatrix* defgrd,
      const int gp,
      ParameterList&  params,
      const string action)
{
  RefCountPtr<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());

      stvk->Evaluate(glstrain,cmat,stress);

      *density = stvk->Density();

      break;
    }
    case m_hyper_polyconvex:
    {
      MAT::HyperPolyconvex* hypo = static_cast <MAT::HyperPolyconvex*>(mat.get());

      const double time = params.get("total time",-1.0);
      hypo->Evaluate(glstrain,defgrd,gp,Id(),time,cmat,stress);

      *density = hypo->Density();

      break;
    }
    case m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());

      const double time = params.get("total time",-1.0);
      anba->Evaluate(glstrain,defgrd,gp,Id(),time,cmat,stress);

      *density = anba->Density();

      break;
    }
    case m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();

      break;
    }
    case m_visconeohooke: /*----------------- Viscous NeoHookean Material */
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      /* Initialization moved to element input. So we can be sure, that material is initialized. */ 
      //if (!visco->Initialized())
      //  visco->Setup(NUMGPT_SOH8);
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();

      break;
    }
    case m_viscoanisotropic: /*------- Viscous Anisotropic Fiber Material */
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Evaluate(glstrain,gp,params,cmat,stress);
      *density = visco->Density();

      break;
    }
    case m_contchainnetw: /*------------ Continuum Chain Network Material */
    {
      MAT::ContChainNetw* chain = static_cast <MAT::ContChainNetw*>(mat.get());
      if (!chain->Initialized())
        chain->Initialize(NUMGPT_SOH8, this->Id());
      chain->Evaluate(glstrain,gp,params,cmat,stress,this->Id());
      *density = chain->Density();

      break;
    }
    case m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
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
//        LINALG::SerialDenseMatrix C_enh(NUMDIM_SOH8,NUMDIM_SOH8);
//        for (int i = 0; i < NUMDIM_SOH8; ++i) C_enh(i,i) = 2.0 * (*glstrain)(i) + 1.0;
//        // off-diagonal terms are already twice in the Voigt-GLstrain-vector
//        C_enh(0,1) =  (*glstrain)(3);  C_enh(1,0) =  (*glstrain)(3);
//        C_enh(1,2) =  (*glstrain)(4);  C_enh(2,1) =  (*glstrain)(4);
//        C_enh(0,2) =  (*glstrain)(5);  C_enh(2,0) =  (*glstrain)(5);
//
//        // -> polar decomposition of (U^mod)^2
//        LINALG::SerialDenseMatrix Q(NUMDIM_SOH8,NUMDIM_SOH8);
//        LINALG::SerialDenseMatrix S(NUMDIM_SOH8,NUMDIM_SOH8);
//        LINALG::SerialDenseMatrix VT(NUMDIM_SOH8,NUMDIM_SOH8);
//        SVD(C_enh,Q,S,VT); // Singular Value Decomposition
//        LINALG::SerialDenseMatrix U_enh(NUMDIM_SOH8,NUMDIM_SOH8);
//        LINALG::SerialDenseMatrix temp(NUMDIM_SOH8,NUMDIM_SOH8);
//        for (int i = 0; i < NUMDIM_SOH8; ++i) S(i,i) = sqrt(S(i,i));
//        temp.Multiply('N','N',1.0,Q,S,0.0);
//        U_enh.Multiply('N','N',1.0,temp,VT,0.0);
//
//        // Second step: determine rotation tensor R from F (F=R*U)
//        // -> polar decomposition of displacement based F
//        SVD(*defgrd,Q,S,VT); // Singular Value Decomposition
//        LINALG::SerialDenseMatrix R(NUMDIM_SOH8,NUMDIM_SOH8);
//        R.Multiply('N','N',1.0,Q,VT,0.0);
//
//        // Third step: determine "enhanced" deformation gradient (F_enh=R*U_enh)
//        defgrd->Multiply('N','N',1.0,R,U_enh,0.0);
//      }

      remo->Evaluate(glstrain,*defgrd,gp,params,cmat,stress,this->Id());
      *density = remo->Density();

      break;
    }
    case m_struct_multiscale: /*------------------- multiscale approach */
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

      // Check if we use EAS on this (macro-)scale
      if (eastype_ != soh8_easnone)
      {
        // In this case, we have to calculate the "enhanced" deformation gradient
        // from the enhanced GL strains with the help of two polar decompositions

        // First step: determine enhanced material stretch tensor U_enh from C_enh=U_enh^T*U_enh
        // -> get C_enh from enhanced GL strains
        LINALG::SerialDenseMatrix C_enh(NUMDIM_SOH8,NUMDIM_SOH8);
        for (int i = 0; i < NUMDIM_SOH8; ++i) C_enh(i,i) = 2.0 * (*glstrain)(i) + 1.0;
        // off-diagonal terms are already twice in the Voigt-GLstrain-vector
        C_enh(0,1) =  (*glstrain)(3);  C_enh(1,0) =  (*glstrain)(3);
        C_enh(1,2) =  (*glstrain)(4);  C_enh(2,1) =  (*glstrain)(4);
        C_enh(0,2) =  (*glstrain)(5);  C_enh(2,0) =  (*glstrain)(5);

        // -> polar decomposition of (U^mod)^2
        LINALG::SerialDenseMatrix Q(NUMDIM_SOH8,NUMDIM_SOH8);
        LINALG::SerialDenseMatrix S(NUMDIM_SOH8,NUMDIM_SOH8);
        LINALG::SerialDenseMatrix VT(NUMDIM_SOH8,NUMDIM_SOH8);
        SVD(C_enh,Q,S,VT); // Singular Value Decomposition
        LINALG::SerialDenseMatrix U_enh(NUMDIM_SOH8,NUMDIM_SOH8);
        LINALG::SerialDenseMatrix temp(NUMDIM_SOH8,NUMDIM_SOH8);
        for (int i = 0; i < NUMDIM_SOH8; ++i) S(i,i) = sqrt(S(i,i));
        temp.Multiply('N','N',1.0,Q,S,0.0);
        U_enh.Multiply('N','N',1.0,temp,VT,0.0);

        // Second step: determine rotation tensor R from F (F=R*U)
        // -> polar decomposition of displacement based F
        SVD(*defgrd,Q,S,VT); // Singular Value Decomposition
        LINALG::SerialDenseMatrix R(NUMDIM_SOH8,NUMDIM_SOH8);
        R.Multiply('N','N',1.0,Q,VT,0.0);

        // Third step: determine "enhanced" deformation gradient (F_enh=R*U_enh)
        defgrd->Multiply('N','N',1.0,R,U_enh,0.0);
      }

      const double time = params.get("total time",-1.0);
      micro->Evaluate(defgrd, cmat, stress, density, gp, Id(), time, action);

      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain,cmat,stress);
      *density = neo->Density();

      break;
    }
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(glstrain,cmat,stress);
      *density = aaa->Density();
      break;
    }
    case m_hyperpolyogden: /*-- slightly compressible hyperelastic polyconvex material for alveoli*/
    {
      MAT::HyperPolyOgden* hpo = static_cast <MAT::HyperPolyOgden*>(mat.get());
      hpo->Evaluate(glstrain,cmat,stress);
      *density = hpo->Density();
      break;
    }
    case m_biocell: /*----------------- Biological Cell Material */
    {
      MAT::BioCell* biocell = static_cast <MAT::BioCell*>(mat.get());
      biocell->Evaluate(glstrain,cmat,stress);
      *density = biocell->Density();

      break;
    }
    default:
      dserror("Illegal type %d of material for element solid3 hex8", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of soh8_mat_sel

/*----------------------------------------------------------------------*
 | material laws for So_weg6                                   gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_mat_sel(
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* stress,
                    LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                    double* density,
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                    LINALG::FixedSizeSerialDenseMatrix<3,3>* defgrd,
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

  // All materials that have a pure LINALG::FixedSizeSerialDenseMatrix
  // interface go to the material law here, all others go through to
  // the old interface
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    default:
    break;
  } // switch (mat->MaterialType())


  // This is a wrapper for the Epetra style material interface
  Epetra_SerialDenseVector stress_e(View,stress->A(),stress->Rows());
  Epetra_SerialDenseMatrix cmat_e(View,cmat->A(),cmat->Rows(),cmat->Rows(),cmat->Columns());
  const Epetra_SerialDenseVector glstrain_e(View,glstrain->A(),glstrain->Rows());
  Epetra_SerialDenseMatrix defgrd_e(View,defgrd->A(),defgrd->Rows(),defgrd->Rows(),defgrd->Columns());
  sow6_mat_sel(&stress_e,&cmat_e,density,&glstrain_e,&defgrd_e,gp,params);

}

/*----------------------------------------------------------------------*
 | material laws for So_weg6                                   maf 08/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      Epetra_SerialDenseMatrix* defgrd,
      const int gp,
      ParameterList&            params)         // algorithmic parameters e.g. time
{
  RefCountPtr<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());

      stvk->Evaluate(glstrain,cmat,stress);

      *density = stvk->Density();

      break;
    }
    case m_hyper_polyconvex:
    {
      MAT::HyperPolyconvex* hypo = static_cast <MAT::HyperPolyconvex*>(mat.get());

      hypo->Evaluate(glstrain,cmat,stress);

      *density = hypo->Density();

      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain, cmat, stress);
      *density = neo->Density();

    break;
    }
    case m_artwallremod: /*-Arterial Wall (Holzapfel) with remodeling (Hariton) */
    {
      MAT::ArtWallRemod* remo = static_cast <MAT::ArtWallRemod*>(mat.get());
      remo->Evaluate(glstrain,*defgrd,gp,params,cmat,stress);
      *density = remo->Density();

      break;
    }
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(glstrain,cmat,stress);
      *density = aaa->Density();
      break;
    }
    default:
      dserror("Illegal type %d of material for element solid3 weg 6", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of sow6_mat_sel

/*----------------------------------------------------------------------*
 | material laws for SoDisp                                   maf 08/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::SoDisp::sodisp_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      ParameterList&            params)         // algorithmic parameters e.g. time
{
  RefCountPtr<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());

      stvk->Evaluate(glstrain,cmat,stress);

      *density = stvk->Density();

      break;
    }
    case m_hyper_polyconvex:
    {
      MAT::HyperPolyconvex* hypo = static_cast <MAT::HyperPolyconvex*>(mat.get());

      hypo->Evaluate(glstrain,cmat,stress);

      *density = hypo->Density();

      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain, cmat, stress);
      *density = neo->Density();

    break;
    }
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(glstrain,cmat,stress);
      *density = aaa->Density();
      break;
    }
    default:
      dserror("Illegal type %d of material for element solid3 disp", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of sow6_mat_sel

/*----------------------------------------------------------------------*
 | material laws for So_tet4                                  gee 10/08|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_mat_sel(
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* stress,
                    LINALG::FixedSizeSerialDenseMatrix<6,6>* cmat,
                    double* density,
                    LINALG::FixedSizeSerialDenseMatrix<6,1>* glstrain,
                    LINALG::FixedSizeSerialDenseMatrix<3,3>* defgrd,
                    const int gp)
{
#ifdef DEBUG
  // I'm not sure whether all of these are always supplied, we'll see....
  if (!stress) dserror("No stress vector supplied");
  if (!cmat) dserror("No material tangent matrix supplied");
  if (!glstrain) dserror("No GL strains supplied");
  if (!defgrd) dserror("No defgrd supplied");
#endif

  // All materials that have a pure LINALG::FixedSizeSerialDenseMatrix
  // interface go to the material law here, all others go through to
  // the old interface
  RCP<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
      stvk->Evaluate(*glstrain,*cmat,*stress);
      *density = stvk->Density();
      return;
      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(*glstrain,*cmat,*stress);
      *density = neo->Density();
      return;
      break;
    }
    case m_aaaneohooke: /*-- special case of generalised NeoHookean material see Raghavan, Vorp */
    {
      MAT::AAAneohooke* aaa = static_cast <MAT::AAAneohooke*>(mat.get());
      aaa->Evaluate(*glstrain,*cmat,*stress);
      *density = aaa->Density();
      return;
      break;
    }
    default:
      dserror("Unknown material to tet4 element");
    break;
  } // switch (mat->MaterialType())

}

/*----------------------------------------------------------------------*
 | material laws for So_tet4                                  vlf 04/07|
 | added as a fast solution by cloning soh8_mat_sel (which is inside a  |
 | different class, and therefore cannot be used in So_tet10)           |
 | one should think about other solutions for this like making this a   |
 | member of a upper class the will be inherited by all so3 classes     |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet4::so_tet4_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      const Epetra_SerialDenseMatrix* defgrd,
      int gp)
{
  RefCountPtr<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());

      stvk->Evaluate(glstrain,cmat,stress);

      *density = stvk->Density();

      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain, cmat, stress);
      *density = neo->Density();

    break;
    }
    default:
      dserror("Illegal type %d of material for element solid3 tet4", mat->MaterialType());
    break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of so_tet4_mat_sel

/*----------------------------------------------------------------------* !!!!
 | material laws for So_tet10                                  vlf 04/07|
 | added as a fast solution by cloning soh8_mat_sel (which is inside a  |
 | different class, and therefore cannot be used in So_tet10)           |
 | one should think about other solutions for this like making this a   |
 | member of a upper class the will be inherited by all so3 classes     |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_tet10::so_tet10_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      const Epetra_SerialDenseMatrix* defgrd,
      int gp)
{
  RefCountPtr<MAT::Material> mat = Material();
  switch (mat->MaterialType())
  {
    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
    {
      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());

      stvk->Evaluate(glstrain,cmat,stress);

      *density = stvk->Density();

      break;
    }
    case m_mooneyrivlin: /*----------------- Mooney-Rivlin Material */
    {
      MAT::MooneyRivlin* moon = static_cast <MAT::MooneyRivlin*>(mat.get());
      moon->Evaluate(glstrain,cmat,stress);
      *density = moon->Density();

      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain, cmat, stress);
      *density = neo->Density();

    break;
    }
    default:
      dserror("Illegal type %d of material for element solid3 tet10", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of so_tet10_mat_sel



#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOLID3
