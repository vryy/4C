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
#ifdef D_SOH8
#ifdef CCADISCRET

// This is just here to get the c++ mpi header, otherwise it would
// use the c version included inside standardtypes.h
#ifdef PARALLEL
#include "mpi.h"
#endif
#include "so_hex8.H"
#include "so_tet4.H"
#include "so_tet10.H"
#include "so_ctet10.H"
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

using namespace std; // cout etc.
using namespace LINALG; // our linear algebra

/*----------------------------------------------------------------------*
 | material laws for So_hex8                                   maf 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::soh8_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      const Epetra_SerialDenseMatrix* defgrd,
      const int gp,
      const int ele_ID,
      const double time,
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

      hypo->Evaluate(glstrain,defgrd,gp,ele_ID,time,cmat,stress);

      *density = hypo->Density();

      break;
    }
    case m_anisotropic_balzani:
    {
      MAT::AnisotropicBalzani* anba = static_cast <MAT::AnisotropicBalzani*>(mat.get());
      
      double avec[3]= {0.0};
      if (this->Type() != DRT::Element::element_sosh8){
        // fiber direction for z-cylinder, calculated via cross-product and beta=45°
        avec[0] = getthicknessvector()[1];
        avec[1] = -1.0 * getthicknessvector()[0];
        avec[2] = sqrt(avec[0]*avec[0] + avec[1]*avec[1]);
      }

      anba->Evaluate(glstrain,defgrd,gp,ele_ID,time,cmat,stress,avec);

      *density = anba->Density();

      break;
    }
    case m_struct_multiscale: /*------------------- multiscale approach */
    {
      // Here macro-micro transition (localization) will take place

      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());

      micro->Evaluate(defgrd, cmat, stress, density, gp, ele_ID, time, action);

      // test case
      //micromat->CalcStressStiffDens(stress, cmat, density, glstrain);

      // perform microscale simulation
      //micromat->PerformMicroSimulation();


      break;
    }
    case m_neohooke: /*----------------- NeoHookean Material */
    {
      MAT::NeoHooke* neo = static_cast <MAT::NeoHooke*>(mat.get());
      neo->Evaluate(glstrain,cmat,stress);
      *density = neo->Density();

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
 | material laws for So_weg6                                   maf 08/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_weg6::sow6_mat_sel(
      Epetra_SerialDenseVector* stress,
      Epetra_SerialDenseMatrix* cmat,
      double* density,
      const Epetra_SerialDenseVector* glstrain,
      const double time)
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
      const double time)
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
    default:
      dserror("Illegal type %d of material for element solid3 disp", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of sow6_mat_sel

#ifdef D_SOTET
/*----------------------------------------------------------------------* !!!!
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
    default:
      dserror("Illegal type %d of material for element solid3 hex8", mat->MaterialType());
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
    default:
      dserror("Illegal type %d of material for element solid3 hex8", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of so_tet10_mat_sel

/*----------------------------------------------------------------------* !!!!
 | material laws for So_ctet10                                  vlf 04/07|
 | added as a fast solution by cloning soh8_mat_sel (which is inside a  |
 | different class, and therefore cannot be used in So_tet10)           |
 | one should think about other solutions for this like making this a   |
 | member of a upper class the will be inherited by all so3 classes     |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_ctet10::so_ctet10_mat_sel(
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
    default:
      dserror("Illegal type %d of material for element solid3 hex8", mat->MaterialType());
      break;
  }

  /*--------------------------------------------------------------------*/
  return;
}  // of so_ctet10_mat_sel

#endif //SO_TET

#endif  // #ifdef CCADISCRET
#endif  // #ifdef D_SOH8
