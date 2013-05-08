/*!----------------------------------------------------------------------
\file so_hex8_thermo_evaluate.cpp
\brief

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*----------------------------------------------------------------------*/

#include "so_hex8.H"

#include "../drt_lib/drt_utils.H"

// include the headers of temperature-dependent materials with history only
#include "../drt_mat/thermostvenantkirchhoff.H"
#include "../drt_mat/robinson.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                           dano 05/10 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::linstiffmass(
  std::vector<int>& lm,  // location matrix
  std::vector<double>& disp,  // current displacements
  std::vector<double>& residual,  // current residual displacements or displacement increment
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* stiffmatrix,  // element stiffness matrix
  LINALG::Matrix<NUMDOF_SOH8,NUMDOF_SOH8>* massmatrix,  // element mass matrix
  LINALG::Matrix<NUMDOF_SOH8,1>* force,  // element internal force vector
  LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* elestress,  // stresses at GP
  LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* elestrain,  // strains at GP
  LINALG::Matrix<NUMGPT_SOH8,MAT::NUM_STRESS_3D>* eleplstrain, // plastic strains at GP
  Teuchos::ParameterList& params,  // algorithmic parameters e.g. time
  const INPAR::STR::StressType iostress,  // stress output option
  const INPAR::STR::StrainType iostrain,  // strain output option
  const INPAR::STR::StrainType ioplstrain  // plastic strain output option
  )
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_8 with 8 GAUSS POINTS*
** ============================================================================*/
  const static std::vector<LINALG::Matrix<NUMNOD_SOH8,1> > shapefcts = soh8_shapefcts();
  const static std::vector<LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> > derivs = soh8_derivs();
  const static std::vector<double> gpweights = soh8_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xcurr;  // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH8,NUMDIM_SOH8> xdisp;

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH8; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH8+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH8+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH8+2];
  }

  LINALG::Matrix<NUMDOF_SOH8,1> nodaldisp;
  // in case of Robinson's material, the (residual) displacements are required
  // residual displacements correspond to current displacement increment
  LINALG::Matrix<NUMDOF_SOH8,1> res_d;
  for (int i=0; i<NUMDOF_SOH8; ++i)
  {
    nodaldisp(i,0) = disp[i];
    res_d(i) = residual[i];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH8,NUMNOD_SOH8> N_XYZ;
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> defgrd(true);
  for (int gp=0; gp<NUMGPT_SOH8; ++gp)
  {
    /* get the inverse of the Jacobian matrix which looks like:
    **            [ x_,r  y_,r  z_,r ]^-1
    **     J^-1 = [ x_,s  y_,s  z_,s ]
    **            [ x_,t  y_,t  z_,t ]
    */
    // compute derivatives N_XYZ at gp w.r.t. material coordinates
    // by N_XYZ = J^-1 * N_rst
    N_XYZ.Multiply(invJ_[gp],derivs[gp]);
    double detJ = detJ_[gp];

    // WATCH OUT: here is the difference to the non-linear method nlnstiffmass()
    // in geometrically linear analysis the deformation gradient is equal to identity
    // no difference between reference and current state
    for (int i=0; i<3; ++i) defgrd(i,i) = 1.0;

    // linear B-operator B = N_XYZ
    // disperse global derivatives to bop-lines
    // bop is arranged as usual (refer to script FE or elsewhere):
    //
    // [ N1,X  0  0  | N2,X  0  0  | ... | Ni,X  0  0  ]
    // [ 0  N1,Y  0  | 0  N2,Y  0  | ... | 0  Ni,Y  0  ]
    // [ 0  0  N1,Z  | 0  0  N2,Z  | ... | 0  0  Ni,Z  ]
    // [ N1,Y N1,X 0 | N2,Y N2,X 0 | ... | Ni,Y Ni,X 0 ]
    // [ 0 N1,Z N1,Y | 0 N2,Z N2,Y | ... | 0 Ni,Z Ni,Y ]
    // [ N1,Z 0 N1,X | N2,Z 0 N2,X | ... | Ni,Z 0 Ni,X ]
    LINALG::Matrix<MAT::NUM_STRESS_3D,NUMDOF_SOH8> boplin;
    for (int i=0; i<NUMNOD_SOH8; ++i)
    {
      boplin(0,NODDOF_SOH8*i+0) = N_XYZ(0,i);
      boplin(0,NODDOF_SOH8*i+1) = 0.0;
      boplin(0,NODDOF_SOH8*i+2) = 0.0;
      boplin(1,NODDOF_SOH8*i+0) = 0.0;
      boplin(1,NODDOF_SOH8*i+1) = N_XYZ(1,i);
      boplin(1,NODDOF_SOH8*i+2) = 0.0;
      boplin(2,NODDOF_SOH8*i+0) = 0.0;
      boplin(2,NODDOF_SOH8*i+1) = 0.0;
      boplin(2,NODDOF_SOH8*i+2) = N_XYZ(2,i);
      /* ~~~ */
      boplin(3,NODDOF_SOH8*i+0) = N_XYZ(1,i);
      boplin(3,NODDOF_SOH8*i+1) = N_XYZ(0,i);
      boplin(3,NODDOF_SOH8*i+2) = 0.0;
      boplin(4,NODDOF_SOH8*i+0) = 0.0;
      boplin(4,NODDOF_SOH8*i+1) = N_XYZ(2,i);
      boplin(4,NODDOF_SOH8*i+2) = N_XYZ(1,i);
      boplin(5,NODDOF_SOH8*i+0) = N_XYZ(2,i);
      boplin(5,NODDOF_SOH8*i+1) = 0.0;
      boplin(5,NODDOF_SOH8*i+2) = N_XYZ(0,i);
    }

    // approximate linearised strain tensor using common naming of strain vector
    // glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(MAT::NUM_STRESS_3D);
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> glstrain(glstrain_epetra.A(),true);
    // E = epsilon_GL == epsilon_1
    // build the linearised strain epsilon = B . d
    glstrain.Multiply(boplin,nodaldisp);

    // return gp strains (only in case of stress/strain output)
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
      if (elestrain == NULL) dserror("strain data not available");
      for (int i = 0; i < 3; ++i)
        (*elestrain)(gp,i) = glstrain(i);
      for (int i = 3; i < 6; ++i)
        (*elestrain)(gp,i) = 0.5 * glstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
      if (elestrain == NULL) dserror("strain data not available");

      // e = F^{T-1} . E . F^{-1}
      LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
      GLtoEA(&glstrain,&defgrd,&euler_almansi);

      (*elestrain)(gp,0) = euler_almansi(0,0);
      (*elestrain)(gp,1) = euler_almansi(1,1);
      (*elestrain)(gp,2) = euler_almansi(2,2);
      (*elestrain)(gp,3) = euler_almansi(0,1);
      (*elestrain)(gp,4) = euler_almansi(1,2);
      (*elestrain)(gp,5) = euler_almansi(0,2);
    }
    break;
    case INPAR::STR::strain_none:
      break;
    default:
      dserror("requested strain type not available");
      break;
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    double scalartemp = 0.0;
    LINALG::Matrix<MAT::NUM_STRESS_3D,MAT::NUM_STRESS_3D> cmat(true);
    LINALG::Matrix<MAT::NUM_STRESS_3D,1> stress(true);

    // default: material call in structural function is purely deformation dependent
    bool young_temp = (params.get<int>("young_temp") == 1);
    if ( (Material()->MaterialType() == INPAR::MAT::m_thermostvenant) && (young_temp==true) )
    {
      Teuchos::RCP<std::vector<double> > temperature_vector
        = params.get<Teuchos::RCP<std::vector<double> > >("nodal_tempnp",Teuchos::null);

      double scalartemp = 0.0;
      // in StructureBaseAlgorithm() temperature not yet available, i.e. ==null
      if (temperature_vector==Teuchos::null)
      {
        MAT::ThermoStVenantKirchhoff* thrstvenant
          = static_cast <MAT::ThermoStVenantKirchhoff*>(Material().get());
        // initialise the temperature field
        scalartemp = thrstvenant->InitTemp();
      }
      // temperature vector is available
      else  // (temperature_vector!=Teuchos::null)
      {
        // get the temperature vector by extraction from parameter list
        LINALG::Matrix<NUMNOD_SOH8,1> etemp(true);
        for (int i=0; i<NUMNOD_SOH8; ++i)
        {
          etemp(i,0) = (*temperature_vector)[i];
        }
        // copy structural shape functions needed for the thermo field
        // identical shapefunctions for the displacements and the temperatures
        scalartemp  = (shapefcts[gp]).Dot(etemp);
      }

      // now set the current temperature vector in the parameter list
      params.set<double>("scalartemp",scalartemp);
    }
    // if Robinson's material --> pass the current temperature to the material
    else if (Material()->MaterialType() == INPAR::MAT::m_vp_robinson)
    {
      // scalar-valued temperature: T = shapefunctions . element temperatures
      // T = N_T^(e) . T^(e)
      // get the temperature vector by extraction from parameter list
      LINALG::Matrix<NUMNOD_SOH8,1> etemp(true);
      LINALG::Matrix<1,1> Ntemp(false);
      LINALG::Matrix<MAT::NUM_STRESS_3D,1> ctemp(true);

      Teuchos::RCP<std::vector<double> > temperature_vector
        = params.get<Teuchos::RCP<std::vector<double> > >("nodal_tempnp",Teuchos::null);
      // in StructureBaseAlgorithm() temperature not yet available, i.e. ==null
      if (temperature_vector==Teuchos::null)
      {
        MAT::Robinson* robinson
          = static_cast <MAT::Robinson*>(Material().get());
        // initialise the temperature field
        scalartemp = robinson->InitTemp();
      }
      // temperature vector is available
      else  // (temperature_vector!=Teuchos::null)
      {
        for (int i=0; i<NUMNOD_SOH8; ++i)
        {
          etemp(i,0) = (*temperature_vector)[i];
        }
        // copy structural shape functions needed for the thermo field
        // identical shapefunctions for the displacements and the temperatures
        scalartemp  = (shapefcts[gp]).Dot(etemp);
      }
      // now set the current temperature vector in the parameter list
      params.set<double>("scalartemp",scalartemp);

      // robinson material is solved using incremental strains
      // calculate incremental strains: Delta strain = B . Delta disp
      LINALG::Matrix<MAT::NUM_STRESS_3D,1> straininc(true);
      straininc.Multiply(boplin,res_d);

      params.set<LINALG::Matrix<MAT::NUM_STRESS_3D,1> >("straininc", straininc);

    } // end Robinson's material
    // default: material call in structural function is purely deformation dependent
    params.set<int>("gp",gp);
    params.set<int>("eleID",Id());
    Teuchos::RCP<MAT::So3Material> so3mat = Teuchos::rcp_dynamic_cast<MAT::So3Material>(Material());
    so3mat->Evaluate(&defgrd,&glstrain,params,&stress,&cmat);
    density = Material()->Density();

    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp plastic strains (only in case of plastic strain output)
    switch (ioplstrain)
    {
    case INPAR::STR::strain_gl:
    {
     if (eleplstrain == NULL) dserror("plastic strain data not available");
     LINALG::Matrix<MAT::NUM_STRESS_3D,1> plglstrain = params.get<LINALG::Matrix<MAT::NUM_STRESS_3D,1> >("plglstrain");
     for (int i = 0; i < 3; ++i)
       (*eleplstrain)(gp,i) = plglstrain(i);
     for (int i = 3; i < 6; ++i)
       (*eleplstrain)(gp,i) = 0.5 * plglstrain(i);
    }
    break;
    case INPAR::STR::strain_ea:
    {
     if (eleplstrain == NULL) dserror("plastic strain data not available");
     LINALG::Matrix<MAT::NUM_STRESS_3D,1> plglstrain = params.get<LINALG::Matrix<MAT::NUM_STRESS_3D,1> >("plglstrain");

     // e = F^{T-1} . E . F^{-1}
     LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> euler_almansi;
     GLtoEA(&plglstrain,&defgrd,&euler_almansi);

     (*eleplstrain)(gp,0) = euler_almansi(0,0);
     (*eleplstrain)(gp,1) = euler_almansi(1,1);
     (*eleplstrain)(gp,2) = euler_almansi(2,2);
     (*eleplstrain)(gp,3) = euler_almansi(0,1);
     (*eleplstrain)(gp,4) = euler_almansi(1,2);
     (*eleplstrain)(gp,5) = euler_almansi(0,2);
    }
    break;
    case INPAR::STR::strain_none:
     break;

    default:
     dserror("requested plastic strain type not available");
    }

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < MAT::NUM_STRESS_3D; ++i)
        (*elestress)(gp,i) = stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");

      // push forward of material stress to the spatial configuration
      LINALG::Matrix<3,3> cauchystress;
      PK2toCauchy(&stress,&defgrd,&cauchystress);

      (*elestress)(gp,0) = cauchystress(0,0);
      (*elestress)(gp,1) = cauchystress(1,1);
      (*elestress)(gp,2) = cauchystress(2,2);
      (*elestress)(gp,3) = cauchystress(0,1);
      (*elestress)(gp,4) = cauchystress(1,2);
      (*elestress)(gp,5) = cauchystress(0,2);
    }
    break;
    case INPAR::STR::stress_none:
      break;
    default:
      dserror("requested stress type not available");
    }

    double detJ_w = detJ*gpweights[gp];

    // update/integrate internal force vector
    if (force != NULL)
    {
      // f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, boplin, stress, 1.0);
    }

    // update/integrate `elastic' and `initial-displacement' stiffness matrix
    if (stiffmatrix != NULL)
    {
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH8> cb;
      cb.Multiply(cmat,boplin);
      stiffmatrix->MultiplyTN(detJ_w,boplin,cb,1.0);
    }

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH8; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH8; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH8*inod+0,NUMDIM_SOH8*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+1,NUMDIM_SOH8*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH8*inod+2,NUMDIM_SOH8*jnod+2) += massfactor;
        }
      }

    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex8::linstiffmass


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::GLtoEA(
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* glstrain,
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* defgrd,
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* euler_almansi
  )
{
  // e = F^{T-1} . E . F^{-1}

  // rewrite Green-Lagrange strain in tensor notation
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> gl;
  gl(0,0) = (*glstrain)(0);
  gl(0,1) = 0.5*(*glstrain)(3);
  gl(0,2) = 0.5*(*glstrain)(5);
  gl(1,0) = gl(0,1);
  gl(1,1) = (*glstrain)(1);
  gl(1,2) = 0.5*(*glstrain)(4);
  gl(2,0) = gl(0,2);
  gl(2,1) = gl(1,2);
  gl(2,2) = (*glstrain)(2);

  // inverse of deformation gradient
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> invdefgrd;
  invdefgrd.Invert((*defgrd));

  // (3x3) = (3x3) (3x3) (3x3)
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
  temp.Multiply(gl,invdefgrd);
  (*euler_almansi).MultiplyTN(invdefgrd,temp);

}  // GLtoEA()


/*----------------------------------------------------------------------*
 | push forward of material to spatial stresses              dano 11/12 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex8::PK2toCauchy(
  LINALG::Matrix<MAT::NUM_STRESS_3D,1>* stress,
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* defgrd,
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8>* cauchystress
  )
{
  // calculate the Jacobi-deterinant
  const double detF = (*defgrd).Determinant();

  // sigma = 1/J . F . S . F^T
  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> pkstress;
  pkstress(0,0) = (*stress)(0);
  pkstress(0,1) = (*stress)(3);
  pkstress(0,2) = (*stress)(5);
  pkstress(1,0) = pkstress(0,1);
  pkstress(1,1) = (*stress)(1);
  pkstress(1,2) = (*stress)(4);
  pkstress(2,0) = pkstress(0,2);
  pkstress(2,1) = pkstress(1,2);
  pkstress(2,2) = (*stress)(2);

  LINALG::Matrix<NUMDIM_SOH8,NUMDIM_SOH8> temp;
  temp.Multiply((1.0/detF),(*defgrd),pkstress);
  (*cauchystress).MultiplyNT(temp,(*defgrd));

}  // PK2toCauchy()


/*----------------------------------------------------------------------*/
