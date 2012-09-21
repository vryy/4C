/*!----------------------------------------------------------------------
\file so_hex20_evaluate.cpp
\brief

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#include "so_hex20.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensevector.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_mat/visconeohooke.H"
#include "../drt_mat/viscoanisotropic.H"
#include "../drt_mat/viscogenmax.H"
#include "../drt_mat/elasthyper.H"
#include "../drt_mortar/mortar_analytical.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/micromaterial.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex20::Evaluate(ParameterList& params,
                                    DRT::Discretization&      discretization,
                                    vector<int>&              lm,
                                    Epetra_SerialDenseMatrix& elemat1_epetra,
                                    Epetra_SerialDenseMatrix& elemat2_epetra,
                                    Epetra_SerialDenseVector& elevec1_epetra,
                                    Epetra_SerialDenseVector& elevec2_epetra,
                                    Epetra_SerialDenseVector& elevec3_epetra)
{
  LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20> elemat1(elemat1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20> elemat2(elemat2_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH20,1> elevec1(elevec1_epetra.A(),true);
  LINALG::Matrix<NUMDOF_SOH20,1> elevec2(elevec2_epetra.A(),true);
  // elevec3 is not used anyway

  // start with "none"
  DRT::ELEMENTS::So_hex20::ActionType act = So_hex20::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff")                        act = So_hex20::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff")                        act = So_hex20::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce")                   act = So_hex20::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass")                    act = So_hex20::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass")                    act = So_hex20::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass")                   act = So_hex20::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress")                          act = So_hex20::calc_struct_stress;
  else if (action=="calc_struct_eleload")                         act = So_hex20::calc_struct_eleload;
  else if (action=="calc_struct_fsiload")                         act = So_hex20::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep")                    act = So_hex20::calc_struct_update_istep;
  else if (action=="calc_struct_update_imrlike")                  act = So_hex20::calc_struct_update_imrlike;
  else if (action=="calc_struct_reset_istep")                     act = So_hex20::calc_struct_reset_istep;
  else if (action=="calc_struct_errornorms")                      act = So_hex20::calc_struct_errornorms;
  else if (action=="postprocess_stress")                          act = So_hex20::postprocess_stress;
  else if (action=="multi_readrestart")                           act = So_hex20::multi_readrestart;
  else if (action=="multi_calc_dens")                             act = So_hex20::multi_calc_dens;
  else dserror("Unknown type of action for So_hex20");
  // what should the element do
  switch(act)
  {
    // linear stiffness
    case calc_struct_linstiff:
    {
      // need current displacement and residual forces
      vector<double> mydisp(lm.size());
      for (unsigned i=0; i<mydisp.size(); ++i) mydisp[i] = 0.0;
      vector<double> myres(lm.size());
      for (unsigned i=0; i<myres.size(); ++i) myres[i] = 0.0;
      soh20_nlnstiffmass(lm,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
    }
    break;

    // nonlinear stiffness and internal force vector
    case calc_struct_nlnstiff:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20>* matptr = NULL;
      if (elemat1.IsInitialized()) matptr = &elemat1;

      // special case: geometrically linear
      if (kintype_ == DRT::ELEMENTS::So_hex20::soh20_linear)
      {
        soh20_linstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                          INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        soh20_nlnstiffmass(lm,mydisp,myres,matptr,NULL,&elevec1,NULL,NULL,params,
                          INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
    }
    break;

    // internal force vector only
    case calc_struct_internalforce:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix to apply linearised EAS-stuff onto
      LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20> myemat(true);

      // special case: geometrically linear
      if (kintype_ == DRT::ELEMENTS::So_hex20::soh20_linear)
      {
        soh20_linstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
                           INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        soh20_nlnstiffmass(lm,mydisp,myres,&myemat,NULL,&elevec1,NULL,NULL,params,
                          INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
    }
    break;

    // linear stiffness and consistent mass matrix
    case calc_struct_linstiffmass:
      dserror("Case 'calc_struct_linstiffmass' not yet implemented");
    break;

    // nonlinear stiffness, internal force vector, and consistent mass matrix
    case calc_struct_nlnstiffmass:
    case calc_struct_nlnstifflmass:
    {
      // need current displacement and residual forces
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==null || res==null) dserror("Cannot get state vectors 'displacement' and/or residual");
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);

      // special case: geometrically linear
      if (kintype_ == DRT::ELEMENTS::So_hex20::soh20_linear)
      {
        soh20_linstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
                           INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        soh20_nlnstiffmass(lm,mydisp,myres,&elemat1,&elemat2,&elevec1,NULL,NULL,params,
                           INPAR::STR::stress_none,INPAR::STR::strain_none);
      }

      if (act==calc_struct_nlnstifflmass) soh20_lumpmass(&elemat2);
    }
    break;

    // evaluate stresses and strains at gauss points
    case calc_struct_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        RefCountPtr<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        RCP<vector<char> > stressdata = params.get<RCP<vector<char> > >("stress", null);
        RCP<vector<char> > straindata = params.get<RCP<vector<char> > >("strain", null);
        if (disp==null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==null) dserror("Cannot get 'stress' data");
        if (straindata==null) dserror("Cannot get 'strain' data");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20> stress;
        LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20> strain;
        INPAR::STR::StressType iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
        INPAR::STR::StrainType iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);

        // special case: geometrically linear
        if (kintype_ == DRT::ELEMENTS::So_hex20::soh20_linear)
        {
          soh20_linstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);
        }
        // standard is: geometrically non-linear with Total Lagrangean approach
        else
        {
          soh20_nlnstiffmass(lm,mydisp,myres,NULL,NULL,NULL,&stress,&strain,params,iostress,iostrain);
        }

        {
          DRT::PackBuffer data;
          AddtoPack(data, stress);
          data.StartPacking();
          AddtoPack(data, stress);
          std::copy(data().begin(),data().end(),std::back_inserter(*stressdata));
        }

        {
          DRT::PackBuffer data;
          AddtoPack(data, strain);
          data.StartPacking();
          AddtoPack(data, strain);
          std::copy(data().begin(),data().end(),std::back_inserter(*straindata));
        }
      }
    }
    break;

    // postprocess stresses/strains at gauss points

    // note that in the following, quantities are always referred to as
    // "stresses" etc. although they might also apply to strains
    // (depending on what this routine is called for from the post filter)
    case postprocess_stress:
    {
      const RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
        params.get<RCP<map<int,RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",null);
      if (gpstressmap==null)
        dserror("no gp stress/strain map available for postprocessing");
      string stresstype = params.get<string>("stresstype","ndxyz");
      int gid = Id();
      LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20> gpstress(((*gpstressmap)[gid])->A(),true);

      Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",null);
      if (poststress==Teuchos::null)
        dserror("No element stress/strain vector available");

      if (stresstype=="ndxyz")
      {
        // extrapolate stresses/strains at Gauss points to nodes
        soh20_expol(gpstress,*poststress);
      }
      else if (stresstype=="cxyz")
      {
        const Epetra_BlockMap& elemap = poststress->Map();
        int lid = elemap.LID(Id());
        if (lid!=-1)
        {
          for (int i = 0; i < NUMSTR_SOH20; ++i)
          {
            double& s = (*((*poststress)(i)))[lid]; // resolve pointer for faster access
            s = 0.;
            for (int j = 0; j < NUMGPT_SOH20; ++j)
            {
              s += gpstress(j,i);
            }
            s *= 1.0/NUMGPT_SOH20;
          }
        }
      }
      else
      {
        dserror("unknown type of stress/strain output on element level");
      }
    }
    break;

    case calc_struct_eleload:
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
    break;

    case calc_struct_fsiload:
      dserror("Case not yet implemented");
    break;

    case calc_struct_update_istep:
    {
      // Update of history for visco material
      RCP<MAT::Material> mat = Material();
      if (mat->MaterialType() == INPAR::MAT::m_visconeohooke)
      {
        MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
        visco->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_viscoanisotropic)
      {
        MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
        visco->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_viscogenmax)
      {
        MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
        viscogenmax->Update();
      }
      else if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
      {
        MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
        micro->Update();
      }
    }
    break;

  case calc_struct_update_imrlike:
  {
    // Update of history for visco material
    RefCountPtr<MAT::Material> mat = Material();
    if (mat->MaterialType() == INPAR::MAT::m_visconeohooke)
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Update();
    }
    else if (mat->MaterialType() == INPAR::MAT::m_viscoanisotropic)
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Update();
    }
    else if (mat->MaterialType() == INPAR::MAT::m_viscogenmax)
    {
      MAT::ViscoGenMax* viscogenmax = static_cast <MAT::ViscoGenMax*>(mat.get());
      viscogenmax->Update();
    }
    else if (mat->MaterialType() == INPAR::MAT::m_struct_multiscale)
    {
      MAT::MicroMaterial* micro = static_cast <MAT::MicroMaterial*>(mat.get());
      micro->Update();
    }
  }
  break;

  case calc_struct_reset_istep:
  {
    // Reset of history for visco material
    RefCountPtr<MAT::Material> mat = Material();
    if (mat->MaterialType() == INPAR::MAT::m_visconeohooke)
    {
      MAT::ViscoNeoHooke* visco = static_cast <MAT::ViscoNeoHooke*>(mat.get());
      visco->Reset();
    }
    else if (mat->MaterialType() == INPAR::MAT::m_viscoanisotropic)
    {
      MAT::ViscoAnisotropic* visco = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      visco->Reset();
    }
    else if (mat->MaterialType() == INPAR::MAT::m_viscogenmax)
    {
      MAT::ViscoAnisotropic* viscogenmax = static_cast <MAT::ViscoAnisotropic*>(mat.get());
      viscogenmax->Reset();
    }
  }
  break;

  //==================================================================================
  case calc_struct_errornorms:
  {
    // IMPORTANT NOTES (popp 10/2010):
    // - error norms are based on a small deformation assumption (linear elasticity)
    // - extension to finite deformations would be possible without difficulties,
    //   however analytical solutions are extremely rare in the nonlinear realm
    // - only implemented for SVK material (relevant for energy norm only, L2 and
    //   H1 norms are of course valid for arbitrary materials)
    // - analytical solutions are currently stored in a repository in the MORTAR
    //   namespace, however they could (should?) be moved to a more general location

    // check length of elevec1
    if (elevec1_epetra.Length() < 3) dserror("The given result vector is too short.");

    // check material law
    RCP<MAT::Material> mat = Material();

    //******************************************************************
        // only for St.Venant Kirchhoff material
        //******************************************************************
            if (mat->MaterialType() == INPAR::MAT::m_stvenant)
          {
            // declaration of variables
            double l2norm = 0.0;
            double h1norm = 0.0;
            double energynorm = 0.0;

            // shape functions, derivatives and integration weights
            const static vector<LINALG::Matrix<NUMNOD_SOH20,1> > vals = soh20_shapefcts();
            const static vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > derivs = soh20_derivs();
            const static std::vector<double> weights = soh20_weights();

            // get displacements and extract values of this element
            RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
            if (disp==null) dserror("Cannot get state displacement vector");
            vector<double> mydisp(lm.size());
            DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

            // nodal displacement vector
            LINALG::Matrix<NUMDOF_SOH20,1> nodaldisp;
            for (int i=0; i<NUMDOF_SOH20; ++i) nodaldisp(i,0) = mydisp[i];

            // reference geometry (nodal positions)
            LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xrefe;
            DRT::Node** nodes = Nodes();
            for (int i=0; i<NUMNOD_SOH20; ++i)
            {
              xrefe(i,0) = nodes[i]->X()[0];
              xrefe(i,1) = nodes[i]->X()[1];
              xrefe(i,2) = nodes[i]->X()[2];
            }

            // deformation gradient = identity tensor (geometrically linear case!)
            LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> defgrd(true);
            for (int i=0;i<NUMDIM_SOH20;++i) defgrd(i,i) = 1;

            //----------------------------------------------------------------
            // loop over all Gauss points
            //----------------------------------------------------------------
            for (int gp=0; gp<NUMGPT_SOH20; gp++)
            {
              // Gauss weights and Jacobian determinant
              double fac = detJ_[gp] * weights[gp];

              // Gauss point in reference configuration
              LINALG::Matrix<NUMDIM_SOH20,1> xgp(true);
              for (int k=0;k<NUMDIM_SOH20;++k)
                for (int n=0;n<NUMNOD_SOH20;++n)
                  xgp(k,0) += (vals[gp])(n) * xrefe(n,k);

              //**************************************************************
                  // get analytical solution
                  LINALG::Matrix<NUMDIM_SOH20,1> uanalyt(true);
              LINALG::Matrix<NUMSTR_SOH20,1> strainanalyt(true);
              LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> derivanalyt(true);

              MORTAR::AnalyticalSolutions3D(xgp,uanalyt,strainanalyt,derivanalyt);
              //**************************************************************

                  //--------------------------------------------------------------
                  // (1) L2 norm
                  //--------------------------------------------------------------

                  // compute displacements at GP
                  LINALG::Matrix<NUMDIM_SOH20,1> ugp(true);
              for (int k=0;k<NUMDIM_SOH20;++k)
                for (int n=0;n<NUMNOD_SOH20;++n)
                  ugp(k,0) += (vals[gp])(n) * nodaldisp(NODDOF_SOH20*n+k,0);

              // displacement error
              LINALG::Matrix<NUMDIM_SOH20,1> uerror(true);
              for (int k=0;k<NUMDIM_SOH20;++k)
                uerror(k,0) = uanalyt(k,0) - ugp(k,0);

              // compute GP contribution to L2 error norm
              l2norm += fac * uerror.Dot(uerror);

              //--------------------------------------------------------------
              // (2) H1 norm
              //--------------------------------------------------------------

              // compute derivatives N_XYZ at GP w.r.t. material coordinates
              // by N_XYZ = J^-1 * N_rst
              LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> N_XYZ(true);
              N_XYZ.Multiply(invJ_[gp],derivs[gp]);

              // compute partial derivatives at GP
              LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> derivgp(true);
              for (int l=0;l<NUMDIM_SOH20;++l)
                for (int m=0;m<NUMDIM_SOH20;++m)
                  for (int k=0;k<NUMNOD_SOH20;++k)
                    derivgp(l,m) += N_XYZ(m,k) * nodaldisp(NODDOF_SOH20*k+l,0);

              // derivative error
              LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> deriverror(true);
              for (int k=0;k<NUMDIM_SOH20;++k)
                for (int m=0;m<NUMDIM_SOH20;++m)
                  deriverror(k,m) = derivanalyt(k,m) - derivgp(k,m);

              // compute GP contribution to H1 error norm
              h1norm += fac * deriverror.Dot(deriverror);
              h1norm += fac * uerror.Dot(uerror);

              //--------------------------------------------------------------
              // (3) Energy norm
              //--------------------------------------------------------------

              // compute linear B-operator
              LINALG::Matrix<NUMSTR_SOH20,NUMDOF_SOH20> bop;
              for (int i=0; i<NUMNOD_SOH20; ++i)
              {
                bop(0,NODDOF_SOH20*i+0) = N_XYZ(0,i);
                bop(0,NODDOF_SOH20*i+1) = 0.0;
                bop(0,NODDOF_SOH20*i+2) = 0.0;
                bop(1,NODDOF_SOH20*i+0) = 0.0;
                bop(1,NODDOF_SOH20*i+1) = N_XYZ(1,i);
                bop(1,NODDOF_SOH20*i+2) = 0.0;
                bop(2,NODDOF_SOH20*i+0) = 0.0;
                bop(2,NODDOF_SOH20*i+1) = 0.0;
                bop(2,NODDOF_SOH20*i+2) = N_XYZ(2,i);

                bop(3,NODDOF_SOH20*i+0) = N_XYZ(1,i);
                bop(3,NODDOF_SOH20*i+1) = N_XYZ(0,i);
                bop(3,NODDOF_SOH20*i+2) = 0.0;
                bop(4,NODDOF_SOH20*i+0) = 0.0;
                bop(4,NODDOF_SOH20*i+1) = N_XYZ(2,i);
                bop(4,NODDOF_SOH20*i+2) = N_XYZ(1,i);
                bop(5,NODDOF_SOH20*i+0) = N_XYZ(2,i);
                bop(5,NODDOF_SOH20*i+1) = 0.0;
                bop(5,NODDOF_SOH20*i+2) = N_XYZ(0,i);
              }

              // compute linear strain at GP
              LINALG::Matrix<NUMSTR_SOH20,1> straingp(true);
              straingp.Multiply(bop,nodaldisp);

              // strain error
              LINALG::Matrix<NUMSTR_SOH20,1> strainerror(true);
              for (int k=0;k<NUMSTR_SOH20;++k)
                strainerror(k,0) = strainanalyt(k,0) - straingp(k,0);

              // compute stress vector and constitutive matrix
              double density = 0.0;
              LINALG::Matrix<NUMSTR_SOH20,NUMSTR_SOH20> cmat(true);
              LINALG::Matrix<NUMSTR_SOH20,1> stress(true);
              soh20_mat_sel(&stress,&cmat,&density,&strainerror,&defgrd,gp,params);

              // compute GP contribution to energy error norm
              energynorm += fac * stress.Dot(strainerror);

              //cout << "UAnalytical:      " << ugp << endl;
              //cout << "UDiscrete:        " << uanalyt << endl;
              //cout << "StrainAnalytical: " << strainanalyt << endl;
              //cout << "StrainDiscrete:   " << straingp << endl;
              //cout << "DerivAnalytical:  " << derivanalyt << endl;
              //cout << "DerivDiscrete:    " << derivgp << endl;
            }
            //----------------------------------------------------------------

            // return results
            elevec1_epetra(0) = l2norm;
            elevec1_epetra(1) = h1norm;
            elevec1_epetra(2) = energynorm;
          }
    else
      dserror("ERROR: Error norms only implemented for SVK material");
  }
  break;


  case multi_calc_dens:
  {
    soh20_homog(params);
  }
  break;


  // read restart of microscale
  case multi_readrestart:
  {
    soh20_read_restart_multi();
  }
  break;

  default:
    dserror("Unknown type of action for So_hex20");
  }
  return 0;
}



/*----------------------------------------------------------------------*
 |  Integrate a Volume Neumann boundary condition (public)               |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex20::EvaluateNeumann(ParameterList& params,
                                           DRT::Discretization&      discretization,
                                           DRT::Condition&           condition,
                                           vector<int>&              lm,
                                           Epetra_SerialDenseVector& elevec1,
                                           Epetra_SerialDenseMatrix* elemat1)
{
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >   ("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val"  );

  /*
  **    TIME CURVE BUSINESS
  */
  // find out whether we will use a time curve
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // find out whether we will use a time curve and get the factor
  const vector<int>* curve  = condition.Get<vector<int> >("curve");
  int curvenum = -1;
  if (curve) curvenum = (*curve)[0];
  double curvefac = 1.0;
  if (curvenum>=0 && usetime)
    curvefac = DRT::Problem::Instance()->Curve(curvenum).f(time);
  // **

/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_20 with 20 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH20,1> > shapefcts = soh20_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > derivs = soh20_derivs();
  const static vector<double> gpweights = soh20_weights();
/* ============================================================================*/

  // update element geometry
   LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xrefe;  // material coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH20; ++i){
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];
  }
  /* ================================================= Loop over Gauss Points */
  for (int gp=0; gp<NUMGPT_SOH20; ++gp) {

    // compute the Jacobian matrix
    LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> jac;
    jac.Multiply(derivs[gp],xrefe);

    // compute determinant of Jacobian
    const double detJ = jac.Determinant();
    if (detJ == 0.0) dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ < 0.0) dserror("NEGATIVE JACOBIAN DETERMINANT");

    double fac = gpweights[gp] * curvefac * detJ;          // integration factor
    // distribute/add over element load vector
      for(int dim=0; dim<NUMDIM_SOH20; dim++) {
      double dim_fac = (*onoff)[dim] * (*val)[dim] * fac;
      for (int nodid=0; nodid<NUMNOD_SOH20; ++nodid) {
        elevec1[nodid*NUMDIM_SOH20+dim] += shapefcts[gp](nodid) * dim_fac;
      }
    }

  }/* ==================================================== end of Loop over GP */

  return 0;
} // DRT::ELEMENTS::So_hex20::EvaluateNeumann


/*----------------------------------------------------------------------*
 |  init the element jacobian mapping (protected)                       |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::InitJacobianMapping()
{
  const static vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > derivs = soh20_derivs();
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xrefe;
  for (int i=0; i<NUMNOD_SOH20; ++i)
  {
    xrefe(i,0) = Nodes()[i]->X()[0];
    xrefe(i,1) = Nodes()[i]->X()[1];
    xrefe(i,2) = Nodes()[i]->X()[2];
  }
  invJ_.resize(NUMGPT_SOH20);
  detJ_.resize(NUMGPT_SOH20);
  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
  {
    //invJ_[gp].Shape(NUMDIM_SOH20,NUMDIM_SOH20);
    invJ_[gp].Multiply(derivs[gp],xrefe);
    detJ_[gp] = invJ_[gp].Invert();
    if (detJ_[gp] == 0.0)
      dserror("ZERO JACOBIAN DETERMINANT");
    else if (detJ_[gp] < 0.0)
      dserror("NEGATIVE JACOBIAN DETERMINANT");

  }
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                           popp 09/11 |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_linstiffmass(
    vector<int>&              lm,             // location matrix
    vector<double>&           disp,           // current displacements
    vector<double>&           residual,       // current residual displ
    LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20>* stiffmatrix, // element stiffness matrix
    LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20>* massmatrix,  // element mass matrix
    LINALG::Matrix<NUMDOF_SOH20,1>* force,                 // element internal force vector
    LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20>* elestress,   // stresses at GP
    LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20>* elestrain,   // strains at GP
    ParameterList&            params,         // algorithmic parameters e.g. time
    const INPAR::STR::StressType   iostress,  // stress output option
    const INPAR::STR::StrainType   iostrain)  // strain output option
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_20 with 27 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH20,1> > shapefcts = soh20_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > derivs = soh20_derivs();
  const static vector<double> gpweights = soh20_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xcurr;  // current  coord. of element
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xdisp;

  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH20; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH20+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH20+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH20+2];
  }

  LINALG::Matrix<NUMDOF_SOH20,1> nodaldisp;
  for (int i=0; i<NUMDOF_SOH20; ++i)
  {
    nodaldisp(i,0) = disp[i];
  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  // CAUTION: defgrd(true): filled with zeros!
  LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> defgrd(true);
  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
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

    // set to initial state as test to receive a linear solution
    for (int i=0; i<3; ++i) defgrd(i,i) = 1.0;

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<NUMSTR_SOH20,NUMDOF_SOH20> bop;
    for (int i=0; i<NUMNOD_SOH20; ++i)
    {
      bop(0,NODDOF_SOH20*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH20*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH20*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_SOH20*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH20*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH20*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_SOH20*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH20*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH20*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_SOH20*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH20*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH20*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_SOH20*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH20*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH20*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_SOH20*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH20*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH20*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    // now build the linear strain
    LINALG::Matrix<NUMSTR_SOH20,1> strainlin(true);
    strainlin.Multiply(bop,nodaldisp);

    // and rename it as glstrain to use the common methods further on

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH20);
    LINALG::Matrix<NUMSTR_SOH20,1> glstrain(glstrain_epetra.A(),true);
    glstrain.Update(1.0,strainlin);

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
      // rewriting Green-Lagrange strains in matrix format
      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> gl;
      gl(0,0) = glstrain(0);
      gl(0,1) = 0.5*glstrain(3);
      gl(0,2) = 0.5*glstrain(5);
      gl(1,0) = gl(0,1);
      gl(1,1) = glstrain(1);
      gl(1,2) = 0.5*glstrain(4);
      gl(2,0) = gl(0,2);
      gl(2,1) = gl(1,2);
      gl(2,2) = glstrain(2);

      // inverse of deformation gradient
      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> invdefgrd;
      invdefgrd.Invert(defgrd);

      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> temp;
      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> euler_almansi;
      temp.Multiply(gl,invdefgrd);
      euler_almansi.MultiplyTN(invdefgrd,temp);

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
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOH20,NUMSTR_SOH20> cmat(true);
    LINALG::Matrix<NUMSTR_SOH20,1> stress(true);
    LINALG::Matrix<NUMSTR_SOH20,1> plglstrain(true);
    soh20_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH20; ++i)
        (*elestress)(gp,i) = stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = stress(0);
      pkstress(0,1) = stress(3);
      pkstress(0,2) = stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = stress(1);
      pkstress(1,2) = stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = stress(2);

      LINALG::Matrix<3,3> temp;
      LINALG::Matrix<3,3> cauchystress;
      temp.Multiply(1.0/detF,defgrd,pkstress);
      cauchystress.MultiplyNT(temp,defgrd);

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
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }
    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH20> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);
    }

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH20; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH20; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH20*inod+0,NUMDIM_SOH20*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH20*inod+1,NUMDIM_SOH20*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH20*inod+2,NUMDIM_SOH20*jnod+2) += massfactor;
        }
      }

    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++
   /* =========================================================================*/
  }/* ==================================================== end of Loop over GP */
   /* =========================================================================*/

  return;
} // DRT::ELEMENTS::So_hex20::soh20_linstiffmass


/*----------------------------------------------------------------------*
 |  evaluate the element (private)                                      |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_nlnstiffmass(
      vector<int>&              lm,             // location matrix
      vector<double>&           disp,           // current displacements
      vector<double>&           residual,       // current residual displ
      LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20>* stiffmatrix, // element stiffness matrix
      LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20>* massmatrix,  // element mass matrix
      LINALG::Matrix<NUMDOF_SOH20,1>* force,                 // element internal force vector
      LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20>* elestress,   // stresses at GP
      LINALG::Matrix<NUMGPT_SOH20,NUMSTR_SOH20>* elestrain,   // strains at GP
      ParameterList&            params,         // algorithmic parameters e.g. time
      const INPAR::STR::StressType   iostress,  // stress output option
      const INPAR::STR::StrainType   iostrain)  // strain output option
{
/* ============================================================================*
** CONST SHAPE FUNCTIONS, DERIVATIVES and WEIGHTS for HEX_20 with 27 GAUSS POINTS*
** ============================================================================*/
  const static vector<LINALG::Matrix<NUMNOD_SOH20,1> > shapefcts = soh20_shapefcts();
  const static vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > derivs = soh20_derivs();
  const static vector<double> gpweights = soh20_weights();
/* ============================================================================*/

  // update element geometry
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xrefe;  // material coord. of element
  LINALG::Matrix<NUMNOD_SOH20,NUMDIM_SOH20> xcurr;  // current  coord. of element
  DRT::Node** nodes = Nodes();
  for (int i=0; i<NUMNOD_SOH20; ++i)
  {
    const double* x = nodes[i]->X();
    xrefe(i,0) = x[0];
    xrefe(i,1) = x[1];
    xrefe(i,2) = x[2];

    xcurr(i,0) = xrefe(i,0) + disp[i*NODDOF_SOH20+0];
    xcurr(i,1) = xrefe(i,1) + disp[i*NODDOF_SOH20+1];
    xcurr(i,2) = xrefe(i,2) + disp[i*NODDOF_SOH20+2];

  }

  /* =========================================================================*/
  /* ================================================= Loop over Gauss Points */
  /* =========================================================================*/
  LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> N_XYZ;
  // build deformation gradient wrt to material configuration
  // in case of prestressing, build defgrd wrt to last stored configuration
  LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> defgrd(false);
  for (int gp=0; gp<NUMGPT_SOH20; ++gp)
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

    // (material) deformation gradient F = d xcurr / d xrefe = xcurr^T * N_XYZ^T
    defgrd.MultiplyTT(xcurr,N_XYZ);

    // Right Cauchy-Green tensor = F^T * F
    LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> cauchygreen;
    cauchygreen.MultiplyTN(defgrd,defgrd);

    // Green-Lagrange strains matrix E = 0.5 * (Cauchygreen - Identity)
    // GL strain vector glstrain={E11,E22,E33,2*E12,2*E23,2*E31}
    Epetra_SerialDenseVector glstrain_epetra(NUMSTR_SOH20);
    LINALG::Matrix<NUMSTR_SOH20,1> glstrain(glstrain_epetra.A(),true);
    glstrain(0) = 0.5 * (cauchygreen(0,0) - 1.0);
    glstrain(1) = 0.5 * (cauchygreen(1,1) - 1.0);
    glstrain(2) = 0.5 * (cauchygreen(2,2) - 1.0);
    glstrain(3) = cauchygreen(0,1);
    glstrain(4) = cauchygreen(1,2);
    glstrain(5) = cauchygreen(2,0);

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
      // rewriting Green-Lagrange strains in matrix format
      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> gl;
      gl(0,0) = glstrain(0);
      gl(0,1) = 0.5*glstrain(3);
      gl(0,2) = 0.5*glstrain(5);
      gl(1,0) = gl(0,1);
      gl(1,1) = glstrain(1);
      gl(1,2) = 0.5*glstrain(4);
      gl(2,0) = gl(0,2);
      gl(2,1) = gl(1,2);
      gl(2,2) = glstrain(2);

      // inverse of deformation gradient
      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> invdefgrd;
      invdefgrd.Invert(defgrd);

      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> temp;
      LINALG::Matrix<NUMDIM_SOH20,NUMDIM_SOH20> euler_almansi;
      temp.Multiply(gl,invdefgrd);
      euler_almansi.MultiplyTN(invdefgrd,temp);

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
    }

    /* non-linear B-operator (may so be called, meaning
    ** of B-operator is not so sharp in the non-linear realm) *
    ** B = F . Bl *
    **
    **      [ ... | F_11*N_{,1}^k  F_21*N_{,1}^k  F_31*N_{,1}^k | ... ]
    **      [ ... | F_12*N_{,2}^k  F_22*N_{,2}^k  F_32*N_{,2}^k | ... ]
    **      [ ... | F_13*N_{,3}^k  F_23*N_{,3}^k  F_33*N_{,3}^k | ... ]
    ** B =  [ ~~~   ~~~~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~~   ~~~ ]
    **      [       F_11*N_{,2}^k+F_12*N_{,1}^k                       ]
    **      [ ... |          F_21*N_{,2}^k+F_22*N_{,1}^k        | ... ]
    **      [                       F_31*N_{,2}^k+F_32*N_{,1}^k       ]
    **      [                                                         ]
    **      [       F_12*N_{,3}^k+F_13*N_{,2}^k                       ]
    **      [ ... |          F_22*N_{,3}^k+F_23*N_{,2}^k        | ... ]
    **      [                       F_32*N_{,3}^k+F_33*N_{,2}^k       ]
    **      [                                                         ]
    **      [       F_13*N_{,1}^k+F_11*N_{,3}^k                       ]
    **      [ ... |          F_23*N_{,1}^k+F_21*N_{,3}^k        | ... ]
    **      [                       F_33*N_{,1}^k+F_31*N_{,3}^k       ]
    */
    LINALG::Matrix<NUMSTR_SOH20,NUMDOF_SOH20> bop;
    for (int i=0; i<NUMNOD_SOH20; ++i)
    {
      bop(0,NODDOF_SOH20*i+0) = defgrd(0,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH20*i+1) = defgrd(1,0)*N_XYZ(0,i);
      bop(0,NODDOF_SOH20*i+2) = defgrd(2,0)*N_XYZ(0,i);
      bop(1,NODDOF_SOH20*i+0) = defgrd(0,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH20*i+1) = defgrd(1,1)*N_XYZ(1,i);
      bop(1,NODDOF_SOH20*i+2) = defgrd(2,1)*N_XYZ(1,i);
      bop(2,NODDOF_SOH20*i+0) = defgrd(0,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH20*i+1) = defgrd(1,2)*N_XYZ(2,i);
      bop(2,NODDOF_SOH20*i+2) = defgrd(2,2)*N_XYZ(2,i);
      /* ~~~ */
      bop(3,NODDOF_SOH20*i+0) = defgrd(0,0)*N_XYZ(1,i) + defgrd(0,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH20*i+1) = defgrd(1,0)*N_XYZ(1,i) + defgrd(1,1)*N_XYZ(0,i);
      bop(3,NODDOF_SOH20*i+2) = defgrd(2,0)*N_XYZ(1,i) + defgrd(2,1)*N_XYZ(0,i);
      bop(4,NODDOF_SOH20*i+0) = defgrd(0,1)*N_XYZ(2,i) + defgrd(0,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH20*i+1) = defgrd(1,1)*N_XYZ(2,i) + defgrd(1,2)*N_XYZ(1,i);
      bop(4,NODDOF_SOH20*i+2) = defgrd(2,1)*N_XYZ(2,i) + defgrd(2,2)*N_XYZ(1,i);
      bop(5,NODDOF_SOH20*i+0) = defgrd(0,2)*N_XYZ(0,i) + defgrd(0,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH20*i+1) = defgrd(1,2)*N_XYZ(0,i) + defgrd(1,0)*N_XYZ(2,i);
      bop(5,NODDOF_SOH20*i+2) = defgrd(2,2)*N_XYZ(0,i) + defgrd(2,0)*N_XYZ(2,i);
    }

    /* call material law cccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ** Here all possible material laws need to be incorporated,
    ** the stress vector, a C-matrix, and a density must be retrieved,
    ** every necessary data must be passed.
    */
    double density = 0.0;
    LINALG::Matrix<NUMSTR_SOH20,NUMSTR_SOH20> cmat(true);
    LINALG::Matrix<NUMSTR_SOH20,1> stress(true);
    soh20_mat_sel(&stress,&cmat,&density,&glstrain,&defgrd,gp,params);
    // end of call material law ccccccccccccccccccccccccccccccccccccccccccccccc

    // return gp stresses
    switch (iostress)
    {
    case INPAR::STR::stress_2pk:
    {
      if (elestress == NULL) dserror("stress data not available");
      for (int i = 0; i < NUMSTR_SOH20; ++i)
        (*elestress)(gp,i) = stress(i);
    }
    break;
    case INPAR::STR::stress_cauchy:
    {
      if (elestress == NULL) dserror("stress data not available");
      const double detF = defgrd.Determinant();

      LINALG::Matrix<3,3> pkstress;
      pkstress(0,0) = stress(0);
      pkstress(0,1) = stress(3);
      pkstress(0,2) = stress(5);
      pkstress(1,0) = pkstress(0,1);
      pkstress(1,1) = stress(1);
      pkstress(1,2) = stress(4);
      pkstress(2,0) = pkstress(0,2);
      pkstress(2,1) = pkstress(1,2);
      pkstress(2,2) = stress(2);

      LINALG::Matrix<3,3> temp;
      LINALG::Matrix<3,3> cauchystress;
      temp.Multiply(1.0/detF,defgrd,pkstress);
      cauchystress.MultiplyNT(temp,defgrd);

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
    // update internal force vector
    if (force != NULL)
    {
      // integrate internal force vector f = f + (B^T . sigma) * detJ * w(gp)
      force->MultiplyTN(detJ_w, bop, stress, 1.0);
    }
    // update stiffness matrix
    if (stiffmatrix != NULL)
    {
      // integrate `elastic' and `initial-displacement' stiffness matrix
      // keu = keu + (B^T . C . B) * detJ * w(gp)
      LINALG::Matrix<6,NUMDOF_SOH20> cb;
      cb.Multiply(cmat,bop);
      stiffmatrix->MultiplyTN(detJ_w,bop,cb,1.0);

      // integrate `geometric' stiffness matrix and add to keu *****************
      LINALG::Matrix<6,1> sfac(stress); // auxiliary integrated stress
      sfac.Scale(detJ_w); // detJ*w(gp)*[S11,S22,S33,S12=S21,S23=S32,S13=S31]
      vector<double> SmB_L(3); // intermediate Sm.B_L
      // kgeo += (B_L^T . sigma . B_L) * detJ * w(gp)  with B_L = Ni,Xj see NiliFEM-Skript
      for (int inod=0; inod<NUMNOD_SOH20; ++inod) {
        SmB_L[0] = sfac(0) * N_XYZ(0, inod) + sfac(3) * N_XYZ(1, inod)
            + sfac(5) * N_XYZ(2, inod);
        SmB_L[1] = sfac(3) * N_XYZ(0, inod) + sfac(1) * N_XYZ(1, inod)
            + sfac(4) * N_XYZ(2, inod);
        SmB_L[2] = sfac(5) * N_XYZ(0, inod) + sfac(4) * N_XYZ(1, inod)
            + sfac(2) * N_XYZ(2, inod);
        for (int jnod=0; jnod<NUMNOD_SOH20; ++jnod) {
          double bopstrbop = 0.0; // intermediate value
          for (int idim=0; idim<NUMDIM_SOH20; ++idim)
            bopstrbop += N_XYZ(idim, jnod) * SmB_L[idim];
          (*stiffmatrix)(3*inod+0,3*jnod+0) += bopstrbop;
          (*stiffmatrix)(3*inod+1,3*jnod+1) += bopstrbop;
          (*stiffmatrix)(3*inod+2,3*jnod+2) += bopstrbop;
        }
      } // end of integrate `geometric' stiffness******************************
    }

    if (massmatrix != NULL) // evaluate mass matrix +++++++++++++++++++++++++
    {
      // integrate consistent mass matrix
      const double factor = detJ_w * density;
      double ifactor, massfactor;
      for (int inod=0; inod<NUMNOD_SOH20; ++inod)
      {
        ifactor = shapefcts[gp](inod) * factor;
        for (int jnod=0; jnod<NUMNOD_SOH20; ++jnod)
        {
          massfactor = shapefcts[gp](jnod) * ifactor;     // intermediate factor
          (*massmatrix)(NUMDIM_SOH20*inod+0,NUMDIM_SOH20*jnod+0) += massfactor;
          (*massmatrix)(NUMDIM_SOH20*inod+1,NUMDIM_SOH20*jnod+1) += massfactor;
          (*massmatrix)(NUMDIM_SOH20*inod+2,NUMDIM_SOH20*jnod+2) += massfactor;
        }
      }

    } // end of mass matrix +++++++++++++++++++++++++++++++++++++++++++++++++++

  }/* ==================================================== end of Loop over GP */

  return;
} // DRT::ELEMENTS::So_hex20::soh20_nlnstiffmass

/*----------------------------------------------------------------------*
 |  lump mass matrix (private)                                          |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_lumpmass(LINALG::Matrix<NUMDOF_SOH20,NUMDOF_SOH20>* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (unsigned int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (unsigned int r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate Hex20 Shape fcts at all 20 Gauss Points                     |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMNOD_SOH20,1> > DRT::ELEMENTS::So_hex20::soh20_shapefcts()
{
  vector<LINALG::Matrix<NUMNOD_SOH20,1> > shapefcts(NUMGPT_SOH20);
  // (r,s,t) gp-locations of fully integrated quadratic Hex 20
  // fill up nodal f at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp)
  {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    DRT::UTILS::shape_function_3D(shapefcts[igp], r, s, t, hex20);
  }
  return shapefcts;
}


/*----------------------------------------------------------------------*
 |  Evaluate Hex20 Shape fct derivs at all 20 Gauss Points              |
 *----------------------------------------------------------------------*/
const vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > DRT::ELEMENTS::So_hex20::soh20_derivs()
{
  vector<LINALG::Matrix<NUMDIM_SOH20,NUMNOD_SOH20> > derivs(NUMGPT_SOH20);
  // (r,s,t) gp-locations of fully integrated quadratic Hex 20
  // fill up df w.r.t. rst directions (NUMDIM) at each gp
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int igp = 0; igp < intpoints.nquad; ++igp)
  {
    const double r = intpoints.qxg[igp][0];
    const double s = intpoints.qxg[igp][1];
    const double t = intpoints.qxg[igp][2];

    DRT::UTILS::shape_function_3D_deriv1(derivs[igp], r, s, t, hex20);
  }
  return derivs;
}

/*----------------------------------------------------------------------*
 |  Evaluate Hex20 Weights at all 20 Gauss Points                       |
 *----------------------------------------------------------------------*/
const vector<double> DRT::ELEMENTS::So_hex20::soh20_weights()
{
  vector<double> weights(NUMGPT_SOH20);
  const DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule_hex_27point;
  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);
  for (int i = 0; i < NUMGPT_SOH20; ++i)
  {
    weights[i] = intpoints.qwgt[i];
  }
  return weights;
}

/*----------------------------------------------------------------------*
 |  shape functions and derivatives for So_hex20                         |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::So_hex20::soh20_shapederiv(
      LINALG::Matrix<NUMNOD_SOH20,NUMGPT_SOH20>** shapefct,   // pointer to pointer of shapefct
      LINALG::Matrix<NUMDOF_SOH20,NUMNOD_SOH20>** deriv,     // pointer to pointer of derivs
      LINALG::Matrix<NUMGPT_SOH20,1>** weights)   // pointer to pointer of weights
{
  // static matrix objects, kept in memory
  static LINALG::Matrix<NUMNOD_SOH20,NUMGPT_SOH20>  f;  // shape functions
  static LINALG::Matrix<NUMDOF_SOH20,NUMNOD_SOH20> df;  // derivatives
  static LINALG::Matrix<NUMGPT_SOH20,1> weightfactors;   // weights for each gp
  static bool fdf_eval;                      // flag for re-evaluate everything

  if (fdf_eval==true) // if true f,df already evaluated
  {
    *shapefct = &f; // return adress of static object to target of pointer
    *deriv = &df; // return adress of static object to target of pointer
    *weights = &weightfactors; // return adress of static object to target of pointer
    return;
  }
  else
  {
    // (r,s,t) gp-locations of fully integrated quadratic Hex 20
    // fill up nodal f at each gp
    // fill up df w.r.t. rst directions (NUMDIM) at each gp
    const DRT::UTILS::GaussRule3D gaussrule_ = DRT::UTILS::intrule_hex_27point;
    const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule_);
    for (int igp = 0; igp < intpoints.nquad; ++igp)
    {
      const double r = intpoints.qxg[igp][0];
      const double s = intpoints.qxg[igp][1];
      const double t = intpoints.qxg[igp][2];

      LINALG::Matrix<NUMNOD_SOH20,1> funct;
      LINALG::Matrix<NUMDIM_SOH20, NUMNOD_SOH20> deriv;
      DRT::UTILS::shape_function_3D(funct, r, s, t, hex20);
      DRT::UTILS::shape_function_3D_deriv1(deriv, r, s, t, hex20);
      for (int inode = 0; inode < NUMNOD_SOH20; ++inode)
      {
        f(inode, igp) = funct(inode);
        df(igp*NUMDIM_SOH20+0, inode) = deriv(0, inode);
        df(igp*NUMDIM_SOH20+1, inode) = deriv(1, inode);
        df(igp*NUMDIM_SOH20+2, inode) = deriv(2, inode);
        weightfactors(igp) = intpoints.qwgt[igp];
      }
    }
    // return adresses of just evaluated matrices
    *shapefct = &f; // return adress of static object to target of pointer
    *deriv = &df; // return adress of static object to target of pointer
    *weights = &weightfactors; // return adress of static object to target of pointer
    fdf_eval = true; // now all arrays are filled statically
  }
  return;
}  // of soh20_shapederiv


/*----------------------------------------------------------------------*
 |  init the element (public)                                           |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::So_hex20Type::Initialize(DRT::Discretization& dis)
{
  for (int i=0; i<dis.NumMyColElements(); ++i)
  {
    if (dis.lColElement(i)->ElementType() != *this) continue;
    DRT::ELEMENTS::So_hex20* actele = dynamic_cast<DRT::ELEMENTS::So_hex20*>(dis.lColElement(i));
    if (!actele) dserror("cast to So_hex20* failed");
    actele->InitJacobianMapping();
  }
  return 0;
}



