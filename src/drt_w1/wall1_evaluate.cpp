/*----------------------------------------------------------------------*/
/*!
\file wall1_evaluate.cpp

\brief This file contains the main evaluate methods of a wall1 element.

\level 1

\maintainer Michael Hiermeier

*/
/*----------------------------------------------------------------------*/
// macros


/*----------------------------------------------------------------------*/
// headers
#include "wall1.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "Epetra_SerialDenseSolver.h"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_contact/contact_analytical.H"
#include "../drt_mat/stvenantkirchhoff.H"

#include "../drt_structure_new/str_elements_paramsinterface.H"

/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            mwgee 12/06|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Wall1::Evaluate(Teuchos::ParameterList&   params,
                                   DRT::Discretization&      discretization,
                                   std::vector<int>&         lm,
                                   Epetra_SerialDenseMatrix& elemat1,
                                   Epetra_SerialDenseMatrix& elemat2,
                                   Epetra_SerialDenseVector& elevec1,
                                   Epetra_SerialDenseVector& elevec2,
                                   Epetra_SerialDenseVector& elevec3)
{
  SetParamsInterfacePtr(params);
  ELEMENTS::ActionType act = ELEMENTS::none;

  if (IsParamsInterface())
  {
    act = ParamsInterface().GetActionType();
  }
  else
  {
    // get the action required
    std::string action = params.get<std::string>("action","calc_none");
    if (action == "calc_none") dserror("No action supplied");
    else if (action=="calc_struct_linstiff")      act = ELEMENTS::struct_calc_linstiff;
    else if (action=="calc_struct_nlnstiff")      act = ELEMENTS::struct_calc_nlnstiff;
    else if (action=="calc_struct_internalforce") act = ELEMENTS::struct_calc_internalforce;
    else if (action=="calc_struct_linstiffmass")  act = ELEMENTS::struct_calc_linstiffmass;
    else if (action=="calc_struct_nlnstiffmass")  act = ELEMENTS::struct_calc_nlnstiffmass;
    else if (action=="calc_struct_nlnstifflmass") act = ELEMENTS::struct_calc_nlnstifflmass;
    else if (action=="calc_struct_nlnstiff_gemm") act = ELEMENTS::struct_calc_nlnstiff_gemm;
    else if (action=="calc_struct_stress")        act = ELEMENTS::struct_calc_stress;
    else if (action=="postprocess_stress")        act = ELEMENTS::struct_postprocess_stress;
    else if (action=="calc_struct_eleload")       act = ELEMENTS::struct_calc_eleload;
    else if (action=="calc_struct_fsiload")       act = ELEMENTS::struct_calc_fsiload;
    else if (action=="calc_struct_update_istep")  act = ELEMENTS::struct_calc_update_istep;
    else if (action=="calc_struct_reset_istep")   act = ELEMENTS::struct_calc_reset_istep;
    else if (action=="calc_struct_energy")        act = ELEMENTS::struct_calc_energy;
    else if (action=="calc_struct_errornorms")    act = ELEMENTS::struct_calc_errornorms;
    else if (action=="calc_struct_mass_volume")   act = ELEMENTS::struct_calc_mass_volume;
    else dserror("Unknown type of action %s for Wall1", action.c_str());
  }
  // get the material law
  Teuchos::RCP<const MAT::Material> actmat = Material();

  // --------------------------------------------------
  // Now do the nurbs specific stuff
  std::vector<Epetra_SerialDenseVector> myknots(2);

  if(Shape()==DRT::Element::nurbs4 or
     Shape()==DRT::Element::nurbs9)
  {
    switch(act)
    {
      case ELEMENTS::struct_calc_linstiff:
      case ELEMENTS::struct_calc_nlnstiffmass:
      case ELEMENTS::struct_calc_nlnstifflmass:
      case ELEMENTS::struct_calc_nlnstiff:
      case ELEMENTS::struct_calc_internalforce:
      case ELEMENTS::struct_calc_stress:
      case ELEMENTS::struct_calc_errornorms:
      case ELEMENTS::struct_calc_mass_volume:
      {
        DRT::NURBS::NurbsDiscretization* nurbsdis =
            dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

        bool zero_sized=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());

        // skip zero sized elements in knot span --- they correspond to interpolated nodes
        if(zero_sized)
          return(0);

        break;
      }
      default :
        myknots.clear();
        break;
    }
  }

  switch(act)
  {
    //==================================================================================
    case ELEMENTS::struct_calc_linstiff:
    {
      // need current displacement and residual forces
      std::vector<double> mydisp(lm.size());
      for (int i=0; i<(int)mydisp.size(); ++i) mydisp[i] = 0.0;
      std::vector<double> myres(lm.size());
      for (int i=0; i<(int)myres.size(); ++i) myres[i] = 0.0;
      std::vector<double> mydispmat(lm.size());
      for (int i=0; i<(int)mydispmat.size(); ++i) mydispmat[i] = 0.0;

      // special case: geometrically linear
      if (kintype_ == INPAR::STR::kinem_linear)
      {
        w1_linstiffmass(lm,mydisp,myres,mydispmat,myknots,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        w1_nlnstiffmass(lm,mydisp,myres,mydispmat,myknots,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_nlnstiffmass:
    case ELEMENTS::struct_calc_nlnstifflmass:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      std::vector<double> mydispmat(lm.size());
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat = discretization.GetState("material_displacement");
        DRT::UTILS::ExtractMyValues(*dispmat,mydispmat,lm);
      }

      // special case: geometrically linear
      if (kintype_ == INPAR::STR::kinem_linear)
      {
        w1_linstiffmass(lm,mydisp,mydispmat,myres,myknots,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        w1_nlnstiffmass(lm,mydisp,mydispmat,myres,myknots,&elemat1,&elemat2,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }

      if (act==ELEMENTS::struct_calc_nlnstifflmass) w1_lumpmass(&elemat2);
      break;
    }
    //==================================================================================
    // NULL-pointer for mass matrix in case of calculating only stiff matrix
    case ELEMENTS::struct_calc_nlnstiff:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      std::vector<double> mydispmat(lm.size());
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat = discretization.GetState("material_displacement");;
        DRT::UTILS::ExtractMyValues(*dispmat,mydispmat,lm);
      }

      // special case: geometrically linear
      if (kintype_ == INPAR::STR::kinem_linear)
      {
        w1_linstiffmass(lm,mydisp,myres,mydispmat,myknots,&elemat1,NULL,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        w1_nlnstiffmass(lm,mydisp,myres,mydispmat,myknots,&elemat1,NULL,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_internalforce:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null) dserror("Cannot get state vectors 'displacement' and/or residual");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      // create a dummy element matrix (initialised to zero)
      // This matrix is not utterly useless. It is used to apply EAS-stuff in a linearised manner
      // onto the internal force vector.
      Epetra_SerialDenseMatrix myemat(lm.size(),lm.size());
      std::vector<double> mydispmat(lm.size());
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat = discretization.GetState("material_displacement");;
        DRT::UTILS::ExtractMyValues(*dispmat,mydispmat,lm);
      }

      // special case: geometrically linear
      if (kintype_ == INPAR::STR::kinem_linear)
      {
        w1_linstiffmass(lm,mydisp,myres,mydispmat,myknots,&myemat,NULL,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      // standard is: geometrically non-linear with Total Lagrangean approach
      else
      {
        w1_nlnstiffmass(lm,mydisp,myres,mydispmat,myknots,&myemat,NULL,&elevec1,NULL,NULL,actmat,params,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      }
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_nlnstiff_gemm:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> dispo = discretization.GetState("old displacement");
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
      if (dispo==Teuchos::null or disp==Teuchos::null or res==Teuchos::null)
        dserror("Cannot get state vectors");
      std::vector<double> mydispo(lm.size());
      DRT::UTILS::ExtractMyValues(*dispo,mydispo,lm);
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      FintStiffMassGEMM(params,lm,mydispo,mydisp,myres,&elemat1,NULL,&elevec1,NULL,NULL,actmat,
                        INPAR::STR::stress_none,INPAR::STR::strain_none);
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_recover:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp =
          discretization.GetState("displacement");
      Teuchos::RCP<const Epetra_Vector> res  =
          discretization.GetState("residual displacement");
      if (disp==Teuchos::null || res==Teuchos::null)
        dserror("Cannot get state vectors \"displacement\" "
          "and/or \"residual displacement\"");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      std::vector<double> myres(lm.size());
      DRT::UTILS::ExtractMyValues(*res,myres,lm);
      w1_recover(lm,mydisp,myres);
      /* ToDo Probably we have to recover the history information of some special
       * materials as well.                                 hiermeier 04/2016  */
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_update_istep:
    {
      // do something with internal EAS, etc parameters
      if (iseas_)
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS blas;  // BLAS front-end dummy
        blas.COPY((*alphao).M()*(*alphao).N(), (*alpha).A(), (*alphao).A());  // alphao := alpha
      }
      SolidMaterial()->Update();
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_reset_istep:
    {
      // do something with internal EAS, etc parameters
      if (iseas_)
      {
        Epetra_SerialDenseMatrix* alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");  // Alpha_{n+1}
        Epetra_SerialDenseMatrix* alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");  // Alpha_n
        Epetra_BLAS blas;
        blas.COPY((*alphao).M()*(*alphao).N(), (*alphao).A(), (*alpha).A());  // alpha := alphao
      }
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_stress:
    {
      // nothing to do for ghost elements
      if (discretization.Comm().MyPID()==Owner())
      {
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
        Teuchos::RCP<std::vector<char> > stressdata = Teuchos::null;
        Teuchos::RCP<std::vector<char> > straindata = Teuchos::null;
        INPAR::STR::StressType iostress = INPAR::STR::stress_none;
        INPAR::STR::StrainType iostrain = INPAR::STR::strain_none;
        if (IsParamsInterface())
        {
          stressdata   = StrParamsInterface().MutableStressDataPtr();
          straindata   = StrParamsInterface().MutableStrainDataPtr();

          iostress   = StrParamsInterface().GetStressOutputType();
          iostrain   = StrParamsInterface().GetStrainOutputType();
        }
        else
        {
          stressdata = params.get<Teuchos::RCP<std::vector<char> > >("stress",Teuchos::null);
          straindata = params.get<Teuchos::RCP<std::vector<char> > >("strain",Teuchos::null);
          iostress = DRT::INPUT::get<INPAR::STR::StressType>(params, "iostress", INPAR::STR::stress_none);
          iostrain = DRT::INPUT::get<INPAR::STR::StrainType>(params, "iostrain", INPAR::STR::strain_none);
        }
        if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
        if (stressdata==Teuchos::null) dserror("Cannot get stress 'data'");
        if (straindata==Teuchos::null) dserror("Cannot get strain 'data'");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        std::vector<double> myres(lm.size());
        DRT::UTILS::ExtractMyValues(*res,myres,lm);
        std::vector<double> mydispmat(lm.size());
        if (structale_)
        {
          Teuchos::RCP<const Epetra_Vector> dispmat = discretization.GetState("material_displacement");;
          DRT::UTILS::ExtractMyValues(*dispmat,mydispmat,lm);
        }
        const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
        Epetra_SerialDenseMatrix stress(intpoints.nquad,Wall1::numstr_);
        Epetra_SerialDenseMatrix strain(intpoints.nquad,Wall1::numstr_);

        // special case: geometrically linear
        if (kintype_ == INPAR::STR::kinem_linear)
        {
          w1_linstiffmass(lm,mydisp,myres,mydispmat,myknots,NULL,NULL,NULL,&stress,&strain,actmat,params,iostress,iostrain);
        }
        // standard is: geometrically non-linear with Total Lagrangean approach
        else
        {
          w1_nlnstiffmass(lm,mydisp,myres,mydispmat,myknots,NULL,NULL,NULL,&stress,&strain,actmat,params,iostress,iostrain);
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
      break;
    }
    //==================================================================================
    // postprocess stresses/strains at gauss points
    /* note that in the following, quantities are always referred to as
     * "stresses" etc. although they might also apply to strains
     * (depending on what this routine is called for from the post filter) */
    case ELEMENTS::struct_postprocess_stress:
    {
      std::string groupname = params.get<std::string>("groupname","gauss_2PK_stresses_xyz");

      //coupling stress are evaluated in the respective coupling element
      //if( groupname != "gauss_2PK_coupling_stresses_xyz" and groupname!= "gauss_cauchy_coupling_stresses_xyz")
      {
        const Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > gpstressmap=
          params.get<Teuchos::RCP<std::map<int,Teuchos::RCP<Epetra_SerialDenseMatrix> > > >("gpstressmap",Teuchos::null);
        if (gpstressmap==Teuchos::null)
          dserror("no gp stress/strain map available for postprocessing");
        std::string stresstype = params.get<std::string>("stresstype","ndxyz");
        int gid = Id();
        Teuchos::RCP<Epetra_SerialDenseMatrix> gpstress = (*gpstressmap)[gid];
        Teuchos::RCP<Epetra_MultiVector> poststress=params.get<Teuchos::RCP<Epetra_MultiVector> >("poststress",Teuchos::null);
        if (poststress==Teuchos::null)
          dserror("No element stress/strain vector available");

        if (stresstype=="ndxyz")
        {
          // extrapolate stresses/strains at Gauss points to nodes
          w1_expol(*gpstress, *poststress);
        }
        else if (stresstype=="cxyz")
        {
          const Epetra_BlockMap& elemap = poststress->Map();
          int lid = elemap.LID(Id());
          const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);
          if (lid!=-1)
          {
            // maximum 4 independent stresses exist in 2D
            for (int i = 0; i < Wall1::numstr_; ++i)
            {
              (*((*poststress)(i)))[lid] = 0.;
              for (int j = 0; j < intpoints.nquad; ++j)
              {
                (*((*poststress)(i)))[lid] += 1.0/intpoints.nquad * (*gpstress)(j,i);
              }
            }
          }
        }
        else
        {
          dserror("unknown type of stress/strain output on element level");
        }
      }
      break;
    }
    case ELEMENTS::struct_calc_energy:
    {
      // need current displacement and residual forces
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state vectors");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      // check if length suffices
      if (elevec1.Length() < 1) dserror("Result vector too short");
      // determine energies
      Energy(params,lm,mydisp,&elevec1,actmat);
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_errornorms:
    {
      // IMPORTANT NOTES (popp 10/2010):
      // - error norms are based on a small deformation assumption (linear elasticity)
      // - extension to finite deformations would be possible without difficulties,
      //   however analytical solutions are extremely rare in the nonlinear realm
      // - only implemented for SVK material (relevant for energy norm only, L2 and
      //   H1 norms are of course valid for arbitrary materials)
      // - analytical solutions are currently stored in a repository in the CONTACT
      //   namespace, however they could (should?) be moved to a more general location

      // check length of elevec1
      if (elevec1.Length() < 3) dserror("The given result vector is too short.");

      // check material law
      Teuchos::RCP<MAT::Material> mat = Material();

      //******************************************************************
      // only for St.Venant Kirchhoff material
      //******************************************************************
      if (mat->MaterialType() == INPAR::MAT::m_stvenant)
      {
        // declaration of variables
        double l2norm = 0.0;
        double h1norm = 0.0;
        double energynorm = 0.0;

        // some definitions
        const int numnode = NumNode();
        const int numdf   = 2;
        const int nd      = numnode*numdf;
        const int numeps  = 4;
        Epetra_SerialDenseMatrix xjm;
        xjm.Shape(2,2);
        double det = 0.0;
        Epetra_SerialDenseMatrix boplin;
        boplin.Shape(numeps,nd);
        Epetra_SerialDenseVector F;
        F.Size(numeps);
        Epetra_SerialDenseVector strain;
        strain.Size(numeps);

        // shape functions, derivatives and integration rule
        Epetra_SerialDenseVector funct(numnode);
        Epetra_SerialDenseMatrix deriv;
        deriv.Shape(2,numnode);
        const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

        // get displacements and extract values of this element
        Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==Teuchos::null) dserror("Cannot get state displacement vector");
        std::vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

        // reference and current geometry (nodal positions)
        Epetra_SerialDenseMatrix xrefe(2,numnode);
        Epetra_SerialDenseMatrix xcure(2,numnode);
        for (int k=0; k<numnode; ++k)
        {
          xrefe(0,k) = Nodes()[k]->X()[0];
          xrefe(1,k) = Nodes()[k]->X()[1];
          xcure(0,k) = xrefe(0,k) + mydisp[k*numdf+0];
          xcure(1,k) = xrefe(1,k) + mydisp[k*numdf+1];
        }

        /*------------------------- get node weights for nurbs elements */
        const DiscretizationType distype = Shape();
        Epetra_SerialDenseVector weights(numnode);
        if (distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
        {
          for (int inode=0; inode<numnode; ++inode)
          {
            DRT::NURBS::ControlPoint* cp =
              dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);
            weights(inode) = cp->W();
          }
        }

        //----------------------------------------------------------------
        // loop over all Gauss points
        //----------------------------------------------------------------
        for (int ip=0; ip<intpoints.nquad; ++ip)
        {
          const double e1 = intpoints.qxg[ip][0];
          const double e2 = intpoints.qxg[ip][1];
          const double wgt = intpoints.qwgt[ip];

          // get values of shape functions and derivatives in the gausspoint
          if (distype != DRT::Element::nurbs4 &&
              distype != DRT::Element::nurbs9)
          {
            // shape functions and their derivatives for polynomials
            DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
            DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
          }
          else
          {
            // nurbs version
            Epetra_SerialDenseVector gp(2);
            gp(0)=e1;
            gp(1)=e2;

            DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
            (funct  ,
             deriv  ,
             gp     ,
             myknots,
             weights,
             distype);
          }

          /*--------------------------------------- compute jacobian Matrix */
          w1_jacobianmatrix(xrefe,deriv,xjm,&det,numnode);

          /*------------------------------------ integration factor  -------*/
          double fac = wgt * det * thickness_;

          /*----------------------------------- calculate operator Blin  ---*/
          w1_boplin(boplin,deriv,xjm,det,numnode);

          /*-------------------------deformation gradient and GL strains ---*/
          w1_defgrad(F,strain,xrefe,xcure,boplin,numnode);

          // Gauss point in reference configuration
          LINALG::Matrix<2,1> xgp(true);
          for (int k=0;k<numdf;++k)
            for (int n=0;n<numnode;++n)
              xgp(k,0) += funct[n] * xrefe(k,n);

          //**************************************************************
          // get analytical solution
          LINALG::Matrix<2,1> uanalyt(true);
          LINALG::Matrix<4,1> strainanalyt(true);
          LINALG::Matrix<2,2> derivanalyt(true);

          CONTACT::AnalyticalSolutions2D(xgp,uanalyt,strainanalyt,derivanalyt);
          //**************************************************************

          //--------------------------------------------------------------
          // (1) L2 norm
          //--------------------------------------------------------------

          // compute displacements at GP
          LINALG::Matrix<2,1> ugp(true);
          for (int k=0;k<numdf;++k)
            for (int n=0;n<numnode;++n)
              ugp(k,0) += funct[n] * (xcure(k,n) - xrefe(k,n));

          // displacement error
          LINALG::Matrix<2,1> uerror(true);
          for (int k=0;k<numdf;++k)
            uerror(k,0) = uanalyt(k,0) - ugp(k,0);

          // compute GP contribution to L2 error norm
          l2norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (2) H1 norm
          //--------------------------------------------------------------

          // compute partial derivatives at GP
          LINALG::Matrix<2,2> derivgp(true);
          derivgp(0,0) = F[0] - 1.0;
          derivgp(0,1) = F[2];
          derivgp(1,0) = F[3];
          derivgp(1,1) = F[1] - 1.0;

          // derivative error
          LINALG::Matrix<2,2> deriverror(true);
          for (int k=0;k<numdf;++k)
            for (int m=0;m<numdf;++m)
              deriverror(k,m) = derivanalyt(k,m) - derivgp(k,m);

          // compute GP contribution to H1 error norm
          h1norm += fac * deriverror.Dot(deriverror);
          h1norm += fac * uerror.Dot(uerror);

          //--------------------------------------------------------------
          // (3) Energy norm
          //--------------------------------------------------------------

          // compute linear strain at GP
          LINALG::Matrix<4,1> straingp(true);
          straingp(0,0) = 0.5 * (F[0] + F[0]) - 1.0;
          straingp(1,0) = 0.5 * (F[1] + F[1]) - 1.0;
          straingp(2,0) = 0.5 * (F[2] + F[3]);
          straingp(3,0) = straingp(2,0);

          // strain error
          LINALG::Matrix<4,1> strainerror(true);
          for (int k=0;k<numeps;++k)
            strainerror(k,0) = strainanalyt(k,0) - straingp(k,0);

          // compute stress vector and constitutive matrix
          Epetra_SerialDenseMatrix C;
          C.Shape(4,4);
          Epetra_SerialDenseMatrix tempstress;
          tempstress.Shape(4,4);
          Epetra_SerialDenseVector tempstrainerror(4);
          tempstrainerror[0] = strainerror(0,0);
          tempstrainerror[1] = strainerror(1,0);
          tempstrainerror[2] = strainerror(2,0);
          tempstrainerror[3] = strainerror(3,0);
          Teuchos::RCP<const MAT::Material> material = Material();
          w1_call_matgeononl(tempstrainerror,tempstress,C,numeps,material,params);
          LINALG::Matrix<4,1> stress(true);
          stress(0,0) = tempstress(0,0);
          stress(1,0) = tempstress(1,1);
          stress(2,0) = tempstress(0,2);
          stress(3,0) = tempstress(0,2);


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
        elevec1[0] = l2norm;
        elevec1[1] = h1norm;
        elevec1[2] = energynorm;
      }
      else
        dserror("ERROR: Error norms only implemented for SVK material");
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_mass_volume:
    {
      // check length of elevec1
      if (elevec1.Length() < 6)
        dserror("The given result vector is too short.");

      // declaration of variables
      double volume_ref = 0.0;
      double volume_mat = 0.0;
      double volume_cur = 0.0;
      double mass_ref = 0.0;
      double mass_mat = 0.0;
      double mass_cur = 0.0;
      double density  = actmat->Density();

      // some definitions
      const int numnode = NumNode();
      const int numdf   = 2;

      Epetra_SerialDenseMatrix xjm;
      Epetra_SerialDenseMatrix xjmmat;
      xjm.Shape(2,2);
      xjmmat.Shape(2,2);
      double det     = 0.0;
      double detmat  = 0.0;
      double detcur  = 0.0;
      double detFmat = 0.0;//F[0]*F[1]-F[2]*F[3];

      // shape functions, derivatives and integration rule
      Epetra_SerialDenseVector funct(numnode);
      Epetra_SerialDenseMatrix deriv;
      deriv.Shape(2,numnode);
      const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

      // get displacements and extract values of this element
      Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==Teuchos::null) dserror("Cannot get state displacement vector");
      std::vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

      std::vector<double> mydispmat(lm.size());
      if (structale_)
      {
        Teuchos::RCP<const Epetra_Vector> dispmat = discretization.GetState("material_displacement");
        DRT::UTILS::ExtractMyValues(*dispmat,mydispmat,lm);
      }

      // reference and current geometry (nodal positions)
      Epetra_SerialDenseMatrix xrefe(2,numnode);
      Epetra_SerialDenseMatrix xcure(2,numnode);
      Epetra_SerialDenseMatrix xmat(2,numnode);
      Epetra_SerialDenseVector strain;
      strain.Size(4);
      Epetra_SerialDenseMatrix boplin;
      boplin.Shape(4,2*numnode);
      Epetra_SerialDenseVector F;
      F.Size(4);

      for (int k=0; k<numnode; ++k)
      {
        xrefe(0,k) = Nodes()[k]->X()[0];
        xrefe(1,k) = Nodes()[k]->X()[1];
        xcure(0,k) = xrefe(0,k) + mydisp[k*numdf+0];
        xcure(1,k) = xrefe(1,k) + mydisp[k*numdf+1];

        // material displacements for structure with ale
        if(structale_ == true)
        {
          xmat(0,k)  = xrefe(0,k) + mydispmat[k*numdf+0];
          xmat(1,k)  = xrefe(1,k) + mydispmat[k*numdf+1];
        }
      }

      /*------------------------- get node weights for nurbs elements */
      const DiscretizationType distype = Shape();
      Epetra_SerialDenseVector weights(numnode);
      if (distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
      {
        for (int inode=0; inode<numnode; ++inode)
        {
          DRT::NURBS::ControlPoint* cp =
            dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);
          weights(inode) = cp->W();
        }
      }

      //----------------------------------------------------------------
      // loop over all Gauss points
      //----------------------------------------------------------------
      for (int ip=0; ip<intpoints.nquad; ++ip)
      {
        const double e1  = intpoints.qxg[ip][0];
        const double e2  = intpoints.qxg[ip][1];
        const double wgt = intpoints.qwgt[ip];

        // get values of shape functions and derivatives in the gausspoint
        if (distype != DRT::Element::nurbs4 &&
            distype != DRT::Element::nurbs9)
        {
          // shape functions and their derivatives for polynomials
          DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
          DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
        }
        else
        {
          // nurbs version
          Epetra_SerialDenseVector gp(2);
          gp(0)=e1;
          gp(1)=e2;

          DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
          (funct  ,
           deriv  ,
           gp     ,
           myknots,
           weights,
           distype);
        }

        //REF ------------------------
        /*--------------------------------------- compute jacobian Matrix */
        w1_jacobianmatrix(xrefe,deriv,xjm,&det,numnode);

        /*------------------------------------ integration factor  -------*/
        double fac = wgt * det * thickness_ ;
        volume_ref+=fac;
        fac = wgt * det * thickness_ * density;
        mass_ref+=fac;

        //MAT ------------------------
        if(structale_)
        {
          w1_jacobianmatrix(xmat,deriv,xjmmat,&detmat,numnode);
          fac = wgt * detmat * thickness_ ;
          volume_mat+=fac;
          fac = wgt * detmat * thickness_ * density;
          mass_mat+=fac;

          w1_boplin(boplin,deriv,xjmmat,detmat,numnode);
          w1_defgrad(F,strain,xmat,xcure,boplin,numnode);
          detFmat = F[0]*F[1]-F[2]*F[3];

          //CUR ------------------------
          /*--------------------------------------- compute jacobian Matrix */
          w1_jacobianmatrix(xcure,deriv,xjm,&detcur,numnode);

          /*------------------------------------ integration factor  -------*/
          fac = wgt * detcur * thickness_ ;
          volume_cur+=fac;
          fac = wgt * detcur * thickness_ * density * 1/detFmat;
          mass_cur+=fac;
        }
        else
        {
          w1_boplin(boplin,deriv,xjm,det,numnode);
          w1_defgrad(F,strain,xrefe,xcure,boplin,numnode);
          detFmat = F[0]*F[1]-F[2]*F[3];

          //CUR ------------------------
          /*--------------------------------------- compute jacobian Matrix */
          w1_jacobianmatrix(xcure,deriv,xjm,&detcur,numnode);

          /*------------------------------------ integration factor  -------*/
          fac = wgt * detcur * thickness_ ;
          volume_cur+=fac;
          fac = wgt * detcur * thickness_ * density * 1/detFmat;
          mass_cur+=fac;
        }
      }
      //----------------------------------------------------------------

      // return results
      if(!structale_)
      {
        volume_mat = volume_ref;
        mass_mat   = mass_ref;
      }

      elevec1[0] = volume_ref;
      elevec1[1] = volume_mat;
      elevec1[2] = volume_cur;
      elevec1[3] = mass_ref;
      elevec1[4] = mass_mat;
      elevec1[5] = mass_cur;
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_eleload:
    {
      dserror("this method is not supposed to evaluate a load, use EvaluateNeumann(...)");
      break;
    }
    //==================================================================================
    case ELEMENTS::struct_calc_predict:
      break;
    //==================================================================================
    default:
    {
      dserror("Unknown type of action for Wall1 %d", act);
      break;
    }
  }
  return 0;

}

/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)  mgit 05/07|
 *----------------------------------------------------------------------*/

int DRT::ELEMENTS::Wall1::EvaluateNeumann(Teuchos::ParameterList&   params,
                                          DRT::Discretization&      discretization,
                                          DRT::Condition&           condition,
                                          std::vector<int>&         lm,
                                          Epetra_SerialDenseVector& elevec1,
                                          Epetra_SerialDenseMatrix* elemat1)
{
  SetParamsInterfacePtr(params);
  // get values and switches from the condition
  const std::vector<int>*    onoff = condition.Get<std::vector<int> >   ("onoff");
  const std::vector<double>* val   = condition.Get<std::vector<double> >("val"  );
  const std::vector<int>*    curve = condition.Get<std::vector<int> >   ("curve");
  const std::vector<int>*    funct = condition.Get<std::vector<int> >   ("funct");

  // check total time
  bool usetime = true;
  double time = -1.0;
  if (IsParamsInterface())
    time = ParamsInterface().GetTotalTime();
  else
    time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  // ensure that at least as many curves/functs as dofs are available
  if (int(onoff->size()) < noddof_)
    dserror("Fewer functions or curves defined than the element has dofs.");

  // factor given by time curves
  std::vector<double> curvefacs(noddof_, 1.0);
  for (int i = 0; i < noddof_; ++i)
  {
    const int curvenum = (curve) ? (*curve)[i] : -1;
    if (curvenum >= 0 && usetime)
      curvefacs[i] = DRT::Problem::Instance()->Curve(curvenum).f(time);
  }

  // no. of nodes on this surface
  const int iel = NumNode();

  // do the isogeometric extras --- get knots and weights
  std::vector<Epetra_SerialDenseVector> myknots(numdim_);
  Epetra_SerialDenseVector weights(iel);

  if(Shape()==DRT::Element::nurbs4
     ||
     Shape()==DRT::Element::nurbs9)
  {
    DRT::NURBS::NurbsDiscretization* nurbsdis =
        dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    bool zero_sized=(*((*nurbsdis).GetKnotVector())).GetEleKnots(myknots,Id());

    // skip zero sized elements in knot span --- they correspond to interpolated nodes
    if(zero_sized)
    {
      return(0);
    }

    for (int inode=0; inode<iel; ++inode)
    {
      DRT::NURBS::ControlPoint* cp =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  // general arrays
  Epetra_SerialDenseMatrix xjm(numdim_,numdim_);  // iso-parametric Jacobian
  double det=0.0;  // determinant of iso-parametric Jacobian

  // quad, tri, etc
  const DiscretizationType distype = Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

  //  std::vector<double>* thick = data_.Get<std::vector<double> >("thick");
  //  if (!thick) dserror("Cannot find vector of nodal thickness");

  // shape functions
  Epetra_SerialDenseVector shapefcts(iel);
  // natural derivatives of shape funcions
  Epetra_SerialDenseMatrix deriv(numdim_,iel);

  // reference co-ordinates of element nodes
  Epetra_SerialDenseMatrix xrefe(numdim_,iel);


  /*----------------------------------------------------- geometry update */
  for (int k=0; k<iel; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];
  }

  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    /*================================== gaussian point and weight at it */
    const double e1 = intpoints.qxg[ip][0];
    const double e2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    /*-------------------- shape functions at gp e1,e2 on mid surface */
    if(distype != DRT::Element::nurbs4
       &&
       distype != DRT::Element::nurbs9)
    {
      // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (shapefcts,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
      (shapefcts,
       deriv  ,
       gp     ,
       myknots,
       weights,
       distype);
    }

    /*--------------------------------------- compute jacobian Matrix */
    w1_jacobianmatrix(xrefe,deriv,xjm,&det,iel);

    /*------------------------------------ integration factor  -------*/
    double fac = wgt * det;

    // load vector ar
    double ar[noddof_];
    // loop the dofs of a node
    // ar[i] = ar[i] * facr * ds * onoff[i] * val[i]
    for (int i=0; i<noddof_; ++i)
    {
      // factor given by spatial function
      const int functnum = (funct) ? (*funct)[i] : -1;
      double functfac = 1.0;
      if (functnum > 0)
      {
        // calculate reference position of GP
        LINALG::SerialDenseMatrix gp_coord(1, numdim_);
        gp_coord.Multiply('T', 'T', 1.0, shapefcts, xrefe, 0.0);

        // write coordinates in another datatype
        double gp_coord2[3]; // the position vector has to be given in 3D!!!
        for (int k = 0; k < numdim_; k++)
          gp_coord2[k] = gp_coord(0, k);
        for (int k = numdim_; k < 3; k++) // set a zero value for the remaining spatial directions
          gp_coord2[k] = 0.0;
        const double* coordgpref = &gp_coord2[0]; // needed for function evaluation

        // evaluate function at current gauss point
        functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(i,coordgpref,time,NULL);
      }

      ar[i] = fac * (*onoff)[i] * (*val)[i] * curvefacs[i] * functfac;
    }

    // add load components
    for (int node=0; node<iel; ++node)
      for (int dof=0; dof<noddof_; ++dof)
         elevec1[node*noddof_+dof] += shapefcts[node] * ar[dof];

  } // for (int ip=0; ip<totngp; ++ip)

  // finished
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_recover(
    const std::vector<int>&         lm,
    const std::vector<double>&      disp,
    const std::vector<double>&      residual)
{
  // for eas
  Epetra_SerialDenseMatrix* alpha      = NULL;
  Epetra_SerialDenseMatrix* eas_inc    = NULL;
  // get access to the interface parameters
  const double step_length   = StrParamsInterface().GetStepLength();

  // have eas?
  if (iseas_)
  {
    // access general eas history stuff stored in element
    // get alpha of previous iteration
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");
    // get the old eas increment
    eas_inc = data_.GetMutable<Epetra_SerialDenseMatrix>("eas_inc");
    if (!alpha || !eas_inc)
      dserror("Missing EAS history data (eas_inc and/or alpha)");
  }

  /* if it is a default step, we have to recover the condensed
   * solution vectors */
  if (StrParamsInterface().IsDefaultStep())
  {
    /* recovery of the enhanced assumed strain increment and
     * update of the eas dofs. */
    if (iseas_)
    {
      // first, store the eas state of the previous accepted Newton step
      StrParamsInterface().SumIntoMyPreviousSolNorm(NOX::NLN::StatusTest::quantity_eas,
          w1_neas(),(*alpha)[0],Owner());

      // get stored EAS history
      Epetra_SerialDenseMatrix* oldfeas =
          data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
      Epetra_SerialDenseMatrix* oldKaainv =
          data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
      Epetra_SerialDenseMatrix* oldKda =
          data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
      if (!oldKaainv or !oldKda or !oldfeas)
        dserror("Missing EAS history-data");

      // we need the (residual) displacement at the previous step
      const int numnode = NumNode();
      Epetra_SerialDenseVector res_d(2*numnode);
      for (int i = 0; i < (2*numnode); ++i) {
        res_d(i) = residual[i];
      }

      // add Kda . res_d to feas
      (*oldfeas).Multiply('T','N',1.0,(*oldKda),res_d,1.0);
      // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
      (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);
    } // if (iseas)
  } // if (*isdefault_step_ptr_)
  /* if it is no default step, we can correct the update and the current eas
   * state without the need for any matrix-vector products. */
  else
  {
    // The first step has to be a default step!
    if (old_step_length_<0.0)
      dserror("The old step length was not defined!");
    /* if this is no full step, we have to adjust the length of the
     * enhanced assumed strain incremental step. */
    if (iseas_)
    {
      /* undo the previous step:
       *            alpha_new = alpha_old - old_step * alpha_inc
       * and update the solution variable with the new step length:
       *            alpha_new = alpha_new + new_step * alpha_inc */
      for (int i=0;i<Wall1::neas_;++i)
        (*alpha)(i,0) += (step_length-old_step_length_)*(*eas_inc)(i,0);
    } // if (nhyb_)
  } // else
  // save the old step length
  old_step_length_ = step_length;

  // Check if the eas incr is tested and if yes, calculate the element
  // contribution to the norm
  if (iseas_)
    StrParamsInterface().SumIntoMyUpdateNorm(NOX::NLN::StatusTest::quantity_eas,
        w1_neas(),(*eas_inc)[0],(*alpha)[0],step_length,Owner());

  // the element internal stuff should be up-to-date for now...
  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            mgit 03/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_nlnstiffmass(
  const std::vector<int>                    & lm         ,
  const std::vector<double>                 & disp       ,
  const std::vector<double>                 & residual   ,
  const std::vector<double>                 & dispmat   ,
  std::vector<Epetra_SerialDenseVector>& myknots    ,
  Epetra_SerialDenseMatrix             * stiffmatrix,
  Epetra_SerialDenseMatrix             * massmatrix ,
  Epetra_SerialDenseVector             * force      ,
  Epetra_SerialDenseMatrix             * elestress  ,
  Epetra_SerialDenseMatrix             * elestrain  ,
  Teuchos::RCP<const MAT::Material>      material   ,
  Teuchos::ParameterList& params,  ///< algorithmic parameters e.g. time
  const INPAR::STR::StressType           iostress   ,
  const INPAR::STR::StrainType           iostrain   )
{
  const int numnode = NumNode();
  const int numdf   = 2;
  const int nd      = numnode*numdf;


   // general arrays
  Epetra_SerialDenseVector funct(numnode);
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  Epetra_SerialDenseMatrix boplin;
  boplin.Shape(4,2*numnode);
  Epetra_SerialDenseVector F;
  F.Size(4);
  Epetra_SerialDenseVector strain;
  strain.Size(4);
  double det;
  Epetra_SerialDenseMatrix xrefe(2,numnode);
  Epetra_SerialDenseMatrix xcure(2,numnode);
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure;
  b_cure.Shape(numeps,nd);
  Epetra_SerialDenseMatrix stress;
  stress.Shape(4,4);
  Epetra_SerialDenseMatrix C;
  C.Shape(4,4);

  // for EAS, in any case declare variables, sizes etc. only in eascase
  Epetra_SerialDenseMatrix* alpha=NULL;  // EAS alphas
  Epetra_SerialDenseMatrix F_enh;  // EAS matrix F_enh
  Epetra_SerialDenseMatrix F_tot;  // EAS vector F_tot
  Epetra_SerialDenseMatrix p_stress;  // first piola-kirchhoff stress vector
  Epetra_SerialDenseMatrix xjm0;  // Jacobian Matrix (origin)
  Epetra_SerialDenseVector F0;  // Deformation Gradient (origin)
  Epetra_SerialDenseMatrix boplin0; // B operator (origin)
  Epetra_SerialDenseMatrix W0;  // W operator (origin)
  Epetra_SerialDenseMatrix G;  // G operator
  Epetra_SerialDenseMatrix Z;  // Z operator
  Epetra_SerialDenseMatrix FCF;  // FCF^T
  Epetra_SerialDenseMatrix Kda;  // EAS matrix Kda
  Epetra_SerialDenseMatrix Kaa;  // EAS matrix Kaa
  Epetra_SerialDenseVector feas; // EAS portion of internal forces
  double detJ0;  // detJ(origin)
  Epetra_SerialDenseMatrix* oldfeas  =NULL;   // EAS history
  Epetra_SerialDenseMatrix* oldKaainv=NULL; // EAS history
  Epetra_SerialDenseMatrix* oldKda   =NULL;    // EAS history

  // arrays for structure with ale (fractional step strategy)
  Epetra_SerialDenseMatrix xmat;
  Epetra_SerialDenseMatrix xjmmat;
  Epetra_SerialDenseMatrix boplinmat;
  Epetra_SerialDenseVector Fmat;
  Epetra_SerialDenseVector FFmatinv;
  double detmat;

  if(structale_ == true)
  {
    xmat.Shape(2,numnode);
    xjmmat.Shape(2,2);
    boplinmat.Shape(4,2*numnode);
    Fmat.Size(4);
    FFmatinv.Size(4);
  }

  // ------------------------------------ check calculation of mass matrix
  double density = 0.0;
  if (massmatrix) density = material->Density();

  /*------- get integraton data ---------------------------------------- */
  const DiscretizationType distype = Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<numnode; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];
    xcure(0,k) = xrefe(0,k) + disp[k*numdf+0];
    xcure(1,k) = xrefe(1,k) + disp[k*numdf+1];

    // material displacements for structure with ale
    if(structale_ == true)
    {
      xmat(0,k)  = xrefe(0,k) + dispmat[k*numdf+0];
      xmat(1,k)  = xrefe(1,k) + dispmat[k*numdf+1];
    }
  }

  /*--------------------------------- get node weights for nurbs elements */
  Epetra_SerialDenseVector weights(numnode);
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnode; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  if (iseas_)
  {
    // allocate EAS quantities
    F_enh.Shape(4,1);
    F_tot.Shape(4,3);
    p_stress.Shape(4,1);
    xjm0.Shape(2,2);
    F0.Size(4);
    boplin0.Shape(4,2*numnode);
    W0.Shape(4,2*numnode);
    G.Shape(4,Wall1::neas_);
    Z.Shape(2*numnode,Wall1::neas_);
    FCF.Shape(4,4);
    Kda.Shape(2*numnode,Wall1::neas_);
    Kaa.Shape(Wall1::neas_,Wall1::neas_);
    feas.Size(Wall1::neas_);

    /*
    ** EAS Update of alphas:
    ** the current alphas are (re-)evaluated out of
    ** Kaa and Kda of previous step to avoid additional element call.
    ** This corresponds to the (innermost) element update loop
    ** in the nonlinear FE-Skript page 120 (load-control alg. with EAS)
    */
    alpha = data_.GetMutable<Epetra_SerialDenseMatrix>("alpha");   // get alpha of previous iteration

    // get stored EAS history
    oldfeas = data_.GetMutable<Epetra_SerialDenseMatrix>("feas");
    oldKaainv = data_.GetMutable<Epetra_SerialDenseMatrix>("invKaa");
    oldKda = data_.GetMutable<Epetra_SerialDenseMatrix>("Kda");
    if (!alpha || !oldKaainv || !oldKda || !oldfeas) dserror("Missing EAS history-data");
    // FixMe deprecated implementation
    if (not IsParamsInterface())
    {
      // we need the (residual) displacement at the previous step
      Epetra_SerialDenseVector res_d(2*numnode);
      for (int i = 0; i < (2*numnode); ++i) {
        res_d(i) = residual[i];
      }

      // add Kda . res_d to feas
      (*oldfeas).Multiply('T','N',1.0,(*oldKda),res_d,1.0);
      // new alpha is: - Kaa^-1 . (feas + Kda . old_d), here: - Kaa^-1 . feas
      (*alpha).Multiply('N','N',-1.0,(*oldKaainv),(*oldfeas),1.0);
    }// if (not IsInterface())
    /* end of EAS Update ******************/

    /* evaluation of EAS variables (which are constant for the following):
    ** -> M defining interpolation of enhanced strains alpha, evaluated at GPs
    ** -> determinant of Jacobi matrix at element origin (r=s=t=0.0)
    ** -> T0^{-T}
    */
    w1_eassetup(boplin0,F0,xjm0,detJ0,xrefe,xcure,distype);
  }

  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    /*================================== gaussian point and weight at it */
    const double e1 = intpoints.qxg[ip][0];
    const double e2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs4
       &&
       distype != DRT::Element::nurbs9)
    {
    // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
     (funct  ,
      deriv  ,
      gp     ,
      myknots,
      weights,
      distype);
    }

    /*--------------------------------------- compute jacobian Matrix */
    w1_jacobianmatrix(xrefe,deriv,xjm,&det,numnode);

    /*------------------------------------ integration factor  -------*/
    double fac = wgt * det * thickness_;

    /*------------------------------compute mass matrix if imass-----*/
    if (massmatrix)
    {
      double facm = fac * density;
      for (int a=0; a<numnode; a++)
      {
        for (int b=0; b<numnode; b++)
        {
          (*massmatrix)(2*a,2*b)     += facm * funct(a) * funct(b); /* a,b even */
          (*massmatrix)(2*a+1,2*b+1) += facm * funct(a) * funct(b); /* a,b odd  */
        }
      }
    }

    /*----------------------------------- calculate operator Blin  ---*/
    w1_boplin(boplin,deriv,xjm,det,numnode);
    //cout.precision(16);
    /*------------ calculate defgrad F^u, Green-Lagrange-strain E^u --*/
    w1_defgrad(F,strain,xrefe,xcure,boplin,numnode);

    // modifications for structural approch with ale
    if(structale_ == true)
    {
      /* -calculate defgrad F^mat, correct Green-Lagrange-strain E^u -*/
      w1_defgradmat(F,Fmat,FFmatinv,strain,xrefe,xmat,boplin,numnode);

      /*---------- compute jacobian Matrix (material configuration) --*/
      w1_jacobianmatrix(xmat,deriv,xjmmat,&detmat,numnode);

      /*---------- calculate operator Blin (material configuration) --*/
      w1_boplin(boplinmat,deriv,xjmmat,detmat,numnode);

      /* -----------------------------replace factors and operators --*/
      fac = wgt * detmat * thickness_;
      boplin = boplinmat;
      F = FFmatinv;
    }

    /*-calculate defgrad F in matrix notation and Blin in current conf.*/
    w1_boplin_cure(b_cure,boplin,F,numeps,nd);

    // EAS technology: "enhance the deformation gradient"  ---- --- EAS
    if (iseas_ == true)
    {
      /*-----calculate the enhanced deformation gradient and--------------------
      -----alsoe the operators G, W0 and Z------------------------------------*/

      w1_call_defgrad_enh(F_enh,xjm0,xjm,detJ0,det,F0,*alpha,e1,e2,G,W0,boplin0,Z);

      /*-----total deformation gradient, Green-Lagrange-strain E^F -----------*/
      w1_call_defgrad_tot(F_enh,F_tot,F,strain);
      /* call material law----------------------------------------------------*/
      w1_call_matgeononl(strain,stress,C,numeps,material,params);

      // return gp strains (only in case of strain output)
      switch (iostrain)
      {
      case INPAR::STR::strain_gl:
      {
        if (elestrain == NULL) dserror("no strain data available");
        (*elestrain)(ip,0) = strain(0);
        (*elestrain)(ip,1) = strain(1);
        (*elestrain)(ip,2) = 0.0;
        (*elestrain)(ip,3) = strain(3);
      }
      break;
      case INPAR::STR::strain_none:
        break;
      case INPAR::STR::strain_ea:
      default:
        dserror("requested strain type not supported");
        break;
      }

      // return gp stresses (only in case of stress output)
      switch (iostress)
      {
      case INPAR::STR::stress_2pk:
      {
        if (elestress == NULL) dserror("no stress data available");
        (*elestress)(ip,0) = stress(0,0);
        (*elestress)(ip,1) = stress(1,1);
        (*elestress)(ip,2) = 0.0;
        (*elestress)(ip,3) = stress(0,2);
      }
      break;
      case INPAR::STR::stress_cauchy:
      {
        if (elestress == NULL) dserror("no stress data available");
        StressCauchy(ip, F_tot(0,0), F_tot(1,1), F_tot(0,2), F_tot(1,2), stress, elestress);
      }
      break;
      case INPAR::STR::stress_none:
        break;
      default:
        dserror("requested stress type not supported");
        break;
      }

      /*-----first piola-kirchhoff stress vector------------------------------*/
      w1_stress_eas(stress,F_tot,p_stress);

      /*-----stiffness matrix kdd---------------------------------------------*/
      if (stiffmatrix) w1_kdd(boplin,W0,F_tot,C,stress,FCF,*stiffmatrix,fac);
      /*-----matrix kda-------------------------------------------------------*/
      w1_kda(FCF,W0,boplin,stress,G,Z,Kda,p_stress,fac);
      /*-----matrix kaa-------------------------------------------------------*/
      w1_kaa(FCF,stress,G,Kaa,fac);
      /*-----nodal forces ----------------------------------------------------*/
      if (force) w1_fint_eas(W0,boplin,G,p_stress,*force,feas,fac);

   }
   else
   {
     params.set<int>("gp",ip);
     w1_call_matgeononl(strain,stress,C,numeps,material,params);

     // return gp strains (only in case of strain output)
     switch (iostrain)
     {
     case INPAR::STR::strain_gl:
     {
       if (elestrain == NULL) dserror("no strain data available");
       (*elestrain)(ip,0) = strain(0);
       (*elestrain)(ip,1) = strain(1);
       (*elestrain)(ip,2) = 0.0;
       (*elestrain)(ip,3) = strain(3);;
     }
     break;
     case INPAR::STR::strain_none:
        break;
     case INPAR::STR::strain_ea:
     default:
        dserror("requested strain type not supported");
        break;
     }

     // return gp stresses (only in case of stress output)
     switch (iostress)
     {
     case INPAR::STR::stress_2pk:
     {
       if (elestress == NULL) dserror("no stress data available");
       (*elestress)(ip,0) = stress(0,0);
       (*elestress)(ip,1) = stress(1,1);
       (*elestress)(ip,2) = 0.0;
       (*elestress)(ip,3) = stress(0,2);
     }
     break;
     case INPAR::STR::stress_cauchy:
     {
       if (elestress == NULL) dserror("no stress data available");
       StressCauchy(ip, F[0], F[1], F[2], F[3], stress, elestress);
     }
     break;
     case INPAR::STR::stress_none:
       break;
     default:
       dserror("requested stress type not supported");
       break;
     }

     /*---------------------- geometric part of stiffness matrix kg ---*/
     if (stiffmatrix) w1_kg(*stiffmatrix,boplin,stress,fac,nd,numeps);

     /*------------------ elastic+displacement stiffness matrix keu ---*/
     if (stiffmatrix) w1_keu(*stiffmatrix,b_cure,C,fac,nd,numeps);

     /*--------------- nodal forces fi from integration of stresses ---*/
     if (force) w1_fint(stress,b_cure,*force,fac,nd);
   }

  } // for (int ip=0; ip<totngp; ++ip)


  // EAS technology: ------------------------------------------------------ EAS
  // subtract EAS matrices from disp-based Kdd to "soften" element

  if (force != NULL && stiffmatrix != NULL)
  {
    if (iseas_ == true)
    {
      // we need the inverse of Kaa
      Epetra_SerialDenseSolver solve_for_inverseKaa;
      solve_for_inverseKaa.SetMatrix(Kaa);
      solve_for_inverseKaa.Invert();


      Epetra_SerialDenseMatrix KdaKaa(2*NumNode(),Wall1::neas_); // temporary Kda.Kaa^{-1}
      KdaKaa.Multiply('N', 'N', 1.0, Kda, Kaa, 1.0);


      // EAS-stiffness matrix is: Kdd - Kda^T . Kaa^-1 . Kad  with Kad=Kda^T
      if (stiffmatrix) (*stiffmatrix).Multiply('N', 'T', -1.0, KdaKaa, Kda, 1.0);

      // EAS-internal force is: fint - Kda^T . Kaa^-1 . feas
      if (force) (*force).Multiply('N', 'N', -1.0, KdaKaa, feas, 1.0);

      // store current EAS data in history
      for (int i=0; i<Wall1::neas_; ++i)
        for (int j=0; j<Wall1::neas_; ++j)
          (*oldKaainv)(i,j) = Kaa(i,j);

      for (int i=0; i<(2*NumNode()); ++i)
        for (int j=0; j<Wall1::neas_; ++j)
        {
          (*oldKda)(i,j) = Kda(i,j);
          (*oldfeas)(j,0) = feas(j);
        }

    }
  }
  // -------------------------------------------------------------------- EAS

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (private)                            popp 09/11|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_linstiffmass(
  const std::vector<int>                    & lm         ,
  const std::vector<double>                 & disp       ,
  const std::vector<double>                 & residual   ,
  const std::vector<double>                 & dispmat   ,
  std::vector<Epetra_SerialDenseVector>& myknots    ,
  Epetra_SerialDenseMatrix             * stiffmatrix,
  Epetra_SerialDenseMatrix             * massmatrix ,
  Epetra_SerialDenseVector             * force      ,
  Epetra_SerialDenseMatrix             * elestress  ,
  Epetra_SerialDenseMatrix             * elestrain  ,
  Teuchos::RCP<const MAT::Material>      material   ,
  Teuchos::ParameterList&                params     ,
  const INPAR::STR::StressType           iostress   ,
  const INPAR::STR::StrainType           iostrain   )
{
  const int numnode = NumNode();
  const int numdf   = 2;
  const int nd      = numnode*numdf;

   // general arrays
  Epetra_SerialDenseVector funct(numnode);
  Epetra_SerialDenseMatrix deriv;
  deriv.Shape(2,numnode);
  Epetra_SerialDenseMatrix xjm;
  xjm.Shape(2,2);
  Epetra_SerialDenseMatrix boplin;
  boplin.Shape(4,2*numnode);
  Epetra_SerialDenseVector F;
  F.Size(4);
  Epetra_SerialDenseVector strain;
  strain.Size(4);
  double det;
  Epetra_SerialDenseMatrix xrefe(2,numnode);
  Epetra_SerialDenseMatrix xcure(2,numnode);
  const int numeps = 4;
  Epetra_SerialDenseMatrix b_cure;
  b_cure.Shape(numeps,nd);
  Epetra_SerialDenseMatrix stress;
  stress.Shape(4,4);
  Epetra_SerialDenseMatrix C;
  C.Shape(4,4);

  // ------------------------------------ check calculation of mass matrix
  double density = 0.0;
  if (massmatrix) density = material->Density();

  /*------- get integraton data ---------------------------------------- */
  const DiscretizationType distype = Shape();

  // gaussian points
  const DRT::UTILS::IntegrationPoints2D  intpoints(gaussrule_);

  /*----------------------------------------------------- geometry update */
  for (int k=0; k<numnode; ++k)
  {
    xrefe(0,k) = Nodes()[k]->X()[0];
    xrefe(1,k) = Nodes()[k]->X()[1];
    xcure(0,k) = xrefe(0,k) + disp[k*numdf+0];
    xcure(1,k) = xrefe(1,k) + disp[k*numdf+1];
  }

  /*--------------------------------- get node weights for nurbs elements */
  Epetra_SerialDenseVector weights(numnode);
  if(distype==DRT::Element::nurbs4 || distype==DRT::Element::nurbs9)
  {
    for (int inode=0; inode<numnode; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
        =
        dynamic_cast<DRT::NURBS::ControlPoint* > (Nodes()[inode]);

      weights(inode) = cp->W();
    }
  }

  /*=================================================== integration loops */
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    /*================================== gaussian point and weight at it */
    const double e1 = intpoints.qxg[ip][0];
    const double e2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // get values of shape functions and derivatives in the gausspoint
    if(distype != DRT::Element::nurbs4
       &&
       distype != DRT::Element::nurbs9)
    {
    // shape functions and their derivatives for polynomials
      DRT::UTILS::shape_function_2D       (funct,e1,e2,distype);
      DRT::UTILS::shape_function_2D_deriv1(deriv,e1,e2,distype);
    }
    else
    {
      // nurbs version
      Epetra_SerialDenseVector gp(2);
      gp(0)=e1;
      gp(1)=e2;

      DRT::NURBS::UTILS::nurbs_get_2D_funct_deriv
  (funct  ,
   deriv  ,
   gp     ,
   myknots,
   weights,
   distype);
    }

    /*--------------------------------------- compute jacobian Matrix */
    w1_jacobianmatrix(xrefe,deriv,xjm,&det,numnode);

    /*------------------------------------ integration factor  -------*/
    double fac = wgt * det * thickness_;

    /*------------------------------compute mass matrix if imass-----*/
    if (massmatrix)
    {
      double facm = fac * density;
      for (int a=0; a<numnode; a++)
      {
        for (int b=0; b<numnode; b++)
        {
          (*massmatrix)(2*a,2*b)     += facm * funct(a) * funct(b); /* a,b even */
          (*massmatrix)(2*a+1,2*b+1) += facm * funct(a) * funct(b); /* a,b odd  */
        }
      }
    }

    /*----------------------------------- calculate operator Blin  ---*/
    w1_boplin(boplin,deriv,xjm,det,numnode);

    /*-------------------------deformation gradient and GL strains ---*/
    w1_defgrad(F,strain,xrefe,xcure,boplin,numnode);

    /*--------------redefine strains -> linear engineering strains ---*/
    strain[0] = 0.5 * (F[0] + F[0]) - 1.0;
    strain[1] = 0.5 * (F[1] + F[1]) - 1.0;
    strain[2] = 0.5 * (F[2] + F[3]);
    strain[3] = strain[2];

    // material call
    w1_call_matgeononl(strain,stress,C,numeps,material,params);

    // return gp strains (only in case of strain output)
    switch (iostrain)
    {
    case INPAR::STR::strain_gl:
    {
     if (elestrain == NULL) dserror("no strain data available");
     (*elestrain)(ip,0) = strain(0);
     (*elestrain)(ip,1) = strain(1);
     (*elestrain)(ip,2) = 0.0;
     (*elestrain)(ip,3) = strain(3);;
    }
    break;
    case INPAR::STR::strain_none:
      break;
    case INPAR::STR::strain_ea:
    default:
      dserror("requested strain type not supported");
      break;
   }

   // return gp stresses (only in case of stress output)
   switch (iostress)
   {
   case INPAR::STR::stress_2pk:
   {
     if (elestress == NULL) dserror("no stress data available");
     (*elestress)(ip,0) = stress(0,0);
     (*elestress)(ip,1) = stress(1,1);
     (*elestress)(ip,2) = 0.0;
     (*elestress)(ip,3) = stress(0,2);
   }
   break;
   case INPAR::STR::stress_cauchy:
   {
     if (elestress == NULL) dserror("no stress data available");
     StressCauchy(ip, F[0], F[1], F[2], F[3], stress, elestress);
   }
   break;
   case INPAR::STR::stress_none:
     break;
   default:
     dserror("requested stress type not supported");
     break;
   }

   /*-------------------------------- linear stiffness matrix keu ---*/
   if (stiffmatrix) w1_keu(*stiffmatrix,boplin,C,fac,nd,numeps);

   /*--------------- nodal forces fi from integration of stresses ---*/
   if (force) w1_fint(stress,boplin,*force,fac,nd);

  } // for (int ip=0; ip<totngp; ++ip)

  return;
}

/*----------------------------------------------------------------------*
 |  jacobian matrix (private)                                  mgit 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_jacobianmatrix(
  const Epetra_SerialDenseMatrix& xrefe,
  const Epetra_SerialDenseMatrix& deriv,
  Epetra_SerialDenseMatrix& xjm,
  double* det,
  const int iel
)
{
  memset(xjm.A(),0,xjm.N()*xjm.M()*sizeof(double));

  for (int k=0; k<iel; k++)
  {
    xjm(0,0) += deriv(0,k) * xrefe(0,k);
    xjm(0,1) += deriv(0,k) * xrefe(1,k);
    xjm(1,0) += deriv(1,k) * xrefe(0,k);
    xjm(1,1) += deriv(1,k) * xrefe(1,k);
  }

/*------------------------------------------ determinant of jacobian ---*/
  *det = xjm[0][0]* xjm[1][1] - xjm[1][0]* xjm[0][1];

  if (*det<0.0) dserror("NEGATIVE JACOBIAN DETERMINANT %8.5f in ELEMENT %d\n",*det,Id());
/*----------------------------------------------------------------------*/

  return;
} // DRT::ELEMENTS::Wall1::w1_jacobianmatrix

/*----------------------------------------------------------------------*
 |  Matrix boplin in reference configuration (private)         mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_boplin(Epetra_SerialDenseMatrix& boplin,
                                     Epetra_SerialDenseMatrix& deriv,
                                     Epetra_SerialDenseMatrix& xjm,
                                     double& det,
                                     const int iel)
{

  double dum;
  double xji[2][2];
  /*---------------------------------------------- inverse of jacobian ---*/
  dum = 1.0/det;
  xji[0][0] = xjm(1,1)* dum;
  xji[0][1] =-xjm(0,1)* dum;
  xji[1][0] =-xjm(1,0)* dum;
  xji[1][1] = xjm(0,0)* dum;
  /*----------------------------- get operator boplin of global derivatives -*/
  /*-------------- some comments, so that even fluid people are able to
   understand this quickly :-)
   the Boplin looks like
       | Nk,x    0   |
       |   0    Nk,y |
       | Nk,y    0   |
       |  0     Nk,x |
  */
  for (int inode=0; inode<iel; inode++)
  {
    int dnode = inode*2;

    boplin(0,dnode+0) = deriv(0,inode)*xji[0][0] + deriv(1,inode)*xji[0][1];
    boplin(1,dnode+1) = deriv(0,inode)*xji[1][0] + deriv(1,inode)*xji[1][1];
    boplin(2,dnode+0) = boplin(1,dnode+1);
    boplin(3,dnode+1) = boplin(0,dnode+0);
  } /* end of loop over nodes */
  return;
}

/* DRT::ELEMENTS::Wall1::w1_boplin */

/*----------------------------------------------------------------------*
 | Deformation gradient F and Green-Langrange strain (private)  mgit 04/07|
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_defgrad(Epetra_SerialDenseVector& F,
                           Epetra_SerialDenseVector& strain,
                           const Epetra_SerialDenseMatrix& xrefe,
                           const Epetra_SerialDenseMatrix& xcure,
                           Epetra_SerialDenseMatrix& boplin,
                           const int iel)
{
  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,X  |
        |  1 + Uy,Y  |
        |      Ux,Y  |
        |      Uy,X  |
  */

  memset(F.A(),0,F.N()*F.M()*sizeof(double));

  F[0] = 1;
  F[1] = 1;
  for (int inode=0; inode<iel; inode++)
  {
     F[0] += boplin(0,2*inode)   * (xcure(0,inode) - xrefe(0,inode));  // F_11
     F[1] += boplin(1,2*inode+1) * (xcure(1,inode) - xrefe(1,inode));  // F_22
     F[2] += boplin(2,2*inode)   * (xcure(0,inode) - xrefe(0,inode));  // F_12
     F[3] += boplin(3,2*inode+1) * (xcure(1,inode) - xrefe(1,inode));  // F_21
  } /* end of loop over nodes */

  /*-----------------------calculate Green-Lagrange strain E -------------*/
  strain[0] = 0.5 * (F[0] * F[0] + F[3] * F[3] - 1.0);  // E_11
  strain[1] = 0.5 * (F[2] * F[2] + F[1] * F[1] - 1.0);  // E_22
  strain[2] = 0.5 * (F[0] * F[2] + F[3] * F[1]);        // E_12
  strain[3] = strain[2];                                // E_21

  /*-----------------------linear engineering strain eps -----------------*/
  /* (choose 2PK stresses for stress output, when using linear strains!)  */
  //strain[0] = 0.5 * (F[0] + F[0]) - 1.0;
  //strain[1] = 0.5 * (F[1] + F[1]) - 1.0;
  //strain[2] = 0.5 * (F[2] + F[3]);
  //strain[3] = strain[2];

  return;
}

/* DRT::ELEMENTS::Wall1::w1_defgrad */

/*----------------------------------------------------------------------*
 | Deformation gradient Fmat and Green-Langrange strain       mgit 04/11|
 | due to structure with ale approach (fractional step method)
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_defgradmat(Epetra_SerialDenseVector& F,
                           Epetra_SerialDenseVector& Fmat,
                           Epetra_SerialDenseVector& FFmatinv,
                           Epetra_SerialDenseVector& strain,
                           const Epetra_SerialDenseMatrix& xrefe,
                           const Epetra_SerialDenseMatrix& xmat,
                           Epetra_SerialDenseMatrix& boplin,
                           const int iel)
{
  /*------------------calculate defgrad --------- (Summenschleife->+=) ---*
  defgrad looks like:

        |  1 + Ux,X  |
        |  1 + Uy,Y  |
        |      Ux,Y  |
        |      Uy,X  |
  */

  memset(Fmat.A(),0,Fmat.N()*Fmat.M()*sizeof(double));

  Fmat[0] = 1;
  Fmat[1] = 1;

  for (int inode=0; inode<iel; inode++)
  {
     Fmat[0] += boplin(0,2*inode)   * (xmat(0,inode) - xrefe(0,inode));  // F_11
     Fmat[1] += boplin(1,2*inode+1) * (xmat(1,inode) - xrefe(1,inode));  // F_22
     Fmat[2] += boplin(2,2*inode)   * (xmat(0,inode) - xrefe(0,inode));  // F_12
     Fmat[3] += boplin(3,2*inode+1) * (xmat(1,inode) - xrefe(1,inode));  // F_21
  } /* end of loop over nodes */

  // determinant of deformation gradient Fmat
  double detFmat = Fmat[0]*Fmat[1]-Fmat[2]*Fmat[3];

  Epetra_SerialDenseVector Fmatinv;
  Fmatinv.Size(4);

  // inverse of Fmat
  Fmatinv[0]=1/detFmat*Fmat[1];
  Fmatinv[1]=1/detFmat*Fmat[0];
  Fmatinv[2]=-1/detFmat*Fmat[2];
  Fmatinv[3]=-1/detFmat*Fmat[3];

  // F.Fmatinv
  FFmatinv[0]=F[0]*Fmatinv[0]+F[2]*Fmatinv[3];
  FFmatinv[1]=F[3]*Fmatinv[2]+F[1]*Fmatinv[1];
  FFmatinv[2]=F[0]*Fmatinv[2]+F[2]*Fmatinv[1];
  FFmatinv[3]=F[3]*Fmatinv[0]+F[1]*Fmatinv[3];

  /*-----------------------calculate Green-Lagrange strain E -------------*/

  strain[0] = 0.5 * (FFmatinv[0] * FFmatinv[0] + FFmatinv[3] * FFmatinv[3] - 1.0);  // E_11
  strain[1] = 0.5 * (FFmatinv[2] * FFmatinv[2] + FFmatinv[1] * FFmatinv[1] - 1.0);  // E_22
  strain[2] = 0.5 * (FFmatinv[0] * FFmatinv[2] + FFmatinv[3] * FFmatinv[1]);        // E_12
  strain[3] = strain[2];                                                            // E_21

  return;
}
/* DRT::ELEMENTS::Wall1::w1_defgradmat */

/*----------------------------------------------------------------------*
 | Deformation gradient F in matrix notation and B in
 reference configuration (private)                             mgit 04/07|
 *----------------------------------------------------------------------*/

void DRT::ELEMENTS::Wall1::w1_boplin_cure(Epetra_SerialDenseMatrix& b_cure,
                                          const Epetra_SerialDenseMatrix& boplin,
                                          const Epetra_SerialDenseVector& F,
                                          const int numeps,
                                          const int nd)
{


     Epetra_SerialDenseMatrix Fmatrix;
     Fmatrix.Shape(4,4);


  /*---------------------------write Vector F as a matrix Fmatrix*/

     Fmatrix(0,0) = F[0];
     Fmatrix(0,2) = 0.5 * F[2];
     Fmatrix(0,3) = 0.5 * F[2];
     Fmatrix(1,1) = F[1];
     Fmatrix(1,2) = 0.5 * F[3];
     Fmatrix(1,3) = 0.5 * F[3];
     Fmatrix(2,1) = F[2];
     Fmatrix(2,2) = 0.5 * F[0];
     Fmatrix(2,3) = 0.5 * F[0];
     Fmatrix(3,0) = F[3];
     Fmatrix(3,2) = 0.5 * F[1];
     Fmatrix(3,3) = 0.5 * F[1];

    /*-------------------------------------------------int_b_cure operator*/
      memset(b_cure.A(),0,b_cure.N()*b_cure.M()*sizeof(double));
      for(int i=0; i<numeps; i++)
        for(int j=0; j<nd; j++)
          for(int k=0; k<numeps; k++)
            b_cure(i,j) += Fmatrix(k,i)*boplin(k,j);
    /*----------------------------------------------------------------*/

  return;
}

/* DRT::ELEMENTS::Wall1::w1_boplin_cure */


//{
//  Teuchos::RCP<MAT::Material> mat = Material();
//  Epetra_SerialDenseMatrix cmat;
//
//  switch(material->mattyp)
//  {
//    case m_stvenant: /*------------------ st.venant-kirchhoff-material */
//    {
//      MAT::StVenantKirchhoff* stvk = static_cast <MAT::StVenantKirchhoff*>(mat.get());
//
//      stvk->Evaluate(glstrain,cmat,stress);
//
//      *density = stvk->Density();
//
//      break;
//    }
//    default:
//      dserror("Illegal type %d of material for wall1 element ", mat->MaterialType());
//      break;
//  }
//
//  /*--------------------------------------------------------------------*/
//  return;
//}  // of w1_mat_sel

/*----------------------------------------------------------------------*
| geometric stiffness part (total lagrange)                   mgit 05/07|
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_kg(Epetra_SerialDenseMatrix& estif,
                                 const Epetra_SerialDenseMatrix& boplin,
                                 const Epetra_SerialDenseMatrix& stress,
                                 const double fac,
                                 const int nd,
                                 const int numeps)
{
  /*---------------------------------------------- perform B^T * SIGMA * B*/
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
      for(int r=0; r<numeps; r++)
         for(int m=0; m<numeps; m++)
            estif(i,j) += boplin(r,i)*stress(r,m)*boplin(m,j)*fac;

  return;

}  // DRT::ELEMENTS::Wall1::w1_kg

/*----------------------------------------------------------------------*
| elastic and initial displacement stiffness (total lagrange)  mgit 05/07
*----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_keu(Epetra_SerialDenseMatrix& estif,
                                  const Epetra_SerialDenseMatrix& b_cure,
                                  const Epetra_SerialDenseMatrix& C,
                                  const double fac,
                                  const int nd,
                                  const int numeps)
{

  /*------------- perform B_cure^T * D * B_cure, whereas B_cure = F^T * B */
  for(int i=0; i<nd; i++)
     for(int j=0; j<nd; j++)
        for(int k=0; k<numeps; k++)
           for(int m=0; m<numeps; m++)
             estif(i,j) +=  b_cure(k,i)*C(k,m)*b_cure(m,j)*fac;

  return;
}  // DRT::ELEMENTS::Wall1::w1_keu


/*----------------------------------------------------------------------*
 | evaluate internal element forces for large def (total Lagr) mgit 05/07  |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_fint(const Epetra_SerialDenseMatrix& stress,
                                   const Epetra_SerialDenseMatrix& b_cure,
                                   Epetra_SerialDenseVector& intforce,
                                   const double fac,
                                   const int nd)

{
  Epetra_SerialDenseVector st;
  st.Size(4);

  st[0] = fac * stress(0,0);
  st[1] = fac * stress(1,1);
  st[2] = fac * stress(0,2);
  st[3] = fac * stress(0,2);

  for(int i=0; i<nd; i++)
    for(int j=0; j<4; j++)
      intforce[i] += b_cure(j,i)*st[j];

  return;
}  // DRT::ELEMENTS::Wall1::w1_fint


/*-----------------------------------------------------------------------------*
| lump mass matrix                                                  bborn 07/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::w1_lumpmass(Epetra_SerialDenseMatrix* emass)
{
  // lump mass matrix
  if (emass != NULL)
  {
    // we assume #elemat2 is a square matrix
    for (int c=0; c<(*emass).N(); ++c)  // parse columns
    {
      double d = 0.0;
      for (int r=0; r<(*emass).M(); ++r)  // parse rows
      {
        d += (*emass)(r,c);  // accumulate row entries
        (*emass)(r,c) = 0.0;
      }
      (*emass)(c,c) = d;  // apply sum of row entries on diagonal
    }
  }
}  // w1_lumpmass

/*-----------------------------------------------------------------------------*
| deliver Cauchy stress                                             bborn 08/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::StressCauchy(
  const int ip,
  const double& F11,
  const double& F22,
  const double& F12,
  const double& F21,
  const Epetra_SerialDenseMatrix& stress,
  Epetra_SerialDenseMatrix* elestress
)
{
  // Question: Is this true for plane stress and/or plane strain mode?

  double detf = F11*F22 - F12*F21;
  // Def.grad. tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix defgrad(2,2);
  defgrad(0,0) = F11;
  defgrad(0,1) = F12;
  defgrad(1,0) = F21;
  defgrad(1,1) = F22;
  // PK2 stress tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix pk2stress(2,2);
  pk2stress(0,0) = stress(0,0);
  pk2stress(0,1) = stress(0,2);
  pk2stress(1,0) = stress(0,2);
  pk2stress(1,1) = stress(1,1);

  // PK1 stress tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix pk1stress(2,2);
  pk1stress.Multiply('N','T',1.0/detf,pk2stress,defgrad,0.0);

  // Cauchy stress tensor in Cartesian matrix notation
  Epetra_SerialDenseMatrix cauchystress(2,2);
  cauchystress.Multiply('N','N',1.0,defgrad,pk1stress,0.0);

  // copy results to array for output
  (*elestress)(ip,0) = cauchystress(0,0);
  (*elestress)(ip,1) = cauchystress(1,1);
  (*elestress)(ip,2) = 0.0;
  (*elestress)(ip,3) = cauchystress(0,1);
}  // StressCauchy


/*-----------------------------------------------------------------------------*
| deliver Cauchy stress                                             bborn 08/08|
*-----------------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1::Energy(
  Teuchos::ParameterList& params,
  const std::vector<int>& lm,
  const std::vector<double>& dis,
  Epetra_SerialDenseVector* energies,
  Teuchos::RCP<const MAT::Material> material
)
{
  // constants
  // element porperties
  const int numnode = NumNode();
  const int edof = numnode * Wall1::noddof_;
  const DiscretizationType distype = Shape();
  // Gaussian points
  const DRT::UTILS::IntegrationPoints2D intpoints(gaussrule_);

  // general arrays
  Epetra_SerialDenseVector shpfct(numnode);  // shape functions at Gauss point
  Epetra_SerialDenseMatrix shpdrv(Wall1::numdim_,numnode);  // parametric derivatives of shape funct. at Gauss point
  Epetra_SerialDenseMatrix Xjm(Wall1::numdim_,Wall1::numdim_);  // material-to-parameter-space Jacobian
  double Xjdet;  // determinant of #Xjm
  Epetra_SerialDenseMatrix boplin(4,edof);
  Epetra_SerialDenseVector Fuv(4);  // disp-based def.grad. vector at t_{n}
  Epetra_SerialDenseVector Ev(4);  // Green-Lagrange strain vector at t_{n}
  Epetra_SerialDenseMatrix Xe(Wall1::numdim_,numnode);  // material/initial element co-ordinates
  Epetra_SerialDenseMatrix xe(Wall1::numdim_,numnode);  // spatial/current element co-ordinates at t_{n}
  Epetra_SerialDenseMatrix bop(Wall1::numstr_,edof);  // non-linear B-op at t_{n}

  Epetra_SerialDenseMatrix massmatrix(lm.size(),lm.size());

  // for EAS, in any case declare variables, sizes etc. only allocated in EAS version
  Epetra_SerialDenseMatrix* alphao=NULL;  // EAS alphas at t_{n}
  Epetra_SerialDenseMatrix Fenhv;  // EAS matrix Fenhv
  Epetra_SerialDenseMatrix Fm;  // total def.grad. matrix at t_{n}
  Epetra_SerialDenseMatrix Xjm0;  // Jacobian Matrix (origin)
  double Xjdet0;  // determinant of #Xjm0
  Epetra_SerialDenseVector Fuv0;  // deformation gradient at origin at t_{n}
  Epetra_SerialDenseMatrix boplin0; // B-operator (origin)
  Epetra_SerialDenseMatrix W0;  // W-operator (origin) at t_{n}
  Epetra_SerialDenseMatrix G;  // G-operator at t_{n}
  Epetra_SerialDenseMatrix Z;  // Z-operator

  // element co-ordinates
  for (int k=0; k<numnode; ++k)
  {
    Xe(0,k) = Nodes()[k]->X()[0];
    Xe(1,k) = Nodes()[k]->X()[1];
    xe(0,k) = Xe(0,k) + dis[k*Wall1::noddof_+0];
    xe(1,k) = Xe(1,k) + dis[k*Wall1::noddof_+1];
  }

  // set-up EAS parameters
  if (iseas_)
  {
    // allocate EAS quantities
    Fenhv.Shape(4,1);
    Fm.Shape(4,3);
    Xjm0.Shape(2,2);
    Fuv0.Size(4);
    boplin0.Shape(4,edof);
    W0.Shape(4,edof);
    G.Shape(4,Wall1::neas_);
    Z.Shape(edof,Wall1::neas_);

    // get alpha of last converged state
    alphao = data_.GetMutable<Epetra_SerialDenseMatrix>("alphao");

    // derivatives at origin
    DRT::UTILS::shape_function_2D_deriv1(shpdrv, 0.0, 0.0, distype);
    // material-to-parameter space Jacobian at origin
    w1_jacobianmatrix(Xe, shpdrv, Xjm0, &Xjdet0, numnode);
    // calculate linear B-operator at origin
    w1_boplin(boplin0, shpdrv, Xjm0, Xjdet0, numnode);
    // displ.-based def.grad. at origin
    w1_defgrad(Fuv0, Ev, Xe, xe, boplin0, numnode);  // at t_{n}
  }

  // integration loops over element domain
  for (int ip=0; ip<intpoints.nquad; ++ip)
  {
    // Gaussian point and weight at it
    const double xi1 = intpoints.qxg[ip][0];
    const double xi2 = intpoints.qxg[ip][1];
    const double wgt = intpoints.qwgt[ip];

    // shape functions and their derivatives
    DRT::UTILS::shape_function_2D(shpfct, xi1, xi2, distype);
    DRT::UTILS::shape_function_2D_deriv1(shpdrv, xi1, xi2, distype);

    // compute Jacobian matrix
    w1_jacobianmatrix(Xe, shpdrv, Xjm, &Xjdet, numnode);

    // integration factor
    double fac = wgt * Xjdet * thickness_;

    // calculate linear B-operator
    w1_boplin(boplin, shpdrv, Xjm, Xjdet, numnode);

    // calculate defgrad F^u, Green-Lagrange-strain E^u
    w1_defgrad(Fuv, Ev, Xe, xe, boplin, numnode);  // at t_{n}

    // calculate non-linear B-operator in current configuration
    w1_boplin_cure(bop, boplin, Fuv, Wall1::numstr_, edof);  // at t_{n} // CHECK THIS: NOT SURE IF bopo NEEDED

    // EAS: The deformation gradient is enhanced
    if (iseas_)
    {
      // calculate the enhanced deformation gradient and
      // also the operators G, W0 and Z
      w1_call_defgrad_enh(Fenhv, Xjm0, Xjm, Xjdet0, Xjdet, Fuv0, *alphao, xi1, xi2, G, W0, boplin0, Z); // at t_{n}

      // total deformation gradient F, and total Green-Lagrange-strain E
      w1_call_defgrad_tot(Fenhv, Fm, Fuv, Ev);  // at t_{n}
    }

    // internal/strain energy
    if (energies) (*energies)(0) += fac * EnergyInternal(material, params, Ev);
  }  // end loop Gauss points

  // bye
  return;
}
