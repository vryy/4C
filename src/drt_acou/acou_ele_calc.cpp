/*!
\file acou_ele_calc.cpp
\brief

<pre>
Maintainer: Svenja Schoeder
            schoeder@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15271
</pre>
*/
/*--------------------------------------------------------------------------*/

#include "acou_ele_calc.H"
#include "acou_ele_action.H"

#include "../drt_fem_general/drt_utils_boundary_integration.H"

#include "../drt_geometry/position_array.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"

#include "../drt_mat/acoustic.H"
#include "../drt_mat/scatra_mat.H"

#include <Teuchos_TimeMonitor.hpp>

namespace
{
  void zeroMatrix (Epetra_SerialDenseMatrix &mat)
  {
    std::memset(mat.A(), 0, sizeof(double)*mat.M()*mat.N());
  }
}


/*----------------------------------------------------------------------*
 * Constructor
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::AcouEleCalc():
    localSolver_(shapes_)
{
  interiorGradnp_.Shape(ndofs_*nsd_,1);
  interiorVelnp_.Shape(ndofs_,1);
  interiorGradn_.Shape(ndofs_*nsd_,1);
  interiorVeln_.Shape(ndofs_,1);
  interiorGradnm_.Shape(ndofs_*nsd_,1);
  interiorVelnm_.Shape(ndofs_,1);
  interiorGradnmm_.Shape(ndofs_*nsd_,1);
  interiorVelnmm_.Shape(ndofs_,1);
  interiorGradnmmm_.Shape(ndofs_*nsd_,1);
  interiorVelnmmm_.Shape(ndofs_,1);
}


/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elemat1_epetra,
                                                 Epetra_SerialDenseMatrix&  elemat2_epetra,
                                                 Epetra_SerialDenseVector&  elevec1_epetra,
                                                 Epetra_SerialDenseVector&  elevec2_epetra,
                                                 Epetra_SerialDenseVector&  elevec3_epetra,
                                                 const DRT::UTILS::GaussIntegration & ,
                                                 bool                       offdiag)
{
  return this->Evaluate( ele, discretization, lm, params, mat,
                         elemat1_epetra, elemat2_epetra,
                         elevec1_epetra, elevec2_epetra, elevec3_epetra,
                         offdiag );
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::Evaluate(DRT::ELEMENTS::Acou* ele,
                                                 DRT::Discretization    &   discretization,
                                                 const std::vector<int> &   lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elemat1,
                                                 Epetra_SerialDenseMatrix&  ,
                                                 Epetra_SerialDenseVector&  elevec1,
                                                 Epetra_SerialDenseVector&  elevec2,
                                                 Epetra_SerialDenseVector&  ,
                                                 bool                       offdiag)
{

  const ACOU::Action action = DRT::INPUT::get<ACOU::Action>(params,"action");

  switch(action)
  {
  case ACOU::project_field:
  {
    ProjectField(ele,params,mat,discretization,lm,elevec1,elevec2);
    break;
  }
  case ACOU::project_optical_field:
  {
    ProjectOpticalField(ele,params,mat,discretization,lm,elevec1,elevec2);
    break;
  }
  case ACOU::interpolate_hdg_to_node:
  {
    ReadGlobalVectors(*ele,discretization,lm);
    NodeBasedValues(ele,discretization,lm,elevec1);
    break;
  }
  case ACOU::update_secondary_solution:
  {
    bool adjoint = params.get<bool>("adjoint");
    double dt = params.get<double>("dt");

    ReadGlobalVectors(*ele,discretization,lm);
    shapes_.Evaluate(*ele);

    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    zeroMatrix(elemat1);
    zeroMatrix(localSolver_.Cmat);
    zeroMatrix(localSolver_.Emat);
    zeroMatrix(localSolver_.Gmat);
    localSolver_.ComputeInteriorMatrices(mat,dt);
    for (unsigned int face=0; face<nfaces_; ++face)
    {
      shapes_.EvaluateFace(*ele, face);
      localSolver_.ComputeFaceMatrices(face,mat,dt);
    }

    if(!adjoint)
    {
      if (dyna_ == INPAR::ACOU::acou_trapezoidal)
        UpdateInteriorVariablesTrap(discretization,*ele);
      else
        UpdateInteriorVariables(discretization,*ele);
    } // if(!adjoint)
    else
    {
      if (dyna_ == INPAR::ACOU::acou_trapezoidal)
        UpdateInteriorVariablesTrapAdjoint(discretization,*ele);
      else
        UpdateInteriorVariablesAdjoint(discretization,*ele);
    }
    break;
  }
  case ACOU::calc_abc:
  {
    int face = params.get<int>("face");
    shapes_.Evaluate(*ele);
    shapes_.EvaluateFace(*ele, face);
    // note: absorbing bcs are treated fully implicitly!
    localSolver_.ComputeAbsorbingBC(ele,params,mat,face,elemat1,elevec1);
    break;
  }
  case ACOU::calc_integr_objf:
  {
    if(ele->Owner() == discretization.Comm().MyPID())
    {
      int face = params.get<int>("face");
      double dt = params.get<double>("dt");
      shapes_.Evaluate(*ele);
      shapes_.EvaluateFace(*ele, face);
      localSolver_.ComputeObjfIntegral(ele,params,mat,face,dt);
    }
    break;
  }
  case ACOU::calc_systemmat_and_residual:
  {
    const bool resonly = params.get<bool>("resonly");
    bool adjoint = params.get<bool>("adjoint");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    shapes_.Evaluate(*ele);

    zeroMatrix(elevec1);

    zeroMatrix(elemat1);
    zeroMatrix(localSolver_.Cmat);
    zeroMatrix(localSolver_.Emat);
    zeroMatrix(localSolver_.Gmat);
    localSolver_.ComputeInteriorMatrices(mat,dt);
    for (unsigned int face=0; face<nfaces_; ++face)
    {
      shapes_.EvaluateFace(*ele, face);
      localSolver_.ComputeFaceMatrices(face,mat,dt);
    }

    ReadGlobalVectors(*ele,discretization,lm);

    if(!adjoint)
    {
      if(!resonly)
      {
        if (dyna_ == INPAR::ACOU::acou_trapezoidal)
          localSolver_.CondenseLocalPartTrap(elemat1,mat,dt);
        else
          localSolver_.CondenseLocalPart(elemat1,mat);
      }

      if (dyna_ == INPAR::ACOU::acou_trapezoidal)
        localSolver_.ComputeResidualTrap(elevec1,mat,interiorGradnp_,interiorVelnp_,traceVal_,dt);
      else
        localSolver_.ComputeResidual(elevec1,mat,interiorGradnp_,interiorVelnp_);

    } // if(!adjoint)
    else
    {
      if(!resonly)
      {
        if (dyna_ == INPAR::ACOU::acou_trapezoidal)
          localSolver_.CondenseLocalPartTrapAdjoint(elemat1,mat);
        else
          localSolver_.CondenseLocalPartAdjoint(elemat1,mat);
      }

      if (dyna_ == INPAR::ACOU::acou_trapezoidal)
        localSolver_.ComputeResidualTrapAdjoint(elevec1,mat,interiorGradnp_,interiorVelnp_,traceVal_,dt);
      else
        localSolver_.ComputeResidualAdjoint(elevec1,mat,interiorGradnp_,interiorVelnp_);
    }
    break;
  }
  case ACOU::update_secondary_solution_and_calc_residual:
  {
    bool adjoint = params.get<bool>("adjoint");
    bool errormaps = params.get<bool>("errormaps");
    double dt = params.get<double>("dt");
    dyna_ = params.get<INPAR::ACOU::DynamicType>("dynamic type");

    ReadGlobalVectors(*ele,discretization,lm);
    shapes_.Evaluate(*ele);

    zeroMatrix(elevec1);
    zeroMatrix(localSolver_.Cmat);
    zeroMatrix(localSolver_.Emat);
    zeroMatrix(localSolver_.Gmat);
    localSolver_.ComputeInteriorMatrices(mat,dt);
    for (unsigned int face=0; face<nfaces_; ++face)
    {
      shapes_.EvaluateFace(*ele, face);
      localSolver_.ComputeFaceMatrices(face,mat,dt);
    }

    if(!adjoint)
    {
      if (dyna_ == INPAR::ACOU::acou_trapezoidal)
        UpdateInteriorVariablesAndComputeResidualTrap(discretization,params,*ele,dt,elevec1,mat,errormaps);
      else
        UpdateInteriorVariablesAndComputeResidual(discretization,params,*ele,elevec1,mat,dt,errormaps); // the standard case
    } // if(!adjoint)
    else
    {
      if (dyna_ == INPAR::ACOU::acou_trapezoidal)
        UpdateInteriorVariablesAndComputeResidualTrapAdjoint(discretization,*ele,elevec1,mat); // the standard case
      else
        UpdateInteriorVariablesAndComputeResidualAdjoint(discretization,*ele,elevec1,mat); // the standard case
    } // else ** if(!adjoint)
    break;
  }
  default:
    dserror("unknown action supplied");
    break;
  } // switch(action)
  return 0;
}

/*----------------------------------------------------------------------*
 * ReadGlobalVectors
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
ReadGlobalVectors(const DRT::Element     & ele,
                  DRT::Discretization    & discretization,
                  const std::vector<int> & lm)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ReadGlobalVectors");

  // read the HDG solution vector (for traces)
  traceVal_.resize(nfaces_*nfdofs_);
  dsassert(lm.size() == traceVal_.size(), "Internal error");
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState("trace");
  DRT::UTILS::ExtractMyValues(*matrix_state,traceVal_,lm);

  traceValm_.resize(nfaces_*nfdofs_);
  if(discretization.HasState("trace_m"))
  {
    Teuchos::RCP<const Epetra_Vector> matrix_state_m = discretization.GetState("trace_m");
    DRT::UTILS::ExtractMyValues(*matrix_state_m,traceValm_,lm);
  }

  {
    interiorValnp_.resize(ndofs_*(nsd_+1));
    Teuchos::RCP<const Epetra_Vector> intvel = discretization.GetState(1,"intvel");
    std::vector<int> localDofs1 = discretization.Dof(1, &ele);
    DRT::UTILS::ExtractMyValues(*intvel,interiorValnp_,localDofs1);

    // now write this in corresponding interiorGradnp_ and interiorVelnp_
    for(unsigned int i=0; i<interiorValnp_.size(); ++i)
    {
      if ((i+1)%(nsd_+1) == 0)
      {
        interiorVelnp_((i+1)/(nsd_+1)-1) = interiorValnp_[i];
      }
      else
      {
        // careful: interiorValnp is sorted as [ Q_1x Q_1y Q_1z V_1 Q_2x Q_2y Q_2z V_2 ... V_ndofs ]
        // interiorVelnp is sorted as [ V_1 V_2 V_3 ... V_ndofs]
        // interiorGradnp should be sorted as [ Q_1x Q_2x Q_3x Q_4x ... Q_ndofsx Q_1y Q_2y ... ]
        int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
        interiorGradnp_(xyz*ndofs_+i/(nsd_+1)) = interiorValnp_[i];
      }
    }

    if(discretization.HasState(1,"intvelm")) // bdf2 and bdf3 need this
    {
      interiorValnm_.resize(ndofs_*(nsd_+1));
      Teuchos::RCP<const Epetra_Vector> intvelm = discretization.GetState(1,"intvelm");
      std::vector<int> localDofs1 = discretization.Dof(1, &ele);
      DRT::UTILS::ExtractMyValues(*intvelm,interiorValnm_,localDofs1);

      // now write this in corresponding interiorGradnp_ and interiorVelnp_
      for(unsigned int i=0; i<interiorValnm_.size(); ++i)
      {
        if ((i+1)%(nsd_+1) == 0)
        {
          interiorVelnm_((i+1)/(nsd_+1)-1) = interiorValnm_[i];
        }
        else
        {
          int xyz = i % (nsd_+1);
          interiorGradnm_(xyz*ndofs_+i/(nsd_+1)) = interiorValnm_[i];
        }
      }

      if(discretization.HasState(1,"intvelmm")) // bdf3 needs this
      {
        interiorValnmm_.resize(ndofs_*(nsd_+1));
        Teuchos::RCP<const Epetra_Vector> intvelmm = discretization.GetState(1,"intvelmm");
        std::vector<int> localDofs1 = discretization.Dof(1, &ele);
        DRT::UTILS::ExtractMyValues(*intvelmm,interiorValnmm_,localDofs1);

        // now write this in corresponding interiorGradnp_ and interiorVelnp_
        for(unsigned int i=0; i<interiorValnmm_.size(); ++i)
        {
          if ((i+1)%(nsd_+1) == 0)
          {
            interiorVelnmm_((i+1)/(nsd_+1)-1) = interiorValnmm_[i];
          }
          else
          {
            int xyz = i % (nsd_+1);
            interiorGradnmm_(xyz*ndofs_+i/(nsd_+1)) = interiorValnmm_[i];
          }
        }
        if(discretization.HasState(1,"intvelmmm")) // bdf4 needs this
        {
          interiorValnmmm_.resize(ndofs_*(nsd_+1));
          Teuchos::RCP<const Epetra_Vector> intvelmmm = discretization.GetState(1,"intvelmmm");
          std::vector<int> localDofs1 = discretization.Dof(1, &ele);
          DRT::UTILS::ExtractMyValues(*intvelmmm,interiorValnmmm_,localDofs1);

          // now write this in corresponding interiorGradnp_ and interiorVelnp_
          for(unsigned int i=0; i<interiorValnmmm_.size(); ++i)
          {
            if ((i+1)%(nsd_+1) == 0)
            {
              interiorVelnmmm_((i+1)/(nsd_+1)-1) = interiorValnmmm_[i];
            }
            else
            {
              int xyz = i % (nsd_+1);
              interiorGradnmmm_(xyz*ndofs_+i/(nsd_+1)) = interiorValnmmm_[i];
            }
          } // if(discretization.HasState(1,"intvelmmm"))
        } // if(discretization.HasState("intvelmm"))
      } // if(discretization.HasState("intvelmm"))
    } // if(discretization.HasState("intvelm"))
  }

  return;
} // ReadGlobalVectors


/*----------------------------------------------------------------------*
 * ProjectField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::ProjectField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Teuchos::RCP<MAT::Material>&         mat,
    DRT::Discretization&                 discretization,
    const std::vector<int>&              lm,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{

  shapes_.Evaluate(*ele);

  // reshape elevec2 as matrix
  dsassert(elevec2.M() == 0 ||
           elevec2.M() == (nsd_+1)*ndofs_, "Wrong size in project vector 2");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  // get function
  const int *start_func = params.getPtr<int>("startfuncno");
  double pulse = params.get<double>("pulse");

  // internal variables
  if (elevec2.M() > 0)
  {
    LINALG::Matrix<ndofs_,nsd_+1> localMat(true);
    localMat.PutScalar(0.);
    Epetra_SerialDenseMatrix massPart;
    massPart.Shape(ndofs_,ndofs_);

    for (unsigned int q=0; q<ndofs_; ++q )
    {
      const double fac = shapes_.jfac(q);
      const double sqrtfac = std::sqrt(fac);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_.xyzreal(d,q); // coordinates of quadrature point in real coordinates
      double u[nsd_];
      double p;
      dsassert(start_func != NULL,"startfuncno not set for initial value");
      EvaluateAll(*start_func, xyz, u, p, rho/pulse); // u and p at quadrature point

      // now fill the components in the one-sided mass matrix and the right hand side
      for (unsigned int i=0; i<ndofs_; ++i)
      {
        massPart(i,q) = shapes_.shfunct(i,q) * sqrtfac;
        for (unsigned int d=0; d<nsd_; ++d)
          localMat(i,d) += shapes_.shfunct(i,q) * u[d] * fac;
        localMat(i,nsd_) += shapes_.shfunct(i,q) * p * fac;
      }
    }
    Epetra_SerialDenseMatrix massMat;
    massMat.Shape(ndofs_,ndofs_);
    massMat.Multiply('N','T',1.,massPart,massPart,0.);
    {
      LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_,nsd_+1> inverseMass;
      LINALG::Matrix<ndofs_,ndofs_> mass(massMat,true);
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(localMat,localMat);
      inverseMass.Solve();
    }

    for (unsigned int r=0; r<ndofs_; ++r )
    {
      elevec2[r*(nsd_+1)+nsd_] += localMat(r,nsd_); // pressure
      for (unsigned int i=0;i<nsd_;++i)
      {
        elevec2[r*(nsd_+1)+i] += localMat(r,i); // velocity
      }
    }
  }

  // trace variable
  LINALG::Matrix<nfdofs_,nfdofs_> mass;
  LINALG::Matrix<nfdofs_,1> trVec;
  dsassert(elevec1.M() == nfaces_*nfdofs_, "Wrong size in project vector 1");

  for (unsigned int face=0; face<nfaces_; ++face)
  {
    shapes_.EvaluateFace(*ele, face);
    mass.PutScalar(0.);
    trVec.PutScalar(0.);

    for (unsigned int q=0; q<nfdofs_; ++q)
    {
      const double fac = shapes_.jfacF(q);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_.xyzFreal(d,q);
      double u[nsd_];
      double p;
      EvaluateAll(*start_func, xyz, u, p, rho/pulse);

      // now fill the components in the mass matrix and the right hand side
      for (unsigned int i=0; i<nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j=0; j<nfdofs_; ++j)
          mass(i,j) += shapes_.shfunctF(i,q) * shapes_.shfunctF(j,q) * fac;
        trVec(i,0) += shapes_.shfunctF(i,q) * p * fac;
      }
    }

    LINALG::FixedSizeSerialDenseSolver<nfdofs_,nfdofs_,1> inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec,trVec);
    inverseMass.Solve();
    for (unsigned int i=0; i<nfdofs_; ++i)
      elevec1(face*nfdofs_+i) = trVec(i,0);
  }

  return 0;
}

/*----------------------------------------------------------------------*
 * ProjectOpticalField
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::AcouEleCalc<distype>::ProjectOpticalField(
    DRT::ELEMENTS::Acou*                 ele,
    Teuchos::ParameterList&              params,
    Teuchos::RCP<MAT::Material>&         mat,
    DRT::Discretization&                 discretization,
    const std::vector<int>&              lm,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2)
{

  bool meshconform = params.get<bool>("mesh conform");

  if(meshconform)
  {
    // this is the easy case: the corresponding optical element has exactly the same nodes
    // and hence we can easily get the scatra dofs of the nodes and substitute the pressure
    // values

    // get the absorption coefficient
    double absorptioncoeff = params.get<double>("absorption");

    Teuchos::RCP<std::vector<double> > nodevals = params.get<Teuchos::RCP<std::vector<double> > >("nodevals");
    int numlightnode = (*nodevals).size()/(nsd_+1);

    double lightxyz[numlightnode][nsd_];
    double values[numlightnode];

    for(int i=0; i<numlightnode; ++i)
    {
      for(unsigned int j=0; j<nsd_; ++j)
        lightxyz[i][j] = (*nodevals)[i*(nsd_+1)+j];
      values[i] = (*nodevals)[i*(nsd_+1)+nsd_];
    }

    // now we have locations "lightxyz" and corresponding phis in "values"
    // what we want to do now is to evaluate the field in the locations of the
    // Gauss points in this element, to assign an initial distribution
    // so we do exactly the same as in ProjectInitalField BUT call an other
    // evaluate function which will interpolate from given lightxyz and values

    shapes_.Evaluate(*ele);
    const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
    double rho = actmat->Density();

    double pulse = params.get<double>("pulse");

    // internal variables
    if (elevec2.M() > 0)
    {
      LINALG::Matrix<ndofs_,nsd_+1> localMat(true);
      localMat.PutScalar(0.);
      Epetra_SerialDenseMatrix massPart;
      massPart.Shape(ndofs_,ndofs_);

      for (unsigned int q=0; q<ndofs_; ++q )
      {
        const double fac = shapes_.jfac(q);
        const double sqrtfac = std::sqrt(fac);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapes_.xyzreal(d,q); // coordinates of quadrature point in real coordinates
        double p=0.0;
        EvaluateLight(lightxyz,values,numlightnode, xyz, p, rho/pulse, absorptioncoeff); // p at quadrature point

        for (unsigned int i=0; i<ndofs_; ++i)
        {
          massPart(i,q) = shapes_.shfunct(i,q) * sqrtfac;
          localMat(i,nsd_) += shapes_.shfunct(i,q) * p * fac;
        }
      }
      Epetra_SerialDenseMatrix massMat;
      massMat.Shape(ndofs_,ndofs_);
      massMat.Multiply('N','T',1.,massPart,massPart,0.);
      {
        LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_,nsd_+1> inverseMass;
        LINALG::Matrix<ndofs_,ndofs_> mass(massMat,true);
        inverseMass.SetMatrix(mass);
        inverseMass.SetVectors(localMat,localMat);
        inverseMass.Solve();
      }

      for (unsigned int r=0; r<ndofs_; ++r )
      {
        elevec2[r*(nsd_+1)+nsd_] += localMat(r,nsd_); // pressure
      }
    }

    // trace variable
    LINALG::Matrix<nfdofs_,nfdofs_> mass;
    LINALG::Matrix<nfdofs_,1> trVec;
    dsassert(elevec1.M() == nfaces_*nfdofs_, "Wrong size in project vector 1");

    for (unsigned int face=0; face<nfaces_; ++face)
    {
      shapes_.EvaluateFace(*ele, face);

      mass.PutScalar(0.);
      trVec.PutScalar(0.);

      for (unsigned int q=0; q<nfdofs_; ++q)
      {
        const double fac = shapes_.jfacF(q);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapes_.xyzFreal(d,q);
        double p = 0.0;

        EvaluateLight(lightxyz,values,numlightnode, xyz, p, rho/pulse, absorptioncoeff); // u and p at quadrature point

        for (unsigned int i=0; i<nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j=0; j<nfdofs_; ++j)
            mass(i,j) += shapes_.shfunctF(i,q) * shapes_.shfunctF(j,q) * fac;
          trVec(i,0) += shapes_.shfunctF(i,q) * p * fac;
        }
      }

      LINALG::FixedSizeSerialDenseSolver<nfdofs_,nfdofs_,1> inverseMass;
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(trVec,trVec);
      inverseMass.Solve();

      for (unsigned int i=0; i<nfdofs_; ++i)
        elevec1(face*nfdofs_+i) = trVec(i,0);
    }
  }
  else
    dserror("non conforming meshes not yet implemented");

  return 0;
}

/*----------------------------------------------------------------------*
 * EvaluateAll
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::EvaluateAll(const int start_func,
    const double (&xyz)[nsd_],
    double (&u)[nsd_],
    double  &p,
    double rho) const
{
  p = DRT::Problem::Instance()->Funct(start_func-1).Evaluate(0,xyz,0.0,NULL);

  // how to determine the spatial derivative for initial velocity field?
  /*
  double eps=1e-5;
  double xyz_eps[nsd_];
  for(unsigned int i=0; i<nsd_; ++i)
  {
    for(unsigned int j=0; j<nsd_; ++j)
      xyz_eps[j]=xyz[j];
    xyz_eps[i]+=eps;
    double p_eps = DRT::Problem::Instance()->Funct(start_func-1).Evaluate(0,xyz_eps,0.0,NULL);
    ///u[i] = 0.0;
    u[i]=-1/rho*(p-p_eps)/eps;
  }
  */
  return;
}

/*----------------------------------------------------------------------*
 * EvaluateLight
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::EvaluateLight(
                   double lightxyz[][nsd_],
                   double values[],
                   int    numnode,
                   const double (&xyz)[nsd_],
                   double  &p,
                   double rho,
                   double absorptioncoeff
) const
{

  // interpolation from nodes
  if (distype == DRT::Element::quad4)
  {
    if(numnode!=4) dserror("wrong number of nodes given");

    //*******************************************************************
    LINALG::Matrix<4,4> coeff(true);
    for(int i=0; i<4; ++i)
    {
      coeff(i,0)=1.0;
      coeff(i,1)=lightxyz[i][0];
      coeff(i,2)=lightxyz[i][1];
      coeff(i,3)=lightxyz[i][0]*lightxyz[i][1];
    }
    LINALG::Matrix<4,4> coeff_N(true);
    coeff_N(0,0)=1.0;
    coeff_N(1,1)=1.0;
    coeff_N(2,2)=1.0;
    coeff_N(3,3)=1.0;
    {
      LINALG::FixedSizeSerialDenseSolver<4,4,4> inverseCoeff;
      inverseCoeff.SetMatrix(coeff);
      inverseCoeff.SetVectors(coeff_N,coeff_N);
      inverseCoeff.Solve();
    }

    p = ( coeff_N(0,0) + coeff_N(1,0) * xyz[0] + coeff_N(2,0) * xyz[1] + coeff_N(3,0) * xyz[0] * xyz[1] ) * values[0]
      + ( coeff_N(0,1) + coeff_N(1,1) * xyz[0] + coeff_N(2,1) * xyz[1] + coeff_N(3,1) * xyz[0] * xyz[1] ) * values[1]
      + ( coeff_N(0,2) + coeff_N(1,2) * xyz[0] + coeff_N(2,2) * xyz[1] + coeff_N(3,2) * xyz[0] * xyz[1] ) * values[2]
      + ( coeff_N(0,3) + coeff_N(1,3) * xyz[0] + coeff_N(2,3) * xyz[1] + coeff_N(3,3) * xyz[0] * xyz[1] ) * values[3];
    p *= -absorptioncoeff;

  }
  else
    dserror("not yet implemented"); // TODO

  return;
}

/*----------------------------------------------------------------------*
 * NodeBasedValues
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::NodeBasedValues(
                          DRT::ELEMENTS::Acou*                 ele,
                          DRT::Discretization&                 discretization,
                          const std::vector<int>&              lm,
                          Epetra_SerialDenseVector&            elevec1)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::NodeBasedValues");

  dsassert(elevec1.M() == (int)nen_*(nsd_+2)+1, "Vector does not have correct size");
  elevec1.Scale(0.0);
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);
  LINALG::Matrix<1,ndofs_> values;

  double cellpress = 0.0;
  for (unsigned int i=0; i<ndofs_; ++i)
    cellpress += interiorVelnp_(i);
  cellpress /= double(ndofs_);
  elevec1((nsd_+2)*nen_) = cellpress;

  for (unsigned int i=0; i<nen_; ++i)
  {
    // evaluate shape polynomials in node
    for (unsigned int idim=0;idim<nsd_;idim++)
      shapes_.xsi(idim) = locations(idim,i);
    shapes_.polySpace_.Evaluate(shapes_.xsi,values);

    // compute values for velocity and pressure by summing over all basis functions
    double sum = 0;
    for (unsigned int k=0; k<ndofs_; ++k)
      sum += values(k) * interiorVelnp_(k);
    elevec1(nsd_*nen_+i) = sum;

    for (unsigned int d=0; d<nsd_; ++d)
    {
      sum = 0.0;
      for (unsigned int k=0; k<ndofs_; ++k)
        sum += values(k) * interiorGradnp_(d*ndofs_+k);
      elevec1(d*nen_+i) = sum;
    }
  }

  LINALG::Matrix<1,nfdofs_> fvalues;
  for (unsigned int face=0; face<nfaces_; ++face)
  {
    // const int * fnodeIds = ele->Faces()[face]->NodeIds();
    for (int i=0; i<DRT::UTILS::DisTypeToNumNodePerFace<distype>::numNodePerFace; ++i)
    {
      // evaluate shape polynomials in node
      for (unsigned int idim=0;idim<nsd_-1;idim++)
        shapes_.xsiF(idim) = locations(idim,i);
      shapes_.polySpaceFace_.Evaluate(shapes_.xsiF,fvalues); // TODO: fix face orientation here

      // compute values for velocity and pressure by summing over all basis functions
      double sum = 0;
      for (unsigned int k=0; k<nfdofs_; ++k)
        sum += fvalues(k) * traceVal_[face*nfdofs_+k];
      if(elevec1((nsd_+1)*nen_+shapes_.faceNodeOrder[face][i]) == 0.0)
        elevec1((nsd_+1)*nen_+shapes_.faceNodeOrder[face][i]) = sum;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * Instance
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype> *
DRT::ELEMENTS::AcouEleCalc<distype>::Instance( bool create )
{
  static AcouEleCalc<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new AcouEleCalc<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 * Done
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}

/*----------------------------------------------------------------------*
 * Constructor ShapeValues
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::ShapeValues::
ShapeValues()
:
polySpace_(DRT::ELEMENTS::Acou::degree),
polySpaceFace_(DRT::ELEMENTS::Acou::degree)
{
  dsassert(polySpace_.Size() == ndofs_, "Wrong polynomial space constructed");
  dsassert(polySpaceFace_.Size() == nfdofs_, "Wrong polynomial space constructed");

  const int degree = DRT::ELEMENTS::Acou::degree;
  LINALG::Matrix<1,ndofs_> values;
  LINALG::Matrix<nsd_,ndofs_> derivs;
  LINALG::Matrix<1,nfdofs_> faceValues;

  // this causes an error for singleton Done()
  // DRT::UTILS::GaussPointCache cache;
  // quadrature_ = cache.Create(distype, degree*2);
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
  // fquadrature_ = cache.Create(facedis, degree*2);
  quadrature_ = DRT::UTILS::GaussPointCache::Instance().Create(distype, degree*2);
  fquadrature_ = DRT::UTILS::GaussPointCache::Instance().Create(facedis, degree*2);

  shfunct.Shape(ndofs_,ndofs_);
  shfunctAvg.Resize(ndofs_);
  shderiv.Shape(ndofs_*nsd_,ndofs_);
  shderxy.Shape(ndofs_*nsd_,ndofs_);
  jfac.Resize(ndofs_);

  dsassert(static_cast<unsigned int>(quadrature_->NumPoints())==ndofs_, "Internal error - not implemented");
  for (unsigned int q=0; q<ndofs_; ++q )
  {
    // gauss point in real coordinates
    const double* gpcoord = quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    polySpace_.Evaluate(xsi,values);
    polySpace_.Evaluate_deriv1(xsi,derivs);

    for (unsigned int i=0; i<ndofs_; ++i)
    {
      shfunct(i,q) = values(i);
      for (unsigned int d=0; d<nsd_; ++d)
        shderiv(i*nsd_+d,q) = derivs(d,i);
    }

    LINALG::Matrix<nen_,1> myfunct(funct.A()+q*nen_,true);
    DRT::UTILS::shape_function<distype>(xsi,myfunct);
  }

  shfunctFNoPermute.Shape(nfdofs_, nfdofs_);
  shfunctF.Shape(nfdofs_, nfdofs_);
  shfunctI.resize(nfaces_);
  for (unsigned int f=0;f<nfaces_; ++f)
    shfunctI[f].Shape(ndofs_, nfdofs_);
  normals.Shape(nsd_, nfdofs_);
  jfacF.Resize(nfdofs_);

  for (unsigned int q=0; q<nfdofs_; ++q )
  {
    const double* gpcoord = fquadrature_->Point(q);

    const unsigned int codim = nsd_-1;
    for (unsigned int idim=0;idim<codim;idim++)
      xsiF(idim) = gpcoord[idim];

    polySpaceFace_.Evaluate(xsiF,faceValues);
    for (unsigned int i=0; i<nfdofs_; ++i)
      shfunctFNoPermute(i,q) = faceValues(i);

    LINALG::Matrix<nfn_,1> myfunct(functF.A()+q*nfn_,true);
    DRT::UTILS::shape_function<facedis>(xsiF,myfunct);
  }

  Epetra_SerialDenseMatrix quadrature(nfdofs_,nsd_,false);
  Epetra_SerialDenseMatrix trafo(nsd_,nsd_,false);
  for (unsigned int f=0; f<nfaces_; ++f)
  {
    DRT::UTILS::BoundaryGPToParentGP<nsd_>(quadrature,trafo,*fquadrature_,distype,
                                           facedis, f);
    for (unsigned int q=0; q<nfdofs_; ++q)
    {
      for (unsigned int d=0; d<nsd_; ++d)
        xsi(d) = quadrature(q,d);
      polySpace_.Evaluate(xsi,values);
      for (unsigned int i=0; i<ndofs_; ++i)
        shfunctI[f](i,q) = values(i);
    }
  }

  if (nsd_ == 2)
    faceNodeOrder = DRT::UTILS::getEleNodeNumberingLines(distype);
  else if (nsd_ == 3)
    faceNodeOrder = DRT::UTILS::getEleNodeNumberingSurfaces(distype);
  else
    dserror("Not implemented for dim != 2, 3");

}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::ShapeValues::Evaluate (const DRT::Element &ele)
{
  dsassert(ele.Shape() == distype, "Internal error");
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(&ele,xyze);

  for (unsigned int i=0; i<ndofs_; ++i)
    shfunctAvg(i) = 0.;
  double faceVol = 0.;

  // evaluate geometry
  for (unsigned int q=0; q<ndofs_; ++q) {
    const double* gpcoord = quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,xyze);
    jfac(q) = xji.Invert(xjm) * quadrature_->Weight(q);

    LINALG::Matrix<nen_,1> myfunct(funct.A()+q*nen_,true);
    LINALG::Matrix<nsd_,1> mypoint(xyzreal.A()+q*nsd_,true);
    mypoint.MultiplyNN(xyze,myfunct);

    // transform shape functions
    for (unsigned int i=0; i<ndofs_; ++i)
      for (unsigned int d=0; d<nsd_; ++d) {
        shderxy(i*nsd_+d,q) = xji(d,0) * shderiv(i*nsd_,q);
        for (unsigned int e=1; e<nsd_; ++e)
          shderxy(i*nsd_+d,q) += xji(d,e) * shderiv(i*nsd_+e,q);
      }

    for (unsigned int i=0; i<ndofs_; ++i)
      shfunctAvg(i) += shfunct(i,q) * jfac(q);
    faceVol += jfac(q);
  }
  faceVol = 1./faceVol;
  for (unsigned int i=0; i<ndofs_; ++i)
    shfunctAvg(i) *= faceVol;
}


/*----------------------------------------------------------------------*
 * EvaluateFace
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void
DRT::ELEMENTS::AcouEleCalc<distype>::ShapeValues::
EvaluateFace (const DRT::Element &ele,
              const unsigned int  face)
{
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;

  // get face position array from element position array
  dsassert(faceNodeOrder[face].size() == nfn_,
           "Internal error");
  for (unsigned int i=0; i<nfn_; ++i)
    for (unsigned int d=0; d<nsd_; ++d)
    {
      xyzeF(d,i) = xyze(d,faceNodeOrder[face][i]);
    }

  // evaluate geometry
  for (unsigned int q=0; q<nfdofs_; ++q)
  {
    const double* gpcoord = fquadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_-1;idim++)
      xsiF(idim) = gpcoord[idim];

    DRT::UTILS::shape_function_deriv1<facedis>(xsiF,derivF);
    double jacdet = 0;
    DRT::UTILS::ComputeMetricTensorForBoundaryEle<facedis>(xyzeF,derivF,metricTensor,jacdet,&normal);
    for (unsigned int d=0; d<nsd_; ++d)
      normals(d,q) = normal(d);
    jfacF(q) = jacdet * fquadrature_->Weight(q);

    LINALG::Matrix<nfn_,1> myfunct(functF.A()+q*nfn_,true);
    LINALG::Matrix<nsd_,1> mypoint(xyzFreal.A()+q*nsd_,true);
    mypoint.MultiplyNN(xyzeF,myfunct);
  }

  // figure out how to permute face indices by checking permutation of nodes
  const int * nodeIds = ele.NodeIds();
  const int * fnodeIds = ele.Faces()[face]->NodeIds();
  const int ndofs1d = DRT::ELEMENTS::Acou::degree+1;
  // easy case: standard orientation
  bool standard = true;
  for (unsigned int i=0; i<nfn_; ++i)
    if (nodeIds[faceNodeOrder[face][i]] != fnodeIds[i])
      standard = false;

  if (standard)
  {
    //std::cout << "standard orientation" << std::endl;
    for (unsigned int q=0; q<nfdofs_; ++q)
      for (unsigned int i=0; i<nfdofs_; ++i)
        shfunctF(i,q) = shfunctFNoPermute(i,q);
  }
  // OK, the orientation is different from what I expect. see if we can find it
  else switch (nsd_)
  {
  case 2:
    // face flipped is the only case
    {
      dsassert(nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
               nodeIds[faceNodeOrder[face][0]] == fnodeIds[1], "Unknown face orientation in 2D");
      for (unsigned int i=0; i<nfdofs_; ++i) {
        //std::cout << nfdofs_-1-i << " ";
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(nfdofs_-1-i,q);
        }
      //std::cout << std::endl;
    }
    break;
  case 3:
    dsassert(distype == DRT::Element::hex8 ||
             distype == DRT::Element::hex20 ||
             distype == DRT::Element::hex27 ||
             distype == DRT::Element::nurbs8 ||
             distype == DRT::Element::nurbs27,
             "Not implemented for given shape");
    if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[1] &&
        nodeIds[faceNodeOrder[face][1]] == fnodeIds[0] &&
        nodeIds[faceNodeOrder[face][2]] == fnodeIds[3] &&
        nodeIds[faceNodeOrder[face][3]] == fnodeIds[2])    // x-direction mirrored
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = ndofs1d-1-ax + ay * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[0] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[3] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[1])    // permute x and y
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = ay + ax * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[3] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[1] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[0])    // y mirrored
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = ax + (ndofs1d-1-ay) * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[3] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[0] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[1])    // x and y mirrored
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = (ndofs1d-1-ax) + (ndofs1d-1-ay) * ndofs1d;
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else if (nodeIds[faceNodeOrder[face][0]] == fnodeIds[2] &&
             nodeIds[faceNodeOrder[face][1]] == fnodeIds[1] &&
             nodeIds[faceNodeOrder[face][2]] == fnodeIds[0] &&
             nodeIds[faceNodeOrder[face][3]] == fnodeIds[3])    // x and y mirrored and permuted
    {
      for (unsigned int i=0; i<nfdofs_; ++i) {
        const int ax = i%ndofs1d;
        const int ay = i/ndofs1d;
        int permute = (ndofs1d-1-ay) + (ndofs1d-1-ax) * ndofs1d; // note that this is for lexicographic ordering
        for (unsigned int q=0; q<nfdofs_; ++q)
          shfunctF(i,q) = shfunctFNoPermute(permute,q);
      }
    }
    else
    {
      for (unsigned int i=0; i<4; ++i)
        std::cout << nodeIds[faceNodeOrder[face][i]] << " " << fnodeIds[i] << "   ";
      std::cout << std::endl << std::flush;
      dserror("Unknown face orientation in 3D");
    }
    break;
  default:
    dserror("Only implemented in 2D and 3D");
    break;
  }
} // EvaluateFace

/*----------------------------------------------------------------------*
 * Constructor LocalSolver
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
LocalSolver(const ShapeValues &shapeValues)
:
shapes_(shapeValues)
{
  // shape all matrices
  Amat.Shape(nsd_*ndofs_,nsd_*ndofs_);
  invAmat.Shape(nsd_*ndofs_,nsd_*ndofs_);
  Bmat.Shape(ndofs_,nsd_*ndofs_);
  Hmat.Shape(nsd_*ndofs_,ndofs_);
  Dmat.Shape(ndofs_,ndofs_);
  Mmat.Shape(ndofs_,ndofs_);
  DmMmat.Shape(ndofs_,ndofs_);

  Cmat.Shape(nsd_*ndofs_,nfaces_*nfdofs_);
  Imat.Shape(nfaces_*nfdofs_,nsd_*ndofs_);
  Emat.Shape(ndofs_,nfaces_*nfdofs_);
  Jmat.Shape(nfaces_*nfdofs_,ndofs_);
  Gmat.Shape(nfaces_*nfdofs_,nfaces_*nfdofs_);
}

/*----------------------------------------------------------------------*
 * UpdateInteriorVariables
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariables(DRT::Discretization &                discretization,
                        DRT::ELEMENTS::Acou &                ele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariables");

  interiorGradn_ = interiorGradnp_;
  interiorVeln_ = interiorVelnp_;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Teuchos::RCP<MAT::Material> mat = ele.Material();
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector tempVec1(nsd_*ndofs_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradn_,0.0);
  tempVec1.Multiply('N','N',1.0,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2;
  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVeln_,0.0);
  tempVec2.Multiply('N','N',-1.0,localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',1.0/rho,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2.Multiply('N','T',1.0,tempMat1,localSolver_.Bmat,0.0);
  tempMat2 += localSolver_.Dmat;

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',1.0,localSolver_.Bmat,interiorVelnp_,1.0);
  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  return;
} // UpdateInteriorVariables


/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesTrap
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesTrap(DRT::Discretization &                discretization,
                            DRT::ELEMENTS::Acou &                ele)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariablesTrap");

  interiorGradn_ = interiorGradnp_;
  interiorVeln_ = interiorVelnp_;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Epetra_SerialDenseVector traceVal_SDV_m(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV_m(i) = traceValm_[i];

  double theta = 0.66;

  Teuchos::RCP<MAT::Material> mat = ele.Material();
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector tempVec1(ndofs_*nsd_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradnp_,0.0);
  tempVec1.Multiply('T','N',1.0-theta,localSolver_.Bmat,interiorVelnp_,1.0);
  tempVec1.Multiply('N','N',1.0-theta,localSolver_.Cmat,traceVal_SDV_m,1.0);
  tempVec1.Multiply('N','N',theta,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec2.Multiply('N','N',-(1.0-theta)/rho,localSolver_.Bmat,interiorGradnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.DmMmat,interiorVelnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.Emat,traceVal_SDV_m,1.0);
  tempVec2.Multiply('N','N',-(theta),localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',theta/rho,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 += localSolver_.DmMmat;
  tempMat2.Scale(theta);
  tempMat2.Multiply('N','T',theta,tempMat1,localSolver_.Bmat,1.0);
  tempMat2 += localSolver_.Mmat;

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',theta,localSolver_.Bmat,interiorVelnp_,1.0);
  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  return;
} // UpdateInteriorVariablesTrap

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesAdjoint(DRT::Discretization &                discretization,
                               DRT::ELEMENTS::Acou &                ele)
{
  interiorGradn_ = interiorGradnp_;
  interiorVeln_ = interiorVelnp_;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Teuchos::RCP<MAT::Material> mat = ele.Material();
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector tempVec1(nsd_*ndofs_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradn_,0.0);
  tempVec1.Multiply('N','N',-1.0/rho,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2;
  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVeln_,0.0);
  tempVec2.Multiply('N','N',-1.0,localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',-1.0,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2.Multiply('N','T',-1.0/rho,tempMat1,localSolver_.Bmat,0.0);
  tempMat2 += localSolver_.Dmat;

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',-1.0/rho,localSolver_.Bmat,interiorVelnp_,1.0);
  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  return;
} // UpdateInteriorVariablesAdjoint

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesTrapAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesTrapAdjoint(DRT::Discretization &            discretization,
                                   DRT::ELEMENTS::Acou &            ele)
{
  interiorGradn_ = interiorGradnp_;
  interiorVeln_ = interiorVelnp_;

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  Epetra_SerialDenseVector traceVal_SDV_m(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV_m(i) = traceValm_[i];

  double theta = 0.66;

  Teuchos::RCP<MAT::Material> mat = ele.Material();
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector tempVec1(ndofs_*nsd_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradnp_,0.0);
  tempVec1.Multiply('T','N',-(1.0-theta)/rho,localSolver_.Bmat,interiorVelnp_,1.0);
  tempVec1.Multiply('N','N',-(1.0-theta)/rho,localSolver_.Cmat,traceVal_SDV_m,1.0);
  tempVec1.Multiply('N','N',-theta/rho,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec2.Multiply('N','N',(1.0-theta),localSolver_.Bmat,interiorGradnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.DmMmat,interiorVelnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.Emat,traceVal_SDV_m,1.0);
  tempVec2.Multiply('N','N',-(theta),localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',-theta,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 += localSolver_.DmMmat;
  tempMat2.Scale(theta);
  tempMat2.Multiply('N','T',-theta/rho,tempMat1,localSolver_.Bmat,1.0);
  tempMat2 += localSolver_.Mmat;

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',-theta/rho,localSolver_.Bmat,interiorVelnp_,1.0);
  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  return;
} // UpdateInteriorVariablesTrapAdjoint

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesAndComputeResidual(DRT::Discretization &     discretization,
                                          Teuchos::ParameterList&    params,
                                          DRT::ELEMENTS::Acou &                ele,
                                          Epetra_SerialDenseVector          & elevec,
                                          const Teuchos::RCP<MAT::Material> &mat,
                                          double dt,
                                          bool errormaps)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariablesAndComputeResidual");

  Epetra_SerialDenseVector tempGradnp;
  Epetra_SerialDenseVector tempVelnp;
  if(dyna_ == INPAR::ACOU::acou_bdf2)
  {
    tempGradnp.Shape(ndofs_*nsd_,1); tempGradnp = interiorGradnp_;
    tempVelnp.Shape(ndofs_,1);       tempVelnp  = interiorVelnp_;
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradn_[i] = interiorGradnp_[i] * 4.0 / 3.0 - interiorGradnm_[i] / 3.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVeln_[i]  = interiorVelnp_[i]  * 4.0 / 3.0 - interiorVelnm_[i]  / 3.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf3)
  {
    tempGradnp.Shape(ndofs_*nsd_,1); tempGradnp = interiorGradnp_;
    tempVelnp.Shape(ndofs_,1);       tempVelnp  = interiorVelnp_;
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradn_[i] = interiorGradnp_[i] * 18.0 / 11.0 - interiorGradnm_[i] * 9.0 / 11.0 + interiorGradnmm_[i] * 2.0 / 11.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVeln_[i]  = interiorVelnp_[i]  * 18.0 / 11.0 - interiorVelnm_[i]  * 9.0 / 11.0 + interiorVelnmm_[i]  * 2.0 / 11.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf4)
  {
    tempGradnp.Shape(ndofs_*nsd_,1); tempGradnp = interiorGradnp_;
    tempVelnp.Shape(ndofs_,1);       tempVelnp  = interiorVelnp_;
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradn_[i] = interiorGradnp_[i] * 48.0 / 25.0 - interiorGradnm_[i] * 36.0 / 25.0 + interiorGradnmm_[i] * 16.0 / 25.0 - interiorGradnmmm_[i] * 3.0 / 25.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVeln_[i]  = interiorVelnp_[i]  * 48.0 / 25.0 - interiorVelnm_[i]  * 36.0 / 25.0 + interiorVelnmm_[i]  * 16.0 / 25.0 - interiorVelnmmm_[i]  * 3.0 / 25.0;
  }
  else
  {
    interiorGradn_ = interiorGradnp_;
    interiorVeln_ = interiorVelnp_;
  }
  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(nsd_*ndofs_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradn_,0.0);
  tempVec1.Multiply('N','N',1.0,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2;
  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVeln_,0.0);
  tempVec2.Multiply('N','N',-1.0,localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',1.0/rho,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2.Multiply('N','T',1.0,tempMat1,localSolver_.Bmat,0.0);
  tempMat2 += localSolver_.Dmat;

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',1.0,localSolver_.Bmat,interiorVelnp_,1.0);

  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  // sort this back to the interior values vector
  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  if (errormaps)
  {
    // first step: calculate the gradient we need (bold p in paper by Nguyen)
    Epetra_SerialDenseVector temp(ndofs_*nsd_);
    temp.Multiply('T','N',1.0,localSolver_.Bmat,interiorVelnp_,0.0);
    temp.Multiply('N','N',1.0,localSolver_.Cmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector p(ndofs_*nsd_);
    p.Multiply('N','N',1.0/dt,localSolver_.invAmat,temp,0.0);

    // second step: postprocess the pressure field: therefore we need some additional matrices!
    const unsigned int ndofspost = DRT::UTILS::FixedPower<DRT::ELEMENTS::Acou::degree+2,nsd_>::value;
    Epetra_SerialDenseMatrix H(ndofspost,ndofspost);
    Epetra_SerialDenseVector R(ndofspost);
    double err_p = CalculateError(ele,H,R,p);

    Teuchos::RCP<std::vector<double> > values = params.get<Teuchos::RCP<std::vector<double> > >("elevals");
    (*values)[ele.Id()] = err_p;
  }


  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  if(dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradnp_[i] = 4.0 / 3.0 * interiorGradnp_[i] - 1.0 / 3.0 * tempGradnp[i];
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVelnp_[i]  = 4.0 / 3.0 * interiorVelnp_[i]  - 1.0 / 3.0 * tempVelnp[i];
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradnp_[i] = 18.0 / 11.0 * interiorGradnp_[i] - 9.0 / 11.0 * tempGradnp[i] + 2.0 / 11.0 * interiorGradnm_[i];
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVelnp_[i]  = 18.0 / 11.0 * interiorVelnp_[i]  - 9.0 / 11.0 * tempVelnp[i]  + 2.0 / 11.0 * interiorVelnm_[i];
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradnp_[i] = 48.0 / 25.0 * interiorGradnp_[i] - 36.0 / 25.0 * tempGradnp[i] + 16.0 / 25.0 * interiorGradnm_[i] - 3.0 / 25.0 * interiorGradnmm_[i];
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVelnp_[i]  = 48.0 / 25.0 * interiorVelnp_[i]  - 36.0 / 25.0 * tempVelnp[i]  + 16.0 / 25.0 * interiorVelnm_[i]  - 3.0 / 25.0 * interiorVelnmm_[i];
  }

  tempVec1.Shape(ndofs_,1);
  tempVec1.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec1.Multiply('N','N',-1.0/rho,localSolver_.Bmat,interiorGradnp_,1.0);

  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0);

  elevec.Multiply('T','N',-1.0,localSolver_.Emat,tempVec2,0.0);

  tempVec1.Shape(ndofs_*nsd_,1);
  tempVec1.Multiply('T','N',1.0,localSolver_.Bmat,tempVec2,0.0);

  tempVec2.Shape(ndofs_*nsd_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);
  tempVec2 += interiorGradnp_;

  elevec.Multiply('T','N',-1.0/rho,localSolver_.Cmat,tempVec2,1.0);

  return;
} // UpdateInteriorVariablesAndComputeResidual

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidualTrap
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesAndComputeResidualTrap(DRT::Discretization &     discretization,
                                              Teuchos::ParameterList&    params,
                                              DRT::ELEMENTS::Acou &     ele,
                                              double                    dt,
                                              Epetra_SerialDenseVector          & elevec,
                                              const Teuchos::RCP<MAT::Material> &mat,
                                              bool errormaps)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariablesAndComputeResidualTrap");

  interiorGradn_ = interiorGradnp_;
  interiorVeln_ = interiorVelnp_;
  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];
  Epetra_SerialDenseVector traceVal_SDV_m(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV_m(i) = traceValm_[i];

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  double theta = 0.66;

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(ndofs_*nsd_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradnp_,0.0);
  tempVec1.Multiply('T','N',1.0-theta,localSolver_.Bmat,interiorVelnp_,1.0);
  tempVec1.Multiply('N','N',1.0-theta,localSolver_.Cmat,traceVal_SDV_m,1.0);
  tempVec1.Multiply('N','N',theta,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec2.Multiply('N','N',-(1.0-theta)/rho,localSolver_.Bmat,interiorGradnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.DmMmat,interiorVelnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.Emat,traceVal_SDV_m,1.0);
  tempVec2.Multiply('N','N',-(theta),localSolver_.Emat,traceVal_SDV,1.0);

  // now, we have to do the Schur complement thing, just as in "CondenseLocalPart" but without C and E
  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',theta/rho,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);


  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 += localSolver_.DmMmat;
  tempMat2.Scale(theta);
  tempMat2.Multiply('N','T',theta,tempMat1,localSolver_.Bmat,1.0);
  tempMat2 += localSolver_.Mmat;


  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',theta,localSolver_.Bmat,interiorVelnp_,1.0);

  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  // sort this back to the interior values vector
  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  // *****************************************************
  // local postprocessing to calculate error maps
  // *****************************************************

  if (errormaps)
  {
    // first step: calculate the gradient we need (bold p in paper by Nguyen)
    Epetra_SerialDenseVector temp(ndofs_*nsd_);
    temp.Multiply('T','N',1.0,localSolver_.Bmat,interiorVelnp_,0.0);
    temp.Multiply('N','N',1.0,localSolver_.Cmat,traceVal_SDV,1.0);
    Epetra_SerialDenseVector p(ndofs_*nsd_);
    p.Multiply('N','N',1.0/dt,localSolver_.invAmat,temp,0.0);

    // second step: postprocess the pressure field: therefore we need some additional matrices!
    const unsigned int ndofspost = DRT::UTILS::FixedPower<DRT::ELEMENTS::Acou::degree+2,nsd_>::value;
    Epetra_SerialDenseMatrix H(ndofspost,ndofspost);
    Epetra_SerialDenseVector R(ndofspost);
    double err_p = CalculateError(ele,H,R,p);

    Teuchos::RCP<std::vector<double> > values = params.get<Teuchos::RCP<std::vector<double> > >("elevals");
    (*values)[ele.Id()] = err_p;
  }

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  tempVec1.Shape(ndofs_,1);
  tempVec1.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec1.Multiply('N','N',-(1.0-theta)/rho,localSolver_.Bmat,interiorGradnp_,1.0);
  tempVec1.Multiply('N','N',-(1.0-theta),localSolver_.DmMmat,interiorVelnp_,1.0);
  tempVec1.Multiply('N','N',-(1.0-theta),localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector f(ndofs_*nsd_);
  f.Multiply('N','N',1.0,localSolver_.Amat,interiorGradnp_,0.0);
  f.Multiply('T','N',1.0-theta,localSolver_.Bmat,interiorVelnp_,1.0);
  f.Multiply('N','N',1.0-theta,localSolver_.Cmat,traceVal_SDV,1.0);

  tempVec1.Multiply('N','N',-1.0,tempMat1,f,1.0);

  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0);

  elevec.Multiply('T','N',-theta,localSolver_.Emat,tempVec2,0.0);

  tempVec1.Shape(ndofs_*nsd_,1);
  tempVec1.Multiply('T','N',theta,localSolver_.Bmat,tempVec2,0.0);
  tempVec1 += f;
  tempVec2.Shape(ndofs_*nsd_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  elevec.Multiply('T','N',-theta/rho,localSolver_.Cmat,tempVec2,1.0);

  elevec.Multiply('T','N',-(1.0-theta)/rho,localSolver_.Cmat,interiorGradnp_,1.0);
  elevec.Multiply('T','N',-(1.0-theta),localSolver_.Emat,interiorVelnp_,1.0);
  elevec.Multiply('N','N',-(1.0-theta),localSolver_.Gmat,traceVal_SDV,1.0);

  return;
} // UpdateInteriorVariablesAndComputeResidualTrap

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidualAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesAndComputeResidualAdjoint(DRT::Discretization &     discretization,
                                DRT::ELEMENTS::Acou &                ele,
                                Epetra_SerialDenseVector          & elevec,
                                const Teuchos::RCP<MAT::Material> &mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::UpdateInteriorVariablesAndComputeResidualAdjoint");

  Epetra_SerialDenseVector tempGradnp;
  Epetra_SerialDenseVector tempVelnp;
  if(dyna_ == INPAR::ACOU::acou_bdf2)
  {
    tempGradnp.Shape(ndofs_*nsd_,1); tempGradnp = interiorGradnp_;
    tempVelnp.Shape(ndofs_,1);       tempVelnp  = interiorVelnp_;
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradn_[i] = interiorGradnp_[i] * 4.0 / 3.0 - interiorGradnm_[i] / 3.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVeln_[i]  = interiorVelnp_[i]  * 4.0 / 3.0 - interiorVelnm_[i]  / 3.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf3)
  {
    tempGradnp.Shape(ndofs_*nsd_,1); tempGradnp = interiorGradnp_;
    tempVelnp.Shape(ndofs_,1);       tempVelnp  = interiorVelnp_;
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradn_[i] = interiorGradnp_[i] * 18.0 / 11.0 - interiorGradnm_[i] * 9.0 / 11.0 + interiorGradnmm_[i] * 2.0 / 11.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVeln_[i]  = interiorVelnp_[i]  * 18.0 / 11.0 - interiorVelnm_[i]  * 9.0 / 11.0 + interiorVelnmm_[i]  * 2.0 / 11.0;
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf4)
  {
    tempGradnp.Shape(ndofs_*nsd_,1); tempGradnp = interiorGradnp_;
    tempVelnp.Shape(ndofs_,1);       tempVelnp  = interiorVelnp_;
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradn_[i] = interiorGradnp_[i] * 48.0 / 25.0 - interiorGradnm_[i] * 36.0 / 25.0 + interiorGradnmm_[i] * 16.0 / 25.0 - interiorGradnmmm_[i] * 3.0 / 25.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVeln_[i]  = interiorVelnp_[i]  * 48.0 / 25.0 - interiorVelnm_[i]  * 36.0 / 25.0 + interiorVelnmm_[i]  * 16.0 / 25.0 - interiorVelnmmm_[i]  * 3.0 / 25.0;
  }
  else
  {
    interiorGradn_ = interiorGradnp_;
    interiorVeln_ = interiorVelnp_;
  }
  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(nsd_*ndofs_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradn_,0.0);
  tempVec1.Multiply('N','N',-1.0/rho,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2;
  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVeln_,0.0);
  tempVec2.Multiply('N','N',-1.0,localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',-1.0,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2.Multiply('N','T',-1.0/rho,tempMat1,localSolver_.Bmat,0.0);
  tempMat2 += localSolver_.Dmat;

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',-1.0/rho,localSolver_.Bmat,interiorVelnp_,1.0);

  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  // sort this back to the interiorvalues vector
  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  if(dyna_ == INPAR::ACOU::acou_bdf2)
  {
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradnp_[i] = 4.0 / 3.0 * interiorGradnp_[i] - 1.0 / 3.0 * tempGradnp[i];
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVelnp_[i]  = 4.0 / 3.0 * interiorVelnp_[i]  - 1.0 / 3.0 * tempVelnp[i];
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf3)
  {
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradnp_[i] = 18.0 / 11.0 * interiorGradnp_[i] - 9.0 / 11.0 * tempGradnp[i] + 2.0 / 11.0 * interiorGradnm_[i];
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVelnp_[i]  = 18.0 / 11.0 * interiorVelnp_[i]  - 9.0 / 11.0 * tempVelnp[i]  + 2.0 / 11.0 * interiorVelnm_[i];
  }
  else if(dyna_ == INPAR::ACOU::acou_bdf4)
  {
    for(unsigned int i=0; i<ndofs_*nsd_; ++i)
      interiorGradnp_[i] = interiorGradnp_[i] * 48.0 / 25.0 - tempGradnp[i] * 36.0 / 25.0 + interiorGradnm_[i] * 16.0 / 25.0 - interiorGradnmm_[i] * 3.0 / 25.0;
    for(unsigned int i=0; i<ndofs_; ++i)
      interiorVelnp_[i]  = interiorVelnp_[i]  * 48.0 / 25.0 - tempVelnp[i]  * 36.0 / 25.0 + interiorVelnm_[i]  * 16.0 / 25.0 - interiorVelnmm_[i]  * 3.0 / 25.0;
  }

  tempVec1.Shape(ndofs_,1);
  tempVec1.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec1.Multiply('N','N',1.0,localSolver_.Bmat,interiorGradnp_,1.0);

  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0); // = W = ( D + 1/rho B A^{-1} B^T )^{-1} ( M V^{n-1} - 1/rho B Q^{n-1} )

  elevec.Multiply('T','N',-1.0,localSolver_.Emat,tempVec2,0.0);

  tempVec1.Shape(ndofs_*nsd_,1);
  tempVec1.Multiply('T','N',-1.0/rho,localSolver_.Bmat,tempVec2,0.0);

  tempVec2.Shape(ndofs_*nsd_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);
  tempVec2 += interiorGradnp_;

  elevec.Multiply('T','N',1.0,localSolver_.Cmat,tempVec2,1.0);

  return;
} // UpdateInteriorVariablesAndComputeResidualAdjoint

/*----------------------------------------------------------------------*
 * UpdateInteriorVariablesAndComputeResidualTrapAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::
UpdateInteriorVariablesAndComputeResidualTrapAdjoint(DRT::Discretization &     discretization,
                                                     DRT::ELEMENTS::Acou &                ele,
                                                     Epetra_SerialDenseVector          & elevec,
                                                     const Teuchos::RCP<MAT::Material> &mat)
{
  interiorGradn_ = interiorGradnp_;
  interiorVeln_ = interiorVelnp_;
  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal_[i];
  Epetra_SerialDenseVector traceVal_SDV_m(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV_m(i) = traceValm_[i];

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  double theta = 0.66;

  // *****************************************************
  // update interior variables first
  // *****************************************************

  Epetra_SerialDenseVector tempVec1(ndofs_*nsd_);
  tempVec1.Multiply('N','N',1.0,localSolver_.Amat,interiorGradnp_,0.0);
  tempVec1.Multiply('T','N',-(1.0-theta)/rho,localSolver_.Bmat,interiorVelnp_,1.0);
  tempVec1.Multiply('N','N',-(1.0-theta)/rho,localSolver_.Cmat,traceVal_SDV_m,1.0);
  tempVec1.Multiply('N','N',-theta/rho,localSolver_.Cmat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec2.Multiply('N','N',(1.0-theta),localSolver_.Bmat,interiorGradnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.DmMmat,interiorVelnp_,1.0);
  tempVec2.Multiply('N','N',-(1.0-theta),localSolver_.Emat,traceVal_SDV_m,1.0);
  tempVec2.Multiply('N','N',-(theta),localSolver_.Emat,traceVal_SDV,1.0);

  // now, we have to do the Schur complement thing, just as in "CondenseLocalPart" but without C and E
  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,nsd_*ndofs_);
  tempMat1.Multiply('N','N',-theta,localSolver_.Bmat,localSolver_.invAmat,0.0);
  tempVec2.Multiply('N','N',-1.0,tempMat1,tempVec1,1.0);


  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 += localSolver_.DmMmat;
  tempMat2.Scale(theta);
  tempMat2.Multiply('N','T',-theta/rho,tempMat1,localSolver_.Bmat,1.0);
  tempMat2 += localSolver_.Mmat;


  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseDmat;
    LINALG::Matrix<ndofs_,ndofs_> D(tempMat2,true);
    inverseDmat.SetMatrix(D);
    inverseDmat.Invert();
  }

  interiorVelnp_.Multiply('N','N',1.0,tempMat2,tempVec2,0.0);

  tempVec1.Multiply('T','N',-theta/rho,localSolver_.Bmat,interiorVelnp_,1.0);

  interiorGradnp_.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  // sort this back to the interior values vector
  for(unsigned int i=0; i<interiorValnp_.size(); ++i)
  {
    if ((i+1)%(nsd_+1) == 0)
    {
      interiorValnp_[i] = interiorVelnp_((i+1)/(nsd_+1)-1);
    }
    else
    {
      int xyz = i % (nsd_+1); // 0 for x, 1 for y and 2 for z (for 3D)
      interiorValnp_[i] = interiorGradnp_(xyz*ndofs_+i/(nsd_+1));
    }
  }

  // tell this change in the interior variables the discretization
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(1,"intvel");
  std::vector<int> localDofs = discretization.Dof(1, &ele);
  Epetra_Vector& secondary = const_cast<Epetra_Vector&>(*matrix_state);

  const Epetra_Map* intdofcolmap = discretization.DofColMap(1);
  for (unsigned int i=0; i<localDofs.size(); ++i)
  {
    const int lid = intdofcolmap->LID(localDofs[i]);
    secondary[lid] = interiorValnp_[i];
  }

  // *****************************************************
  // compute residual second (reuse intermediate matrices)
  // *****************************************************

  tempVec1.Shape(ndofs_,1);
  tempVec1.Multiply('N','N',1.0,localSolver_.Mmat,interiorVelnp_,0.0);
  tempVec1.Multiply('N','N',(1.0-theta),localSolver_.Bmat,interiorGradnp_,1.0);
  tempVec1.Multiply('N','N',-(1.0-theta),localSolver_.DmMmat,interiorVelnp_,1.0);
  tempVec1.Multiply('N','N',-(1.0-theta),localSolver_.Emat,traceVal_SDV,1.0);

  Epetra_SerialDenseVector f(ndofs_*nsd_);
  f.Multiply('N','N',1.0,localSolver_.Amat,interiorGradnp_,0.0);
  f.Multiply('T','N',-(1.0-theta)/rho,localSolver_.Bmat,interiorVelnp_,1.0);
  f.Multiply('N','N',-(1.0-theta)/rho,localSolver_.Cmat,traceVal_SDV,1.0);

  tempVec1.Multiply('N','N',-1.0,tempMat1,f,1.0);

  tempVec2.Shape(ndofs_,1);
  tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0);

  elevec.Multiply('T','N',-theta,localSolver_.Emat,tempVec2,0.0);

  tempVec1.Shape(ndofs_*nsd_,1);
  tempVec1.Multiply('T','N',-theta/rho,localSolver_.Bmat,tempVec2,0.0);
  tempVec1 += f;
  tempVec2.Shape(ndofs_*nsd_,1);
  tempVec2.Multiply('N','N',1.0,localSolver_.invAmat,tempVec1,0.0);

  elevec.Multiply('T','N',theta,localSolver_.Cmat,tempVec2,1.0);

  elevec.Multiply('T','N',1.0-theta,localSolver_.Cmat,interiorGradnp_,1.0);
  elevec.Multiply('T','N',-(1.0-theta),localSolver_.Emat,interiorVelnp_,1.0);
  elevec.Multiply('N','N',-(1.0-theta),localSolver_.Gmat,traceVal_SDV,1.0);

  return;
}

template <DRT::Element::DiscretizationType distype>
double DRT::ELEMENTS::AcouEleCalc<distype>::
CalculateError(DRT::ELEMENTS::Acou & ele,Epetra_SerialDenseMatrix & h, Epetra_SerialDenseVector & rhs, Epetra_SerialDenseVector & p)
{
  DRT::UTILS::LagrangeBasis<nsd_> postpoly(DRT::ELEMENTS::Acou::degree+1);
  LINALG::Matrix<nsd_,1> xsi;
  const unsigned int ndofspost = DRT::UTILS::FixedPower<DRT::ELEMENTS::Acou::degree+2,nsd_>::value;

  LINALG::Matrix<nsd_,ndofspost> derivs;
  LINALG::Matrix<1,ndofs_> myvalues;
  LINALG::Matrix<nsd_,nen_> deriv;
  LINALG::Matrix<nsd_,nsd_> xjm, xji;

  LINALG::Matrix<1,ndofspost> values;

  Teuchos::RCP<DRT::UTILS::GaussPoints> postquad;

  postquad = DRT::UTILS::GaussPointCache::Instance().Create(distype, (DRT::ELEMENTS::Acou::degree+1)*2);

  LINALG::Matrix<nsd_,nen_> xyze;
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(&ele,xyze);

  for (int q=0; q<postquad->NumPoints(); ++q)
  {
    const double* gpcoord = postquad->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi,values);
    postpoly.Evaluate_deriv1(xsi,derivs);

    DRT::UTILS::shape_function_deriv1<distype>(xsi,deriv);
    xjm.MultiplyNT(deriv,xyze);
    const double jfac = xji.Invert(xjm) * postquad->Weight(q);

    // transform shape functions derivatives
    for (unsigned int i=0; i<ndofspost; ++i)
    {
      double res[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
      {
        res[d] = xji(d,0) * derivs(0,i);
        for (unsigned int e=1; e<nsd_; ++e)
          res[d] += xji(d,e) * derivs(e,i);
      }
      for (unsigned int d=0; d<nsd_; ++d)
        derivs(d,i) = res[d];
    }

    shapes_.polySpace_.Evaluate(xsi,myvalues);

    for (unsigned int j=0; j<ndofspost; ++j)
      h(0,j) += values(j) * jfac;
    for (unsigned int j=0; j<ndofs_; ++j)
      rhs(0) += myvalues(j) * jfac * interiorVelnp_[j];

    for (unsigned int i=1; i<ndofspost; ++i)
    {
      for (unsigned int j=0; j<ndofspost; ++j)
      {
        double t = 0;
        for (unsigned int d=0; d<nsd_; ++d)
          t += derivs(d,i) * derivs(d,j);
        h(i,j) += t * jfac;
      }
    }
    double ugrad[nsd_];
    for (unsigned int d=0; d<nsd_; ++d)
      ugrad[d]=0.0;
    for (unsigned int d=0; d<nsd_; ++d)
      for (unsigned int j=0; j<ndofs_; ++j)
        ugrad[d] += myvalues(j) * p(j+d*ndofs_);
    for (unsigned int i=1; i<ndofspost; ++i)
      for (unsigned int d=0; d<nsd_; ++d)
        rhs(i) += ugrad[d] * derivs(d,i) * jfac;

  } // for (int q=0; q<postquad->NumPoints(); ++q)


  LINALG::FixedSizeSerialDenseSolver<ndofspost,ndofspost,1> inverseH;
  LINALG::Matrix<ndofspost,ndofspost> Hmat(h,true);
  LINALG::Matrix<ndofspost,1> Rvec(rhs,true);
  inverseH.SetMatrix(Hmat);
  inverseH.SetVectors(Rvec,Rvec);
  inverseH.Solve();

  double err_p = 0.0;
  double numerical_post = 0.0;
  double numerical = 0.0;
  for (unsigned int q=0; q<ndofs_; ++q)
  {
    const double* gpcoord = shapes_.quadrature_->Point(q);
    for (unsigned int idim=0;idim<nsd_;idim++)
      xsi(idim) = gpcoord[idim];

    postpoly.Evaluate(xsi,values);
    for (unsigned int i=0; i<ndofspost; ++i)
      numerical_post += values(i) * rhs(i);
    for (unsigned int i=0; i<ndofs_; ++i)
      numerical += shapes_.shfunct(i,q) * interiorVelnp_(i);

    err_p +=  (numerical_post - numerical) * (numerical_post - numerical) * shapes_.jfac(q);
  } // for (int q=0; q<postquad->NumPoints(); ++q)

  return err_p;
} // FillMatrices


/*----------------------------------------------------------------------*
 * ComputeAbsorbingBC
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeAbsorbingBC(DRT::ELEMENTS::Acou*        ele,
                   Teuchos::ParameterList&     params,
                   Teuchos::RCP<MAT::Material> & mat,
                   int                         face,
                   Epetra_SerialDenseMatrix    &elemat,
                   Epetra_SerialDenseVector    &elevec)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeAbsorbingBC");

  bool resonly = params.get<bool>("resonly");
  bool adjoint = params.get<bool>("adjoint");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();

  if(!resonly)
  {

    // loop over number of shape functions
    for (unsigned int p=0; p<nfdofs_; ++p)
    {
      // loop over number of shape functions
      for (unsigned int q=0; q<=p; ++q)
      {
        double tempG = 0.0;
        for (unsigned int i=0; i<nfdofs_; ++i)
        {
          tempG += shapes_.jfacF(i) * shapes_.shfunctF(p,i) * shapes_.shfunctF(q,i);
        }

        elemat(face*nfdofs_+p,face*nfdofs_+q) =
        elemat(face*nfdofs_+q,face*nfdofs_+p) += tempG / rho / c;

      } // for (unsigned int q=0; q<nfdofs_; ++q)
    } // for (unsigned int p=0; p<nfdofs_; ++p)
  } // if(!resonly)

  // for the adjoint problem, we have a right hand side for this boundary condition
  //if(adjoint && resonly)
  if(adjoint)
  {
    // get the values for the source term!
    Teuchos::RCP<Epetra_MultiVector> adjointrhs = params.get<Teuchos::RCP<Epetra_MultiVector> >("adjointrhs");
    int step = params.get<int>("step");
    int stepmax = adjointrhs->NumVectors();

    if(step>0)
    {
      // adjointrhs is a nodebased vector. Which nodes do belong to this element and which values are associated to those?
      const int * fnodeIds = ele->Faces()[face]->NodeIds();
      int numfnode = ele->Faces()[face]->NumNode();
      double fnodexyz[numfnode][nsd_];
      double values[numfnode];

      for(int i=0; i<numfnode; ++i)
      {
        int localnodeid = adjointrhs->Map().LID(fnodeIds[i]);
        values[i] = adjointrhs->operator ()(stepmax-step)->operator [](localnodeid); // in inverse order -> we're integrating backwards in time
        for(unsigned int d=0; d<nsd_; ++d)
          fnodexyz[i][d] = shapes_.xyzeF(i,d);
//          fnodexyz[i][d] = ele->Faces()[face]->Nodes()[i]->X()[d];
      }// for(int i=0; i<numfnode; ++i)

      // now get the nodevalues to the dof values just as in ProjectOpticalField function

      LINALG::Matrix<nfdofs_,nfdofs_> mass(true);
      LINALG::Matrix<nfdofs_,1> trVec(true);

      for(unsigned int q=0; q<nfdofs_; ++q)
      {
        const double fac = shapes_.jfacF(q);
        double xyz[nsd_];
        for (unsigned int d=0; d<nsd_; ++d)
          xyz[d] = shapes_.xyzFreal(d,q);
        double val = 0.0;

        EvaluateFaceAdjoint(fnodexyz,values,numfnode,xyz,val);

        for (unsigned int i=0; i<nfdofs_; ++i)
        {
          // mass matrix
          for (unsigned int j=0; j<nfdofs_; ++j)
            mass(i,j) += shapes_.shfunctF(i,q) * shapes_.shfunctF(j,q) * fac;
          trVec(i,0) += shapes_.shfunctF(i,q) * val * fac;
        }
      } // for(unsigned int q=0; q<nfdofs_; ++q)

      LINALG::FixedSizeSerialDenseSolver<nfdofs_,nfdofs_,1> inverseMass;
      inverseMass.SetMatrix(mass);
      inverseMass.SetVectors(trVec,trVec);
      inverseMass.Solve();

      for (unsigned int i=0; i<nfdofs_; ++i) // integration points
      {
        const double fac = shapes_.jfacF(i);
        double p = 0.0;
        for(unsigned int q=0; q<nfdofs_; ++q) // shape functions
        {
          p += shapes_.shfunctF(q,i) * trVec(q,0);
        }
        for(unsigned int r=0; r<nfdofs_; ++r)
        {
          elevec(face*nfdofs_+r) +=  p * fac * shapes_.shfunctF(r,i);
        }

      } // for all integration points
    } // if(step>0)
  } // if(adjoint)
  return;
}

/*----------------------------------------------------------------------*
 * ComputeObjfIntegral
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeObjfIntegral(DRT::ELEMENTS::Acou*        ele,
                    Teuchos::ParameterList&     params,
                    Teuchos::RCP<MAT::Material> & mat,
                    int                         face,
                    double                      dt)
{
  double integr_objf = params.get<double>("integr_objf");

  Teuchos::RCP<Epetra_MultiVector> rhsvec = params.get<Teuchos::RCP<Epetra_MultiVector> >("rhsvec");
  const int * fnodeIds = ele->Faces()[face]->NodeIds();
  int numfnode = ele->Faces()[face]->NumNode();
  double fnodexyz[numfnode][nsd_];
  double values[numfnode];

  for(int t=0; t<rhsvec->NumVectors(); ++t)
  {
    for(int i=0; i<numfnode; ++i)
    {
      int localnodeid = rhsvec->Map().LID(fnodeIds[i]);
      double temp = rhsvec->operator ()(t)->operator [](localnodeid);
      values[i] =  temp * temp;
      for(unsigned int d=0; d<nsd_; ++d)
        fnodexyz[i][d] = ele->Faces()[face]->Nodes()[i]->X()[d];

    }// for(int i=0; i<numfnode; ++i)

    LINALG::Matrix<nfdofs_,nfdofs_> mass(true);
    LINALG::Matrix<nfdofs_,1> trVec(true);

    for(unsigned int q=0; q<nfdofs_; ++q)
    {
      const double fac = shapes_.jfacF(q);
      double xyz[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        xyz[d] = shapes_.xyzFreal(d,q);
      double val = 0.0;

      EvaluateFaceAdjoint(fnodexyz,values,numfnode,xyz,val);

      for (unsigned int i=0; i<nfdofs_; ++i)
      {
        // mass matrix
        for (unsigned int j=0; j<nfdofs_; ++j)
          mass(i,j) += shapes_.shfunctF(i,q) * shapes_.shfunctF(j,q) * fac;
        trVec(i,0) += shapes_.shfunctF(i,q) * val * fac;
      }
    } // for(unsigned int q=0; q<nfdofs_; ++q)

    LINALG::FixedSizeSerialDenseSolver<nfdofs_,nfdofs_,1> inverseMass;
    inverseMass.SetMatrix(mass);
    inverseMass.SetVectors(trVec,trVec);
    inverseMass.Solve();

    // trVec now holds the dofbased values for p-p_m
    // integrate this
    double integralvalue = 0.0;
    for(unsigned int i=0; i<nfdofs_; ++i) // sum integration points
    {
      const double fac = shapes_.jfacF(i);
      for (unsigned int j=0; j<nfdofs_; ++j) // sum shape functions
        integralvalue += fac * shapes_.shfunctF(j,i) * trVec(j,0);
    }
    integr_objf += integralvalue * dt;
  } // for(int t=0; t<rhsvec->NumVectors(); ++t)

  params.set<double>("integr_objf",integr_objf);

  return;
} // ComputeObjfIntegral

template<DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::EvaluateFaceAdjoint(
                    double fnodexyz[][nsd_],
                    double values[],
                    int numfnode,
                    const double (&xyz)[nsd_],
                    double &val) const
{
  const DRT::Element::DiscretizationType facedis = DRT::UTILS::DisTypeToFaceShapeType<distype>::shape;
  if (facedis == DRT::Element::line2)
  {
    if(numfnode!=2) dserror("NEIN");
    val = values[1] * sqrt( (fnodexyz[0][0]-xyz[0])*(fnodexyz[0][0]-xyz[0]) + (fnodexyz[0][1]-xyz[1])*(fnodexyz[0][1]-xyz[1]) )
        + values[0] * sqrt( (fnodexyz[1][0]-xyz[0])*(fnodexyz[1][0]-xyz[0]) + (fnodexyz[1][1]-xyz[1])*(fnodexyz[1][1]-xyz[1]) );
    double dist = sqrt( (fnodexyz[0][0]-fnodexyz[1][0])*(fnodexyz[0][0]-fnodexyz[1][0]) + (fnodexyz[0][1]-fnodexyz[1][1])*(fnodexyz[0][1]-fnodexyz[1][1]) );
    val /= dist;
  }
  else
    dserror("not yet implemented"); // TODO

  return;
}


/*----------------------------------------------------------------------*
 * ComputeInteriorMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeInteriorMatrices(const Teuchos::RCP<MAT::Material> &mat, double dt)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeInteriorMatrices");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();

  zeroMatrix(Amat);
  zeroMatrix(Mmat);
  zeroMatrix(Dmat);
  zeroMatrix(Bmat);
  zeroMatrix(DmMmat);

  Epetra_SerialDenseMatrix  massPart(ndofs_,ndofs_);
  Epetra_SerialDenseMatrix  gradPart(ndofs_*nsd_,ndofs_);

  // loop quadrature points
  for (unsigned int q=0; q<ndofs_; ++q)
  {
    const double sqrtfac = std::sqrt(shapes_.jfac(q));
    // loop shape functions
    for (unsigned int i=0; i<ndofs_; ++i)
    {
      const double valf = shapes_.shfunct(i,q) * sqrtfac;
      massPart(i,q) = valf;
      for (unsigned int d=0; d<nsd_; ++d)
      {
        const double vald = shapes_.shderxy(i*nsd_+d,q) * sqrtfac;
        gradPart (d*ndofs_+i,q) = vald;
      }
    }
  }

  // multiply matrices to perform summation over quadrature points
  for (unsigned int i=0; i<ndofs_; ++i)
    for (unsigned int j=0; j<ndofs_; ++j)
    {
      double sum = 0;
      double sums[nsd_];
      for (unsigned int d=0; d<nsd_; ++d)
        sums[d] = 0;
      for (unsigned int k=0; k<ndofs_; ++k)
      {
        sum += massPart(i,k) * massPart(j,k);
        for (unsigned int d=0; d<nsd_; ++d)
          sums[d] += gradPart(d*ndofs_+i,k) * massPart(j,k);
      }
      Dmat(i,j) = Mmat(i,j) = sum;
      for (unsigned int d=0; d<nsd_; ++d)
      {
        Amat(d*ndofs_+i,d*ndofs_+j) = sum;
        Bmat(j,d*ndofs_+i) = -sums[d];
      }
    }

  Amat.Scale(1.0/dt);
  Mmat.Scale(1.0/dt/rho/c/c);
  Dmat.Scale(1.0/dt/rho/c/c);

  invAmat = Amat;

  // we will need the inverse of A
  {
    LINALG::FixedSizeSerialDenseSolver<nsd_*ndofs_,nsd_*ndofs_> inverseAmat;
    LINALG::Matrix<nsd_*ndofs_,nsd_*ndofs_> A(invAmat,true);
    inverseAmat.SetMatrix(A);
    inverseAmat.Invert();
  }

  return;
}


/*----------------------------------------------------------------------*
 * ComputeResidual
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeResidual(Epetra_SerialDenseVector          & elevec,
                const Teuchos::RCP<MAT::Material> &mat,
                Epetra_SerialDenseVector          & interiorGradnp,
                Epetra_SerialDenseVector          & interiorVelnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeResidual");

  /*
                              -1
                   +---------+    +---------+
                   |       T |    |    n-1  |
    n              | A   -B  |    | A Q     |
   R  = - [ C  E ] |         |    |    n-1  |
                   | B    D  |    | M V     |
                   +---------+    +---------+

  */

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector tempVec1(ndofs_);
  tempVec1.Multiply('N','N',1.0,Mmat,interiorVelnp,0.0);
  tempVec1.Multiply('N','N',-1.0/rho,Bmat,interiorGradnp,1.0); // = M V^{n-1} - 1/rho B Q^{n-1}

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,ndofs_*nsd_);
  tempMat1.Multiply('N','N',1.0/rho,Bmat,invAmat,0.0); // = 1/rho B A^{-1}
  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 = Dmat;
  tempMat2.Multiply('N','T',1.0,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
    LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
    inverseinW.SetMatrix(inv);
    inverseinW.Invert();
  }
  // tempMat2 = ( D + 1/rho B A^{-1} B^T )^{-1}

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0); // = W = ( D + 1/rho B A^{-1} B^T )^{-1} ( M V^{n-1} - 1/rho B Q^{n-1} )

  elevec.Multiply('T','N',-1.0,Emat,tempVec2,0.0);

  tempVec1.Shape(ndofs_*nsd_,1);
  tempVec1.Multiply('T','N',1.0,Bmat,tempVec2,0.0);

  tempVec2.Shape(ndofs_*nsd_,1);
  tempVec2.Multiply('N','N',1.0,invAmat,tempVec1,0.0);
  tempVec2 += interiorGradnp;

  tempVec1.Multiply('N','N',1.0,Amat,interiorGradnp,0.0);

  elevec.Multiply('T','N',-1.0/rho,Cmat,tempVec2,1.0);

  return;
} // ComputeResidual

/*----------------------------------------------------------------------*
 * ComputeResidualTrap
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeResidualTrap(Epetra_SerialDenseVector          & elevec,
                const Teuchos::RCP<MAT::Material> &mat,
                Epetra_SerialDenseVector          & interiorGradnp,
                Epetra_SerialDenseVector          & interiorVelnp,
                std::vector<double>               traceVal,
                double                            dt)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeResidualTrap");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal[i];

  double theta = 0.66;

  {
    Epetra_SerialDenseVector tempVec1(ndofs_);
    tempVec1.Multiply('N','N',1.0,Mmat,interiorVelnp,0.0);
    tempVec1.Multiply('N','N',-(1.0-theta)/rho,Bmat,interiorGradnp,1.0);
    tempVec1.Multiply('N','N',-(1.0-theta),DmMmat,interiorVelnp,1.0);
    tempVec1.Multiply('N','N',-(1.0-theta),Emat,traceVal_SDV,1.0);

    Epetra_SerialDenseVector f(ndofs_*nsd_);
    f.Multiply('N','N',1.0,Amat,interiorGradnp,0.0);
    f.Multiply('T','N',1.0-theta,Bmat,interiorVelnp,1.0);
    f.Multiply('N','N',1.0-theta,Cmat,traceVal_SDV,1.0);

    Epetra_SerialDenseMatrix tempMat1;
    tempMat1.Shape(ndofs_,ndofs_*nsd_);
    tempMat1.Multiply('N','N',theta/rho,Bmat,invAmat,0.0); // = 1/rho B A^{-1}

    tempVec1.Multiply('N','N',-1.0,tempMat1,f,1.0); // = M V^{n-1} - 1/rho B Q^{n-1}

    Epetra_SerialDenseMatrix tempMat2;
    tempMat2.Shape(ndofs_,ndofs_);
    tempMat2 = DmMmat;
    tempMat2.Scale(theta);
    tempMat2 += Mmat;
    tempMat2.Multiply('N','T',theta,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
    {
      LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
      LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
      inverseinW.SetMatrix(inv);
      inverseinW.Invert();
    }

    Epetra_SerialDenseVector tempVec2(ndofs_);
    tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0); // = W = ( D + 1/rho B A^{-1} B^T )^{-1} ( M V^{n-1} - 1/rho B Q^{n-1} )

    //Epetra_SerialDenseVector testelevec(nfaces_*nfdofs_);

    elevec.Multiply('T','N',-theta,Emat,tempVec2,0.0);

    tempVec1.Shape(ndofs_*nsd_,1);
    tempVec1.Multiply('T','N',theta,Bmat,tempVec2,0.0);
    tempVec1 += f;
    tempVec2.Shape(ndofs_*nsd_,1);
    tempVec2.Multiply('N','N',1.0,invAmat,tempVec1,0.0);

    elevec.Multiply('T','N',-theta/rho,Cmat,tempVec2,1.0);

    elevec.Multiply('T','N',-(1.0-theta)/rho,Cmat,interiorGradnp,1.0);
    elevec.Multiply('T','N',-(1.0-theta),Emat,interiorVelnp,1.0);
    elevec.Multiply('N','N',-(1.0-theta),Gmat,traceVal_SDV,1.0);
  }

  return;
} // ComputeResidualTrap

/*----------------------------------------------------------------------*
 * ComputeResidualAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeResidualAdjoint(Epetra_SerialDenseVector          & elevec,
                       const Teuchos::RCP<MAT::Material> &mat,
                       Epetra_SerialDenseVector          & interiorGradnp,
                       Epetra_SerialDenseVector          & interiorVelnp)
{
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector tempVec1(ndofs_);
  tempVec1.Multiply('N','N',1.0,Mmat,interiorVelnp,0.0);
  tempVec1.Multiply('N','N',1.0,Bmat,interiorGradnp,1.0); // = M V^{n-1} - 1/rho B Q^{n-1}

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,ndofs_*nsd_);
  tempMat1.Multiply('N','N',-1.0,Bmat,invAmat,0.0); // = 1/rho B A^{-1}
  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 = Dmat;
  tempMat2.Multiply('N','T',-1.0/rho,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
    LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
    inverseinW.SetMatrix(inv);
    inverseinW.Invert();
  }

  Epetra_SerialDenseVector tempVec2(ndofs_);
  tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0); // = W = ( D + 1/rho B A^{-1} B^T )^{-1} ( M V^{n-1} - 1/rho B Q^{n-1} )

  elevec.Multiply('T','N',-1.0,Emat,tempVec2,0.0);

  tempVec1.Shape(ndofs_*nsd_,1);
  tempVec1.Multiply('T','N',-1.0/rho,Bmat,tempVec2,0.0);

  tempVec2.Shape(ndofs_*nsd_,1);
  tempVec2.Multiply('N','N',1.0,invAmat,tempVec1,0.0);
  tempVec2 += interiorGradnp;

  tempVec1.Multiply('N','N',1.0,Amat,interiorGradnp,0.0);

  elevec.Multiply('T','N',1.0,Cmat,tempVec2,1.0);

  return;
} // ComputeResidualAdjoint

/*----------------------------------------------------------------------*
 * ComputeResidualTrapAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeResidualTrapAdjoint(Epetra_SerialDenseVector          & elevec,
                           const Teuchos::RCP<MAT::Material> &mat,
                           Epetra_SerialDenseVector          & interiorGradnp,
                           Epetra_SerialDenseVector          & interiorVelnp,
                           std::vector<double>               traceVal,
                           double                            dt)
{

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  Epetra_SerialDenseVector traceVal_SDV(nfaces_*nfdofs_);
  for(unsigned i=0; i<nfaces_*nfdofs_; ++i)
    traceVal_SDV(i) = traceVal[i];

  double theta = 0.66;

  {
    Epetra_SerialDenseVector tempVec1(ndofs_);
    tempVec1.Multiply('N','N',1.0,Mmat,interiorVelnp,0.0);
    tempVec1.Multiply('N','N',(1.0-theta),Bmat,interiorGradnp,1.0);
    tempVec1.Multiply('N','N',-(1.0-theta),DmMmat,interiorVelnp,1.0);
    tempVec1.Multiply('N','N',-(1.0-theta),Emat,traceVal_SDV,1.0);

    Epetra_SerialDenseVector f(ndofs_*nsd_);
    f.Multiply('N','N',1.0,Amat,interiorGradnp,0.0);
    f.Multiply('T','N',-(1.0-theta)/rho,Bmat,interiorVelnp,1.0);
    f.Multiply('N','N',-(1.0-theta)/rho,Cmat,traceVal_SDV,1.0);

    Epetra_SerialDenseMatrix tempMat1;
    tempMat1.Shape(ndofs_,ndofs_*nsd_);
    tempMat1.Multiply('N','N',-theta,Bmat,invAmat,0.0); // = 1/rho B A^{-1}

    tempVec1.Multiply('N','N',-1.0,tempMat1,f,1.0); // = M V^{n-1} - 1/rho B Q^{n-1}

    Epetra_SerialDenseMatrix tempMat2;
    tempMat2.Shape(ndofs_,ndofs_);
    tempMat2 = DmMmat;
    tempMat2.Scale(theta);
    tempMat2 += Mmat;
    tempMat2.Multiply('N','T',-theta/rho,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
    {
      LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
      LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
      inverseinW.SetMatrix(inv);
      inverseinW.Invert();
    }

    Epetra_SerialDenseVector tempVec2(ndofs_);
    tempVec2.Multiply('N','N',1.0,tempMat2,tempVec1,0.0); // = W = ( D + 1/rho B A^{-1} B^T )^{-1} ( M V^{n-1} - 1/rho B Q^{n-1} )

    elevec.Multiply('T','N',-theta,Emat,tempVec2,0.0);

    tempVec1.Shape(ndofs_*nsd_,1);
    tempVec1.Multiply('T','N',-theta/rho,Bmat,tempVec2,0.0);
    tempVec1 += f;
    tempVec2.Shape(ndofs_*nsd_,1);
    tempVec2.Multiply('N','N',1.0,invAmat,tempVec1,0.0);

    elevec.Multiply('T','N',theta,Cmat,tempVec2,1.0);

    elevec.Multiply('T','N',(1.0-theta),Cmat,interiorGradnp,1.0);
    elevec.Multiply('T','N',-(1.0-theta),Emat,interiorVelnp,1.0);
    elevec.Multiply('N','N',-(1.0-theta),Gmat,traceVal_SDV,1.0);
  }

  return;
} // ComputeResidualTrapAdjoint

/*----------------------------------------------------------------------*
 * ComputeFaceMatrices
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
ComputeFaceMatrices(const int                          face,
                    const Teuchos::RCP<MAT::Material> &mat,
                    double dt)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::ComputeFaceMatrices");

  // Compute the matrices C, E and G
  // Here, we don't consider material properties! Don't forget this during condensation

  // tau = rho / t_c where t_c is the characteristic time scale or tau = rho * omega_c with the characteristic frequency omega_c
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();
  double c = actmat->SpeedofSound();
  double tau = 1.0 / c / c / rho / dt;// / 1.0e6;

  // loop over number of shape functions
  for (unsigned int p=0; p<nfdofs_; ++p)
  {
    // loop over number of shape functions
    for (unsigned int q=0; q<ndofs_; ++q)
    {
      // C and E

      // numerical integration: sum over quadrature points
      double tempE = 0.0;
      for (unsigned int i=0; i<nfdofs_; ++i)
      {
        double temp = shapes_.jfacF(i) * shapes_.shfunctF(p,i) * shapes_.shfunctI[face](q,i);
        tempE += temp;
        for(unsigned int j=0; j<nsd_; ++j)
        {
          double temp_d = temp*shapes_.normals(j,i);
          Cmat(j*ndofs_+q,face*nfdofs_+p) += temp_d;
        }
      }
      Emat(q,face*nfdofs_+p) -= tau * tempE;

    } // for (unsigned int q=0; q<ndofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // G

  // loop over number of shape functions
  for (unsigned int p=0; p<nfdofs_; ++p)
  {
    // loop over number of shape functions
    for (unsigned int q=0; q<=p; ++q)
    {
      double tempG = 0.0;
      for (unsigned int i=0; i<nfdofs_; ++i)
      {
        tempG += shapes_.jfacF(i) * shapes_.shfunctF(p,i) * shapes_.shfunctF(q,i);
      }
      Gmat(face*nfdofs_+p,face*nfdofs_+q) =
      Gmat(face*nfdofs_+q,face*nfdofs_+p) += tau * tempG;

    } // for (unsigned int q=0; q<nfdofs_; ++q)
  } // for (unsigned int p=0; p<nfdofs_; ++p)

  // one term is still missing in D!!
  for (unsigned int p=0; p<ndofs_; ++p)
  {
    for (unsigned int q=0; q<=p; ++q)
    {
      double tempD = 0.0;
      for (unsigned int i=0; i<nfdofs_; ++i)
      {
        tempD += shapes_.jfacF(i) * shapes_.shfunctI[face](p,i) * shapes_.shfunctI[face](q,i);
      }
      Dmat(p,q) =
      Dmat(q,p) += tau * tempD;
      DmMmat(p,q) = DmMmat(q,p) += tau * tempD;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 * CondenseLocalPart
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
CondenseLocalPart(Epetra_SerialDenseMatrix &eleMat,
                  const Teuchos::RCP<MAT::Material>& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::CondenseLocalPart");

  /*
   THE MATRIX
                              -1
                  +---------+    +-----+
                  |       T |    |   T |
                  | A   -B  |    | -C  |
   K = - [ C  E ] |         |    |     | + G  = - C V - E W + G
                  | B    D  |    |  E  |
                  +---------+    +-----+

                 -1  T  -1           -1  T
   W = [ D + B  A   B  ]   [ E + B  A   C  ]

        -1      T    T
   V = A   [ - C  + B  W ]

  */

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();


  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,ndofs_*nsd_);
  tempMat1.Multiply('N','N',1.0/rho,Bmat,invAmat,0.0); // = 1/rho B A^{-1}

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 = Dmat;
  tempMat2.Multiply('N','T',1.0,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
  Epetra_SerialDenseMatrix tempMat3;
  tempMat3.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat3 = Emat;
  tempMat3.Multiply('N','N',1.0,tempMat1,Cmat,1.0); // = E + 1/rho B A^{-1} C^T

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
    LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
    inverseinW.SetMatrix(inv);
    inverseinW.Invert();
  }
  // tempMat2 = ( D + 1/rho B A^{-1} B^T )^{-1}

  eleMat = Gmat; // = G
  tempMat1.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat1.Multiply('N','N',1.0,tempMat2,tempMat3,0.0); // = W
  eleMat.Multiply('T','N',-1.0,Emat,tempMat1,1.0); // = - E W + G

  tempMat2.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat2 = Cmat;
  tempMat2.Multiply('T','N',1.0,Bmat,tempMat1,-1.0); // = - C^T + B^T W

  tempMat3.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat3.Multiply('N','N',1.0,invAmat,tempMat2,0.0); // = V = A^{-1} ( - C^T + B^T W )

  eleMat.Multiply('T','N',-1.0/rho,Cmat,tempMat3,1.0); // = K = - C V - E W + G

  return;
} // CondenseLocalPart

/*----------------------------------------------------------------------*
 * CondenseLocalPartTrap
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
CondenseLocalPartTrap(Epetra_SerialDenseMatrix &eleMat,
                  const Teuchos::RCP<MAT::Material>& mat,
                  double dt)
{
  TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::AcouEleCalc::CondenseLocalPartTrap");

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  double theta = 0.66;

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,ndofs_*nsd_);
  tempMat1.Multiply('N','N',theta/rho,Bmat,invAmat,0.0); // = 1/rho B A^{-1}

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);

  tempMat2 = DmMmat;
  tempMat2.Scale(theta);
  tempMat2 += Mmat;

  tempMat2.Multiply('N','T',theta,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
  Epetra_SerialDenseMatrix tempMat3;
  tempMat3.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat3 = Emat;
  tempMat3.Multiply('N','N',theta,tempMat1,Cmat,theta); // = E + 1/rho B A^{-1} C^T

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
    LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
    inverseinW.SetMatrix(inv);
    inverseinW.Invert();
  }
  // tempMat2 = ( D + 1/rho B A^{-1} B^T )^{-1}

  eleMat = Gmat; // = G
  eleMat.Scale(theta);
  tempMat1.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat1.Multiply('N','N',1.0,tempMat2,tempMat3,0.0); // = W
  eleMat.Multiply('T','N',-theta,Emat,tempMat1,1.0); // = - E W + G

  tempMat2.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat2 = Cmat;
  tempMat2.Multiply('T','N',theta,Bmat,tempMat1,-theta); // = - C^T + B^T W

  tempMat3.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat3.Multiply('N','N',1.0,invAmat,tempMat2,0.0); // = V = A^{-1} ( - C^T + B^T W )

  eleMat.Multiply('T','N',-theta/rho,Cmat,tempMat3,1.0); // = K = - C V - E W + G

  return;
}

/*----------------------------------------------------------------------*
 * CondenseLocalPartAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
CondenseLocalPartAdjoint(Epetra_SerialDenseMatrix &eleMat,
                         const Teuchos::RCP<MAT::Material>& mat)
{

  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();


  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,ndofs_*nsd_);
  tempMat1.Multiply('N','N',-1.0,Bmat,invAmat,0.0); // = 1/rho B A^{-1}

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);
  tempMat2 = Dmat;
  tempMat2.Multiply('N','T',-1.0/rho,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
  Epetra_SerialDenseMatrix tempMat3;
  tempMat3.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat3 = Emat;
  tempMat3.Multiply('N','N',-1.0/rho,tempMat1,Cmat,1.0); // = E + 1/rho B A^{-1} C^T

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
    LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
    inverseinW.SetMatrix(inv);
    inverseinW.Invert();
  }
  // tempMat2 = ( D + 1/rho B A^{-1} B^T )^{-1}

  eleMat = Gmat; // = G
  tempMat1.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat1.Multiply('N','N',1.0,tempMat2,tempMat3,0.0); // = W
  eleMat.Multiply('T','N',-1.0,Emat,tempMat1,1.0); // = - E W + G

  tempMat2.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat2 = Cmat;
  tempMat2.Multiply('T','N',-1.0/rho,Bmat,tempMat1,1.0/rho); // = - C^T + B^T W

  tempMat3.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat3.Multiply('N','N',1.0,invAmat,tempMat2,0.0); // = V = A^{-1} ( - C^T + B^T W )

  eleMat.Multiply('T','N',1.0,Cmat,tempMat3,1.0); // = K = - C V - E W + G

  return;
}

/*----------------------------------------------------------------------*
 * CondenseLocalPartTrapAdjoint
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::AcouEleCalc<distype>::LocalSolver::
CondenseLocalPartTrapAdjoint(Epetra_SerialDenseMatrix &eleMat,
                             const Teuchos::RCP<MAT::Material>& mat)
{
  const MAT::AcousticMat* actmat = static_cast<const MAT::AcousticMat*>(mat.get());
  double rho = actmat->Density();

  double theta = 0.66;

  Epetra_SerialDenseMatrix tempMat1;
  tempMat1.Shape(ndofs_,ndofs_*nsd_);
  tempMat1.Multiply('N','N',-theta,Bmat,invAmat,0.0); // = 1/rho B A^{-1}

  Epetra_SerialDenseMatrix tempMat2;
  tempMat2.Shape(ndofs_,ndofs_);

  tempMat2 = DmMmat;
  tempMat2.Scale(theta);
  tempMat2 += Mmat;

  tempMat2.Multiply('N','T',-theta/rho,tempMat1,Bmat,1.0); // = D + 1/rho B A^{-1} B^T
  Epetra_SerialDenseMatrix tempMat3;
  tempMat3.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat3 = Emat;
  tempMat3.Multiply('N','N',-theta/rho,tempMat1,Cmat,theta); // = E + 1/rho B A^{-1} C^T

  {
    LINALG::FixedSizeSerialDenseSolver<ndofs_,ndofs_> inverseinW;
    LINALG::Matrix<ndofs_,ndofs_> inv(tempMat2,true);
    inverseinW.SetMatrix(inv);
    inverseinW.Invert();
  }
  // tempMat2 = ( D + 1/rho B A^{-1} B^T )^{-1}

  eleMat = Gmat; // = G
  eleMat.Scale(theta);
  tempMat1.Shape(ndofs_,nfaces_*nfdofs_);
  tempMat1.Multiply('N','N',1.0,tempMat2,tempMat3,0.0); // = W
  eleMat.Multiply('T','N',-theta,Emat,tempMat1,1.0); // = - E W + G

  tempMat2.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat2 = Cmat;
  tempMat2.Multiply('T','N',-theta/rho,Bmat,tempMat1,theta/rho); // = - C^T + B^T W

  tempMat3.Shape(ndofs_*nsd_,nfaces_*nfdofs_);
  tempMat3.Multiply('N','N',1.0,invAmat,tempMat2,0.0); // = V = A^{-1} ( - C^T + B^T W )

  eleMat.Multiply('T','N',theta,Cmat,tempMat3,1.0); // = K = - C V - E W + G

  return;
}

// template classes
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::AcouEleCalc<DRT::Element::nurbs27>;
