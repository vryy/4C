/*!-----------------------------------------------------------------------------------------------------------
 \file spring3_evaluate.cpp
 \brief three dimensional spring element (can be connected to beam3eb elements)

<pre>
Maintainer: Dhrubajyoti Mukherjee
            mukherjee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15270
</pre>

 *-----------------------------------------------------------------------------------------------------------*/

#include "spring3.H"
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_beamcontact/beam3contact_utils.H"
#include "../drt_statmech/statmech_manager.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_mat/stvenantkirchhoff.H"
#include "../drt_inpar/inpar_structure.H"
#include "../drt_inpar/inpar_statmech.H"
#include "../drt_lib/standardtypes_cpp.H"

#include "Sacado.hpp"
#include "../headers/FAD_utils.H"
typedef Sacado::Fad::DFad<double> FAD;

/*-----------------------------------------------------------------------------------------------------------*
 |  evaluate the element (public)                                                             mukherjee 04/15|
 *----------------------------------------------------------------------------------------------------------*/
int DRT::ELEMENTS::Spring3::Evaluate(Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  DRT::ELEMENTS::Spring3::ActionType act = Spring3::calc_none;
  // get the action required
  std::string action = params.get<std::string>("action","calc_none");
  if (action == "calc_none") dserror("No action supplied");
  else if (action=="calc_struct_linstiff") act = Spring3::calc_struct_linstiff;
  else if (action=="calc_struct_nlnstiff") act = Spring3::calc_struct_nlnstiff;
  else if (action=="calc_struct_internalforce") act = Spring3::calc_struct_internalforce;
  else if (action=="calc_struct_linstiffmass") act = Spring3::calc_struct_linstiffmass;
  else if (action=="calc_struct_nlnstiffmass") act = Spring3::calc_struct_nlnstiffmass;
  else if (action=="calc_struct_nlnstifflmass") act = Spring3::calc_struct_nlnstifflmass;
  else if (action=="calc_struct_stress") act = Spring3::calc_struct_stress;
  else if (action=="calc_struct_eleload") act = Spring3::calc_struct_eleload;
  else if (action=="calc_struct_fsiload") act = Spring3::calc_struct_fsiload;
  else if (action=="calc_struct_update_istep") act = Spring3::calc_struct_update_istep;
  else if (action=="calc_struct_reset_istep") act = Spring3::calc_struct_reset_istep;
  else if (action=="postprocess_stress") act = Spring3::postprocess_stress;
  else if (action=="calc_struct_ptcstiff") act = Spring3::calc_struct_ptcstiff;
  else if (action=="calc_struct_energy") act = Spring3::calc_struct_energy;
  else
  {
    std::cout<<action<<std::endl;
    dserror("Unknown type of action for Spring3");
  }

  switch(act)
  {
  case Spring3::calc_struct_ptcstiff:
  {
    EvaluatePTC<2,3,3>(params,elemat1);
  }
  break;
  /*in case that only linear stiffness matrix is required b3_nlstiffmass is called with zero displacement and
     residual values*/
  case Spring3::calc_struct_linstiff:
  {
    //only nonlinear case implemented!
    dserror("linear stiffness matrix called, but not implemented");

  }
  break;
  //calculate internal energy
  case Spring3::calc_struct_energy:
  {
    // need current global displacement and get them from discretization
    // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

    t3_energy(params,mydisp,&elevec1);
  }
  break;

  //nonlinear stiffness and mass matrix are calculated even if only nonlinear stiffness matrix is required
  case Spring3::calc_struct_nlnstiffmass:
  case Spring3::calc_struct_nlnstifflmass:
  case Spring3::calc_struct_nlnstiff:
  case Spring3::calc_struct_internalforce:
  {
    // need current global displacement and residual forces and get them from discretization
    // making use of the local-to-global map lm one can extract current displacement and residual values for each degree of freedom
    //
    // get element displacements
    Teuchos::RCP<const Epetra_Vector> disp = discretization.GetState("displacement");
    if (disp==Teuchos::null) dserror("Cannot get state vectors 'displacement'");
    std::vector<double> mydisp(lm.size());
    DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

    // get residual displacements
    Teuchos::RCP<const Epetra_Vector> res  = discretization.GetState("residual displacement");
    if (res==Teuchos::null) dserror("Cannot get state vectors 'residual displacement'");
    std::vector<double> myres(lm.size());
    DRT::UTILS::ExtractMyValues(*res,myres,lm);

    /*first displacement vector is modified for proper element evaluation in case of periodic boundary conditions; in case that
     *no periodic boundary conditions are to be applied the following code line may be ignored or deleted*/
    if( params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
      NodeShift<2,3>(params,mydisp);


    //only if random numbers for Brownian dynamics are passed to element, get element velocities
    std::vector<double> myvel(lm.size());
    if( params.get<  Teuchos::RCP<Epetra_MultiVector> >("RandomNumbers",Teuchos::null) != Teuchos::null)
    {
      Teuchos::RCP<const Epetra_Vector> vel  = discretization.GetState("velocity");
      DRT::UTILS::ExtractMyValues(*vel,myvel,lm);
    }

    // for engineering strains instead of total lagrange use t3_nlnstiffmass2
    if (act == Spring3::calc_struct_nlnstiffmass)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
    else if (act == Spring3::calc_struct_nlnstifflmass)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,&elemat2,&elevec1);
    else if (act == Spring3::calc_struct_nlnstiff)
      t3_nlnstiffmass(params,myvel,mydisp,&elemat1,NULL,&elevec1);
    else if (act == Spring3::calc_struct_internalforce)
      t3_nlnstiffmass(params,myvel,mydisp,NULL,NULL,&elevec1);

    /*at the end of an iteration step the geometric configuration has to be updated: the starting point for the
     * next iteration step is the configuration at the end of the current step */
    Qold_ = Qnew_;
    dispthetaold_= dispthetanew_;

    /*
      //the following code block can be used to check quickly whether the nonlinear stiffness matrix is calculated
      //correctly or not by means of a numerically approximated stiffness matrix
      //The code block will work for all higher order elements.
      if(Id() == 3) //limiting the following tests to certain element numbers
      {
        //assuming the same number of DOF for all nodes
        int numdof = NumDofPerNode(*(Nodes()[0]));
        int nnode  = NumNode();

        //variable to store numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_approx;
        stiff_approx.Shape(numdof*nnode,numdof*nnode);


        //relative error of numerically approximated stiffness matrix
        Epetra_SerialDenseMatrix stiff_relerr;
        stiff_relerr.Shape(numdof*nnode,numdof*nnode);

        //characteristic length for numerical approximation of stiffness
        double h_rel = 1e-9;

        //flag indicating whether approximation leads to significant relative error
        int outputflag = 0;

        //calculating strains in new configuration
        for(int i=0; i<numdof; i++) //for all dof
        {
          for(int k=0; k<nnode; k++)//for all nodes
          {

            Epetra_SerialDenseVector force_aux;
            force_aux.Size(numdof*nnode);

            //create new displacement and velocity vectors in order to store artificially modified displacements
            std::vector<double> vel_aux(myvel);
            std::vector<double> disp_aux(mydisp);

            //modifying displacement artificially (for numerical derivative of internal forces):
            disp_aux[numdof*k + i] += h_rel;
            vel_aux[numdof*k + i]  += h_rel / params.get<double>("delta time",0.01);

            t3_nlnstiffmass(params,vel_aux,disp_aux,NULL,NULL,&force_aux);

            //computing derivative d(fint)/du numerically by finite difference
            for(int u = 0 ; u < numdof*nnode ; u++ )
              stiff_approx(u,k*numdof+i) = ( pow(force_aux[u],2) - pow(elevec1(u),2) )/ (h_rel * (force_aux[u] + elevec1(u) ) );

          } //for(int k=0; k<nnode; k++)//for all nodes
        } //for(int i=0; i<numdof; i++) //for all dof


        for(int line=0; line<numdof*nnode; line++)
        {
          for(int col=0; col<numdof*nnode; col++)
          {
            stiff_relerr(line,col)= fabs( ( pow(elemat1(line,col),2) - pow(stiff_approx(line,col),2) )/ ( (elemat1(line,col) + stiff_approx(line,col)) * elemat1(line,col) ));

            //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
            if ( fabs( stiff_relerr(line,col) ) < h_rel*1000 || isnan( stiff_relerr(line,col)) || elemat1(line,col) == 0) //isnan = is not a number
              stiff_relerr(line,col) = 0;

            if ( stiff_relerr(line,col) > 0)
              outputflag = 1;

          } //for(int col=0; col<numdof*nnode; col++)
        } //for(int line=0; line<numdof*nnode; line++)

        if(outputflag ==1)
        {
          std::cout<<"\n\n acutally calculated stiffness matrix in Element "<<Id()<<": "<< elemat1;
          std::cout<<"\n\n approximated stiffness matrix in Element "<<Id()<<": "<< stiff_approx;
          std::cout<<"\n\n rel error stiffness matrix in Element "<<Id()<<": "<< stiff_relerr;
        }

      } //end of section in which numerical approximation for stiffness matrix is computed
     */

  }
  break;
  case calc_struct_update_istep:
  {
    Theta0_=Theta_;
  }
  break;
  case calc_struct_reset_istep:
  {
    //nothing to do
  }
  break;
  case calc_struct_stress:
  {
    //no stress calculation implemented! Do not crash simulation and just keep quiet!
  }
  break;
  case postprocess_stress:
  {
    //no stress calculation for postprocess. Does not really make sense!
    dserror("No stress output for Spring3!");
  }
  break;
  default:
    dserror("Unknown type of action for Spring3 %d", act);
    break;
  }
  return 0;

}

/*-----------------------------------------------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition (public)                                   mukherjee 04/15|
 *----------------------------------------------------------------------------------------------------------*/

int DRT::ELEMENTS::Spring3::EvaluateNeumann(Teuchos::ParameterList&  params,
    DRT::Discretization&     discretization,
    DRT::Condition&          condition,
    std::vector<int>&        lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  dserror("This method needs to be modified for bio-polymer networks!");
  return (0);
}


/*-----------------------------------------------------------------------------------------------------------*
 | Evaluate PTC damping (public)                                                              mukherjee 04/15|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim, int dof> //number of nodes, number of dimensions of embedding space, number of degrees of freedom per node
void DRT::ELEMENTS::Spring3::EvaluatePTC(Teuchos::ParameterList&   params,
    Epetra_SerialDenseMatrix& elemat1)
{

  dserror("PTC is not implemented for Spring3 elements!");

  return;
} //DRT::ELEMENTS::Spring3::EvaluatePTC

/*--------------------------------------------------------------------------------------*
 | calculation of elastic energy                                         mukherjee 04/15|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::t3_energy(Teuchos::ParameterList&   params,
    std::vector<double>&      disp,
    Epetra_SerialDenseVector* intenergy)
{
  dserror("This method is yet to be configured for bio-polymer networks!");

  return;
}

/*--------------------------------------------------------------------------------------*
 | Nonlinear stiffness                                                   mukherjee 04/15|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::t3_nlnstiffmass(Teuchos::ParameterList&   params,
    std::vector<double>&      vel,
    std::vector<double>&      disp,
    Epetra_SerialDenseMatrix* stiffmatrix,
    Epetra_SerialDenseMatrix* massmatrix,
    Epetra_SerialDenseVector* force)
{
  /*
   * It is observed that for a mixed problems, such is the case for biopolymer network simulations (),
   * the method "Evaluate" hands in the larger matrices and vectors of size of element described in the
   *  input file. For example, if the computational volume contains both Beam and Truss elements. The Evaluate
   *  hand into the method a 12x12 matrix. However, for truss element we need only 6x6. Therefore, an
   *  appropriate mapping needs to be established to ensure proper assemblies for corresponding DOFs.
   *  The algorithm implemented here is valid only for linear elements i.e. element containing two nodes.
   */
  //6x6 Stiffness Matrix of the Truss
  Epetra_SerialDenseMatrix DummyStiffMatrix;
  DummyStiffMatrix.Shape(6,6);
  DummyStiffMatrix.Scale(0);
  //6x6 force vector of the Truss
  Epetra_SerialDenseVector DummyForce;
  DummyForce.Size(6);
  DummyForce.Scale(0);
  //1x6 velocity vector
  LINALG::Matrix<1,6> DummyVel;
  DummyVel.Clear();
  //1x6 displacement vector
  LINALG::Matrix<1,6> DummyDisp;
  DummyDisp.Clear();
  // Map velocity global level into element level
  if (vel.size()>12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (vel.size()==6)
  {
    for (int i=0; i<6; i++)
      DummyVel(i)+=vel[i];
  }
  else if (vel.size()==12)
  {
    for (int i=0; i<3; i++)
    {
      DummyVel(i)+=vel[i];
      DummyVel(i+3)+=vel[i+6];
    }

  }
  // Map displacement global level into element level
  if (disp.size()>12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (disp.size()==6)
  {
    for (int i=0; i<6; i++)
      DummyDisp(i)+=disp[i];
  }
  else if (disp.size()==12)
  {
    for (int i=0; i<3; i++)
    {
      DummyDisp(i)+=disp[i];
      DummyDisp(i+3)+=disp[i+6];
    }
  }

  t3_nlnstiffmass_spring(DummyDisp,DummyStiffMatrix,DummyForce);

  if(params.get<std::string>("internalforces","no")=="yes")
    f_ = Teuchos::rcp(new Epetra_SerialDenseVector(DummyForce));


  /*the following function call applies statistical forces and damping matrix according to the fluctuation dissipation theorem;
   * it is dedicated to the application of truss3 elements in the frame of statistical mechanics problems; for these problems a
   * special vector has to be passed to the element packed in the params parameter list; in case that the control routine calling
   * the element does not attach this special vector to params the following method is just doing nothing, which means that for
   * any ordinary problem of structural mechanics it may be ignored*/


  // Map element level into global 12 by 12 element
  if (force->Length()>12)
    dserror("Vector is larger than 12. Please use different mapping strategy!");
  else if (force->Length()==6)
  {
    for (int i=0; i<6; i++)
      (*force)(i)+=DummyForce(i);
  }
  else if (force->Length()==12)
  {
    for (int i=0; i<3; i++)
    {
      (*force)(i)+=DummyForce(i);
      (*force)(i+6)+=DummyForce(i+3);
    }
  }

  // Map element level into global 12 by 12 element
  if (stiffmatrix->RowDim()>12)
    dserror("Matrix is larger than 12. Please use different mapping strategy!");
  else if(stiffmatrix->RowDim()==6)
  {
    for (int i=0; i<6; i++)
      for (int j=0; j<6; j++)
        (*stiffmatrix)(i,j)+=DummyStiffMatrix(i,j);
  }
  else if(stiffmatrix->RowDim()==12)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
        (*stiffmatrix)(i,j)+=DummyStiffMatrix(i,j);
        (*stiffmatrix)(i,j+6)+=DummyStiffMatrix(i,j+3);
        (*stiffmatrix)(i+6,j+6)+=DummyStiffMatrix(i+3,j+3);
        (*stiffmatrix)(i+6,j)+=DummyStiffMatrix(i+3,j);
      }
  }

  Teuchos::ParameterList sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  if(sdyn.get<std::string> ("DYNAMICTYP", "none")=="StatMech")
    torsion_stiffmass(params, disp, stiffmatrix, force);

  /*****************Compare magnitude of calculated forces********************/
//  LINALG::Matrix<6,1> axialforce;
//  axialforce.Clear();
//  LINALG::Matrix<6,1> moment;
//  moment.Clear();
//  for (int i=0; i<3; i++)
//  {
//    axialforce(i)=(*force)(i);
//    axialforce(i+3)=(*force)(i+6);
//    moment(i)=(*force)(i+3);
//    moment(i+3)=(*force)(i+9);
//  }
//
//  NormForce=axialforce.Norm2();
//  NormMoment=moment.Norm2();
//  if (NormMoment!=0)
//    RatioNormForceMoment=NormForce/NormMoment;
//  else
//    RatioNormForceMoment=0;
//  std::cout<<"Spring3 Element ID="<<this->Id()<<std::endl;
//  std::cout<<"Norm Axial Force="<<NormForce<<std::endl;
//  std::cout<<"Norm Moment="<<NormMoment<<std::endl;
  /*****************Compare magnitude of calculated forces********************/

  return;
}




/*------------------------------------------------------------------------------------------------------------*
 | nonlinear stiffness matrix (private)                                                       mukherjee 04/15|
 | engineering strain measure, in case of spring                                                             |
 *-----------------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::t3_nlnstiffmass_spring(const LINALG::Matrix<1,6>&      DummyDisp,
                                                  Epetra_SerialDenseMatrix& DummyStiffMatrix,
                                                  Epetra_SerialDenseVector& DummyForce)
{
  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // Spring stiffness
  double spring =statmechparams.get<double> ("KTOR1_LINK", 0.0);

//  double spring= 1e-1;

  //current node position (first entries 0 .. 2 for first node, 3 ..5 for second node)
  LINALG::Matrix<6,1> xcurr;

  //auxiliary vector for both internal force and stiffness matrix: N^T_(,xi)*N_(,xi)*xcurr
  LINALG::Matrix<6,1> aux;

  //current nodal position (first
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + DummyDisp(j); //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + DummyDisp(3+j); //second node
  }

  //computing auxiliary vector aux = 4.0*N^T_{,xi} * N_{,xi} * xcurr
  aux(0) = (xcurr(0) - xcurr(3));
  aux(1) = (xcurr(1) - xcurr(4));
  aux(2) = (xcurr(2) - xcurr(5));
  aux(3) = (xcurr(3) - xcurr(0));
  aux(4) = (xcurr(4) - xcurr(1));
  aux(5) = (xcurr(5) - xcurr(2));


  // resulting force scaled by current length
  double forcescalar=spring;

  //computing global internal forces
  for (int i=0; i<6; ++i)
    DummyForce(i) = forcescalar * aux(i);


  //computing linear stiffness matrix
  for (int i=0; i<3; ++i)
  {
    //stiffness entries for first node
    DummyStiffMatrix(i,i)    =  forcescalar;
    DummyStiffMatrix(i,3+i)  = -forcescalar;
    //stiffness entries for second node
    DummyStiffMatrix(i+3,i+3)=  forcescalar;
    DummyStiffMatrix(i+3,i)  = -forcescalar;
  }

  lcurr_ = sqrt(pow(aux(0),2)+pow(aux(1),2)+pow(aux(2),2));

  return;
} // DRT::ELEMENTS::Spring3::bt_nlnstiffmass3




/*-------------------------------------------------------------------------------------------*
 | Calculate change in angle from reference configuration                     mukherjee 09/14|
 *-------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::CalcDeltaTheta(std::vector<double>& disp,
                                                 LINALG::Matrix<1,3>& thetacurr)
{
  //current tangential vector
  std::vector<LINALG::Matrix<3,1> >  tcurrNode;
  tcurrNode.resize(3);
  //current node position (first entries 0,1,2 for first node, 3,4,5 for second node)
  LINALG::Matrix<1,6> xcurr(true);


  // Calculate current directional vector of the truss element (v_1 in derivation)
  // v_1 direction 1---->2 i.e. v1=d2-d1;
  LINALG::Matrix<1,3> diff_disp_curr(true);

  //current nodal position and directional vector of the truss
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[    j];  //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[6 + j];  //second node
    diff_disp_curr(j) = (xcurr(j+3) - xcurr(j));      // v1 = d2 - d1
  }

  //current tangent vector
  std::vector<LINALG::Matrix<3,1> > Tcurr;
  Tcurr.resize(2);
  //current tangent vector
  if (FilamentIsReissner_)
  {
    this->GetCurrTangents(disp, Tcurr);
    for (int node=0; node<2; node++)
    for(int i=0; i<3; i++)
    {
      tcurrNode[node](i)=Tcurr[node](i);
    }
  }
  else
    for (int node=0; node<2; node++)
    {
      tcurrNode[node].Clear();
      for (int j=0; j<3; ++j)
      {
        tcurrNode[node](j)   = trefNode_[node](j) + disp[6*node +3 + j];
      }
    }

  for (int location=0; location<3; location++) // Location of torsional spring. There are three locations
  {
    double dotprod=0.0;
    LINALG::Matrix  <1,3> crossprod(true);
    double CosTheta=0.0;
    double SinTheta=0.0;

    if (location==0)
    {
      double norm_v_curr = diff_disp_curr.Norm2();
      double norm_t1_curr=tcurrNode[location].Norm2();

      for (int j=0; j<3; ++j)
      {
        dotprod +=  tcurrNode[location](j) * diff_disp_curr(j);
      }

      CosTheta = dotprod/(norm_v_curr*norm_t1_curr);

      //Cross Product
      crossprod(0) = tcurrNode[location](1)*diff_disp_curr(2) - tcurrNode[location](2)*diff_disp_curr(1);
      crossprod(1) = tcurrNode[location](2)*diff_disp_curr(0) - tcurrNode[location](0)*diff_disp_curr(2);
      crossprod(2) = tcurrNode[location](0)*diff_disp_curr(1) - tcurrNode[location](1)*diff_disp_curr(0);

      double norm= crossprod.Norm2();
      SinTheta= norm/(norm_v_curr*norm_t1_curr);

    }
    else if (location==1)
    {
      double norm_v_curr = diff_disp_curr.Norm2();
      double norm_t2_curr= tcurrNode[location].Norm2();
      for (int j=0; j<3; ++j)
      {
        dotprod +=  tcurrNode[location](j) * diff_disp_curr(j); // From the opposite direction v_2 =-v_1
      }

      CosTheta = dotprod/(norm_v_curr*norm_t2_curr);

      // cross product
      crossprod(0) = tcurrNode[location](1)*diff_disp_curr(2) - tcurrNode[location](2)*diff_disp_curr(1);
      crossprod(1) = tcurrNode[location](2)*diff_disp_curr(0) - tcurrNode[location](0)*diff_disp_curr(2);
      crossprod(2) = tcurrNode[location](0)*diff_disp_curr(1) - tcurrNode[location](1)*diff_disp_curr(0);
      double norm= crossprod.Norm2();
      SinTheta= norm/(norm_v_curr*norm_t2_curr);
    }
    else // i.e. for calculation of reference angle between t1 & t2
    {
      double norm_t1_curr = tcurrNode[location-2].Norm2();
      double norm_t2_curr=tcurrNode[location-1].Norm2();
      for (int j=0; j<3; ++j)
        dotprod +=  tcurrNode[location-1](j) * tcurrNode[location-2](j);

      CosTheta = dotprod/(norm_t1_curr*norm_t2_curr);

      // cross product
      crossprod(0) = tcurrNode[location-2](1)*tcurrNode[location-1](2) - tcurrNode[location-2](2)*tcurrNode[location-1](1);
      crossprod(1) = tcurrNode[location-2](2)*tcurrNode[location-1](0) - tcurrNode[location-2](0)*tcurrNode[location-1](2);
      crossprod(2) = tcurrNode[location-2](0)*tcurrNode[location-1](1) - tcurrNode[location-2](1)*tcurrNode[location-1](0);

      double norm= crossprod.Norm2();
      SinTheta= norm/(norm_t1_curr*norm_t2_curr);
    }

    double ThetaBoundary1=M_PI/4;
    double ThetaBoundary2=3*M_PI/4;

    if (SinTheta >=0)
    {
      if (CosTheta >= cos(ThetaBoundary1))
        thetacurr(location)=asin(SinTheta);
      else if (CosTheta <= cos(ThetaBoundary2))
        thetacurr(location)=M_PI-asin(SinTheta);
      else
        thetacurr(location)=acos(CosTheta);
    }
    else
      dserror("Angle more than 180 degrees!");

    deltatheta_(location)=thetacurr(location)-ThetaRef_[location];
  }

  return;
}

/*--------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices                                mukherjee 09/14|
 *--------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::torsion_stiffmass(Teuchos::ParameterList&   params,
                                              std::vector<double>&      disp,
                                              Epetra_SerialDenseMatrix* stiffmatrix,
                                              Epetra_SerialDenseVector* force)
{

  // Calculate current directional vector of the truss element (v_1 in derivation)
  // v_1 direction 1---->2 i.e. v1=d2-d1;
  LINALG::Matrix<1,3> diff_disp_curr(true);

  //current node position (first entries 0,1,2 for first node, 3,4,5 for second node)
  LINALG::Matrix<1,6> xcurr(true);
  LINALG::Matrix<1,6> theta_0(true);

  //current tangential vector (first entries 0,1,2 for first node, 3,4,5 for second node)
  LINALG::Matrix<1,3> tcurrNode1(true);
  LINALG::Matrix<1,3> tcurrNode2(true);

  //current nodal position and directional vector of the truss
  for (int j=0; j<3; ++j)
  {
    xcurr(j  )   = Nodes()[0]->X()[j] + disp[    j];  //first node
    xcurr(j+3)   = Nodes()[1]->X()[j] + disp[6 + j];  //second node
    diff_disp_curr(j) = (xcurr(j+3) - xcurr(j));        // v1 = d2 - d1
  }
  std::vector<LINALG::Matrix<3,1> > Tcurr;
  Tcurr.resize(2);
  //current tangent vector
  if (FilamentIsReissner_)
  {
    LINALG::Matrix<3,1> Tcurr1(true);
    LINALG::Matrix<3,1> Tcurr2(true);

    this->GetCurrTangents(disp, Tcurr);
    this->TcurrBeam3r(Tcurr1,Tcurr2);
    for(int i=0; i<3; i++)
    {
      tcurrNode1(i)=Tcurr[0](i);
      tcurrNode2(i)=Tcurr[1](i);
    }
  }
  else
  {
    for (int j=0; j<3; ++j)
    {
      tcurrNode1(j)   = trefNode_[0](j) + disp[3 + j];  //first node
      tcurrNode2(j)   = trefNode_[1](j) + disp[9 + j];  //second node
    }
  }


  deltatheta_.Clear();
  LINALG::Matrix<1,3> thetacurr(true);
  CalcDeltaTheta(disp, thetacurr);
  double ThetaBoundary1=M_PI/4;
  double ThetaBoundary2=3*M_PI/4;

  /********************Check if work potential function is unique*************************************/
//  //Energy potential
//  FAD W = 0.0;
//  // Current tangents
//  LINALG::TMatrix<FAD,3,1> t1(true);
//    LINALG::TMatrix<FAD,3,1> t2(true);
//    LINALG::TMatrix<FAD,3,1> t1_unit(true);
//    LINALG::TMatrix<FAD,3,1> t2_unit(true);
//
//    // Reference tangents
//  LINALG::TMatrix<FAD,3,1> t10(true);
//    LINALG::TMatrix<FAD,3,1> t20(true);
//    LINALG::TMatrix<FAD,3,1> t10_unit(true);
//    LINALG::TMatrix<FAD,3,1> t20_unit(true);
//
//    t1(0)=1.0;
//    t2(1)=1.0;
//    // Compute terms at reference config.
//    t10(0)=1.0;
//    t20(1)=1.0;
//    FAD norm_t20 = pow(FADUTILS::ScalarProduct(t20,t20),0.5);
//    t20_unit.Update(1.0/norm_t20,t20,0.0);
//
//    FAD norm_t10 = pow(FADUTILS::ScalarProduct(t10,t10),0.5);
//    t10_unit.Update(1.0/norm_t10,t10,0.0);
//    double delta=1.0e-6;
//
//
//    for(int i=0;i<=360;i++)
//    {
//      double theta = 360-1.0*i;
//      t2(0)=cos(theta/180.0*M_PI);
//      t2(1)=sin(theta/180.0*M_PI);
//
//      for(int j=0;j<3;j++)
//      {
//        t1(j).diff(j,6);
//        t2(j).diff(j+3,6);
//      }
//
//      FAD norm_t2 = pow(FADUTILS::ScalarProduct(t2,t2),0.5);
//      t2_unit.Update(1.0/norm_t2,t2,0.0);
//
//      FAD norm_t1 = pow(FADUTILS::ScalarProduct(t1,t1),0.5);
//      t1_unit.Update(1.0/norm_t1,t1,0.0);
//
//      W=pow((FADUTILS::ScalarProduct(t1_unit,t2_unit)-FADUTILS::ScalarProduct(t10_unit,t20_unit)),2);
//
//      std::cout<<"Energy Potential FAD="<<W<<std::endl;
//
//      std::cout << "theta: " << theta << std::endl;
//      std::cout << "t1: " << t1 << std::endl;
//      std::cout << "t2: " << t2 << std::endl;
//      std::cout << "t1_unit: " << t1_unit << std::endl;
//      std::cout << "t2_unit: " << t2_unit << std::endl;

      /*%%%%%%%%%%%%%Check with  analytical results %%%%%%%%%%%%%%%%%%%%%*/
//
//      // Reference tangents
//      LINALG::Matrix<3,1> tangent10(true);
//      LINALG::Matrix<3,1> tangent20(true);
//      LINALG::Matrix<3,1> tangent10_unit(true);
//      LINALG::Matrix<3,1> tangent20_unit(true);
//
//      // Compute terms at reference config.
//      tangent10(0)=1.0;
//      tangent20(1)=1.0;
//      tangent20_unit.Scale(1.0/tangent20.Norm2());
//      tangent10_unit.Scale(1.0/tangent10.Norm2());
//
//
//      LINALG::Matrix<3,1> t1_m(true);
//      LINALG::Matrix<3,1> t1_l(true);
//      LINALG::Matrix<3,1> t1_r(true);
//      LINALG::Matrix<3,1> t2_m(true);
//
//      t1_m(0)=1.0;
//      t1_l(0)=1.0;
//      t1_r(0)=1.0;
//
//      for (int j=0;j<3;j++)
//      {
//        t2_m(j)=t2(j).val();
//      }
//
//      for(int k=0;k<3;k++)
//      {
//        t1_l(k)-=delta;
//        t1_r(k)+=delta;
//
//        t1_l.Scale(1.0/t1_l.Norm2());
//        t1_r.Scale(1.0/t1_r.Norm2());
//
//        double W_m = pow((FADUTILS::ScalarProduct(t1_m,t2_m)-FADUTILS::ScalarProduct(tangent20_unit,tangent10_unit)),2);
//        double W_l = pow((FADUTILS::ScalarProduct(t1_l,t2_m)-FADUTILS::ScalarProduct(tangent20_unit,tangent10_unit)),2);
//        double W_r = pow((FADUTILS::ScalarProduct(t1_r,t2_m)-FADUTILS::ScalarProduct(tangent20_unit,tangent10_unit)),2);
//
//        std::cout << "W_m: " << W_m << std::endl;
//        std::cout << "W_l: " << W_l << std::endl;
//        std::cout << "W_r: " << W_r << std::endl;
//
//        std::cout << "left difference: " << (W_m-W_l)/delta << std::endl;
//        std::cout << "right difference: " << (W_r-W_m)/delta << std::endl;
//
//        t1_l(k)+=delta;
//        t1_r(k)-=delta;
//      }

//    }
 /*******************End check work potential****************************************/

  /*%%%%%% Calculate torsional stiffness matrices and forces between tangents at node 1 & node 2 %%%%%%%*/
  //6x6 Stiffness Matrix between tangents at node 1 & node 2
  Epetra_SerialDenseMatrix TorStiffmatrixNode3;
  TorStiffmatrixNode3.Shape(6,6);
  TorStiffmatrixNode3.Scale(0);
  //6x6 force vector between tangents at node 1 & node 2. Contributions to vector {t1, t2} in respective order
  Epetra_SerialDenseVector TorForceNode3;
  TorForceNode3.Size(6);
  TorForceNode3.Scale(0);

  // Calculate torsional stiffness matrices and forces between tangents at node 1 & node 2
  if (thetacurr(2)>=ThetaBoundary1 && thetacurr(2)<=ThetaBoundary2)
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffTangentDot(params, tcurrNode1, tcurrNode2,TorStiffmatrixNode3, TorForceNode3);

  // Uncomment this part to verify the stiffness matrix with FAD type of variable
//     FADMyTorsionalStiffTangentCos(params, double(ThetaRef_[2]), tcurrNode1, tcurrNode2, TorStiffmatrixNode3, TorForceNode3);
  }
  else if ((thetacurr(2)>=0 && thetacurr(2)<ThetaBoundary1) || (thetacurr(2)>ThetaBoundary2 && thetacurr(2)<=M_PI))
  {
    // Calculation based on dot product of vectors
    MyTorsionalStiffTangentDot(params, tcurrNode1, tcurrNode2,TorStiffmatrixNode3, TorForceNode3);

  }
  else
   dserror("Angle out of range!");
  // Map element level into global 12 by 12 element
  if (force->Length()!=12)
    dserror("This element does not need torsional element!");
  else if (force->Length()==12)
  {
    for (int i=0; i<3; i++)
    {
      (*force)(i+3)+=TorForceNode3(i);
      (*force)(i+9)+=TorForceNode3(i+3);
    }
  }
  // Map element level into global 12 by 12 element
  if (stiffmatrix->RowDim()!=12)
    dserror("This element does not require torsional element!");
  else if(stiffmatrix->RowDim()==12)
  {
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
      {
        (*stiffmatrix)(i+3,j+3)+=TorStiffmatrixNode3(i,j);
        (*stiffmatrix)(i+9,j+3)+=TorStiffmatrixNode3(i+3,j);
        (*stiffmatrix)(i+3,j+9)+=TorStiffmatrixNode3(i,j+3);
        (*stiffmatrix)(i+9,j+9)+=TorStiffmatrixNode3(i+3,j+3);
      }
  }
  return;
}

/*-------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices and forces between tangents of beam       mukherjee 09/14|
 *-------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::MyTorsionalStiffTangentCos(Teuchos::ParameterList&   params,
                                                    double theta,
                                                    double deltatheta,
                                                    LINALG::Matrix<1,3> & tcurr1,
                                                    LINALG::Matrix<1,3> & tcurr2,
                                                    Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                    Epetra_SerialDenseVector & TorForce)
{
  // Norms of the tangential vectors and directional displacement vector
  double norm_t2= tcurr2.Norm2();
  double norm_t1= tcurr1.Norm2();
  double s=cos(theta);

  // Calculate gradient of theta w.r.t. unknown dofs , in this case {t1,t2}
  LINALG::Matrix<1,3> A(true);  // Identical to "A" in derivation
  LINALG::Matrix<1,3> B(true);  // Identical to "B" in derivation

  for(int i=0; i<3; i++)
  {
    A(i)+=tcurr1(i)*s/(pow(norm_t1,2)*sin(theta)) - tcurr2(i)/(norm_t1*norm_t2*sin(theta));
    B(i)+=tcurr2(i)*s/(pow(norm_t2,2)*sin(theta)) - tcurr1(i)/(norm_t1*norm_t2*sin(theta));
  }

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::Matrix<1,6> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_c(j)= A(j);
    aux_c(j+3)= B(j);
  }

  // Spring stiffness
  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  double spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);

  // Calculate torsional forces
  for(int j=0; j<6; j++)
    TorForce(j)=spring*deltatheta*aux_c(j);

  //%%%%%%% Calculation of stiffness matrix %%%%%%%%

  // Calculate auxiliary vectors for computation of stiffnessmatrices
  LINALG::Matrix<3,3> dadt1(true);
  LINALG::Matrix<3,3> dadt2(true);  // Identical to "C" in derivation
  LINALG::Matrix<3,3> dbdt1(true);
  LINALG::Matrix<3,3> dbdt2(true);  // Identical to "D" in derivation

  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      dadt1(i,j)+= (i==j)* s/(pow(norm_t1,2)* sin(theta)) - 2*s*tcurr1(i)*tcurr1(j)/(sin(theta)*pow(norm_t1,4)) - tcurr1(i)*A(j)/(pow(norm_t1,2)*pow(sin(theta),2)) + tcurr2(i)*tcurr1(j)/(pow(norm_t1,3)*norm_t2*sin(theta)) + tcurr2(i)*A(j)*s/(norm_t1*norm_t2*pow(sin(theta),2));

      dadt2(i,j)+= tcurr2(i)*tcurr2(j)/(norm_t1*pow(norm_t2,3)*sin(theta)) + tcurr2(i)*B(j)*s/(norm_t1*norm_t2*pow(sin(theta),2)) -(i==j)/(norm_t1*norm_t2*sin(theta)) -tcurr1(i)*B(j)/(pow(norm_t1,2)*pow(sin(theta),2));

      dbdt1(i,j)+= -(i==j)/(norm_t1*norm_t2*sin(theta))+ tcurr1(i)*tcurr1(j)/(pow(norm_t1,3)*norm_t2*sin(theta))+ tcurr1(i)*A(j)*s/(norm_t1*norm_t2*pow(sin(theta),2)) -tcurr2(i)*A(j)/(pow(norm_t2,2)*pow(sin(theta),2));

      dbdt2(i,j)+= (i==j)*s/(pow(norm_t2,2)*sin(theta))- 2*s*tcurr2(i)*tcurr2(j)/(sin(theta)*pow(norm_t2,4))- tcurr2(i)*B(j)/(pow(norm_t2,2)*pow(sin(theta),2))+ tcurr1(i)*tcurr2(j)/(norm_t1*pow(norm_t2,3)*sin(theta))+ tcurr1(i)*B(j)*s/(norm_t1*norm_t2*pow(sin(theta),2));
    }


  // Create torsional siffness matrix
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
    {
      TorStiffmatrix(i,j)+= spring*aux_c(i)*aux_c(j) + spring* deltatheta * dadt1(i,j);
      TorStiffmatrix(i,j+3)+= spring*aux_c(i)*aux_c(j+3) + spring* deltatheta * dadt2(i,j);
      TorStiffmatrix(i+3,j)+= spring*aux_c(i+3)*aux_c(j) + spring* deltatheta * dbdt1(i,j);
      TorStiffmatrix(i+3,j+3)+= spring*aux_c(i+3)*aux_c(j+3) + spring* deltatheta * dbdt2(i,j);
    }
  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::MyTorsionalStiffTangentDot(Teuchos::ParameterList&   params,
                                                      LINALG::Matrix<1,3> & tangentcurr1,
                                                      LINALG::Matrix<1,3> & tangentcurr2,
                                                      Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                      Epetra_SerialDenseVector & TorForce)
{

  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //total no of dofs
  const int totdof = 6; // (t1,t2)

  //FAD calculated stiff matrix for validation purposes
  LINALG::TMatrix<FAD,totdof,totdof> stiffmatrix_check(true);
  LINALG::TMatrix<FAD,1,totdof> force_check(true);

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // current tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tcurr1(true);
  // current tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tcurr2(true);
  // reference tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tref1(true);
  // reference tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tref2(true);
  // unit current tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tcurr1_unit(true);
  // unit current tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tcurr2_unit(true);
  // unit reference tangent at Node 1
  LINALG::TMatrix<FAD,3,1> tref1_unit(true);
  // unit reference tangent at Node 2
  LINALG::TMatrix<FAD,3,1> tref2_unit(true);

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tangentcurr1(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= tangentcurr2(i);
    disp_tot[i+3].diff(i+3, totdof);
    tcurr1(i)=disp_tot[i];
    tcurr2(i)=disp_tot[i+3];
    tref1(i)=trefNode_[0](i); // trefNode_ nodal index varies from [0,1]
    tref2(i)=trefNode_[1](i);
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_tref1=pow(FADUTILS::ScalarProduct(tref1,tref1),0.5);
  FAD norm_tref2=pow(FADUTILS::ScalarProduct(tref2,tref2),0.5);
  FAD norm_tcurr1=pow(FADUTILS::ScalarProduct(tcurr1,tcurr1),0.5);
  FAD norm_tcurr2=pow(FADUTILS::ScalarProduct(tcurr2,tcurr2),0.5);

  // Recalculate unit vectors
  tref1_unit.Update(1.0/norm_tref1,tref1,0.0);
  tref2_unit.Update(1.0/norm_tref2,tref2,0.0);
  tcurr1_unit.Update(1.0/norm_tcurr1,tcurr1,0.0);
  tcurr2_unit.Update(1.0/norm_tcurr2,tcurr2,0.0);

  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();
  // Spring stiffness
  FAD spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);

  //computing energy potential (equation 3.3)
  FAD W=0.5*spring*pow((FADUTILS::ScalarProduct(tcurr1_unit,tcurr2_unit)-FADUTILS::ScalarProduct(tref1_unit,tref2_unit)),2);

  //Compute linearization with FAD for checking
  for (int i=0; i<totdof; i++)
    force_check(i)=W.dx(i);

  // Analytical matrix for computation of forces
  LINALG::TMatrix<FAD,1,3> A(true);
  LINALG::TMatrix<FAD,1,3> B(true);
  // Auxiliary matrix
  LINALG::TMatrix<FAD,3,3> t1_aux(true);
  LINALG::TMatrix<FAD,3,3> t2_aux(true);

  // Calculate dot product
  FAD t1t2curr_unit=FADUTILS::ScalarProduct(tcurr1_unit,tcurr2_unit);
  FAD t1t2ref_unit=FADUTILS::ScalarProduct(tref1_unit,tref2_unit);
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
    {
      t1_aux(i,j)= (i==j)/norm_tcurr1 -tcurr1(i)*tcurr1(j)/pow(norm_tcurr1,3);
      t2_aux(i,j)= (i==j)/norm_tcurr2 -tcurr2(i)*tcurr2(j)/pow(norm_tcurr2,3);
      A(i)+=spring*(t1t2curr_unit-t1t2ref_unit)*t1_aux(i,j)*tcurr2_unit(j);
      B(i)+=spring*(t1t2curr_unit-t1t2ref_unit)*t2_aux(i,j)*tcurr1_unit(j);
    }

  // FAD force vector for computation of forces
  LINALG::TMatrix<FAD,1,totdof> FADforce(true);
  for (int i=0; i<3; i++)
  {
    FADforce(i)=A(i);
    TorForce(i)=A(i).val();
    FADforce(i+3)=B(i);
    TorForce(i+3)=B(i).val();
  }

  // FAD force vector for computation of forces
  LINALG::TMatrix<FAD,totdof,totdof> FADstiffMatrix(true);

  for (int i=0; i<totdof; i++)
    for (int j=0; j<totdof; j++)
      TorStiffmatrix(i,j)=FADforce(i).dx(j);

  return;
}



/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::FADMyTorsionalStiffTangentCos(Teuchos::ParameterList&   params,
                                                      double theta_0,
                                                      LINALG::Matrix<1,3> & tcurrNode1,
                                                      LINALG::Matrix<1,3> & tcurrNode2,
                                                      Epetra_SerialDenseMatrix & TorStiffmatrix,
                                                      Epetra_SerialDenseVector & TorForce)
{
  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //FAD calculated stiff matrix for validation purposes
  LINALG::TMatrix<FAD,6,6> stiffmatrix_check;
  stiffmatrix_check.Clear();
  LINALG::TMatrix<FAD,1,6> force_check;

  //total no of dofs
  const int totdof = 6;

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // diff_disp = d2 -d1
  LINALG::TMatrix<FAD,1,3> tangent2;
  // t1 = tangent
  LINALG::TMatrix<FAD,1,3> tangent1;

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tcurrNode1(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= tcurrNode2(i);
    disp_tot[i+3].diff(i+3, totdof);
    tangent1(0,i)=disp_tot[i];
    tangent2(0,i)=disp_tot[i+3];
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_t1= 0.0;
  FAD norm_t2= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_t1+=pow(tangent1(0,i),2);
    norm_t2+=pow(tangent2(0,i),2);
  }
  norm_t1=pow(norm_t1,0.5);
  norm_t2=pow(norm_t2,0.5);
  if (norm_t1==0.0)
    norm_t1=1.0e-14;
  if (norm_t2==0.0)
    norm_t2=1.0e-14;

  //computing the change of angle theta (equation 3.3)
  FAD deltatheta=0.0;
  FAD theta=0.0;
  FAD dotprod=0.0;
  FAD s=0.0;

  for (int j=0; j<3; ++j)
    dotprod +=  tangent1(0,j) * tangent2(0,j);

  s = dotprod/(norm_t1*norm_t2);

  // Owing to round-off errors the variable s can be slightly
  // outside the admissible range [-1.0;1.0]. We take care for this
  // preventing potential floating point exceptions in acos(s)
  if (s>1.0)
  {
    if ((s-1.0)>1.0e-14)
      dserror("s out of admissible range [-1.0;1.0]");
    else // tiny adaptation of s accounting for round-off errors
      s = 1.0-1.0e-14;
  }
  if (s<-1.0)
  {
    if ((s+1.0)<-1.0e-14)
      dserror("s out of admissible range [-1.0;1.0]");
    else // tiny adaptation of s accounting for round-off errors
      s = -1.0+1.0e-14;
  }
  if (s==0.0)
    s = 1.0e-14;
  else if (s==1.0)
    s = 1-1.0e-14;
  else if (s==-1.0)
    s = -1+1.0e-14;

  theta=acos(s);

  deltatheta=theta-theta_0;

  Teuchos::ParameterList statmechparams = DRT::Problem::Instance()->StatisticalMechanicsParams();

  // Calculate auxiliary vectors for computation of forces and stiffnessmatrix
  LINALG::TMatrix<FAD,1,3> aux_a(true);  // Identical to "a" in derivation
  LINALG::TMatrix<FAD,1,3> aux_b(true);  // Identical to "b" in derivation
  LINALG::TMatrix<FAD,1,6> aux_c(true);  // Identical to "c" in derivation

  for(int j=0; j<3; j++)
  {
    aux_a(0,j)= tangent1(0,j)*s/(pow(norm_t1,2)*sin(theta)) - tangent2(0,j)/(norm_t2*norm_t1*sin(theta));
    aux_b(0,j)= tangent2(0,j)*s/(pow(norm_t2,2)*sin(theta)) - tangent1(0,j)/(norm_t2*norm_t1*sin(theta));
    aux_c(0,j)= aux_a(0,j);
    aux_c(0,j+3)= aux_b(0,j);
  }

  // Spring stiffness
  FAD spring =statmechparams.get<double> ("KTOR2_LINK", 0.0);

  // Calculate torsional forces
  for(int j=0; j<6; j++)
    force_check(0,j)=spring*deltatheta*aux_c(0,j);

  //shifting values from fixed size matrix to epetra matrix *stiffmatrix
  for(int i = 0; i < totdof; i++)
  {
    for(int j = 0; j < totdof; j++)
    {
      stiffmatrix_check(i,j) = force_check(0,i).dx(j);
    }
  } //for(int i = 0; i < dofpn*nnode; i++)

  LINALG::TMatrix<FAD,6,6> stiff_relerr;
  stiff_relerr.Clear();
  for(int row=0; row<6; row++)
  {
    for(int col=0; col<6; col++)
    {
      if (fabs(stiffmatrix_check(row,col))<1.0e-30 || fabs(TorStiffmatrix(row,col))<1.0e-30 || fabs(stiffmatrix_check(row,col))==0.0 || fabs(TorStiffmatrix(row,col))==0 )
        stiff_relerr(row,col) = 0;
      else
        stiff_relerr(row,col)= fabs(stiffmatrix_check(row,col) - TorStiffmatrix(row,col));

      //suppressing small entries whose effect is only confusing and NaN entires (which arise due to zero entries)
      if ( fabs(stiff_relerr(row,col)) < 1.0e-15 ) //isnan = is not a number
        stiff_relerr(row,col) = 0;
    } //for(int col=0; col<3*nnode; col++)
  } //for(int line=0; line<3*nnode; line++)


  std::cout<<"\n\n original stiffness matrix corresponding to tangential dofs: "<< std::endl;
  for(int i = 0; i< 6; i++)
  {
    for(int j = 0; j< 6; j++)
    {
      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << TorStiffmatrix(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }

  std::cout<<"\n\n analytical stiffness matrix corresponding to tangential dofs: "<< std::endl;
  for(int i = 0; i< 6; i++)
  {
    for(int j = 0; j< 6; j++)
    {
      std::cout << std::setw(9) << std::setprecision(4) << std::scientific << stiffmatrix_check(i,j)<<"\t";
    }
    std::cout<<std::endl;
  }

  //std::cout<<"\n\n FAD stiffness matrix"<< stiffmatrix_check;
  std::cout<<"\n\n rel error of stiffness matrix corresponding to tangential dofs"<< stiff_relerr;
  std::cout<<"Analytical Force corresponding to tangential dofs= "<< force_check << std::endl;
  std::cout<<"Original Force corresponding to tangential dofs="<<TorForce<<std::endl;
  return;
}

/*--------------------------------------------------------------------------------------------------*
 | Calculate torsional stiffness matrices between tangent of beam and directional vector of truss   |
 | using automatic differentiation                                                  mukherjee 09/14 |
 *--------------------------------------------------------------------------------------------------*/
void DRT::ELEMENTS::Spring3::FADThetaLinearisation(LINALG::Matrix<1,3> & diff_disp,
                                                  LINALG::Matrix<1,3> & tcurr,
                                                  LINALG::Matrix<1,3> & A,
                                                  LINALG::Matrix<1,3> & B)
{
  //1x6 Auxilary vector to verify linearisation of angle w.r.t. variables. Contributions to vector {t1, v1} in respective order
  LINALG::Matrix<1,6> AuxilaryVector(true);
  for (int i=0; i<3; i++)
  {
    AuxilaryVector(i)=A(i);
    AuxilaryVector(i+3)=B(i);
  }

  //see also so_nstet_nodalstrain.cpp, so_nstet.H, autodiff.cpp and autodiff.H
  //FAD calculated linearisation for validation purposes
  LINALG::TMatrix<FAD,1,6> vector_check;

  //total no of dofs
  const int totdof = 6;

  //matrix for current positions and tangents
  std::vector<FAD> disp_tot(totdof,0.0);

  // diff_disp = d2 -d1
  LINALG::TMatrix<FAD,1,3> FADdiff_disp;
  // t1 = tangent
  LINALG::TMatrix<FAD,1,3> tangent;

  for (int i=0; i<3; i++)
  {
    disp_tot[i]=tcurr(i);
    disp_tot[i].diff(i,totdof);
    disp_tot[i+3]= diff_disp(i);
    disp_tot[i+3].diff(i+3, totdof);
    tangent(i)=disp_tot[i];
    FADdiff_disp(i)=disp_tot[i+3];
  }

  // Norms of the tangential vectors and directional displacement vector
  FAD norm_v= 0.0;
  FAD norm_t= 0.0;
  for (int i=0; i<3; i++)
  {
    norm_v+=pow(FADdiff_disp(i),2);
    norm_t+=pow(tangent(i),2);
  }
  norm_v=pow(norm_v,0.5);
  norm_t=pow(norm_t,0.5);

  //computing the change of angle theta (equation 3.3)
  FAD Theta=0.0;
  FAD dotprod=0.0;
  FAD CosTheta=0.0;

  for (int j=0; j<3; ++j)
    dotprod +=  tangent(j) * FADdiff_disp(j);

  CosTheta = dotprod/(norm_v*norm_t); // Cosine of angle

  Theta=acos(CosTheta);

  for(int i = 0; i < totdof; i++)
  {
    vector_check(i) = Theta.dx(i);
  }


  LINALG::TMatrix<FAD,1,totdof> vector_relerr;
  vector_relerr.Clear();

  for(int index=0; index<totdof; index++)
  {
    if (fabs(vector_check(index))<1.0e-30 || fabs(AuxilaryVector(index))<1.0e-30 )
      vector_relerr(index) = 0;
    else
      vector_relerr(index)= fabs(AuxilaryVector(index) - vector_check(index));
  }

  //  std::cout<<"\n\n rel error of vector"<< vector_relerr;
  //  std::cout<<"Analytical vector= "<< vector_check<< std::endl;
  //  std::cout<<"Numerical Force="<<AuxilaryVector<<std::endl;
  return;
}

/*-----------------------------------------------------------------------------------------------------------*
 | shifts nodes so that proper evaluation is possible even in case of periodic boundary conditions; if two   |
 | nodes within one element are separated by a periodic boundary, one of them is shifted such that the final |
 | distance in R^3 is the same as the initial distance in the periodic space; the shift affects computation  |
 | on element level within that very iteration step, only (no change in global variables performed)          |                                 |
 |                                                                                   (public) mukherjee 04/14|
 *----------------------------------------------------------------------------------------------------------*/
template<int nnode, int ndim> //number of nodes, number of dimensions
inline void DRT::ELEMENTS::Spring3::NodeShift(Teuchos::ParameterList& params,  //!<parameter list
                                             std::vector<double>&    disp) //!<element disp vector
{
  /* get number of degrees of freedom per node; note: the following function assumes the same number of degrees
   * of freedom for each element node*/
  int numdof = NumDofPerNode(*(Nodes()[0]));
  if (nnode==2 && disp.size()==12)
    numdof = 6;
  double time = params.get<double>("total time",0.0);
  double starttime = params.get<double>("STARTTIMEACT",0.0);
  double dt = params.get<double>("delta time");
  double shearamplitude = params.get<double> ("SHEARAMPLITUDE", 0.0);
  int curvenumber = params.get<int> ("CURVENUMBER", -1)-1;
  int dbcdispdir = params.get<int> ("DBCDISPDIR", -1)-1;
  Teuchos::RCP<std::vector<double> > defvalues = Teuchos::rcp(new std::vector<double>(3,0.0));
  Teuchos::RCP<std::vector<double> > periodlength = params.get("PERIODLENGTH", defvalues);
  INPAR::STATMECH::DBCType dbctype = params.get<INPAR::STATMECH::DBCType>("DBCTYPE", INPAR::STATMECH::dbctype_std);
  bool shearflow = false;
  if(dbctype==INPAR::STATMECH::dbctype_shearfixed || dbctype==INPAR::STATMECH::dbctype_sheartrans || dbctype==INPAR::STATMECH::dbctype_affineshear)
    shearflow = true;
  /*only if periodic boundary conditions are in use, i.e. params.get<double>("PeriodLength",0.0) > 0.0, this
   * method has to change the displacement variables*/
  if(periodlength->at(0) > 0.0)
    //loop through all nodes except for the first node which remains fixed as reference node
    for(int i=1;i<nnode;i++)
    {
      for(int dof= ndim - 1; dof > -1; dof--)
      {
        /*if the distance in some coordinate direction between some node and the first node becomes smaller by adding or subtracting
         * the period length, the respective node has obviously been shifted due to periodic boundary conditions and should be shifted
         * back for evaluation of element matrices and vectors; this way of detecting shifted nodes works as long as the element length
         * is smaller than half the periodic length*/
        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) + periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] += periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the upper z-surface and reenter through the lower
           *z-surface, the shear has to be substracted from nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+dbcdispdir] += shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }

        if( fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - periodlength->at(dof) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) < fabs( (Nodes()[i]->X()[dof]+disp[numdof*i+dof]) - (Nodes()[0]->X()[dof]+disp[numdof*0+dof]) ) )
        {
          disp[numdof*i+dof] -= periodlength->at(dof);

          /*the upper domain surface orthogonal to the z-direction may be subject to shear Dirichlet boundary condition; the lower surface
           *may be fixed by DBC. To avoid problmes when nodes exit the domain through the lower z-surface and reenter through the upper
           *z-surface, the shear has to be added to nodal coordinates in that case */
          if(shearflow && dof == 2 && curvenumber >= 0 && time>starttime && fabs(time-starttime)>dt/1e4)
            disp[numdof*i+dbcdispdir] -= shearamplitude*DRT::Problem::Instance()->Curve(curvenumber).f(time);
        }
      }
    }


  return;

}//DRT::ELEMENTS::Spring3::NodeShift

