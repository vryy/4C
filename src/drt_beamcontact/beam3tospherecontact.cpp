/*----------------------------------------------------------------------------*/
/*!
\file beam3tospherecontact.cpp

\brief class to handle contact between a 3D beam element and a rigid sphere

\level 3

\maintainer Maximilian Grill
*/
/*----------------------------------------------------------------------------*/

#include "beam3tospherecontact.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3/beam3r.H"
#include "../drt_beam3/beam3eb.H"
#include "../drt_beaminteraction/beam3contact_defines.H"
#include "../drt_rigidsphere/rigidsphere.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                     grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::Beam3tospherecontact(const DRT::Discretization& pdiscret,
                                    const DRT::Discretization& cdiscret,
                                    const std::map<int,int>& dofoffsetmap,
                                    DRT::Element* element1,
                                    DRT::Element* element2):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffsetmap_(dofoffsetmap),
element1_(element1),
element2_(element2),
contactflag_(false)
{
  ele1pos_.Clear();
  ele2pos_.Clear();

  fc1_.Clear();
  fc2_.Clear();

  gap_=0.0;

  // initialize class variables for contact point coordinates
  x1_.Clear();
  x2_.Clear();
  normal_.Clear();

//  firstcall_ = true;

  const DRT::ElementType & eot1 = element1_->ElementType();
  const DRT::ElementType & eot2 = element2_->ElementType();

  if(eot1 != DRT::ELEMENTS::Beam3Type::Instance() and eot1 != DRT::ELEMENTS::Beam3rType::Instance()
     and eot1 != DRT::ELEMENTS::Beam3ebType::Instance())
  {
    dserror("How did you get here? element1_ has to be of type beam3, beam3r or beam3eb!!!");
  }

  if(eot2 != DRT::ELEMENTS::RigidsphereType::Instance())
    dserror("How did you get here? element2_ has to be of type rigidsphere!!!");

  nodalcontactflag_.assign(element1_->NumNode(),false);

  return;
}

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                                grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::Beam3tospherecontact(const Beam3tospherecontact& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffsetmap_(old.dofoffsetmap_),
element1_(old.element1_),
element2_(old.element2_),
ele1pos_(old.ele1pos_),
ele2pos_(old.ele2pos_)
{
  dserror("ERROR: Copy constructor incomplete");
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                            grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, double& pp)
{
  //**********************************************************************
  // Evaluation of contact forces and stiffness
  //**********************************************************************
  // (1) Closest Point Projection
  // (2) Compute some auxiliary quantities
  //     -> normal vector, gap, shape functions, contact flag,
  //     -> linearizations of all geometric quantities
  // (3) Compute contact forces and stiffness
  //     -> stiffness terms are directly assembled to global matrix
  //     -> contact forces are only returned as global vector
  // (4) Perform some finite difference checks
  //     -> only if the flag BEAMCONTACTFDCHECKS is defined
  // set class variable for status of gapfunction

  //Note: element1 is always the beam, element2 is the sphere

  ClosestPointProjection();

  // vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues> N1_i;
  LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues> N1_i_xi;
  LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues> N1_i_xixi;

  // coords and derivatives of the two contact points
  LINALG::TMatrix<TYPE,3,1> x1;                            // = x1
  LINALG::TMatrix<TYPE,3,1> x2;                            // = x2
  LINALG::TMatrix<TYPE,3,1> dx1;                            // = x1,xi
  LINALG::TMatrix<TYPE,3,1> ddx1;                          // = x1,xixi

  // initialize
  LINALG::TMatrix<TYPE,3,1> normal;
  TYPE gap= 0.0;
  TYPE norm = 0.0;

  bool validclosestpointprojection = true;

  if (abs(FADUTILS::CastToDouble(xicontact_))< (1.0 + XIETATOL))    // ToDo when to reset nodalcontactflag_?
  {
    nodalcontactflag_[0]=false;
    nodalcontactflag_[1]=false;
  }
  else
  {
    contactflag_ = false;
    validclosestpointprojection = false;
  }

  if(validclosestpointprojection)
  {
    //**********************************************************************
    // (1) Compute some auxiliary quantities
    //**********************************************************************

    // call function to fill variables for shape functions and their derivatives
    GetShapeFunctions(N1_i,N1_i_xi,N1_i_xixi,xicontact_);

    // call function to fill variables with coords and derivs of the contact point
    ComputeCoordsAndDerivs(x1,x2,dx1,ddx1,N1_i,N1_i_xi,N1_i_xixi);

    // call function to compute scaled normal and gap in possible contact point
    ComputeNormal(normal,gap,norm,x1,x2);

    // call function to evaluate contact status
    CheckContactStatus(pp);

    //**********************************************************************
    // (2) Compute contact forces and stiffness
    //**********************************************************************

    // call function to evaluate and assemble contact forces
    EvaluateFcContact(pp,gap,normal,fint,N1_i,contactflag_);
    // call function to evaluate and assemble contact stiffness
    EvaluateStiffcContact(pp,gap,normal,norm,stiffmatrix,x1,x2,dx1,
                          ddx1,N1_i,N1_i_xi,N1_i_xixi,contactflag_);

  }
  else    // TODO do this only for the nodes at beam end points in case of C1-continuous beam centerline representation
  {
    dsassert(numnodes==2,"Beamtospherecontact: Check for nodal contact only implemented for 2-noded beam elements!");
    //************************************************************************
    // NEW: NODAL-BEAM-TO-SPHERE CONTACT
    //************************************************************************
    // distance of beam nodes and center of sphere is computed
    // -> penalty contribution if distance is smaller than sum of radii
    // (point-to-point contact: xi=+-1, hence no contribution from linearization of xi)
    // this also has consequences for AUTOMATICDIFF: do not use xicontact_ here!
    TYPE XiContact = 0.0;
    for (int inode=0; inode<2; ++inode)
    {
      /* if inode=0 is active contact (gap<0), we do NOT want to run the loop for inode=1
       * because the case validclosestpointprojection=false and nodalcontactflag_[0]=true and nodalcontactflag_[1]=true
       * is highly unlikely and class variables like gap_ are overwritten if we run this loop again
       * this leads to wrong output for gap_ of active BTSPH pairs in Beam3contactmanager */
      if (nodalcontactflag_[0]!=true)
      {
        // Set XiContact to +-1.0
        switch(inode)
        {
          case 0:
            XiContact=-1.0;   // first node
            break;
          case 1:
            XiContact=1.0;    // second node
            break;
          default:
            dserror("This must not happen!");
            break;
        }

        // Now do exactly the same as for "normal" contact based on closest point projection
        // (except for contribution from linearization of xi)

        // TODO: make more efficient: use knowledge that contact point is a node -> matrix of shape functions, etc. simplifies

        //**********************************************************************
        // (1) Compute some auxiliary quantities
        //**********************************************************************

        // call function to fill variables for shape functions and their derivatives
        GetShapeFunctions(N1_i,N1_i_xi,N1_i_xixi,XiContact);

        // call function to fill variables with coords and derivs of the contact point
        ComputeCoordsAndDerivs(x1,x2,dx1,ddx1,N1_i,N1_i_xi,N1_i_xixi);

        // call function to compute scaled normal and gap in possible contact point
        ComputeNormal(normal,gap,norm,x1,x2);

        // evaluate nodal contact status
        if ( gap_ < 0)
        {
          nodalcontactflag_[inode]=true;
        }
        else nodalcontactflag_[inode] = false;

        //**********************************************************************
        // (2) Compute contact forces and stiffness
        //**********************************************************************

        // call function to evaluate and assemble contact forces
        EvaluateFcContact(pp,gap,normal,fint,N1_i,nodalcontactflag_[inode]);

        // call function to evaluate and assemble contact stiffness
        EvaluateStiffcContact(pp,gap,normal,norm,stiffmatrix,x1,x2,dx1,
                              ddx1,N1_i,N1_i_xi,N1_i_xixi,nodalcontactflag_[inode],false);
      }

    }
  }

  return GetContactFlag();
}

/*----------------------------------------------------------------------*
 |  Closest point projection                                 grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ClosestPointProjection()
{
  // local variable for beam element coordinate
  TYPE eta=0.0;

  // vectors for shape functions and their derivatives
  LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues> N1_i;
  LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues> N1_i_xi;
  LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues> N1_i_xixi;

  // coords and derivatives of the two contact points
  LINALG::TMatrix<TYPE,3,1> x1;                            // = x1
  LINALG::TMatrix<TYPE,3,1> x2;                            // = x2
  LINALG::TMatrix<TYPE,3,1> dx1;                            // = x1,xi
  LINALG::TMatrix<TYPE,3,1> ddx1;                          // = x1,xixi
  LINALG::TMatrix<TYPE,3,1> delta_x;                       // = x1 - x2

  // initialize function f and Jacobian df for Newton iteration
  TYPE f;
  TYPE df;

  // initial scalar residual (L2-norm of f)
  double residual = 0.0;

  int iter=0;

  //**********************************************************************
  // local Newton iteration
  //**********************************************************************
  for (int i=0;i<BEAMCONTACTMAXITER;++i)
  {
    iter++;

    // reset shape function variables to zero
    N1_i.Clear();
    N1_i_xi.Clear();
    N1_i_xixi.Clear();

    // update shape functions and their derivatives
    GetShapeFunctions(N1_i,N1_i_xi,N1_i_xixi,eta);

    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(x1,x2,dx1,ddx1,N1_i,N1_i_xi,N1_i_xixi);

    // compute delta_x = x1-x2
    for (int j=0; j<3; ++j)
      delta_x(j) = x1(j) - x2(j);

    // compute norm of difference vector to scale the equations
    // (this yields better conditioning)
    // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
    // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
    double norm_delta_x = FADUTILS::CastToDouble(FADUTILS::VectorNorm<3>(delta_x));

    // the closer the beams get, the smaller is norm
    // norm is not allowed to be too small, else numerical problems occur
    if (norm_delta_x < NORMTOL)
    {
      dserror("Contact points x1 and x2 are identical. Choose smaller time step!");
    }

    // evaluate f at current eta
    EvaluateOrthogonalityCondition(f,delta_x,norm_delta_x,dx1);

    // compute the scalar residuum
    residual = abs(FADUTILS::CastToDouble(f));

    // check if Newton iteration has converged
    if (residual < BEAMCONTACTTOL) break;

    // evaluate Jacobian of f at current eta
    EvaluateLinOrthogonalityCondition(df,delta_x,norm_delta_x,dx1,ddx1);

    // singular df
    if (abs(df)<COLINEARTOL)
    {
      dserror("No solution for Closest Point Projection!");
    }
    // regular df (inversion possible)
    else
    {
      // update element coordinates of contact point
      eta += -f/df;
    }
  }
  //**********************************************************************

  // Newton iteration unconverged
  if (residual > BEAMCONTACTTOL)
  {
    eta = 1e+12;
  }

  // store final result and return
  xicontact_=eta;

  #ifdef AUTOMATICDIFF
  // Set xi1_ as (additional) primary variable for automatic differentiation
  // The dependence between the infinitesimal changes delta xi1_ and the
  // the increments of the primary displacement variables delta disp have to be given explicitly, since
  // no explicit relation between the finite quantities xi1_ and disp exists.
  // The latter would have been necessary if the full linearization had to be computed directly with Sacado!!!

  // The 3*numnodes*numnodalvalues+3 primary DoFs are the components of the nodal positions / tangents. The additional
  // degree of freedom (+1) represents the dependency on the beam parameter coordinate xi, which is necessary in beam contact.
  xicontact_.diff((3*numnodes*numnodalvalues+3+1)-1,3*numnodes*numnodalvalues+3+1);
  #endif

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate function f in CPP                               grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::EvaluateOrthogonalityCondition(TYPE& f, const LINALG::TMatrix<TYPE,3,1>& delta_x, const double norm_delta_x, const LINALG::TMatrix<TYPE,3,1>& dx1)
{
  // reset f
  f=0.0;

  // evaluate f
  // see Wriggers, Computational Contact Mechanics, equation (12.5)

  for (int i=0;i<3;i++)
  {
    f += delta_x(i)*dx1(i) / norm_delta_x;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df in CPP                              grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::EvaluateLinOrthogonalityCondition(TYPE& df,
                                                         LINALG::TMatrix<TYPE,3,1>& delta_x,
                                                         const double norm_delta_x,
                                                         const LINALG::TMatrix<TYPE,3,1>& dx1,
                                                         const LINALG::TMatrix<TYPE,3,1>& ddx1)

{
  // reset df
  df=0;

  // evaluate df
  // see Wriggers, Computational Contact Mechanics, equation (12.7)
  for(int i=0;i<3;i++)
  {
    df += (dx1(i)*dx1(i) + delta_x(i)*ddx1(i)) / norm_delta_x;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Check if contact is active or inactive                    grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::CheckContactStatus(double& pp)
{
  // check contact condition
   contactflag_ = (gap_ < 0) ? true : false;

  return;
}

/*----------------------------------------------------------------------*
 |  Get global dofs of a node                                grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
std::vector<int> CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::GetGlobalDofs(DRT::Node* node)
{
  // get dofs in beam contact discretization
  std::vector<int> cdofs = ContactDiscret().Dof(node);

  // get dofs in problem discretization via offsetmap
  std::vector<int> pdofs((int)(cdofs.size()));
  for (int k=0;k<(int)(cdofs.size());++k)
  {
    pdofs[k]=(dofoffsetmap_.find(cdofs[k]))->second;
  }

  return pdofs;
}

/*----------------------------------------------------------------------*
 |  Compute contact forces                                   grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::EvaluateFcContact(const double& pp,
                                                               const TYPE& gap,
                                                               const LINALG::TMatrix<TYPE,3,1>& normal,
                                                               Epetra_Vector& fint,
                                                               const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i,
                                                               const bool contactactive)
{
  fc1_.Clear();
  fc2_.Clear();

  // get dimensions for vectors fc1 and fc2
  const int dim1 = 3*numnodes*numnodalvalues;

  // temporary vectors for contact forces, DOF-GIDs and owning procs
  Epetra_SerialDenseVector fc1_copy(dim1);
  Epetra_SerialDenseVector fc2_copy(3);
  std::vector<int> lm1(dim1);
  std::vector<int> lm2(3);
  std::vector<int> lmowner1(dim1);
  std::vector<int> lmowner2(3);

  // flag indicating assembly
  bool DoNotAssemble = false;

  //**********************************************************************
  // evaluate contact forces for active pairs
  //**********************************************************************
  if (contactactive)
  {
    // node ids of both elements
    const int* node_ids1 = element1_->NodeIds();
    const int* node_ids2 = element2_->NodeIds();

    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    // prepare assembly
    for (int i=0;i<numnodes;++i)
    {
      // get node pointer and dof ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lm1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
        lmowner1[3*numnodalvalues*i+j] = node->Owner();
      }
    }

    // compute force vector Fc1
    for (int i=0;i<numnodes*numnodalvalues;++i)
    {
      for (int j=0;j<3;++j)
      {
        fc1_(3*i+j) = (- pp*gap) * normal(j) * N1_i(i);
      }
    }

    //********************************************************************
    // Compute Fc2 (force acting on second element)
    //********************************************************************
    // get node pointer and dof ids
    DRT::Node* node = ContactDiscret().gNode(node_ids2[0]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    // compute force vector Fc2 and prepare assembly
    for (int j=0;j<3;++j)
    {
      fc2_(j) = pp*gap * normal(j);
      lm2[j] = NodeDofGIDs[j];
      lmowner2[j] = node->Owner();
    }
  }
  //**********************************************************************
  // no forces for inactive pairs
  //**********************************************************************
  else
  {
    // set flag to avoid assembly
    DoNotAssemble = true;
  }

  //**********************************************************************
  // assemble contact forces
  //**********************************************************************
  if (!DoNotAssemble)
  {
    for (int i=0; i<dim1; ++i)
      fc1_copy[i] = FADUTILS::CastToDouble(fc1_(i));

    for (int i=0; i<3; ++i)
      fc2_copy[i] = FADUTILS::CastToDouble(fc2_(i));

    // assemble fc1 and fc2 into global contact force vector
    LINALG::Assemble(fint,fc1_copy,lm1,lmowner1);
    LINALG::Assemble(fint,fc2_copy,lm2,lmowner2);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate contact stiffness                               grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::EvaluateStiffcContact(const double& pp,
                                                            const TYPE& gap,
                                                            const LINALG::TMatrix<TYPE,3,1>& normal,
                                                            const TYPE& norm,
                                                            LINALG::SparseMatrix& stiffmatrix,
                                                            const LINALG::TMatrix<TYPE,3,1>& x1,
                                                            const LINALG::TMatrix<TYPE,3,1>& x2,
                                                            const LINALG::TMatrix<TYPE,3,1>& dx1,
                                                            const LINALG::TMatrix<TYPE,3,1>& ddx1,
                                                            const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i,
                                                            const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i_xi,
                                                            const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i_xixi,
                                                            bool activecontact,
                                                            bool linxi)
{
  // get dimensions for vector fc1 and fc2
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3;

  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
  LINALG::TMatrix<TYPE, dim1, dim1+dim2> stiffc1;
  LINALG::TMatrix<TYPE, dim2, dim1+dim2> stiffc2;
  Epetra_SerialDenseMatrix stiffc1_copy(dim1,dim1+dim2);
  Epetra_SerialDenseMatrix stiffc2_copy(dim2,dim1+dim2);
  std::vector<int> lmrow1(dim1);
  std::vector<int> lmrow2(dim2);
  std::vector<int> lmrowowner1(dim1);
  std::vector<int> lmrowowner2(dim2);
  std::vector<int> lmcol1(dim1+dim2);
  std::vector<int> lmcol2(dim1+dim2);

  // flag indicating assembly
  bool DoNotAssemble = false;

  // initialize stiffness to zero
  for (int i=0;i<dim1;i++)
    for (int j=0;j<(dim1+dim2);j++)
      stiffc1(i,j) = 0.0;
  for (int i=0;i<dim2;i++)
    for (int j=0;j<(dim1+dim2);j++)
      stiffc2(i,j) = 0.0;

  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  // initialize delta_xi here because we need it outside the if-environment for FAD Check
      LINALG::TMatrix<TYPE,dim1+dim2, 1> delta_xi;

  if (activecontact)
  {
    // auxiliary stiffmatrix for part III of linearization to avoid tensor notation
    LINALG::TMatrix<TYPE,dim1+dim2,dim1+dim2> stiffc_III;

    // node ids of both elements
    const int* node_ids1 = element1_->NodeIds();
    const int* node_ids2 = element2_->NodeIds();

    // initialize storage for linearizations
    LINALG::TMatrix<TYPE,3,1> distance;
    TYPE normdist = 0.0;
    LINALG::TMatrix<TYPE,dim1+dim2, 1> delta_gap;
    LINALG::TMatrix<TYPE,3, dim1+dim2> delta_x1_minus_x2;
    LINALG::TMatrix<TYPE,3, dim1+dim2> delta_n;

    //********************************************************************
    // evaluate linearizations and distance
    //********************************************************************
    // linearization of contact point
    // (not needed in case of check for nodal points in contact)
    if(linxi) ComputeLinXi(delta_xi,x1,x2,dx1,ddx1,N1_i,N1_i_xi);
    else
    {
      delta_xi.Clear();
    }

    // evaluation of distance
    ComputeDistance(distance, normdist, normal, norm);

    // linearization of gap function which is equal to delta d
    ComputeLinGap(delta_gap,delta_xi,x1,x2,dx1,N1_i,normdist,normal,norm,gap,delta_x1_minus_x2);

    // linearization of normal vector
    ComputeLinNormal(delta_n,delta_xi,normal,norm,dx1,N1_i);

    //********************************************************************
    // prepare assembly
    //********************************************************************
    // fill lmrow1 and lmrowowner1
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and dof ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmrow1[3*numnodalvalues*i+j]=NodeDofGIDs[j];
        lmrowowner1[3*numnodalvalues*i+j]=node->Owner();
      }
    }

    // fill lmrow2 and lmrowowner2
    // get pointer and node ids
    DRT::Node* node = ContactDiscret().gNode(node_ids2[0]);
    std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3;++j)
    {
      lmrow2[j]=NodeDofGIDs[j];
      lmrowowner2[j]=node->Owner();
    }

    // fill lmcol1 and lmcol2
    for (int i=0;i<numnodes;++i)
    {
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);

      for (int j=0;j<3*numnodalvalues;++j)
      {
        lmcol1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
        lmcol2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
      }
    }

    // fill lmcol1 and lmcol2
    // get pointer and node ids
    node = ContactDiscret().gNode(node_ids2[0]);
    NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3;++j)
    {
      lmcol1[3*numnodes*numnodalvalues+j] = NodeDofGIDs[j];
      lmcol2[3*numnodes*numnodalvalues+j] = NodeDofGIDs[j];
    }

    //********************************************************************
    // evaluate contact stiffness
    // (1) stiffc1 of first element
    //********************************************************************

    //********************************************************************
    // part I
    //********************************************************************
    LINALG::TMatrix<TYPE,dim1,1> N1T_normal;
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<numnodes*numnodalvalues;j++)
      {
        N1T_normal(3*j+i) += normal(i) * N1_i(j);
      }
    }

    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<(dim1+dim2);j++)
      {
        stiffc1(i,j) = -pp * N1T_normal(i) * delta_gap(j);
      }
    }

    //********************************************************************
    // part II
    //********************************************************************
    for  (int i=0;i<3;i++)
    {
      for (int j=0;j<numnodes*numnodalvalues;j++)
      {
        for (int k=0;k<dim1+dim2;k++)
        {
            stiffc1(3*j+i,k) -= pp*gap * N1_i(j) * delta_n(i,k);
        }
      }
    }

    //********************************************************************
    // part III
    //********************************************************************
    LINALG::TMatrix<TYPE,dim1,1> N1xiT_normal;
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<numnodes*numnodalvalues;j++)
      {
        N1xiT_normal(3*j+i) += normal(i) * N1_i_xi(j);
      }
    }

    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        stiffc1(i,j) -= pp*gap * N1xiT_normal(i) * delta_xi(j);
      }
    }

    //********************************************************************
    // evaluate contact stiffness
    // (2) stiffc2 of second element
    //********************************************************************

    //********************************************************************
    // part I
    //********************************************************************

    for (int i=0;i<3;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        stiffc2(i,j) = pp * normal(i) * delta_gap(j);
      }
    }

    //********************************************************************
    // part II
    //********************************************************************
    for (int i=0;i<3;i++)
    {
      for (int k=0;k<dim1+dim2;k++)
      {
        stiffc2(i,k) -= pp*gap * delta_n(i,k);
      }
    }
    //********************************************************************
    // part III
    //********************************************************************

    // There is no third part for spheres, since no shape functions are applied!
  }
  //**********************************************************************
  // no stiffness for inactive pairs
  //**********************************************************************
  else
  {
     // set flag to avoid assembly
    DoNotAssemble = true;
  }

  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  if (!DoNotAssemble)
  {
#ifndef AUTOMATICDIFF
    // change sign of stiffc1 and stiffc2 due to time integration.
    // according to analytical derivation there is no minus sign, but for
    // our time integration methods the negative stiffness must be assembled.
    for (int j=0;j<dim1+dim2;j++)
    {
      for (int i=0;i<dim1;i++)
        stiffc1_copy(i,j) = - FADUTILS::CastToDouble(stiffc1(i,j));
      for (int i=0;i<dim2;i++)
        stiffc2_copy(i,j) = - FADUTILS::CastToDouble(stiffc2(i,j));
    }
#else
    for (int j=0;j<dim1+dim2;j++)
    {
      for (int i=0;i<dim1;i++)
        stiffc1_copy(i,j) = - FADUTILS::CastToDouble(fc1_(i).dx(j) + fc1_(i).dx(dim1+dim2) * delta_xi(j) );
      for (int i=0;i<dim2;i++)
        stiffc2_copy(i,j) = - FADUTILS::CastToDouble(fc2_(i).dx(j) + fc2_(i).dx(dim1+dim2) * delta_xi(j) );
    }

#ifdef FADCHECKS
    std::cout << "BTSPH Contact Pair: " << element1_->Id() << " / " << element2_->Id() << std::endl;

    std::cout << "stiffc1: " << std::endl;
    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        std::cout << std::setw(14) << - stiffc1(i,j).val();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "stiffc1_FAD: " << std::endl;
    for (int i=0;i<dim1;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        std::cout << std::setw(14) << stiffc1_copy(i,j);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "stiffc2: " << std::endl;
    for (int i=0;i<dim2;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        std::cout << std::setw(14) << - stiffc2(i,j).val();
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    std::cout << "stiffc2_FAD: " << std::endl;
    for (int i=0;i<dim2;i++)
    {
      for (int j=0;j<dim1+dim2;j++)
      {
        std::cout << std::setw(14) << stiffc2_copy(i,j);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

#endif // FADCHECKS

#endif // AUTOMATICDIFF

   // now finally assemble stiffc1 and stiffc2
    stiffmatrix.Assemble(0,stiffc1_copy,lmrow1,lmrowowner1,lmcol1);
    stiffmatrix.Assemble(0,stiffc2_copy,lmrow2,lmrowowner2,lmcol2);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute normal vector in contact point                   grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeNormal(LINALG::TMatrix<TYPE,3,1>& normal,
                                                          TYPE& gap,
                                                          TYPE& norm,
                                                          const LINALG::TMatrix<TYPE,3,1>& x1,
                                                          const LINALG::TMatrix<TYPE,3,1>& x2)
{
  // compute non-unit normal
  for (int i=0;i<3;i++)
    normal(i) = x1(i)-x2(i);

  // compute length of normal
  norm = FADUTILS::VectorNorm<3>(normal);
  if (norm < NORMTOL) dserror("ERROR: Normal of length zero! --> change time step!");

  // compute unit normal and store it in class variable
  for (int i=0;i<3;i++)
    {
    normal(i) /= norm;
    normal_(i)=normal(i);
    }

  // evaluate scalar gap function
  ComputeGap(gap,norm);

  return;
}

/*----------------------------------------------------------------------*
 |  Compute scalar gap function                               grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeGap(TYPE& gap,
                                                                         const TYPE& norm)
{
  const DRT::ELEMENTS::Beam3Base* beamele = static_cast<const DRT::ELEMENTS::Beam3Base*>(element1_);

  if (beamele == NULL)
    dserror("cast to beam base failed!");

  // compute radii of both elements
  double radius_ele1 = MANIPULATERADIUS * beamele->GetCircularCrossSectionRadiusForInteractions();

  double radius_ele2 = (static_cast<DRT::ELEMENTS::Rigidsphere*>(element2_))->Radius();

  // compute gap to be returned
  gap = norm - radius_ele1 - radius_ele2;
  gap_ = gap;

  return;
}


/*----------------------------------------------------------------------*
 | compute contact point coordinates and their derivatives    grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeCoordsAndDerivs(LINALG::TMatrix<TYPE,3,1>& x1,
                                                             LINALG::TMatrix<TYPE,3,1>& x2,
                                                             LINALG::TMatrix<TYPE,3,1>& dx1,
                                                             LINALG::TMatrix<TYPE,3,1>& ddx1,
                                                             const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i,
                                                             const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i_xi,
                                                             const LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i_xixi)
{
  // reset input variables
  x1.Clear();
  x2.Clear();
  dx1.Clear();
  ddx1.Clear();

#ifdef AUTOMATICDIFF
  // The 3*numnodes*numnodalvalues+3 primary DoFs are the components of the nodal positions / tangents. The additional
  // degree of freedom (+1) represents the dependency on the parameter coordinates xi, which is necessary in beam contact.
  for (int i=0;i<3*numnodes*numnodalvalues;i++)
    ele1pos_(i).diff(i,3*numnodes*numnodalvalues+3+1);

  for (int i=0;i<3;i++)
    ele2pos_(i).diff(3*numnodes*numnodalvalues+i,3*numnodes*numnodalvalues+3+1);
#endif // AUTOMATICDIFF

  // compute output variables
  for (int i=0;i<3;i++)
      x2(i)=ele2pos_(i);

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<numnodes*numnodalvalues;j++)
    {
      x1(i) += N1_i(j) * ele1pos_(3*j+i);              // x1 = N1 * x~1
      dx1(i) += N1_i_xi(j) * ele1pos_(3*j+i);          // dx1 = N1,xi * x~1
      ddx1(i) += N1_i_xixi(j) * ele1pos_(3*j+i);         // ddx1 = N1,xixi * x~1
    }
  }

  // store coordinates of contact point into class variables
  for(int i=0;i<3;i++)
  {
    x1_(i)=x1(i);
    x2_(i)=x2(i);
  }

  return;
}

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::GetShapeFunctions(LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i,
                                                      LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i_xi,
                                                      LINALG::TMatrix<TYPE,1,numnodes*numnodalvalues>& N1_i_xixi,
                                                      const TYPE& eta)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype1 = element1_->Shape();

  if (numnodalvalues==1)
  {
    // get values and derivatives of shape functions
    DRT::UTILS::shape_function_1D(N1_i, eta, distype1);
    DRT::UTILS::shape_function_1D_deriv1(N1_i_xi, eta, distype1);
    DRT::UTILS::shape_function_1D_deriv2(N1_i_xixi, eta, distype1);
  }
  else if (numnodalvalues==2)
  {
    if ( element1_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
       dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");

     double length1 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->jacobi();

     // get values and derivatives of shape functions
     DRT::UTILS::shape_function_hermite_1D(N1_i, eta, length1, distype1);
     DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi, eta, length1, distype1);
     DRT::UTILS::shape_function_hermite_1D_deriv2(N1_i_xixi, eta, length1, distype1);
  }
  else
    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");

  return;
}

/*----------------------------------------------------------------------*
 |  Linearizations of contact point                          grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeLinXi(LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues+3, 1>& delta_xi,
                                                          const LINALG::TMatrix<TYPE,3,1>& x1,
                                                          const LINALG::TMatrix<TYPE,3,1>& x2,
                                                          const LINALG::TMatrix<TYPE,3,1>& dx1,
                                                          const LINALG::TMatrix<TYPE,3,1>& ddx1,
                                                          const LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues>& N1_i,
                                                          const LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues>& N1_i_xi)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3;

  // matrices to compute Lin_Xi (in case of beam-sphere contact: scalar and vector)
  TYPE L=0.0;
  LINALG::TMatrix<TYPE,1,dim1+dim2> B;

  LINALG::TMatrix<TYPE,3,1> delta_x;

  // compute L
  for (int i=0;i<3;i++)
  {
    delta_x(i) = x1(i) - x2(i);
    L +=  dx1(i)*dx1(i) + delta_x(i)*ddx1(i);
  }

  if (L == 0) dserror("ERROR: L = 0");


  for (int i=0; i<3; ++i)
  {
    for (int j=0; j<numnodes*numnodalvalues; ++j)
    {
      B(3*j+i) += -delta_x(i)*N1_i_xi(j) - dx1(i) * N1_i(j);
    }
  }

  for (int i=0; i<dim2; ++i)
    B(dim1+i) += dx1(i);

  // finally the linearizations / directional derivatives
  for (int i=0;i<dim1+dim2;i++)
  {
    delta_xi(i) = B(i)/L;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  Compute distance vector                                  grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeDistance(LINALG::TMatrix<TYPE,3,1>& distance,
                                                                 TYPE& normdist,
                                                                 const LINALG::TMatrix<TYPE,3,1>& normal,
                                                                 const TYPE& norm)
{
  // compute distance vector
  for (int i=0;i<3;i++)
    distance(i) = normal(i) * norm;

  // compute scalar distance
  normdist = norm;

  return;
}

/*----------------------------------------------------------------------*
 | Compute linearization of gap                              grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeLinGap(LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues+3, 1>& delta_gap,
                                                LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues+3, 1>& delta_xi,
                                                const LINALG::TMatrix<TYPE,3,1>& x1,
                                                const LINALG::TMatrix<TYPE,3,1>& x2,
                                                const LINALG::TMatrix<TYPE,3,1>& dx1,
                                                const LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues>& N1_i,
                                                const TYPE& normdist,
                                                const LINALG::TMatrix<TYPE,3,1>& normal,
                                                const TYPE& norm,
                                                const TYPE& gap,
                                                LINALG::TMatrix<TYPE,3, 3*numnodes*numnodalvalues+3>& delta_x1_minus_x2)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3;


  // delta g := delta_r/||delta_r||*auxiliary_matri1 delta d, with auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  LINALG::TMatrix<TYPE,3,dim1+dim2>  auxiliary_matrix1(true);

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      auxiliary_matrix1(i,j)+=dx1(i)*delta_xi(j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<numnodes*numnodalvalues;j++)
    {
      auxiliary_matrix1(i,3*j+i)+= N1_i(j);
    }
  }

  for (int i=0;i<3;i++)
  {
      auxiliary_matrix1(i,i+dim1)+= -1.0;
  }

  // compute linearization of gap
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      delta_gap(j) +=  (x1(i)-x2(i)) * auxiliary_matrix1(i,j)/normdist;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | Compute linearization of normal vector                    grill 09/14|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::ComputeLinNormal(LINALG::TMatrix<TYPE,3, 3*numnodes*numnodalvalues+3>& delta_normal,
                                                            const LINALG::TMatrix<TYPE,3*numnodes*numnodalvalues+3, 1>&  delta_xi,
                                                            const LINALG::TMatrix<TYPE,3,1>& normal,
                                                            const TYPE& norm_delta_x,
                                                            const LINALG::TMatrix<TYPE,3, 1>& x1_xi,
                                                            const LINALG::TMatrix<TYPE,1, numnodes*numnodalvalues>& N1_i)
{
  const int dim1 = 3*numnodes*numnodalvalues;
  const int dim2 = 3;

  //delta n := auxiliary_matri2*auxiliary_matrix1* delta d, with auxiliary_matri2 = (I-nxn)/||r1-r2||
  //and auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))

  LINALG::TMatrix<TYPE,3,dim1+dim2>  auxiliary_matrix1(true);
  LINALG::TMatrix<TYPE,3,3>  auxiliary_matrix2(true);

  //compute auxiliary_matrix1
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<dim1+dim2;j++)
    {
      auxiliary_matrix1(i,j)+=x1_xi(i)*delta_xi(j);
    }
  }

  for (int i=0;i<3;i++)
  {
    for (int j=0;j<numnodes*numnodalvalues;j++)
    {
      auxiliary_matrix1(i,3*j+i)+= N1_i(j);
    }
  }

  for (int i=0;i<3;i++)
  {
      auxiliary_matrix1(i,dim1+i)+= -1.0;
  }

  //compute auxiliary_matrix2
  for (int i=0;i<3;i++)
  {
    auxiliary_matrix2(i,i)+= 1.0/norm_delta_x;
    for (int j=0;j<3;j++)
    {
      auxiliary_matrix2(i,j)+= -normal(i)*normal(j)/norm_delta_x;
    }
  }

  // compute linearization of normal vector
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<dim1+dim2;k++)
        delta_normal(i,k) +=  auxiliary_matrix2(i,j) * auxiliary_matrix1(j,k);

  return;
}

/*----------------------------------------------------------------------*
 |  Update nodal coordinates (public)                           popp 04/10|
 *----------------------------------------------------------------------*/
template<const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                                         Epetra_SerialDenseMatrix& newele2pos)
{
  for (int i=0;i<3*numnodalvalues;i++)
  {
    for (int j=0;j<numnodes;j++)
    {
      ele1pos_(3*numnodalvalues*j+i)=newele1pos(i,j);
    }
  }

  for (int i=0; i<3; ++i)
    ele2pos_(i)=newele2pos(i,0);

  return;
}

template<const int numnodes, const int numnodalvalues>
void CONTACT::Beam3tospherecontact<numnodes, numnodalvalues>::Print() const
{
  // ToDo add further information here
  printf("\nInstance of Beam3tospherecontact: element GIDs %i (beam) and %i (sphere)",element1_->Id(),element2_->Id());
  std::cout << "\ncontactflag_: " << contactflag_ << "\tnodalcontactflag: " << nodalcontactflag_[0] << nodalcontactflag_[1];
  std::cout << "\ngap_: " << gap_;
  std::cout << "\nxicontact_: " << xicontact_;
  std::cout << "\nnormal_: " << normal_;
  std::cout << "\nx1_: " << x1_;
  std::cout << "\nx2_: " << x2_;
  std::cout << "\nele1pos_: " << ele1pos_;
  std::cout << "\nele2pos_: " << ele2pos_;

  return;
}


Teuchos::RCP<CONTACT::Beam3tospherecontactinterface> CONTACT::Beam3tospherecontactinterface::Impl( const int numnodes,
                                                                      const int numnodalvalues,
                                                                      const DRT::Discretization& pdiscret,
                                                                      const DRT::Discretization& cdiscret,
                                                                      const std::map<int,int>& dofoffsetmap,
                                                                      DRT::Element* element1,
                                                                      DRT::Element* element2)
{

  switch (numnodalvalues)
  {
    case 1:
    {
      switch (numnodes)
      {
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tospherecontact<2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2));
        }
        case 3:
        {
          return Teuchos::rcp (new CONTACT::Beam3tospherecontact<3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2));
        }
        case 4:
        {
          return Teuchos::rcp (new CONTACT::Beam3tospherecontact<4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2));
        }
        case 5:
        {
          return Teuchos::rcp (new CONTACT::Beam3tospherecontact<5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2));
        }
        default:
          dserror("No valid template parameter for the number of nodes (numnodes = 2,3,4,5 for Reissner beams) available!");
          break;
      }
      break;
    }
    case 2:
    {
      switch (numnodes)
      {
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tospherecontact<2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2));
        }
        default:
          dserror("No valid template parameter for the number of nodes (only numnodes = 2 for Kirchhoff beams valid so far) available!");
          break;
      }
      break;
    }
    default:
      dserror("No valid template parameter for the Number of nodal values (numnodalvalues = 1 for Reissner beams, numnodalvalues = 2 for Kirchhoff beams) available!");
      break;
  }
  return Teuchos::null;
}
