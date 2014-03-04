/*!-----------------------------------------------------------------------------------------------------------
\file Beam3tospherecontact.cpp
\brief One beam contact pair (two beam elements)

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "beam3tospherecontact.H"
#include "beam3contact_defines.H"
#include "../drt_inpar/inpar_beamcontact.H"
#include "../drt_inpar/inpar_contact.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_exporter.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_beam3/beam3.H"
#include "../drt_beam3ii/beam3ii.H"
#include "../drt_rigidsphere/rigidsphere.H"
#include "../drt_inpar/inpar_statmech.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                      popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3tospherecontact::Beam3tospherecontact(const DRT::Discretization& pdiscret,
                                    const DRT::Discretization& cdiscret,
                                    const int& dofoffset,
                                    DRT::Element* element1,
                                    DRT::Element* element2):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffset_(dofoffset),
element1_(element1),
element2_(element2),
contactflag_(false)
{

  // initialize augmented lagrange multiplier to zero
  lmuzawa_ = 0.0;
  gap_=0.0;

  ele1pos_.Reshape(3,element1->NumNode());
  ele2pos_.Reshape(3,1);


  // initialize class variables for contact point coordinates
  x1_.Size(3);
  x2_.Size(3);
  normal_.Size(3);
  for(int i=0;i<3;i++)
  {
    x1_[i]=0.0;
    x2_[i]=0.0;
    normal_[i]=0.0;
    firstcall_ = true;

    const DRT::ElementType & eot1 = element1_->ElementType();
    const DRT::ElementType & eot2 = element2_->ElementType();

    if(eot1 != DRT::ELEMENTS::Beam3Type::Instance() and eot1 != DRT::ELEMENTS::Beam3iiType::Instance())
      dserror("How did you get here? element1_ has to be of type beam3 or beam3ii!!!");

    if(eot2 != DRT::ELEMENTS::RigidsphereType::Instance())
      dserror("How did you get here? element2_ has to be of type rigidsphere!!!");

  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: constructor
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  copy-constructor (public)                                 popp 04/10|
 *----------------------------------------------------------------------*/
CONTACT::Beam3tospherecontact::Beam3tospherecontact(const Beam3tospherecontact& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffset_(old.dofoffset_),
element1_(old.element1_),
element2_(old.element2_),
ele1pos_(old.ele1pos_),
ele2pos_(old.ele2pos_)
{
  dserror("ERROR: Copy constructor incomplete");
  return;
}
/*----------------------------------------------------------------------*
 |  end: copy-constructor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                             popp 04/10|
 *----------------------------------------------------------------------*/
bool CONTACT::Beam3tospherecontact::Evaluate(LINALG::SparseMatrix& stiffmatrix, Epetra_Vector& fint, double& pp)
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

  // initialize the two element coordinates xi and eta of contact point
  double XiContact = xicontact_;

  // number of nodes of each element
  const int numnode1 = element1_->NumNode();

  // vectors for shape functions and their derivatives
  Epetra_SerialDenseVector funct1(numnode1);          // = N1
  Epetra_SerialDenseMatrix deriv1(1,numnode1);        // = N1,xi
  Epetra_SerialDenseMatrix secondderiv1(1,numnode1);  // = N1,xixi

  // coords and derivatives of the two contact points
  std::vector<double> x1(3);                            // = x1
  std::vector<double> x2(3);                            // = x2
  std::vector<double> dx1(3);                            // = x1,xi
  std::vector<double> ddx1(3);                          // = x1,xixi

  // initialize
  std::vector<double> normal(3);
  double gap= 0.0;
  double norm = 0.0;

  if (abs(XiContact)< (1.0 + XIETATOL))
  {
    std::cout << "Auswertung von Paar:" << element1_->Id() << "/" << element2_->Id() << std::endl;
  }
  else
  {
    contactflag_ = false;
    return false;
  }

  //**********************************************************************
  // (1) Compute some auxiliary quantities
  //**********************************************************************

  // call function to fill variables for shape functions and their derivatives
  GetShapeFunctions(funct1,deriv1,secondderiv1,XiContact);

  // call function to fill variables with coords and derivs of the contact point
  ComputeCoordsAndDerivs(x1,x2,dx1,ddx1,funct1,deriv1,secondderiv1,numnode1);

  // call function to compute scaled normal and gap in possible contact point
  ComputeNormal(normal,gap,norm,x1,x2);

  // call function to evaluate contact status
  CheckContactStatus(pp);
  
  // store coordinates of contact point into class variables
  for(int i=0;i<3;i++)
  {
    x1_[i]=x1[i];
    x2_[i]=x2[i];
  }

  //**********************************************************************
  // (2) Compute contact forces and stiffness
  //**********************************************************************

  // call function to evaluate and assemble contact forces
  EvaluateFcContact(pp,gap,normal,fint,funct1,numnode1);
  // call function to evaluate and assemble contact stiffness
  EvaluateStiffcContact(pp,gap,normal,norm,stiffmatrix,x1,x2,dx1,
                        ddx1,funct1,deriv1,secondderiv1,numnode1,XiContact);
  return true;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate the element
  *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  Closest point projection                                  popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ClosestPointProjection()
{


  // local variables for element coordinates
  double eta=0.0;

  // number of nodes of each element
  const int numnode1 = Element1()->NumNode();

  // vectors for shape functions and their derivatives
  Epetra_SerialDenseVector funct1(numnode1);          // = N1
  Epetra_SerialDenseMatrix deriv1(1,numnode1);        // = N1,xi
  Epetra_SerialDenseMatrix secondderiv1(1,numnode1);  // = N1,xixi

  // coords and derivatives of the two contacting points
  std::vector<double> x1(3);                            // = x1
  std::vector<double> x2(3);                            // = x2
  std::vector<double> dx1(3);                           // = x1,xi
  std::vector<double> ddx1(3);                          // = x1,xixi

  // initialize function f and Jacobian df for Newton iteration
  double f;
  double df;

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
    for (int j=0;j<numnode1;j++) funct1[j]=0;
    for (int j=0;j<numnode1;j++) deriv1(0,j)=0;
    for (int j=0;j<numnode1;j++) secondderiv1(0,j)=0;

    // update shape functions and their derivatives
    GetShapeFunctions(funct1,deriv1,secondderiv1,eta);

    // update coordinates and derivatives of contact points
    ComputeCoordsAndDerivs(x1,x2,dx1,ddx1,funct1,deriv1,secondderiv1,numnode1);

    // compute norm of difference vector to scale the equations
    // (this yields better conditioning)
    double norm = sqrt((x1[0]-x2[0])*(x1[0]-x2[0])
                     + (x1[1]-x2[1])*(x1[1]-x2[1])
                     + (x1[2]-x2[2])*(x1[2]-x2[2]));


    // the closer the beams get, the smaller is norm
    // norm is not allowed to be too small, else numerical problems occur
    if (norm < NORMTOL)
    {
      dserror("Contact points x1 and x2 are identical. Choose smaller time step!");
    }

    // evaluate f at current eta
    EvaluateNewtonF(f,x1,x2,dx1,norm);

    // compute the scalar residuum
    residual = abs(f);

    // check if Newton iteration has converged
    if (residual < BEAMCONTACTTOL) break;

    // evaluate Jacobian of f at current eta
    EvaluateNewtonGradF(df,x1,x2,dx1,ddx1,norm);

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

  // store and return final result
  xicontact_=eta;

  return;
}
/*----------------------------------------------------------------------*
|  end: Closest point projection
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate function f in CPP                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::EvaluateNewtonF(double& f, const std::vector<double>& x1,
     const std::vector<double>& x2,  const std::vector<double>& dx1, const double& norm)
{
  // reset f
  f=0;
  
  // evaluate f
  // see Wriggers, Computational Contact Mechanics, equation (12.5)

  for (int i=0;i<3;i++)
  {
    f += (x1[i]-x2[i])*dx1[i] / norm;
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate function f in CPP
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df in CPP                               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::EvaluateNewtonGradF( double& df, const std::vector<double>& x1,
                                                         const std::vector<double>& x2, const std::vector<double>& dx1,
                                                         const std::vector<double>& ddx1,const double& norm)

{
  // reset df
  df=0;
  
  // evaluate df
  // see Wriggers, Computational Contact Mechanics, equation (12.7)
  for(int i=0;i<3;i++)
  {
    df += (dx1[i]*dx1[i] + (x1[i]-x2[i])*ddx1[i]) / norm;
  }

  return;  
}
/*----------------------------------------------------------------------*
 |  end: Evaluate Jacobian df in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Check if conact is active or inactive                     popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::CheckContactStatus(double& pp)
{
  // check contact condition valid for both pure penalty
  // (lmuzawa = 0) and augmented lagrange (lmuzawa != 0)
  if (lmuzawa_ - pp * gap_ > 0) contactflag_ = true;
  else                          contactflag_ = false;

  return;
}
/*----------------------------------------------------------------------*
 |  end: Check if conact is active or inactive
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  Get global dofs of a node                                 popp 12/10|
 *----------------------------------------------------------------------*/
std::vector<int> CONTACT::Beam3tospherecontact::GetGlobalDofs(DRT::Node* node)
{
  // get dofs in beam contact discretization
  std::vector<int> cdofs = ContactDiscret().Dof(node);

  // get dofs in problem discretization via offset
  std::vector<int> pdofs((int)(cdofs.size()));
  for (int k=0;k<(int)(cdofs.size());++k)
    pdofs[k] = cdofs[k]-DofOffset();

  return pdofs;
}
/*----------------------------------------------------------------------*
 |  end: Get global dofs of a node
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute contact forces                                    popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::EvaluateFcContact(const double& pp,
     const double& gap, const std::vector<double>& normal, Epetra_Vector& fint,
     Epetra_SerialDenseVector funct1, const int numnode1)
{
  // get dimensions for vectors fc1 and fc2
  const int dim1 = 3*numnode1;
  
  // temporary vectors for contact forces, DOF-GIDs and owning procs
  Epetra_SerialDenseVector fc1(dim1);
  Epetra_SerialDenseVector fc2(3);
  std::vector<int> lm1(dim1);
  std::vector<int> lm2(3);
  std::vector<int> lmowner1(dim1);
  std::vector<int> lmowner2(3);

  // flag indicating assembly
  bool DoNotAssemble = false;
  
  //**********************************************************************
  // evaluate contact forces for active pairs
  //**********************************************************************
  if (contactflag_)
  {
    // node ids of both elements
    const int* node_ids1 = element1_->NodeIds();
    const int* node_ids2 = element2_->NodeIds();

    //********************************************************************
    // Compute Fc1 (force acting on first element)
    //********************************************************************
    for (int i=0;i<numnode1;++i)
    {
      // get node pointer and dof ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs = GetGlobalDofs(node);

      // compute force vector Fc1 and prepare assembly
      for (int j=0;j<3;++j)
      {
        fc1[3*i+j] = (lmuzawa_ - pp*gap) * normal[j] * funct1[i];
        lm1[3*i+j] = NodeDofGIDs[j];
        lmowner1[3*i+j] = node->Owner();
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
        fc2[j] = -(lmuzawa_ - pp*gap) * normal[j];
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
    
    // compute forces
    for (int i=0;i<3*numnode1;++i) fc1[i]=0;
    for (int i=0;i<3;++i) fc2[i]=0;
  }
  
  //**********************************************************************
  // assemble contact forces
  //**********************************************************************
  if (!DoNotAssemble)
  {  
    // assemble fc1 and fc2 into global contact force vector
    LINALG::Assemble(fint,fc1,lm1,lmowner1);
    LINALG::Assemble(fint,fc2,lm2,lmowner2);
  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute contact forces
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate contact stiffness                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::EvaluateStiffcContact(const double& pp,
     const double& gap, const std::vector<double>& normal, const double& norm, 
     LINALG::SparseMatrix& stiffmatrix,const std::vector<double>& x1,const std::vector<double>& x2,
     const std::vector<double>& dx1, const std::vector<double>& ddx1,
     const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseMatrix& deriv1,
     const Epetra_SerialDenseMatrix& secondderiv1, const int numnode1, double& XiContact)
{
  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs 
  Epetra_SerialDenseMatrix stiffc1(3*numnode1,3*(numnode1+1));
  Epetra_SerialDenseMatrix stiffc2(3,3*(numnode1+1));
  std::vector<int> lmrow1(3*numnode1);
  std::vector<int> lmrow2(3);
  std::vector<int> lmrowowner1(3*numnode1);
  std::vector<int> lmrowowner2(3);
  std::vector<int> lmcol1(3*(numnode1+1));
  std::vector<int> lmcol2(3*(numnode1+1));
  
  // flag indicating assembly
  bool DoNotAssemble = false;
  
  //**********************************************************************
  // evaluate contact stiffness for active pairs
  //**********************************************************************
  if (contactflag_)
  {

    // auxiliary stiffmatrix for part III of linearization to avoid tensor notation
    Epetra_SerialDenseMatrix stiffc_III(3*(numnode1+1),3*(numnode1+1));
    
    // node ids of both elements
    const int* node_ids1 = element1_->NodeIds();
    const int* node_ids2 = element2_->NodeIds();
    
    // initialize storage for linearizations
    std::vector<double> delta_xi(3*(numnode1+1));
    std::vector<double> distance(3);
    double normdist = 0.0;
    std::vector<double> delta_gap(3*(numnode1+1));
    Epetra_SerialDenseMatrix delta_x1_minus_x2(3,3*(numnode1+1));
    Epetra_SerialDenseMatrix delta_n(3, 3*(numnode1+1));
    
    //********************************************************************
    // evaluate linearizations and distance
    //********************************************************************
    // linearization of contact point
    ComputeLinXiAndLinEta(delta_xi,x1,x2,dx1,ddx1,funct1,deriv1,normal,norm,numnode1,XiContact);
    
    // evaluation of distance
    ComputeDistance(distance, normdist, normal, norm);
    
    // linearization of gap function which is equal to delta d
    ComputeLinGap(delta_gap,delta_xi,x1,x2,dx1,funct1,normdist,numnode1,normal,norm,gap,delta_x1_minus_x2);

    // linearization of normal vector
    ComputeLinNormal(delta_n,x1,x2,norm,numnode1,delta_x1_minus_x2,normal,XiContact);

    //********************************************************************
    // prepare assembly
    //********************************************************************
    // fill lmrow1 and lmrowowner1
    for (int i=0;i<numnode1;++i)
    {
      // get pointer and dof ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);
      
      for (int j=0;j<3;++j)
      {
        lmrow1[3*i+j]=NodeDofGIDs[j];
        lmrowowner1[3*i+j]=node->Owner();
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
    for (int i=0;i<numnode1;++i)
    {  
      // get pointer and node ids
      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);
      
      for (int j=0;j<3;++j)
      {
        lmcol1[3*i+j] = NodeDofGIDs[j];
        lmcol2[3*i+j] = NodeDofGIDs[j];
      }
    }
    
    // fill lmcol1 and lmcol2
    // get pointer and node ids
    node = ContactDiscret().gNode(node_ids2[0]);
    NodeDofGIDs =  GetGlobalDofs(node);

    for (int j=0;j<3;++j)
    {
      lmcol1[3*numnode1+j] = NodeDofGIDs[j];
      lmcol2[3*numnode1+j] = NodeDofGIDs[j];
    }
        
    //********************************************************************
    // index vectors for access to shape function vectors and matrices
    //********************************************************************
    // In literature shape functions are stored in matrices but here only
    // a vector is used. To access an element of a shape function matrix,
    // use index vectors, that relate the vector index to the row index of
    // the shape function matrices. Use an if-contruction to select only
    // the entries of the shape function matrices, that are != 0.
    // --> if (j%3 == i): only quasi-diagonal entries are considered
    //********************************************************************
    // intialize index vectors
    std::vector<int> index1(3*(numnode1+1));
    std::vector<int> index2(3*(numnode1+1));
    
    // fill the index vectors
    for (int i=0;i<3*numnode1;++i)
    {
      index1[i]=(int)floor(i/3);
      index2[i]=(int)floor(i/3);
    }
    for (int i=3*numnode1;i<3*(numnode1+1);++i)
    {
      index1[i]=(int)floor((i-3*numnode1)/3);
      index2[i]=(int)floor((i-3*numnode1)/3);
    }
    
    //********************************************************************
    // evaluate contact stiffness
    // (1) stiffc1 of first element
    //********************************************************************
    
    //********************************************************************
    // part I
    //********************************************************************
    std::vector<double> normal_t_N1(3*numnode1);
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnode1;j++)
      {
        if (j%3 == i)
        {
          int row = index1[j];
          normal_t_N1[j] += normal[i] * funct1[row];
        }
      }
    }

    for (int i=0;i<3*numnode1;i++)
    {
      for (int j=0;j<3*(numnode1+1);j++)
      {
        stiffc1(i,j) = -pp * normal_t_N1[i] * delta_gap[j];
      }
    }

    //********************************************************************
    // part II
    //********************************************************************
    for  (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnode1;j++)
      {
        for (int k=0;k<3*(numnode1+1);k++)
        {
          if (j%3 == i)
          {
            int row = index1[j];
            stiffc1(j,k) += (lmuzawa_ - pp*gap) * funct1[row] * delta_n(i,k);
          }
        }
      }
    }

    //********************************************************************
    // part III
    //********************************************************************
    std::vector<double> normal_t_N1_xi(3*numnode1);
    for (int i=0;i<3;i++)
    {
      for (int j=0;j<3*numnode1;j++)
      {
        if (j%3 == i)
        {
          int row = index1[j];
          normal_t_N1_xi[j] += normal[i] * deriv1(0,row);
        }
      }
    }

    for (int i=0;i<3*numnode1;i++)
    {
      for (int j=0;j<3*(numnode1+1);j++)
      {
        stiffc1(i,j) += (lmuzawa_ - pp*gap) * normal_t_N1_xi[i] * delta_xi[j];
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
      for (int j=0;j<3*(numnode1+1);j++)
      {
        stiffc2(i,j) = pp * normal[i] * delta_gap[j];
      }
    }

    //********************************************************************
    // part II
    //********************************************************************
    for (int i=0;i<3;i++)
    {
      for (int k=0;k<3*(numnode1+1);k++)
      {
        stiffc2(i,k) += -(lmuzawa_  - pp*gap) * delta_n(i,k);
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
    
    // compute stiffness
    for (int i=0;i<3*numnode1;i++)
      for (int j=0;j<3*(numnode1+1);j++)
        stiffc1(i,j) = 0;
    for (int i=0;i<3;i++)
      for (int j=0;j<3*(numnode1+1);j++)
        stiffc2(i,j) = 0;
  }
  
  //**********************************************************************
  // assemble contact stiffness
  //**********************************************************************
  // change sign of stiffc1 and stiffc2 due to time integration.
  // according to analytical derivation there is no minus sign, but for
  // our time integration methods the negative stiffness must be assembled.
  for (int j=0;j<3*(numnode1+1);j++)
  {
    for (int i=0;i<3*numnode1;i++)
      stiffc1(i,j) = stiffc1(i,j) * (-1);
    for (int i=0;i<3;i++)
      stiffc2(i,j) = stiffc2(i,j) * (-1);
  }

  // now finally assemble stiffc1 and stiffc2
  if (!DoNotAssemble)
  {
    stiffmatrix.Assemble(0,stiffc1,lmrow1,lmrowowner1,lmcol1);
    stiffmatrix.Assemble(0,stiffc2,lmrow2,lmrowowner2,lmcol2);
  }
  
  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate contact stiffness
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute normal vector in contact point                    popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeNormal(std::vector<double>& normal, double& gap,
    double& norm, const std::vector<double>& x1, const std::vector<double>& x2)
{
  // compute non-unit normal
  for (int i=0;i<3;i++) normal[i] = x1[i]-x2[i];

  // compute length of normal
  norm = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
  if (norm < NORMTOL) dserror("ERROR: Normal of length zero! --> change time step!");
  
  // compute unit normal
  for (int i=0;i<3;i++)
    {
    normal[i] /= norm;
    normal_[i]=normal[i];
    }

  // evaluate scalar gap function
  ComputeGap(gap,norm);
  
  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute normal vector in contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute scalar gap function                                popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeGap(double& gap, const double& norm)
{
  // get moments of inertia of both elements
  // NOTE: here Iyy_ = Izz_ due to the circular cross section
  double MomentOfInertia_ele1 = 0;
  double radius_ele2 = 0;

  const DRT::ElementType & eot1 = element1_->ElementType();

  if ( eot1 == DRT::ELEMENTS::Beam3Type::Instance() )
    MomentOfInertia_ele1 = (static_cast<DRT::ELEMENTS::Beam3*>(element1_))->Iyy();

  if ( eot1 == DRT::ELEMENTS::Beam3iiType::Instance() )
    MomentOfInertia_ele1 = (static_cast<DRT::ELEMENTS::Beam3ii*>(element1_))->Iyy();
  
  radius_ele2 = (static_cast<DRT::ELEMENTS::Rigidsphere*>(element2_))->Radius();

  // compute radii of both elements
  double radius_ele1=0;
  ComputeEleRadius(radius_ele1,MomentOfInertia_ele1);
  
  // comute gap to be returned
  gap = norm - radius_ele1 - radius_ele2;
  gap_ = gap;

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute scalar gap function
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  compute radius from moment of inertia                     popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeEleRadius(double& radius, const double& moi)
{
  // fixed formula for circular cross sections: the factor f can be used to manipulate the geometrical radius if not equal to 1
  radius = MANIPULATERADIUS * sqrt(sqrt(4 * moi / M_PI));
  
  return;
}

/*----------------------------------------------------------------------*
 | compute contact point coordinates and their derivatives     popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeCoordsAndDerivs(std::vector<double>& x1,
     std::vector<double>& x2, std::vector<double>& dx1,std::vector<double>& ddx1,
     const Epetra_SerialDenseVector& funct1,const Epetra_SerialDenseMatrix& deriv1,
     const Epetra_SerialDenseMatrix& secondderiv1,const int& numnode1)
{
  // reset input variables
  for (int i=0;i<3;i++)
  {
    x1[i]   = 0.0;
    x2[i]   = 0.0;
    dx1[i]  = 0.0;
    ddx1[i] = 0.0;
  }
    
  // auxialiary variables for the nodal coordinates
  Epetra_SerialDenseMatrix coord1(3,numnode1);
  
  // full coord1 and coord2
  for (int i=0;i<3;i++)
    for (int j=0;j<numnode1;j++)      
      coord1(i,j)=ele1pos_(i,j);      
  for (int i=0;i<3;i++)
      x2[i]=ele2pos_(i,0);
  
  // compute output variable
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<numnode1;j++)
      x1[i] += funct1[j] * coord1(i,j);              // x1 = N1 * x~1
    for (int j=0;j<numnode1;j++)
      dx1[i] += deriv1(0,j) * coord1(i,j);          // dx1 = N1,xi * x~1
    for (int j=0;j<numnode1;j++)
      ddx1[i] += secondderiv1(0,j) * coord1(i,j);    // ddx1 = N1,xixi * x~1
  }
  
  return;
}
/*----------------------------------------------------------------------*
 | end: compute contact point coordinates and their derivatives         |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                  popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::GetShapeFunctions(Epetra_SerialDenseVector& funct1,
                                                      Epetra_SerialDenseMatrix& deriv1,
                                                      Epetra_SerialDenseMatrix& secondderiv1,
                                                      const double& eta)
{
  // get both discretization types
  const DRT::Element::DiscretizationType distype1 = element1_->Shape();
      
  // get values and derivatives of shape functions
  DRT::UTILS::shape_function_1D(funct1, eta, distype1);
  DRT::UTILS::shape_function_1D_deriv1(deriv1, eta, distype1);
  DRT::UTILS::shape_function_1D_deriv2(secondderiv1, eta, distype1);
    
  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Linearizations of contact point                           popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeLinXiAndLinEta(
    std::vector<double>& delta_xi, const std::vector<double>& x1, const std::vector<double>& x2,
    const std::vector<double>& dx1, const std::vector<double>& ddx1,
    const Epetra_SerialDenseVector& funct1, const Epetra_SerialDenseMatrix& deriv1,
    const std::vector<double>& normal, const double& norm, const int numnode1, const double& XiContact)
{
  
  // matrices to compute Lin_Xi and Lin_Eta
  double L=0.0;
  Epetra_SerialDenseMatrix B(3*(numnode1+1),1);
  
  // compute L elementwise
  for (int i=0;i<3;i++)
  {
    L +=  dx1[i]*dx1[i] + (x1[i]-x2[i])*ddx1[i];
  }

  if (L == 0) dserror("ERROR: L = 0");

  
  // index vectors for access to shape function vectors and matrices
  std::vector<int> index1(3*(numnode1));
  
  // fill the index vectors
  for (int i=0;i<3*numnode1;i++)
  {
    index1[i]=(int) floor(i/3);
  }
  
  // compute B elementwise
  for (int i=0;i<3*(numnode1+1);i++)
  {
    // first block
    if (i<3*numnode1)
    {
      int j=i%3;
      int row1 = index1[i];

      B(i,0) = -(x1[j]-x2[j])*deriv1(0,row1) - dx1[j]*funct1[row1];
    }
    // second block
    else
    {
      B(i,0) = dx1[i-3*numnode1];
    }
  }

  // finally the linearizations / directional derivatives
  for (int i=0;i<3*(numnode1+1);i++)
  {
    delta_xi[i] = B(i,0)/L;
  }
  
  return;
  
}
/*----------------------------------------------------------------------*
 |  end: Linearizations of contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute distance vector                                   popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeDistance(std::vector<double>& distance,
     double& normdist, const std::vector<double>& normal, const double& norm)
{
  // compute distance vector
  for (int i=0;i<(int)normal.size();i++)
    distance[i] = normal[i] * norm;
  
  // compute scalar distance
  normdist = norm;
  
  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute distance vector
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of gap                               popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeLinGap(std::vector<double>& delta_gap,
      std::vector<double>& delta_xi,
      const std::vector<double>& x1, const std::vector<double>& x2, 
      const std::vector<double>& dx1,
      const Epetra_SerialDenseVector& funct1,
      const double& normdist, const int& numnode1,
      const std::vector<double>& normal, const double& norm, const double& gap,
      Epetra_SerialDenseMatrix& delta_x1_minus_x2)
{
  // index vectors for access to shape function vectors and matrices
  std::vector<int> index1(3*(numnode1));
  
  // fill the index vectors
  for (int i=0;i<3*numnode1;i++)
  {
    index1[i]=(int) floor(i/3);
  }

  // compute linearization of disctance vector
  for (int i=0;i<3;i++)
  {
    for (int j=0;j<3*(numnode1+1);j++)
    {
      // standard part for each j
      delta_x1_minus_x2(i,j) = dx1[i]*delta_xi[j];

      // only for first block
      if (j<3*numnode1 && j%3 == i)
      {
        int row = index1[j];
        delta_x1_minus_x2(i,j) += funct1[row];
      }
      // only for second block
      else if(j%3 == i)
      {
        delta_x1_minus_x2(i,j) += -1.0;
      }
    }
  }
  
  // compute linearization of gap
  for (int i=0;i<3;i++)
    for (int j=0;j<3*(numnode1+1);j++)
      delta_gap[j] +=  (x1[i]-x2[i])/normdist * delta_x1_minus_x2(i,j);
  
  return;  
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of gap
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute linearization of normal                           popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::ComputeLinNormal(Epetra_SerialDenseMatrix& delta_n,
     const std::vector<double>& x1, const std::vector<double>& x2,
     const double& norm,const int& numnode1,
     const Epetra_SerialDenseMatrix& delta_x1_minus_x2, const std::vector<double>& normal, 
     const double& XiContact)
{
  // local vectors for shape functions and their derivatives
  Epetra_SerialDenseVector funct1(numnode1);
  Epetra_SerialDenseMatrix deriv1(1,numnode1);
  Epetra_SerialDenseMatrix secondderiv1(1,numnode1);

  // get shape functions and their derivatives
  GetShapeFunctions(funct1,deriv1,secondderiv1,XiContact);
    
  // tensor product of normal (x) normal
  Epetra_SerialDenseMatrix n_tp_n(3,3);
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      n_tp_n(i,j) = normal[i] * normal[j];

  // build a 3x3-identity matrix
  Epetra_SerialDenseMatrix identity(3,3);
  for (int i=0;i<3;i++) identity(i,i) = 1;
      
  // compute linearization of normal
  for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
      for (int k=0;k<3*(numnode1+1);k++)
        delta_n(i,k) += (identity(i,j) - n_tp_n(i,j)) * delta_x1_minus_x2(j,k) / norm;
  
  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute linearization of normal
 *----------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  Reset Uzawa-based Lagrange mutliplier                   popp 04/2010|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::Resetlmuzawa()
{
  lmuzawa_ = 0.0;
  return;
}
/*----------------------------------------------------------------------*
 |  end: Reset Uzawa-based Lagrange mutliplier
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update Uzawa-based Lagrange mutliplier                  popp 04/2010|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::Updatelmuzawa(const double& currentpp)
{
  // only update for active pairs, else reset
  if (contactflag_) lmuzawa_ -=  currentpp * GetGap();
  else              lmuzawa_ = 0.0;
  
  return;
}
/*----------------------------------------------------------------------*
 |  end: Update Uzawa-based Lagrange mutliplier
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update nodal coordinates (public)                           popp 04/10|
 *----------------------------------------------------------------------*/
void CONTACT::Beam3tospherecontact::UpdateElePos(Epetra_SerialDenseMatrix& newele1pos,
                                         Epetra_SerialDenseMatrix& newele2pos)
{
  ele1pos_ = newele1pos;
  ele2pos_ = newele2pos;
  
  return;
}
/*----------------------------------------------------------------------*
 |  end: Update nodal coordinates (public)
 *----------------------------------------------------------------------*/
