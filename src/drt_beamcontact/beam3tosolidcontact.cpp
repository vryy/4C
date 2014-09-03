/*!-----------------------------------------------------------------------------------------------------------
\file beam3tosolidcontact.cpp
\brief One beam and solid contact pair (two elements)

<pre>
Maintainer: Christoph Meier
            meier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15262
</pre>

*-----------------------------------------------------------------------------------------------------------*/

#include "beam3tosolidcontact.H"
#include "beam3contact_defines.H"
#include "beam3contact_utils.H"
#include "beam3contact_tangentsmoothing.H"
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
#include "../drt_beam3eb/beam3eb.H"
#include "../drt_lib/drt_element.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_inpar/inpar_statmech.H"

/*----------------------------------------------------------------------*
 |  constructor (public)                                     meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol ,const int numnodes , const int numnodalvalues>
CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Beam3tosolidcontact(
                                                                     const DRT::Discretization& pdiscret,
                                                                     const DRT::Discretization& cdiscret,
                                                                     const std::map<int,int>& dofoffsetmap,
                                                                     DRT::Element* element1,
                                                                     DRT::Element* element2,
                                                                     Teuchos::ParameterList beamcontactparams):
pdiscret_(pdiscret),
cdiscret_(cdiscret),
dofoffsetmap_(dofoffsetmap),
element1_(element1),
element2_(element2),
sgn_(1.0),
firstcall_(true),
lmuzawa_(0.0),
gap_(0.0),
gap_original_(0.0),
contactflag_(false),
elementscolinear_(false),
elementscrossing_(false),
shiftnodalvalues_(false),
xi1_(0.0),
xi2_(0.0)
{

  std::cout << "Warning: Currently only a dummy constructor is implemented!" << std::endl;

//  for (int i=0;i<3;i++)
//  {
//    r1_(i)=0.0;
//    r2_(i)=0.0;
//    normal_(i)=0.0;
//    normal_old_(i)=0.0;
//  }
//  for (int i=0;i<3*numnodes*numnodalvalues;i++)
//  {
//    ele1pos_(i)=0.0;
//    ele2pos_(i)=0.0;
//  }
//  for (int i=0;i<3*numnodes;i++)
//  {
//    nodaltangentssmooth1_(i)=0.0;
//    nodaltangentssmooth2_(i)=0.0;
//  }
//
//  ngf_= DRT::INPUT::IntegralValue<int>(beamcontactparams,"BEAMS_NEWGAP");
//  smoothing_ = DRT::INPUT::IntegralValue<INPAR::BEAMCONTACT::Smoothing>(beamcontactparams,"BEAMS_SMOOTHING");
//
//  const DRT::ElementType & eot1 = element1_->ElementType();
//
//  if (eot1 == DRT::ELEMENTS::Beam3Type::Instance())
//
//  if (smoothing_ == INPAR::BEAMCONTACT::bsm_cpp and eot1 != DRT::ELEMENTS::Beam3Type::Instance() and eot1 != DRT::ELEMENTS::Beam3iiType::Instance() )
//    dserror("Tangent smoothing only implemented for beams of type beam3 and beam3ii!");
//
//  //In case of tangent smoothing for both elements the 2 direct neighbor elements are determined and saved in the B3CNeighbor-Class
//  //variables neighbors1_ and neighbors2_
//  if (smoothing_ == INPAR::BEAMCONTACT::bsm_cpp)
//  {
//    neighbors1_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element1);
//    neighbors2_ = CONTACT::B3TANGENTSMOOTHING::DetermineNeigbors(element2);
//  }
//  else
//  {
//    neighbors1_=Teuchos::null;
//    neighbors2_=Teuchos::null;
//  }
//
//  if (element1->ElementType() != element2->ElementType())
//    dserror("The class beam3contact only works for contact pairs of the same beam element type!");

  return;
}
/*----------------------------------------------------------------------*
 |  end: constructor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  copy-constructor (public)                                meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Beam3tosolidcontact(const Beam3tosolidcontact& old):
pdiscret_(old.pdiscret_),
cdiscret_(old.cdiscret_),
dofoffsetmap_(old.dofoffsetmap_)
{
  dserror("ERROR: Copy constructor incomplete");
  return;
}
/*----------------------------------------------------------------------*
 |  end: copy-constructor
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate the element (public)                             meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Evaluate(
     LINALG::SparseMatrix& stiffmatrix,
     Epetra_Vector& fint,
     const double& pp)
{

  std::cout << "Warning: Currently only a dummy version of the method Evaluate is implemented!" << std::endl;

  std::cout << "element1_->Id(): " << element1_->Id() << std::endl;
  std::cout << "element2_->Id(): " << element2_->Id() << std::endl;

//  //**********************************************************************
//  // Evaluation of contact forces and stiffness
//  //**********************************************************************
//  // (1) Closest Point Projection (CPP)
//  //     -> find closest point where contact forces are evaluated
//  // (2) Compute some auxiliary quantities
//  //     -> normal vector, gap, shape functions, contact flag,
//  //     -> linearizations of all geometric quantities
//  // (3) Compute contact forces and stiffness
//  //     -> stiffness terms are directly assembled to global matrix
//  //     -> contact forces are only returned as global vector
//  // (4) Perform some finite difference checks
//  //     -> only if the flag BEAMCONTACTFDCHECKS is defined
//
//  //**********************************************************************
//  // (1) Closest Point Projection (CPP)
//  //**********************************************************************
//
//  ClosestPointProjection();
//
//  //**********************************************************************
//  // (2) Compute some auxiliary quantities
//  //**********************************************************************
//
//  // vectors for shape functions and their derivatives
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta
//
//  // coords and derivatives of the two contacting points
//  LINALG::TMatrix<BTSTYPE, 3, 1> r1(true);                               // = r1
//  LINALG::TMatrix<BTSTYPE, 3, 1> r2(true);                               // = r2
//  LINALG::TMatrix<BTSTYPE, 3, 1> r1_xi(true);                            // = r1,xi
//  LINALG::TMatrix<BTSTYPE, 3, 1> r2_xi(true);                            // = r2,eta
//  LINALG::TMatrix<BTSTYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
//  LINALG::TMatrix<BTSTYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
//  LINALG::TMatrix<BTSTYPE, 3, 1> delta_r(true);                          // = r1-r2
//  BTSTYPE norm_delta_r= 0.0;                                             // = g
//
//  // Calculate normal vector for all neighbor elements in order to have the quantity normal_old_ available
//  // in the next time step. This is necessary in order to avoid an undetected crossing of the beams when a beam
//  // element of one physical beam slides across the boundary node of two adjacent beam elements on the second physical beam
//
//  // this exludes pairs with IDs i and i+2, i.e. contact with the next but one element: It has proven useful that NEIGHBORTOL ca. 3,
//  //since the sum of |xi1| and |xi2| has to be larger than three for a contact pair consisting of element i and element i+2 of the same physical beam!
//  if ((BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(xi1_)) + BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(xi2_))) < NEIGHBORTOL && !elementscolinear_ && !elementscrossing_)
//  {
//    // update shape functions and their derivatives
//    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, xi1_, xi2_);
//    // update coordinates and derivatives of contact points
//    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);
//    // store coordinates of contact point into class variables
//    for(int i=0;i<3;i++)
//    {
//      r1_(i)=r1(i);
//      r2_(i)=r2(i);
//    }
//
//    // call function to compute scaled normal and gap of contact point and store this quantities in the corresponding
//    // class variables. The auxiliary variables delta_r and norm_delta_r might be useful for later applications.
//    ComputeNormal(delta_r, norm_delta_r);
//  }
//  else
//  {
//    contactflag_ = false;
//    return false;
//  }
//
//  // Check if the CPP found for this contact pair is really on the considered element, i.e. xi \in [-1;1]
//  if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(xi1_)) < (1.0 + XIETATOL) && BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(xi2_)) < (1.0 + XIETATOL))
//  {
//    //std::cout << "Auswertung von Paar:" << element1_->Id() << "/" << element2_->Id() << std::endl;
//  }
//  else
//  {
//    contactflag_ = false;
//    return false;
//  }
//
//  // call function to evaluate contact status
//  CheckContactStatus(pp);
//
//  //**********************************************************************
//  // (3) Compute contact forces and stiffness
//  //**********************************************************************
//
//  // call function to evaluate and assemble contact forces
//  EvaluateFcContact(pp, &fint, N1, N2);
//
//  // call function to evaluate and assemble contact stiffness
//  EvaluateStiffcContact(pp,norm_delta_r,delta_r,stiffmatrix,r1,r2,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi,N1_xixi,N2_xixi);

  return true;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate the element
  *---------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute contact forces                                   meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateFcContact(
     const double& pp,
     Epetra_Vector* fint,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2,
     LINALG::TMatrix<BTSTYPE, 3*numnodes*numnodalvalues, 1>* fc1_FAD,
     LINALG::TMatrix<BTSTYPE, 3*numnodes*numnodalvalues, 1>* fc2_FAD)
{
//  // get dimensions for vectors fc1 and fc2
//  const int dim1 = 3*numnodes*numnodalvalues;
//  const int dim2 = 3*numnodes*numnodalvalues;
//
//  // temporary vectors for contact forces, DOF-GIDs and owning procs
//  LINALG::TMatrix<BTSTYPE, dim1, 1> fc1(true);
//  LINALG::TMatrix<BTSTYPE, dim2, 1> fc2(true);
//  Epetra_SerialDenseVector fcontact1(dim1);
//  Epetra_SerialDenseVector fcontact2(dim2);
//  std::vector<int>  lm1(dim1);
//  std::vector<int>  lm2(dim2);
//  std::vector<int>  lmowner1(dim1);
//  std::vector<int>  lmowner2(dim2);
//
//  // flag indicating assembly
//  bool DoNotAssemble = false;
//
//  //**********************************************************************
//  // evaluate contact forces for active pairs
//  //**********************************************************************
//  if (contactflag_)
//  {
//    // node ids of both elements
//    const int* node_ids1 = element1_->NodeIds();
//    const int* node_ids2 = element2_->NodeIds();
//
//    //********************************************************************
//    // Compute Fc1 (force acting on first element)
//    //********************************************************************
//    for (int i=0;i<dim1;++i)
//    {
//      for (int j=0;j<3;++j)
//      {
//        fc1(i) +=  sgn_*N1(j,i)*normal_(j)*(lmuzawa_ - pp*gap_);
//      }
//    }
//
//    for (int i=0;i<numnodes;++i)
//    {
//      // get node pointer and dof ids
//      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
//      std::vector<int> NodeDofGIDs = GetGlobalDofs(node);
//
//      // compute force vector Fc1 and prepare assembly
//      for (int j=0;j<3*numnodalvalues;++j)
//      {
//        lm1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
//        lmowner1[3*numnodalvalues*i+j] = node->Owner();
//      }
//    }
//    //********************************************************************
//    // Compute Fc2 (force acting on second element)
//    //********************************************************************
//    for (int i=0;i<dim1;++i)
//    {
//      for (int j=0;j<3;++j)
//      {
//        fc2(i) +=  -sgn_*N2(j,i)*normal_(j)*(lmuzawa_ - pp*gap_);
//      }
//    }
//
//    for (int i=0;i<numnodes;++i)
//    {
//      // get node pointer and dof ids
//      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
//      std::vector<int> NodeDofGIDs = GetGlobalDofs(node);
//
//      // compute force vector Fc1 and prepare assembly
//      for (int j=0;j<3*numnodalvalues;++j)
//      {
//        lm2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
//        lmowner2[3*numnodalvalues*i+j] = node->Owner();
//      }
//    }
//    #ifdef AUTOMATICDIFF
//    if (fc1_FAD != NULL and fc2_FAD != NULL)
//    {
//      for (int i=0;i<dim1;++i)
//      {
//          (*fc1_FAD)(i) = fc1(i);
//      }
//      for (int i=0;i<dim2;++i)
//      {
//          (*fc2_FAD)(i) = fc2(i);
//      }
//    }
//    #endif
//  }
//  //**********************************************************************
//  // no forces for inactive pairs
//  //**********************************************************************
//  else
//  {
//    // set flag to avoid assembly
//    DoNotAssemble = true;
//
//    // compute forces
//    for (int i=0;i<3*numnodes*numnodalvalues;++i) fc1(i)=0;
//    for (int i=0;i<3*numnodes*numnodalvalues;++i) fc2(i)=0;
//  }
//  //**********************************************************************
//  // assemble contact forces
//  //**********************************************************************
//  if (!DoNotAssemble and fint != NULL)
//  {
//    for (int i=0;i<dim1;++i)
//    {
//        fcontact1[i] = BEAMCONTACT::CastToDouble(fc1(i));
//    }
//    for (int i=0;i<dim2;++i)
//    {
//        fcontact2[i] = BEAMCONTACT::CastToDouble(fc2(i));
//    }
//    // assemble fc1 and fc2 into global contact force vector
//    LINALG::Assemble(*fint,fcontact1,lm1,lmowner1);
//    LINALG::Assemble(*fint,fcontact2,lm2,lmowner2);
//
//  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute contact forces
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate contact stiffness                               meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateStiffcContact(
     const double& pp,
     const BTSTYPE& norm_delta_r,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& delta_r,
     LINALG::SparseMatrix& stiffmatrix,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& r1,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& r2,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xi,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xi,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xixi,
     const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xixi,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1_xixi,
     const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi)
{
//  // get dimensions for vectors fc1 and fc2
//  const int dim1 = 3*numnodes*numnodalvalues;
//  const int dim2 = 3*numnodes*numnodalvalues;
//
//  // temporary matrices for stiffness and vectors for DOF-GIDs and owning procs
//  LINALG::TMatrix<BTSTYPE, dim1, dim1+dim2> stiffc1(true);
//  LINALG::TMatrix<BTSTYPE, dim2, dim1+dim2> stiffc2(true);
//  Epetra_SerialDenseMatrix stiffcontact1(dim1,dim1+dim2);
//  Epetra_SerialDenseMatrix stiffcontact2(dim2,dim1+dim2);
//  std::vector<int>  lmrow1(dim1);
//  std::vector<int>  lmrow2(dim2);
//  std::vector<int>  lmrowowner1(dim1);
//  std::vector<int>  lmrowowner2(dim2);
//  std::vector<int>  lmcol1(dim1+dim2);
//  std::vector<int>  lmcol2(dim1+dim2);
//
//  // flag indicating assembly
//  bool DoNotAssemble = false;
//
//  //**********************************************************************
//  // evaluate contact stiffness for active pairs
//  //**********************************************************************
//  if (contactflag_)
//  {
//
//    // node ids of both elements
//    const int* node_ids1 = element1_->NodeIds();
//    const int* node_ids2 = element2_->NodeIds();
//
//    // initialize storage for linearizations
//    LINALG::TMatrix<BTSTYPE, dim1+dim2, 1>  delta_xi(true);
//    LINALG::TMatrix<BTSTYPE, dim1+dim2, 1> delta_eta(true);
//    LINALG::TMatrix<BTSTYPE, dim1+dim2, 1> delta_gap(true);
//    LINALG::TMatrix<BTSTYPE, 3, dim1+dim2> delta_x1_minus_x2(true);
//    LINALG::TMatrix<BTSTYPE, 3, dim1+dim2> delta_n(true);
//
//    //********************************************************************
//    // evaluate linearizations and distance
//    //********************************************************************
//    // linearization of contact point
//    ComputeLinXiAndLinEta(delta_xi,delta_eta,delta_r,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi);
//
//    #ifdef FADCHECKS
//      cout << "delta_xi: " << endl;
//      for (int i=0;i<dim1+dim2;i++)
//        cout << delta_xi(i).val() << "  ";
//      cout << endl << "delta_eta: " << endl;
//      for (int i=0;i<dim1+dim2;i++)
//        cout << delta_eta(i).val() << "  ";
//      cout << endl;
//      FADCheckLinXiAndLinEta(delta_r,r1_xi,r2_xi,r1_xixi,r2_xixi,N1,N2,N1_xi,N2_xi);
//    #endif
//
//    // linearization of gap function which is equal to delta d
//    ComputeLinGap(delta_gap,delta_xi,delta_eta,delta_r,norm_delta_r,r1_xi,r2_xi,N1,N2);
//
//    // linearization of normal vector
//    ComputeLinNormal(delta_n,delta_xi,delta_eta,norm_delta_r,r1_xi,r2_xi,N1,N2);
//
//    //********************************************************************
//    // prepare assembly
//    //********************************************************************
//    // fill lmrow1 and lmrowowner1
//    for (int i=0;i<numnodes;++i)
//    {
//      // get pointer and dof ids
//      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
//      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);
//
//      for (int j=0;j<3*numnodalvalues;++j)
//      {
//        lmrow1[3*numnodalvalues*i+j]=NodeDofGIDs[j];
//        lmrowowner1[3*numnodalvalues*i+j]=node->Owner();
//      }
//    }
//
//    // fill lmrow2 and lmrowowner2
//    for (int i=0;i<numnodes;++i)
//    {
//      // get pointer and node ids
//      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
//      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);
//
//      for (int j=0;j<3*numnodalvalues;++j)
//      {
//        lmrow2[3*numnodalvalues*i+j]=NodeDofGIDs[j];
//        lmrowowner2[3*numnodalvalues*i+j]=node->Owner();
//      }
//    }
//
//    // fill lmcol1 and lmcol2
//    for (int i=0;i<numnodes;++i)
//    {
//      // get pointer and node ids
//      DRT::Node* node = ContactDiscret().gNode(node_ids1[i]);
//      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);
//
//      for (int j=0;j<3*numnodalvalues;++j)
//      {
//        lmcol1[3*numnodalvalues*i+j] = NodeDofGIDs[j];
//        lmcol2[3*numnodalvalues*i+j] = NodeDofGIDs[j];
//      }
//    }
//
//    // fill lmcol1 and lmcol2
//    for (int i=0;i<numnodes;++i)
//    {
//      // get pointer and node ids
//      DRT::Node* node = ContactDiscret().gNode(node_ids2[i]);
//      std::vector<int> NodeDofGIDs =  GetGlobalDofs(node);
//
//      for (int j=0;j<3*numnodalvalues;++j)
//      {
//        lmcol1[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
//        lmcol2[3*numnodalvalues*numnodes+3*numnodalvalues*i+j] = NodeDofGIDs[j];
//      }
//    }
//
//    //********************************************************************
//    // evaluate contact stiffness
//    // (1) stiffc1 of first element
//    //********************************************************************
//
//    //********************************************************************
//    // part I
//    //********************************************************************
//    LINALG::TMatrix<BTSTYPE,dim1,1> N1T_normal(true);
//    for (int i=0;i<3;i++)
//    {
//      for (int j=0;j<dim1;j++)
//      {
//        N1T_normal(j)+=N1(i,j)*normal_(i);
//      }
//    }
//
//    for (int i=0;i<dim1;i++)
//    {
//      for (int j=0;j<dim1+dim2;j++)
//      {
//        stiffc1(i,j) = -sgn_* pp * N1T_normal(i) * delta_gap(j);
//      }
//    }
//    //********************************************************************
//    // part II
//    //********************************************************************
//    for  (int i=0;i<3;i++)
//    {
//      for (int j=0;j<dim1;j++)
//      {
//        for (int k=0;k<dim1+dim2;k++)
//        {
//            stiffc1(j,k) += sgn_*(lmuzawa_ - pp*gap_)*N1(i,j)*delta_n(i,k);
//        }
//      }
//    }
//    //********************************************************************
//    // part III
//    //********************************************************************
//    LINALG::TMatrix<BTSTYPE,dim1,1> N1xiT_normal(true);
//    for (int i=0;i<3;i++)
//    {
//      for (int j=0;j<dim1;j++)
//      {
//        N1xiT_normal(j) += N1_xi(i,j)*normal_(i);
//      }
//    }
//
//    for (int i=0;i<dim1;i++)
//    {
//      for (int j=0;j<dim1+dim2;j++)
//      {
//        stiffc1(i,j) += sgn_*(lmuzawa_ - pp*gap_)*N1xiT_normal(i)*delta_xi(j);
//      }
//    }
//    //********************************************************************
//    // evaluate contact stiffness
//    // (2) stiffc2 of second element
//    //********************************************************************
//
//    //********************************************************************
//    // part I
//    //********************************************************************
//    LINALG::TMatrix<BTSTYPE,dim2,1> N2T_normal(true);
//    for (int i=0;i<3;i++)
//    {
//      for (int j=0;j<dim2;j++)
//      {
//        N2T_normal(j)+=N2(i,j)*normal_(i);
//      }
//    }
//
//    for (int i=0;i<dim2;i++)
//    {
//      for (int j=0;j<dim1+dim2;j++)
//      {
//        stiffc2(i,j) = sgn_* pp * N2T_normal(i) * delta_gap(j);
//      }
//    }
//    //********************************************************************
//    // part II
//    //********************************************************************
//    for  (int i=0;i<3;i++)
//    {
//      for (int j=0;j<dim2;j++)
//      {
//        for (int k=0;k<dim1+dim2;k++)
//        {
//            stiffc2(j,k) += -sgn_*(lmuzawa_ - pp*gap_)*N2(i,j)*delta_n(i,k);
//        }
//      }
//    }
//    //********************************************************************
//    // part III
//    //********************************************************************
//    LINALG::TMatrix<BTSTYPE,dim1,1> N2xiT_normal(true);
//    for (int i=0;i<3;i++)
//    {
//      for (int j=0;j<dim2;j++)
//      {
//        N2xiT_normal(j) += N2_xi(i,j)*normal_(i);
//      }
//    }
//
//    for (int i=0;i<dim2;i++)
//    {
//      for (int j=0;j<dim1+dim2;j++)
//      {
//        stiffc2(i,j) += -sgn_*(lmuzawa_ - pp*gap_)*N2xiT_normal(i)*delta_eta(j);
//      }
//    }
//
//    #ifdef AUTOMATICDIFF
//      LINALG::TMatrix<BTSTYPE, dim1, 1> fc1_FAD(true);
//      LINALG::TMatrix<BTSTYPE, dim2, 1> fc2_FAD(true);
//
//      EvaluateFcContact(pp, NULL, N1, N2, &fc1_FAD, &fc2_FAD);
//      LINALG::TMatrix<BTSTYPE, dim1, dim1+dim2> stiffc1_FAD(true);
//      LINALG::TMatrix<BTSTYPE, dim2, dim1+dim2> stiffc2_FAD(true);
//      for (int j=0;j<dim1+dim2;j++)
//      {
//        for (int i=0;i<dim1;i++)
//          stiffc1_FAD(i,j) = -(fc1_FAD(i).dx(j)+fc1_FAD(i).dx(dim1+dim2)*delta_xi(j)+fc1_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
//        for (int i=0;i<dim2;i++)
//          stiffc2_FAD(i,j) = -(fc2_FAD(i).dx(j)+fc2_FAD(i).dx(dim1+dim2)*delta_xi(j)+fc2_FAD(i).dx(dim1+dim2+1)*delta_eta(j));
//      }
//
//    #endif
//
//  } //if(contactflag_)
//
//  //**********************************************************************
//  // no stiffness for inactive pairs
//  //**********************************************************************
//  else
//  {
//     // set flag to avoid assembly
//    DoNotAssemble = true;
//
//    // compute stiffness
//    for (int i=0;i<dim1;i++)
//      for (int j=0;j<dim1+dim2;j++)
//        stiffc1(i,j) = 0;
//    for (int i=0;i<dim2;i++)
//      for (int j=0;j<dim1+dim2;j++)
//        stiffc2(i,j) = 0;
//  }
//
//  //**********************************************************************
//  // assemble contact stiffness
//  //**********************************************************************
//  // change sign of stiffc1 and stiffc2 due to time integration.
//  // according to analytical derivation there is no minus sign, but for
//  // our time integration methods the negative stiffness must be assembled.
//
//  // now finally assemble stiffc1 and stiffc2
//  if (!DoNotAssemble)
//  {
//    for (int j=0;j<dim1+dim2;j++)
//    {
//      for (int i=0;i<dim1;i++)
//        stiffcontact1(i,j) = -BEAMCONTACT::CastToDouble(stiffc1(i,j));
//      for (int i=0;i<dim2;i++)
//        stiffcontact2(i,j) = -BEAMCONTACT::CastToDouble(stiffc2(i,j));
//    }
//
//    stiffmatrix.Assemble(0,stiffcontact1,lmrow1,lmrowowner1,lmcol1);
//    stiffmatrix.Assemble(0,stiffcontact2,lmrow2,lmrowowner2,lmcol2);
//
////    cout << "Steifigkeitsmatrix: " << (*(stiffmatrix.EpetraMatrix())) << endl;
////    cout << "test: " << (*(stiffmatrix.EpetraMatrix()))[0,0] << endl;
//  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate contact stiffness
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Linearizations of contact point                          meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeLinXiAndLinEta(
                                                                          LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                                                                          LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& delta_r,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xixi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xixi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2_xi)
{
//  //**********************************************************************
//  // we have to solve the following system of equations:
//  //  _              _       _      _       _              _      _       _
//  // | L(1,1)  L(1,2) |    | Lin_Xi  |    |  B(1,1)  B(1,2) |   | Lin_d1 |
//  // |                | *  |         | =  |                 | * |        |
//  // |_L(2,1)  L(2,2)_|    |_Lin_Eta_|    |_B(2,1)  B(2,2)_ |   |_Lin_d2_|
//  //
//  // this can be done easily because it is a linear 2x2-system.
//  // we obtain the solution by inverting matrix L:
//  //
//  // [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
//  //
//  //**********************************************************************
//
//  const int dim1 = 3*numnodes*numnodalvalues;
//  const int dim2 = 3*numnodes*numnodalvalues;
//
//  // matrices to compute Lin_Xi and Lin_Eta
//  LINALG::TMatrix<BTSTYPE,2,2> L(true);
//  LINALG::TMatrix<BTSTYPE,2,2> L_inv(true);
//  LINALG::TMatrix<BTSTYPE,2,dim1+dim2> B(true);
//  LINALG::TMatrix<BTSTYPE,2,dim1+dim2> D(true);
//
//  // compute L elementwise
//  L(0,0)=::BEAMCONTACT::ScalarProduct(r1_xi, r1_xi) + ::BEAMCONTACT::ScalarProduct(delta_r, r1_xixi);
//  L(1,1)=-::BEAMCONTACT::ScalarProduct(r2_xi, r2_xi) + ::BEAMCONTACT::ScalarProduct(delta_r, r2_xixi);
//  L(0,1)=-::BEAMCONTACT::ScalarProduct(r2_xi, r1_xi);
//  L(1,0)=-L(0,1);
//
//  // invert L by hand
//  BTSTYPE det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
//  if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(det_L)) < DETERMINANTTOL) dserror("ERROR: Determinant of L = 0");
//  L_inv(0,0) =  L(1,1) / det_L;
//  L_inv(0,1) = -L(0,1) / det_L;
//  L_inv(1,0) = -L(1,0) / det_L;
//  L_inv(1,1) =  L(0,0) / det_L;
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim1;j++)
//    {
//      B(0,j)+= -delta_r(i)*N1_xi(i,j) - r1_xi(i)*N1(i,j);
//      B(1,j)+= - r2_xi(i)*N1(i,j);
//    }
//  }
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim2;j++)
//    {
//      B(0,j+dim1)+= r1_xi(i)*N2(i,j);
//      B(1,j+dim1)+= -delta_r(i)*N2_xi(i,j) + r2_xi(i)*N2(i,j);
//    }
//  }
//
//  // compute D = L^-1 * B
//  D.Multiply(L_inv, B);
//
//  // finally the linearizations / directional derivatives
//  for (int i=0;i<dim1+dim2;i++)
//  {
//    delta_xi(i) = D(0,i);
//    delta_eta(i) = D(1,i);
//  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Linearizations of contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of gap                              meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeLinGap(
                                                                       LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>& delta_gap,
                                                                       const LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                                                                       const LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                       const LINALG::TMatrix<BTSTYPE, 3, 1>& delta_r,
                                                                       const BTSTYPE& norm_delta_r,
                                                                       const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xi,
                                                                       const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xi,
                                                                       const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                       const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2)
{
//  const int dim1 = 3*numnodes*numnodalvalues;
//  const int dim2 = 3*numnodes*numnodalvalues;
//
//  // delta g := delta_r/||delta_r||*auxiliary_matri1 delta d, with auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))
//
//  LINALG::TMatrix<BTSTYPE,3,dim1+dim2>  auxiliary_matrix1(true);
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim1+dim2;j++)
//    {
//      auxiliary_matrix1(i,j)+=r1_xi(i)*delta_xi(j)-r2_xi(i)*delta_eta(j);
//    }
//  }
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim1;j++)
//    {
//      auxiliary_matrix1(i,j)+= N1(i,j);
//    }
//  }
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim2;j++)
//    {
//      auxiliary_matrix1(i,j+dim1)+= -N2(i,j);
//    }
//  }
//
//  // compute linearization of gap
//  for (int i=0;i<3;i++)
//    for (int j=0;j<dim1+dim2;j++)
//      delta_gap(j) +=  sgn_*delta_r(i) * auxiliary_matrix1(i,j)/norm_delta_r;

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of gap
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 | Compute linearization of normal vector                    meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeLinNormal(
                                                                          LINALG::TMatrix<BTSTYPE, 3, 2*3*numnodes*numnodalvalues>& delta_normal,
                                                                          const LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>&  delta_xi,
                                                                          const LINALG::TMatrix<BTSTYPE, 2*3*numnodes*numnodalvalues, 1>& delta_eta,
                                                                          const BTSTYPE& norm_delta_r,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xi,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                          const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2)
{
//  const int dim1 = 3*numnodes*numnodalvalues;
//  const int dim2 = 3*numnodes*numnodalvalues;
//
//  //delta n := auxiliary_matri2*auxiliary_matrix1* delta d, with auxiliary_matri2 = (I-nxn)/||r1-r2||
//  //and auxiliary_matri1 = (r1_xi*delta_xi-r2_xi*delta_eta + (N1, -N2))
//
//  LINALG::TMatrix<BTSTYPE,3,dim1+dim2>  auxiliary_matrix1(true);
//  LINALG::TMatrix<BTSTYPE,3,3>  auxiliary_matrix2(true);
//
//  //compute auxiliary_matrix1
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim1+dim2;j++)
//    {
//      auxiliary_matrix1(i,j)+=r1_xi(i)*delta_xi(j)-r2_xi(i)*delta_eta(j);
//    }
//  }
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim1;j++)
//    {
//      auxiliary_matrix1(i,j)+= N1(i,j);
//    }
//  }
//
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<dim2;j++)
//    {
//      auxiliary_matrix1(i,j+dim1)+= -N2(i,j);
//    }
//  }
//
//  //compute auxiliary_matrix2
//  for (int i=0;i<3;i++)
//  {
//    auxiliary_matrix2(i,i)+= 1.0/norm_delta_r;
//    for (int j=0;j<3;j++)
//    {
//      auxiliary_matrix2(i,j)+= -normal_(i)*normal_(j)/norm_delta_r;
//    }
//  }
//
//  // compute linearization of normal vector
//  for (int i=0;i<3;i++)
//    for (int j=0;j<3;j++)
//      for (int k=0;k<dim1+dim2;k++)
//        delta_normal(i,k) +=  auxiliary_matrix2(i,j) * auxiliary_matrix1(j,k);

  return;
}
/*----------------------------------------------------------------------*
 | end: Compute linearization of normal vector
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Closest point projection                                  meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ClosestPointProjection()
{
//
//  // local variables for element coordinates
//  BTSTYPE eta1=0.0;
//  BTSTYPE eta2=0.0;
//
//  // vectors for shape functions and their derivatives
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N1(true);        // = N1
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N2(true);        // = N2
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N1_xi(true);     // = N1,xi
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N2_xi(true);     // = N2,eta
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N1_xixi(true);   // = N1,xixi
//  LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues> N2_xixi(true);   // = N2,etaeta
//
//  // coords and derivatives of the two contacting points
//  LINALG::TMatrix<BTSTYPE, 3, 1> r1(true);                               // = r1
//  LINALG::TMatrix<BTSTYPE, 3, 1> r2(true);                               // = r2
//  LINALG::TMatrix<BTSTYPE, 3, 1> r1_xi(true);                            // = r1,xi
//  LINALG::TMatrix<BTSTYPE, 3, 1> r2_xi(true);                            // = r2,eta
//  LINALG::TMatrix<BTSTYPE, 3, 1> r1_xixi(true);                          // = r1,xixi
//  LINALG::TMatrix<BTSTYPE, 3, 1> r2_xixi(true);                          // = r2,etaeta
//  LINALG::TMatrix<BTSTYPE, 3, 1> delta_r(true);                          // = r1-r2
//
//  //Tangent and derivatives for tangent field smoothing (only for Reissner beams)
//  LINALG::TMatrix<BTSTYPE,3,1> t1(true);
//  LINALG::TMatrix<BTSTYPE,3,1> t1_xi(true);
//  LINALG::TMatrix<BTSTYPE,3,1> t2(true);
//  LINALG::TMatrix<BTSTYPE,3,1> t2_xi(true);
//
//  // initialize function f and Jacobian df for Newton iteration
//  LINALG::TMatrix<BTSTYPE,2,1> f(true);
//  LINALG::TMatrix<BTSTYPE,2,2> df(true);
//  LINALG::TMatrix<BTSTYPE,2,2> dfinv(true);
//
//  // initial scalar residual (L2-norm of f)
//  BTSTYPE residual = 0.0;
//
//  int iter=0;
//
//  //set these excluding criteria to false in the default case
//  elementscrossing_ = false;
//  shiftnodalvalues_ = false;
//
//  //**********************************************************************
//  // local Newton iteration
//  //**********************************************************************
//  for (int i=0;i<BEAMCONTACTMAXITER;++i)
//  {
//    iter++;
//    // reset shape function variables to zero
//    N1.Clear();
//    N2.Clear();
//    N1_xi.Clear();
//    N2_xi.Clear();
//    N1_xixi.Clear();
//    N2_xixi.Clear();
//
//    // update shape functions and their derivatives
//    GetShapeFunctions(N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi, eta1, eta2);
//    // update coordinates and derivatives of contact points
//    ComputeCoordsAndDerivs(r1, r2, r1_xi, r2_xi, r1_xixi, r2_xixi, N1, N2, N1_xi, N2_xi, N1_xixi, N2_xixi);
//    // compute delta_r=r1-r2
//    for(int i=0;i<3;i++)
//      delta_r(i)=r1(i)-r2(i);
//
//    // compute norm of difference vector to scale the equations
//    // (this yields better conditioning)
//    // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
//    // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
//    double norm_delta_r = BEAMCONTACT::CastToDouble(BEAMCONTACT::VectorNorm<3>(delta_r));
//
//    // The closer the beams get, the smaller is norm_delta_r, but
//    // norm_delta_r is not allowed to be too small, else numerical problems occur.
//    // It can happen quite often that the centerlines of two beam elements of the same physical beam
//    // cross in one point and norm_delta_r = 0. Since in this case |eta1|>1 and |eta2|>1 they will be sorted out later anyways.
//    if (norm_delta_r < NORMTOL)
//    {
//      // this exludes pairs with IDs i and i+2, i.e. contact with the next but one element
//      if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(eta1)) + BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(eta2)) < NEIGHBORTOL)
//      {
//        std::cout << "Warning! pair " << element1_->Id() << " / " << element2_->Id() << ": Nodal Values shifted! "<< std::endl;
//
//        //Shift nodal values of contact pair by a small pre-defined value in order to enable contact evaluation for contact pairs with r1=r2
//        //TODO: May this can be done in a more beautiful way in the future
//        // It is checked by the method pairs_[i]->GetShiftStatus() in the function void CONTACT::Beam3cmanager::UpdateConstrNorm(double* cnorm)
//        // that all active contact pairs fulfill shifnodalvalues_ = false in the converged configuration!!!
//        ShiftNodalPositions();
//        shiftnodalvalues_=true;
//        continue;
//
//      }
//      else
//      {
//        elementscrossing_ = true;
//        break;
//      }
//    }
//
//    // Evaluate nodal tangents in each case. However, they are used only if smoothing_=INPAR::BEAMCONTACT::bsm_cpp
//    CONTACT::B3TANGENTSMOOTHING::ComputeTangentsAndDerivs<numnodes,numnodalvalues>(t1, t1_xi, nodaltangentssmooth1_, N1, N1_xi);
//    CONTACT::B3TANGENTSMOOTHING::ComputeTangentsAndDerivs<numnodes,numnodalvalues>(t2, t2_xi, nodaltangentssmooth2_, N2, N2_xi);
//
//    // evaluate f at current eta1, eta2
//    EvaluateOrthogonalityCondition(f, delta_r, norm_delta_r, r1_xi, r2_xi, t1, t2);
//
//    // compute the scalar residuum
//    residual = BEAMCONTACT::VectorNorm<2>(f);
//
//    // check if Newton iteration has converged
//    if (BEAMCONTACT::CastToDouble(residual) < BEAMCONTACTTOL) break;
//
//    // evaluate Jacobian of f at current eta1, eta2
//    // Note: Parallel elements can not be handled with this beam contact formulation; such pairs
//    // are sorted out within ComputeCoordsAndDerivs(...) and the local Newton loop is terminated!
//    EvaluateLinOrthogonalityCondition(df, dfinv, delta_r, norm_delta_r, r1_xi, r2_xi, r1_xixi, r2_xixi, t1, t2, t1_xi, t2_xi);
//
//    if (elementscolinear_) break;
//
//    #ifdef FADCHECKS
//      cout << "df: " << df << endl;
//      BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(eta1, eta2);
//      FADCheckLinOrthogonalityCondition(delta_r,r1_xi,r2_xi);
//    #endif
//
//    // update element coordinates of contact point
//    eta1 += -dfinv(0,0)*f(0) - dfinv(0,1)*f(1);
//    eta2 += -dfinv(1,0)*f(0) - dfinv(1,1)*f(1);
//  }
//  //**********************************************************************
//
//  // Newton iteration unconverged after BEAMCONTACTMAXITER
//  if (residual > BEAMCONTACTTOL)
//  {
//    eta1 = 1e+12;
//    eta2 = 1e+12;
//  }
//
//  // store and return final result
//  xi1_=eta1;
//  xi2_=eta2;
//
//  // Set xi1_ and xi2_ as primary variables for automatic differentiation
//  // The dependence between the infinitesimal changes delta xi1_ and delta xi2_ and the
//  // the increments of the primary displacement variables delta disp have to be given explicitly, since
//  // no explicit relation between the finite quantities xi1_, xi2_ and disp exists.
//  // The latter would have been necessary if the full linearization had to be computed directly with Sacado!!!
//  #ifdef AUTOMATICDIFF
//    BEAMCONTACT::SetFADParCoordDofs<numnodes, numnodalvalues>(xi1_, xi2_);
//  #endif

  return;
}
/*----------------------------------------------------------------------*
|  end: Closest point projection
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  evaluate shape functions and derivatives                 meier 01/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetShapeFunctions(
                                                                            LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                            LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                            LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                            LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2_xi,
                                                                            LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1_xixi,
                                                                            LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2_xixi,
                                                                            const BTSTYPE& eta1,
                                                                            const BTSTYPE& eta2)
{
//  // get both discretization types
//  const DRT::Element::DiscretizationType distype1 = element1_->Shape();
//  const DRT::Element::DiscretizationType distype2 = element2_->Shape();
//
//  LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues> N1_i(true);
//  LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues> N1_i_xi(true);
//  LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues> N1_i_xixi(true);
//  LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues> N2_i(true);
//  LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues> N2_i_xi(true);
//  LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues> N2_i_xixi(true);
//
//  if (numnodalvalues==1)
//  {
//    // get values and derivatives of shape functions
//    DRT::UTILS::shape_function_1D(N1_i, eta1, distype1);
//    DRT::UTILS::shape_function_1D(N2_i, eta2, distype2);
//    DRT::UTILS::shape_function_1D_deriv1(N1_i_xi, eta1, distype1);
//    DRT::UTILS::shape_function_1D_deriv1(N2_i_xi, eta2, distype2);
//    DRT::UTILS::shape_function_1D_deriv2(N1_i_xixi, eta1, distype1);
//    DRT::UTILS::shape_function_1D_deriv2(N2_i_xixi, eta2, distype2);
//  }
//  else if (numnodalvalues==2)
//  {
//
//    if ( element1_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
//      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");
//
//    if ( element2_->ElementType() != DRT::ELEMENTS::Beam3ebType::Instance() )
//      dserror("Only elements of type Beam3eb are valid for the case numnodalvalues=2!");
//
//    double length1 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->jacobi();
//    double length2 = 2*(static_cast<DRT::ELEMENTS::Beam3eb*>(element2_))->jacobi();
//
//    // get values and derivatives of shape functions
//    DRT::UTILS::shape_function_hermite_1D(N1_i, eta1, length1, distype1);
//    DRT::UTILS::shape_function_hermite_1D(N2_i, eta2, length2, distype2);
//    DRT::UTILS::shape_function_hermite_1D_deriv1(N1_i_xi, eta1, length1, distype1);
//    DRT::UTILS::shape_function_hermite_1D_deriv1(N2_i_xi, eta2, length2, distype2);
//    DRT::UTILS::shape_function_hermite_1D_deriv2(N1_i_xixi, eta1, length1, distype1);
//    DRT::UTILS::shape_function_hermite_1D_deriv2(N2_i_xixi, eta2, length2, distype2);
//  }
//  else
//    dserror("Only beam elements with one (nodal positions) or two (nodal positions + nodal tangents) values are valid!");
//
//  //Assemble the individual shape functions in matrices, such that: r1=N1*d1, r1_xi=N1_xi*d1, r1_xixi=N1_xixi*d1, r2=N2*d2, r2_xi=N2_xi*d2, r2_xixi=N2_xixi*d2
//  AssembleShapefunctions(N1_i, N1_i_xi, N1_i_xixi, N1, N1_xi, N1_xixi);
//  AssembleShapefunctions(N2_i, N2_i_xi, N2_i_xixi, N2, N2_xi, N2_xixi);

  return;
}
/*----------------------------------------------------------------------*
 |  end: evaluate shape functions and derivatives
 *----------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------------------------*
 |  Assemble all shape functions                                                                  meier 01/14|
 *-----------------------------------------------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::AssembleShapefunctions(
                                                                                const LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues>& N_i,
                                                                                const LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues>& N_i_xi,
                                                                                const LINALG::TMatrix<BTSTYPE,1,numnodes*numnodalvalues>& N_i_xixi,
                                                                                LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N,
                                                                                LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N_xi,
                                                                                LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N_xixi)
{
//  //assembly_N is just an array to help assemble the matrices of the shape functions
//  //it determines, which shape function is used in which column of N
//  int assembly_N[3][3*numnodes*numnodalvalues];
//
//  //Initialize to zero
//  for (int i=0;i<3*numnodes*numnodalvalues;i++)
//  {
//    for (int j=0;j<3; j++)
//    {
//      assembly_N[j][i]=0.0;
//    }
//  }
//
//  /*
//  Set number of shape functions for each 3*3 block:
//  e.g. second order Reissner beam (numnodes=3, numnodalvalues=1)
//  int assembly_N[3][9]=  { {1,0,0,2,0,0,3,0,0},
//                           {0,1,0,0,2,0,0,3,0},
//                           {0,0,1,0,0,2,0,0,3}};
//
//  e.g. Kirchhoff beam (numnodes=2, numnodalvalues=2)
//  int assembly_N[3][12]=  {{1,0,0,2,0,0,3,0,0,4,0,0},
//                           {0,1,0,0,2,0,0,3,0,0,4,0},
//                           {0,0,1,0,0,2,0,0,3,0,0,4}};
//  */
//
//  for (int i=0;i<numnodes*numnodalvalues;i++)
//  {
//    assembly_N[0][3*i]=i+1;
//    assembly_N[1][3*i+1]=i+1;
//    assembly_N[2][3*i+2]=i+1;
//  }
//
//  //Assemble the matrices of the shape functions
//  for (int i=0; i<3*numnodes*numnodalvalues; i++)
//  {
//    for (int j=0; j<3; j++)
//    {
//      if(assembly_N[j][i]==0)
//      {
//        N(j,i)=0;
//        N_xi(j,i)=0;
//        N_xixi(j,i)=0;
//      }
//      else
//      {
//        N(j,i)=N_i(assembly_N[j][i]-1);
//        N_xi(j,i)=N_i_xi(assembly_N[j][i]-1);
//        N_xixi(j,i)=N_i_xixi(assembly_N[j][i]-1);
//      }
//    }
//  }

  return;
}

/*----------------------------------------------------------------------*
 | compute contact point coordinates and their derivatives   meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeCoordsAndDerivs(
                                                                                LINALG::TMatrix<BTSTYPE,3,1>& r1,
                                                                                LINALG::TMatrix<BTSTYPE,3,1>& r2,
                                                                                LINALG::TMatrix<BTSTYPE,3,1>& r1_xi,
                                                                                LINALG::TMatrix<BTSTYPE,3,1>& r2_xi,
                                                                                LINALG::TMatrix<BTSTYPE,3,1>& r1_xixi,
                                                                                LINALG::TMatrix<BTSTYPE,3,1>& r2_xixi,
                                                                                const LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N1,
                                                                                const LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N2,
                                                                                const LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N1_xi,
                                                                                const LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N2_xi,
                                                                                const LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N1_xixi,
                                                                                const LINALG::TMatrix<BTSTYPE,3,3*numnodes*numnodalvalues>& N2_xixi)
{
//  r1.Clear();
//  r2.Clear();
//  r1_xi.Clear();
//  r2_xi.Clear();
//  r1_xixi.Clear();
//  r2_xixi.Clear();
//
//#ifdef AUTOMATICDIFF
//  BEAMCONTACT::SetFADDispDofs<numnodes, numnodalvalues>(ele1pos_,ele2pos_);
//#endif
//
//  // compute output variable
//  for (int i=0;i<3;i++)
//  {
//    for (int j=0;j<3*numnodes*numnodalvalues;j++)
//    {
//      r1(i)+=N1(i,j)*ele1pos_(j);
//      r2(i)+=N2(i,j)*ele2pos_(j);
//      r1_xi(i)+=N1_xi(i,j)*ele1pos_(j);
//      r2_xi(i)+=N2_xi(i,j)*ele2pos_(j);
//      r1_xixi(i)+=N1_xixi(i,j)*ele1pos_(j);
//      r2_xixi(i)+=N2_xixi(i,j)*ele2pos_(j);
//    }
//  }

  return;
}
/*----------------------------------------------------------------------*
 | end: compute contact point coordinates and their derivatives         |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate function f in CPP                               meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateOrthogonalityCondition(
                                                                                        LINALG::TMatrix<BTSTYPE,2,1>& f,
                                                                                        const LINALG::TMatrix<BTSTYPE,3,1>& delta_r,
                                                                                        const double norm_delta_r,
                                                                                        const LINALG::TMatrix<BTSTYPE,3,1>& r1_xi,
                                                                                        const LINALG::TMatrix<BTSTYPE,3,1>& r2_xi,
                                                                                        const LINALG::TMatrix<BTSTYPE,3,1>& t1,
                                                                                        const LINALG::TMatrix<BTSTYPE,3,1>& t2)
{
//  // reset f
//  f.Clear();
//
//  // evaluate f
//  // see Wriggers, Computational Contact Mechanics, equation (12.5)
//  if (smoothing_ == INPAR::BEAMCONTACT::bsm_none) //non-smoothed
//  {
//    for (int i=0;i<3;i++)
//    {
//      f(0) += delta_r(i)*r1_xi(i) / norm_delta_r;
//      f(1) += -delta_r(i)*r2_xi(i) / norm_delta_r;
//    }
//  }
//  else //smoothed
//  {
//    for (int i=0;i<3;i++)
//    {
//      f(0) += delta_r(i)*t1(i) / norm_delta_r;
//      f(1) += -delta_r(i)*t2(i) / norm_delta_r;
//    }
//  }


  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate function f in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Evaluate Jacobian df in CPP                              meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::EvaluateLinOrthogonalityCondition(
                                                                                            LINALG::TMatrix<BTSTYPE,2,2>& df,
                                                                                            LINALG::TMatrix<BTSTYPE,2,2>& dfinv,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& delta_r,
                                                                                            const double norm_delta_r,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& r1_xi,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& r2_xi,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& r1_xixi,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& r2_xixi,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& t1,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& t2,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& t1_xi,
                                                                                            const LINALG::TMatrix<BTSTYPE,3,1>& t2_xi)

{
//  // reset df and dfinv
//  df.Clear();
//  dfinv.Clear();
//
//  // evaluate df
//  // see Wriggers, Computational Contact Mechanics, equation (12.7)
//
//  if (smoothing_ == INPAR::BEAMCONTACT::bsm_none) //non-smoothed
//  {
//    for(int i=0;i<3;i++)
//    {
//      df(0,0) += (r1_xi(i)*r1_xi(i) + delta_r(i)*r1_xixi(i)) / norm_delta_r;
//      df(0,1) += -r1_xi(i)*r2_xi(i) / norm_delta_r;
//      df(1,0) += -r2_xi(i)*r1_xi(i) / norm_delta_r;
//      df(1,1) += (r2_xi(i)*r2_xi(i) - delta_r(i)*r2_xixi(i)) / norm_delta_r;
//    }
//  }
//  else //smoothed
//  {
//    for(int i=0;i<3;i++)
//    {
//      df(0,0) += (r1_xi(i)*t1(i) + delta_r(i)*t1_xi(i)) / norm_delta_r;
//      df(0,1) += -t1(i)*r2_xi(i) / norm_delta_r;
//      df(1,0) += -t2(i)*t1_xi(i) / norm_delta_r;
//      df(1,1) += (r2_xi(i)*t2(i) - delta_r(i)*t2_xi(i)) / norm_delta_r;
//    }
//  }
//
//  // Inverting (2x2) matrix df by hard coded formula, so that it is
//  // possible to handle colinear vectors, because they lead to det(df) =0
//  BTSTYPE det_df = df(0,0)*df(1,1)-df(1,0)*df(0,1);
//
//  //********************************************************************
//  // ASSUMPTION:
//  // If det_df=0 we assume, that the two elements have an identical
//  // neutral axis. These contact objects will be rejected. The outcome
//  // of this physically rare phenomenon is that handling of line contact
//  // is not possible with this approach.
//  //********************************************************************
//
//  // singular df
//  if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(det_df))<COLINEARTOL)
//  {
//    // sort out
//    elementscolinear_ = true;
//  }
//  // regular df (inversion possible)
//  else
//  {
//    // do not sort out
//    elementscolinear_ = false;
//
//    // invert df
//    dfinv(0,0)=df(1,1)/det_df;
//    dfinv(0,1)=-df(0,1)/det_df;
//    dfinv(1,0)=-df(1,0)/det_df;
//    dfinv(1,1)=df(0,0)/det_df;
//  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Evaluate Jacobian df in CPP
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Compute normal vector in contact point                   meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ComputeNormal(LINALG::TMatrix<BTSTYPE,3,1>& delta_r,
                                                                                        BTSTYPE& norm_delta_r)
{
//  // compute non-unit normal
//  for (int i=0;i<3;i++) delta_r(i) = r1_(i)-r2_(i);
//
//  // compute length of normal
//  norm_delta_r = BEAMCONTACT::VectorNorm<3>(delta_r);
//
//  if (BEAMCONTACT::CastToDouble(norm_delta_r) < NORMTOL)
//    dserror("ERROR: Normal of length zero! --> change time step!");
//
//  // compute unit normal
//  for (int i=0;i<3;i++)
//  {
//    normal_(i) = delta_r(i) / norm_delta_r;
//  }
//
//  //In the first call after the contact pair element has been created the vector normal_old is initialized with the current normal vector
//  if (firstcall_)
//  {
//    for (int i=0;i<3;i++)
//    {
//      normal_old_(i) = normal_(i);
//    }
//    firstcall_ = false;
//  }
//
//  //Calculation of element radius via moment of inertia (only valid for beams with circular cross section!!!)
//  const DRT::ElementType & eot1 = element1_->ElementType();
//
//  double Iyy1 = 0.0;
//  double Iyy2 = 0.0;
//  if (eot1 == DRT::ELEMENTS::Beam3Type::Instance())
//  {
//    Iyy1 = (static_cast<DRT::ELEMENTS::Beam3*>(element1_))->Iyy();
//    Iyy2 = (static_cast<DRT::ELEMENTS::Beam3*>(element2_))->Iyy();
//  }
//  else if (eot1 == DRT::ELEMENTS::Beam3iiType::Instance())
//  {
//    Iyy1 = (static_cast<DRT::ELEMENTS::Beam3ii*>(element1_))->Iyy();
//    Iyy2 = (static_cast<DRT::ELEMENTS::Beam3ii*>(element2_))->Iyy();
//  }
//  else if (eot1 == DRT::ELEMENTS::Beam3ebType::Instance())
//  {
//    Iyy1 = (static_cast<DRT::ELEMENTS::Beam3eb*>(element1_))->Iyy();
//    Iyy2 = (static_cast<DRT::ELEMENTS::Beam3eb*>(element2_))->Iyy();
//  }
//  double radius1 = MANIPULATERADIUS * sqrt(sqrt(4 * Iyy1 / M_PI));
//  double radius2 = MANIPULATERADIUS * sqrt(sqrt(4 * Iyy2 / M_PI));
//
//  BTSTYPE gap=0.0;
//  sgn_=1.0;
//  if (ngf_)
//  {
//    if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(BEAMCONTACT::ScalarProduct(normal_, normal_old_))) < NORMALTOL)
//      dserror("ERROR: Rotation too large! --> Choose smaller Time step!");
//
//    gap = BEAMCONTACT::Signum(BEAMCONTACT::ScalarProduct(normal_, normal_old_))*norm_delta_r - radius1 - radius2;
//    sgn_= BEAMCONTACT::CastToDouble(BEAMCONTACT::Signum(BEAMCONTACT::ScalarProduct(normal_, normal_old_)));
//  }
//  else
//  {
//    gap = norm_delta_r - radius1 - radius2;
//  }
//
//  // also set class variable
//  gap_ = gap;
//
//  // for comparison reasons we calculate in each case additionally the original gap function defintion,
//  // thus gap_original==gap_ for ngf_==false
//  gap_original_ = norm_delta_r - radius1 - radius2;

  return;
}
/*----------------------------------------------------------------------*
 |  end: Compute normal vector in contact point
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Check if conact is active or inactive                    meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::CheckContactStatus(const double& pp)
{
//  // check contact condition valid for both pure penalty
//  // (lmuzawa = 0) and augmented lagrange (lmuzawa != 0)
//  if (lmuzawa_ - pp * gap_ > 0) contactflag_ = true;
//  else                          contactflag_ = false;

  return;
}
/*----------------------------------------------------------------------*
 |  end: Check if conact is active or inactive
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Get global dofs of a node                                 meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
std::vector<int> CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetGlobalDofs(const DRT::Node* node)
{
  // get dofs in beam contact discretization
  const std::vector<int> cdofs = ContactDiscret().Dof(node);

  // get dofs in problem discretization via offset
  std::vector<int> pdofs((int)(cdofs.size()));
  for (int k=0;k<(int)(cdofs.size());++k)
  {
    pdofs[k]=(dofoffsetmap_.find(cdofs[k]))->second;
  }

  return pdofs;
}
/*----------------------------------------------------------------------*
 |  end: Get global dofs of a node
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Reset Uzawa-based Lagrange mutliplier                  meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Resetlmuzawa()
{
  //lmuzawa_ = 0.0;

  return;
}
/*----------------------------------------------------------------------*
 |  end: Reset Uzawa-based Lagrange mutliplier
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Change the sign of the normal vector                   meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::InvertNormal()
{
//  for (int i=0; i<3;i++)
//    normal_(i) = -normal_(i);
}
/*----------------------------------------------------------------------*
 |  end: Change the sign of the old normal vector                       |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Shift current normal vector to old normal vector       meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::UpdateClassVariablesStep()
{
//  for (int j=0;j<3;j++)
//        normal_old_(j) = normal_(j);

}
/*----------------------------------------------------------------------*
 |  end: Shift current normal vector to old normal vector
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Check if there is a difference of old and new gap      meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
bool CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::GetNewGapStatus()
{
//  BTSTYPE gap_diff = gap_-gap_original_;
//
//  if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(gap_diff)) < GAPTOL)
//    return false;
//
//  else
    return true;


}
/*----------------------------------------------------------------------*
 |  end: Check if there is a difference of old and new gap               |
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update Uzawa-based Lagrange mutliplier                 meier 02/2014|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::Updatelmuzawa(const double& currentpp)
{
//  // only update for active pairs, else reset
//  if (contactflag_) lmuzawa_ -=  currentpp * GetGap();
//  else              lmuzawa_ = 0.0;

  return;
}
/*----------------------------------------------------------------------*
 |  end: Update Uzawa-based Lagrange mutliplier
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update nodal coordinates (public)                        meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::UpdateElePos( Epetra_SerialDenseMatrix& newele1pos,
                                                                                        Epetra_SerialDenseMatrix& newele2pos)
{
//  for (int i=0;i<3*numnodalvalues;i++)
//  {
//    for (int j=0;j<numnodes;j++)
//    {
//      ele1pos_(3*numnodalvalues*j+i)=newele1pos(i,j);
//      ele2pos_(3*numnodalvalues*j+i)=newele2pos(i,j);
//    }
//  }

  return;
}
/*----------------------------------------------------------------------*
 |  end: Update nodal coordinates (public)
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Update nodal tangents for tangent smoothing (public)      meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::UpdateEleSmoothTangents(std::map<int,LINALG::Matrix<3,1> >& currentpositions)
{

//  //Tangent smoothing is only possible for Reissner beam elements --> dserror() otherwise
//  if (numnodalvalues>1)
//    dserror("Tangent smoothing only possible for Reissner beam elements (numnodalvalues=1)!!!");
//
//  LINALG::Matrix<3*numnodes,1> elepos_aux(true);
//  //Tangent smoothing only possible with data type double (not with Sacado FAD)
//  for (int i=0;i<3*numnodes;i++)
//    elepos_aux(i)=BEAMCONTACT::CastToDouble(ele1pos_(i));
//
//  nodaltangentssmooth1_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element1_,neighbors1_);
//
//  elepos_aux.Clear();
//  //Tangent smoothing only possible with data type double (not with Sacado FAD)
//  for (int i=0;i<3*numnodes;i++)
//    elepos_aux(i)=BEAMCONTACT::CastToDouble(ele2pos_(i));
//
//  nodaltangentssmooth2_=CONTACT::B3TANGENTSMOOTHING::CalculateNodalTangents<numnodes>(currentpositions,elepos_aux ,element2_,neighbors2_);

}
/*----------------------------------------------------------------------*
 |  end: Update nodal coordinates (public)
 *----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  Shift Nodal positions (public)                           meier 02/14|
 *----------------------------------------------------------------------*/
template<const int numnodessol, const int numnodes , const int numnodalvalues>
void CONTACT::Beam3tosolidcontact<numnodessol, numnodes, numnodalvalues>::ShiftNodalPositions()
{
//  //Reissner beams
//  if (numnodalvalues == 1)
//  {
//    for (int i=0; i<numnodes; i++)
//    {
//      for (int j=0;j<3;j++)
//      {
//        ele1pos_(3*i + j) = ele1pos_(3*i + j) + SHIFTVALUE * normal_old_(j);
//      }
//    }
//  }
//  //Kirchhoff beams
//  else if (numnodalvalues == 2)
//  {
//    if (numnodes == 2)
//    {
//      for (int j=0;j<3;j++)
//      {
//        ele1pos_(j) = ele1pos_(j) + SHIFTVALUE * normal_old_(j);
//        ele1pos_(6+j) = ele1pos_(6+j) + SHIFTVALUE * normal_old_(j);
//      }
//    }
//    else
//    {
//      dserror("Only numnodes = 2 possible for Kirchhoff beams!!!");
//    }
//  }
//  else
//  {
//    dserror("The parameter numnodalvalues can only have the values 1 or 2!!!");
//  }

  return;
}


Teuchos::RCP<CONTACT::Beam3tosolidcontactinterface> CONTACT::Beam3tosolidcontactinterface::Impl( const int numnodessol,
                                                                                                 const int numnodes,
                                                                                                 const int numnodalvalues,
                                                                                                 const DRT::Discretization& pdiscret,
                                                                                                 const DRT::Discretization& cdiscret,
                                                                                                 const std::map<int,int>& dofoffsetmap,
                                                                                                 DRT::Element* element1,
                                                                                                 DRT::Element* element2,
                                                                                                 Teuchos::ParameterList beamcontactparams)
{

  if (numnodalvalues!=1 and numnodalvalues!=2)
    dserror("Only the values 1 and 2 are valid for numnodalvalues!");

  if (numnodalvalues!=2 and numnodes!=2)
    dserror("Only the values numnodes=2 is possible for Kirchhoff beams, i.e. if numnodalvalues=2!");

  if (numnodes!=2 and numnodes!=3 and numnodes!=4 and numnodes!=5)
    dserror("Only the values 2, 3, 4 and 5 are valid for numnodes!");

  if (numnodessol!=3 and numnodessol!=6 and numnodessol!=4 and numnodessol!=8 and numnodessol!=9)
    dserror("Only the values 3, 4, 6, 8 and 9 are valid for numnodessol!");


  switch (numnodessol)
  {
    case 3:
    {
      switch (numnodalvalues)
      {
        case 1:
        {
          switch (numnodes)
          {
            case 2:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 3:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 4:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 5:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
          }
          break;
        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<3,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
    case 4:
    {
      switch (numnodalvalues)
      {
        case 1:
        {
          switch (numnodes)
          {
            case 2:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 3:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 4:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 5:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
          }
          break;
        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<4,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
    case 6:
    {
      switch (numnodalvalues)
      {
        case 1:
        {
          switch (numnodes)
          {
            case 2:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 3:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 4:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 5:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
          }
          break;
        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<6,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
    case 8:
    {
      switch (numnodalvalues)
      {
        case 1:
        {
          switch (numnodes)
          {
            case 2:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 3:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 4:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 5:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
          }
          break;
        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<8,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
    case 9:
    {
      switch (numnodalvalues)
      {
        case 1:
        {
          switch (numnodes)
          {
            case 2:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,2,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 3:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,3,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 4:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,4,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
            case 5:
            {
              return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,5,1>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
            }
          }
          break;
        }
        case 2:
        {
          return Teuchos::rcp (new CONTACT::Beam3tosolidcontact<9,2,2>(pdiscret,cdiscret,dofoffsetmap,element1,element2,beamcontactparams));
        }
      }
      break;
    }
  }

  return Teuchos::null;

}

#ifdef BTSFADCHECKS
  /*----------------------------------------------------------------------*
   |  FAD-Check for Linearizations of contact point            meier 02/14|
   *----------------------------------------------------------------------*/
  template<const int numnodessol, const int numnodes , const int numnodalvalues>
  void CONTACT::Beam3tosolcontact<numnodessol, numnodes, numnodalvalues>::FADCheckLinXiAndLinEta(const LINALG::TMatrix<BTSTYPE, 3, 1>& delta_r,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xi,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xi,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xixi,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xixi,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N1_xi,
                                                                                  const LINALG::TMatrix<BTSTYPE, 3, 3*numnodes*numnodalvalues>& N2_xi)
  {
//    LINALG::TMatrix<BTSTYPE,2,1>f(true);
//
//    // compute norm of difference vector to scale the equations
//    // (this yields better conditioning)
//    // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
//    // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
//    double norm_delta_r = BEAMCONTACT::CastToDouble(BEAMCONTACT::VectorNorm<3>(delta_r));
//
//    // evaluate f of CCP condition
//    // see Wriggers, Computational Contact Mechanics, equation (12.5)
//    for (int i=0;i<3;i++)
//    {
//      f(0) += delta_r(i)*r1_xi(i) / norm_delta_r;
//      f(1) += -delta_r(i)*r2_xi(i) / norm_delta_r;
//    }
//
//    //**********************************************************************
//    // we have to solve the following system of equations:
//    //  _              _       _      _       _              _      _       _
//    // | L(1,1)  L(1,2) |    | Lin_Xi  |    |  B(1,1)  B(1,2) |   | Lin_d1 |
//    // |                | *  |         | =  |                 | * |        |
//    // |_L(2,1)  L(2,2)_|    |_Lin_Eta_|    |_B(2,1)  B(2,2)_ |   |_Lin_d2_|
//    //
//    // this can be done easily because it is a linear 2x2-system.
//    // we obtain the solution by inverting matrix L:
//    //
//    // [Lin_Xi; Lin_Eta] = L^-1 * B * [Lin_d1; Lin_d2] = D * [Lin_d1; Lin_d2]
//    //
//    //**********************************************************************
//
//    const int dim1 = 3*numnodes*numnodalvalues;
//    const int dim2 = 3*numnodes*numnodalvalues;
//
//    // matrices to compute Lin_Xi and Lin_Eta
//    LINALG::TMatrix<BTSTYPE,2,2> L(true);
//    LINALG::TMatrix<BTSTYPE,2,2> L_inv(true);
//    LINALG::TMatrix<BTSTYPE,2,dim1+dim2> B(true);
//    LINALG::TMatrix<BTSTYPE,2,dim1+dim2> D(true);
//
//    // compute L elementwise
//    L(0,0)= f(0).dx(2*3*numnodes*numnodalvalues);
//    L(0,1)= f(0).dx(2*3*numnodes*numnodalvalues+1);
//    L(1,0)= f(1).dx(2*3*numnodes*numnodalvalues);
//    L(1,1)= f(1).dx(2*3*numnodes*numnodalvalues+1);
//
//    // invert L by hand
//    BTSTYPE det_L = L(0,0)*L(1,1) - L(0,1)*L(1,0);
//    if (BEAMCONTACT::CastToDouble(BEAMCONTACT::Norm(det_L)) < DETERMINANTTOL) dserror("ERROR: Determinant of L = 0");
//    L_inv(0,0) =  L(1,1) / det_L;
//    L_inv(0,1) = -L(0,1) / det_L;
//    L_inv(1,0) = -L(1,0) / det_L;
//    L_inv(1,1) =  L(0,0) / det_L;
//
//    for (int j=0;j<dim1+dim2;j++)
//    {
//      B(0,j)= -f(0).dx(j);
//      B(1,j)= -f(1).dx(j);
//    }
//
//    // compute D = L^-1 * B
//    D.Multiply(L_inv, B);
//
//    cout << "linxi and lineta: " << endl;
//
//    cout << D << endl;

    return;
  }
  /*----------------------------------------------------------------------*
   |  End: FAD-Check for Linearizations of contact point
   *----------------------------------------------------------------------*/

  /*----------------------------------------------------------------------*
   |  FAD-Check for Linearizations of CCP                      meier 02/14|
   *----------------------------------------------------------------------*/
  template<const int numnodessol, const int numnodes , const int numnodalvalues>
  void CONTACT::Beam3tosolcontact<numnodessol, numnodes, numnodalvalues>::FADCheckLinOrthogonalityCondition( const LINALG::TMatrix<BTSTYPE, 3, 1>& delta_r,
                                                                                                             const LINALG::TMatrix<BTSTYPE, 3, 1>& r1_xi,
                                                                                                             const LINALG::TMatrix<BTSTYPE, 3, 1>& r2_xi)
  {
//    LINALG::TMatrix<BTSTYPE,2,1>f(true);
//
//    // compute norm of difference vector to scale the equations
//    // (this yields better conditioning)
//    // Note: Even if automatic differentiation via FAD is applied, norm_delta_r has to be of type double
//    // since this factor is needed for a pure scaling of the nonlinear CCP and has not to be linearized!
//    double norm_delta_r = BEAMCONTACT::CastToDouble(BEAMCONTACT::VectorNorm<3>(delta_r));
//
//    // evaluate f of CCP condition
//    // see Wriggers, Computational Contact Mechanics, equation (12.5)
//    for (int i=0;i<3;i++)
//    {
//      f(0) += delta_r(i)*r1_xi(i) / norm_delta_r;
//      f(1) += -delta_r(i)*r2_xi(i) / norm_delta_r;
//    }
//
//    LINALG::TMatrix<BTSTYPE,2,2>df(true);
//
//    for (int i=0;i<2;i++)
//    {
//      for (int j=0;j<2;j++)
//      {
//        df(i,j)=f(i).dx(2*3*numnodes*numnodalvalues+j);
//      }
//    }
//
//    cout << "df: " << endl;
//
//    cout << df << endl;

    return;
  }
  /*----------------------------------------------------------------------*
   |  End: FAD-Check for Linearizations of CPP
   *----------------------------------------------------------------------*/
#endif //#ifdef BTSFADCHECKS

//Possible template cases: this is necessary for the compiler
template class CONTACT::Beam3tosolidcontact<3,2,1>;
template class CONTACT::Beam3tosolidcontact<3,3,1>;
template class CONTACT::Beam3tosolidcontact<3,4,1>;
template class CONTACT::Beam3tosolidcontact<3,5,1>;
template class CONTACT::Beam3tosolidcontact<3,2,2>;

template class CONTACT::Beam3tosolidcontact<4,2,1>;
template class CONTACT::Beam3tosolidcontact<4,3,1>;
template class CONTACT::Beam3tosolidcontact<4,4,1>;
template class CONTACT::Beam3tosolidcontact<4,5,1>;
template class CONTACT::Beam3tosolidcontact<4,2,2>;

template class CONTACT::Beam3tosolidcontact<6,2,1>;
template class CONTACT::Beam3tosolidcontact<6,3,1>;
template class CONTACT::Beam3tosolidcontact<6,4,1>;
template class CONTACT::Beam3tosolidcontact<6,5,1>;
template class CONTACT::Beam3tosolidcontact<6,2,2>;

template class CONTACT::Beam3tosolidcontact<8,2,1>;
template class CONTACT::Beam3tosolidcontact<8,3,1>;
template class CONTACT::Beam3tosolidcontact<8,4,1>;
template class CONTACT::Beam3tosolidcontact<8,5,1>;
template class CONTACT::Beam3tosolidcontact<8,2,2>;

template class CONTACT::Beam3tosolidcontact<9,2,1>;
template class CONTACT::Beam3tosolidcontact<9,3,1>;
template class CONTACT::Beam3tosolidcontact<9,4,1>;
template class CONTACT::Beam3tosolidcontact<9,5,1>;
template class CONTACT::Beam3tosolidcontact<9,2,2>;
