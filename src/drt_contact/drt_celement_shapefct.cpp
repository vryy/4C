/*!----------------------------------------------------------------------
\file drt_celement_shapefct.cpp
\brief

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15264
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include "drt_celement.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/linalg_utils.H"
#include "drt_contact_integrator.H"
#include "contactdefines.H"


/*----------------------------------------------------------------------*
 |  1D/2D shape function repository                           popp 04/08|
 *----------------------------------------------------------------------*/
void CONTACT::CElement::ShapeFunctions(CElement::ShapeType shape,
                                       const double* xi,
                                       vector<double>& val,
                                       vector<double>& deriv)
{
  switch(shape)
  {
  // *********************************************************************
  // 1D standard linear shape functions
  // (used for interpolation of displacemt field)
  // *********************************************************************
  case CElement::lin1D:
  {
    val[0] = 0.5*(1-xi[0]);
    val[1] = 0.5*(1+xi[0]); 
    deriv[0] = -0.5;
    deriv[1] =  0.5;
    break;
  }
  // *********************************************************************
  // 1D dual linear shape functions
  // (used for interpolation of Lagrange mutliplier field)
  // *********************************************************************
  case CElement::lindual1D:
  {
    val[0] = 0.5*(1-3*xi[0]);
    val[1] = 0.5*(1+3*xi[0]);
    deriv[0] = -1.5;
    deriv[1] =  1.5;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (const replacing linear)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case CElement::lindual1D_edge0:
  {
    val[0] = 0;
    val[1] = 1;
    deriv[0] = 0;
    deriv[1] = 0;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (const replacing linear)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // *********************************************************************
  case CElement::lindual1D_edge1:
  {
    val[0] = 1;
    val[1] = 0;
    deriv[0] = 0;
    deriv[1] = 0;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear replacing quad)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (only form a basis and have to be adapted for distorted elements)
  // *********************************************************************
  case CElement::dual1D_base_for_edge0:
  {
    val[0] = xi[0];
    val[1] = 1-xi[0]; 
    deriv[0] =  1;
    deriv[1] = -1;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear replacing quad)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (only form a basis and have to be adapted for distorted elements)
  // *********************************************************************
  case CElement::dual1D_base_for_edge1:
  {
    val[0] = -xi[0];
    val[1] = 1+xi[0]; 
    deriv[0] = -1;
    deriv[1] =  1;
    break;
  }
  // *********************************************************************
  // 1D standard quadratic shape functions
  // (used for interpolation of displacemt field)
  // *********************************************************************
  case CElement::quad1D:
  {
    val[0] = 0.5*xi[0]*(xi[0]-1);
    val[1] = 0.5*xi[0]*(xi[0]+1);
    val[2] = (1-xi[0])*(1+xi[0]);       
    deriv[0] = xi[0]-0.5;
    deriv[1] = xi[0]+0.5;
    deriv[2] = -2*xi[0];
    break;
  }
  // *********************************************************************
  // 1D dual quadratic shape functions
  // (used for interpolation of displacemt field)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case CElement::quaddual1D:
  {
    // establish fundamental data  
    double detg = 0.0;
    int nnodes = NumNode();
  
    LINALG::SerialDenseMatrix coord(3,nnodes);
    coord = GetNodalCoords();
    
    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    CONTACT::Integrator integrator(CONTACTNGP,true);
    
    Epetra_SerialDenseMatrix me(nnodes,nnodes);
    Epetra_SerialDenseMatrix de(nnodes,nnodes);
    
    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i), 0.0};
      EvaluateShape1D(gpc, val, deriv, nnodes);
      detg = Jacobian1D(val,deriv,coord);
      
      for (int j=0;j<nnodes;++j)
        for (int k=0;k<nnodes;++k)
        {
          me(j,k)+=integrator.Weight(i)*val[j]*val[k]*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*val[j]*detg;
        }  
    }
    
    // invert bi-ortho matrix me
    LINALG::SymmetricInverse(me,nnodes);
    
    // get solution matrix with dual parameters
    Epetra_SerialDenseMatrix ae(nnodes,nnodes);
    ae.Multiply('N','N',1.0,de,me,0.0);
    
    // evaluate dual shape functions at loc. coord. xi
    // need standard shape functions at xi first
    EvaluateShape1D(xi, val, deriv, nnodes);
    
    vector<double> valtemp(nnodes);
    vector<double> derivtemp(nnodes);
    for (int i=0;i<nnodes;++i)
    {
      valtemp[i]=0.0;
      derivtemp[i]=0.0;
      for (int j=0;j<nnodes;++j)
      {
        valtemp[i]+=ae(i,j)*val[j];
        derivtemp[i]+=ae(i,j)*deriv[j];
      }
    }
    val=valtemp;
    deriv=derivtemp;
    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case CElement::quaddual1D_edge0:
  {
    // establish fundamental data  
    double detg = 0.0;
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    coord = GetNodalCoords();
    
    // empty shape function vals + derivs
    vector<double> valquad(nnodes);
    vector<double> derivquad(nnodes);
    vector<double> vallin(nnodes-1);
    vector<double> derivlin(nnodes-1);
    vector<double> valtemp(nnodes);
    vector<double> derivtemp(nnodes);
    
    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    CONTACT::Integrator integrator(CONTACTNGP,true);
    
    Epetra_SerialDenseMatrix me(nnodes-1,nnodes-1);
    Epetra_SerialDenseMatrix de(nnodes-1,nnodes-1);
    
    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i), 0.0};
      ShapeFunctions(CElement::quad1D,gpc,valquad,derivquad);
      ShapeFunctions(CElement::dual1D_base_for_edge0,gpc,vallin,derivlin);
      detg = Jacobian1D(valquad,derivquad,coord);
      
      for (int j=1;j<nnodes;++j)
        for (int k=1;k<nnodes;++k)
        {
          me(j-1,k-1)+=integrator.Weight(i)*vallin[j-1]*valquad[k]*detg;
          de(j-1,k-1)+=(j==k)*integrator.Weight(i)*valquad[k]*detg;
        }  
    }

    // invert bi-ortho matrix me
    // CAUTION: This is a non-symmetric inverse operation!
    double detme = me(0,0)*me(1,1)-me(0,1)*me(1,0);
    Epetra_SerialDenseMatrix meold(nnodes-1,nnodes-1);
    meold=me;
    me(0,0) =  1/detme*meold(1,1);
    me(0,1) = -1/detme*meold(0,1);
    me(1,0) = -1/detme*meold(1,0);
    me(1,1) =  1/detme*meold(0,0);

    // get solution matrix with dual parameters
    Epetra_SerialDenseMatrix ae(nnodes-1,nnodes-1);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    ShapeFunctions(CElement::dual1D_base_for_edge0,xi,vallin,derivlin);
    for (int i=1;i<nnodes;++i)
      for (int j=1;j<nnodes;++j)
      {
        valtemp[i]+=ae(i-1,j-1)*vallin[j-1];
        derivtemp[i]+=ae(i-1,j-1)*derivlin[j-1];
      }
    
    val[0] = 0;
    val[1] = valtemp[1];
    val[2] = valtemp[2];
    deriv[0] =  0;
    deriv[1] = derivtemp[1];
    deriv[2] = derivtemp[2];

    break;
  }
  // *********************************************************************
  // 1D modified dual shape functions (linear)
  // (used for interpolation of Lagrange mult. field near boundaries)
  // (including adaption process for distorted elements)
  // *********************************************************************
  case CElement::quaddual1D_edge1:
  {
    // establish fundamental data  
    double detg = 0.0;
    int nnodes = NumNode();
    LINALG::SerialDenseMatrix coord(3,nnodes);
    coord = GetNodalCoords();
    
    // empty shape function vals + derivs
    vector<double> valquad(nnodes);
    vector<double> derivquad(nnodes);
    vector<double> vallin(nnodes-1);
    vector<double> derivlin(nnodes-1);
    vector<double> valtemp(nnodes);
    vector<double> derivtemp(nnodes);
    
    // compute entries to bi-ortho matrices me/de with Gauss quadrature
    CONTACT::Integrator integrator(CONTACTNGP,true);
    
    Epetra_SerialDenseMatrix me(nnodes-1,nnodes-1);
    Epetra_SerialDenseMatrix de(nnodes-1,nnodes-1);
    
    for (int i=0;i<integrator.nGP();++i)
    {
      double gpc[2] = {integrator.Coordinate(i), 0.0};
      ShapeFunctions(CElement::quad1D,gpc,valquad,derivquad);
      ShapeFunctions(CElement::dual1D_base_for_edge1,gpc,vallin,derivlin);
      detg = Jacobian1D(valquad,derivquad,coord);
      
      for (int j=0;j<nnodes-1;++j)
        for (int k=0;k<nnodes-1;++k)
        {
          me(j,k)+=integrator.Weight(i)*vallin[j]*valquad[2*k]*detg;
          de(j,k)+=(j==k)*integrator.Weight(i)*valquad[2*k]*detg;
        }  
    }
    
    // invert bi-ortho matrix me
    // CAUTION: This is a non-symmetric inverse operation!
    double detme = me(0,0)*me(1,1)-me(0,1)*me(1,0);
    Epetra_SerialDenseMatrix meold(nnodes-1,nnodes-1);
    meold=me;
    me(0,0) =  1/detme*meold(1,1);
    me(0,1) = -1/detme*meold(0,1);
    me(1,0) = -1/detme*meold(1,0);
    me(1,1) =  1/detme*meold(0,0);

    // get solution matrix with dual parameters
    Epetra_SerialDenseMatrix ae(nnodes-1,nnodes-1);
    ae.Multiply('N','N',1.0,de,me,0.0);

    // evaluate dual shape functions at loc. coord. xi
    ShapeFunctions(CElement::dual1D_base_for_edge1,xi,vallin,derivlin);
    for (int i=0;i<nnodes-1;++i)
      for (int j=0;j<nnodes-1;++j)
      {
        valtemp[2*i]+=ae(i,j)*vallin[j];
        derivtemp[2*i]+=ae(i,j)*derivlin[j];
      }
    
    val[0] = valtemp[0];
    val[1] = 0;
    val[2] = valtemp[2];
    deriv[0] = derivtemp[0];
    deriv[1] = 0;
    deriv[2] = derivtemp[2];
    
    break;
  }
  // *********************************************************************
  // Unkown shape function type
  // *********************************************************************
  default:
    dserror("ERROR: Unknown shape function type identifier");
  }
  
  return;
}

/*----------------------------------------------------------------------*
 |  Evaluate shape functions - LINEAR / QUAD 1D               popp 01/08|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::EvaluateShape1D(const double* xi, vector<double>& val,
                                        vector<double>& deriv, const int valdim)
{
  if (!xi)
    dserror("ERROR: EvaluateShape1D called with xi=NULL");
  
  // 2D linear case (2noded line element)
  if ((valdim==2)&& (Shape()==line2))
    ShapeFunctions(CElement::lin1D,xi,val,deriv);

  // 2D quadratic case (3noded line element)
  else if ((valdim==3) && (Shape()==line3))
    ShapeFunctions(CElement::quad1D,xi,val,deriv);
  
  // unknown case
  else
    dserror("ERROR: EvaluateShape1D called for unknown CElement type");

  return true;
}

/*----------------------------------------------------------------------*
 |  Evaluate dual shape functions - LINEAR / QUAD 1D          popp 12/07|
 *----------------------------------------------------------------------*/
bool CONTACT::CElement::EvaluateShapeDual1D(const double* xi, vector<double>& val,
                                            vector<double>& deriv, const int valdim)
{
  if (!xi)
    dserror("ERROR: EvaluateShapeDual1D called with xi=NULL");
  if (!IsSlave())
    dserror("ERROR: EvaluateShapeDual1D called for master element");
  
  // get node number and node pointers
  DRT::Node** mynodes = Nodes();
  if (!mynodes) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
   
  // 2D linear case (2noded line element)
  if ((valdim==2) && (Shape()==line2))
  {
    // check for boundary nodes
    CNode* mycnode0 = static_cast<CNode*> (mynodes[0]);
    CNode* mycnode1 = static_cast<CNode*> (mynodes[1]);
    if (!mycnode0) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    if (!mycnode1) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    bool isonbound0 = mycnode0->IsOnBound();
    bool isonbound1 = mycnode1->IsOnBound();
    
    // both nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1)
      ShapeFunctions(CElement::lindual1D,xi,val,deriv);
    
    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1)
      ShapeFunctions(CElement::lindual1D_edge0,xi,val,deriv);
     
    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1)
      ShapeFunctions(CElement::lindual1D_edge1,xi,val,deriv);
      
    // both nodes are on boundary: infeasible case
    else
      dserror("ERROR: EvaluateShapeDual1D: Element with 2 boundary nodes!");
  }
  
  // 2D quadratic case (3noded line element)
  else if ((valdim==3) && (Shape()==line3))
  {
    // check for boundary nodes
    CNode* mycnode0 = static_cast<CNode*> (mynodes[0]);
    CNode* mycnode1 = static_cast<CNode*> (mynodes[1]);
    CNode* mycnode2 = static_cast<CNode*> (mynodes[2]);
    if (!mycnode0) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    if (!mycnode1) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    if (!mycnode2) dserror("ERROR: EvaluateShapeDual1D: Null pointer!");
    bool isonbound0 = mycnode0->IsOnBound();
    bool isonbound1 = mycnode1->IsOnBound();
    bool isonbound2 = mycnode2->IsOnBound();
    
    // all 3 nodes are interior: use unmodified dual shape functions
    if (!isonbound0 && !isonbound1 && !isonbound2)
      ShapeFunctions(CElement::quaddual1D,xi,val,deriv);
    
    // node 0 is on boundary: modify dual shape functions
    else if (isonbound0 && !isonbound1 && !isonbound2)
      ShapeFunctions(CElement::quaddual1D_edge0,xi,val,deriv);
     
    // node 1 is on boundary: modify dual shape functions
    else if (!isonbound0 && isonbound1 && !isonbound2)
      ShapeFunctions(CElement::quaddual1D_edge1,xi,val,deriv);
    
    // node 2 is on boundary: infeasible case
    else if (isonbound2)
      dserror("ERROR: EvlautaeShapeDual1D: Middle boundary node");
    
    // nodes 0 and 1 are on boundary: infeasible case
    else 
      dserror("ERROR: EvaluateShapeDual1D: Element with 2 boundary nodes");
  }
  
  // unknown case
  else
    dserror("ERROR: EvaluateShapeDual1D called for unknown element type");
  
  return true;
}

#endif  // #ifdef CCADISCRET
