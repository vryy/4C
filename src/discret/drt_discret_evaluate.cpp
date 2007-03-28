/*!----------------------------------------------------------------------
\file drt_discret_fillcomplete.cpp
\brief

<pre>
Maintainer: Michael Gee
            gee@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15239
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "drt_discret.H"
#include "drt_dserror.H"
#include "linalg_utils.H"
#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_SerialDenseVector.h"




/*----------------------------------------------------------------------*
 |  evaluate (public)                                        mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::Evaluate(
                              ParameterList&                params, 
                              RefCountPtr<Epetra_CrsMatrix> systemmatrix1,
                              RefCountPtr<Epetra_CrsMatrix> systemmatrix2,
                              RefCountPtr<Epetra_Vector>    systemvector1,
                              RefCountPtr<Epetra_Vector>    systemvector2,
                              RefCountPtr<Epetra_Vector>    systemvector3)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  
  // see what we have for input
  bool havesysmatrix1 = false;
  bool havesysmatrix2 = false;
  bool havesysvector1 = false;
  bool havesysvector2 = false;
  bool havesysvector3 = false;
  if (systemmatrix1!=null) havesysmatrix1 = true;
  if (systemmatrix2!=null) havesysmatrix2 = true;
  if (systemvector1!=null) havesysvector1 = true;
  if (systemvector2!=null) havesysvector2 = true;
  if (systemvector3!=null) havesysvector3 = true;
  
  // see what we want to assemble (default is no assembly)
  const bool assemblemat1 = params.get("assemble matrix 1",false);
  const bool assemblemat2 = params.get("assemble matrix 2",false);
  const bool assemblevec1 = params.get("assemble vector 1",false);
  const bool assemblevec2 = params.get("assemble vector 2",false);
  const bool assemblevec3 = params.get("assemble vector 3",false);
  // check whether we have system matrices and vectors supplied to do this
  if (assemblemat1 && !havesysmatrix1) dserror("Do not have system matrix 1 for assembly");
  if (assemblemat2 && !havesysmatrix2) dserror("Do not have system matrix 2 for assembly");
  if (assemblevec1 && !havesysvector1) dserror("Do not have system vector 1 for assembly");
  if (assemblevec2 && !havesysvector2) dserror("Do not have system vector 2 for assembly");
  if (assemblevec3 && !havesysvector3) dserror("Do not have system vector 3 for assembly");

  // define element matrices and vectors
  Epetra_SerialDenseMatrix elematrix1;
  Epetra_SerialDenseMatrix elematrix2;
  Epetra_SerialDenseVector elevector1;
  Epetra_SerialDenseVector elevector2;
  Epetra_SerialDenseVector elevector3;
  
  // loop over column elements
  const int numcolele = NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* actele = lColElement(i);
    
    // get element location vector, dirichlet flags and ownerships
    vector<int> lm;
    vector<int> lmowner;
    actele->LocationVector(lm,lmowner);
    
    // get dimension of element matrices and vectors
    // Reshape element matrices and vectors and init to zero
    const int eledim = (int)lm.size();
    elematrix1.Shape(eledim,eledim); 
    elematrix2.Shape(eledim,eledim); 
    elevector1.Size(eledim);
    elevector2.Size(eledim);
    elevector3.Size(eledim);
    
    // call the element evaluate method
    int err = actele->Evaluate(params,*this,lm,elematrix1,elematrix2,
                               elevector1,elevector2,elevector3);
    if (err) dserror("Proc %d: Element %d returned err=%d",Comm().MyPID(),actele->Id(),err);

    if (assemblemat1) LINALG::Assemble(*systemmatrix1,elematrix1,lm,lmowner);
    if (assemblemat2) LINALG::Assemble(*systemmatrix2,elematrix2,lm,lmowner);
    if (assemblevec1) LINALG::Assemble(*systemvector1,elevector1,lm,lmowner);
    if (assemblevec2) LINALG::Assemble(*systemvector2,elevector2,lm,lmowner);
    if (assemblevec3) LINALG::Assemble(*systemvector3,elevector3,lm,lmowner);

    
  } // for (int i=0; i<numcolele; ++i)
  return;
}

extern "C"
{
  void dyn_facfromcurve(int actcurve,double T,double *fac); 
}
/*----------------------------------------------------------------------*
 |  evaluate Neumann conditions (public)                     mwgee 12/06|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateNeumann(ParameterList& params, Epetra_Vector& systemvector)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  
  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Point Neumann conditions and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"PointNeumann") continue;
    DRT::Condition& cond = *(fool->second);
    const vector<int>* nodeids = cond.Get<vector<int> >("Node Ids");
    if (!nodeids) dserror("PointNeumann condition does not have nodal cloud");
    const int nnode = (*nodeids).size();
    vector<int>*    curve  = cond.Get<vector<int> >("curve");
    vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
    vector<double>* val    = cond.Get<vector<double> >("val");
    // Neumann BCs for some historic reason only have one curve
    int curvenum = -1;
    if (curve) curvenum = (*curve)[0]; 
    double curvefac = 1.0;
      if (curvenum>=0 && usetime)
        dyn_facfromcurve(curvenum,time,&curvefac); 
    for (int i=0; i<nnode; ++i)
    {
      // do only nodes in my row map
      if (!NodeRowMap()->MyGID((*nodeids)[i])) continue;
      DRT::Node* actnode = gNode((*nodeids)[i]);
      if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
      const int  numdf = actnode->Dof().NumDof();
      const int* dofs  = actnode->Dof().Dofs();
      for (int j=0; j<numdf; ++j)
      {
        if ((*onoff)[j]==0) continue;
        const int gid = dofs[j];
        double value  = (*val)[j];
        value *= curvefac;
        const int lid = systemvector.Map().LID(gid);
        if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
        systemvector[lid] += value;
      }
    }
  }
  //--------------------------------------------------------
  // loop through line/surface/volume Neumann BCs and evaluate them
  //--------------------------------------------------------
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
    if (fool->first == (string)"LineNeumann" ||
        fool->first == (string)"SurfaceNeumann" ||
        fool->first == (string)"VolumeNeumann"
       )
    {
      DRT::Condition& cond = *(fool->second);
      map<int,RefCountPtr<DRT::Element> >& geom = cond.Geometry();
      map<int,RefCountPtr<DRT::Element> >::iterator curr;
      Epetra_SerialDenseVector elevector;
      for (curr=geom.begin(); curr!=geom.end(); ++curr)
      {
        // get element location vector, dirichlet flags and ownerships
        vector<int> lm;
        vector<int> lmowner;
        curr->second->LocationVector(lm,lmowner);
        elevector.Size((int)lm.size());
        curr->second->EvaluateNeumann(params,*this,cond,lm,elevector);
        LINALG::Assemble(systemvector,elevector,lm,lmowner);
      }
    }
  return;
}


static void DoDirichletCondition(DRT::Condition&      cond,
                                 DRT::Discretization& dis,
                                 const bool           usetime,
                                 const double         time,
                                 Epetra_Vector&       systemvector,
                                 Epetra_Vector&       toggle);
				 
static double EvaluateFunction(DRT::Node*        node,
		               int               index,
			       int		 funct_num);				 

				
/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DRT::Discretization::EvaluateDirichlet(ParameterList& params, 
                                            Epetra_Vector& systemvector,
                                            Epetra_Vector& toggle)
{
  if (!Filled()) dserror("FillComplete() was not called");
  if (!HaveDofs()) dserror("AssignDegreesOfFreedom() was not called");
  
  // get the current time
  bool usetime = true;
  const double time = params.get("total time",-1.0);
  if (time<0.0) usetime = false;
  
  // make temp. copy of system vector
  Epetra_Vector backup(systemvector);

  multimap<string,RefCountPtr<Condition> >::iterator fool;
  //--------------------------------------------------------
  // loop through Dirichlet conditions and evaluate them
  //--------------------------------------------------------
  // Note that this method does not sum up but 'sets' values in systemvector.
  // For this reason, Dirichlet BCs are evaluated hierarchical meaning
  // in this order:
  //                VolumeDirichlet
  //                SurfaceDirichlet
  //                LineDirichlet
  //                PointDirichlet
  // This way, lower entities override higher ones which is
  // equivalent to inheritance of dirichlet BCs as done in the old
  // ccarat discretization with design          (mgee 1/07)
  
  // Do VolumeDirichlet first
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::VolumeDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do SurfaceDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::SurfaceDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do LineDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::LineDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  // Do PointDirichlet
  for (fool=condition_.begin(); fool!=condition_.end(); ++fool)
  {
    if (fool->first != (string)"Dirichlet") continue;
    if (fool->second->Type() != DRT::Condition::PointDirichlet) continue;
    DoDirichletCondition(*(fool->second),*this,usetime,time,systemvector,toggle);
  }
  
  // copy all values not marked as Dirichlet in toggle from
  // temporary copy back to systemvector
  const int mylength = systemvector.MyLength();
  for (int i=0; i<mylength; ++i)
    if (toggle[i]==0.0)
      systemvector[i] = backup[i];
  
  return;
}


/*----------------------------------------------------------------------*
 |  evaluate Dirichlet conditions (public)                   mwgee 01/07|
 *----------------------------------------------------------------------*/
void DoDirichletCondition(DRT::Condition&      cond,
                          DRT::Discretization& dis,
                          const bool           usetime,
                          const double         time,
                          Epetra_Vector&       systemvector,
                          Epetra_Vector&       toggle)
{
  const vector<int>* nodeids = cond.Get<vector<int> >("Node Ids");
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");
  const int nnode = (*nodeids).size();
  vector<int>*    curve  = cond.Get<vector<int> >("curve");
  vector<int>*    funct  = cond.Get<vector<int> >("funct");
  vector<int>*    onoff  = cond.Get<vector<int> >("onoff");
  vector<double>* val    = cond.Get<vector<double> >("val");

  for (int i=0; i<nnode; ++i)
  {
    // do only nodes in my row map
    if (!dis.NodeRowMap()->MyGID((*nodeids)[i])) continue;
    DRT::Node* actnode = dis.gNode((*nodeids)[i]);
    if (!actnode) dserror("Cannot find global node %d",(*nodeids)[i]);
    const int  numdf = actnode->Dof().NumDof();
    const int* dofs  = actnode->Dof().Dofs();
    for (int j=0; j<numdf; ++j)
    {
      if ((*onoff)[j]==0) 
      {
        const int lid = systemvector.Map().LID(dofs[j]);
        if (lid<0) dserror("Global id %d not on this proc in system vector",dofs[j]);
        toggle[lid] = 0.0;
        continue;
      }
      const int gid = dofs[j];
      double value  = (*val)[j];
      
      // factor given by time curve
      double curvefac = 1.0;
      int    curvenum = -1;
      if (curve) curvenum = (*curve)[j];
      if (curvenum>=0 && usetime)
        dyn_facfromcurve(curvenum,time,&curvefac);
      //cout << "Dirichlet value " << value << " curvefac " <<  curvefac << endl;
            
      // factor given by spatial function
      double functfac = 1.0;
      int funct_num = -1;
      if (funct) funct_num = (*funct)[j];
       {
         if (funct_num>0)
         functfac = EvaluateFunction(actnode,j,funct_num);
       }
      //cout << "Dirichlet value " << value << " functfac " <<  functfac << endl;
	 
      //apply factors to dirichlet value
      value *= (curvefac*functfac);
      
      const int lid = systemvector.Map().LID(gid);
      if (lid<0) dserror("Global id %d not on this proc in system vector",gid);
      systemvector[lid] = value;
      toggle[lid] = 1.0;
    }
  }
  return;
}


// some extern declarations necessary for routine 'EvaluateFunction'
extern "C"
{ 
double pss_evaluate_funct(struct _ST_NODE* funct, double x, double y, double z);
}
#include "../headers/standardtypes.h"
#include "../pss_full/pss_parser.h"

extern struct _FUNCT* funct;
extern int            numfunct; 


/*----------------------------------------------------------------------*
 |  evaluate spatial function (public)                       g.bau 03/07|
 *----------------------------------------------------------------------*/
double EvaluateFunction(DRT::Node*      node,
		        int             index,
			int             funct_num)
{
  double  xi;
  double  length,length_1,length_2;
  double  xp[3],x1[3],x2[3];
  double  fac;
  double  a,d;
  double  h,um;
  
  if (funct_num < 1 || funct_num > numfunct)
    dserror("Unknown function: FUNCT%d",funct_num); 
  
  // function number funct_num has now to be reduced by one
  funct_num -= 1;  
  
  // get coordinates of current node
  xp[0]     = node->X()[0];
  xp[1]     = node->X()[1];
  xp[2]     = node->X()[2];

  //switch to the correct funtion 
  switch (funct[funct_num].functtyp)
  {
    case funct_line_lin:  // linear function
      
      FUNCT_LINE_LIN     *f_line_lin;
       
      f_line_lin = funct[funct_num].typ.funct_line_lin;
      if (f_line_lin==NULL)
	dserror("failed to read function %d",index);

      length = f_line_lin->length;

      // x1: vector along the line 
      x1[0] = f_line_lin->x2[0] - f_line_lin->x1[0];
      x1[1] = f_line_lin->x2[1] - f_line_lin->x1[1];
      x1[2] = f_line_lin->x2[2] - f_line_lin->x1[2];

      // x2: vector from the beginning of the line the the point 
      x2[0] = xp[0] - f_line_lin->x1[0];
      x2[1] = xp[1] - f_line_lin->x1[1];
      x2[2] = xp[2] - f_line_lin->x1[2];

      // length_1 = projection of x2 onto x1 
      length_1 = ( x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] ) / length;
      // length_2 = length of the vector x2 
      length_2 = sqrt( x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2] );

      // check for a point not on the line 
      if ( FABS(length_1 - length_2) > 10e-6 )
        //dswarning(1,6);
	;
	
      // calculate xi and check for a point outside the range of the funct 
      xi =  length_1/length;
      if (xi < 0.0 || xi > 1.0)
        //dswarning(1,5);
	;
      // calculate function value at point p 
      fac = f_line_lin->b + xi * f_line_lin->m;
      break;

  //----------------------------------------------------------------

    case funct_radius_lin:  // linear function 
    
      FUNCT_RADIUS_LIN   *f_radius_lin;
    
      f_radius_lin = funct[funct_num].typ.funct_radius_lin;

      length = f_radius_lin->length;

      // x2: vector from the beginning of the line the the point 
      x1[0] = xp[0] - f_radius_lin->x1[0];
      x1[1] = xp[1] - f_radius_lin->x1[1];
      x1[2] = xp[2] - f_radius_lin->x1[2];

      // length_1 = length of the vector x1 
      length_1 = sqrt( x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2] );

      // calculate xi and check for a point outside the range of the funct 
      xi =  length_1/length;
      if (xi < 0.0 || xi > 1.0)
        //dswarning(1,5);
	;

      // calculate function value at point p 
      fac = f_radius_lin->b + xi * f_radius_lin->m;
      break;

  //----------------------------------------------------------------

    case funct_line_quad: //quadratic parabola 
    
      FUNCT_LINE_QUAD    *f_line_quad;
    
      f_line_quad = funct[funct_num].typ.funct_line_quad;

      length = f_line_quad->length;

      // x1: vector along the line 
      x1[0] = f_line_quad->x2[0] - f_line_quad->x1[0];
      x1[1] = f_line_quad->x2[1] - f_line_quad->x1[1];
      x1[2] = f_line_quad->x2[2] - f_line_quad->x1[2];

      // x2: vector from the beginning of the line the the point 
      x2[0] = xp[0] - f_line_quad->x1[0];
      x2[1] = xp[1] - f_line_quad->x1[1];
      x2[2] = xp[2] - f_line_quad->x1[2];

      // length_1 = projection of x2 onto x1 
      length_1 = ( x1[0]*x2[0] + x1[1]*x2[1] + x1[2]*x2[2] ) / length;
      // length_2 = length of the vector x2 
      length_2 = sqrt( x2[0]*x2[0] + x2[1]*x2[1] + x2[2]*x2[2] );

      // check for a point not on the line 
      if ( FABS(length_1 - length_2) > 10e-6 )
        //dswarning(1,6);
	;

      // calculate xi and check for a point outside the range of the funct 
      xi =  length_1/length;
      if (xi < 0.0 || xi > 1.0)
        //dswarning(1,5);
	;

      // calculate function value at point p 
      fac = 1.0 - 4 * (xi - 1.0/2.0)*(xi - 1.0/2.0);
      break;

  //----------------------------------------------------------------

    case funct_radius_quad: //quadratic parabola 
    
      FUNCT_RADIUS_QUAD  *f_radius_quad;
    
      f_radius_quad = funct[funct_num].typ.funct_radius_quad;

      length = f_radius_quad->length;

      // x1: vector from the beginning of the line the the point 
      x1[0] = xp[0] - f_radius_quad->x1[0];
      x1[1] = xp[1] - f_radius_quad->x1[1];
      x1[2] = xp[2] - f_radius_quad->x1[2];

      // length_1 = length of the vector x1 */
      length_1 = sqrt( x1[0]*x1[0] + x1[1]*x1[1] + x1[2]*x1[2] );

      // calculate xi and check for a point outside the range of the funct 
      xi =  length_1/length;
      if ( xi > 1.0)
        //dswarning(1,5);
	;

      // calculate function value at point p 
      fac = 1.0 - xi * xi ;
      break;

  //----------------------------------------------------------------

    case funct_bel:  // spatial function for beltrami flow 
      // set some constants 
      a    = PI/4.0;
      d    = PI/2.0;
      // calculate values 
      switch (index)
      {
        case 0:
          fac  = -a * ( exp(a*xp[0]) * sin(a*xp[1] + d*xp[2]) +
              exp(a*xp[2]) * cos(a*xp[0] + d*xp[1]) );
          break;
        case 1:
          fac  = -a * ( exp(a*xp[1]) * sin(a*xp[2] + d*xp[0]) +
              exp(a*xp[0]) * cos(a*xp[1] + d*xp[2]) );
          break;
        case 2:
          fac  = -a * ( exp(a*xp[2]) * sin(a*xp[0] + d*xp[1]) +
              exp(a*xp[1]) * cos(a*xp[2] + d*xp[0]) );
          break;
        case 3:
          fac  = -a*a/2 * ( exp(2*a*xp[0]) + exp(2*a*xp[1]) + exp(2*a*xp[2])
              + 2* sin(a*xp[0]+d*xp[1]) * cos(a*xp[2]+d*xp[0]) * exp(a*(xp[1]+xp[2]))
              + 2* sin(a*xp[1]+d*xp[2]) * cos(a*xp[0]+d*xp[1]) * exp(a*(xp[2]+xp[0]))
              + 2* sin(a*xp[2]+d*xp[0]) * cos(a*xp[1]+d*xp[2]) * exp(a*(xp[0]+xp[1])));
          break;
        default:
          fac = 1.0;
          break;
      }
      break;

  //----------------------------------------------------------------

    case funct_kim:  // spatial function for kim-moin flow 
      // set some constants 
      a    = 2.0;
      // calculate values 
      switch (index)
      {
        case 0:
          fac  = - cos(a*PI*xp[0]) * sin(a*PI*xp[1]);
          break;
        case 1:
          fac  = + sin(a*PI*xp[0]) * cos(a*PI*xp[1]);
          break;
        case 2:
          fac  = -1.0/4.0 * ( cos(2.0*a*PI*xp[0]) + cos(2.0*a*PI*xp[1]) );
          break;
        default:
          fac = 1.0;
          break;
      }
      break;

  //----------------------------------------------------------------

    case funct_cyl:
    
      FUNCT_CYL          *f_cyl;  
      f_cyl = funct[funct_num].typ.funct_cyl;

      // set some constants 
      h    = 0.41;
      um = f_cyl->um;

      // calculate values 
      fac = 16*um*xp[1]*xp[2]*(h-xp[1])*(h-xp[2]) / (h*h*h*h);
      break;

  //----------------------------------------------------------------

  case funct_explicit:
  {
    FUNCT_EXPLICIT* f = funct[funct_num].typ.funct_explicit;
    fac = pss_evaluate_funct(f->funct,
		       xp[0] - f->x[0],
		       xp[1] - f->x[1],
		       xp[2] - f->x[2]);
    break;
  }


    default:  // default: no function 
      fac = 1.0;
      break;
  } // end of switch (funct[funct_num].functtyp)


  return fac;

} // end of EvaluateFunction	




#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
