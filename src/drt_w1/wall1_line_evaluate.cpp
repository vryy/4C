/*!----------------------------------------------------------------------
\file wall1_line_evaluate.cpp
\brief

<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/
#ifdef D_WALL1
#ifdef CCADISCRET

#include "wall1.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_dserror.H"
#include "../drt_lib/drt_timecurve.H"


/*----------------------------------------------------------------------*
 |  Integrate a Line Neumann boundary condition (public)     mgit 03/07|
 *----------------------------------------------------------------------*/

int DRT::ELEMENTS::Wall1Line::EvaluateNeumann(ParameterList& params,
                              DRT::Discretization&      discretization,
                              DRT::Condition&           condition,
                              vector<int>&              lm,
                              Epetra_SerialDenseVector& elevec1)
{
  RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
  if (disp==null) dserror("Cannot get state vector 'displacement'");
  vector<double> mydisp(lm.size());
  DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);

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
    curvefac = DRT::UTILS::TimeCurveManager::Instance().Curve(curvenum).f(time);

   // set number of nodes
  const int iel   = NumNode();
  
  const DiscretizationType distype = Shape();
  
  // gaussian points 
  const DRT::UTILS::GaussRule1D gaussrule = getOptimalGaussrule(distype); 
  const DRT::UTILS::IntegrationPoints1D  intpoints = getIntegrationPoints1D(gaussrule);
  
  // get values and switches from the condition
  const vector<int>*    onoff = condition.Get<vector<int> >("onoff");
  const vector<double>* val   = condition.Get<vector<double> >("val");
   
  // allocate vector for shape functions and for derivatives
  Epetra_SerialDenseVector    funct(iel);
  Epetra_SerialDenseMatrix    deriv(1,iel);
  
  // node coordinates
  Epetra_SerialDenseMatrix xye(2,iel);
    
  // get node coordinates
  for(int i=0;i<iel;++i)
  {
    xye(0,i)=Nodes()[i]->X()[0]; 
    xye(1,i)=Nodes()[i]->X()[1];
  }
   
  // loop over integration points //new
  for (int gpid=0;gpid<intpoints.nquad;gpid++)
  {  
     const double e1 = intpoints.qxg[gpid];
     
  // get shape functions and derivatives in the line
     DRT::UTILS::shape_function_1D(funct,e1,distype);
     DRT::UTILS::shape_function_1D_deriv1(deriv,e1,distype);
  
  // compute infinitesimal line element dr for integration along the line
     const double dr = w1_substitution(xye,deriv,iel);
    
  // load vector ar
     double ar[NODDOF_W1];
  // loop the dofs of a node
  // ar[i] = ar[i] * facr * ds * onoff[i] * val[i]*curvefac
     for (int i=0; i<NODDOF_W1; ++i)
     {
        ar[i] = intpoints.qwgt[gpid] * dr * (*onoff)[i]*(*val)[i] * curvefac;
     }

   // add load components
     for (int node=0; node<iel; ++node)
     {	 
       for (int j=0; j<NODDOF_W1; ++j)
       {
          elevec1[node*NODDOF_W1+j] += funct[node]*ar[j];
       }
     }
  } 
  return 0;
}


DRT::UTILS::GaussRule1D DRT::ELEMENTS::Wall1Line::getOptimalGaussrule(const DiscretizationType& distype)
{
  DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
  switch (distype)
    {
    case line2:
      rule = DRT::UTILS::intrule_line_2point;
      break;
    case line3:
      rule = DRT::UTILS::intrule_line_3point;
      break;
    default: 
    dserror("unknown number of nodes for gaussrule initialization");
    }
  return rule;
}

// determinant of jacobian matrix

double  DRT::ELEMENTS::Wall1Line::w1_substitution(const Epetra_SerialDenseMatrix& xye,
                                                  const Epetra_SerialDenseMatrix& deriv,
                                                  const int iel)
{
	  /*
	  |                                            0 1 
	  |                                           +-+-+
	  |       0 1              0...iel-1          | | | 0
	  |      +-+-+             +-+-+-+-+          +-+-+
	  |      | | | 1     =     | | | | | 0        | | | .
	  |      +-+-+             +-+-+-+-+       *  +-+-+ .
	  |                                           | | | .
	  |                                           +-+-+
	  |                                           | | | iel-1
	  |		           	     	                  +-+-+
	  |
	  |       dxyzdrs             deriv^T          xye^T
	  |
	  |
	  |                       +-        -+
	  |  	   	    	    	    | dx   dy  |
	  |  	  yields   dxydr =	| --   --  |
	  | 	                    | dr   dr  |
	  |                       +-        -+
	  |
	  */
  // compute derivative of parametrization
  double dr = 0.0;
  Epetra_SerialDenseMatrix der_par (1,2);
  int err = der_par.Multiply('N','T',1.0,deriv,xye,0.0);
  if (err!=0)
    dserror("Multiply failed");
  dr=sqrt(der_par(0,0)*der_par(0,0)+der_par(0,1)*der_par(0,1));
  return dr;
}


int DRT::ELEMENTS::Wall1Line::Evaluate(ParameterList& params,
                                DRT::Discretization&      discretization,
                                vector<int>&              lm,
                                Epetra_SerialDenseMatrix& elematrix1,
                                Epetra_SerialDenseMatrix& elematrix2,
                                Epetra_SerialDenseVector& elevector1,
                                Epetra_SerialDenseVector& elevector2,
                                Epetra_SerialDenseVector& elevector3)
{
  const DiscretizationType distype = Shape();

  // start with "none"
  DRT::ELEMENTS::Wall1Line::ActionType act = Wall1Line::none;

  // get the required action
  string action = params.get<string>("action","none");
  if (action == "none") dserror("No action supplied");
  else if (action=="calc_struct_constrarea")       act = Wall1Line::calc_struct_constrarea;
  else if (action=="calc_struct_areaconstrstiff")  act= Wall1Line::calc_struct_areaconstrstiff;
  else dserror("Unknown type of action for Soh8Surface");
  //create communicator
  const Epetra_Comm& Comm = discretization.Comm();
  // what the element has to do
  switch(act)
  {
    //just compute the enclosed volume (e.g. for initialization)
    case calc_struct_constrarea:
    {
      if (distype!=line2)
      {
        dserror("Area Constraint only works for line2 curves!");
      }
      //We are not interested in volume of ghosted elements
      if(Comm.MyPID()==Owner())
      {
        // element geometry update
        RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
        if (disp==null) dserror("Cannot get state vector 'displacement'");
        vector<double> mydisp(lm.size());
        DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
        const int numnod = NumNode();
        Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_W1);  // material coord. of element
        Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_W1);  // material coord. of element
        for (int i=0; i<numnod; ++i)
        {
          xsrefe(i,0) = Nodes()[i]->X()[0];
          xsrefe(i,1) = Nodes()[i]->X()[1];
          
          xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_W1];
          xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_W1+1];
        }
        //compute area between line and x-Axis
        double areaele =  0.5*(xscurr(0,1)+xscurr(1,1))*(xscurr(1,0)-xscurr(0,0));
        
        // get RIGHT area out of parameterlist and maximum ConditionID
        char areaname[30];
        const int ID =params.get("ConditionID",-1);
        const int maxID=params.get("MaxID",0);
        const int minID=params.get("MinID",1000000);
        if (ID<0)
        {
          dserror("Condition ID for area constraint missing!");
        }
        if (maxID<ID)
        {
          params.set("MaxID",ID);
        }
        if (minID>ID)
        {
          params.set("MinID",ID);
        }
        sprintf(areaname,"computed area %d",ID);
        double areacondval = params.get(areaname,0.0);
        //update volume in parameter list
        params.set(areaname, areacondval+areaele);
      }
      
    }
    break;
    case calc_struct_areaconstrstiff:
    {
      if (distype!=line2)
      {
        dserror("Area Constraint only works for line2 curves!");
      }  // element geometry update
      RefCountPtr<const Epetra_Vector> disp = discretization.GetState("displacement");
      if (disp==null) 
      {
        dserror("Cannot get state vector 'displacement'");
      }
      vector<double> mydisp(lm.size());
      DRT::UTILS::ExtractMyValues(*disp,mydisp,lm);
      const int numnod = NumNode();
      Epetra_SerialDenseMatrix xsrefe(numnod,NUMDIM_W1);  // material coord. of element
      Epetra_SerialDenseMatrix xscurr(numnod,NUMDIM_W1);  // material coord. of element
      for (int i=0; i<numnod; ++i)
      {
        xsrefe(i,0) = Nodes()[i]->X()[0];
        xsrefe(i,1) = Nodes()[i]->X()[1];
        
        xscurr(i,0) = xsrefe(i,0) + mydisp[i*NODDOF_W1];
        xscurr(i,1) = xsrefe(i,1) + mydisp[i*NODDOF_W1+1];
      }
      //call submethods
      ComputeAreaConstrStiff(xscurr,elematrix1);
      ComputeAreaConstrDeriv(xscurr,elevector1);
      //apply the right lagrange multiplier and right signs to matrix and vectors
      const int ID =params.get("ConditionID",-1);
      RCP<Epetra_Vector> lambdav=rcp(new Epetra_Vector(*(params.get<RCP<Epetra_Vector> >("LagrMultVector"))));
      if (ID<0)
      {
        dserror("Condition ID for area constraint missing!");
      }
      const int minID =params.get("MinID",0);
      //update corresponding column in "constraint" matrix
      elevector2=elevector1;
      elevector1.Scale(1*(*lambdav)[ID-minID]);
      elematrix1.Scale(1*(*lambdav)[ID-minID]);
      //compute area between line and x-Axis
      if(Comm.MyPID()==Owner())
      {
        double areaele =  0.5*(xscurr(0,1)+xscurr(1,1))*(xscurr(1,0)-xscurr(0,0));

        // get RIGHT area out of parameterlist
        char areaname[30];
        sprintf(areaname,"computed area %d",ID);
        double areacondval = params.get(areaname,0.0);
        //update volume in parameter list
        params.set(areaname, areacondval+areaele);
      }
    }
    break;
    default:
      dserror("Unimplemented type of action for Soh8Surface");

   }
   return 0;
}

/*----------------------------------------------------------------------*
 * Compute first derivatives of area                            tk 10/07*
 * with respect to the displacements                                    *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::ComputeAreaConstrDeriv(Epetra_SerialDenseMatrix xscurr,
    Epetra_SerialDenseVector& elevector)
{
  if (elevector.Length()!=4)
  {
    cout<<"Length of element Vector: "<<elevector.Length()<<endl;
    dserror("That is not the right size!");
  }
  //implementation of simple analytic solution
  elevector[0]=-xscurr(0,1)-xscurr(1,1);
  elevector[1]=xscurr(1,0)-xscurr(0,0);
  elevector[2]=xscurr(0,1)+xscurr(1,1);
  elevector[3]=xscurr(1,0)-xscurr(0,0);
  elevector.Scale(-0.5);
  return ;
}

/*----------------------------------------------------------------------*
 * Compute influence of area constraint on stiffness matrix.    tk 10/07*
 * Second derivatives of areas with respect to the displacements        *
 * ---------------------------------------------------------------------*/
void DRT::ELEMENTS::Wall1Line::ComputeAreaConstrStiff(Epetra_SerialDenseMatrix xscurr,
    Epetra_SerialDenseMatrix& elematrix)
{
  elematrix(0,0)=0.0;
  elematrix(0,1)=-0.5;
  elematrix(0,2)=0.0;
  elematrix(0,3)=-0.5;

  elematrix(1,0)=-0.5;
  elematrix(1,1)=0.0;
  elematrix(1,2)=0.5;
  elematrix(1,3)=0.0;
  
  elematrix(2,0)=0.0;
  elematrix(2,1)=0.5;
  elematrix(2,2)=0.0;
  elematrix(2,3)=0.5;

  elematrix(3,0)=-0.5;
  elematrix(3,1)=0.0;
  elematrix(3,2)=0.5;
  elematrix(3,3)=0.0;

  elematrix.Scale(-1.0);
  return ;
}




#endif  // #ifdef CCADISCRET
#endif // #ifdef D_WALL1
