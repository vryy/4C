/*!----------------------------------------------------------------------
\file drt_constraint_manager.H

\class DRT::ConstrManager

\brief Class controlling constraints and containing the necessary data

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_CONSTRAINT_MANAGER_H
#define DRT_CONSTRAINT_MANAGER_H

#include "Teuchos_RefCountPtr.hpp"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"


#include "Teuchos_ParameterList.hpp"
#include "./drt_discret.H"
#include "./linalg_solver.H"
#include "./linalg_utils.H"


using namespace std;
using namespace Teuchos;



class ConstrManager
{
public:

	class StruGenAlpha;
	//! @name Constructors and destructors and related methods

	  /*!
	  \brief Standard Constructor
	  
	  \param disc (in): Discretization
	  \param disp (in): current displacement to compute initial constraint volume 

	  */
	
	  ConstrManager(DRT::Discretization& disc,
			  RCP<Epetra_Vector> disp);
	  
	  /*!
	  	  \brief Destructor

	 */
	  ~ConstrManager(){return;};

	  //@}
	  
	  //! @name Public Member Functions
	  /*!
	    \brief Change stiffness matrix and force vector according to the volume constraint.
	    Values of lagrange multiplier are taken from intern variable. 
	    Difference between current and prescribed volume is calculated and stored as well.
	    
	    \param time (in):        time at end of time step
	    \param disp (in):        Displacement at end of time step
	    \param fint (in/out):	 Vector of internal forces
	    \param stiff (in/out):   Stiffness matrix
	  */  
	  void StiffnessAndInternalForces(const double time, 
			  RCP<Epetra_Vector> disp,
			  RefCountPtr<Epetra_Vector> fint,
			  RefCountPtr<Epetra_CrsMatrix> stiff);
	  
	  /*!
	   \brief Return norm of difference between actual and prescribed volumes
	  */
	  double GetErrorNorm()
	  {
		  return constrainterr_->Norm2();
	  };
	  	  
	  /*!
	  	   \brief Return number of constrained volumes
	  */ 
	  int GetNumberOfConstraints()
	  {
		  return numConstrID_; 
	  };
	  
	  /*!
	  	   \brief Return DofRowMap of volume i
	  */
	  Epetra_Map GetDofMap(int i)
	  {
		  return *(voldofrowmaps_[i]);
	  };
	  
	  /*!
	  	   \brief Scale increment of lagrange multiplier increment by a double d
	  */
	  void ScaleLagrIncr(double d)
	  {
		  lagrMultInc_->Scale(d);
		  return;
	  };
	  
	  /*!
	  	  	   \brief Scale lagrange multiplier increment by a double d
	  */
	  void ScaleLagrMult(double d)
	  	  {
	  		  lagrMultVec_->Scale(d);
	  		  return;
	  	  };
  
	  /*!
	  	   \brief Update increment of lagrange multiplier Delta lambda_{n+1}=Delta lambda_{n}+factor*(vect+volerr)
	  */
	  void UpdateLagrIncr(double factor, Epetra_SerialDenseVector vect);
	  
	  /*!
	  	   \brief Update lagrange multiplier lambda_{n+1}=lambda_{n}+factor*(volerr)
	  */	  
	  void UpdateLagrMult(double factor);
	  /*!
	  	   \brief Return lagrange multiplier increment
	  */
	  double GetLagrIncr(int i)
	  {
		  return (*lagrMultInc_)[i];
	  };
	  
	  /*!
	  	  	   \brief Compute difference between current and prescribed volume at a given time and a given displacement 
	  */
	  void ComputeError(double time,RCP<Epetra_Vector> disp);
	  
	  /*!
	  	   \brief Return differences between prescribed and actual value of volume number i
	  */
	  double GetError(int i)
	  {
		  return (*constrainterr_)[i];
	  }
	  
	  /*!
	  	   \brief Add Increment to lagrange multiplier
	  */
	  void UpdateLagrMult();
	  
	  
	  /*!
	   \brief Return matrix containing first derivatives of constraints wrt displacements 
	   */
	  RCP<Epetra_MultiVector> GetConstraintMatrix()
	  {
		  return constrVec_;
	  }
	  
	  /*!
	  	  	   \brief Return lagrange multiplier for volume i
	  */
	  double GetLagrMult(int i)
	  {
		  return (*lagrMultVec_)[i];		  
	  }	  
	  
	  bool HaveConstraint()
	  {
		  return haveconstraint_;
	  }
	  
	  /*!
	   \brief Return current value of volume i
	  */
	  double GetCurrValue(int i)
	  {
		  return (*actvalues_)[i];
	  }
	  
	  /*!
	   \brief Computes inner product of constr multivector times lagrange multiplier increment 
	   vector as required for the uzawa algorithm to compute new displacement increment 
	   
	   \param dotprod (out):		resulting value of B^T * Delta lambda 
	   */
	  void ComputeConstrTimesLagrIncr(RCP<Epetra_Vector> dotprod);
	  
	  /*!
	    \brief Return inner product of constr multivector times displacement increment
	    
	    \param disi (in)		:		displacement increment Delta u
	    \param dotprod (out) 	:		resulting value of B * Delta u
	    */
	  void ComputeConstrTimesDisi(Epetra_Vector disi,RCP<Epetra_SerialDenseVector> dotprod);
	  //@}
private:
	
	  //! name@ Private Member Functions 	
	  //methods used by constructor and other public functions
	  /*!
	    \brief synchronize processors after evaluating constraint volume 

	  */  
	
	  void SynchronizeConstraint(ParameterList& params,
				RCP<Epetra_SerialDenseVector>& vect,const char* resultstr);
	  
	  /*!
	      \brief Setup separate dofrowmaps for any condition ID and initialize lagrange multiplier

	  */  
	  
	  void ManageIDs(ParameterList& params);
	  //@}
	
private:
	
	// don't want = operator, cctor and destructor
	
	ConstrManager operator = (const ConstrManager& old);
	ConstrManager(const ConstrManager& old);    

	
	//! name@ Private Attributes
	DRT::Discretization& actdisc_;	
	RefCountPtr<Epetra_SerialDenseVector> referencevalues_;
	RefCountPtr<Epetra_SerialDenseVector> initialvalues_;
	RefCountPtr<Epetra_SerialDenseVector> actvalues_;
	RefCountPtr<Epetra_SerialDenseVector> constrainterr_;
	int minConstrID_;
	int maxConstrID_;
	int numConstrID_;
	double fact_; 
	RefCountPtr<Epetra_SerialDenseVector> lagrMultVec_;
	RefCountPtr<Epetra_SerialDenseVector> lagrMultInc_;
	map<int, RCP<Epetra_Map> > voldofrowmaps_;
	RefCountPtr<Epetra_MultiVector> constrVec_;
	bool haveconstraint_;
	//@}
    

}; //class

#endif /*CONSTRAINTMANAGER_H*/
