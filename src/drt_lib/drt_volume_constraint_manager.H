/*!----------------------------------------------------------------------
\file drt_volume_constraint_manager.H

\class DRT::VolConstrManager

\brief Class controlling volume constraint and containing the necessary data

<pre>
Maintainer: Thomas Kloeppel
            kloeppel@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_VOLUME_CONSTRAINT_MANAGER_H
#define DRT_VOLUME_CONSTRAINT_MANAGER_H

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


namespace DRT
{
class VolConstrManager
{
public:

	class StruGenAlpha;
	//! @name Constructors and destructors and related methods

	  /*!
	  \brief Standard Constructor
	  
	  \param disc (in): Discretization
	  \param disp (in): current displacement to compute initial constraint volume 

	  */
	
	  VolConstrManager(DRT::Discretization& disc,
			  RCP<Epetra_Vector> disp);
	  
	  /*!
	  	  \brief Destructor

	 */
	  ~VolConstrManager(){return;};

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
	  double GetVolumeErrorNorm()
	  {
		  return volerr_->Norm2();
	  };
	  	  
	  /*!
	  	   \brief Return number of constrained volumes
	  */ 
	  int GetNumberOfVolumes()
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
	  void ComputeVolumeError(double time,RCP<Epetra_Vector> disp);
	  
	  /*!
	  	   \brief Return differences between prescribed and actual value of volume number i
	  */
	  double GetVolumeError(int i)
	  {
		  return (*volerr_)[i];
	  }
	  
	  /*!
	  	   \brief Return additional vector containing constraints
	  */
	  RCP<Epetra_Vector> GetConstrVec()
	  {
		  return constrVec_;
	  };
	  
	  /*!
	  	   \brief Add Increment to lagrange multiplier
	  */
	  void UpdateLagrMult();
	  
	  /*!
	  	  	   \brief Return lagrange multiplier for volume i
	  */
	  double GetLagrMult(int i)
	  {
		  return (*lagrMultVec_)[i];		  
	  }	  
	  
	  /*!
	   \brief Return current value of volume i
	  */
	  double GetCurrVol(int i)
	  {
		  return (*actvol_)[i];
	  }
	  //@}
private:
	
	  //! name@ Private Member Functions 	
	  //methods used by constructor and other public functions
	  /*!
	    \brief synchronize processors after evaluating constraint volume 

	  */  
	
	  void SynchronizeVolConstraint(ParameterList& params,
				RCP<Epetra_SerialDenseVector>& vect);
	  
	  /*!
	      \brief Setup separate dofrowmaps for any condition ID and initialize lagrange multiplier

	  */  
	  void SetupVolDofrowmaps();
	  //@}
	
private:
	
	// don't want = operator, cctor and destructor
	
	VolConstrManager operator = (const VolConstrManager& old);
	VolConstrManager(const DRT::VolConstrManager& old);    

	
	//! name@ Private Attributes
	DRT::Discretization& actdisc_;	
	RefCountPtr<Epetra_SerialDenseVector> referencevolumes_;
	RefCountPtr<Epetra_SerialDenseVector> initialvolumes_;
	RefCountPtr<Epetra_SerialDenseVector> actvol_;
	RefCountPtr<Epetra_SerialDenseVector> volerr_;
	int minConstrID_;
	int maxConstrID_;
	int numConstrID_;
	double fact_; 
	RefCountPtr<Epetra_SerialDenseVector> lagrMultVec_;
	RefCountPtr<Epetra_SerialDenseVector> lagrMultInc_;
	map<int, RCP<Epetra_Map> > voldofrowmaps_;
	RefCountPtr<Epetra_Vector> constrVec_;
	//@}
    

}; //class
}//namespace DRT

#endif /*VOLUMECONSTRAINTMANAGER_H*/
