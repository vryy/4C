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
	//! @name Constructors and destructors and related methods

	  /*!
	  \brief Standard Constructor

	  */
	
	  VolConstrManager(const double time, DRT::Discretization& disc,
			  RCP<Epetra_Vector> disp,RCP<Epetra_Vector> resdisp);
	  
	  ~VolConstrManager(){return;};

	  //@}
	  
	  /*!
	    \brief Compute volume difference and change stiffness matrix and force vector accordingly
	    
	    Volume difference will be stored in volerr_.
	  */  
	  void StiffnessAndInternalForces(const double time, 
			  RCP<Epetra_Vector> disp,
			  RefCountPtr<Epetra_Vector> fint,
			  RefCountPtr<Epetra_CrsMatrix> stiff);
	  
	  /*
	   \brief Return norm of difference between actual and prescribed volumes
	  */
	  double GetVolumeErrorNorm()
	  {
		  return volerr_->Norm2();
	  };
	  	  
	  /*
	  	   \brief Return number of constrained volumes
	  */ 
	  int GetNumberOfVolumes()
	  {
		  return numConstrID_; 
	  };
	  
	  /*
	  	   \brief Return DofRowMap of volume i
	  */
	  Epetra_Map GetDofMap(int i)
	  {
		  return *(voldofrowmaps_[i]);
	  };
	  
	  /*
	  	   \brief Scale increment of lagrange multiplier increment by a double d
	  */
	  void ScaleLagrIncr(double d)
	  {
		  lagrMultInc_->Scale(d);
		  return;
	  };
	  
	  /*
	  	   \brief Update lagrange multiplier as needed in Uzawa algorithm
	  */
	  void UpdateLagrIncr(double factor, Epetra_SerialDenseVector vect);
	  
	  /*
	  	   \brief Return lagrange multiplier increment
	  */
	  double GetLagrIncr(int i)
	  {
		  return (*lagrMultInc_)[i];
	  };
	  
	  /*
	  	   \brief Return differences between prescribed and actual value of volume number i
	  */
	  double GetVolumeError(int i)
	  {
		  return (*volerr_)[i];
	  }
	  
	  /*
	  	   \brief Return additional vector containing constraints
	  */
	  RCP<Epetra_Vector> GetConstrVec()
	  {
		  return constrVec_;
	  };
	  
	  /*
	  	   \brief Add Increment to lagrange multiplier
	  */
	  void UpdateLagrMult();
	  
	  
private:
	
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

	
private:
	
	// don't want = operator, cctor and destructor
	
	VolConstrManager operator = (const VolConstrManager& old);
	VolConstrManager(const DRT::VolConstrManager& old);    

	
	//variables needed for volume constraint boundary condition
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
	
    

}; //class
}//namespace DRT

#endif /*VOLUMECONSTRAINTMANAGER_H*/
