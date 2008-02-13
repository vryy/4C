

#ifdef CCADISCRET
#ifndef INV_ANALYSIS_H
#define INV_ANALYSIS_H

#include "Epetra_Vector.h"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_condition.H"
#include "../drt_lib/drt_element.H"
#include "../drt_structure/strugenalpha.H"
#include <math.h>

class Inv_analysis : public StruGenAlpha {
	
	private:
		vector<int> surfneum_nodes_;
		vector<double> final_disp_;
		vector<double> measured_disp_;
		vector<double> residual_disp_;
		vector<double> residual_disp_o_;                
		vector<double> p_ ;    // parameter vector
		vector<double> p_o_;   // old parameter vector
		double aim_function_;
		RefCountPtr<Epetra_CrsMatrix> jacobian_;
		vector<double> parameters_;
		RefCountPtr<Epetra_Vector>  ref_disp_;
		const MATERIAL& mat_;
		MATERIAL& nonconstmat_;
		int count_;
		int numb_run_;
		
	public:
			
	  Inv_analysis(ParameterList& params,
	                          DRT::Discretization& dis,
	                          LINALG::Solver& solver,
	                          IO::DiscretizationWriter& output,
	                          bool init = true);
	  
		void get_surfneum_nodes(vector<DRT::Condition*> surfneum);
		void get_final_displacment(const RefCountPtr<Epetra_Vector> disp);
		void get_measured_displacment (const double measured_displacment_value);
		void get_old_parameters();
		bool compare_displacment ();
		void calculate_new_parameters();
		void evaluate();
		void Integrate();
		
		// virtual void Integrate();

};
		




#endif /*INV_ANALYSIS_H_*/
#endif  // #ifdef CCADISCRET
