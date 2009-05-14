/*!----------------------------------------------------------------------
\file drt_potential_volume.H

\brief  Class controlling volume stresses due to potential forces
        between mesoscopic structures

<pre>
Maintainer: Ursula Mayer
            mayer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>

*----------------------------------------------------------------------*/

#ifndef DRT_POTENTIAL_VOLUME_H_
#define DRT_POTENTIAL_VOLUME_H_


#include "../drt_surfstress/drt_potential.H"


using namespace std;
using namespace Teuchos;


namespace POTENTIAL
{
/*!
\brief  Class controlling volume stresses due to potential forces

*/
class VolumePotential : public Potential
{
  
public:

  /*!
    \brief Standard Constructor

    \param discret (in): Discretization
   */
  VolumePotential(
      Teuchos::RCP<DRT::Discretization>     discretRCP,
      DRT::Discretization&                  discret
      );

  /*!
    \brief Destructor
   */
  ~VolumePotential(){};

  /*!
    \brief Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
   */
  void EvaluatePotential( ParameterList& p,
      RefCountPtr<Epetra_Vector> disp,
      RefCountPtr<Epetra_Vector> fint,
      RefCountPtr<LINALG::SparseMatrix> stiff);

  
  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena for surface elements

    \param element    (in)  : element
    \param gaussrule  (in)  : gaussrule of element
    \param params     (in)  : first partial derivatives of surface area
    \param K_surf     (out) : element stiffness matrix
    \param fint       (out) : element internal force vector
   */
  void StiffnessAndInternalForcesPotential(
      const DRT::Element*               element,
      const DRT::UTILS::GaussRule2D&    gaussrule,
      ParameterList&                    params,
      vector<int>&                      lm,
      Epetra_SerialDenseMatrix&         K_surf,
      Epetra_SerialDenseVector&         fint);
  
  
  /*!
    \brief Routine to calculate additional internal forces and
    stiffnesses due to interfacial phenomena for line elements

    \param element    (in)  : element
    \param gaussrule  (in)  : gaussrule of element
    \param params     (in)  : first partial derivatives of surface area
    \param K_surf     (out) : element stiffness matrix
    \param fint       (out) : element internal force vector
   */
  void StiffnessAndInternalForcesPotential(
      const DRT::Element*               element,
      const DRT::UTILS::GaussRule1D&    gaussrule,
      ParameterList&                    params,
      vector<int>&                      lm,
      Epetra_SerialDenseMatrix&         K_surf,
      Epetra_SerialDenseVector&         fint);


private:

  //don't want = operator and cctor
  VolumePotential operator = (const POTENTIAL::VolumePotential& old);
  VolumePotential(const POTENTIAL::VolumePotential& old);


  /*!
  \brief Update displacements
   */
  void UpdateDisplacementsOfPotentialDiscretization(
    Teuchos::RCP<Epetra_Vector>   idisp_solid
    );


  /*!
  \brief  compute internal force and stiffness matrix
          for surface elements
   */
  void computeFandK(
      const DRT::Element*              actEle,
      const DRT::UTILS::GaussRule2D&   gaussrule,
      std::map<int,std::set<int> >&    potElements,
      vector<int>&                     lm,
      Epetra_SerialDenseMatrix&        K_surf,
      Epetra_SerialDenseVector&        F_int,
      RefCountPtr<DRT::Condition>      cond,
      const int                        label,
      const double                     curvefac
      );
  
  
  /*!
  \brief  compute internal force and stiffness matrix
          for line elements
   */
  void computeFandK(
      const DRT::Element*              actEle,
      const DRT::UTILS::GaussRule1D&   gaussrule,
      std::map<int,std::set<int> >&    potElements,
      vector<int>&                     lm,
      Epetra_SerialDenseMatrix&        K_surf,
      Epetra_SerialDenseVector&        F_int,
      RefCountPtr<DRT::Condition>      cond,
      const int                        label,
      const double                     curvefac
      );
  
  
  
  //! compute integration factors for surface elements
  double ComputeFactor(
      const DRT::Element*                     element, 
      LINALG::SerialDenseVector&              funct, 
      LINALG::SerialDenseMatrix&              deriv, 
      const DRT::UTILS::IntegrationPoints2D&  intpoints,
      const int                               gp,
      LINALG::Matrix<3,1>&                    x_gp,
      const double                            curve_fac
      );
 
  
  //! compute integration factors for line elements
  double ComputeFactor(
      const DRT::Element*                     element, 
      LINALG::SerialDenseVector&              funct, 
      LINALG::SerialDenseMatrix&              deriv, 
      const DRT::UTILS::IntegrationPoints1D&  intpoints,
      const int                               gp,
      LINALG::Matrix<3,1>&                    x_gp,
      const double                            curve_fac
      );
  
  
  //! returns gauss rule for given distype
  void GetGaussRule2D(
      const DRT::Element::DiscretizationType& distype, 
      DRT::UTILS::GaussRule2D&                rule_pot
      );
  
  
  //! compute surface normal in gaussian point
  LINALG::Matrix<3,1> ComputeNormalInGP(
      const DRT::Element*                     element, 
      const LINALG::Matrix<3,1>&              x_gp
      );
  
  
  //! determine valid contribution
  bool DetermineValidContribution(
      const LINALG::Matrix<3,1>&              x_gp,
      const LINALG::Matrix<3,1>&              x_pot_gp,
      const LINALG::Matrix<3,1>&              normal_gp,
      const LINALG::Matrix<3,1>&              normal_pot_gp
      );
 
  
  //! name@ Private Attributes

  //! current position of boundary elements
  std::map<int,LINALG::Matrix<3,1> >          currentpositions_;

  //! for each potential condition label there is a set of boundary elements
 std::map<int,std::set<int> >                 elementsByLabel_;

  //! potential boundary
  LINALG::MapExtractor                        potboundary_;

  //! @name local copies of input parameters
  Teuchos::RCP<DRT::Discretization>           potentialdis_;
 
  //! @name displacement vector of redundant boundary discretization
  Teuchos::RCP<Epetra_Vector>                 idisp_total_;
  
  //! @name displacement vector of boundary discretization, which is part of the local solid dis
  Teuchos::RCP<Epetra_Vector>                 idisp_onproc_;
  
  //! @name importer
  Teuchos::RCP<Epetra_Import>                 importer_;


}; //class

}//namespace UTILS

#endif /*DRT_POTENTIAL_VOLUME_H_*/
