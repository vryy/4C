/*!----------------------------------------------------------------------
\file drt_potential_volume.H

\brief  Class controlling volume stresses due to 
        intermolecular interaction forces
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



#include "../drt_potential/drt_potential.H"

namespace POTENTIAL
{
/*!
\brief  Class controlling volume stresses due to 
        intermolecular interaction forces
        between mesoscopic structures

*/
class VolumePotential : public Potential
{

public:

  /*!
    \brief Standard Constructor

    \param discret (in): Discretization
    \param discret (in): Discretization
    \param discret (in): treetype
   */
  VolumePotential(
      Teuchos::RCP<DRT::Discretization>     discretRCP,
      DRT::Discretization&                  discret,
      const GEO::TreeType&                  treetype
      );

  
  /*!
    \brief Destructor
   */
  ~VolumePotential(){};

  
  /*!
    \brief Call discretization to calculate additional contributions
    due to interfacial energy

    \param p (in): parameter list
    \param disp (in): displacement
    \param fint (out): internal forces
    \param stiff (out): stiffness matrix
   */
  void EvaluatePotential( 
		  ParameterList& p,
		  RCP<Epetra_Vector> 			disp,
		  RCP<Epetra_Vector> 			fint,
		  RCP<LINALG::SparseMatrix> 	stiff
		  );
  
  
  /*!
    \brief calculates internal forces and
    stiffnesses due to interfacial phenomena for volume elements
    in 3D

    \param element    (in)  : element
    \param gaussrule  (in)  : gaussrule of element
    \param params     (in)  : first partial derivatives of surface area
    \param K_surf     (out) : element stiffness matrix
    \param fint       (out) : element internal force vector
   */
  void StiffnessAndInternalForcesPotential(
      const DRT::Element*             element,
      const DRT::UTILS::GaussRule3D&  gaussrule,
      ParameterList&                  params,
      vector<int>&                    lm,
      Epetra_SerialDenseMatrix&       K_stiff,
      Epetra_SerialDenseVector&       F_int
      );
  
  
  /*!
    \brief calculates internal forces and
    stiffnesses due to interfacial phenomena for surface elements
    in 2D

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
      Epetra_SerialDenseMatrix&         K_stiff,
      Epetra_SerialDenseVector&         fint
      );
  
  
  // TEST Van Der Waals potential
  void TestEvaluatePotential(
    ParameterList&                      p,
    RCP<Epetra_Vector>          disp,
    RCP<Epetra_Vector>          fint,
    RCP<LINALG::SparseMatrix>   stiff,
    const double                        time,
    const int                           step
    );


  
private:

  //don't want = operator and cctor
  VolumePotential operator = (const POTENTIAL::VolumePotential& old);
  VolumePotential(const POTENTIAL::VolumePotential& old);


  /*!
  \brief Update displacements
   */
  void UpdateDisplacementsOfPotentialDiscretization(
    Teuchos::RCP<const Epetra_Vector>   idisp_solid
    );
  
  
  void EvaluateVolumePotentialCondition(
      ParameterList&                          params,
      RCP<LINALG::SparseMatrix>       systemmatrix1,
      RCP<LINALG::SparseMatrix>       systemmatrix2,
      RCP<Epetra_Vector>              systemvector1,
      Teuchos::RCP<Epetra_Vector>             systemvector2,
      Teuchos::RCP<Epetra_Vector>             systemvector3,
      const string&                           condstring
      );
  

  /*!
    \brief parallel volume tree search
     */
  void TreeSearch(
      const DRT::Element*      element,
      const double             radius,
      const int                label
      );
  
  
  /*!
    \brief parallel volume dummy tree search
     */
  void TreeSearchDummy();
  
  
  /*!
  \brief volume tree search
   */
  void TreeSearchElement(
      const LINALG::Matrix<3,2>&                              eleXAABB,
      const int                                               label,
      const bool                                              dummy
      );

  
  /*!
    \brief  compute internal force and stiffness matrix
            for local volume elements and nonlocal volume pecs
   */
  void ComputeFandK(
      const DRT::Element*                                             actEle,
      const DRT::UTILS::GaussRule3D&                                  gaussrule,
      vector<int>&                                                    lm,
      Epetra_SerialDenseMatrix&                                       K_stiff,
      Epetra_SerialDenseVector&                                       F_int,
      RCP<DRT::Condition>                                     cond,
      const double                                                    curvefac);


  /*!
      \brief  compute internal force and stiffness matrix
              for local surface elements and nonlocal surface pecs
   */
  void ComputeFandK(
      const DRT::Element*                                         actEle,
      const DRT::UTILS::GaussRule2D&                              gaussrule,
      vector<int>&                                                lm,
      Epetra_SerialDenseMatrix&                                   K_stiff,
      Epetra_SerialDenseVector&                                   F_int,
      RCP<DRT::Condition>                                 cond,
      const double                                                curvefac
  );

  
  //! compute integration factors for volume elements
  double ComputeFactor(
      const DRT::Element*                     element,
      Epetra_SerialDenseVector&               funct,
      Epetra_SerialDenseMatrix&               deriv,
      const DRT::UTILS::IntegrationPoints3D&  intpoints,
      const int                               gp,
      LINALG::Matrix<3,1>&                    x_gp,
      const double                            curve_fac
      );
  
  
  //! compute integration factors for volume pecs
  double ComputeFactor(
      PotentialElementContainer&              pec,
      Epetra_SerialDenseVector&               funct,
      Epetra_SerialDenseMatrix&               deriv,
      const DRT::UTILS::IntegrationPoints3D&  intpoints,
      const int                               gp,
      LINALG::Matrix<3,1>&                    x_gp,
      const double                            curve_fac
      );


  //! compute integration factors for surface elements 
  double ComputeFactor(
      const DRT::Element*                     element,
      Epetra_SerialDenseVector&               funct,
      Epetra_SerialDenseMatrix&               deriv,
      const DRT::UTILS::IntegrationPoints2D&  intpoints,
      const int                               gp,
      LINALG::Matrix<2,1>&                    x_gp,
      const double                            curve_fac
      );
  
  
  //! compute integration factors for surface pecs
  double ComputeFactor(
      PotentialElementContainer&              pec,
      Epetra_SerialDenseVector&               funct,
      Epetra_SerialDenseMatrix&               deriv,
      const DRT::UTILS::IntegrationPoints2D&  intpoints,
      const int                               gp,
      LINALG::Matrix<2,1>&                    x_gp,
      const double                            curve_fac
      );


  //! returns gauss rule for given distype
  void GetGaussRule3D(
      const DRT::Element::DiscretizationType& distype,
      DRT::UTILS::GaussRule3D&                rule_pot
      );
  
  
  //! returns gauss rule for given distype
  void GetGaussRule2D(
      const DRT::Element::DiscretizationType& distype,
      DRT::UTILS::GaussRule2D&                rule_pot
      );
 


  //! name@ Private Attributes
  //! search tree type
  GEO::TreeType                                           treetype_;
  
  //! for each potential condition label there is a set of elements
  std::map<int,std::set<int> >                            elementsByLabel_;
  
  //! for each element there is a potential condition 
  std::map<int,int >                                      labelByElement_;
  
  //! column vector of displacements 
  Teuchos::RCP<const Epetra_Vector>                       disp_col_;

  //! current position of boundary elements
  std::map<int,LINALG::Matrix<3,1> >                      currentpositions_;

  //! list of current XAABBs 
  std::map<int,LINALG::Matrix<3,2> >                      elemXAABBList_;
  
  //! local potential element ids
  std::map<int,std::set<int> >                            localEleIds_;
  
  //! nonlocal potential element containers
  std::map<int,std::vector<PotentialElementContainer> >   nonlocalPecs_;
             

}; //class

}//namespace UTILS

#endif /*DRT_POTENTIAL_VOLUME_H_*/
