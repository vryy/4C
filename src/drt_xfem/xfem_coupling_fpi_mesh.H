/*----------------------------------------------------------------------*/
/*!
\file xfem_coupling_fpi_mesh.H

\brief manages mesh based coupling of fluid and porous media
xfluid class and the cut-library

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef XFEM_COUPLING_FPI_MESH_H_
#define XFEM_COUPLING_FPI_MESH_H_

#include "xfem_coupling_mesh.H"

namespace XFEM
{

/*!
specialized class for coupling of fluid with an porous media (Darcy Flow)
  \author ager
  \date 06/16
 */
class MeshCouplingFPI : public MeshCoupling
{
public:

  enum coupled_field
  {
    ps_ps,
    ps_pf,
    pf_ps,
    pf_pf
  };

  //! constructor
  explicit MeshCouplingFPI(
      Teuchos::RCP<DRT::Discretization>&  bg_dis,   ///< background discretization
      const std::string &                 cond_name,///< name of the condition, by which the derived cutter discretization is identified
      Teuchos::RCP<DRT::Discretization>&  cond_dis,  ///< discretization from which cutter discretization can be derived
      const int                           coupling_id,///< id of composite of coupling conditions
      const double                        time,      ///< time
      const int                           step,       ///< time step
      MeshCouplingFPI::coupled_field      field      ///< which field is coupled to the fluid
  );

  //! destructor
  virtual ~MeshCouplingFPI()
  {
    return;
  }

  //! cutter dis should be loaded into the cut?
  virtual bool CutGeometry()
  {
    return (coupled_field_ == MeshCouplingFPI::ps_ps);
  }

  void SetFullDispnp(Teuchos::RCP<const Epetra_Vector> dispnp, Teuchos::RCP<const Epetra_Vector> pres)
  {
    fulldispnp_ = dispnp;
    fullpres_ = pres;
  }

  void InitializeStrucPresMap(Teuchos::RCP<const Epetra_Map> pfmap, Teuchos::RCP<const Epetra_Map> psmap )
  {
    if (pfmap->NumMyElements() != psmap->NumMyElements())
      dserror("InitializeStrucPresMap: (pfmap->NumGlobalElements() != psmap->NumGlobalElements())!");

    for (int lid = 0; lid < pfmap->NumMyElements(); ++lid)
      lm_struct_x_lm_pres_[psmap->GID(lid)] = pfmap->GID(lid)+1; //z-component of structure --> pressure
  }

  /*!
   Return prescribed velocities and traction vectors for a GNBC boundary condition.
   Also returns the projection matrix (to the plane of the surface) needed for the GNBC condition.
   */
  template <DRT::Element::DiscretizationType DISTYPE, class T1, class M3>
  void EvaluateCouplingConditions(
      T1&       projection_matrix,                ///< Projection matrix
      M3&       normal                         ///< surface normal of cut element
  )
{
  EvalProjectionMatrix<DISTYPE>(projection_matrix,normal);
  return;
};

  // finalize the interface true residual vector
  virtual void CompleteStateVectors();

  virtual void ZeroStateVectors_FPI();

  virtual void GmshOutput(
      const std::string & filename_base,
      const int step,
      const int gmsh_step_diff,
      const bool gmsh_debug_out_screen
  );

  virtual void GmshOutputDiscretization(
    std::ostream& gmshfilecontent
  );

  virtual void LiftDrag(
      const int step,
      const double time
  ) const;

  virtual void ReadRestart(
      const int step
  );

  // interface foces
  Teuchos::RCP<Epetra_Vector> ITrueResidual(){ return itrueresidual_; }

  // for assembly of fluid interface forces
  Teuchos::RCP<Epetra_Vector> IForcecol(){ return iforcecol_; }

private:

  void Output(
      const int step,
      const double time,
      const bool write_restart_data
  );

protected:

  //! Initialize
  virtual void Init();

  //! Initializes configurationmap
  virtual void InitConfigurationMap();

  //! set the name of the coupling object based on the field coupling
  virtual void SetCouplingName();

  //! Updates configurationmap for specific Gausspoint
  virtual void UpdateConfigurationMap_GP( double& kappa_m, double& visc_m, double& visc_s,
      double& visc_stab, double& full_stab, const LINALG::Matrix<3,1>& x, const DRT::Condition* cond);

  virtual void InitStateVectors();

  virtual bool HasMovingInterface(){ return true; }

  virtual void SetConditionSpecificParameters();

#if 0
  template <DRT::Element::DiscretizationType DISTYPE>
  double CalcPorosity(
      const Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DISTYPE> > & si,
      int coup_sid)
  {
    double porosity;
    double J = 0.0;
    double pres = 0.0;
    J = TDetDeformationGradient<DISTYPE>(si,coup_sid,pres);
  //  std::cout << "J: " << J << " pres: " << pres << std::endl;
    Teuchos::ParameterList params; //empty parameter list;

    DRT::Element* ele = cond_dis_->gElement(GetEmbeddedElementId(coup_sid));
    {
      //access second material in structure element
      if (ele->NumMaterial() > 1)
      {
        structmat_ = Teuchos::rcp_dynamic_cast<MAT::StructPoro>(ele->Material(0));
        if(structmat_->MaterialType() != INPAR::MAT::m_structporo and
           structmat_->MaterialType() != INPAR::MAT::m_structpororeaction and
           structmat_->MaterialType() != INPAR::MAT::m_structpororeactionECM)
          dserror("invalid structure material for poroelasticity");
      }
      else
        dserror("no second material defined for element %i",ele->Id());
    }

    structmat_->ComputePorosity(
        params,
        pres,
        J,
        1, //finally check what to do here Todo:)
        porosity,
        false);
    return porosity;
  }

  /*!
  \brief Templated Calculate Determinate of the Deformation Gradient at GP

  */
  template <DRT::Element::DiscretizationType DISTYPE>
  double TDetDeformationGradient(
      const Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementInterface<DISTYPE> > & si,
      int coup_sid,
      double& pres)
{
DRT::Element* ele = cond_dis_->gElement(GetEmbeddedElementId(coup_sid));
DRT::Element* fele = cutter_dis_->gElement(coup_sid);

if (fele->Shape() == DRT::Element::quad4)
{
  pres = 0.0;

  const DRT::Element::DiscretizationType SLAVE_DISTYPE = DRT::Element::quad4;
  const size_t PARENT_NEN = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;
  const unsigned int SLAVE_NUMDOF = 3;

  Teuchos::RCP<DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DISTYPE,SLAVE_DISTYPE,SLAVE_NUMDOF> > ser =
    Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::XFLUID::SlaveElementRepresentation<DISTYPE,SLAVE_DISTYPE,SLAVE_NUMDOF> >(si);
  if (ser == Teuchos::null) dserror("Cast of SER failed!");

  LINALG::Matrix<3,1> slave_xsi(true);
  ser->GetSlaveXsi(slave_xsi);


  DRT::UTILS::CollectedGaussPoints intpoints = DRT::UTILS::CollectedGaussPoints(1); //reserve just for 1 entry ...
  intpoints.Append(slave_xsi(0,0), slave_xsi(1,0),0.0, 1.0);

  // get coordinates of gauss point w.r.t. local parent coordinate system
  LINALG::SerialDenseMatrix pqxg(1,SLAVE_NUMDOF);
  LINALG::Matrix<SLAVE_NUMDOF,SLAVE_NUMDOF>  derivtrafo(true);

  DRT::UTILS::BoundaryGPToParentGP<SLAVE_NUMDOF>( pqxg,
      derivtrafo,
      intpoints ,
      DISTYPE,
      fele->Shape(),
      GetFaceLidOfEmbeddedElement(coup_sid));

  LINALG::Matrix<SLAVE_NUMDOF,1> pxsi        (true);

  // coordinates of the current integration point in parent coordinate system
    for (uint idim=0;idim<SLAVE_NUMDOF ;idim++)
  {
    pxsi(idim) = pqxg(0,idim);
  }

  LINALG::Matrix<PARENT_NEN,1> pfunc_loc(true); // derivatives of parent element shape functions in parent element coordinate system
  LINALG::Matrix<SLAVE_NUMDOF,PARENT_NEN> pderiv_loc(true); // derivatives of parent element shape functions in parent element coordinate system

  // evaluate derivatives of parent element shape functions at current integration point in parent coordinate system
  DRT::UTILS::shape_function<DISTYPE>(pxsi,pfunc_loc);
  DRT::UTILS::shape_function_deriv1<DISTYPE>(pxsi,pderiv_loc);
  //
  // get Jacobian matrix and determinant w.r.t. spatial configuration
  //
  // |J| = det(xjm) * det(Jmat^-1) = det(xjm) * 1/det(Jmat)
  //
  //    _                     _
  //   |  x_1,1  x_2,1  x_3,1  |           d x_i
  //   |  x_1,2  x_2,2  x_3,2  | = xjm  = --------
  //   |_ x_1,3  x_2,3  x_3,3 _|           d s_j
  //    _                     _
  //   |  X_1,1  X_2,1  X_3,1  |           d X_i
  //   |  X_1,2  X_2,2  X_3,2  | = Jmat = --------
  //   |_ X_1,3  X_2,3  X_3,3 _|           d s_j
  //
  LINALG::Matrix<SLAVE_NUMDOF,SLAVE_NUMDOF>    xjm;
  LINALG::Matrix<SLAVE_NUMDOF,SLAVE_NUMDOF>   Jmat;

  LINALG::Matrix<SLAVE_NUMDOF,PARENT_NEN>  xrefe (true);   // material coord. of parent element
  LINALG::Matrix<SLAVE_NUMDOF,PARENT_NEN>  xcurr (true);   // current  coord. of parent element

  // update element geometry of parent element
  {
    DRT::Node** nodes = ele->Nodes();
    for (uint inode=0;inode<PARENT_NEN;++inode)
    {
      for (unsigned int idof=0;idof<SLAVE_NUMDOF;++idof)
      {
        int lid = fulldispnp_->Map().LID(GetCondDis()->Dof(0,ele->Nodes()[inode],idof));

        const double* x = nodes[inode]->X();
        xrefe(idof,inode)   = x[idof];
        xcurr(idof,inode)   = xrefe(idof,inode) + fulldispnp_->operator [](lid);
      }
      int lidp = fullpres_->Map().LID(lm_struct_x_lm_pres_.operator [](GetCondDis()->Dof(0,ele->Nodes()[inode],2)));
      pres += fullpres_->operator [](lidp)*pfunc_loc(inode);
    }
  }
  xjm.MultiplyNT (pderiv_loc,xcurr);
  Jmat.MultiplyNT(pderiv_loc,xrefe);
  double det  = xjm.Determinant();
  double detJ = Jmat.Determinant();
//  std::cout << std::setprecision(30) << "det: " << det << " detJ: " << detJ << std::endl;
  const double J = det/detJ;

//  //Linearisation of Jacobian (atm missing!)
//  // D J[d] = J div(d) = J div(N_i d_i) = J dN_i/dx_j d_ij = J dN_i/dzeta_k dzeta_k/dx_j d_ij
//
////  LINALG::Matrix<dim,numnodes>  auxJLin (true);   // matrix to be filled with linearization scalars
//  xjm.Invert();
//  {
////    DRT::Node** nodes = sele.ParentElement()->Nodes();
//    for (int inode=0;inode<numnodes;++inode)
//    {
//      for (unsigned int i=0;i<dim;++i)
//      {
//        for (unsigned int j=0;j<dim;++j)
//        {
//          JLin[sele.MoData().ParentDof()[inode*dim+i]] += J*pderiv_loc(j,inode)*xjm(i,j);
////          auxJLin(idof,inode)+=J*pderiv_loc(j,inode)*xjm(i,j); //matrix that could be filled with linearizations
//        }
//      }
//    }
//  }
  return J;
}
else
  dserror("TDetDeformationGradient for type %s not yet implemented", (DRT::DistypeToString(fele->Shape())).c_str());
return 0.0;
}
#endif

  //! Access to map from cutting surface element id to embedded element id
  int GetEmbeddedElementId(int sid)
  {
    std::map<int,int>::iterator it = cutting_emb_gid_map_.find(sid);

    if (it == cutting_emb_gid_map_.end())
      dserror("Could not find parent element ID for side %d.", sid);

    return it->second;
  }

  //! Access to map from cutting surface element id to coupling element's face lid
  int GetFaceLidOfEmbeddedElement(int sid)
  {
    std::map<int,int>::iterator it = cutting_emb_face_lid_map_.find(sid);

    if (it == cutting_emb_face_lid_map_.end())
      dserror("No local id of the face-element % d w.r.t embedded element stored!", sid);

    return it->second;
  }

  //! create map from cutting surface element id to embedded element id
  void CreateCuttingToEmbeddedElementMap();

  //------------------------------- vectors -----------------------------
  //! @name cutter-dis state vectors
  Teuchos::RCP<Epetra_Vector>    itrueresidual_; //! interface forces acting on the structural surface (= -iforcenp)
  Teuchos::RCP<Epetra_Vector>    iforcecol_;     //! interface forces acting on the fluid surface (column vector assembly)
  //@}

  //---------------------------------configuration--------------------------------
  //! @name type of poro field coupled to xfluid handled by this mesh coupling object!
  MeshCouplingFPI::coupled_field coupled_field_;

  //! Full BJ Variant or BJSaffmann?
  bool full_BJ_;


  //! corresponding poro structure material
  //Teuchos::RCP< MAT::StructPoro > structmat_;

  Teuchos::RCP<const Epetra_Vector> fulldispnp_;
  Teuchos::RCP<const Epetra_Vector> fullpres_;

  //! map from cutting surface element id to embedded element id
  std::map<int,int> cutting_emb_gid_map_;

  //! map from cutting surface element id to embedded element's face lid
  std::map<int,int> cutting_emb_face_lid_map_;

  //! map from structural x dof to pres dof of a node!
  std::map<int,int> lm_struct_x_lm_pres_;

  double BJ_coeff_;

    //@}
};

} // namespace XFEM

#endif /*XFEM_COUPLING_FPI_MESH_H_*/
