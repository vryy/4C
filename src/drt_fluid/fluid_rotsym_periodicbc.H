/*-----------------------------------------------------------*/
/*!

\brief Methods to apply rotationally symmetric periodic boundary
       conditions for fluid problems on element level

\maintainer Martin Kronbichler

\level 1

*/
/*-----------------------------------------------------------*/

#ifndef FLUID_ROTSYM_PERIODICBC_H
#define FLUID_ROTSYM_PERIODICBC_H

#include <vector>


#include "fluid_rotsym_periodicbc_utils.H"

namespace FLD
{
  /// numdofpernode = number of dofs per node for fluid problem
  template <DRT::Element::DiscretizationType distype, int numdofpernode,
      DRT::ELEMENTS::Fluid::EnrichmentType enrtype>

  /*!
  \brief  This class manages local transformations(rotation) of velocity fields

          o used for rotationally symmetric boundary conditions
          o used for rotationally symmetric boundary conditions

  \author gjb

   */
  class RotationallySymmetricPeriodicBC
  {
   public:
    /// number of nodes for this element type including virtual nodes
    static const int elenumnode = DRT::ELEMENTS::MultipleNumNode<enrtype>::multipleNode *
                                  DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    /// number of nodes for this element type (only real nodes)
    static const int elenumnodereal =
        DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

    /// standard constructor
    explicit RotationallySymmetricPeriodicBC() : rotangle_(0.0), rotmat_(true) { return; };

    /// empty destructor
    virtual ~RotationallySymmetricPeriodicBC(){};


    /// is any rotation needed on the element level?
    bool HasRotSymmPBC()
    {
      if (slavenodelids_.empty())
      {
        return (false);
      }
      else
      {
        return (true);
      }
    }

    /// prepare the class for this element
    void Setup(DRT::Element* ele)
    {
      // clean everything
      rotangle_ = 0.0;
      slavenodelids_.clear();
      rotmat_.Clear();

      DRT::Node** nodes = ele->Nodes();
      slavenodelids_.reserve(elenumnodereal);

      for (int inode = 0; inode < elenumnodereal; inode++)
      {
        if (FLD::IsSlaveNodeOfRotSymPBC(nodes[inode], rotangle_))
        {
          // store the lid of this node
          slavenodelids_.push_back(inode);
        }
      }

      // yes, we prepare the rotation matrix (identity matrix when there are no slave node)
      for (int i = 0; i < elenumnode * numdofpernode; i++)
      {
        rotmat_(i, i) = 1.0;
      }
      const double c = cos(rotangle_);
      const double s = sin(rotangle_);
      for (unsigned int slnode = 0; slnode < slavenodelids_.size(); slnode++)
      {
        // get local id of current slave node
        int inode = slavenodelids_[slnode];
        // velocity x- and y-components have to be rotated
        rotmat_(inode * numdofpernode, inode * numdofpernode) = c;
        rotmat_(inode * numdofpernode, inode * numdofpernode + 1) = s;
        rotmat_(inode * numdofpernode + 1, inode * numdofpernode) = -s;
        rotmat_(inode * numdofpernode + 1, inode * numdofpernode + 1) = c;
        // corresponding z- component and pressure entries remain unchanged!
      }

      return;
    }

    /// rotate velocity vector used in element routine if necessary
    void RotateMyValuesIfNecessary(std::vector<double>& myvalues)
    {
      if (HasRotSymmPBC())
      {
        // rotate velocity vectors to right position (use transposed rotation matrix!!)
        const double c = cos(rotangle_);
        const double s = sin(rotangle_);
        for (unsigned int slnodes = 0; slnodes < slavenodelids_.size(); slnodes++)
        {
          int i = slavenodelids_[slnodes];
          double xvalue = myvalues[0 + (i * numdofpernode)];
          double yvalue = myvalues[1 + (i * numdofpernode)];
          myvalues[0 + (i * numdofpernode)] = c * xvalue - s * yvalue;
          myvalues[1 + (i * numdofpernode)] = s * xvalue + c * yvalue;
        }
      }
      return;
    }

    /// rotate velocity vector used in element routine if necessary
    void RotateMyValuesIfNecessary(LINALG::Matrix<numdofpernode - 1, elenumnode>& myvalues)
    {
      if (HasRotSymmPBC())
      {
        // rotate velocity vectors to right position (use transposed rotation matrix!!)
        const double c = cos(rotangle_);
        const double s = sin(rotangle_);
        for (unsigned int slnodes = 0; slnodes < slavenodelids_.size(); slnodes++)
        {
          const int i = slavenodelids_[slnodes];
          const double xvalue = myvalues(0, i);
          const double yvalue = myvalues(1, i);
          myvalues(0, i) = c * xvalue - s * yvalue;
          myvalues(1, i) = s * xvalue + c * yvalue;
        }
      }
      return;
    }

    /// rotate velocity vector used in element routine if necessary
    template <int rows, int cols>
    void RotateMyValuesIfNecessary(LINALG::Matrix<rows, cols>& myvalues)
    {
      if (HasRotSymmPBC())
      {
        // rotate velocity vectors to right position (use transposed rotation matrix!!)
        const double c = cos(rotangle_);
        const double s = sin(rotangle_);
        for (unsigned int slnodes = 0; slnodes < slavenodelids_.size(); slnodes++)
        {
          const int i = slavenodelids_[slnodes];
          const double xvalue = myvalues(0, i);
          const double yvalue = myvalues(1, i);
          myvalues(0, i) = c * xvalue - s * yvalue;
          myvalues(1, i) = s * xvalue + c * yvalue;
        }
      }
      return;
    }

    /// rotate element matrix and vectors if necessary (first version)
    void RotateMatandVecIfNecessary(
        LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode>& elemat1,
        LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode>& elemat2,
        LINALG::Matrix<numdofpernode * elenumnode, 1>& elevec1)
    {
      if (HasRotSymmPBC())
      {
        // execute the rotation
        /*
        K_rot = Q* K * Q^T
        b_rot = Q * b
        with Q as defined in the Setup() function
         */
        if (elemat1.IsInitialized())  // do not try to access an uninitialized matrix!
        {
          LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> elematold1(
              elemat1);
          LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> tempmatrix(true);
          tempmatrix.MultiplyNT(elematold1, rotmat_);
          elemat1.Multiply(rotmat_, tempmatrix);
        }
        if (elemat2.IsInitialized())
        {
          LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> elematold2(
              elemat2);
          LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> tempmatrix(true);
          tempmatrix.MultiplyNT(elematold2, rotmat_);
          elemat2.Multiply(rotmat_, tempmatrix);
        }
        if (elevec1.IsInitialized())
        {
          LINALG::Matrix<numdofpernode * elenumnode, 1> elevec1old(elevec1);
          elevec1.Multiply(rotmat_, elevec1old);
        }
      }
      return;
    }

    /// rotate element matrix and vectors if necessary (second version)
    void RotateMatandVecIfNecessary(
        LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode>& elemat1,
        LINALG::Matrix<numdofpernode * elenumnode, 1>& elevec1)
    {
      if (HasRotSymmPBC())
      {
        // execute the rotation
        /*
        K_rot = Q* K * Q^T
        b_rot = Q * b
        with Q as defined in the Setup() function
         */
        if (elemat1.IsInitialized())  // do not try to access an uninitialized matrix!
        {
          LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> elemat1old(
              elemat1);
          LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> tempmatrix(true);
          tempmatrix.MultiplyNT(elemat1old, rotmat_);
          elemat1.Multiply(rotmat_, tempmatrix);
        }

        if (elevec1.IsInitialized())
        {
          LINALG::Matrix<numdofpernode * elenumnode, 1> elevec1old(elevec1);
          elevec1.Multiply(rotmat_, elevec1old);
        }
      }
      return;
    }


   protected:
   private:
    //! vector of local slave node ids of the applied periodic surface boundary conditions
    std::vector<int> slavenodelids_;

    //! angle of rotation (slave plane <--> master plane)
    double rotangle_;

    //! rotation matrix
    LINALG::Matrix<numdofpernode * elenumnode, numdofpernode * elenumnode> rotmat_;
  };

}  // namespace FLD

#endif  // #ifndef FLUID_ROTSYM_PERIODICBC_H
