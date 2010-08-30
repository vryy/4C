#ifdef STKADAPTIVE

#ifndef STK_LINEARSYSTEM_H
#define STK_LINEARSYSTEM_H

#include <fei_Factory_Trilinos.hpp>
#include <fei_LinearSystem.hpp>

#include <Intrepid_FieldContainer.hpp>

#include <stk_linsys/DofMapper.hpp>
#include <stk_linsys/LinearSystem.hpp>
#include <stk_linsys/LinsysFunctions.hpp>
#include <stk_linsys/ImplDetails.hpp>

#include "stk_mesh.H"

namespace STK
{

class LinearSystem
{
public:

  LinearSystem( Mesh & mesh, stk::mesh::FieldBase * field );

  LinearSystem( Mesh & mesh, stk::mesh::FieldBase * field1, stk::mesh::FieldBase * field2 );

  void Commit();

  int Solve( int & status, Teuchos::ParameterList & params );

  void Done();

  void Assemble( stk::mesh::Bucket & bucket,
                 Intrepid::FieldContainer<double> & localStiffMatrix,
                 stk::mesh::Selector & sel,
                 stk::mesh::Part & active,
                 stk::mesh::Part & dbc );

  void Assemble( stk::mesh::Bucket & bucket,
                 Intrepid::FieldContainer<double> & localRHS,
                 stk::mesh::Selector & sel,
                 stk::mesh::Part & dbc );

private:

  void Setup();

  Mesh & m_mesh;
  std::vector<stk::mesh::FieldBase*> m_fields;
  fei::SharedPtr<fei::Factory> m_factory;
  stk::linsys::LinearSystem m_ls;
  stk::linsys::DofMapper * m_dof_mapper;

  fei::SharedPtr<fei::LinearSystem> m_fei_ls;
  fei::SharedPtr<fei::Matrix> m_fei_matrix;
  fei::SharedPtr<fei::Vector> m_fei_rhs;
  fei::SharedPtr<fei::Vector> m_fei_x;
};

}

#endif
#endif
