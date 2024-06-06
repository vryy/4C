/*----------------------------------------------------------------------*/
/*! \file
\brief boundary calc base routines
\level 2

*/
/*--------------------------------------------------------------------------*/


#include "4C_elemag_ele_boundary_calc.hpp"

#include "4C_discretization_fem_general_node.hpp"
#include "4C_elemag_ele.hpp"
#include "4C_elemag_ele_action.hpp"
#include "4C_global_data.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Discret::ELEMENTS::ElemagBoundaryImplInterface*
Discret::ELEMENTS::ElemagBoundaryImplInterface::Impl(const Core::Elements::Element* ele)
{
  switch (ele->Shape())
  {
    case Core::FE::CellType::quad4:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::quad4>::Instance();
    }
    case Core::FE::CellType::quad8:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::quad8>::Instance();
    }
    case Core::FE::CellType::quad9:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::quad9>::Instance();
    }
    case Core::FE::CellType::tri3:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::tri3>::Instance();
    }
    case Core::FE::CellType::tri6:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::tri6>::Instance();
    }
    case Core::FE::CellType::line2:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::line2>::Instance();
    }
    case Core::FE::CellType::line3:
    {
      return ElemagBoundaryImpl<Core::FE::CellType::line3>::Instance();
    }
    case Core::FE::CellType::nurbs2:  // 1D nurbs boundary element
    {
      return ElemagBoundaryImpl<Core::FE::CellType::nurbs2>::Instance();
    }
    case Core::FE::CellType::nurbs3:  // 1D nurbs boundary element
    {
      return ElemagBoundaryImpl<Core::FE::CellType::nurbs3>::Instance();
    }
    case Core::FE::CellType::nurbs4:  // 2D nurbs boundary element
    {
      return ElemagBoundaryImpl<Core::FE::CellType::nurbs4>::Instance();
    }
    case Core::FE::CellType::nurbs9:  // 2D nurbs boundary element
    {
      return ElemagBoundaryImpl<Core::FE::CellType::nurbs9>::Instance();
    }
    default:
      FOUR_C_THROW(
          "Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->num_node());
      break;
  }
}

template <Core::FE::CellType distype>
Discret::ELEMENTS::ElemagBoundaryImpl<distype>*
Discret::ELEMENTS::ElemagBoundaryImpl<distype>::Instance(Core::UTILS::SingletonAction action)
{
  static auto singleton_owner = Core::UTILS::MakeSingletonOwner(
      []()
      {
        return std::unique_ptr<Discret::ELEMENTS::ElemagBoundaryImpl<distype>>(
            new Discret::ELEMENTS::ElemagBoundaryImpl<distype>());
      });

  return singleton_owner.Instance(action);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
Discret::ELEMENTS::ElemagBoundaryImpl<distype>::ElemagBoundaryImpl()
    : xyze_(true),
      funct_(true),
      deriv_(true),
      unitnormal_(true),
      velint_(true),
      drs_(0.0),
      fac_(0.0)
{
  return;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ElemagBoundaryImpl<distype>::evaluate_neumann(
    Discret::ELEMENTS::ElemagBoundary* ele, Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseMatrix* elemat1_epetra)
{
  return 0;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype>
int Discret::ELEMENTS::ElemagBoundaryImpl<distype>::Evaluate(Discret::ELEMENTS::ElemagBoundary* ele,
    Teuchos::ParameterList& params, Discret::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  /* the term representing absorbing first order boundary conditions for the
   * here given problem looks like < lambda, mu > over Gamma_ext, hence it belongs
   * to the matrix Gmat evaluated at Gamma_ext. When condensing the local
   * unknowns we build K with G as summand. Hence, we can just add the terms
   * resulting from this boundary condition to K (and hence G)
   */
  const EleMag::Action action = Core::UTILS::GetAsEnum<EleMag::Action>(params, "action");
  switch (action)
  {
    case EleMag::calc_abc:
    {
      const int* nodeids = ele->NodeIds();

      Core::Elements::Element* parent = ele->parent_element();
      Teuchos::RCP<Core::Elements::FaceElement>* faces = parent->Faces();
      bool same = false;
      for (int i = 0; i < parent->NumFace(); ++i)
      {
        const int* nodeidsfaces = faces[i]->NodeIds();

        if (faces[i]->num_node() != ele->num_node()) break;

        for (int j = 0; j < ele->num_node(); ++j)
        {
          if (nodeidsfaces[j] == nodeids[j])
            same = true;
          else
          {
            same = false;
            break;
          }
        }
        if (same == true)
        {
          // i is the number we were searching for!!!!
          params.set<int>("face", i);
          ele->parent_element()->Evaluate(params, discretization, lm, elemat1_epetra,
              elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
          // break;
        }
      }
      if (same == false && (faces[0]->num_node() != ele->num_node()))
      {
        // in this case we have a three dimensional problem and the absorbing boundary condition on
        // a line and not on a surface. hence, we have to evaluate the abc term only at a part of
        // the face. here, we want to figure, which part!
        int elenode = ele->num_node();
        int face = -1;
        if (elenode != 2)
          FOUR_C_THROW(
              "absorbing line in 3d not implemented for higher order geometry approximation");
        // find the first face which contains the line!
        for (int i = 0; i < parent->NumFace(); ++i)
        {
          const int* nodeidsfaces = faces[i]->NodeIds();

          int count = 0;
          for (int j = 0; j < faces[i]->num_node(); ++j)
          {
            for (int n = 0; n < elenode; ++n)
            {
              count += (nodeidsfaces[j] == nodeids[n]);
            }
          }
          if (count == elenode)
          {
            same = true;
            face = i;
            params.set<int>("face", i);

            const int* nodeidsface = faces[face]->NodeIds();
            Teuchos::RCP<std::vector<int>> indices = Teuchos::rcp(new std::vector<int>(elenode));
            for (int j = 0; j < faces[face]->num_node(); ++j)
            {
              for (int n = 0; n < elenode; ++n)
                if (nodeids[n] == nodeidsface[j]) (*indices)[n] = j;
            }
            params.set<Teuchos::RCP<std::vector<int>>>("nodeindices", indices);
            ele->parent_element()->Evaluate(params, discretization, lm, elemat1_epetra,
                elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
          }
        }
        if (same == false) FOUR_C_THROW("no face contains absorbing line");
        // now, we know which face contains the line, but we have to tell the element, which nodes
        // we are talking about! therefore, we create a vector of ints and this vector stores the
        // relevant nodes, for example: the line has two nodes, we store which position these nodes
        // have in the face element
      }
      // if (same == false)
      //  FOUR_C_THROW("either nodeids are sorted differently or a boundary element does not know to
      //  whom it belongs");
      break;
    }
    case EleMag::bd_integrate:
    {
      const int* nodeids = ele->NodeIds();

      Core::Elements::Element* parent = ele->parent_element();
      Teuchos::RCP<Core::Elements::FaceElement>* faces = parent->Faces();
      bool same = false;
      for (int i = 0; i < parent->NumFace(); ++i)
      {
        const int* nodeidsfaces = faces[i]->NodeIds();

        if (faces[i]->num_node() != ele->num_node()) break;
        /*
        for(int j=0; j<ele->num_node(); ++j)
        {
          if(nodeidsfaces[j]==nodeids[j])
            same = true;
          else
          {
            same = false;
            break;
          }
        }
        */
        int count = 0;
        for (int j = 0; j < ele->num_node(); ++j)
        {
          for (int k = 0; k < ele->num_node(); ++k)
          {
            if (nodeidsfaces[j] == nodeids[k]) count++;
          }
        }
        if (count == ele->num_node()) same = true;


        if (same == true)
        {
          // i is the number we were searching for!!!!
          params.set<int>("face", i);
          ele->parent_element()->Evaluate(params, discretization, lm, elemat1_epetra,
              elemat2_epetra, elevec1_epetra, elevec2_epetra, elevec3_epetra);
          break;
        }
      }

      break;
    }
    default:
      FOUR_C_THROW("unknown action %d provided to ElemagBoundaryImpl", action);
      break;
  }
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
