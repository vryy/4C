/*----------------------------------------------------------------------*/
/*! \file
\brief Declaration

\level 2
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/
#ifndef LINALG_FIXEDSIZEBLOCKMATRIX_H
#define LINALG_FIXEDSIZEBLOCKMATRIX_H

#include "linalg_fixedsizematrix.H"

namespace LINALG
{
  template <class value_type, unsigned int brows, unsigned int bcols>
  class BlockMatrix
  {
   public:
    typedef typename value_type::scalar_type scalar_type;

    BlockMatrix() { std::fill(blocks_, blocks_ + brows * bcols, static_cast<value_type*>(NULL)); }

    ~BlockMatrix()
    {
      for (unsigned i = 0; i < brows * bcols; ++i)
      {
        delete blocks_[i];
      }
    }

    bool IsUsed(unsigned row, unsigned col) const { return blocks_[Position(row, col)] != NULL; }

    bool IsUsed(unsigned pos) const { return blocks_[pos] != NULL; }

    void Clear(unsigned row, unsigned col)
    {
      int p = Position(row, col);
      delete blocks_[p];
      blocks_[p] = NULL;
    }

    void AddView(unsigned row, unsigned col, value_type& matrix)
    {
      Clear(row, col);
      int p = Position(row, col);
      blocks_[p] = new value_type(matrix, true);
    }

    value_type* operator()(unsigned row, unsigned col)
    {
      int p = Position(row, col);
      value_type* b = blocks_[p];
      if (b == NULL)
      {
        b = new value_type(true);
        blocks_[p] = b;
      }
      return b;
    }

    const value_type* operator()(unsigned row, unsigned col) const
    {
      const value_type* b = blocks_[Position(row, col)];
#ifdef DEBUG
      if (b == NULL) dserror("null block access");
#endif
      return b;
    }

    value_type* operator()(unsigned pos)
    {
      value_type* b = blocks_[pos];
      if (b == NULL)
      {
        b = new value_type(true);
        blocks_[pos] = b;
      }
      return b;
    }

    const value_type* operator()(unsigned pos) const
    {
#ifdef DEBUG
      if (pos >= brows * bcols) dserror("block index out of range");
#endif
      const value_type* b = blocks_[pos];
#ifdef DEBUG
      if (b == NULL) dserror("null block access");
#endif
      return b;
    }

    template <class lhs, class rhs, unsigned int inner>
    inline void Multiply(
        const BlockMatrix<lhs, brows, inner>& left, const BlockMatrix<rhs, inner, bcols>& right)
    {
      PutScalar(0.);
      for (unsigned int ic = 0; ic < bcols; ++ic)
      {
        for (unsigned int i = 0; i < inner; ++i)
        {
          if (right.IsUsed(i, ic))
          {
            for (unsigned int ir = 0; ir < brows; ++ir)
            {
              if (left.IsUsed(ir, i))
              {
                value_type* b = (*this)(ir, ic);
                b->Multiply(1, *left(ir, i), *right(i, ic), 1);
              }
            }
          }
        }
      }
    }

    void PutScalar(scalar_type scalar)
    {
      for (unsigned i = 0; i < brows * bcols; ++i)
      {
        if (blocks_[i] != NULL)
        {
          blocks_[i]->PutScalar(scalar);
        }
      }
    }

    template <int rows, int cols>
    void AssembleTo(LINALG::TMatrix<scalar_type, rows, cols>& dest, double f)
    {
      const int local_rows = rows / brows;
      const int local_cols = cols / bcols;
      const unsigned row_blocks = brows;
      const unsigned col_blocks = bcols;
      for (unsigned icb = 0; icb < col_blocks; ++icb)
      {
        for (unsigned irb = 0; irb < row_blocks; ++irb)
        {
          if (IsUsed(irb, icb))
          {
            LINALG::Matrix<local_rows, local_cols>& local = *(*this)(irb, icb);
            for (int ic = 0; ic < local_cols; ++ic)
            {
              unsigned c = col_blocks * ic + icb;
              for (int ir = 0; ir < local_rows; ++ir)
              {
                unsigned r = row_blocks * ir + irb;
                dest(r, c) += f * local(ir, ic);
              }
            }
          }
        }
      }
    }

    friend std::ostream& operator<<(
        std::ostream& stream, const BlockMatrix<value_type, brows, bcols>& matrix)
    {
      stream << "BlockMatrix<" << brows << "," << bcols << ">\n";
      for (unsigned int ir = 0; ir < brows; ++ir)
      {
        for (unsigned int ic = 0; ic < bcols; ++ic)
        {
          int p = matrix.Position(ir, ic);
          value_type* b = matrix.blocks_[p];
          stream << "[" << ir << "," << ic << "] = ";
          if (b == NULL)
          {
            stream << "NULL\n";
          }
          else
          {
            stream << (*b);
          }
        }
      }
      return stream;
    }

   private:
    int Position(unsigned row, unsigned col) const
    {
#ifdef DEBUG
      if (row >= brows or col >= bcols) dserror("block index out of range");
#endif
      return row + brows * col;
    }

    value_type* blocks_[brows * bcols];
  };

}  // namespace LINALG

#endif
