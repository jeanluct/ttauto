// <LICENSE
//   ttauto: a C++ library for building train track automata
//
//   https://github.com/jeanluct/ttauto
//
//   Copyright (C) 2010-2026  Jean-Luc Thiffeault   <jeanluc@math.wisc.edu>
//                            Erwan Lanneau <erwan.lanneau@ujf-grenoble.fr>
//
//   This file is part of ttauto.
//
//   ttauto is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   ttauto is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with ttauto.  If not, see <http://www.gnu.org/licenses/>.
// LICENSE>

#ifndef TRAINTRACKS_MATHMATRIX_PERMPLUS1_HPP
#define TRAINTRACKS_MATHMATRIX_PERMPLUS1_HPP

#include <iostream>
#include <jlt/vector.hpp>
#include <jlt/mathmatrix.hpp>

// Sparse matrix type for matrices of the form permutation or permutation+1.

namespace traintracks {

class mathmatrix_permplus1;

// Left-multiplication by a permutation or permutation+1 matrix.
jlt::mathmatrix<int> operator*(const mathmatrix_permplus1& pm,
			       const jlt::mathmatrix<int>& A);

// Right-multiplication by a permutation or permutation+1 matrix.
jlt::mathmatrix<int> operator*(const jlt::mathmatrix<int>& A,
			       const mathmatrix_permplus1& pm);

std::ostream& operator<<(std::ostream& strm, const mathmatrix_permplus1& pm);

std::ostream& printMathematicaForm(std::ostream& strm,
				   const mathmatrix_permplus1& pm);


class mathmatrix_permplus1
{
public:
  typedef int			value_type;
  typedef int*			pointer;
  typedef const int*		const_pointer;
  typedef int&			reference;
  typedef const int&		const_reference;
  typedef size_t		size_type;

private:
  typedef jlt::vector<int>	Vec;
  typedef jlt::mathmatrix<int>	Mat;

  Vec rperm;			// Row permutation.
  Vec cperm;			// Inverse (column) permutation.
  int p1row, p1col;		// The 'plus 1' entry's coordinates.
  				// If both -1, then permutation matrix.

public:
  mathmatrix_permplus1(const Vec& rperm_ = Vec(),
			 const int p1row_ = -1, const int p1col_ = -1)
    : rperm(rperm_), cperm(rperm_.size(),-1), p1row(p1row_), p1col(p1col_)
  {
    for (int i = 0; i < (int)rperm.size(); ++i)
      {
	if (rperm[i] < 0 || rperm[i] >= (int)rperm.size() || cperm[rperm[i]] != -1)
	  {
	    std::cerr << "Bad row permutation in ";
	    std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
	    std::exit(1);
	  }
	cperm[rperm[i]] = i;
      }
  }

  // Construct from a mathmatrix.
  //
  // p1row and p1col are set by the entry where both the column sum
  // and row sum are 2.
  mathmatrix_permplus1(const Mat& M)
    : rperm(M.dim(),-1), cperm(M.dim(),-1), p1row(-1), p1col(-1)
  {
    const int n = M.dim();

    for (int i = 0; i < n; ++i)
      {
	int colsum = 0, rowsum = 0;
	for (int j = 0; j < n; ++j)
	  {
	    if (!(M(i,j) == 0 || M(i,j) == 1))
	      {
		std::cerr << "Matrix should contain only ones or zeros in ";
		std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
		std::exit(1);
	      }
	    colsum += M(i,j);
	    rowsum += M(j,i);
	  }

	if (colsum == 2)
	  {
	    if (p1row != -1)
	      {
		std::cerr << "More than one +1 row in ";
		std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
		std::exit(1);
	      }
	    p1row = i;
	  }
	else if (colsum != 1)
	  {
	    std::cerr << "Bad row sum in ";
	    std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
	    std::exit(1);
	  }

	if (rowsum == 2)
	  {
	    if (p1col != -1)
	      {
		std::cerr << "More than one +1 col in ";
		std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
		std::exit(1);
	      }
	    p1col = i;
	  }
	else if (rowsum != 1)
	  {
	    std::cerr << "Bad col sum in ";
	    std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
	    std::exit(1);
	  }
      }

    if ((p1row == -1) != (p1col == -1))
      {
	std::cerr << "Inconsistent +1 row/col in ";
	std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
	std::exit(1);
      }

    for (int i = 0; i < n; ++i)
      {
	for (int j = 0; j < n; ++j)
	  {
	    if (M(i,j) == 1 && !(i == p1row && j == p1col))
	      {
		if (rperm[i] != -1)
		  {
		    std::cerr << "Non-permutation row in ";
		    std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
		    std::exit(1);
		  }
		rperm[i] = j;
	      }
	  }
	if (rperm[i] == -1)
	  {
	    std::cerr << "Missing permutation row in ";
	    std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
	    std::exit(1);
	  }
	if (cperm[rperm[i]] != -1)
	  {
	    std::cerr << "Duplicate permutation column in ";
	    std::cerr << "traintracks::mathmatrix_permplus1::mathmatrix_permplus1\n";
	    std::exit(1);
	  }
	cperm[rperm[i]] = i;
      }
  }

  // Expand sparse representation back into a dense square matrix.
  Mat full() const
  {
    const int n = this->dim();
    Mat M(n,n);

    for (int j = 0; j < n; ++j) M(j,rperm[j]) = 1;
    if (p1row >= 0 && p1col >= 0) M(p1row,p1col) = 1;

    return M;
  }

  // Return true if permutation matrix (no "+1").
  bool is_perm() const { return (p1row < 0 || p1col < 0); }

  const Vec& row_perm() const { return rperm; }

  int plus1_row() const { return p1row; }
  int plus1_col() const { return p1col; }

  // Return inverse permutation as column permutation indices.
  const Vec& column_perm() const { return cperm; }

  // Find the order of the permutation part.
  // i.e., smallest integer j such that rperm^j == id.
  int order() const
  {
    const int n = this->dim();
    jlt::vector<int> p(n), pnew(n);
    for (int i = 0; i < n; ++i) p[i] = i;

    if (rperm == p) return 1;

    for (int m = 2; m <= n; ++m)
      {
	for (int i = 0; i < n; ++i)
	  {
	    pnew[rperm[i]] = p[i];
	  }
	if (pnew == rperm) return m;
	p = pnew;
      }

    std::cerr << "Bad permutation in traintrack::mathmatrix_permplus1::order().\n";
    exit(1);
  }

  size_type dim() const { return rperm.size(); }// Always a square matrix.
  size_type rows() const { return dim(); }	// Number of rows.
  size_type columns() const { return dim(); }	// Number of columns.

  friend Mat operator*(const mathmatrix_permplus1& pm, const Mat& A);
  friend Mat operator*(const Mat& A, const mathmatrix_permplus1& pm);
  friend std::ostream& operator<<(std::ostream& strm,
				  const mathmatrix_permplus1& pm);
};

inline
jlt::mathmatrix<int> operator*(const mathmatrix_permplus1& pm,
			       const jlt::mathmatrix<int>& A)
{
  const int n = pm.dim();

  JLT_MATRIX_ASSERT(A.isSquare());
  JLT_MATRIX_ASSERT(n == (int)A.dim());

  jlt::mathmatrix<int> pmA(n,n);

  for (int j = 0; j < n; ++j)
    {
      for (int i = 0; i < n; ++i)
	{
	  pmA(i,j) = A(pm.rperm[i],j);
	}
      if (pm.p1row >= 0 && pm.p1col >= 0) pmA(pm.p1row,j) += A(pm.p1col,j);
    }

  return pmA;
}

inline
jlt::mathmatrix<int> operator*(const jlt::mathmatrix<int>& A,
			       const mathmatrix_permplus1& pm)
{
  const int n = pm.dim();

  JLT_MATRIX_ASSERT(A.isSquare());
  JLT_MATRIX_ASSERT(n == (int)A.dim());

  // Find column permutation.
  // This means that right-multiplication by pm is a little slower than left.
  jlt::vector<int> cperm(pm.column_perm());

  jlt::mathmatrix<int> pmA(n,n);

  for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
	{
	  pmA(i,j) = A(i,cperm[j]);
	}
      if (pm.p1row >= 0 && pm.p1col >= 0) pmA(i,pm.p1col) += A(i,pm.p1row);
    }

  return pmA;
}

inline
std::ostream& operator<<(std::ostream& strm, const mathmatrix_permplus1& pm)
{
  const int n = pm.dim();

  strm << "(";
  for (int i = 0; i < n-1; ++i)
    {
      strm << pm.rperm[i];
      if (n >= 10) strm << " ";
    }
  strm << pm.rperm.back();

  if (pm.p1row >= 0 && pm.p1col >= 0)
    {
      strm << ", ";
      // Print in the Ham&Song form:
      // edge p1row -> rperm[p1row] + p1col
      // Only print p1row and p1col.
      strm << pm.p1row << "->" << pm.p1col;
    }
  strm << ")";

  return strm;
}

inline
std::ostream& printMathematicaForm(std::ostream& strm,
				   const mathmatrix_permplus1& pm)
{
  const int n = pm.dim();

  strm << "{{";
  for (int i = 0; i < n-1; ++i)
    {
      strm << pm.row_perm()[i]+1 << ",";
    }
  strm << pm.row_perm().back()+1 << "}";
  if (pm.plus1_row() >= 0 && pm.plus1_col() >= 0)
    {
      strm << ",{";
      strm << pm.plus1_row()+1 << "," << pm.plus1_col()+1 << "}";
    }
  strm << "}";

  return strm;
}

} // namespace traintracks

#endif // TRAINTRACKS_MATHMATRIX_PERMPLUS1_HPP
