#include <iostream>
#include <jlt/stlio.hpp>
#include "mathmatrix_permplus1.hpp"

using namespace traintracks;
using namespace jlt;
using namespace std;


void testM(const mathmatrix<int>& Mpm)
{
  const int n = Mpm.dim();

  cout << "Matrix pm:\n";
  Mpm.printMatrixForm(cout) << endl;

  mathmatrix_permplus1 pm(Mpm);

  cout << "In permutation+1 form:\n" << pm << endl;

  cout << "\nFull again:\n"; pm.full().printMatrixForm(cout);

  cout << "\nPermutation has order " << pm.order() << endl;

  cout << endl;

  mathmatrix<int> M(identity_matrix<int>(n));
  for (int i = 0; i < n; ++i)
    {
      for (int j = 0; j < n; ++j)
	{
	  M(i,j) = i - 3 + n*j + i*j;
	}
    }

  cout << "Matrix M:\n";
  M.printMatrixForm(cout)<< endl;

  cout << "pm * M: (mathmatrix)\n";
  (Mpm * M).printMatrixForm(cout)<< endl;

  cout << "pm * M: (mathmatrix_permplus1)\n";
  (pm * M).printMatrixForm(cout);

  if (Mpm*M == pm*M)
    {
      cout << "The matrices are equal!\n";
    }
  else
    {
      cout << "The matrices are not equal...\n";
    }

  cout << endl;

  cout << "M * pm: (mathmatrix)\n";
  (M * Mpm).printMatrixForm(cout)<< endl;

  cout << "M * pm: (mathmatrix_permplus1)\n";
  (M * pm).printMatrixForm(cout);

  if (M*Mpm == M*pm)
    {
      cout << "The matrices are equal!\n";
    }
  else
    {
      cout << "The matrices are not equal...\n";
    }
}


int main()
{
  const int n = 6;

  mathmatrix<int> Mpm(n,n);
  Mpm(0,4) = 1;
  Mpm(1,0) = 1;
  Mpm(2,1) = 1;
  Mpm(3,2) = 1;
  Mpm(4,3) = 1;
  Mpm(5,5) = 1;

  testM(Mpm); // Test a plain permutation matrix.

  Mpm(5,4) = 1;	// The +1 element.

  testM(Mpm); // Test a permutation+1 matrix.

  return 0;
}
