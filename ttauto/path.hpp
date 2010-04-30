#ifndef PATH_HPP
#define PATH_HPP

#include <iostream>
#include <jlt/vector.hpp>

namespace traintracks {

// A path is a vector of integers with values between m0 and m1.

// operator++ is defined to step through all possible paths.

// TODO:
//
// - Check content between [m0,m1]
// - operator-- to come back.
// - Should operator++ really return false at the end of a path?
//   Breaks convention.  Could reset path, and define == 0?
// - Convert a path into a long unsigned int and vice versa?

class path : public jlt::vector<int>
{
  // Branches in path take values from m0 to m1.
  const int m0;
  const int m1;
public:
  path(const int m0_, const int m1_, const int n_ = 0)
    : jlt::vector<int>(n_), m0(m0_), m1(m1_)
  {
  }

  // Copy constructor.
  path(const path& pp) : jlt::vector<int>(pp), m0(pp.m0), m1(pp.m1)
  {
  }

  path& operator=(const path& pp)
  {
    if (m0 != pp.m0 || m1 != pp.m1)
      {
	std::cerr << "Cannot be equated in path::operator=.\n";
	std::exit(1);
      }

    return (path&)jlt::vector<int>::operator=(pp);
  }

  bool operator==(const path& pp) const
  {
    return (m0 == pp.m0 && m1 == pp.m1 && std::operator==(*this,pp));
  }

  // Cycle through all paths of fixed length, letting values at each
  // slot tun from m0 to m1.
  bool operator++()
  {
    // Increment the path.
    const int pathlength = size();

    bool incr = false;

    for (int b = 0; b < pathlength; ++b)
      {
	// The sequence that branches cycle through is m0 ... m1
	if ((*this)[b] != m1)
	  {
	    incr = true;
	    // Increment the branch.
	    ++((*this)[b]);
	    break;
	  }
	else
	  {
	    // Otherwise reset branch at that position, and let
	    // the loop move on to the next position.
	    (*this)[b] = m0;
	  }
      }
    // If nothing was changed, we're done.
    return incr;
  }

  // Almost the same as the vector < for paths of equal lengths, but
  // shorter paths are smaller than longer ones.
  bool operator<(const path& pp) const
  {
    if (m0 != pp.m0 || m1 != pp.m1)
      {
	std::cerr << "Cannot be compared in path::operator=.\n";
	std::exit(1);
      }
    if (this->size() != pp.size()) return (this->size() < pp.size());

    return std::operator<(*this,pp);
  }

  bool operator>(const path& pp) const { return (pp < *this); }

  bool operator<=(const path& pp) const { return !(pp < *this); }

  bool operator>=(const path& pp) const { return !(pp > *this); }

  friend std::ostream& operator<<(std::ostream& strm, const path& pp);
};

// Print a path.
std::ostream& operator<<(std::ostream& strm, const path& pp)
{
  if (pp.empty()) return strm;

  for (path::const_iterator i = pp.begin(); i != pp.end(); ++i)
    {
      strm << *i;
      // Appropriate if 0 <= m0 < m1 < 10.
      // Otherwise leave a space between path elements.
      if (!(pp.m0 >= 0 && pp.m1 < 10) && i != pp.end()-1) strm << " ";
    }

  return strm;
}

} // namespace traintracks

#endif // PATH_HPP
