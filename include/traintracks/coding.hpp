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

#ifndef TRAINTRACKS_CODING_HPP
#define TRAINTRACKS_CODING_HPP

#include "traintracks/traintrack.hpp"

namespace traintracks {
namespace coding {

traintrack::intVec coding_from_monogon(const traintrack& tt,
                                       const int mono,
                                       const int dir = 1);

traintrack::intVec coding(const traintrack& tt, const int dir = 1);

traintrack::intVec canonical_coding_from_monogon(const traintrack& tt,
                                                 const int mono,
                                                 const int dir = 1);

traintrack::intVec canonical_coding(const traintrack& tt, const int dir = 1);

} // namespace coding
} // namespace traintracks

#endif // TRAINTRACKS_CODING_HPP
