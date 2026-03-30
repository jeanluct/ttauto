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

#ifndef TRAINTRACKS_BUILD_HPP
#define TRAINTRACKS_BUILD_HPP

#include <jlt/vector.hpp>
#include "traintracks/traintrack.hpp"

namespace traintracks {

// Build train tracks for fixed puncture count N and fixed number N2 of
// punctured bigons, enumerating compatible higher-prong strata.
jlt::vector<traintrack> build_traintrack_list(const int N, const int N2 = 0);

// Build train tracks by sweeping punctured-bigon count N2 from 0 to N-3,
// aggregating build_traintrack_list(N-N2,N2) across that range.
jlt::vector<traintrack> build_traintrack_list_sweep_bigons(const int N);

} // namespace traintracks

#endif // TRAINTRACKS_BUILD_HPP
