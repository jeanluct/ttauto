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

#ifndef TRAINTRACKS_GRAPH_ISO_HPP
#define TRAINTRACKS_GRAPH_ISO_HPP

#include <vector>

namespace traintracks {

class traintrack;

namespace graph_iso {

struct canonical_witness
{
  bool valid = false;
  std::vector<int> multigon_order;
  std::vector<int> multigon_rank;
  std::vector<int> prong_shift;
};

// Orientation-preserving isotopy test using blown-up prong graph.
//
// Phase-1 semantics:
// - preserves branch incidence and multiplicity at each multigon prong,
// - preserves cyclic prong order up to rotation (for k>=3),
// - preserves puncture flags and multigon labels.
//
// This is a direct structural matcher, not a coding-sequence generator.
bool is_isotopic_oriented(const traintrack& lhs, const traintrack& rhs);

// Deterministic orientation-preserving witness for canonical relabeling.
canonical_witness canonical_witness_oriented(const traintrack& tt);

} // namespace graph_iso
} // namespace traintracks

#endif // TRAINTRACKS_GRAPH_ISO_HPP
