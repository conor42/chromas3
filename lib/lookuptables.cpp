// Lookup tables for case conversion, complement, base bitflag representation, and char alphabetic index.
// 
// Copyright 2021 Conor N. McCarthy
//
// This file is part of Chromas 3.
//
// Chromas 3 is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Chromas 3 is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Chromas 3. If not, see < https://www.gnu.org/licenses/>.

#include "lookuptables.h"

LookupTables::static_initializer LookupTables::initializer_;

std::array<char, 256> LookupTables::uppercase, LookupTables::lowercase;
std::array<uint8_t, 256> LookupTables::char_index;
std::array<char, 256> LookupTables::complement;
std::array<uint8_t, 256> LookupTables::base_flags;

void LookupTables::Initialize()
{
	const std::string bases = "TCAG";
	const std::array<const std::string, 15> iupac_codes = {
		"A:A",
		"B:CGT",
		"C:C",
		"D:AGT",
		"G:G",
		"H:ACT",
		"K:GT",
		"M:AC",
		"N:ACGT",
		"R:AG",
		"S:CG",
		"T:T",
		"V:ACG",
		"W:AT",
		"Y:CT"
	};

	for (int i = 0; i < uppercase.size(); i++) {
		uppercase[i] = char(std::toupper(i));
		lowercase[i] = char(std::tolower(i));
	}

	char_index.fill(CHAR_UNKNOWN);

	// Complement is undefined for invalid nucleotide codes so leave them unchanged
	for(int i = 0; i < complement.size(); ++i)
		complement[i] = char(i);
	for (size_t i = 0; i < 26; ++i) {
		//                                           ABCDEFGHIJKLMNOPQRSTUVWXYZ
		static const char alphabetic_complement[] = "TVGHEFCDIJMLKNOPQYSAABWXRZ";
		unsigned char c = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i];
		complement[c] = alphabetic_complement[i];
		complement[std::tolower(c)] = std::tolower(alphabetic_complement[i]);
		char_index[c] = (uint8_t)i;
		char_index[std::tolower(c)] = (uint8_t)i;
	}

	base_flags.fill(0);
	for (auto& str : iupac_codes) {
		uint8_t flags = 0;
		for (auto it = str.cbegin() + 2; it != str.cend(); ++it) {
			flags |= uint8_t(1) << bases.find(*it, 0);
		}
		char base = str[0];
		base_flags[(unsigned char)base] = flags;
		base_flags[(unsigned char)Lowercase(base)] = flags;
	}
	base_flags[(unsigned char)'U'] = base_flags[(unsigned char)'T'];
	base_flags[(unsigned char)'u'] = base_flags[(unsigned char)'T'];
}
