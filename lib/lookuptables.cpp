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

const std::array<const std::string, 15> LookupTables::iupac_codes = {
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

std::array<char, 256> LookupTables::uppercase, LookupTables::lowercase;
std::array<uint8_t, 256> LookupTables::char_index;
std::array<uint8_t, 256> LookupTables::iupac_index;
std::array<char, 256> LookupTables::complement;
std::array<uint8_t, 256> LookupTables::base_flags;
std::array<uint32_t, LookupTables::CHAR_UNKNOWN + 1> LookupTables::amino_acid_redundant_matrix;

void LookupTables::Initialize()
{
	const std::array<const std::string, 15> aa_iupac_codes_ = {
		"B:DN",
		"X:ACDEFGHIKLMNPQRSTVWY",
		"Z:EQ"
	};
	const std::string bases = "TCAG";

	for (int i = 0; i < uppercase.size(); ++i) {
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
	iupac_index.fill(IUPAC_UNDEFINED_INDEX);
	for (size_t i = 0; i < iupac_codes.size(); ++i) {
		auto& str = iupac_codes[i];
		uint8_t flags = 0;
		for (auto it = str.cbegin() + 2; it != str.cend(); ++it) {
			flags |= uint8_t(1) << bases.find(*it, 0);
		}
		uint8_t base = str[0];
		base_flags[base] = flags;
		base_flags[lowercase[base]] = flags;
		iupac_index[base] = uint8_t(i);
		iupac_index[lowercase[base]] = uint8_t(i);
	}
	base_flags[(unsigned char)'U'] = base_flags[(unsigned char)'T'];
	base_flags[(unsigned char)'u'] = base_flags[(unsigned char)'T'];
	iupac_index['U'] = iupac_index['T'];
	iupac_index['u'] = iupac_index['T'];

	for (size_t i = 0; i < amino_acid_redundant_matrix.size(); ++i)
		amino_acid_redundant_matrix[i] = (1u << i);
	for (size_t i = 0; i < aa_iupac_codes_.size(); ++i) {
		auto& str = aa_iupac_codes_[i];
		uint32_t flags = 0;
		for (auto it = str.cbegin() + 2; it != str.cend(); ++it) {
			flags |= 1u << CharIndex(*it);
		}
		size_t aa_index = CharIndex(str[0]);
		amino_acid_redundant_matrix[aa_index] |= flags;
	}
	amino_acid_redundant_matrix[CHAR_UNKNOWN] = 0;
}
