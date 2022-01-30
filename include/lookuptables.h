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

#pragma once

#include <array>
#include <string>

class LookupTables
{
public:
	static const size_t IUPAC_UNDEFINED_INDEX = 15;

	static const std::array<const std::string, 15> iupac_codes;

	static inline char Uppercase(char base) {
		return uppercase[(unsigned char)base];
	}

	static inline char Lowercase(char base) {
		return lowercase[(unsigned char)base];
	}

	static inline uint8_t CharIndex(char base) {
		return char_index[(unsigned char)base];
	}

	static inline uint8_t IupacIndex(char base) {
		return iupac_index[(unsigned char)base];
	}

	static inline char Complement(char base) {
		return complement[(unsigned char)base];
	}

	static inline uint8_t BaseFlags(char base) {
		return base_flags[(unsigned char)base];
	}

	static inline uint8_t BaseFlagsComplement(char base) {
		return base_flags[(unsigned char)complement[(unsigned char)base]];
	}

	static inline bool BaseMatch(char base, char query) {
		return (BaseFlags(base) & BaseFlags(query)) == BaseFlags(base);
	}

	static inline bool AminoAcidMatch(char base, char query) {
		// CharIndex() can't handle the stop codon '*' so test identity too.
		return base == query || (amino_acid_redundant_matrix[CharIndex(base)] & (1u << CharIndex(query))) != 0;
	}

private:
	static class static_initializer
	{
	public:
		static_initializer() { Initialize(); }
	} initializer_;

	static void Initialize();

	static const char CHAR_UNKNOWN = 26;

	static std::array<char, 256> uppercase, lowercase;
	static std::array<uint8_t, 256> char_index;
	static std::array<uint8_t, 256> iupac_index;
	static std::array<char, 256> complement;
	static std::array<uint8_t, 256> base_flags;
	static std::array<uint32_t, CHAR_UNKNOWN + 1> amino_acid_redundant_matrix;
};

