// Genetic code tables and conversion methods.
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

#include <algorithm>
#include <array>
#include <cassert>
#include <string>
#include <vector>
#include "lookuptables.h"

class GeneticCodes
{
public:
	struct Codon
	{
		char amino_acid;
		bool can_start;
	};

private:
	// Lookup table size is 1 bit per base for 3 bases
	static const size_t TABLE_SIZE = 1 << (4 * 3);

	struct CodonTable
	{
		std::string name;
		std::array<Codon, TABLE_SIZE> table;

		CodonTable(const char* line);
	};

public:
	static const int STANDARD = 0;
	static const char AA_D_OR_N = 'B';
	static const char AA_E_OR_Q = 'Z';
	static const char AA_UNKNOWN = 'X';

	// Iterator to enumerate names of all available genetic codes
	class Iterator
	{
		using value_type = const std::string;
		using pointer_type = const value_type*;
		using reference_type = const value_type&;

	public:
		Iterator(std::vector<CodonTable>::const_iterator it)
			: iterator_(it) {}

		reference_type operator*() const {
			return iterator_->name;
		}

		pointer_type operator->() const {
			return &iterator_->name;
		}

		Iterator& operator++() {
			++iterator_;
			return *this;
		}

		Iterator operator++(int) {
			auto tmp = *this;
			++iterator_;
			return tmp;
		}

		bool operator==(const Iterator& rval) const {
			return iterator_.operator==(rval.iterator_);
		}

		bool operator!=(const Iterator& rval) const	{
			return iterator_.operator!=(rval.iterator_);
		}

	private:
		std::vector<CodonTable>::const_iterator iterator_;
	};

	static Iterator cbegin() {
		return Iterator(codon_tables_.cbegin());
	}

	static Iterator cend() {
		return Iterator(codon_tables_.cend());
	}

	static int Search(const char* name);

	static Codon TranslateForward(const char* sequence, size_t pos, size_t length, size_t genetic_code)	{
		assert(genetic_code < codon_tables_.size());

		if (genetic_code >= codon_tables_.size())
			return { AA_UNKNOWN, false };

		unsigned int i = LookupTables::BaseFlags(sequence[pos + 2])
			+ (LookupTables::BaseFlags(sequence[pos + 1]) << 4)
			+ (LookupTables::BaseFlags(sequence[pos]) << 8);
		return codon_tables_[genetic_code].table[i];
	}

	static Codon TranslateReverseComplement(const char* sequence, size_t pos, size_t length, size_t genetic_code) {
		assert(genetic_code < codon_tables_.size());

		if (genetic_code >= codon_tables_.size())
			return { AA_UNKNOWN, false };

		unsigned int i = LookupTables::BaseFlagsComplement(sequence[pos])
			+ (LookupTables::BaseFlagsComplement(sequence[pos + 1]) << 4)
			+ (LookupTables::BaseFlagsComplement(sequence[pos + 2]) << 8);
		return codon_tables_[genetic_code].table[i];
	}

private:
	static class static_initializer
	{
	public:
		static_initializer() { Initialize(); }
	} initializer_;

	static void Initialize();

	static std::vector<CodonTable> codon_tables_;
};