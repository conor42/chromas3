// Alignment functions.
// 
// Copyright 2022 Conor N. McCarthy
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

#include "sequence.h"

class Alignment
{
private:
	static const size_t BASE_INDEX_COUNT = 17;
	static const int ALIGN_MATCH = 2;
	static const int ALIGN_MISMATCH = -6;
	static const int ALIGN_GAP_OPEN = -4;
	static const int ALIGN_GAP_EXTEND = ALIGN_GAP_OPEN * 2;

	class ScoringMatrix
	{
	public:
		ScoringMatrix() {
			matrix_[0][0] = 0;
		}

		void Initialize(int match, int mismatch);

		int* operator[](size_t i) {
			return matrix_[i];
		}

#ifdef _DEBUG
		void debug_print();
#endif

	private:
		int matrix_[BASE_INDEX_COUNT][BASE_INDEX_COUNT];
	};

	struct Node
	{
		int score;
		int quality;
		int gap_penalty;
	};

public:
	struct Result
	{
		int score;
		ptrdiff_t start;
	};

	static bool Search(const NucleotideSequence& sequence, size_t start, size_t end, size_t max_result, const char* query, int min_percent, Result* result);
	static bool VectorSearch5(const NucleotideSequence& sequence, const char* query, int min_percent, size_t min_match, Result* result);
	static bool VectorSearch3(const NucleotideSequence& sequence, size_t start, const char* query, int min_percent, size_t min_match, Result* result);

private:
	static void InitializeMatrix() {
		matrix_.Initialize(ALIGN_MATCH, ALIGN_MISMATCH);
	}

	static int* GetSearchMatrix(char c) {
		return matrix_[LookupTables::IupacIndex(c)];
	}

	static int ComputeSearchMinScore(ptrdiff_t length, int min_percent)
	{
		return static_cast<int>(length * ALIGN_GAP_OPEN + length * (ALIGN_MATCH - ALIGN_GAP_OPEN) * min_percent / 100);
	}

	static ScoringMatrix matrix_;
};