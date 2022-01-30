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

#include <bitset>
#include "align.h"

Alignment::ScoringMatrix Alignment::matrix_;

void Alignment::ScoringMatrix::Initialize(int match, int mismatch)
{
	for (size_t i = 0; i < LookupTables::iupac_codes.size(); i++) {
		uint8_t base = LookupTables::iupac_codes[i][0];
		uint8_t flags = LookupTables::BaseFlags(base);
		size_t count = std::bitset<8>(flags).count();
		for (size_t j = 0; j < LookupTables::iupac_codes.size(); j++) {
			uint8_t base2 = LookupTables::iupac_codes[j][0];
			uint8_t flags2 = LookupTables::BaseFlags(base2);
			size_t count2 = std::bitset<8>(flags2).count();
			size_t matches = std::bitset<8>(flags2 & flags).count();
			if (!matches) {
				matrix_[i][j] = mismatch;
			}
			else if (count + count2 > 2) {
				matrix_[i][j] = mismatch + int((match - mismatch) * sqrtf(matches / float(count * count2)) + 0.5f);
			}
			else {
				matrix_[i][j] = match;
			}
		}
		matrix_[i][LookupTables::IUPAC_UNDEFINED_INDEX] = mismatch;
		matrix_[LookupTables::IUPAC_UNDEFINED_INDEX][i] = mismatch;
	}
	matrix_[LookupTables::IUPAC_UNDEFINED_INDEX][LookupTables::IUPAC_UNDEFINED_INDEX] = mismatch;
	size_t n = LookupTables::IupacIndex('N');
}

static const size_t VECTOR_MIN_MATCH = 11;

bool Alignment::Search(const NucleotideSequence& sequence, size_t start, size_t end, size_t max_result, const char* query, int min_percent, Alignment::Result* result_)
{
	InitializeMatrix();

	start -= (start != 0);

	size_t across = strlen(query);
	auto scores = std::make_unique<int[]>(across);
	size_t down = end - start;
	auto seq_down = sequence.cbegin() + start;

	for (size_t x = 1; x <= across; x++)
		scores[across - x] = static_cast<int>(x * ALIGN_GAP_OPEN);

	bool backwards = max_result < sequence.Length();
	int min_score = ComputeSearchMinScore(across, min_percent);

	int prev_score = -1;
	int prev_score_2 = -1;
	ptrdiff_t max_y = max_result;

	Alignment::Result result;
	result.start = -1 - start;
	result.score = -1;

	for (ptrdiff_t y = down - 1; y >= 0; y--) {
		int prev_x = 0;
		int prev_d = 0;
		const int* row_matrix = GetSearchMatrix(seq_down[y]);

		for (ptrdiff_t x = across - 1; x >= 0; x--) {
			int match = prev_d + row_matrix[LookupTables::IupacIndex(query[x])];
			prev_d = scores[x];
			scores[x] = std::max(match, std::max(prev_x, prev_d) + ALIGN_GAP_OPEN);
			prev_x = scores[x];
		}

		if (prev_x <= prev_score && prev_score >= prev_score_2 && prev_score >= min_score && (!backwards || y < max_y)) {
			result.start = y + 1;
			result.score = prev_score;
			if (backwards)
				break;
		}

		prev_score_2 = prev_score;
		prev_score = prev_x;
	}

	if (result.start < 0 && !start && prev_score >= min_score) {
		result.start = 0;
		result.score = prev_score;
	}

	result.start += start;
	*result_ = result;
	return result.score >= min_score;
}

static const int N_QUALITY = 5;

bool Alignment::VectorSearch5(const NucleotideSequence& sequence, const char* query, int min_percent, size_t min_match, Alignment::Result* result_)
{
	InitializeMatrix();

	Result result = { -1, -1 };

	size_t across = strlen(query);
	if (!across || !sequence.Length()) {
		*result_ = result;
		return false;
	}

	auto scores = std::make_unique<Node[]>(across);

	const int* row_matrix = GetSearchMatrix(sequence[0]);
	int quality = sequence.QualityOrDefault(0);
	if (LookupTables::Uppercase(sequence[0]) == 'N') {
		quality = N_QUALITY;
	}
	for (size_t x = 0; x < across; x++) {
		scores[x].score = row_matrix[LookupTables::IupacIndex(query[x])] * quality;
		scores[x].quality = quality;
		scores[x].gap_penalty = ALIGN_GAP_OPEN;
	}

	int min_score = ComputeSearchMinScore(min_match, min_percent);
	size_t first = std::min(across, min_match) - 1;

	for (size_t y = 1; y < sequence.Length(); y++) {
		row_matrix = GetSearchMatrix(sequence[y]);
		quality = sequence.QualityOrDefault(y);
		if (LookupTables::Uppercase(sequence[y]) == 'N')
			quality = N_QUALITY;
		Node diagonal_node = scores[0];
		scores[0].score = row_matrix[LookupTables::IupacIndex(query[0])] * quality;
		scores[0].quality = quality;

		for (size_t x = 1; x < across; x++) {
			int score_diagonal = diagonal_node.score + row_matrix[LookupTables::IupacIndex(query[x])] * quality;
			int match_quality = diagonal_node.quality + quality;
			diagonal_node = scores[x];
			int score_across = scores[x - 1].score + scores[x - 1].gap_penalty * quality;

			scores[x].score += scores[x].gap_penalty * quality;
			scores[x].quality += quality;
			scores[x].gap_penalty = ALIGN_GAP_EXTEND;

			if (scores[x].score < score_across || x == across - 1) {
				scores[x].score = score_across;
				scores[x].quality = scores[x - 1].quality + quality;
			}
			if (scores[x].score < score_diagonal) {
				scores[x].score = score_diagonal;
				scores[x].quality = match_quality;
				scores[x].gap_penalty = ALIGN_GAP_OPEN;
			}
		}

		if (y >= first) {
			int score = scores[across - 1].score / scores[across - 1].quality;
			if (result.score < score && score >= min_score) {
				result.score = score;
				result.start = y;
			}
		}
	}

	*result_ = result;
	return result.score >= min_score;
}

bool Alignment::VectorSearch3(const NucleotideSequence& sequence, size_t start, const char* query, int min_percent, size_t min_match, Alignment::Result* result_)
{
	InitializeMatrix();

	Result result = { -1, -1 };

	size_t across = strlen(query);
	ptrdiff_t down = sequence.Length() - start;

	if (!across || down <= 0) {
		*result_ = result;
		return false;
	}

	auto seq_down = sequence.cbegin() + start;
	auto scores = std::make_unique<Node[]>(across);

	const int* row_matrix = GetSearchMatrix(seq_down[down - 1]);
	int quality = sequence.QualityOrDefault(start + down - 1);
	if (LookupTables::Uppercase(seq_down[down - 1]) == 'N')
		quality = N_QUALITY;
	for (size_t x = 0; x < across; x++) {
		scores[x].score = row_matrix[LookupTables::IupacIndex(query[x])] * quality;
		scores[x].quality = quality;
		scores[x].gap_penalty = ALIGN_GAP_OPEN;
	}

	int min_score = ComputeSearchMinScore(min_match, min_percent);
	ptrdiff_t first = down - std::min(across, min_match);

	for (ptrdiff_t y = down - 2; y >= 0; --y) 	{
		row_matrix = GetSearchMatrix(seq_down[y]);
		quality = sequence.QualityOrDefault(start + y);
		if (LookupTables::Uppercase(seq_down[y]) == 'N')
			quality = N_QUALITY;
		Node diagonal_node = scores[across - 1];
		scores[across - 1].score = row_matrix[LookupTables::IupacIndex(query[across - 1])] * quality;
		scores[across - 1].quality = quality;

		for (ptrdiff_t x = across - 2; x >= 0; --x) {
			int diagonal = diagonal_node.score + row_matrix[LookupTables::IupacIndex(query[x])] * quality;
			int match_quality = diagonal_node.quality + quality;
			diagonal_node = scores[x];
			int iAcross = scores[x + 1].score + scores[x + 1].gap_penalty * quality;

			scores[x].score += scores[x].gap_penalty * quality;
			scores[x].quality += quality;
			scores[x].gap_penalty = ALIGN_GAP_EXTEND;

			if (scores[x].score < iAcross || x == 0) {
				scores[x].score = iAcross;
				scores[x].quality = scores[x + 1].quality + quality;
			}
			if (scores[x].score < diagonal) {
				scores[x].score = diagonal;
				scores[x].quality = match_quality;
				scores[x].gap_penalty = ALIGN_GAP_OPEN;
			}
		}

		if (y <= first) {
			int score = scores[0].score / scores[0].quality;
			if (result.score < score && score >= min_score) {
				result.score = score;
				result.start = y;
			}
		}
	}

	*result_ = result;
	return result.score >= min_score;
}