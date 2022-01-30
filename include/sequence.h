// Nucleotide sequence class.
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

#include <string>
#include "lookuptables.h"
#include "seqcontainer.h"

class NucleotideSequence
{
public:
	static const size_t NOT_FOUND = (size_t)-1;
	static const size_t TRACE_COUNT = 4;
	static const uint8_t DEFAULT_BASE_QUALITY = 40;

	using base_type = char;
	using quality_type = uint8_t;
	using peak_type = int32_t;
	using trace_type = int32_t;

	struct Description
	{
		std::string name;
	};

	struct BaseCounts
	{
		size_t a, c, g, t, other;
	};

	struct Traces
	{
		SeqContainer<peak_type> peaks;
		std::unique_ptr<trace_type[]> heights[4];
		size_t trace_length;
	};

	NucleotideSequence();

	NucleotideSequence(const base_type* sequence, size_t length = (size_t)-1, const char* name = nullptr);

	template<class BaseIterator, class QualityIterator>
	NucleotideSequence(BaseIterator begin_base, BaseIterator end_base, QualityIterator begin_qual, QualityIterator end_qual, const char* name = nullptr)
		: sequence_(begin_base, end_base),
		quality_(begin_qual, end_qual)
	{
		ConstructQuality();
		SetName(name);
	}

	template<class PeakIterator, class TraceIterator>
	void LoadTraces(PeakIterator begin_peak, PeakIterator end_peak, const TraceIterator begin_trace[TRACE_COUNT], const TraceIterator end_trace[TRACE_COUNT])
	{
		traces_ = std::make_unique<Traces>();
		traces_->peaks = SeqContainer<peak_type>(begin_peak, end_peak);
		traces_->trace_length = 0;

		for (size_t i = 0; i < TRACE_COUNT; ++i)
			traces_->trace_length = std::max<size_t>(traces_->trace_length, end_trace[i] - begin_trace[i]);

		for (size_t i = 0; i < TRACE_COUNT; ++i) {
			std::unique_ptr<trace_type[]>& heights = traces_->heights[i];
			heights = std::make_unique<trace_type[]>(traces_->trace_length);
			size_t length = 0;
			for (TraceIterator y = begin_trace[i]; y != end_trace[i]; ++y)
				heights[length++] = *y;
			while(length < traces_->trace_length)
				heights[length++] = 0;
		}
	}

	NucleotideSequence& operator=(const NucleotideSequence& rval) = default;

	NucleotideSequence& operator=(NucleotideSequence&& rval) = default;

	void SetName(const char* name);

	size_t Length() const {
		return sequence_.Length();
	}

	bool Empty() const {
		return sequence_.Length() == 0;
	}

	base_type operator[](size_t pos) const {
		return sequence_[pos];
	}

	const base_type* cbegin() const {
		return sequence_.cbegin();
	}

	const base_type* cend() const {
		return sequence_.cend();
	}

	const quality_type* QualityBegin() const {
		return quality_.cbegin();
	}

	const quality_type* QualityEnd() const {
		return quality_.cend();
	}

	const peak_type* PeakBegin() const {
		return traces_ ? traces_->peaks.cbegin() : nullptr;
	}

	const peak_type* PeakEnd() const {
		return traces_ ? traces_->peaks.cend() : nullptr;
	}

	bool HasTraces() const {
		return traces_.operator bool();
	}

	bool HasValidQuality() const;

	quality_type QualityOrDefault(size_t pos) const {
		return quality_[pos] > 1 ? quality_[pos] : DEFAULT_BASE_QUALITY;
	}

	BaseCounts ComputeBaseCounts() const;

	float ComputePercentGC() const;

	float ComputeSpacing() const;

	void Replace(size_t start_pos, size_t old_length, const NucleotideSequence& source);

	void Replace(size_t pos, size_t old_length, base_type c, quality_type quality, peak_type peak);

	void DeleteSubsequence(size_t start_pos, size_t length);

	void ReverseComplement(peak_type trace_length_max);

	size_t SearchSequenceForward(size_t start_pos, const std::string& query, bool both_strands) const;

	size_t SearchSequenceBackward(size_t start_pos, const std::string& query, bool both_strands) const;

	size_t SearchByAlignmentFwd(size_t start_pos, const std::string& query, int min_percent) const;

	size_t SearchByAlignmentBack(size_t start_pos, const std::string& query, int min_percent) const;

	size_t FindInTranslation(size_t start_pos, bool backwards, const std::string& query, int genetic_code) const;

	size_t FindNextN(size_t start_pos) const;

	size_t FindNextRedundant(size_t start_pos) const;

	bool IsRedundant(size_t pos) const;

	size_t ComputeQualityStart(unsigned int window, int quality) const;

	size_t ComputeQualityEnd(size_t start_pos, unsigned int window, int quality) const;

private:

	static const quality_type MAX_QUALITY = 93;

	void ConstructQuality();

	void InitializeQuality(size_t start_pos, size_t length);

	bool MatchSequence(size_t pos, const std::string& query, bool both_strands) const;

	bool MatchForward(size_t pos, const base_type* query, size_t length) const {
		for (size_t i = 0; i < length; ++pos, ++i) {
			if (!LookupTables::BaseMatch(sequence_[pos], query[i]))
				return false;
		}
		return true;
	}

	bool MatchReverse(size_t pos, const base_type* query, size_t length) const {
		for (ptrdiff_t i = length - 1; i >= 0; ++pos, --i) {
			if (!LookupTables::BaseMatch(sequence_[pos], LookupTables::Complement(query[i])))
				return false;
		}
		return true;
	}

	bool MatchTranslation(size_t pos, const base_type* query, int genetic_code) const;

	base_type UppercaseBase(size_t i) const {
		return LookupTables::Uppercase(sequence_[i]);
	}

	Description description_;
	SeqContainer<base_type> sequence_;
	SeqContainer<quality_type> quality_;
	std::unique_ptr<Traces> traces_;
};

