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

#include "sequence.h"

NucleotideSequence::NucleotideSequence()
{
}

NucleotideSequence::NucleotideSequence(const base_type* seq, size_t length, const char* name)
	: sequence_(seq, length)
{
	SetName(name);
}

void NucleotideSequence::ConstructQuality()
{
	if (quality_.Length() < sequence_.Length()) {
		size_t start_pos = quality_.Length();
		size_t length = sequence_.Length() - quality_.Length();
		quality_.ResizeSubsequence(start_pos, 0, length);
		InitializeQuality(start_pos, length);
	}
	else if (quality_.Length() > sequence_.Length()) {
		quality_.Truncate(sequence_.Length());
	}
}

void NucleotideSequence::SetName(const char* name)
{
	if (name)
		description_.name = name;
}

bool NucleotideSequence::HasValidQuality() const
{
	for (auto q = QualityBegin(); q < QualityEnd(); ++q) {
		if (*q > 1)
			return true;
	}
	return false;
}

NucleotideSequence::BaseCounts NucleotideSequence::ComputeBaseCounts() const
{
	BaseCounts bc = { 0, 0, 0, 0, 0 };
	for (auto base = cbegin(); base != cend(); ++base) {
		switch (LookupTables::Uppercase(*base)) {
		case 'A': ++bc.a; break;
		case 'C': ++bc.c; break;
		case 'G': ++bc.g; break;
		case 'U':
		case 'T': ++bc.t; break;
		default: ++bc.other;
		}
	}
	return bc;
}

float NucleotideSequence::ComputePercentGC() const
{
	size_t gc = 0, at = 0;
	for (auto base = cbegin(); base != cend(); ++base) {
		switch (LookupTables::Uppercase(*base)) {
		case 'A':
		case 'T':
		case 'U':
		case 'W': ++at; break;
		case 'C':
		case 'G':
		case 'S': ++gc; break;
		}
	}
	return gc * 100.0f / (gc + at + !(gc + at));
}

float NucleotideSequence::ComputeSpacing() const
{
	if (Length() > 1 && traces_) {
		peak_type total = 0;
		for (auto peak = PeakBegin() + 1; peak != PeakEnd(); ++peak)
			total += peak[0] - peak[-1];
		return (float)total / float(Length() - 1);
	}
	return 0.0f;
}

void NucleotideSequence::Replace(size_t start_pos, size_t old_length, const NucleotideSequence& source)
{
	sequence_.Replace(start_pos, old_length, source.sequence_);
	quality_.Replace(start_pos, old_length, source.quality_);
	if (traces_ && !source.Empty()) {
		if (source.traces_->peaks) {
			traces_->peaks.Replace(start_pos, old_length, source.traces_->peaks);
		}
		else {
			assert(false);
			traces_->peaks.ResizeSubsequence(start_pos, old_length, source.Length());
			// Zero any extra peaks.
			// If this function is used only for undo/redo the source will always have peaks if the destination
			// has peaks. Editing a peak-containing sequence outside of a trace editor can result in zeroing
			// new peaks, which should normally be avoided if it will be saved back into a trace file.
			if (source.Length() > old_length)
				traces_->peaks.FillSubsequence(start_pos + old_length, source.Length() - old_length, 0);
		}
	}
}

void NucleotideSequence::Replace(size_t pos, size_t old_length, base_type c, quality_type q, peak_type peak)
{
	sequence_.Replace(pos, old_length, c);
	quality_.Replace(pos, old_length, q);
	if (traces_)
		traces_->peaks.Replace(pos, old_length, peak);
}

void NucleotideSequence::DeleteSubsequence(size_t start_pos, size_t length)
{
	sequence_.DeleteSubsequence(start_pos, length);
	quality_.DeleteSubsequence(start_pos, length);
	if (traces_)
		traces_->peaks.DeleteSubsequence(start_pos, length);
}

void NucleotideSequence::InitializeQuality(size_t start_pos, size_t length)
{
	quality_.FillSubsequence(start_pos, length, 0);
}

void NucleotideSequence::ReverseComplement(peak_type trace_length_max)
{
	size_t back = Length() - 1;
	for (size_t fwd = 0; fwd <= back; ++fwd, --back) {
		base_type c = LookupTables::Complement(sequence_[fwd]);
		sequence_[fwd] = LookupTables::Complement(sequence_[back]);
		sequence_[back] = c;
		std::swap(quality_[fwd], quality_[back]);
	}
	if (traces_) {
		SeqContainer<peak_type>& peaks = traces_->peaks;
		back = Length() - 1;
		for (size_t fwd = 0; fwd <= back; ++fwd, --back) {
			peak_type i = trace_length_max - 1 - peaks[fwd];
			peaks[fwd] = trace_length_max - 1 - peaks[back];
			peaks[back] = i;
		}
	}
}

bool NucleotideSequence::MatchSequence(size_t pos, const std::string& query, bool both_strands) const
{
	if (MatchForward(pos, query.c_str(), query.length()))
		return true;
	if (both_strands)
		return MatchReverse(pos, query.c_str(), query.length());
	return false;
}

size_t NucleotideSequence::SearchSequenceForward(size_t start_pos, const std::string& query, bool both_strands) const
{
	assert(start_pos <= Length() && start_pos + query.length() <= Length());
	size_t end = Length() - query.length();
	for (size_t i = start_pos; i <= end; ++i)
		if (MatchSequence(i, query, both_strands))
			return i;
	return NOT_FOUND;
}

size_t NucleotideSequence::SearchSequenceBackward(size_t start_pos, const std::string& query, bool both_strands) const
{
	assert(start_pos <= Length() && start_pos + query.length() <= Length());
	for (ptrdiff_t i = std::min(start_pos, Length() - query.length()); i >= 0; --i)
		if (MatchSequence(i, query, both_strands))
			return i;
	return NOT_FOUND;
}

size_t NucleotideSequence::FindNextN(size_t start_pos) const
{
	for (size_t i = start_pos; i < Length(); ++i)
		if (LookupTables::Uppercase(sequence_[i]) == 'N')
			return i;
	return NOT_FOUND;
}

size_t NucleotideSequence::FindNextRedundant(size_t start_pos) const
{
	for (size_t i = start_pos; i < Length(); ++i)
		if (IsRedundant(i))
			return i;
	return NOT_FOUND;
}

bool NucleotideSequence::IsRedundant(size_t pos) const
{
	auto flags = LookupTables::BaseFlags(sequence_[pos]);
	return (flags & (flags - 1)) != 0;
}
