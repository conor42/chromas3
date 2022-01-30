// Tests for sequence classes.
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

#include "catch_amalgamated.hpp"
#include "sequence.h"
#include "ab1file.h"
#include "scffile.h"

template<typename T>
static int elemcmp(const T* first, const void* second, size_t count)
{
    return memcmp(first, second, count * sizeof(T));
}

TEST_CASE("Sequence modification", "[modification]")
{
    static const char seq[] = "ACTG";
    static const uint8_t qual[] = "\x3F\x3A\x3D\x30";
    static const int32_t peaks[] = { 5,17,29,42 };
    static const int32_t* traces[4] = { nullptr, nullptr, nullptr, nullptr };
    NucleotideSequence sequence(seq + 0, seq + 4, qual + 0, qual + 4, nullptr);
    sequence.LoadTraces(peaks + 0, peaks + 4, traces, traces);
    NucleotideSequence sequence2(seq + 0, seq + 4, qual + 0, qual + 4, nullptr);
    sequence2.LoadTraces(peaks + 0, peaks + 4, traces, traces);
    NucleotideSequence empty_sequence;

    REQUIRE(sequence.Length() == 4);
    REQUIRE(empty_sequence.Length() == 0);

    SECTION("Replace subsequence")
    {
        sequence.Replace(1, 1, 'R', 1, 11);
        REQUIRE(sequence.Length() == 4);
        REQUIRE(elemcmp(sequence.cbegin(), "ARTG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x1\x3D\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,11,29,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Insert subsequence")
    {
        sequence.Replace(1, 0, 'R', 1, 11);
        REQUIRE(sequence.Length() == 5);
        REQUIRE(elemcmp(sequence.cbegin(), "ARCTG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x1\x3A\x3D\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,11,17,29,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Insert subsequence to expand beyond reserve")
    {
        sequence.Replace(3, 0, sequence2);
        REQUIRE(sequence.Length() == 8);
        REQUIRE(elemcmp(sequence.cbegin(), "ACTACTGG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x3A\x3D\x3F\x3A\x3D\x30\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,17,29,5,17,29,42,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Insert subsequence at head")
    {
        sequence.Replace(0, 0, 'R', 1, 3);
        REQUIRE(sequence.Length() == 5);
        REQUIRE(elemcmp(sequence.cbegin(), "RACTG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x1\x3F\x3A\x3D\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 3,5,17,29,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Insert subsequence at tail")
    {
        sequence.Replace(4, 0, 'R', 1, 47);
        REQUIRE(sequence.Length() == 5);
        REQUIRE(elemcmp(sequence.cbegin(), "ACTGR", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x3A\x3D\x30\x1", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,17,29,42,47 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Insert subsequence into empty sequence")
    {
        empty_sequence.Replace(0, 0, sequence);
        REQUIRE(sequence.Length() == 4);
        REQUIRE(elemcmp(sequence.cbegin(), "ACTG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x3A\x3D\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,17,29,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Insert empty sequence")
    {
        sequence.Replace(1, 0, empty_sequence);
        REQUIRE(sequence.Length() == 4);
        REQUIRE(elemcmp(sequence.cbegin(), "ACTG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x3A\x3D\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,17,29,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Delete subsequence from middle")
    {
        sequence.DeleteSubsequence(2, 1);
        REQUIRE(sequence.Length() == 3);
        REQUIRE(elemcmp(sequence.cbegin(), "ACG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3F\x3A\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 5,17,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Delete subsequence at head")
    {
        sequence.DeleteSubsequence(0, 2);
        REQUIRE(sequence.Length() == 2);
        REQUIRE(elemcmp(sequence.cbegin(), "TG", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x3D\x30", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 29,42 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
    SECTION("Delete subsequence at tail")
    {
        sequence.DeleteSubsequence(1, 3);
        REQUIRE(sequence.Length() == 1);
        REQUIRE(*sequence.cbegin() == 'A');
        REQUIRE(*sequence.QualityBegin()== '\x3F');
        REQUIRE(*sequence.PeakBegin() == 5);
    }
    SECTION("Reverse+complement")
    {
        sequence.ReverseComplement(54);
        REQUIRE(sequence.Length() == 4);
        REQUIRE(elemcmp(sequence.cbegin(), "CAGT", sequence.Length()) == 0);
        REQUIRE(elemcmp(sequence.QualityBegin(), "\x30\x3D\x3A\x3F", sequence.Length()) == 0);
        int32_t peaks_cmp[] = { 11,24,36,48 };
        REQUIRE(elemcmp(sequence.PeakBegin(), peaks_cmp, sequence.Length()) == 0);
    }
}

TEST_CASE("ab1 constructor", "[ab1_load]")
{
    Ab1File ab1file("test.ab1");
    Ab1File::Iterator<char> seq_begin, seq_end;
    Ab1File::Iterator<uint8_t> qual_begin, qual_end;
    Ab1File::Iterator<int32_t> peak_begin, peak_end;
    Ab1File::Iterator<int32_t> trace_begin[4], trace_end[4];
    REQUIRE(ab1file.SearchTag("PBAS", 1, seq_begin, seq_end) == Ab1File::SearchResult::SUCCESS);
    REQUIRE(ab1file.SearchTag("PCON", 1, qual_begin, qual_end) == Ab1File::SearchResult::SUCCESS);
    REQUIRE(ab1file.SearchTag("PLOC", 1, peak_begin, peak_end) == Ab1File::SearchResult::SUCCESS);
    for (int32_t i = 0; i < 4; ++i) {
        REQUIRE(ab1file.SearchTag("DATA", 9 + i, trace_begin[i], trace_end[i]) == Ab1File::SearchResult::SUCCESS);
    }
    NucleotideSequence sequence(seq_begin, seq_end, qual_begin, qual_end, "test");
    sequence.LoadTraces(peak_begin, peak_end, trace_begin, trace_end);
    REQUIRE(sequence.Length() != 0);
    REQUIRE(sequence.HasTraces());
    time_t datetime;
    REQUIRE(ab1file.DateTime("RUND", "RUNT", 1, datetime) == Ab1File::SearchResult::SUCCESS);
    REQUIRE(ab1file.DateTime("RUND", "RUNT", 2, datetime) == Ab1File::SearchResult::SUCCESS);
}

TEST_CASE("scf constructor", "[scf_load]")
{
    ScfFile scffile("test.scf");
    ScfFile::TraceIterator<int32_t> trace_begin[4], trace_end[4];
    for (int32_t i = 0; i < 4; ++i) {
        REQUIRE(scffile.Traces("ACGT"[i], trace_begin[i], trace_end[i]));
    }
    int32_t buffer[16384];
    for (int32_t i = 0; i < 4; ++i) {
        ScfFile::TraceIterator<int32_t> end = trace_end[i];
        size_t j = 0;
        for (ScfFile::TraceIterator<int32_t> it = trace_begin[i]; it != end; ++it)
            buffer[j++] = *it;
    }
}

TEST_CASE("search by alignment", "[align_search]")
{
    NucleotideSequence sequence("ACGATCAGACTGCGAAGATTCCATACAGCG");
    REQUIRE(sequence.SearchByAlignmentFwd(0, "CAGACAGCG", 80) == 5);
    REQUIRE(sequence.SearchByAlignmentBack(sequence.Length() - 1, "CAGACAGCG", 80) == 21);
}