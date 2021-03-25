// Tests for codon translation.
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
#include "geneticcodes.h"

// TCAG
static const char standard_code[] = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";
static const char can_start[] =     "---M---------------M---------------M----------------------------";

TEST_CASE("Codon translation", "[translation]")
{
    static char seq[64 * 3 + 1] = "";
    static const char bases[] = "TCAG";
    for (size_t i = 0; i < 64; ++i) {
        size_t pos = i * 3;
        seq[pos] = bases[(i >> 4) & 3];
        seq[pos + 1] = bases[(i >> 2) & 3];
        seq[pos + 2] = bases[i & 3];
    }

    GeneticCodes::Codon codon;
    for (size_t i = 0; i < 64; ++i) {
        codon = GeneticCodes::TranslateForward(seq, i * 3, sizeof(seq), 0);
        REQUIRE(codon.amino_acid == standard_code[i]);
        REQUIRE(codon.can_start == (can_start[i] == 'M'));
    }
	codon = GeneticCodes::TranslateForward("TCN", 0, 3, 0);
	REQUIRE(codon.amino_acid == 'S');
	REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("CTN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'L');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("CCN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'P');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("CGN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'R');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("ACN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'T');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("GTN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'V');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("GCN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'A');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("GGN", 0, 3, 0);
    REQUIRE(codon.amino_acid == 'G');
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("RAY", 0, 3, 0); // (A/G)A(T/C)
    REQUIRE(codon.amino_acid == 'B'); // D or N
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("SAR", 0, 3, 0); // (G/C)A(A/G)
    REQUIRE(codon.amino_acid == 'Z'); // E or Q
    REQUIRE(codon.can_start == false);
    codon = GeneticCodes::TranslateForward("HTG", 0, 3, 0); // (T/C/A)TG
    REQUIRE(codon.amino_acid == 'X');
    REQUIRE(codon.can_start == true);

    char leucine[] = "CTN";
    static const char redundant[] = "BDHKMRSVWY";
    for (const char* p = redundant; *p; ++p) {
		leucine[2] = *p;
		codon = GeneticCodes::TranslateForward(leucine, 0, 3, 0);
		REQUIRE(codon.amino_acid == 'L');
		REQUIRE(codon.can_start == false);
    }
}