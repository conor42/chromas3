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

#include <fstream>
#include <iostream>
#include <string>
#include "geneticcodes.h"
#include "strutil.h"
#include "system.h"

GeneticCodes::static_initializer GeneticCodes::impl_;
static const char standard_code[] = "Standard; FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG ---M---------------M---------------M";

static const size_t CODON_COUNT = 64;

static size_t CodonIndexToFlags(size_t index)
{
	size_t flags = 0;
	for (size_t i = 0; i < 6; i += 2) {
		flags |= size_t(1) << (((index >> i) & 3) + i * 2);
	}
	return flags;
}

static bool CodonIndexMatch(size_t index, size_t flags)
{
	size_t i = CodonIndexToFlags(index);
	return (i & flags) == i;
}

GeneticCodes::CodonTable::CodonTable(const char* line)
{
	const char* p = strchr(line, ';');
	if (!p)
		return;
	p = SkipSpacesBackward(line, p);
	ptrdiff_t i = p - line;
	if (!i)
		return;
	name.assign(line, i);
	p = SkipSpaces(p + 1);
	std::array<Codon, CODON_COUNT> codons;
	for (i = 0; i < CODON_COUNT; ++i) {
		char aa = LookupTables::Uppercase(p[i]);
		if (aa != '*' && (aa < 'A' || aa > 'Z'))
			return;
		codons[i].amino_acid = aa;
	}
	p = strchr(p, ' ');
	if (p) {
		p = SkipSpaces(p);
		for (i = 0; i < 64; ++i) {
			char c = LookupTables::Uppercase(p[i]);
			if (c != '-' && c != 'M')
				break;
			codons[i].can_start = (c == 'M');
		}
		for (; i < 64; ++i)
			codons[i].can_start = false;
	}
	for (size_t i = 0; i < TABLE_SIZE; ++i) {
		unsigned int aa_flags = 0;
		size_t aa_index = 0;
		for (size_t j = 0; j < 64; ++j) {
			if (CodonIndexMatch(j, i)) {
				aa_flags |= 1 << LookupTables::CharIndex(codons[j].amino_acid);
				aa_index = j;
			}
		}
		if ((aa_flags & (aa_flags - 1)) == 0) {
			table[i] = codons[aa_index];
		}
		else if (aa_flags == ((1 << LookupTables::CharIndex('D')) | (1 << LookupTables::CharIndex('N')))) {
			table[i].amino_acid = AA_D_OR_N;
			table[i].can_start = false;
		}
		else if (aa_flags == ((1 << LookupTables::CharIndex('E')) | (1 << LookupTables::CharIndex('Q')))) {
			table[i].amino_acid = AA_E_OR_Q;
			table[i].can_start = false;
		}
		else {
			table[i].amino_acid = AA_UNKNOWN;
			table[i].can_start = false;
		}
	}
}

int GeneticCodes::Search(const char* name)
{
    for (size_t i = 0; i < impl_.codon_tables.size(); ++i)
        if (impl_.codon_tables[i].name.compare(name) == 0)
            return (int)i;
    return -1;
}

void GeneticCodes::Initialize()
{
    std::string path = System::ProgramDataDir();
    System::AppendName(path, "geneticcodes");

    std::ifstream stream(path);
    if (stream.fail()) {
        std::cerr << "Genetic code file not found; loading standard code only." << std::endl;
		impl_.codon_tables.emplace_back(standard_code);
        return;
    }

    do {
        std::string line;
        std::getline(stream, line);
        if (stream.fail()) {
            std::cerr << "Error reading genetic code file; only " << impl_.codon_tables.size() << " table(s) loaded." << std::endl;
            break;
        }
		impl_.codon_tables.emplace_back(line.c_str());
    } while (!stream.eof());

    if (impl_.codon_tables.empty())
		impl_.codon_tables.emplace_back(standard_code);
}
