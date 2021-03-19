// Simple string utilities.
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

#include <cstring>

inline char* SkipSpaces(const char* p)
{
	while (*p == ' ')
		++p;
	return const_cast<char*>(p);
}

inline char* SkipSpacesBackward(const char* begin, const char* p)
{
	while (p > begin && p[-1] == ' ')
		--p;
	return const_cast<char*>(p);
}

inline ptrdiff_t FindSubstring(const char* data, size_t length, const char* query)
{
	size_t query_len = strlen(query);
	if (!query_len || length < query_len)
		return -1;

	const char* end = data + length - query_len + 1;
	for (const char* p = data; p < end; ++p) {
		p = reinterpret_cast<const char*>(memchr(p, query[0], end - p));
		if (!p)
			break;
		if (memcmp(p + 1, query + 1, query_len - 1) == 0)
			return p - data;
	}
	return -1;
}