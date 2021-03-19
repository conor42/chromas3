// System class containing platform-specific methods.
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

#include <Shlobj.h>
#include "system.h"

const char System::data_dir_name[] = "libchromas";

static size_t WideCharToMultiByte(const wchar_t* src, std::string dst)
{
	int src_len = (int)wcslen(src);
	int dst_len = WideCharToMultiByte(CP_UTF8, 0, src, src_len, nullptr, 0, nullptr, nullptr);
	if (dst_len) {
		dst.resize(dst_len);
		// This would break on non-contiguous std::string storage, but Microsoft's are contiguous
		// and the standard is changing to allow it anyway.
		WideCharToMultiByte(CP_UTF8, 0, src, src_len, &dst[0], (int)dst.size(), nullptr, nullptr);
	}
	else {
		dst.clear();
	}
	return dst_len;
}

std::string System::ProgramDataDir()
{
    PWSTR path = nullptr;
    std::string returned_path;
    if (SHGetKnownFolderPath(FOLDERID_RoamingAppData, 0, nullptr, &path) == S_OK) {
        WideCharToMultiByte(path, returned_path);
        AppendName(returned_path, data_dir_name);
    }
    CoTaskMemFree(path);
    return std::move(returned_path);
}

void System::AppendName(std::string& path, const char* name)
{
    if (path.back() != '\\' && path.back() != '/')
        path.push_back('\\');
    path.append(name);
}
