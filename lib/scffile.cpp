// SCF file reader
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
#include "exception.h"
#include "scffile.h"

// Avoid large memory consumption where a very large file (of the wrong type) is specified.
static const size_t MAX_FILE_SIZE = size_t(1) << 25;

const char ScfFile::base_order[ScfFile::TRACE_COUNT] = { 'A','C','G','T' };

static inline const void ThrowCorruptHeader()
{
    throw invalid_file_format("SCF file header is corrupted.");
}

ScfFile::ScfFile(const char* path)
{
    std::ifstream stream;
    stream.exceptions(std::ios_base::failbit | std::ios_base::badbit);
    stream.open(path, std::ios_base::in | std::ios_base::binary | std::ios_base::ate);

    file_size_ = stream.tellg();
    if (file_size_ > MAX_FILE_SIZE)
        throw std::invalid_argument(path);
    stream.seekg(0, stream.beg);

    file_buffer_ = std::make_unique<char[]>(file_size_);
    stream.read(file_buffer_.get(), file_size_);
    stream.close();

    ValidateHeader();
}

void ScfFile::ValidateHeader()
{
    header_ = reinterpret_cast<const ScfHeader*>(file_buffer_.get());
    if (header_->version[0] > '3')
        throw unsupported_file_format("SCF file versions above 3.x are not supported.");

    size_t sample_size = header_->SampleSize();

    if (header_->samples_offset > file_size_ || header_->samples_offset + header_->samples * sample_size * TRACE_COUNT > file_size_)
        ThrowCorruptHeader();

    size_t base_size = (header_->version[0] < '3') ? sizeof(ScfBase) : BASE_SIZE_V3;
    if(header_->bases_offset > file_size_ || header_->bases_offset + header_->bases * base_size > file_size_)
        ThrowCorruptHeader();

    if (header_->comments_offset > file_size_ || header_->comments_offset + header_->comments_size > file_size_)
        ThrowCorruptHeader();
}
