// ab1 file reader
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
#include "ab1file.h"
#include "endian.h"
#include "exception.h"
#include "log.h"
#include "strutil.h"

// Avoid large memory consumption where a very large file (of the wrong type) is specified.
static const size_t MAX_FILE_SIZE = size_t(1) << 26;

static const char ab1_magic[] = "ABIF";

static inline const void ThrowCorruptDirectory()
{
    throw invalid_file_format("ab1 file directory is corrupted.");
}

Ab1File::Ab1File(const char* path)
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

    FindDirectory();
}

Ab1File::SearchResult Ab1File::SearchTag(const char* tag, int32_t number, DirEntry& converted)
{
    DirEntry* entry = SearchTag(tag, number);
    if (!entry) {
        Log(2) << "ab1 tag not found: " << tag << " number " << number << endl;
        return SearchResult::NOT_FOUND;
    }
    converted.element_len = system_endian(entry->element_len);
    converted.elements = system_endian(entry->elements);

    size_t data_pos = DataPos(entry);
    size_t data_length = (size_t)converted.element_len * converted.elements;
    if (data_pos > file_size_ || data_pos + data_length > file_size_
        || data_pos + data_length < data_pos || entry->bytes < data_length)
        ThrowCorruptDirectory();

    converted.data_type = (Type)system_endian((unsigned short)entry->data_type);
    converted.bytes = system_endian(entry->bytes);
    converted.data = (uint32_t)data_pos;

    Log(2) << "Found ab1 tag " << tag << " number " << number << ", data " << data_pos << ", elements " << converted.elements
        << ", size " << converted.element_len << ", type " << (int)converted.data_type << endl;

    return SearchResult::SUCCESS;
}

Ab1File::SearchResult Ab1File::SearchTag(const char* tag, int32_t number, std::string& result)
{
    result.clear();

    DirEntry converted;
    SearchResult found = SearchTag(tag, number, converted);
    if (found != SearchResult::SUCCESS)
        return found;
    if (converted.element_len != 1)
        return SearchResult::INCOMPATIBLE_TYPE;

    if (converted.data_type == PSTRING) {
        if(!converted.elements)
            ThrowCorruptDirectory();
        result.reserve(converted.elements - 1);
        for (size_t i = 1; i < converted.elements; ++i)
            result.push_back(file_buffer_[converted.data + i]);
    }
    else if (converted.data_type == CSTRING) {
        result.reserve(converted.elements);
        for (size_t i = 0; i < converted.elements; ++i)
            result.push_back(file_buffer_[converted.data + i]);
    }
    else {
        return SearchResult::INCOMPATIBLE_TYPE;
    }
    return SearchResult::SUCCESS;
}

Ab1File::SearchResult Ab1File::SearchTag(const char* tag, int32_t number, int32_t& value)
{
    Iterator<int32_t> begin, end;
    SearchResult found = SearchTag(tag, number, begin, end);
    if (found != SearchResult::SUCCESS)
        return found;
    value = *begin;
    return SearchResult::SUCCESS;
}

Ab1File::SearchResult Ab1File::DateTime(const char* date_tag, const char* time_tag, int32_t number, std::time_t& datetime)
{
    Iterator<Date> date_begin, date_end;
    Iterator<Time> time_begin, time_end;
    SearchResult found;
    if ((found = SearchTag(date_tag, number, date_begin, date_end)) != SearchResult::SUCCESS)
        return found;

    std::tm tm;
    tm.tm_year = system_endian(date_begin->year) - 1900;
    tm.tm_mon = date_begin->month;
    tm.tm_mon -= (tm.tm_mon != 0);
    tm.tm_mday = date_begin->day;

    if ((found = SearchTag(time_tag, number, time_begin, time_end)) != SearchResult::SUCCESS) {
        Log(2) << "Date found but time not found, result = " << (int)found << endl;
        tm.tm_hour = 0;
        tm.tm_min = 0;
        tm.tm_sec = 0;
    }
    else {
        tm.tm_hour = time_begin->hour;
        tm.tm_min = time_begin->minute;
        tm.tm_sec = time_begin->second;
    }
    datetime = std::mktime(&tm);

    return SearchResult::SUCCESS;
}

void Ab1File::FindDirectory()
{
    // 1990s legacy files often have a short length of prepended data as a result of file transfer from Mac OS.
    ptrdiff_t start = FindSubstring(file_buffer_.get(), std::min<size_t>(file_size_, 1024), ab1_magic);
    if (start < 0)
        throw invalid_file_format("ab1 file header not found.");
    file_start_ = start;

    Header header(file_buffer_.get(), file_start_);

    Log(1) << "ab1 header found at offset " << file_start_  << ", file version " << system_endian(header.version) << endl;
    if(system_endian(header.dir_pos.element_len) != sizeof(DirEntry))
        throw invalid_file_format("unsupported ab1 directory entry structure.");

    size_t dir_pos = system_endian(header.dir_pos.data);
    size_t dir_end = dir_pos + system_endian(header.dir_pos.elements) * sizeof(DirEntry);
    Log(1) << "ab1 directory located at offset " << dir_pos << " with " << ((dir_end - dir_pos) / sizeof(DirEntry)) << " entries." << endl;
    if (dir_pos > (file_size_ - file_start_) || dir_end > (file_size_ - file_start_))
        throw invalid_file_format("ab1 file header is corrupted.");

    dir_begin_ = reinterpret_cast<DirEntry*>(file_buffer_.get() + file_start_ + dir_pos);
    dir_src_ = dir_begin_;
    dir_end_ = reinterpret_cast<DirEntry*>(file_buffer_.get() + file_start_ + dir_end);
}

Ab1File::DirEntry* Ab1File::SearchTag(const char* tag, int32_t number)
{
    if (dir_end_ == dir_begin_)
        return nullptr;

    DirEntry* end = dir_src_;
    // Using a persistent source pointer allows faster search of sorted tags if queries are also sorted.
    do {
        ++dir_src_;
        if (dir_src_ == dir_end_)
            dir_src_ = dir_begin_;
        if (memcmp(tag, &dir_src_->name, sizeof(dir_src_->name)) == 0 && system_endian((uint32_t)dir_src_->number) == number)
            return dir_src_;
    } while (dir_src_ != end);

    return nullptr;
}

size_t Ab1File::DataPos(DirEntry* entry)
{
    if (system_endian(entry->bytes) <= sizeof(entry->data))
        return reinterpret_cast<const char*>(&entry->data) - file_buffer_.get();
    return system_endian(entry->data) + file_start_;
}
