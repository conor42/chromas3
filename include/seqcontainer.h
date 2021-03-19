// Chromas sequence container class.
//
// Copyright 2021 Conor N. McCarthy
//
// Used by the sequence class to store peaks and quality data in addition to nucleotide sequences.
// A reserve at the end allows for some expansion due to editing.
// No bounds checking is done except for assertions for use with debugging. Sanitisation of function parameters
// must be done by the calling app.
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

#include <memory>
#include <cassert>

template<typename value_type>
class SeqContainer
{
public:
	SeqContainer()
		: length_(0), max_length_(0)
	{
	}

	SeqContainer(const value_type* sequence, size_t length)
		: length_(0), max_length_(0)
	{
		Reallocate(length);
		std::copy_n(sequence, length, sequence_.get());
		length_ = length;
	}

	template<class InputIterator>
	SeqContainer(InputIterator begin, InputIterator end)
		: length_(0), max_length_(0)
	{
		Reallocate(end - begin);

		for (InputIterator it = begin; it != end; ++it)
			sequence_[length_++] = *it;
		assert(length_ == end - begin);
	}

	SeqContainer(const SeqContainer& source) = default;

	SeqContainer(SeqContainer&& source) noexcept {
		return operator=(std::move(source));
	}

	SeqContainer& operator=(const SeqContainer& rval) = default;

	SeqContainer& operator=(SeqContainer&& rval) noexcept {
		sequence_ = std::move(rval.sequence_);
		length_ = rval.length_;
		max_length_ = rval.max_length_;
		// Source is left indeterminate but valid.
		// unique_ptr members can take care of themselves, but lengths should be set to zero.
		rval.length_ = 0;
		rval.max_length_ = 0;
		return *this;
	}

	~SeqContainer() = default;

	// Resize the subsequence of length old_length which begins at start_pos.
	// It will grow or shrink to new_length. Extra elements will be left undefined, or excess elements
	// will be deleted from the end of the subsequence.
	void ResizeSubsequence(size_t start_pos, size_t old_length, size_t new_length) {
		assert(start_pos <= length_ && start_pos + old_length <= length_);

		if (new_length == old_length)
			return;

		if (new_length > old_length) {
			Reallocate(length_ + new_length - old_length);
			std::move_backward(sequence_.get() + start_pos + old_length, sequence_.get() + length_, sequence_.get() + length_ + new_length - old_length);
		}
		else {
			std::move(sequence_.get() + start_pos + old_length, sequence_.get() + length_, sequence_.get() + start_pos + new_length);
		}
		length_ += new_length - old_length;
	}

	// Replace the subsequence of length old_length, which begins at start_pos, with source.
	// The subsequence will be resized to fit source exactly.
	// This is typically used for undo/redo and copy/paste.
	void Replace(size_t start_pos, size_t old_length, const SeqContainer& source) {
		ResizeSubsequence(start_pos, old_length, source.length_);
		OverwriteSubsequence(start_pos, source.cbegin(), source.length_);
	}


	// Replace the subsequence of length old_length, which begins at start_pos, with one new element.
	void Replace(size_t start_pos, size_t old_length, value_type value) {
		ResizeSubsequence(start_pos, old_length, 1);
		OverwriteSubsequence(start_pos, &value, 1);
	}

	// Delete the subsequence of length old_length, which begins at start_pos.
	void DeleteSubsequence(size_t start_pos, size_t length) {
		ResizeSubsequence(start_pos, length, 0);
	}

	// Truncate the sequence to new_length.
	void Truncate(size_t new_length) {
		ResizeSubsequence(new_length, length_ - new_length, 0);
	}

	// Fill the subsequence of length old_length, which begins at start_pos, with a value.
	void FillSubsequence(size_t start_pos, size_t length, value_type value) {
		assert(start_pos <= length_ && start_pos + length <= length_);
		std::fill_n(sequence_.get() + start_pos, length, value);
	}

	const value_type* cbegin() const {
		return sequence_.get();
	}

	const value_type* cend() const {
		return sequence_.get() + length_;
	}

	size_t Length() const {
		return length_;
	}

	operator const value_type*() const {
		return cbegin();
	}

	const value_type& operator[](size_t pos) const {
		return sequence_[pos];
	}

	value_type& operator[](size_t pos) {
		return sequence_[pos];
	}

	operator bool() const {
		return sequence_.operator bool();
	}

private:

	size_t ReserveSize(size_t size) {
		return (size >> 5) + 1;
	}

	void Reallocate(size_t new_length) {
		if (max_length_ < new_length) {
			// The reserve at the end allows for some editing, but if this is exceeded, copy the sequence to a new buffer.
			size_t new_max_length = new_length + ReserveSize(new_length);
			auto new_buffer = std::make_unique<value_type[]>(new_max_length);
			std::copy_n(sequence_.get(), length_, new_buffer.get());
			sequence_ = std::move(new_buffer);
			max_length_ = new_max_length;
		}
	}

	void OverwriteSubsequence(size_t start_pos, const value_type* source, size_t length) {
		assert(start_pos <= length_ && start_pos + length <= length_);
		std::copy_n(source, length, sequence_.get() + start_pos);
	}

	std::unique_ptr<value_type[]> sequence_;
	size_t length_;
	size_t max_length_;
};

