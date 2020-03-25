// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \brief Provides the public interface for search algorithms.
 */

#pragma once

#include <array>

#include <sdsl/suffix_trees.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/slice.hpp>
#include <seqan3/std/ranges>

namespace seqan3 {

//*!TODO quick ugly reimplentation of bi_fm_index_cursor

template <typename index_t>
class bi_fm_index_cursor_ng2
{
public:
    using size_type     = typename index_t::size_type;
    using alphabet_type = typename index_t::alphabet_type;
    using rank_type     = std::decay_t<decltype(std::declval<alphabet_type>().to_rank())>;

private:
public: //!TODO

    index_t const* index {};

    size_type fwd_lb {};
    size_type rev_lb {};
    size_type length {}; // should probably be called something like "count". `fwd_rb = fwd_lb + length-1` and `rev_rb = rev_lb + length -1`
    size_type depth  {}; // If locate is adjusted, this variable is not needed.

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Default constructor. Accessing member functions on a default constructed object is undefined behavior.
    //        Default construction is necessary to make this class semi-regular and e.g., to allow construction of
    //        std::array of cursors.
    bi_fm_index_cursor_ng2() noexcept                                           = default;
    bi_fm_index_cursor_ng2(bi_fm_index_cursor_ng2 const &) noexcept             = default;
    bi_fm_index_cursor_ng2 & operator=(bi_fm_index_cursor_ng2 const &) noexcept = default;
    bi_fm_index_cursor_ng2(bi_fm_index_cursor_ng2 &&) noexcept                  = default;
    bi_fm_index_cursor_ng2 & operator=(bi_fm_index_cursor_ng2 &&) noexcept      = default;
    ~bi_fm_index_cursor_ng2()                                                   = default;

    //! \brief Construct from given index.
    bi_fm_index_cursor_ng2(index_t const & _index) noexcept
      : index  {&_index}
      , length {index->size()}
    {}

protected:
public://!TODO
    bi_fm_index_cursor_ng2(index_t const & _index, size_type _fwd_lb, size_type _rev_lb, size_type _length, size_type _depth) noexcept
      : index  {&_index}
      , fwd_lb {_fwd_lb}
      , rev_lb {_rev_lb}
      , length {_length}
      , depth  {_depth}
    {}

public:

	/* Extends cursor by on character on the right
	 */
    auto extend_right(rank_type const c) const noexcept
    {
        assert(index != nullptr);
        auto&      csa            = index->fwd_fm.index;
        auto const c_begin        = csa.C[c];
        auto const [rank_l, s, b] = csa.wavelet_tree.lex_count(fwd_lb, fwd_lb+length, c);
        return bi_fm_index_cursor_ng2{*index, c_begin + rank_l, rev_lb + s, length -b -s, depth+1};
    }
    auto extend_left(rank_type const c) const noexcept
    {
        assert(index != nullptr);
        auto&      csa            = index->rev_fm.index;
        auto const c_begin        = csa.C[c];
        auto const [rank_l, s, b] = csa.wavelet_tree.lex_count(rev_lb, rev_lb+length, c);
        return bi_fm_index_cursor_ng2{*index, fwd_lb + s, c_begin + rank_l, length -b -s, depth+1};
    }

	/* Loops through all possible character extensions on the right
	 */
    template <typename CB>
    void extend_right_cb(CB const& cb) const noexcept {
        assert(index != nullptr);
		auto& csa = index->fwd_fm.index;
		csa.wavelet_tree.lex_count_fast_cb(fwd_lb, fwd_lb+length, [&](size_t rank_l, size_t s, size_t b, size_t, size_type c) noexcept {
	        size_type const c_begin = csa.C[c];
	        cb(c, bi_fm_index_cursor_ng2{*index, c_begin + rank_l, rev_lb + s, length -b -s, depth+1});
		});
    }

    template <typename CB>
    void extend_left_cb(CB const& cb) const noexcept
    {
        assert(index != nullptr);
        auto& csa = index->rev_fm.index;
		csa.wavelet_tree.lex_count_fast_cb(rev_lb, rev_lb+length, [&](size_t rank_l, size_t s, size_t b, size_t, size_type c) noexcept {
	        size_type const c_begin = csa.C[c];
	        cb(c, bi_fm_index_cursor_ng2{*index, fwd_lb + s, c_begin + rank_l, length -b -s, depth+1});
		});
    }

    /*!\brief Counts the number of occurrences of the searched query in the text.
     * \returns Number of occurrences of the searched query in the text.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exceptions
     *
     * No-throw guarantee.
     */
    size_type count() const noexcept
    {
        assert(index != nullptr);
        return length;
    }


    /*!TODO this should not be based on fwd_lb but on rev_lb, that allows us to drop `depth` variable
     */
    template <typename delegate_t>
    void locate(delegate_t delegate) const
    {
        assert(index != nullptr);
        auto os = index->size() - depth - 1;

        for (size_type i = 0; i < count(); ++i)
        {
            //if (os < index->fwd_fm.index[fwd_lb + i]) continue;
            size_type loc               = os - index->fwd_fm.index[fwd_lb + i];
            size_type sequence_rank     = index->fwd_fm.text_begin_rs.rank(loc + 1);
            size_type sequence_position = loc - index->fwd_fm.text_begin_ss.select(sequence_rank);
            delegate(sequence_rank - 1, sequence_position);
        }
    }

    /** !TODO Don't really think this should be part of the cursor
     */
    template <typename query_alphabet_type>
    auto convert(query_alphabet_type const query_c) const -> size_type {
        static_assert(std::convertible_to<query_alphabet_type, alphabet_type>,
                     "The character must be convertible to the alphabet of the index.");

        auto index_c = to_rank(static_cast<alphabet_type>(query_c)) + 1;
        if constexpr(!std::same_as<alphabet_type, sdsl::plain_byte_alphabet>)
        {
            assert(index->fwd_fm.index.char2comp[index_c] == index->rev_fm.index.char2comp[index_c]);
            assert(not (index_c == 0 && index->fwd_fm.index.char2comp[index_c] > 0));
            return index->fwd_fm.index.char2comp[index_c];
        } else {
            return index_c;
        }
    }

};
}

