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

#include <seqan3/search/fm_index/bi_fm_index_cursor_ng2.hpp>

namespace seqan3
{

template <typename cursor_t, typename query_t, typename search_scheme_t, typename delegate_t>
struct Search_ng2 {
	query_t const&          query;
	search_scheme_t const&  search;
	size_t                  qidx;
	delegate_t const&       delegate;

	using index_alphabet_type = typename cursor_t::alphabet_type;
	using query_alphabet_type = innermost_value_type_t<query_t>;
	using index_rank_type     = std::decay_t<decltype(std::declval<index_alphabet_type>().to_rank())>;

	constexpr static size_t index_sigma = alphabet_size<index_alphabet_type>;

	Search_ng2(cursor_t const& _cursor, query_t const& _query, search_scheme_t const& _search, size_t _qidx, delegate_t const& _delegate)
		: query    {_query}
		, search   {_search}
		, qidx     {_qidx}
		, delegate {_delegate}
	{
		assert(search.pi.size() > 1);
		assert(search.pi.size() == query.size());
		search_next(_cursor, 0, 0);
	}

	bool goToRight(int pos) const noexcept {
		assert(pos >= 0 and pos < search.pi.size());
		return (pos == 0) ? (search.pi[0] < search.pi[1]) : (search.pi[pos-1] < search.pi[pos]);
	}

	auto extend(cursor_t const& cur, int pos, index_rank_type rank) const noexcept {
		if (goToRight(pos)) {
			return cur.extend_right(rank);
		}
		return cur.extend_left(rank);
	}

	template <typename CB>
	void extend_cb(cursor_t const& cur, int pos, CB const& cb) const noexcept {
		if (goToRight(pos) > 0) {
			cur.extend_right_cb(cb);
		} else {
			cur.extend_left_cb(cb);
		}
	}


	// single search step
	void search_next(cursor_t const& cur, int e, size_t pos) const noexcept {
		if (cur.count() == 0) {
			return;
		}

		if (pos == query.size()) {
			delegate(qidx, cur);
			return;
		}
		assert(pos >= 0);
		assert(pos < search.l.size());
		assert(pos < search.u.size());
		assert(pos < search.pi.size());
		assert(search.pi[pos] >= 0);
		assert(search.pi[pos] < query.size());

		auto rank_c = cur.convert(query[search.pi[pos]]);

		std::array<cursor_t, index_sigma+1> cursors;
		if (search.l[pos] <= e+1 and e+1 <= search.u[pos]) {
			extend_cb(cur, pos, [&](auto c, auto newCur) {
				if (c > index_sigma) return;
				cursors[c] = newCur;
			});

			// search for substitute
			for (index_rank_type i{1}; i <= index_sigma; ++i) {
				if (i == rank_c) continue;
				search_next(cursors[i], e+1, pos+1);
			}

			// search for hit (when substitutes are allowed, using performance increase)
			if (search.l[pos] <= e and e <= search.u[pos]) {
				search_next(cursors[rank_c], e, pos+1);
			}

		// search for hit (no substitutes are allowed)
		} else if (search.l[pos] <= e and e <= search.u[pos]) {
			auto newCur = extend(cur, pos, rank_c);
			search_next(newCur, e, pos+1);
		}
	}
};

template <typename index_t, typename queries_t, typename search_schemes_t, typename delegate_t>
void search_ng2(index_t const & index, queries_t && queries, search_schemes_t const & search_scheme, delegate_t && delegate)
{
	auto internal_delegate = [&delegate] (size_t qidx, auto const & it) {
		it.locate([&](auto p1, auto p2) {
			delegate(qidx, p1, p2);
		});
	};

	using query_alphabet_t = innermost_value_type_t<queries_t>;
	auto rootCursor = bi_fm_index_cursor_ng2<index_t>{index};
	for (size_t qidx{0}; qidx < queries.size(); ++qidx) {
		auto const& query = queries[qidx];
		for (size_t j{0}; j < search_scheme.size(); ++j) {
			auto const& search = search_scheme[j];
			Search_ng2{rootCursor, query, search, qidx, internal_delegate};
		}
	}
}

}
