#include <sstream>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/search/fm_index/bi_fm_index.hpp>
#include <seqan3/search/algorithm/search_ng2.hpp>

int main(int argc, char ** argv)
{
    using seqan3::operator""_dna4;
    auto genome = std::vector<std::vector<seqan3::dna4>>{"GAATTAACGAAC"_dna4, "GATTCTCTCTA"_dna4};
    auto index  = seqan3::bi_fm_index{genome};

    index.init();

	// searching for GAAT
	std::cout << "\nsearch GAAT step by step\n";
	{
		auto cursor = index.cursor();
		std::cout << "count: " << cursor.count() << "\n"; // count: 25 (23 pos + 2 sentinels)

		cursor = cursor.extend_right(cursor.convert('G'_dna4));
		std::cout << "count: " << cursor.count() << "\n"; // count: 3

		cursor = cursor.extend_right(cursor.convert('A'_dna4));
		std::cout << "count: " << cursor.count() << "\n"; // count: 3

		cursor = cursor.extend_right(cursor.convert('A'_dna4));
		std::cout << "count: " << cursor.count() << "\n"; // count: 2

		cursor = cursor.extend_right(cursor.convert('T'_dna4));
		std::cout << "count: " << cursor.count() << "\n"; // count: 1
	}

	// searching for GAAT in a loop
	std::cout << "\nsearch GAAT in a loop\n";
	{
		auto cursor = index.cursor();
		for (auto v : "GAAT"_dna4) {
			cursor = cursor.extend_right(cursor.convert(v));
		}
		std::cout << "count: " << cursor.count() << "\n"; // count: 3, 3, 2, 1
	}

	// fast searching for $1, A, C, G, T and $2
	// Notice "$1" is the 0 sentinel from sdsl
	//        "$2" is the sentinel we add for collections
	std::cout << "\nfast search for $1, A, C, G, T and $2\n";
	{
		auto cursor = index.cursor();
		cursor.extend_right_cb([](auto c, auto newCursor) {
			if (c == 0) {
				std::cout << "$1 ";
			} else if (c == seqan3::alphabet_size<seqan3::dna4>+1) {
				std::cout << "$2 ";
			} else {
				seqan3::debug_stream << " " << seqan3::dna4{}.assign_rank(c-1) << " ";
			}
			std::cout << newCursor.count() << "\n";
		});
	}







    return 0;
}
