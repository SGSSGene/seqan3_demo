#include <sstream>
#include <string>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/output.hpp>

int main()
{
    using seqan3::operator""_dna5;

    seqan3::sequence_file_output fout{std::ostringstream{}, seqan3::format_fasta{}};

    for(int i = 0; i < 5; ++i) // some criteria
    {
        std::string id{"test_id"};
        seqan3::dna5_vector seq{"ACGT"_dna5};

        // ...

        fout.emplace_back(seq, id);
    }
}
