# Demo to show possible improvements to sdsl and seqan3

Instructions:
1. clone this repository
2. run `mkdir seqan3_demo/build; cd seqan3_demo/build`
3. run cmake `cmake ..`
4. enjoy `./cursor_demo`


This repo does not use submodules. Everything is included via subrepos.

Check commit '84801d139' to see the additions to sdsl.

Check file 'seqan3/search/fm_index/bi_fm_index_cursor_ng2.hpp' for cursor implementation.
Check file 'src/main.cpp' for usage of cursor.

Also file 'seqan3/search/search_algorithm/search_ng2.hpp' has a sketch of how a search could look like.
This search will not work with "normal" search schemes, but need so called expanded ones. (Code to do this is currently missing)

