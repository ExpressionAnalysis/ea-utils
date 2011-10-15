///A tutorial about local alignments.
#include <iostream>
#include <seqan/align.h>

using namespace seqan;

int main(int argc, char **argv)
{
///Example 1: This program applies the Smith-Waterman algorithm to compute the best local alignment between two given sequences.
	if (argc < 3) {
		fprintf(stderr, "usage: align-sw <seq1> <seq2>\n");
		exit(1);
	}

	Align< String<char> > ali;
	appendValue(rows(ali), argv[1]);
	appendValue(rows(ali), argv[2]);
    ::std::cout << "Score = " << localAlignment(ali, Score<int>(3,-3,-2, -2), SmithWaterman()) << ::std::endl;
	::std::cout << ali;
	::std::cout << "Aligns Seq1[" << clippedBeginPosition(row(ali, 0)) << ":" << (clippedEndPosition(row(ali, 0))-1) << "]";
	::std::cout << " and Seq2[" << clippedBeginPosition(row(ali, 1)) << ":" <<  (clippedEndPosition(row(ali, 1))-1) << "]" << ::std::endl << ::std::endl;
	return 0;
}
