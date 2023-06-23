/**
* @mainpage
*
* Conditional partition function for density-2 RNA pseudoknots
*/
#include "base_types.hh"
#include "cmdline.hh"
#include "cmdline.cc"
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <cstring>
#include <string>
#include <cassert>
extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}

class CParty;

/**
* Methods for the computation of the partition function
* and associated structure classes.
*/
class CParty {
	public:
		std::string seq_;

		short *S_;
		short *S1_;

		paramT *params_;

        bool garbage_collect_;
        CParty(const std::string &seq, bool garbage_collect)
        : seq_(seq), params_(scale_parameters()), garbage_collect_(garbage_collect)
        {
        make_pair_matrix();

        S_ = encode_sequence(seq.c_str(),0);
        S1_ = encode_sequence(seq.c_str(),1);
        }

        

        ~CParty() {
        free(params_);
        free(S_);
        free(S1_);
        }
    };

/**
* @brief Simple driver for CParty
*
* Reads sequence from command line or stdin and calls folding methods.
*/
int
main(int argc,char **argv) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
	std::getline(std::cin,seq);
	}


	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}
	

	bool verbose;
	verbose = args_info.verbose_given;


	CParty cparty(seq,!args_info.noGC_given);

    cmdline_parser_free(&args_info);

	std::cout << seq << std::endl;
	
	if (verbose) {
		std::cout << "Hi"<<std::endl;
	}

	return 0;
}