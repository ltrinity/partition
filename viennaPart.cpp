#include <cstring>
#include <string>
#include <iostream>
extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}

//! type of base pair
typedef int32_t pair_type;
//! type of energy
typedef int32_t energy_t;
//! type of position
typedef int32_t cand_pos_t;

/**
 * From SparseRNAFolD
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
energy_t HairpinE(auto const& seq, auto const& S, auto const& S1, auto const& params, cand_pos_t i, cand_pos_t j) {
	
	const pair_type ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
}

int main(){
    vrna_params_load_defaults();
    //scale_parameters();
    make_pair_matrix();
    // Initialize RNA folding parameters
    vrna_md_t md;
    paramT* params = get_scaled_parameters(37.0,md);
    std::string seq_="CAAAG";
    auto const& S_ = encode_sequence(seq_.c_str(),0);
	auto const& S1_ = encode_sequence(seq_.c_str(),1);
    
    for(int i=0; i < 6;i++){
        for(int j=0; j < 6;j++){
            const pair_type ptype_closing = pair[S_[i]][S_[j]];
            if (ptype_closing!=0){
                energy_t hairpinEnergy = HairpinE(seq_,S_,S1_,params,i,j);
                std::cout << hairpinEnergy << std::endl;
            }
        }
    }
    return 0;
}
