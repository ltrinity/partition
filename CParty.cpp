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
		cand_pos_t n_;

		short *S_;
		short *S1_;

		paramT *params_;

        bool garbage_collect_;

		std::vector<energy_t> V_;

		CParty(const std::string &seq, bool garbage_collect, std::string restricted)
		: seq_(seq),n_(seq.length()),params_(scale_parameters()),garbage_collect_(garbage_collect){

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
 	* @brief Determines the MFE energy for a given sequence
	*/
	energy_t fold(auto const& seq, auto &V, auto const& S, auto const& S1, auto const& params, auto const& n, auto const& garbage_collect,const cand_pos_t*p_table,const cand_pos_t*last_j_array,const cand_pos_t*in_pair_array,const cand_pos_t*up_array) {
		std::cout << "fold" << std::endl;
		// for (cand_pos_t i=n; i>0; --i) {
		// 	for (cand_pos_t j=i+TURN+1; j<=n; j++ ) {
		// 		bool evaluate = evaluate_restriction(i,j,last_j_array,in_pair_array);
		// 		// ------------------------------
		// 		// W: split case
		// 		bool pairedkj = 0;
		// 		energy_t w_split = INF;
		// 		energy_t wm_split = INF;
		// 		energy_t wm2_split = INF;
		// 		for ( auto it=CL[j].begin();CL[j].end() != it;++it ) {
		// 			const cand_pos_t k=it->first;
		// 			// Decode the energies
		// 			const energy_t v_kj = it->third >> 2;
		// 			const energy_t v_kjw = it->fourth >> 2;

		// 			const bool can_pair = up_array[k-1] >= (k-i) ? true: false;
					
		// 			wm_split = std::min( wm_split, WM[k-1] + v_kj );
		// 			if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
		// 			wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
		// 			w_split = std::min( w_split, W[k-1] + v_kjw );		
					
		// 		}
		// 		if(p_table[j]<0) w_split = std::min(w_split,W[j-1]);
		// 		if(p_table[j]<0) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
		// 		if(p_table[j]<0) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
				
		// 		energy_t w  = w_split; // entry of W w/o contribution of V
		// 		energy_t wm = wm_split; // entry of WM w/o contribution of V


		// 		cand_pos_t i_mod= (uint_least32_t)i%(MAXLOOP+1);

		// 		const pair_type ptype_closing = pair[S[i]][S[j]];
		// 		const bool restricted = p_table[i] == -1 || p_table[j] == -1;

		// 		const bool unpaired = (p_table[i]<-1 && p_table[j]<-1);
		// 		const bool paired = (p_table[i] == j && p_table[j] == i);
		// 		energy_t v = INF;
		// 		// ----------------------------------------
		// 		// cases with base pair (i,j)
		// 		if(ptype_closing>0 && evaluate && !restricted) { // if i,j form a canonical base pair
		// 			bool canH = (paired || unpaired);
		// 			if(up_array[j-1]<(j-i-1)) canH=false;
					
		// 			energy_t v_h = canH ? HairpinE(seq,S,S1,params,i,j) : INF;
		// 			// info of best interior loop decomposition (if better than hairpin)
		// 			cand_pos_t best_l=0;
		// 			cand_pos_t best_k=0;
		// 			energy_t best_e;

		// 			energy_t v_iloop=INF;

		// 			cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
		// 			if((p_table[i]<-1 && p_table[j] < -1) || p_table[i] == j) { 
		// 				for ( cand_pos_t k=i+1; k<=max_k; k++) {
		// 					cand_pos_t k_mod= k%(MAXLOOP+1);
							
		// 					energy_t cank = (up_array[k-1]>=(k-i-1)-1);
		// 					cand_pos_t min_l=std::max(k+TURN+1, k+j-i- MAXLOOP-2);
							
		// 					for (size_t l=j-1; l>=min_l; --l) {
		// 						assert(k-i+j-l-2<=MAXLOOP);
		// 						energy_t canl = ((up_array[j-1]>=(j-l-1)-1) | cank);
		// 						energy_t v_iloop_kl = INF & canl;
								
		// 						v_iloop_kl = v_iloop_kl + V(k_mod,l) + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
		// 						// v_iloop_kl = v_iloop_kl + V(k_mod,l) + ILoopE(S,S1,params,ptype_closing,i,j,k,l);

		// 						if ( v_iloop_kl < v_iloop) {
		// 							v_iloop = v_iloop_kl;
		// 							best_l=l;
		// 							best_k=k;
		// 							best_e=V(k_mod,l);
		// 						}
		// 					}
							
		// 				}
		// 			}
		// 			const energy_t v_split = E_MbLoop(dmli1,dmli2,S,params,i,j,p_table);

		// 			v = std::min(v_h,std::min(v_iloop,v_split));
		// 			// register required trace arrows from (i,j)
		// 			if ( v_iloop < std::min(v_h,v_split) ) {
		// 				if ( is_candidate(CL,cand_comp,best_k,best_l) ) {
		// 					avoid_trace_arrow(ta);
		// 				} else {
		// 					register_trace_arrow(ta,i,j,best_k,best_l,best_e);
		// 				}
		// 			}
		// 			V(i_mod,j) = v;
		// 		} else {
		// 			V(i_mod,j) = INF;
		// 		} // end if (i,j form a canonical base pair)
							
		// 		cand_pos_t ip1_mod = (uint_least32_t)(i+1) %(MAXLOOP+1);
		// 		energy_t vi1j = V(ip1_mod,j);
		// 		energy_t vij1 = V(i_mod,j-1);
		// 		energy_t vi1j1 = V(ip1_mod,j-1);
					
	
		// 		w  = std::min(w_v, w_split);
		// 		wm = std::min(wm_v, wm_split);
		// 		if ( w_v < w_split || wm_v < wm_split || paired) {
		// 			cand_pos_t k_mod = (uint_least32_t) k%(MAXLOOP+1);
		// 			// Encode the dangles into the energies
		// 			energy_t w_enc = (w_v << 2) | d;
		// 			energy_t wm_enc = (wm_v << 2) | d;
		// 			register_candidate(CL, i, j,V(i_mod,j), wm_enc,w_enc);
		// 			// always keep arrows starting from candidates
		// 			inc_source_ref_count(ta,i,j);
		// 		}		

		// 		W[j]       = w;
		// 		WM[j]      = wm;
		// 		WM2[j]     = wm2_split;

		// 	} // end loop j
		// }
		// return W[n];
	}


/**
 * @brief Fills the restriction arrays from SPARSEMFEFOLD
 * p_table will contain the index of each base pair
 * X or x tells the program the base cannot pair and . sets it as unpaired but can pair
 * Pseudoknots (denoted by [ ], < >, or { } ) are filled the same way as ( )
 * That is, a structure like this (<)> is not possible.
 * @param structure Input structure
 * @param p_table Restricted array
 * @param last_j_array Restricted Array
 * @param in_pair_array Restricted Array
 */
void detect_restricted_pairs(auto const &structure, cand_pos_t *p_table, cand_pos_t *last_j_array, cand_pos_t *in_pair_array){
	cand_pos_t i, j, count = 0, length = structure.length(),last_j=length;
	std::vector<cand_pos_t>  pairs;
	pairs.push_back(length);

	for (i=length; i >=1; --i){
		if ((structure[i-1] == 'x') || (structure[i-1] == 'X'))
			p_table[i] = -1;
		else if (structure[i-1] == '.')
			p_table[i] = -2;
		if (structure[i-1] == ')' || structure[i-1] == ']' || structure[i-1] == '}' || structure[i-1] == '>'){
			pairs.push_back(i);
			count++;
		}
		last_j_array[i] = pairs[pairs.size()-1];
		in_pair_array[i] = count;
		if (structure[i-1] == '(' || structure[i-1] == '[' || structure[i-1] == '{' || structure[i-1] == '<'){
			j = pairs[pairs.size()-1];
			pairs.erase(pairs.end()-1);
			p_table[i] = j;
			p_table[j] = i;
			count--;
		}
	}
	pairs.pop_back();
	if (pairs.size() != 0)
	{
		fprintf (stderr, "The given structure is not valid: more left parentheses than right parentheses: \n");
		exit (1);
	}
}

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
	cand_pos_t n = seq.length();
	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (n,'.');
	if(restricted != "" && restricted.length() != n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		exit(0);
	}
	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}
	

	bool verbose;
	verbose = args_info.verbose_given;


	CParty cparty(seq,!args_info.noGC_given,restricted);

	// Make replicate mx array in linear space
	cand_pos_t last_j_array[n+1] = {0};
	cand_pos_t in_pair_array[n+1] = {0};
	cand_pos_t p_table[n+1] = {0};
	cand_pos_t up_array[n+1] = {0};
	
	cmdline_parser_free(&args_info);

	std::cout << seq << std::endl;
	
	detect_restricted_pairs(restricted,p_table,last_j_array,in_pair_array);
	
	energy_t mfe = fold(cparty.seq_,cparty.V_,cparty.S_,cparty.S1_,cparty.params_,cparty.n_,cparty.garbage_collect_, p_table,last_j_array,in_pair_array,up_array);		
	

	if (verbose) {
		std::cout << "Hi"<<std::endl;
	}

	return 0;
}