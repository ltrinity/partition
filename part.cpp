#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// matrix printer function
// for future work, print 1d vector as a matrix
void printMatrix(const vector<vector<double>>& M, const string& sequence){
    // print the sequence
    for (int i = 0; i < sequence.length(); i++) {
        cout << "\t" << sequence[i];
    }
    cout << "\n";
    // print the matrix
    for (int i = 0; i < sequence.length(); i++) {
        cout << sequence[i] << "\t" ;
        for (int j = 0; j < sequence.length(); j++) {
            cout << M[i][j] << "\t";
        }
        cout << endl;
    }
}

// canonical pair checker
std::pair<bool, double> canPair(char base1, char base2) {
    // two or three hydrogen bonds for AU or CG respectively
    if ((base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A')){
        return std::make_pair(true, -2.0);
    } else if ((base1 == 'G' && base2 == 'C') || (base1 == 'C' && base2 == 'G')) {
        return std::make_pair(true, -3.0);
    } else {
        return std::make_pair(false, 0.0);
    }
}

// traceback method
void computeTrace(const string& sequence, int n, int minLoop, vector<vector<double>> V, vector<vector<double>> WM1) {
    // iterate through each column and row
    for (int j = minLoop + 1; j < n; j++) {
        for (int i = j - minLoop - 1; i >= 0; i--) {
            // check for pair
            auto hairpin = canPair(sequence[i], sequence[j]);
            bool hairpin_possible = hairpin.first;
            if (hairpin_possible) {
                string db = string(n, '.');
                db[i] = '(';
                db[j] = ')';
                //stacked pair
                auto stackedPair = canPair(sequence[i+1], sequence[j-1]);
                bool stack_possible = stackedPair.first;
                double stack_energy = stackedPair.second;
                if(stack_possible){
                    cout << db << " " << V[i][j]-V[i+1][j-1]-stack_energy  << endl;
                    db[i+1] = '(';
                    db[j-1] = ')';
                    cout << db << " " << V[i][j]-V[i+1][j-1]-0.2 << endl;
                } else {
                    cout << db << " " << V[i][j] << endl;
                }
                //internal bulge
                for(int m = j-2; m > 0; m--){
                    int k = i + 1;
                    int l = m;
                    while (l <= j-1 && k < l){
                        auto internalPair = canPair(sequence[k], sequence[l]);
                        bool internal_possible = internalPair.first;
                        if(internal_possible){
                            db[k] = '(';
                            db[l] = ')';
                            cout << db << " " << V[i][j] << endl;
                        }
                        k++;
                        l++;
                    }
                }
                //multiloop
            }
        }
    }
}

// Compute the pf using simple energy model
void computePartitionFunction(const string& sequence, int n, int minLoop) {

    // Initialize V matrix
    vector<vector<double>> V(n, vector<double>(n, 0.0));
    // Initialize WM1  matrix
    vector<vector<double>> WM1(n, vector<double>(n, 0.0));

    // Compute V matrix
    for (int j = minLoop + 1; j < n; j++) {
        for (int i = j - minLoop - 1; i >= 0; i--) {
            auto hairpin = canPair(sequence[i], sequence[j]);
            bool hairpin_possible = hairpin.first;
            int hairpin_energy = hairpin.second;
            if (hairpin_possible) {
                //hairpin and stack
                int hairpinLoopSize = j-i-1;
                V[i][j] = 0.1*hairpinLoopSize + hairpin_energy;
                //stacked pair
                auto stackedPair = canPair(sequence[i+1], sequence[j-1]);
                bool stack_possible = stackedPair.first;
                int stackedPairEnergy = stackedPair.second;
                if(stack_possible){
                    V[i][j] += V[i+1][j-1] + stackedPairEnergy;
                }
                //internal bulge
                for(int m = j-2; m > 0; m--){
                    int k = i + 1;
                    int l = m;
                    while (l <= j-1 && k < l){
                        auto internalPair = canPair(sequence[k], sequence[l]);
                        bool internal_possible = internalPair.first;
                        int internalEnergy = internalPair.second;
                        if(internal_possible){
                            //unpaired count in between the pairs
                            int internalLoopSize = (k-i-1) + (j-l-1);
                            V[i][j] += V[k][l] + internalEnergy + internalLoopSize*0.1;
                        }
                        k++;
                        l++;
                    }
                }

                //multiloop
                WM1[i][j] += WM1[i][j-1] + V[i][j];
            }
        }
    }

    cout << "V" << endl;
    printMatrix(V, sequence);

    // Print the structures in V from traceback
    cout << sequence << endl;
    computeTrace(sequence, n, minLoop, V, WM1);

    //cout << "WM1 Matrix:" << endl;
    //printMatrix(WM1, sequence);

    // Initialize W matrix
    vector<double>  W(n, 0.0);
    //base case
    W[0] = 1.0;

    //Compute W matrix
    for (int j = 1; j < n ; j++) {
        W[j] += W[j-1];
        for (int r = 0; r < j; r++){
            W[j] += max(W[r-1],1.0)*V[r][j];
        }
    }

    //print W
    cout << "W" <<endl;
    for (int i = 0; i < sequence.length(); i++) {
        cout << "\t" << W[i] ;
    }
    cout << endl;
}

int main() {
    //basic hairpin
    //string sequence = "CAAAG";
    //basic stack
    //string sequence = "CCAAAGG";
    //basic internal
    string sequence = "CCAAAGAG";
    int n = sequence.length();
    //minloopsize
    int minLoop = 3;
    computePartitionFunction(sequence, n, minLoop);
    return 0;
}
