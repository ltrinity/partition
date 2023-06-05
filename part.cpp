#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

//penalties
float unpairedCost = 0.1;
float unpairedMLCost = 0.1;
float bpMLCost = 0.5;
float MLinitCost = 1;

// add brackets for structure output
string addBrackets(string structure, int i, int j){
    structure[i] = '(';
    structure[j] = ')';
    return structure;
}

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

// recursive traceback method
void traceback(const vector<vector<double>>& M, const string& sequence, int i, int j, string dotBracket, bool leftTraceAllowed) {
    // empty, also not possible in V or not allowed
    if (i == -1 or j ==-1) {
        return;
    }
    // i unpaired
    traceback(M, sequence, i-1, j, dotBracket, false);

    // j unpaired with unambiguous trace back
    // avoids left up and up left reaching same cell multiple ways
    if(leftTraceAllowed){
        traceback(M, sequence, i, j-1, dotBracket, true);
    }

    // paired
    if (M[i][j] < 0) {
        // for visualizing structures
        string dotBracketRecursive = addBrackets(dotBracket, i, j);
        cout << dotBracketRecursive << endl;
        traceback(M, sequence, i-1, j-1, dotBracketRecursive, true);
    }
}

// Compute the pf using simple energy model
void computePartitionFunction(const string& sequence, int n, int minLoop) {
    // Initialize V matrix
    vector<vector<double>> V(n, vector<double>(n, 0.0));
    // Initialize WM1  matrix
    vector<vector<double>> WM1(n, vector<double>(n, 0.0));
    // Initialize WM  matrix
    vector<vector<double>> WM(n, vector<double>(n, 0.0));
    // Initialize VM  matrix
    vector<vector<double>> VM(n, vector<double>(n, 0.0));

    // Compute V matrix
    for (int j = minLoop + 1; j < n; j++) {
        for (int i = j - minLoop - 1; i >= 0; i--) {
            auto hairpin = canPair(sequence[i], sequence[j]);
            bool hairpin_possible = hairpin.first;
            int hairpin_energy = hairpin.second;
            if (hairpin_possible) {
                // hairpin and stack
                int hairpinLoopSize = j-i-1;
                V[i][j] = unpairedCost*hairpinLoopSize + hairpin_energy;
                // stacked pair
                auto stackedPair = canPair(sequence[i+1], sequence[j-1]);
                bool stack_possible = stackedPair.first;
                int stackedPairEnergy = stackedPair.second;
                if(stack_possible){
                    float stackSumEnergy = V[i+1][j-1] + stackedPairEnergy;
                    V[i][j] += stackSumEnergy;
                }
                // internal bulge
                for(int m = j-2; m > 0; m--){
                    int k = i + 1;
                    int l = m;
                    while (l <= j-1 && k < l-minLoop){
                        auto internalPair = canPair(sequence[k], sequence[l]);
                        bool internal_possible = internalPair.first;
                        int internalEnergy = internalPair.second;
                        if(internal_possible){
                            // unpaired count in between the pairs
                            int internalLoopSize = (k-i-1) + (j-l-1);
                            float internalSumEnergy = V[k][l] + internalEnergy + internalLoopSize*unpairedCost;
                            V[i][j] += internalSumEnergy;
                        }
                        k++;
                        l++;
                    }
                }
                // multiloop
                if(V[i][j]<0){
                    WM1[i][j] += (V[i][j] + bpMLCost);
                }
                if(WM1[i][j-1]<0){
                    WM1[i][j] += (WM1[i][j-1] + unpairedMLCost);
                }
                if(WM[i][j-1]<0){
                    WM[i][j] += (WM[i][j-1] + unpairedMLCost);
                }         
                for (int r = i; r < j-1; r++){
                    //terminal branch
                    if(V[r][j] < 0){
                        WM[i][j] += ((unpairedMLCost*(r-i)) + V[r][j] + bpMLCost);
                    }
                    //intermediate branch
                    if(WM[i][r] + V[r+1][j] < 0){
                        WM[i][j] += (WM[i][r] + V[r+1][j] + bpMLCost);
                    }
                }
                for (int h = i+2; h <= j-1; h++){
                    if(WM[i+1][h-1] < 0 && WM1[h][j-1] < 0){
                        VM[i][j] += (WM[i+1][h-1] + WM1[h][j-1] + bpMLCost + MLinitCost);
                        V[i][j] += (VM[i][j]);
                    }
                }

            }
        }
    }

    // Initialize W matrix
    vector<vector<double>> W(n, vector<double>(n, 0.0));
    // base case
    for (int i = 0; i < n; i++){
        W[i][i] = 1.0;
    }

    // Compute W matrix
    for (int i = 0; i <= n-1; i++){
        for (int j = 1; j < n ; j++) {
            W[i][j] += W[i][j-1];
            for (int r = i; r < j; r++){
                W[i][j] += max(W[i][r-1],1.0)*V[r][j];
            }
        }
    }

    //cout << "V" << endl;
    //printMatrix(V, sequence);

    //cout << "WM1" << endl;
    //printMatrix(WM1, sequence);

    //cout << "WM" << endl;
    //printMatrix(WM, sequence);

    //cout << "VM" << endl;
    //printMatrix(VM, sequence);

    // print W
    cout << "W" <<endl;
    printMatrix(W, sequence);
    //for (int i = 0; i < sequence.length(); i++) {
    //    cout << "\t" << W[i] ;
    //}
    //cout << endl;
    
    traceback(V, sequence, n-2, n-1, string(n, '.'), true);
}

int main() {
    // basic hairpin
    string sequence = "CAAAG";
    // basic stack
    //string sequence = "CCAAAGG";
    // basic internal
    //string sequence = "CCAAAGCAAAGG";
    int n = sequence.length();
    // minloopsize
    int minLoop = 3;
    computePartitionFunction(sequence, n, minLoop);
    return 0;
}
