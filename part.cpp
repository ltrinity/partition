#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//matrix printer
void printMatrix(const vector<vector<double>>& M, const string& sequence){
    //print the sequence
    cout << "Input Sequence: " << sequence << endl;
    cout  << "\t";
    for (int i = 0; i < sequence.length(); i++) {
        cout << sequence[i] << "\t";
    }
    cout << "\n";
    //print the matrix
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

//recursive traceback method
void traceback(const vector<vector<double>>& M, const string& sequence, int i, int j, string dotBracket, bool leftTraceAllowed) {
    //empty, also not possible in V or not allowed
    if (i == -1 or j ==-1) {
        return;
    }
    // i unpaired
    traceback(M, sequence, i-1, j, dotBracket, false);

    // j unpaired with unambiguous trace back
    //avoids left up and up left reaching same cell multiple ways
    if(leftTraceAllowed){
        traceback(M, sequence, i, j-1, dotBracket, true);
    }

    //paired
    if (M[i][j] > 0) {
        //for visualizing structures
        string newDotBracket = dotBracket;
        newDotBracket[i] = '(';
        newDotBracket[j] = ')';
        cout << newDotBracket << endl;
        traceback(M, sequence, i-1, j-1, newDotBracket, true);
    }
}

// Compute the pf using simple energy model
void computePartitionFunction(const string& sequence) {

    int n = sequence.length();
    int minHairpinUnpaired = 3;

    // Initialize V matrix
    vector<vector<double>> V(n, vector<double>(n, 0.0));
    // Initialize WM1  matrix
    vector<vector<double>> WM1(n, vector<double>(n, 0.0));

    // Compute V matrix
    for (int j = minHairpinUnpaired + 1; j < n; j++) {
        for (int i = j - minHairpinUnpaired - 1; i >= 0; i--) {
            auto result = canPair(sequence[i], sequence[j]);
            bool pair_possible = result.first;
            int pair_energy = result.second;
            if (pair_possible) {
                //hairpin and stack
                int hairpinLoopSize = j-i-1;
                int stackedPairEnergy = V[i + 1][j - 1];
                V[i][j] = 0.1*hairpinLoopSize + pair_energy + stackedPairEnergy;
                
                //internal bulge
                for(int m = j-2; m > 0; m--){
                    int k = i + 1;
                    int l = m;
                    while (l <= j-1 && k < l){
                        if(V[k][l] > 0){
                            //unpaired count in between the pairs
                            int internalLoopSize = (k-i-1) + (j-l-1);
                            V[i][j] += V[k][l] + internalLoopSize*0.1;
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

    cout << "V Matrix:" << endl;
    printMatrix(V, sequence);

    // Print the structures in V from traceback
    cout << sequence << endl;
    traceback(V, sequence, n-2, n-1, string(n, '.'), true);

    //Compute WM1 matrix
    //for(int i =0; i < n -1;i++){
    //    for (int j = i+1; j < n ; j++) {
    //        WM1[i][j] += WM1[i][j-1] + V[i][j];
    //    }
    //}

    //cout << "WM1 Matrix:" << endl;
    //printMatrix(WM1, sequence);


    // Initialize W matrix
    vector<vector<double>> W(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; i++) {
        W[i][i] = 1.0;
    }

    //Compute W matrix
    for (int j = 1; j < n ; j++) {
        W[0][j] += W[0][j-1];
        for (int r = 0; r < j; r++){
            W[0][j] += max(W[0][r-1],1.0)*V[r][j];
        }
    }
    cout << "W Matrix:" << endl;
    printMatrix(W, sequence);
}

int main() {
    //basic hairpin
    string sequence = "CAAAG";
    //basic stack
    //string sequence = "CCAAAGG";
    //basic internal
    //string sequence = "CCAAAGAG";
    computePartitionFunction(sequence);

    return 0;
}
