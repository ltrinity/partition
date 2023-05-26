#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

//matrix printer
void printMatrix(const vector<vector<int>>& M, const string& sequence){
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
bool canPair(char base1, char base2) {
    // if bases form a canonical pair
    if ((base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A') ||
        (base1 == 'G' && base2 == 'C') || (base1 == 'C' && base2 == 'G')) {
        return true;
    } else {
        return false;
    }
}

//recursive traceback method
void traceback(const vector<vector<int>>& M, const string& sequence, int i, int j, string dotBracket, bool leftTraceAllowed) {
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
    vector<vector<int>> V(n, vector<int>(n, 0));

    // Compute V matrix
    for (int j = minHairpinUnpaired + 1; j < n; j++) {
        for (int i = j - minHairpinUnpaired - 1; i >= 0; i--) {
            if (canPair(sequence[i], sequence[j])) {
                //hairpin and stack
                int hairpinLoopSize = j-i-1;
                int stackedPair = V[i + 1][j - 1];
                V[i][j] = hairpinLoopSize + stackedPair;
                
                //internal bulge
                for(int m = j-2; m > 0; m--){
                    int k = i + 1;
                    int l = m;
                    while (l <= j-1 && k < l){
                        if(V[k][l] > 0){
                            //unpaired count in between the pairs
                            int internalLoopSize = (k-i-1) + (j-l-1);
                            V[i][j] += V[k][l] + internalLoopSize;
                        }
                        k++;
                        l++;
                    }
                }
            }
        }
    }

    cout << "V Matrix:" << endl;
    printMatrix(V, sequence);

    // Print the structures in V from traceback
    cout << sequence << endl;
    traceback(V, sequence, n-2, n-1, string(n, '.'), true);

    // Initialize W matrix
    vector<vector<int>> W(n, vector<int>(n, 0));
    for (int i = 0; i < n; i++) {
        W[i][i] = 1;
    }

    //Compute W matrix
    for (int j = 1; j < n ; j++) {
        W[0][j] += W[0][j-1];
        for (int r = 0; r < j; r++){
            W[0][j] += max(W[0][r-1],1)*V[r][j];
        }
    }
    cout << "W Matrix:" << endl;
    printMatrix(W, sequence);
}

int main() {
    //basic hairpin
    //string sequence = "CAAAG";
    //basic stack
    //string sequence = "CCAAAGG";
    //basic internal
    string sequence = "CCAAAGAG";
    computePartitionFunction(sequence);

    return 0;
}
