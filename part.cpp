#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

bool DEBUG = true;

// constraints
int minLoop = 3;

// penalties
float unpairedCost = 0.1;
float unpairedMLCost = 0.1;
float bpMLCost = 0.5;
float MLinitCost = 1;

//constants
float temperature = 37.0;
float gasConstant = 0.0821;

// index vector saves space
// we only need the upper triangle
int* index;

// calculate matrix index vis a vis vector
int matInd(int i, int j){
    return (index[i] + j - i);
}

// boltzmann function
float B(float energy){
    return exp(-energy/(temperature*gasConstant));
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

// vectir printer function
void printVector(const double* M, string sequence){
    // Print the sequence
    cout << "\t";
    for (int i = 0; i < sequence.length(); i++) {
        cout << sequence[i]<< "\t" ;
    }
    cout << "\n";

    int size = sequence.length();
    // Print the matrix
    for (int i = 0; i < size; i++) {
        cout << sequence[i] << "\t";
        //empty spaces for bottom triangle
        int x = i;
        while( x > 0){
            cout << "\t";
            x--;
        }
        for (int j = 0; j < size-i; j++) {
            cout << std::setprecision(2) << M[index[i] + j] << "\t";
        }
        cout << endl;
    }
}



// Compute the pf using simple energy model
void computePartitionFunction(string sequence, int n, int minLoop) {

    //fill index vector
    int vectorLength = (n*(n+1))/2;
    int rowInc = n;
    for(int i = 1; i <= n; i++){
        index[i] += index[i-1]+rowInc;
        rowInc--;
    }
    
    // Initialize V vector
    double* V = new double[vectorLength];
    // Initialize WM1 vector
    double* WM1 = new double[vectorLength];
    // Initialize WM vector
    double* WM = new double[vectorLength];
    // Initialize VM vector
    double* VM = new double[vectorLength];

    int count = 0;

    // Compute V matrix
    for (int i = n-1-minLoop; i >= 0; i--) {
        for (int j = n-1; j >= i + minLoop; j--) {
            auto hairpin = canPair(sequence[i], sequence[j]);
            bool hairpin_possible = hairpin.first;
            int hairpin_energy = hairpin.second;
            if (hairpin_possible) {
                // hairpin
                count += 1;
                int ijLoc = matInd(i,j);
                int hairpinLoopSize = j-i-1;
                V[ijLoc] = B(unpairedCost*hairpinLoopSize + hairpin_energy);
                // stacked pair
                auto stackedPair = canPair(sequence[i+1], sequence[j-1]);
                bool stack_possible = stackedPair.first;
                int stackedPairEnergy = stackedPair.second;
                if(stack_possible){
                    count += 1;
                    float stackSumEnergy = V[matInd(i+1,j-1)] * B(stackedPairEnergy);
                    V[ijLoc] += stackSumEnergy;
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
                            count += 1;
                            // unpaired count in between the pairs
                            int internalLoopSize = (k-i-1) + (j-l-1);
                            float internalSumEnergy = V[matInd(k,l)] * B(internalEnergy + internalLoopSize*unpairedCost);
                            V[ijLoc] += internalSumEnergy;
                        }
                        k++;
                        l++;
                    }
                }
                // multiloop
                // check if terminal (rightmost) branch possible
                WM1[ijLoc] += (V[ijLoc] * B(bpMLCost));
                
                //option for j unpaired and penalty inside of multiloop
                int ijsub1Loc = matInd(i, j-1);

                // move to terminal branch i,j-1 (j unpaired)
                WM1[ijLoc] += (WM1[ijsub1Loc] * B(unpairedMLCost));

                // unpaired base between branches (j unpaired)
                WM[ijLoc] += (WM[ijsub1Loc] * B(unpairedMLCost));

                // option for additional branch 
                for (int r = i; r < j-1; r++){
                    //initial branch
                    WM[ijLoc] += (V[matInd(r,j)] * B(unpairedMLCost*(r-i) + bpMLCost));

                    // intermediate branch
                    WM[ijLoc] += (WM[matInd(i,r)] * V[matInd(r+1,j)] * B(bpMLCost));
                }

                // at least two branches required for multiloop
                for (int h = i+2; h <= j-1; h++){
                        VM[ijLoc] += (WM[matInd(i+1,h-1)] * WM1[matInd(h,j-1)] * B(bpMLCost + MLinitCost));
                        V[ijLoc] += (VM[ijLoc]);
                        if(VM[ijLoc] > 0){
                            count += 1;
                        }
                }
            }
        }
    }

    // Initialize W vector
    double* W = new double[sequence.length()];
    // base case
    W[0] = 1.0;
    // Compute W
    for (int j = 1; j < n; j++) {
        for (int r = 0; r < j; r++){
            W[j] += (max(W[r-1],1.0)*V[matInd(r,j)]);
        }
        W[j] += W[j-1];
    }
    

    if(DEBUG){
        cout << "V" << endl;
        printVector(V, sequence);

         cout << "\nWM1" << endl;
         printVector(WM1, sequence);

         cout << "\nWM" << endl;
         printVector(WM, sequence);

         cout << "\nVM" << endl;
         printVector(VM, sequence);

        cout << "\nW" <<endl;

         // Print the sequence
        cout << "\t";
        for (int i = 0; i < sequence.length(); i++) {
            cout << sequence[i]<< "\t" ;
        }
        cout << "\n\t";
        // Print the matrix
        for (int i = 0; i < sequence.length(); i++) {
            cout << W[i] << "\t" ;
        }
        cout << "\n";
        
    }
    
    cout << "count: " << count << endl;

    cout << "partition function: " << W[index[0]+(n-1)] << endl;
}

int main() {
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    // basic hairpin
    //string sequence = "CAAAG";
    // basic stack
    //string sequence = "CCAAAGG";
    // basic internal
    //string sequence = "CACAAAGG";
    // basic multiloop
    string sequence = "CCAAAGCAAAGG";
    int n = sequence.length();
    index = new int[n];
    cout << sequence << endl;
    computePartitionFunction(sequence, n, minLoop);
    return 0;
}
