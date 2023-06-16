#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstring>

using namespace std;

bool DEBUG = false;

// constraints
int minLoop = 3;

// penalties
float unpairedCost = 0.1;
float unpairedMLCost = 0.1;
float bpMLCost = 0.5;
float MLinitCost = 1;

//constants
float temp = 37.0;
float gasConstant = 0.0821;

// index vector saves space
// we only need the upper triangle
vector<int> indexVec;

// calculate matrix index vis a vis vector
int matInd(int i, int j){
    return (indexVec[i] + j - i);
}

// boltzmann function
float Bolt(float energy){
    return exp(-energy/(temp*gasConstant));
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
            cout << std::setprecision(2) << M[indexVec[i] + j] << "\t";
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
        indexVec[i] += indexVec[i-1]+rowInc;
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
	            V[ijLoc] = Bolt(unpairedCost*hairpinLoopSize + hairpin_energy);
                // stacked pair
                auto stackedPair = canPair(sequence[i+1], sequence[j-1]);
                bool stack_possible = stackedPair.first;
                int stackedPairEnergy = stackedPair.second;
                if(stack_possible){
                    count += 1;
                    float stackSumEnergy = V[matInd(i+1,j-1)] * Bolt(stackedPairEnergy);
                    V[ijLoc] += stackSumEnergy;
                }
                // internal bulge
                for(int m = j-2; m > 0; m--){
                    int k = i + 1;
                    int l = m;
                    // possibly remove this check
                    while (l <= j-1 && k < l-minLoop){
                        auto internalPair = canPair(sequence[k], sequence[l]);
                        bool internal_possible = internalPair.first;
                        int internalEnergy = internalPair.second;
                        if(internal_possible){
                            count += 1;
                            // unpaired count in between the pairs
                            int internalLoopSize = (k-i-1) + (j-l-1);
                            float internalSumEnergy = V[matInd(k,l)] * Bolt(internalEnergy + internalLoopSize*unpairedCost);
                            V[ijLoc] += internalSumEnergy;
                        }
                        k++;
                        l++;
                    }
                }
                // multiloop
                // terminal (rightmost) branch
                WM1[ijLoc] += (V[ijLoc] * Bolt(bpMLCost));
                
                // option for j unpaired and penalty inside of multiloop
                int ijsub1Loc = matInd(i, j-1);

                // move to terminal branch i,j-1 (j unpaired)
                WM1[ijLoc] += (WM1[ijsub1Loc] * Bolt(unpairedMLCost));

                // unpaired base between branches (j unpaired)
                WM[ijLoc] += (WM[ijsub1Loc] * Bolt(unpairedMLCost));

                // option for additional branch 
                for (int r = i; r < j-1; r++){
                    //initial branch
                    WM[ijLoc] += (V[matInd(r,j)] * Bolt(unpairedMLCost*(r-i) + bpMLCost));

                    // intermediate branch
                    WM[ijLoc] += (WM[matInd(i,r)] * V[matInd(r+1,j)] * Bolt(bpMLCost));
                }

                // at least two branches required for multiloop
                for (int h = i+2; h <= j-1; h++){
                        VM[ijLoc] += (WM[matInd(i+1,h-1)] * WM1[matInd(h,j-1)] * Bolt(bpMLCost + MLinitCost));
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
         // Print W
        cout << "\t";
        for (int i = 0; i < sequence.length(); i++) {
            cout << sequence[i]<< "\t" ;
        }
        cout << "\n\t";
        for (int i = 0; i < sequence.length(); i++) {
            cout << W[i] << "\t" ;
        }
        cout << "\n";
    }
    
    cout << "count: " << count << endl;

    cout << "partition function: " << W[indexVec[0]+(n-1)] << endl;

    // PIPE DOTPLOT VIA Gnuplot
    // Open Gnuplot
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    //output png
    fprintf(gnuplotPipe, "set terminal png; set output 'plot.png'\n");
    //set acis
    fprintf(gnuplotPipe, "set yrange [0:%d]; set xrange [0:%d]\n",n,n);
    //set ticks
    fprintf(gnuplotPipe, "set tics scale 0;set xtics format ''\n");
    fprintf(gnuplotPipe, "set x2tics(");
    for(float i = 0;i< n-1;i++){
        fprintf(gnuplotPipe, "'%c' %f, ",sequence.at(i),i+0.5) ;   
    }
    fprintf(gnuplotPipe, "'%c' %f)\n",sequence.at(n-1),n-0.5) ;  
    fprintf(gnuplotPipe, "set ytics(");
    for(float i = 0;i< n-1;i++){
        fprintf(gnuplotPipe, "'%c' %f, ",sequence.at(n-i-1),i+0.5) ;   
    }
    fprintf(gnuplotPipe, "'%c' %f)\n",sequence.at(0),n-0.5) ;  
    fprintf(gnuplotPipe, "set grid ytics x2tics\n"); 
    //init rect from V/W[0,n-1] with sqrt scaling
    int idx = 0;
    for(int i=n-1; i >= 0; i--) {
      for(int j=n-i-1; j < n; j++) {
          double size = sqrt(V[idx++]/W[indexVec[0]+(n-1)]);
          fprintf(gnuplotPipe, "set object %d rect center %f, %f size %f, %f fc black\n", idx, j+0.5, i+0.5, size, size);

      }
    }
    fprintf(gnuplotPipe, "plot 0 notitle\n");

    // Close the pipe
    pclose(gnuplotPipe);
}

int main() {
    std::cout << std::fixed;
    std::cout << std::setprecision(2);
    // basic hairpin
    string seq = "CAAAG";
    // basic stack
    //string sequence = "CCAAAGG";
    // basic internal
    //string sequence = "CACAAAGG";
    // basic multiloop
    //string seq = "CCAAAGCAAAGG";
    int n = seq.length();
    indexVec.resize(n, 0);
    cout << seq << endl;
    computePartitionFunction(seq, n, minLoop);
    return 0;
}
