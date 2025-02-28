#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;
typedef vector<float> vec;  // A vector of floating point numbers
typedef vector<vec> matrix; // A matrix: a vector of vectors

// Type Constants
enum class OptionType {
    PUT,
    CALL
};

class FiniteDifferenceMethod {
public:

    float timeHorizon = 0;                    // The number of years until expiry
    int timesteps = 0;                        // How many timesteps
    float maxS = 0;                           // Maximum stock price price, where the value of the option is considered to be the same as the price of the stock 
    int priceSteps = 0;                       // How many steps to divide the price axis into
    float strikePrice = 0;                    // The strike Price of the option
    float deltaT = 0;                         // The time delta accross each timestep
    float deltaS = 0;                         // The price delta for the price axis
    float volatility = 0;                     // The volatility for the noise (normally distributed)
    float riskFreeR  = 0;                     // The risk free interest rate that is used in the black scholes model
    float dividend = 0; ;                     // Placeholder for dividend percentage reductions
    float startPrice = 0;                     // Initial Price
    matrix grid;                              // The grid we will use for all methods;

    // Constructor
    FiniteDifferenceMethod(float timeHorizon, int timesteps, float maxS, int priceSteps, float strikePrice, float  volatility, float riskFreeR , float startPrice) : 
        timeHorizon(timeHorizon), timesteps(timesteps), maxS(maxS), priceSteps(priceSteps), strikePrice(strikePrice), volatility(volatility), riskFreeR (riskFreeR ), startPrice(startPrice) {
        deltaT = timeHorizon / timesteps;
        deltaS = maxS / priceSteps;
        // Define the grid, a vector of vectors of floats
        grid.resize(priceSteps, vector<float>(timesteps, 0.0f));

    };

    // Print the FDM grid to the console, for debugging
    void printGrid(int maxPriceSteps = 45, int maxTimesteps = 20) {
        for (int s = 0; s < priceSteps; s++) {
            for (int t = 0; t < timesteps; t++) {
                if (t % (timesteps / maxTimesteps) == 0 && s % (priceSteps / maxPriceSteps) == 0)
                    cout << setw(12) << setprecision(4) << grid[s][t];
            }
            if (s % (priceSteps / maxPriceSteps) == 0)
                cout << endl;
        }
        cout << endl;
    };

    // Static method for printing a matrix, for debugging
    static void printMatrix(matrix& m) {
        for (int s = 0; s < m.size(); s++) {
            for (int t = 0; t < m[0].size(); t++) {
                cout << setw(12) << setprecision(4) << m[s][t];
            };
            cout << endl;
        };
        cout << endl;
    };

    // Save grid to csv for use in MATLAB or Python, for debugging
    void saveGridToCSV() {
        cout << "Writing to CSV..." << endl;
        ofstream out;
        const string file = "./surfaceTest.csv";
        out.open(file, ios::out);
        for (int s = 0; s < grid.size(); s++) {
            for (int t = 0; t < grid[s].size(); t++) {
                out << grid[s][t];
                if (t == grid[s].size() - 1)
                    out << endl;
                else
                    out << ',';
            };
        };
        out.close();
    };

    // Get the interpolated price based on the current (start) price and the grid at t = T
    float getPrice() {
        float priceBelow = grid[int(startPrice / deltaS)][0]; // 
        float priceAbove = grid[int(startPrice / deltaS) + 1][0];
        float stockPriceTrunc = int(startPrice / deltaS) * (deltaS);
        float stockPriceCeil = int(startPrice / deltaS) * (deltaS) + deltaS;
        float price = priceBelow + (priceAbove - priceBelow) * (startPrice - stockPriceTrunc) / (stockPriceCeil - stockPriceTrunc);
        cout << "Price: " << price << endl;
        return price;
    };

    // Setup boundary conditions for the different option types
    void setBoundaries(OptionType optionType) {
        if (optionType == OptionType::PUT) {
            // 1. Final timestep price for a put max(strike - price, 0)
            for (int p = 0; p < priceSteps; p++) {
                grid[p][timesteps - 1] = max(strikePrice - (deltaS * p), 0.0f);
            };
            // 2. Top and bottom boundary (over all timesteps)
            for (int t = 0; t < timesteps; t++) {
                grid[0][t] = strikePrice * exp(-1 * riskFreeR * (timesteps - 1 - t) * deltaT);
                grid[priceSteps - 1][timesteps - 1 - t] = 0;
            };
        }
        else if (optionType == OptionType::CALL) {
            // 1. Final timestep price
            for (float p = 0; p < priceSteps; p++) {
                grid[p][timesteps - 1] = max((deltaS * p) - strikePrice, 0.0f);
            };
            // 2. Top and bottom boundary (over all timesteps)
            for (int t = 0; t < timesteps; t++) {
                grid[0][t] = 0;
                grid[priceSteps - 1][t] = maxS - strikePrice * exp(-1.0f * riskFreeR * (timesteps - 1 - t) * deltaT);
            };
        }
        else {
            // Check (for safety)
            throw runtime_error("Incorrect option types");
        }
    }

    // Explicit finite difference method
    void explicitMethod(OptionType optionType) {
        // Check for stability
        if (deltaT > pow(deltaS, 2) / (pow(volatility, 2) * pow(maxS, 2))) {
            cerr << "NOT STABLE" << endl;
            throw runtime_error("The deltas are unstable");
        };

        // Set boundaries
        setBoundaries(optionType);

        // Backwards marching
        // Initialise the coefficients 
        float A_n = 0.0f, B_n = 0.0f, C_n = 0.0f;
        // For each timestep (marchin backwards)
        for (int t = 1; t < timesteps; t++) {
            // For each node (not on the boundary)
            for (int s = 1; s < priceSteps - 1; s++) {

                A_n = 0.5 * (pow(volatility, 2) * pow(s, 2) - s * (riskFreeR  - dividend)) * deltaT;
                B_n = 1 - (riskFreeR  + pow(volatility, 2) * pow(s, 2)) * deltaT;
                C_n = 0.5 * (pow(volatility, 2) * pow(s, 2) + s * (riskFreeR  - dividend)) * deltaT;
                
                grid[s][timesteps - 1 - t] = 
                    A_n * grid[s - 1][timesteps - t] +
                    B_n * grid[s][timesteps - t] + 
                    C_n * grid[s + 1][timesteps - t];

                // Set to zero if less than threshold
                if (grid[s][timesteps - 1 - t] < pow(2, -10))
                    grid[s][timesteps - 1 - t] = 0;
            };
        };
    };

    // Matrix inversion method for use in the fully implicit and Crank Nicolson methods
    static vec gaussSeidelSolve(matrix& M, vec& w, int max_iterations = 1000, double tol = 1e-10) {
        // This function solves M.v = w iteratively using the Gauss-Seidel method.
        // M: The matrix to invert
        // w: the RHS of the equation, the column vector to multiply by the inverse
        int n = w.size();
        vec v(n, 0); // Initialize solution vector with zeros

        for (int iter = 0; iter < max_iterations; iter++) {
            vec v_old = v;  // Store previous iteration values

            for (int j = 0; j < n; j++) {
                float sum = w[j];

                if (j > 0)
                    sum -= M[j][j - 1] * v[j - 1];

                if (j < n - 1)
                    sum -= M[j][j + 1] * v[j + 1];

                // Check for division by zero
                if (abs(M[j][j]) < 1e-8) {
                    throw runtime_error("Zero pivot encountered in Gauss-Seidel method.");
                }

                v[j] = sum / M[j][j];
            }

            float error = 0;
            for (int j = 0; j < n; j++) 
                error = max(error, fabs(v[j] - v_old[j]));

            // Solution has converged
            if (error < tol)
                    break; 
                
        }

        return v;
    }

    // Fully Implicit finite difference method
    void implicitMethod(OptionType optionType) {
        
        // Set boundaries
        setBoundaries(optionType);

        // The matrix (initially sparce) we will populate and then invert
        matrix linearSystem(priceSteps, vector<float>(priceSteps, 0));

        // Create a_n, b_n and c_n for all prices and build matrix
        for (int s = 0; s < priceSteps; s++) {
            // populate a_n
            if (s > 0)
                linearSystem[s][s - 1] = -0.5 * (pow(volatility, 2) * pow(s, 2) - s * (riskFreeR - dividend)) * deltaT;

            // populate b_n
            linearSystem[s][s] = 1 + (riskFreeR + pow(volatility, 2) * pow(s, 2)) * deltaT;

            // populate c_n
            if (s < priceSteps - 1)
                linearSystem[s][s + 1] = -0.5 * (pow(volatility, 2) * pow(s, 2) + s * (riskFreeR - dividend)) * deltaT;
        }
        
        // Print matrix of coefficients
        //printMatrix(linearSystem);
        
        // Loop over all timesteps
        for (int t = 1; t < timesteps; t++) {
            
            // Get prices at timestep
            vec timeStepPrices(priceSteps, 0);
            for (int p = 0; p < priceSteps; p++) {
                timeStepPrices[p] = grid[p][timesteps - t];
            };

            vec newVec = gaussSeidelSolve(linearSystem, timeStepPrices);
            
            // Update grid
            for (int i = 1; i < priceSteps - 1; i++) {
                grid[i][timesteps - 1 - t] = newVec[i];
            }
            
        }
    }

    // Crank Nicolson - Hybrid explicit / implicit
    void crankNicolson(OptionType optionType) {
        
        // Set boundaries
        setBoundaries(optionType);

        // Initialise matrix A and B s.t. A*V_m = B*V_m+1
        matrix A(priceSteps, vector<float>(priceSteps, 0));
        matrix B(priceSteps, vector<float>(priceSteps, 0));

        // Create a_n, b_n and c_n for all prices and build matrix
        for (int s = 0; s < priceSteps; s++) {
            // populate a_n
            if (s > 0) {
                A[s][s - 1] = -0.25 * (pow(volatility, 2) * pow(s, 2) - s * (riskFreeR - dividend)) * deltaT; // a_n
                B[s][s - 1] = 0.25 * (pow(volatility, 2) * pow(s, 2) - s * (riskFreeR - dividend)) * deltaT; // A_n
            };

            // populate b_n
            A[s][s] = 1 + 0.5 * (riskFreeR + pow(volatility, 2) * pow(s, 2)) * deltaT; //b_n
            B[s][s] = 1 - 0.5 * (riskFreeR + pow(volatility, 2) * pow(s, 2)) * deltaT; // B_n

            // populate c_n
            if (s < priceSteps - 1) {
                A[s][s + 1] = -0.25 * (pow(volatility, 2) * pow(s, 2) + s * (riskFreeR - dividend)) * deltaT; // c_n
                B[s][s + 1] = 0.25 * (pow(volatility, 2) * pow(s, 2) + s * (riskFreeR - dividend)) * deltaT; // C_n
            };

        };

        // March backwards through the timesteps to update based on A*V_m-1 = B*V_m 
        // =>  V_m-1 = A^-1 * B * V_m  via Gauss Siedel

        for (int t = 0; t < timesteps - 1; t++) {
            // 1. find B.S_m

            // Vector, the result of B * S_m
            vec BV(priceSteps, 0);
            for (int x = 0; x < priceSteps; x++) {
                for (int y = 0; y < priceSteps; y++) {
                    BV[x] += B[x][y] * grid[y][timesteps - 1 - t];
                };
            };

            vec A_invBV = gaussSeidelSolve(A, BV);

            // FOR DEBUGGING THE LARGE NEGATIVE VALUE IN BS[-1]
            /*if (t == 0) {
                for (int z = 0; z < A_invBV.size(); z++) {
                    cout << setw(12) << A_invBV[z] << setw(12) <<  grid[z][timesteps - 1 - t] << setw(12) << BV[z] << setw(12) << B[z][z] << endl;
                }

                cout << B[B.size() - 1][B.size() - 2] << "    " << B[B.size() - 1][B.size() - 1] << endl;
                
            }*/

            //if (t==0)
                //printMatrix(B);
            //printMatrix(A);

            // Update grid
            for (int i = 1; i < priceSteps - 1; i++) {
                grid[i][timesteps - 2 - t] = A_invBV[i];
            };

        };
    };

};

// Struct object for storing inputs from the  user
struct fdmInputs {
    float timeHorizon = 0;
    float strikePrice = 0;
    float riskFreeR = 0;
    float volatility = 0;
    float maxS = 0;
    float startPrice = 0;
};

// Function to prompt input from the user
fdmInputs promptInputs() {
    
    fdmInputs inputs;

    cout << "Enter time horizon (years): ";
    cin >> inputs.timeHorizon;

    cout << "Enter strike price: ";
    cin >> inputs.strikePrice;

    cout << "Enter risk-free interest rate (as a decimal, e.g., 0.05 for 5%): ";
    cin >> inputs.riskFreeR;

    cout << "Enter volatility (as a decimal, e.g., 0.2 for 20%): ";
    cin >> inputs.volatility;

    cout << "Enter starting price: ";
    cin >> inputs.startPrice;

    // Compute maxS based on the strike price
    inputs.maxS = 3 * inputs.strikePrice;

    return inputs;
}

int main()
{   

    

    //float timeHorizon = 1;          // One year
    //float strikePrice = 100;        // Strike price of the option
    //float startPrice = 100;         // Price when t=0
    //float riskFreeR = 0.05;         // Risk Free interest rate
    //float volatility = 0.2;         // Volatility in rand

    fdmInputs inputs = promptInputs();

    float timeHorizon = inputs.timeHorizon;          // One year
    float strikePrice = inputs.strikePrice;        // Strike price of the option
    float startPrice = inputs.startPrice;         // Price when t=0
    float riskFreeR = inputs.riskFreeR;         // Risk Free interest rate
    float volatility = inputs.volatility;         // Volatility in rand

    float dividends = 0;            // NOT USED

    float maxS = 3 * strikePrice;   // Should be three or four times the excersise price

    long int timesteps = 3000;      // Daily granularity 
    long int priceSteps = 180;      // Price granularity


    cout << "EXPLICIT PUT" << endl;
    FiniteDifferenceMethod fdm(timeHorizon, timesteps, maxS, priceSteps, strikePrice, volatility, riskFreeR, startPrice);
    fdm.explicitMethod(OptionType::PUT);
    //fdm.printGrid(fdm.priceSteps, fdm.timesteps);
    //FiniteDifferenceMethod::printMatrix(fdm.grid);
    //fdm.printGrid();
    fdm.getPrice();    
    //fdm.saveGridToCSV();

    cout << "EXPLICIT CALL" << endl;
    FiniteDifferenceMethod fdm2(timeHorizon, timesteps, maxS, priceSteps, strikePrice, volatility, riskFreeR, startPrice);
    fdm2.explicitMethod(OptionType::CALL);
    //FiniteDifferenceMethod::printMatrix(fdm2.grid);
    //fdm2.saveGridToCSV();
    //fdm2.printGrid(fdm2.priceSteps, fdm2.timesteps);
    fdm2.getPrice();
    

    cout << "IMPLICIT PUT" << endl;
    FiniteDifferenceMethod fdm3(timeHorizon, timesteps, maxS, priceSteps, strikePrice, volatility, riskFreeR, startPrice);
    fdm3.implicitMethod(OptionType::PUT);
    //fdm3.printGrid(fdm3.priceSteps, fdm3.timesteps);
    //FiniteDifferenceMethod::printMatrix(fdm3.grid);
    //fdm3.printGrid();
    fdm3.getPrice();    
    //fdm3.saveGridToCSV();

    cout << "IMPLICIT CALL" << endl;
    FiniteDifferenceMethod fdm4(timeHorizon, timesteps, maxS, priceSteps, strikePrice, volatility, riskFreeR, startPrice);
    fdm4.implicitMethod(OptionType::CALL);
    //FiniteDifferenceMethod::printMatrix(fdm4.grid);
    //fdm4.saveGridToCSV();
    //fdm4.printGrid(fdm4.priceSteps, fdm4.timesteps);
    fdm4.getPrice();


    cout << "CRANK PUT" << endl;
    FiniteDifferenceMethod fdm5(timeHorizon, timesteps, maxS, priceSteps, strikePrice, volatility, riskFreeR, startPrice);
    fdm5.crankNicolson(OptionType::PUT);
    //FiniteDifferenceMethod::printMatrix(fdm5.grid);
    //fdm5.printGrid();
    fdm5.getPrice();
    //fdm5.saveGridToCSV();

    cout << "CRANK CALL" << endl;
    FiniteDifferenceMethod fdm6(timeHorizon, timesteps, maxS, priceSteps, strikePrice, volatility, riskFreeR, startPrice);
    fdm6.crankNicolson(OptionType::CALL);
    //fdm6.saveGridToCSV();
    //fdm6.printGrid(fdm6.priceSteps, fdm6.timesteps);
    fdm6.getPrice();

    cout << "Process Complete" << endl;
    return 0;
}

