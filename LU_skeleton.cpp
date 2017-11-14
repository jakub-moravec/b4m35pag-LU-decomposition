#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <string>
#include <functional>
#include <vector>
#include <thread>
#include <cmath>

using namespace std;
using namespace std::chrono;

class LU {
	public:
		// It reads a matrix from the binary file.
		void readMatrixFromInputFile(const string& inputFile)	{
			ifstream bin(inputFile.c_str(), ifstream::in | ifstream::binary);
			if (bin)	{
				uint64_t n = 0;
				bin.read((char*) &n, sizeof(uint64_t));
				A.resize(n, vector<double>(n, 0.0));
				L = U = A;
				for (uint64_t r = 0; r < n; ++r)
					bin.read((char*) A[r].data(), n*sizeof(double));
			} else {
				throw invalid_argument("Cannot open the input file!");
			}

			bin.close();
		}

    void u_row(int k) {
        for (int j = k+1; j <  A.size(); j++) {
            U[j][k] = A[k][j];
        }
    }

    void l_col(int k) {
        L[k][k] = 1;
        for (int i = k+1; i < A.size(); i++) {
            L[i][k] = A[i][k] / U[k][k];
        }
    }

    void a_block(int k, int start_row, int height) {
        for (int i = start_row; i < start_row + height ; i++) {
            for (int j = k + 1; j < A.size(); j++) {
                A[i][j] = A[i][j] - L[i][k] * U[j][k];
            }
        }
    }


    double decompose()	{
        high_resolution_clock::time_point start = high_resolution_clock::now();

        int K = (int) A.size();
        std::vector<std::thread> thread_list;

        // paralel
        int worker_count = 8;
        cout << "Threads = " << worker_count << endl;

        for (int k = 0; k < K; k++) {

            int step = (int) ceil((K - k) / (double) worker_count);

            for (auto &t : thread_list) {
                t.join();
            }
            thread_list.clear();

            // LU
            U[k][k] = A[k][k];
            int lu_step = step * 2;
            thread_list.emplace_back(&LU::u_row, this, k);
            thread_list.emplace_back(&LU::l_col, this, k);

            for (auto &t : thread_list) {
                t.join();
            }
            thread_list.clear();

            // A
            for (int i = k + 1; i < K; i += step) {
                int height =i + step > K ? K - i : step;
                thread_list.emplace_back(&LU::a_block, this, k, i, height);
            }
        }

        for (auto &t : thread_list) {
            t.join();
        }
        thread_list.clear();


        double runtime = duration_cast<duration<double>>(high_resolution_clock::now()-start).count();

        return runtime;
    }

        void writeResults(const string& outputFile)	{
            ofstream bout(outputFile.c_str(), ofstream::out | ofstream::binary | ofstream::trunc);
            if (bout)	{
                uint64_t n = A.size();
                for (uint64_t r = 0; r < n; ++r)
                    bout.write((char*) L[r].data(), n*sizeof(double));
                for (uint64_t r = 0; r < n; ++r) {
                    std::vector<double> column(n);
                    for (uint64_t q = 0; q < n; ++q) {
                        column[q] = U[q][r];
                    }
                    bout.write((char *) column.data(), n * sizeof(double));
                }
            } else {
                throw invalid_argument("Cannot open the input file!");
            }

            bout.close();
        }
		
	private:

		vector<vector<double>> A, L, U;
		friend ostream& operator<<(ostream&, const LU&);
};

// Print the matrices A, L, and U in an instance of LU class.
ostream& operator<<(ostream& out, const LU& lu)	{
	function<void(const vector<vector<double>>&)> printMatrix = [&](const vector<vector<double>>& M)	{
		int n = M.size();
		for (int i = 0; i < n; ++i)	{
			for (int j = 0; j < n; ++j)
				out<<" "<<setw(10)<<M[i][j];
			out<<endl;
		}
	};

	out<<"Matrix A:"<<endl;
	printMatrix(lu.A);
	out<<endl<<"Lower matrix:"<<endl;
	printMatrix(lu.L);
	out<<endl<<"Upper matrix:"<<endl;
	printMatrix(lu.U);

	return out;
}

int main(int argc, char* argv[])	{
	if (argc <= 1 || argc > 3)	{
		cout<<"LU decomposition of a square matrix."<<endl;
		cout<<endl<<"Usage:"<<endl;
		cout<<"\t"<<argv[0]<<" inputMatrix.bin [output.bin]"<<endl;
		return 1;
	}

	string inputFile = argv[1], outputFile;
	if (argc == 3)
		outputFile = argv[2];

	LU lu;
	lu.readMatrixFromInputFile(inputFile);
	double totalDuration = lu.decompose();
	// Decomposition is printed only if the output file is not written.
	if (outputFile.empty())
		cout<<lu<<endl;
	else
		lu.writeResults(outputFile);

	cout<<"computational time: "<<totalDuration<<" s"<<endl;

	return 0;
}

