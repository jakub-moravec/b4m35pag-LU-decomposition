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
        for (int j = k + 1; j < A.size(); j++) {
            U[k][j] = A[k][j];
        }
    }

    void l_col(int k) {
        L[k][k] = 1;
        for (int i = k + 1; i < A.size(); i++) {
            L[i][k] = A[i][k] / U[k][k];
        }
    }

    void a_block(int k, int start_row, int start_col, int width, int height) { // work
        for (int m = start_row; m < start_row + height ; m++) {
            for (int n = start_col; n < start_col + width ; n++) {
                A[m][n] = A[m][n] - L[m][k] * U[k][n];
            }
        }
    }


    double decompose()	{
			high_resolution_clock::time_point start = high_resolution_clock::now();

			// paralel
            int K = A.size();
			int worker_count = thread::hardware_concurrency();
			cout << "Worker count = " << worker_count << endl;
            for (int k = 0; k < K; k++) {
                U[k][k] = A[k][k];

                std::thread t1(&LU::u_row, this, k);
                std::thread t2(&LU::l_col, this, k);

                t1.join();
                t2.join();

                std::vector<std::thread> thread_list;
                int step = ceil(((A.size() - k) / std::sqrt(worker_count)));
                for (int a = k + 1; a < K; a += step) {
                    for (int b = k + 1; b < K; b += step) {
                        int tmp_width = step;
                        int tmp_height = step;
                        if (a + step > K) {
                            tmp_height = K - a;
                        }
                        if (b + step > K) {
                            tmp_width = K - b;
                        }
                        thread_list.emplace_back(&LU::a_block, this, k, a, b, tmp_width, tmp_height);
                    }
                }

                for (int i = 0; i < thread_list.size(); i++) {
                    thread_list[i].join();
                }
                thread_list.clear();
            }


			double runtime = duration_cast<duration<double>>(high_resolution_clock::now()-start).count();

			return runtime;
		}

		void writeResults(const string& outputFile)	{
			ofstream bout(outputFile.c_str(), ofstream::out | ofstream::binary | ofstream::trunc);
			if (bout)	{
				uint64_t n = A.size();
				for (uint64_t r = 0; r < n; ++r)
					bout.write((char*) L[r].data(), n*sizeof(double));
				for (uint64_t r = 0; r < n; ++r)
					bout.write((char*) U[r].data(), n*sizeof(double));
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

