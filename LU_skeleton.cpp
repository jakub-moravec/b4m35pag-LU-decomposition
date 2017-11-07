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

        ulong K;
        int k = 0;
        int step_width = 0;
        int step_height = 0;
        int I = 0;
        int J = 0;

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

    void u_row() {
        for (int j = k + 1; j < K; j++) {
            U[k][j] = A[k][j];
        }
    }

    void l_col() {
        L[k][k] = 1;
        for (int i = k + 1; i < K; i++) {
            L[i][k] = A[i][k] / U[k][k];
        }
    }

    /**
     *
     * @param k
     * @param I start i index
     * @param J start j index
     */
    void a_block() { // work
        int localI = I;
        int localJ = J;
        for (int i = localI; i < localI + step_height; i++) {
            for (int j = localJ; j < localJ + step_width; j++) {
                A[i][j] = A[i][j] - L[i][k] * U[k][j];
            }
        }
    }


    double decompose()	{
			high_resolution_clock::time_point start = high_resolution_clock::now();

            K = A.size();

			int worker_count = thread::hardware_concurrency();
			cout << "Worker count = " << worker_count << endl;

            for (k = 0; k < K; k++) {

                U[k][k] = A[k][k];

                std::thread t1(&LU::u_row, this);
                std::thread t2(&LU::l_col, this);

                t1.join();
                t2.join();

                std::vector<std::thread> thread_list;
                int step = ceil(((A.size() - k) / std::sqrt(worker_count)));
                for (int localI = k + 1; localI < K; localI += step) {
                    for (int localJ = k + 1; localJ < K; localJ += step) {
                        step_width = localJ + step > K ? K - localJ : step;
                        step_height = localI + step > K? K - localI : step;

                        I = localI; // fixme kritická sekce
                        J = localJ; // fixme kritická sekce
                        thread_list.push_back(std::thread(&LU::a_block, this));
                    }
                }

                // todo barier
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

