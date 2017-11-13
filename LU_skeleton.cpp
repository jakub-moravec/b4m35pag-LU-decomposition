#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <functional>
#include <vector>
#include <thread>
#include <cmath>
#include <atomic>

using namespace std;
using namespace std::chrono;

class LU {

	struct parameter {
		int i;
		int k;
		double stepHeight;
	} ;

	atomic<unsigned long> numberOfFinishedThreads{};
	atomic<int> runThread{};
    atomic<bool> programEnd{};

    unsigned long worker_count;
	std::vector<std::thread> thread_list;
	std::vector<parameter> parameters;

	public:

        ulong K = 0;

		explicit LU() {
			// spawn threads
			worker_count = thread::hardware_concurrency();
            numberOfFinishedThreads = worker_count;
            runThread = -1;
            programEnd = false;

            for (int i = 0; i < worker_count; ++i) {
                struct parameter dummy_par;
                parameters.push_back(dummy_par);
                thread_list.emplace_back(&LU::a_block, this, i);
            }
		}

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

                AA = A;
			} else {
				throw invalid_argument("Cannot open the input file!");
			}

			bin.close();
		}

    void u_row(int k) {
        for (int j = k + 1; j < K; j++) {
            U[k][j] = A[k][j];
        }
    }

    void l_col(int k) {
        L[k][k] = 1;
        for (int i = k + 1; i < K; i++) {
            L[i][k] = A[i][k] / U[k][k];
        }
    }

    void a_block(int thread_index) {
        do {
            // wait for bariere release
            while (thread_index > runThread) {
                // end thread
                if(programEnd.load()) {
                    return;
                }
                this_thread::sleep_for(chrono::milliseconds(1));
            }

            // fetch arguments
			struct parameter par = parameters[thread_index];

//            cout << "T: " << thread_index << ", K: " << par.k << ", I: " << par.i  << ", H: " << par.stepHeight << "\n";

            for (int i = par.i; i < par.i + par.stepHeight; i++) {
                for (int j = par.k + 1; j < K; j++) {
                    A[i][j] = A[i][j] - L[i][par.k] * U[par.k][j];
                }
            }


            // increment number of threads fetched by barier
            numberOfFinishedThreads++;

        } while (!programEnd.load());
	}
//TODO remove - Pavluv check code
	bool checkResult(){
		uint64_t n = A.size();
		bool correct = true;
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				double sum = 0;
				for (int k = 0; k < n; ++k) {
					sum += L[i][k] * U[k][j];
				}

				if(!double_equals(sum, AA[i][j])){
					cout << "Result of L*U [" << i <<", " << j << "] "   << sum << " should be " << AA[i][j] << endl;
					correct = false;
				}
			}
		}
		return correct;
	}

	bool double_equals(double a, double b, double epsilon = 0.001)
	{
		return std::fabs(a - b) < epsilon;
	}


	double decompose()	{
		high_resolution_clock::time_point start = high_resolution_clock::now();

		K = A.size();

		for (int k = 0; k < K; k++) {

			// barier->wait()
            while (numberOfFinishedThreads < worker_count) {
                this_thread::sleep_for(chrono::milliseconds(1));
            }

			U[k][k] = A[k][k];

			std::thread t1(&LU::u_row, this, k);
			std::thread t2(&LU::l_col, this, k);

			t1.join();
			t2.join();

            // recount A
			double step = ceil((A.size() - k) / (double) worker_count);
            int localThreadIndex = 0;
            for (int i = k + 1; i < K; i += step) {
                // předání argumentu vláknu
                struct parameter par;
                par.i = i;
                par.k = k;
                par.stepHeight = i + step > K ? K - i : step;
                parameters[localThreadIndex] = par;

                // release one thread
                ++runThread;
                --numberOfFinishedThreads;

                localThreadIndex++;
			}
            runThread = -1;
        }

        // join all threads
        programEnd = true;
        for (int i = 0; i < thread_list.size(); ++i) {
            thread_list[i].join();
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

		vector<vector<double>> A, L, U, AA; // FIXME remove AA
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

    out<<"Matrix AA:"<<endl;
    printMatrix(lu.AA);
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
	//TODO remove - Pavluv kod
	if(lu.checkResult()){
		cout<<"correct "<<endl;
	}

	return 0;
}

