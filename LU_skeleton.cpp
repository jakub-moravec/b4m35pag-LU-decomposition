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
		int j;
		int k;
		double stepWidth;
		double stepHeight;
	} ;

	atomic<unsigned long> numberOfFinishedThreads{};
	atomic<unsigned long> numberOfSpawnedThreads{};
    atomic<bool> programEnd{};

	int worker_count;
	std::vector<std::thread> thread_list;
    std::vector<bool> thread_readynnes;
	std::vector<parameter> parameters;

	public:

        ulong K = 0;

		explicit LU() {
			// spawn threads
			worker_count = thread::hardware_concurrency() - 1;
            numberOfFinishedThreads = worker_count;
            numberOfSpawnedThreads = 0;
            programEnd = false;
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
            // fetch arguments
			struct parameter par = parameters[thread_index];

            cout << "T: " << thread_index << ", K: " << par.k << ", I: " << par.i << ", J: " << par.j << ", W: " << par.stepWidth << ", H: " << par.stepHeight << "\n";

            // lock thread
            thread_readynnes[thread_index] = false;

            for (int i = par.i; i < par.i + par.stepHeight; i++) {
                for (int j = par.j; j < par.j + par.stepWidth; j++) {
                    A[i][j] = A[i][j] - L[i][par.k] * U[par.k][j];
                }
            }


            // increment number of threads fetched by barier
            numberOfFinishedThreads++;

            // wait for bariere release
            while (!thread_readynnes[thread_index]) {
                // end thread
                if(programEnd.load()) {
                    return;
                }
                this_thread::sleep_for(chrono::milliseconds(1));
            }

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
//					cout << "Result of L*U [" << i <<", " << j << "] "   << sum << " should be " << AA[i][j] << endl;
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
			while (numberOfFinishedThreads.load() < numberOfSpawnedThreads.load()) {
				this_thread::sleep_for(chrono::milliseconds(1)); // wait
                numberOfSpawnedThreads = 0;
			}

			U[k][k] = A[k][k];

			std::thread t1(&LU::u_row, this, k);
			std::thread t2(&LU::l_col, this, k);

			t1.join();
			t2.join();

			// release barier
			numberOfFinishedThreads = 0;

			// recount A
			double step = ceil(((A.size() - k) / std::sqrt(worker_count)));
            int localThreadIndex = 0;
			for (int i = k + 1; i < K; i += step) {
				for (int j = k + 1; j < K; j += step) {

                    // předání argumentu vláknu
					struct parameter par;
					par.i = i;
					par.j = j;
					par.k = k;
					par.stepWidth = j + step > K ? K - j : step;
					par.stepHeight = i + step > K ? K - i : step;

                    if (k == 0) {
                        // spawn thread if needed
                        thread_readynnes.push_back(true);
						parameters.push_back(par);
                        thread_list.emplace_back(&LU::a_block, this, localThreadIndex);
                        ++numberOfSpawnedThreads;
                    } else {
                        // release one thread
                        parameters[localThreadIndex] = par;
                        thread_readynnes[localThreadIndex] = true;
                        ++numberOfSpawnedThreads;
                    }

                    localThreadIndex++;
				}
			}
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

