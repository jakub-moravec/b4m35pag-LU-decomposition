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
#include <atomic>

using namespace std;
using namespace std::chrono;

class LU {

	atomic<unsigned long> numberOfFinishedThreads{};
	atomic<unsigned long> numberOfSpawnedThreads{};
    atomic<bool> parametersReaded{};

	int worker_count;
	std::vector<std::thread> thread_list;
    std::vector<bool> thread_readynnes;

	public:

        ulong K = 0;
        int k = 0;
        double step_width = 0;
        double step_height = 0;
        int I = 0;
        int J = 0;
        int threadIndex = 0;

		explicit LU() {
			// spawn threads
			worker_count = thread::hardware_concurrency();
            numberOfFinishedThreads = worker_count;
            numberOfSpawnedThreads = 0;
            parametersReaded = true;
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

    void a_block() {
        do {
            // fetch arguments
            int localI = I;
            int localJ = J;
            double local_step_widht = step_width;
            double local_step_height = step_height;
            int localThreadIndex = threadIndex;

            parametersReaded = true;

//            cout << "A thread: I: " << localI << ", J: " << localJ << ", T: " << localThreadIndex  << "\n";

            // lock thread
            thread_readynnes[localThreadIndex] = false;

            for (int i = localI; i < localI + local_step_height; i++) {
                for (int j = localJ; j < localJ + local_step_widht; j++) {
                    A[i][j] = A[i][j] - L[i][k] * U[k][j];
                }
            }


            // increment number of threads fetched by barier
            numberOfFinishedThreads++;

            // wait for bariere release
            while (!thread_readynnes[localThreadIndex]) {
                this_thread::sleep_for(chrono::milliseconds(1));
            }
        } while (k < K - 1);
	}


    double decompose()	{
		high_resolution_clock::time_point start = high_resolution_clock::now();

		K = A.size();

		for (k = 0; k < K; k++) {

			// barier->wait()
			while (numberOfFinishedThreads.load() < numberOfSpawnedThreads.load()) {
				this_thread::sleep_for(chrono::milliseconds(1)); // wait
                numberOfSpawnedThreads = 0;
			}

			U[k][k] = A[k][k];

			std::thread t1(&LU::u_row, this);
			std::thread t2(&LU::l_col, this);

			t1.join();
			t2.join();

			// release barier
			numberOfFinishedThreads = 0;

			// recount A
			double step = ceil(((A.size() - k) / std::sqrt(worker_count)));
            int localThreadIndex = 0;
			for (int localI = k + 1; localI < K; localI += step) {
				for (int localJ = k + 1; localJ < K; localJ += step) {
					double local_step_width = localJ + step > K ? K - localJ : step;
					double local_step_height = localI + step > K ? K - localI : step;


                    // kritická sekce // FIXME slow here
                    parametersReaded = false;
					I = localI;
					J = localJ;
                    step_width = local_step_width;
                    step_height = local_step_height;
					threadIndex = localThreadIndex;

                    if (k == 0) {
                        // spawn thread if needed
                        thread_readynnes.push_back(true);
                        thread_list.emplace_back(&LU::a_block, this);
                        ++numberOfSpawnedThreads;
                    } else {
                        // release one thread
                        thread_readynnes[localThreadIndex] = true;
                        ++numberOfSpawnedThreads;
                    }

                    while (!parametersReaded.load()) {
                        this_thread::sleep_for(chrono::milliseconds(1)); // wait
                    }

                    localThreadIndex++;
				}
			}
		}

        // todo join all threads

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

