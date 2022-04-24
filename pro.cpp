
/**
 * PRL - Project 2 - Preorder
 * @file pro.cpp
 * @author Ondřej Krejčí (xkrejc69)
 * @date 4/2022
 */

#include <array>
#include <map>
#include <iostream> // cout, cerr

#include <mpi.h>

enum ERROR_CODE{
	ERR_ARGUMENTS = 1,
	ERR_COMMUNICATION
};

// Abort the execution after an error has occured
void err_exit(const MPI_Comm& comm, const ERROR_CODE code){
	MPI_Abort(comm, code);
}

// Prints error for the MPI error code
void printError(const int ec){
	char estring[MPI_MAX_ERROR_STRING];
	int len;
	MPI_Error_string(ec, estring, &len);
	std::cerr << estring << std::endl;
}

// Retrieves the rank of a processor
int getCommRank(const MPI_Comm& comm){
	int commRank;
	MPI_Comm_rank(comm, &commRank);
	return commRank;
}

// Retrieves the size of a communicator
int getCommSize(const MPI_Comm& comm){
	int commSize;
	MPI_Comm_size(comm, &commSize);
	return commSize;
}

std::string parseArgs(int argc, char** argv){
	if(argc < 2){
		std::cerr << "Requires one argument representing nodes of a binary tree as an array";
		exit(ERR_ARGUMENTS);
	}
	return std::string(argv[1]);
}

int main(int argc, char** argv){
	MPI_Init(&argc, &argv);

	const int rank = getCommRank(MPI_COMM_WORLD);
	const std::string nodes = parseArgs(argc, argv);

	MPI_Finalize();

	return 0;
}
