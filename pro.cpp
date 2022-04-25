
/**
 * PRL - Project 2 - Preorder
 * @file pro.cpp
 * @author Ondřej Krejčí (xkrejc69)
 * @date 4/2022
 */

#include <vector>
#include <iostream> // cout, cerr
#include <limits>

#include <math.h> // ceil, log2

#include <mpi.h>

constexpr int PREORDER_TAG = 1;
constexpr int REQ_TAG = 2;
constexpr int WEIGHT_TAG = 3;
constexpr int SUCC_TAG = 4;

enum ERROR_CODE{
	ERR_ARGUMENTS = 1,
	ERR_COMMUNICATION
};

typedef struct{
	unsigned edge;
	unsigned reverse;
} adjacency;

constexpr int STOP = std::numeric_limits<int>::max();

// Abort the execution after an error has occured
void err_exit(const MPI_Comm& comm, const ERROR_CODE code){
	MPI_Abort(comm, code);
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
		std::cerr << "Requires one argument representing nodes of a binary tree as an array\n";
		err_exit(MPI_COMM_WORLD, ERR_ARGUMENTS);
	}
	return std::string(argv[1]);
}

std::vector<adjacency> buildAdjacencyList(const unsigned vertex, const size_t count){
	std::vector<adjacency> adjacency_list;

	// adding edges for the existing vertices
	if(vertex != 0){ // root does not have a parent
		adjacency_list.push_back({((vertex - 1) * 2) + 1, (vertex - 1) * 2});
	}
	if((2 * (vertex + 1)) <= count){ // left (<= because the index has a +1 offset)
		adjacency_list.push_back({4 * vertex, (4 * vertex) + 1});
	}
	if(((2 * (vertex + 1)) + 1) <= count){ // right
		adjacency_list.push_back({(4 * vertex) + 2, (4 * vertex) + 3});
	}
	return adjacency_list;
}

int main(int argc, char** argv){
	MPI_Init(&argc, &argv);

	const std::string nodes = parseArgs(argc, argv);
	const size_t count = nodes.length();
	const int commSize = getCommSize(MPI_COMM_WORLD);
	const int requiredProcs = (count == 1) ? 1 : (2 * count) - 2;
	if(requiredProcs != commSize){
		std::cerr << "Algorithm requires exaclty " << requiredProcs << " processes (used " << commSize << ")\n";
		err_exit(MPI_COMM_WORLD, ERR_ARGUMENTS);
	}

	const int rank = getCommRank(MPI_COMM_WORLD); // rank determines the rank of the edge
	const bool even = (rank % 2) == 0; // even index == forward edge

	const int to = even ? (rank / 2) + 1 : rank / 4; // index of vertex the edge goes to

	const std::vector<adjacency> adjacencyListTo = buildAdjacencyList(to, count); // relevant part of the adjacency list

	const unsigned reverseEdge = even ? rank + 1 : rank - 1; // index of reverse edge

	std::vector<int> etour(commSize);
	for(unsigned i = 0; i < adjacencyListTo.size(); i++){ // calculate next
		if(adjacencyListTo[i].edge == reverseEdge){
			if((i + 1) < adjacencyListTo.size()){ // next(e_r) <> nil
				etour[rank] = adjacencyListTo[i + 1].edge;
			}
			else{
				etour[rank] = adjacencyListTo[0].edge;
			}
		}
	}

	if(rank == 3){ // root correction
		etour[rank] = rank;
	}

	// prepare suffix sum variables
	int pred;
	int succ;
	unsigned weight = even ? 1 : 0;
	MPI_Request reqs[2];
	// predecessor index
	if(rank == 0){ // first edge has no predecessors
		pred = STOP;
	}
	else{
		MPI_Irecv(&pred, 1, MPI_INT, MPI_ANY_SOURCE, PREORDER_TAG, MPI_COMM_WORLD, &reqs[0]);
	}
	// successor index
	if(rank == 3){ // edge to root was removed, no need to send messages to itself
		succ = STOP;
	}
	else{
		succ = etour[rank];
		MPI_Isend(&rank, 1, MPI_INT, succ, PREORDER_TAG, MPI_COMM_WORLD, &reqs[1]);
	}
	if(rank == 0)
		MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
	else if(rank == 3)
		MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
	else
		MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

	unsigned succWeight;
	int succSucc;
	int predPred;

	if(succ == STOP){
		weight = 0;
	}
	for(unsigned i = 0; i < unsigned(ceil(log2(commSize))); i++){
		MPI_Request lreqs[6];
		unsigned reqCount = 0;

		if(succ != STOP){
			MPI_Isend(&pred, 1, MPI_INT, succ, REQ_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
		}

		if(pred != STOP){
			MPI_Irecv(&predPred, 1, MPI_INT, pred, REQ_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
			MPI_Isend(&weight, 1, MPI_UNSIGNED, pred, WEIGHT_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
			MPI_Isend(&succ, 1, MPI_INT, pred, SUCC_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
		}

		if(succ != STOP){
			MPI_Irecv(&succWeight, 1, MPI_UNSIGNED, succ, WEIGHT_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
			MPI_Irecv(&succSucc, 1, MPI_INT, succ, SUCC_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
		}

		MPI_Waitall(reqCount, lreqs, MPI_STATUSES_IGNORE);

		if(succ != STOP){
			weight += succWeight;
			succ = succSucc;

		}
		if(pred != STOP){
			pred = predPred;
		}
	}
	// no need for correction, last edge always has a weight of 0

	unsigned reqCount = 0;
	std::vector<unsigned> preorder(commSize);
	if(even){ // send the computed weight to the correct index
		const unsigned result = count - weight;
		MPI_Isend(&result, 1, MPI_UNSIGNED, to, PREORDER_TAG, MPI_COMM_WORLD, &reqs[reqCount++]);
	}
	if(rank < static_cast<int>(count) && rank != 0){ // receive the computed weight
		MPI_Irecv(&preorder[rank], 1, MPI_UNSIGNED, MPI_ANY_SOURCE, PREORDER_TAG, MPI_COMM_WORLD, &reqs[reqCount++]);
	}
	if(reqCount > 0){
		MPI_Waitall(reqCount, reqs, MPI_STATUSES_IGNORE);
	}

	if(rank == 0){
		preorder[0] = 0;
	}
	
	// gather the preorder traversal indices to root
	MPI_Gather(&preorder[rank], 1, MPI_UNSIGNED, preorder.data(), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	
	if(rank == 0){
		std::string out;
		out.resize(count + 1);
		for(unsigned i = 0; i < count; i++){
			out[preorder[i]] = nodes[i];
		}
		std::cout << out << std::endl;
	}

	MPI_Finalize();
	return 0;
}
