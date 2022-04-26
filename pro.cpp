
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

constexpr int PRED_TAG = 1;
constexpr int WEIGHT_TAG = 2;
constexpr int SUCC_TAG = 3;
constexpr int PREORDER_TAG = 4;

constexpr int ERR_ARGUMENTS = 1;

typedef struct{
	unsigned edge;
	unsigned reverse;
} adjacency;

constexpr int STOP = std::numeric_limits<int>::max();

// Abort the execution after an error has occured
void err_exit(const MPI_Comm& comm, const int code){
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

// Checks that the required argument exists and returns it
std::string parseArgs(int argc, char** argv, const int rank){
	if(argc < 2){
		if(rank == 0){
			std::cerr << "Requires one argument representing nodes of a binary tree as an array\n";
		}
		err_exit(MPI_COMM_WORLD, ERR_ARGUMENTS);
	}
	return std::string(argv[1]);
}

// Builds part of the adjancency list for the given vertex
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

// Calculates value of the Euler tour for the given edge (parameter is already just the reverse edge)
int getEtourPart(const unsigned reverseEdge, const std::vector<adjacency>& adjacencyListTo, const int rank){
	int etour;
	for(unsigned i = 0; i < adjacencyListTo.size(); i++){ // calculate next
		if(adjacencyListTo[i].edge == reverseEdge){
			if((i + 1) < adjacencyListTo.size()){ // next(e_r) <> nil
				etour = adjacencyListTo[i + 1].edge;
			}
			else{
				etour = adjacencyListTo[0].edge;
			}
		}
	}

	if(etour == 0){ // root correction
		etour = rank;
	}

	return etour;
}

// Calculates suffix sum (modified/simplified for use with the preorder algorithm) and returns the weight for the edge
unsigned preorderSuffixSum(const int rank, const bool even, const int etour, const int commSize){
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
		MPI_Irecv(&pred, 1, MPI_INT, MPI_ANY_SOURCE, PRED_TAG, MPI_COMM_WORLD, &reqs[0]);
	}

	// successor index
	if(etour == rank){ // edge to root was removed, no need to send messages to itself
		succ = STOP;
		weight = 0;
	}
	else{
		succ = etour;
		MPI_Isend(&rank, 1, MPI_INT, succ, PRED_TAG, MPI_COMM_WORLD, &reqs[1]); // informs the successor about its predecessor
	}

	// exchange info about initial predecessors
	if(rank == 0)
		MPI_Wait(&reqs[1], MPI_STATUS_IGNORE);
	else if(etour == rank)
		MPI_Wait(&reqs[0], MPI_STATUS_IGNORE);
	else
		MPI_Waitall(2, reqs, MPI_STATUSES_IGNORE);

	unsigned succWeight; // weight of the successor
	int succSucc; // successor of the successor
	int predPred; // predecessor of the predecessor
	for(unsigned i = 0; i < unsigned(ceil(log2(commSize))); i++){ // internal suffix sum loop
		MPI_Request lreqs[6];
		unsigned reqCount = 0;

		if(pred != STOP){
			MPI_Irecv(&predPred, 1, MPI_INT, pred, PRED_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
			MPI_Isend(&weight, 1, MPI_UNSIGNED, pred, WEIGHT_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
			MPI_Isend(&succ, 1, MPI_INT, pred, SUCC_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
		}

		if(succ != STOP){
			MPI_Isend(&pred, 1, MPI_INT, succ, PRED_TAG, MPI_COMM_WORLD, &lreqs[reqCount++]);
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

		MPI_Barrier(MPI_COMM_WORLD);
	}

	// no need for correction, last edge always has a weight of 0 (is a reverse edge)

	return weight;
}

// Sends preorde indexes to the appropriate process, then gathers them at the root
std::vector<unsigned> exchangePreorder(const int rank, const bool even, const int commSize, 
	const size_t nodeCount, const unsigned weight, const int to
){
	MPI_Request reqs[2];
	unsigned reqCount = 0;
	std::vector<unsigned> preorder(commSize);

	if(even){ // send the computed weight to the correct index
		const unsigned result = nodeCount - weight;
		MPI_Isend(&result, 1, MPI_UNSIGNED, to, PREORDER_TAG, MPI_COMM_WORLD, &reqs[reqCount++]);
	}

	if(rank < static_cast<int>(nodeCount) && rank != 0){ // receive the computed weight
		MPI_Irecv(&preorder[rank], 1, MPI_UNSIGNED, MPI_ANY_SOURCE, PREORDER_TAG, MPI_COMM_WORLD, &reqs[reqCount++]);
	}

	if(reqCount > 0){
		MPI_Waitall(reqCount, reqs, MPI_STATUSES_IGNORE);
	}

	if(rank == 0){
		preorder[0] = 0;
	}
	
	// gather the preorder traversal indices at the root
	MPI_Gather(&preorder[rank], 1, MPI_UNSIGNED, preorder.data(), 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);

	return preorder;
}

// Print the tree's nodes in preorder
void preorderPrint(const std::string& nodes, const std::vector<unsigned>& preorder){
	const size_t count = nodes.size();
	std::string out;
	out.resize(count + 1);
	for(unsigned i = 0; i < count; i++){
		out[preorder[i]] = nodes[i];
	}
	std::cout << out << std::endl;
}

// Calculate preorder traversal and print the tree's nodes
int main(int argc, char** argv){
	MPI_Init(&argc, &argv);

	const int rank = getCommRank(MPI_COMM_WORLD); // rank determines the rank of the edge
	const int commSize = getCommSize(MPI_COMM_WORLD);
	const std::string nodes = parseArgs(argc, argv, rank);
	const size_t count = nodes.length();
	const int requiredProcs = (count == 1) ? 1 : (2 * count) - 2;
	if(requiredProcs != commSize){
		if(rank == 0){
			std::cerr << "Algorithm requires exaclty " << requiredProcs << " processes (used " << commSize << ")\n";
		}
		err_exit(MPI_COMM_WORLD, ERR_ARGUMENTS);
	}

	if(count == 1){ // handle root-only tree as a special case
		std::cout << nodes[0] << std::endl;
		MPI_Finalize();
		return 0;
	}

	const bool even = (rank % 2) == 0; // even index == forward edge

	const int to = even ? (rank / 2) + 1 : rank / 4; // index of vertex the edge goes to

	const std::vector<adjacency> adjacencyListTo = buildAdjacencyList(to, count); // relevant part of the adjacency list
	const unsigned reverseEdge = even ? rank + 1 : rank - 1; // index of the reverse edge
	const int etour = getEtourPart(reverseEdge, adjacencyListTo, rank); // Etour value for the edge

	const unsigned weight = preorderSuffixSum(rank, even, etour, commSize);

	const std::vector<unsigned> preorder = exchangePreorder(rank, even, commSize, count, weight, to);
	
	if(rank == 0){
		preorderPrint(nodes, preorder);
	}

	MPI_Finalize();
	return 0;
}
