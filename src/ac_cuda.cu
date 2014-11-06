#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <queue>

using namespace std;

#include "swsharp/swsharp.h"
#include "swsharp/cuda_utils.h"
#include "table_node.h"
#include "ac_table.h"
#include "ac_cuda_utils.h"

// ***************************************************************************
// PUBLIC
extern void* indicesTableCreateGpu(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer);
// ***************************************************************************

// ***************************************************************************
// PRIVATE
typedef vector<vector<int> > Candidates;
typedef vector<int> Candidate;

static TableGpu* copyTableToGpu(TabNode* table);
static void deleteTableGpu(TableGpu* table);

__global__ static void findCandidates(TableGpu* automata, 
    int automataLen, int* candidates);

// ***************************************************************************
// PUBLIC
extern void* indicesTableCreateGpu(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer) {

    vector<TabNode*>* aut = static_cast<vector<TabNode*>*>(automata);
    vector<TableGpu*> gpuTables;
    gpuTables.reserve(automataLen);

    for (int i = 0; i < automataLen; ++i) {
        TabNode* autH = (*aut)[i];
        gpuTables.push_back(copyTableToGpu(autH));

    }

    dim3 dimGrid(1,1,1);
    dim3 dimBlock(1,1,1);

    TableGpu* gpuTablesD;
    CUDA_SAFE_CALL(cudaMalloc(&gpuTablesD, sizeof(TableGpu*) * automataLen));
    CUDA_SAFE_CALL(cudaMemcpy(gpuTablesD, &gpuTables[0], 
        sizeof(TableGpu*) * automataLen, 
        cudaMemcpyHostToDevice));

    int* candidatesD; 
    CUDA_SAFE_CALL(cudaMalloc(&candidatesD, sizeof(int) * 5001 * automataLen));
    findCandidates<<<dimGrid, dimBlock>>>(gpuTablesD, automataLen, candidatesD);

    // clean up
    CUDA_SAFE_CALL(cudaFree(gpuTablesD));
    for (int i = 0; i < automataLen; ++i) {
        deleteTableGpu(gpuTables[i]);
    }
    gpuTables.clear();
    return NULL;
}
// ***************************************************************************

// ***************************************************************************
// PRIVATE

static TableGpu* copyTableToGpu(TabNode* table) {
    TableGpu* copyAut = (TableGpu*) malloc(sizeof(TableGpu));

    copyAut->numStates = table->numStates;
    copyAut->table = table->table;

    // flatten and copy positions vector
    int start = 0;
    vector<int> positions;
    vector<vector<uint16> > &v = table->positions;

    for (int i = 0; i < copyAut->numStates; ++i) {
        copyAut->table[i * TABLE_WIDTH + POSITIONS_START] = start;

        positions.insert(positions.end(), v[i].begin(), v[i].end());

        start += v[i].size();
    }

    uint16* positionsD;
    CUDA_SAFE_CALL(cudaMalloc(&positionsD, sizeof(uint16) * positions.size()));

    CUDA_SAFE_CALL(cudaMemcpy(positionsD, &positions[0], 
        sizeof(uint16) * positions.size(),
        cudaMemcpyHostToDevice));

    copyAut->positions = positionsD;

    // state table
    int* statesD;
    CUDA_SAFE_CALL(cudaMalloc(&statesD, sizeof(int) * copyAut->numStates * TABLE_WIDTH));

    CUDA_SAFE_CALL(cudaMemcpy(statesD, copyAut->table,
        sizeof(int) * copyAut->numStates * TABLE_WIDTH, 
        cudaMemcpyHostToDevice));

    copyAut->table = statesD;

    TableGpu* autD;
    CUDA_SAFE_CALL(cudaMalloc(&autD, sizeof(TableGpu)));

    CUDA_SAFE_CALL(cudaMemcpy(autD, copyAut, sizeof(TableGpu), 
        cudaMemcpyHostToDevice));

    delete copyAut;
    return autD;
}

static void deleteTableGpu(TableGpu* table) {
    CUDA_SAFE_CALL(cudaFree(table->table));
    CUDA_SAFE_CALL(cudaFree(table->positions));
    CUDA_SAFE_CALL(cudaFree(table));    
}

// ***************************************************************************

// ***************************************************************************
// GPU Modules

__global__ static void findCandidates(TableGpu* automata, 
    int automataLen, int* candidates) {

    printf("heklo worls\n");
    return;
}