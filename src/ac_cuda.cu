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

// TODO: to void
static TableGpu* copyTableToGpu(TabNode* table, TableGpu** hostCopy);
static void deleteTableGpu(TableGpu* table, TableGpu* hostCopy);

static void chainGpuCreate(Chain* chain, ChainGpu** chainD, ChainGpu** chainH);
static void chainGpuDelete(ChainGpu* chainD, ChainGpu* chainH);

__global__ static void findCandidates(TableGpu** automata, 
    int automataLen, int* candidates);

// ***************************************************************************
// PUBLIC
extern void* indicesTableCreateGpu(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer) {

    vector<TabNode*>* aut = static_cast<vector<TabNode*>*>(automata);

    vector<TableGpu*> gpuTables;
    vector<TableGpu*> hostTables;

    gpuTables.reserve(automataLen);
    hostTables.reserve(automataLen);

    vector<ChainGpu*> gpuChains;
    vector<ChainGpu*> hostChains;

    int numTargets = databaseLen - databaseStart;

    gpuChains.reserve(numTargets);
    hostChains.reserve(numTargets);

    //**************************************************************************
    // SEND AUTOMATA TO GPU
    
    for (int i = 0; i < automataLen; ++i) {
        TabNode* autH = (*aut)[i];

        TableGpu* hostCopy;

        gpuTables.push_back(copyTableToGpu(autH, &hostCopy));
        hostTables.push_back(hostCopy);
        // printf("%p\n", gpuTables[i]);
    }

    TableGpu** gpuTablesD;
    CUDA_SAFE_CALL(cudaMalloc(&gpuTablesD, sizeof(TableGpu*) * automataLen));
    CUDA_SAFE_CALL(cudaMemcpy(gpuTablesD, &gpuTables[0], 
        sizeof(TableGpu*) * automataLen, TO_GPU));    

    //**************************************************************************

    //**************************************************************************
    // SEND DATABASE TO GPU

    for (int i = databaseStart; i < databaseLen; ++i) {
        ChainGpu* chainD;
        ChainGpu* chainH;

        chainGpuCreate(database[i], &chainD, &chainH);

        gpuChains.push_back(chainD);
        hostChains.push_back(chainH);
    }

    ChainGpu** chainsGpuD;
    CUDA_SAFE_CALL(cudaMalloc(&chainsGpuD, sizeof(ChainGpu*) * numTargets));
    CUDA_SAFE_CALL(cudaMemcpy(chainsGpuD, &gpuChains[0], 
        sizeof(ChainGpu*) * numTargets, TO_GPU));

    //**************************************************************************


    dim3 dimGrid(1,1,1);
    dim3 dimBlock(1,1,1);

    int* candidatesD; 
    CUDA_SAFE_CALL(cudaMalloc(&candidatesD, sizeof(int) * 5001 * automataLen));

    findCandidates<<<dimGrid, dimBlock>>>(gpuTablesD, automataLen, candidatesD);

    // clean up
    for (int i = 0; i < automataLen; ++i) {
        deleteTableGpu(gpuTables[i], hostTables[i]);
    }

    for (int i = 0; i < numTargets; ++i) {
        chainGpuDelete(gpuChains[i], hostChains[i]);
    }

    CUDA_SAFE_CALL(cudaFree(gpuTablesD));
    CUDA_SAFE_CALL(cudaFree(chainsGpuD));
    CUDA_SAFE_CALL(cudaFree(candidatesD));

    gpuTables.clear();
    return NULL;
}
// ***************************************************************************

// ***************************************************************************
// PRIVATE
static TableGpu* copyTableToGpu(TabNode* table, TableGpu** hostCopy) {
    *hostCopy = (TableGpu*) malloc(sizeof(TableGpu));

    (*hostCopy)->numStates = table->numStates;
    (*hostCopy)->table = table->table;

    // flatten and copy positions vector
    int start = 0;
    vector<int> positions;
    vector<vector<uint16> > &v = table->positions;

    for (int i = 0; i < (*hostCopy)->numStates; ++i) {
        (*hostCopy)->table[i * TABLE_WIDTH + POSITIONS_START] = start;

        positions.insert(positions.end(), v[i].begin(), v[i].end());

        start += v[i].size();
    }

    uint16* positionsD;
    CUDA_SAFE_CALL(cudaMalloc(&positionsD, sizeof(uint16) * positions.size()));

    CUDA_SAFE_CALL(cudaMemcpy(positionsD, &positions[0], 
        sizeof(uint16) * positions.size(),
        cudaMemcpyHostToDevice));

    (*hostCopy)->positions = positionsD;

    // state table
    int* statesD;
    CUDA_SAFE_CALL(cudaMalloc(&statesD, sizeof(int) * (*hostCopy)->numStates * TABLE_WIDTH));

    CUDA_SAFE_CALL(cudaMemcpy(statesD, (*hostCopy)->table,
        sizeof(int) * (*hostCopy)->numStates * TABLE_WIDTH, 
        cudaMemcpyHostToDevice));

    (*hostCopy)->table = statesD;

    TableGpu* autD;
    CUDA_SAFE_CALL(cudaMalloc(&autD, sizeof(TableGpu)));

    CUDA_SAFE_CALL(cudaMemcpy(autD, (*hostCopy), sizeof(TableGpu), 
        cudaMemcpyHostToDevice));

    // printf("autd %p\n", autD);
    return autD;
}

static void deleteTableGpu(TableGpu* table, TableGpu* hostCopy) {
    CUDA_SAFE_CALL(cudaFree(hostCopy->table));
    CUDA_SAFE_CALL(cudaFree(hostCopy->positions));
    CUDA_SAFE_CALL(cudaFree(table));    

    free(hostCopy);
}

static void chainGpuCreate(Chain* chain, ChainGpu** chainD, ChainGpu** chainH) {
    (*chainH) = (ChainGpu*) malloc(sizeof(ChainGpu));

    int len = chainGetLength(chain);
    const char* codesH = chainGetCodes(chain);
    char* codesD;

    CUDA_SAFE_CALL(cudaMalloc(&codesD, sizeof(char) * len));
    CUDA_SAFE_CALL(cudaMemcpy(codesD, codesH, sizeof(char) * len, TO_GPU));

    (*chainH)->len = len;
    (*chainH)->codes = codesD;

    CUDA_SAFE_CALL(cudaMalloc(chainD, sizeof(ChainGpu)));
    CUDA_SAFE_CALL(cudaMemcpy(*chainD, *chainH, sizeof(ChainGpu), TO_GPU));
}

static void chainGpuDelete(ChainGpu* chainD, ChainGpu* chainH) {
    cudaFree(chainH->codes);
    cudaFree(chainD);

    free(chainH);
}
// ***************************************************************************

// ***************************************************************************
// GPU Modules

__global__ static void findCandidates(TableGpu** automata, 
    int automataLen, int* candidates) {

    printf("heklo worls\n");
    return;
}

// ***************************************************************************
