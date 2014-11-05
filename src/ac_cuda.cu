#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <queue>

using namespace std;

#include "swsharp/swsharp.h"
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
static TableGpu* copyTableToGpu(TabNode* table);
static void deleteTableGpu(TableGpu* table);

__global__ static void findCandidates(TableGpu* automata, int automataLen);

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
    cudaMalloc(&gpuTablesD, sizeof(TableGpu*) * automataLen);
    cudaMemcpy(gpuTablesD, &gpuTables[0], 
        sizeof(TableGpu*) * automataLen, 
        cudaMemcpyHostToDevice);

    findCandidates<<<dimGrid, dimBlock>>>(gpuTablesD, automataLen);

    // clean up
    cudaFree(gpuTablesD);
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

    int* states;
    cudaMalloc(&states, sizeof(int) * copyAut->numStates);
    cudaMemcpy(states, copyAut->table, sizeof(int) * copyAut->numStates,
        cudaMemcpyHostToDevice);

    copyAut->table = states;

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
    cudaMalloc(&positionsD, sizeof(uint16) * positions.size());
    cudaMemcpy(positionsD, &positions[0], 
        sizeof(uint16) * positions.size(),
        cudaMemcpyHostToDevice);

    copyAut->positions = positionsD;

    TableGpu* autD;
    cudaMalloc(&autD, sizeof(TableGpu));
    cudaMemcpy(autD, copyAut, sizeof(TableGpu), 
        cudaMemcpyHostToDevice);

    delete copyAut;
    return autD;
}

static void deleteTableGpu(TableGpu* table) {
    cudaFree(table->table);
    cudaFree(table->positions);
    cudaFree(table);    
}

// ***************************************************************************

// ***************************************************************************
// GPU Modules

__global__ static void findCandidates(TableGpu* automata, int automataLen) {
    return;
}

// ***************************************************************************
