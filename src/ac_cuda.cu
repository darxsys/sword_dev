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
    int automataLen, ChainGpu** database, int databaseLen,
    int* candidates);

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
    int* candidatesD;
    CUDA_SAFE_CALL(cudaMalloc(&candidatesD, sizeof(int) * 5001 * automataLen));
    int* candidatesH = (int*) malloc(sizeof(int) * 5001 * automataLen);


    dim3 dimGrid(1,1,1);
    // TODO: change
    dim3 dimBlock(automataLen,1,1);

    findCandidates<<<dimGrid, dimBlock>>>(gpuTablesD, automataLen, chainsGpuD, 
        numTargets, candidatesD);

    CUDA_SAFE_CALL(cudaMemcpy(candidatesH, candidatesD, 
        sizeof(int) * 5001 * automataLen, FROM_GPU));

    //**************************************************************************
    // EXTRACT CANDIDATES

    Candidates* candidates = new Candidates();
    candidates->reserve(automataLen);
    for (int i = 0; i < automataLen; ++i) {
        Candidate queryCandidates;

        int size = candidatesH[i * 5001];
        // printf("Size: %d\n", size);
        for (int j = 0; j < size; ++j) {
            // printf("Kandidat: %d\n", candidatesH[i * 5001 + j + 1]);
            queryCandidates.push_back(candidatesH[i * 5001 + j + 1]);
        }

        candidates->push_back(queryCandidates);
    }

    free(candidatesH);

    //**************************************************************************
    // CLEAN UP
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
    //**************************************************************************

    return static_cast<void*>(candidates);
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
    int automataLen, ChainGpu** database, int databaseLen,
    int* candidates) {

    int index = threadIdx.x;

    int candidatesSize = 0;

    if (index < automataLen) {
        // printf("Inside\n");
        candidates[index * 5001] = 0;

        int* table = automata[index]->table;

        for (int i = 0; i < databaseLen; ++i) {
            // printf("there is a db\n");
            int targetLen = database[i]->len;
            char* codes = database[i]->codes;

            int state = 0;
            int numHits = 0;
            // do transitions on the automaton and fill candidates

            // printf("targetlen %d\n", targetLen);
            for (int k = 0; k < targetLen; ++k) {
                int code = codes[k];

                while(table[state * TABLE_WIDTH + code] == 0 && state != 0) {
                    // fprintf(stderr, "no transition\n");
                    state = table[state * TABLE_WIDTH + FAIL_COL];
                }

                if (state == table[state * TABLE_WIDTH + code]) {
                    // printf("Character %c continuing\n", c + 'A');
                    continue;
                }

                state = table[state * TABLE_WIDTH + code];
                // fprintf(stderr, "Current state: %d\n", state);
                if (table[state * TABLE_WIDTH + FINAL_COL]) {
                    // fprintf(stderr, "A hit. Sign: %c\n", c + 'A');
                    // printf("A hit\n");
                    numHits++;
                }

            }

            if (numHits > 0) {
                candidates[index * 5001 + (candidatesSize + 1) % 5000] = i;
                candidatesSize++;
                candidates[index * 5001] = candidatesSize;
            }
        }

    }

    // printf("heklo worls\n");
    return;
}

// ***************************************************************************
