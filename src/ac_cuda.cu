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

__global__ static void findCandidatesBlocks(TableGpu** automata, 
    int automatonIndex, ChainGpu** database, int databaseLen,
    int numBlocks, int* candidates);

// ***************************************************************************
// PUBLIC
extern void* indicesTableCreateGpu(Chain** database, 
    int databaseStart, int databaseLen, void* automata,
    int automataLen, int seedLen, Scorer* scorer) {

    fprintf(stderr,"Creating indices\n");
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
    fprintf(stderr,"Sending automata to gpu\n");
    for (int i = 0; i < automataLen; ++i) {
        TabNode* autH = (*aut)[i];

        TableGpu* hostCopy;

        gpuTables.push_back(copyTableToGpu(autH, &hostCopy));
        hostTables.push_back(hostCopy);
        // fprintf(stderr,"%p\n", gpuTables[i]);
    }
    fprintf(stderr,"Done\n");

    TableGpu** gpuTablesD;
    CUDA_SAFE_CALL(cudaMalloc(&gpuTablesD, sizeof(TableGpu*) * automataLen));
    CUDA_SAFE_CALL(cudaMemcpy(gpuTablesD, &gpuTables[0], 
        sizeof(TableGpu*) * automataLen, TO_GPU));    

    //**************************************************************************

    //**************************************************************************
    // SEND DATABASE TO GPU

    fprintf(stderr,"Sending database to gpu.\n");
    for (int i = databaseStart; i < databaseLen; ++i) {
        ChainGpu* chainD;
        ChainGpu* chainH;

        chainGpuCreate(database[i], &chainD, &chainH);

        gpuChains.push_back(chainD);
        hostChains.push_back(chainH);
    }
    fprintf(stderr,"DOne\n");

    ChainGpu** chainsGpuD;
    CUDA_SAFE_CALL(cudaMalloc(&chainsGpuD, sizeof(ChainGpu*) * numTargets));
    CUDA_SAFE_CALL(cudaMemcpy(chainsGpuD, &gpuChains[0], 
        sizeof(ChainGpu*) * numTargets, TO_GPU));

    //**************************************************************************
    fprintf(stderr,"Allocating candidates.\n");
    int* candidatesD;
    CUDA_SAFE_CALL(cudaMalloc(&candidatesD, sizeof(int) * 5001 * automataLen));
    int* candidatesH = (int*) malloc(sizeof(int) * 5001 * automataLen);

    //**************************************************************************
    // INVOKE KERNEL
    int grid_x = min(5000, databaseLen);
    dim3 dimGrid(grid_x, 1, 1);
    dim3 dimBlock(1, 1, 1);

    for (int i = 0; i < automataLen; ++i) {
        fprintf(stderr,"Invoking kernel\n");
        
        findCandidatesBlocks<<<dimGrid, dimBlock>>>(gpuTablesD, i, chainsGpuD,
            databaseLen, grid_x, candidatesD);
    }

    // int grid_x = automataLen;
    // int block_x = 1;

    // dim3 dimGrid(grid_x,1,1);
    // // fprintf(stderr,"Automata len: %d\n", automataLen);
    // dim3 dimBlock(block_x,1,1);
    // findCandidates<<<dimGrid, dimBlock>>>(gpuTablesD, automataLen, chainsGpuD, 
    //     numTargets, candidatesD);

    CUDA_SAFE_CALL(cudaMemcpy(candidatesH, candidatesD, 
        sizeof(int) * 5001 * automataLen, FROM_GPU));

    //**************************************************************************
    // EXTRACT CANDIDATES

    fprintf(stderr,"Extracting candidates\n");
    Candidates* candidates = new Candidates();
    candidates->reserve(automataLen);

    for (int i = 0; i < automataLen; ++i) {
        Candidate queryCandidates;

        for (int j = 0; j < grid_x; ++j) {
            if (candidatesH[i * 5001 + j + 1] > -1) {
                queryCandidates.push_back(candidatesH[i * 5001 +j + 1]);
            }
        }

        candidates->push_back(queryCandidates);
    }

    // for (int i = 0; i < automataLen; ++i) {
    //     Candidate queryCandidates;

    //     int size = candidatesH[i * 5001];
    //     // fprintf(stderr,"Size: %d\n", size);
    //     for (int j = 0; j < size; ++j) {
    //         // fprintf(stderr,"Kandidat: %d\n", candidatesH[i * 5001 + j + 1]);
    //         queryCandidates.push_back(candidatesH[i * 5001 + j + 1]);
    //         // fprintf(stderr, "Candidate is: %d\n", candidatesH[i * 5001 + j + 1]);
    //     }

    //     candidates->push_back(queryCandidates);
    // }

    free(candidatesH);
    fprintf(stderr,"Done\n");

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
    fprintf(stderr,"Done cleaning\n");
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
    CUDA_SAFE_CALL(cudaMalloc(&statesD, 
        sizeof(int) * (*hostCopy)->numStates * TABLE_WIDTH));

    CUDA_SAFE_CALL(cudaMemcpy(statesD, (*hostCopy)->table,
        sizeof(int) * (*hostCopy)->numStates * TABLE_WIDTH, 
        cudaMemcpyHostToDevice));

    (*hostCopy)->table = statesD;

    TableGpu* autD;
    CUDA_SAFE_CALL(cudaMalloc(&autD, sizeof(TableGpu)));

    CUDA_SAFE_CALL(cudaMemcpy(autD, (*hostCopy), sizeof(TableGpu), 
        cudaMemcpyHostToDevice));

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
    const char* codes = chainGetCodes(chain);

    int rest = len % 4;
    int mallocLen = rest > 0 ? 1 : 0;
    char4* codesH = (char4*) malloc(sizeof(char4) * (len / 4 + mallocLen));

    for (int i = 0; i < (len - rest) / 4; ++i) {
        codesH[i].x = codes[i * 4];
        codesH[i].y = codes[i * 4 + 1];
        codesH[i].z = codes[i * 4 + 2];
        codesH[i].w = codes[i * 4 + 3];
    }

    if (rest) {
        codesH[len / 4].x = codes[len-rest];

        if (rest >= 2) {
            codesH[len / 4].y = codes[len-rest + 1];
        } else {
            codesH[len / 4].y = 0;
        }

        if (rest == 3) {
            codesH[len / 4].z = codes[len-rest + 2];
        } else {
            codesH[len / 4].z = 0;
        }

        codesH[len / 4].w = 0;
    }

    char4* codesD;

    CUDA_SAFE_CALL(cudaMalloc(&codesD, sizeof(char4) * (len / 4 + mallocLen)));
    CUDA_SAFE_CALL(cudaMemcpy(codesD, codesH, 
        sizeof(char4) * (len / 4 + mallocLen), TO_GPU));

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

    int index = blockIdx.x;
    int candidatesSize = 0;

    if (index < automataLen) {
        candidates[index * 5001] = 0;
        int* table = automata[index]->table;

        for (int i = 0; i < databaseLen; ++i) {
            // printf("there is a db\n");
            int targetLen = database[i]->len;
            char4* codes = database[i]->codes;

            // if (index == 0 && threadIdx.x == 0) {
                // printf("Analyzing target %d\n", i);
                // printf("Target len %d\n", targetLen);
            // }

            // copy codes to shared memory
            // char4 ch4;
            // for (int k = threadIdx.x; k < targetLen / 4; k += 32) {
            //     ch4 = codes[k];

            //     // if (index == 0) {
            //         // printf("Thread %d Current chars: %d\n", threadIdx.x, k);
            //     // }

            //     int j = k * 4;

            //     // if (index == 0) {
            //     //     printf("Setting codesShared to %d\n", j);
            //     // }

            //     codesShared[j] = ch4.x;
            //     codesShared[j+1] = ch4.y;
            //     codesShared[j+2] = ch4.z;
            //     codesShared[j+3] = ch4.w;
            // }

            // // read the rest
            // if ((targetLen % 4) && threadIdx.x == 0) {
            //     int k = targetLen % 4;

            //     codesShared[targetLen - k] = codes[targetLen / 4].x;
            //     codesShared[targetLen - k + 1] = codes[targetLen / 4].y;
            //     codesShared[targetLen - k + 2] = codes[targetLen / 4].z;
            //     codesShared[targetLen - k + 3] = codes[targetLen / 4].w;
            // } else {

            // }

            // __syncthreads();

            // do transitions on the automaton and fill the candidates
            // if (threadIdx.x == 0) {
            int state = 0;
            int numHits = 0;
            char buff[4];
            // printf("targetlen %d\n", targetLen);
            int limit = targetLen % 4 ? targetLen / 4 + 1 : targetLen / 4;
            // printf("limit %d\n", limit);

            for (int k = 0; k < limit; ++k) {
                char4 in = codes[k];
                buff[0] = in.x;
                buff[1] = in.y;
                buff[2] = in.z;
                buff[3] = in.w;
                // printf("Current code %d\n", code);
                for (int j = 0; j < 4; ++j) {
                    char code = buff[j];

                    while(table[state * TABLE_WIDTH + code] == 0 && state != 0) {
                        state = table[state * TABLE_WIDTH + FAIL_COL];
                    }

                    if (state == table[state * TABLE_WIDTH + code]) {
                        continue;
                    }

                    state = table[state * TABLE_WIDTH + code];
                    if (table[state * TABLE_WIDTH + FINAL_COL]) {
                        numHits++;
                    }
                }

            }

            if (numHits > 0) {
                candidates[index * 5001 + (candidatesSize % 5000 + 1)] = i;
                candidatesSize = min(candidatesSize + 1, 5000);
                // candidatesSize = candidatesSize + 1 > 5000 ? 5000 : candidatesSize + 1;
            }
                // printf("Thread id: %d num hits: %d\n", index, numHits);
            // }
            // __syncthreads();
        }

        // printf("Candidates size: %d\n", candidatesSize);
        if (threadIdx.x == 0)
            candidates[index * 5001] = candidatesSize;
    }
}


__global__ static void findCandidatesBlocks(TableGpu** automata, 
    int automataIndex, ChainGpu** database, int databaseLen,
    int numBlocks, int* candidates) {

    int index = automataIndex;

    candidates[index * 5001] = 0;
    candidates[index * 5001 + 1 + blockIdx.x] = -1;
    int* table = automata[index]->table;

    for (int i = blockIdx.x; i < databaseLen; i += numBlocks) {
        // printf("there is a db\n");
        int targetLen = database[i]->len;
        char4* codes = database[i]->codes;

        // do transitions on the automaton and fill the candidates
        int state = 0;
        int numHits = 0;
        char buff[4];
        // printf("targetlen %d\n", targetLen);
        int limit = targetLen % 4 ? targetLen / 4 + 1 : targetLen / 4;
        // printf("limit %d\n", limit);

        for (int k = 0; k < limit; ++k) {
            char4 in = codes[k];
            buff[0] = in.x;
            buff[1] = in.y;
            buff[2] = in.z;
            buff[3] = in.w;
            // printf("Current code %d\n", code);
            for (int j = 0; j < 4; ++j) {
                char code = buff[j];

                while(table[state * TABLE_WIDTH + code] == 0 && state != 0) {
                    state = table[state * TABLE_WIDTH + FAIL_COL];
                }

                if (state == table[state * TABLE_WIDTH + code]) {
                    continue;
                }

                state = table[state * TABLE_WIDTH + code];
                if (table[state * TABLE_WIDTH + FINAL_COL]) {
                    numHits++;
                }
            }

        }

        if (numHits > 0) {
            candidates[index * 5001 + 1 + blockIdx.x] = i;
        }
    }
}

// ***************************************************************************
