#include<stdio.h>
#include<stdlib.h>
#include "submodels.h"

int main(int argc, char *argv[])
{

    if (argc !=3)
    {
        printf("Usage: ./a.out <input file> <output directory>\n");
        abort();
    }
    
    // input file and output directory
    char *inputFile = argv[1];
    char *outputDir = argv[2];
    
    // read input parameters from the input file
    input params;
    read_input_file(&params, inputFile);
    
    // call to solver
    solve(params, outputDir);

    return 0;
}

