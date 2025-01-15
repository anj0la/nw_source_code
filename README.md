# Needleman-Wunsch Algorithm - Serial and Parallel Implementation

This repository contains the serial and parallel source code of the Needleman-Wunsch algorithm. 
When running the program, we assume that the two sequences are of equal length.

## Required Compilers
- GCC compiler
- Nvidia NVC compiler (need a Linux system or virtual OS with a Nvdia GPU)

## To run the serial version of the algorithm:

1. Compile the program using the following line:
```bash
 gcc nw_serial.c -o nw_serial
```
2. Execute the code in your chosen terminal (system terminal or using an IDE for C):
```bash
./[executable_file_name] [first file] [second file] [length] [match] [mismatch] [gap_penalty]
```   
3. With nw_serial, seq1.txt and seq2.txt:
```bash
./nw_serial seq1.txt seq2.txt 4 1 -1 -2
```
## To run the parallel version of the algorithm:

1. Compile the program using the following line:
```bash
nvc -fast -Minfo=all -mp nw_openmp.c -o nw_openmp
```
2. Execute the code in your chosen terminal (system terminal or using an IDE for C):
```bash
./[executable_file_name] [first file] [second file] [length] [thread] [match] [mismatch] [gap_penalty]
```
3. With nw_serial, seq1.txt, seq2.txt and 8 threads:
```bash
./nw_openmp seq1.txt seq2.txt 4 8 1 -1 -2
```
