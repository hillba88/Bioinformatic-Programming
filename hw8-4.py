#! /usr/bin/env python3

import Bio.Seq
import time

# function reads a file, cretes a sequence instance
# for each line, and returns list of sequence instances
def seq_list(file):
    ret_list = []
    try:
        with open(file, "r") as seq_file:
            for line in seq_file.readline():
                ret_list.append(Bio.Seq.Seq(line))
    except IOError:
        print("IO Error!")

    return ret_list

# function runs seq_list a specified number of times and measures time elapsed
def run_func(times, file):
    start = time.time()
    for num in range(times):
        seq_list(file)

    end = time.time()

    return "Ran {0}, took {1}".format(times, end - start)

infile = "hw8/hw-8-4.txt"

# run function with different iteration amounts
print(run_func(100, infile))
print(run_func(1000, infile))
print(run_func(10000, infile))
print(run_func(100000, infile))
print(run_func(1000000, infile))
