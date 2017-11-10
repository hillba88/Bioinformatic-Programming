#! /usr/bin/env python3

import subprocess

# run the Linux command sort to order the data in hw-8-1.txt
result = subprocess.run("sort hw8/hw-8-1.txt", shell=True,
                              stdout=subprocess.PIPE, stderr=subprocess.PIPE)
linux_result = result.stdout.decode('ascii')

# use python command to sort hw-8-1.txt using fileIO and store list
python_result = ""
try:
    with open("hw8/hw-8-1.txt") as my_file:
        lines = my_file.readlines()
        lines.sort()
        for item in lines:
            python_result += item

except IOError:
    print("IO Error!")

# test to see if both methods produce the same results. Need to split strings
# into lists of objects and compare each object.
for item in enumerate(python_result.split("\n")):
        print("{0}=={1}=={2}".format(item[1], linux_result.split("\n")[item[0]],
                                     item[1]==linux_result.split("\n")[item[0]]))
