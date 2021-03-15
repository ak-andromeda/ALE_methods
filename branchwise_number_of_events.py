import os, sys, re
from collections import defaultdict


def main():

    branchwise = defaultdict(list)
    input = [file for file in os.listdir(".") if file.endswith("uml_rec")]

    for file in input:
        inhandle = open(file)
        for line in inhandle:
            line_elements = re.split("\t", line.rstrip())
            if (line_elements[0] == "S_terminal_branch") or (line_elements[0] == "S_internal_branch"):
                if len(branchwise[line_elements[1]]) == 0:
                    branchwise[line_elements[1]].append(float(line_elements[2]))
                    branchwise[line_elements[1]].append(float(line_elements[3]))
                    branchwise[line_elements[1]].append(float(line_elements[4]))
                    branchwise[line_elements[1]].append(float(line_elements[5]))
                    branchwise[line_elements[1]].append(float(line_elements[6]))
                else:
                    branchwise[line_elements[1]][0] += float(line_elements[2])
                    branchwise[line_elements[1]][1] += float(line_elements[3])
                    branchwise[line_elements[1]][2] += float(line_elements[4])
                    branchwise[line_elements[1]][3] += float(line_elements[5])
                    branchwise[line_elements[1]][4] += float(line_elements[6])

    print("Branch\tDuplications\tTransfers\tLosses\tOriginations\tCopyNum")
    for entry in branchwise:
        print(entry + "\t" + str(branchwise[entry][0]) + "\t" + str(branchwise[entry][1]) + "\t" + str(branchwise[entry][2]) + "\t" + str(branchwise[entry][3]) + "\t" + str(branchwise[entry][4]))


if __name__ == "__main__":
    main()
