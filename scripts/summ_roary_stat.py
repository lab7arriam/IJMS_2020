#/usr/bin/python3.7

import sys
import os
import csv
from collections import defaultdict

constant_cols = 14

def main():
    filename=sys.argv[1]
    row_dict=defaultdict(list)
    with open(filename,'r',newline='') as stat:
        cr = csv.reader(stat,delimiter='\t')
        
        for row in cr:
            if row[0]=='Count':
                for name in row[1:]:
                    #create a dictionary with lists where 0-ind stands for the number of genes in less then n strains while the second number has the opposite meaning
                    row_dict[name]=[0,0]
            else:
                if int(row[0])<=51 #same for 35:
                    i=1
                    for name in row_dict.keys():
                        row_dict[name][0]+=int(row[i])
                        i+=1
                else:
                    i=1
                    for name in row_dict.keys():
                        row_dict[name][1]+=int(row[i])
                        i+=1
    print(len([name for name in row_dict.keys() if row_dict[name][0]>=row_dict[name][1]]))
    #print(row_dict)
    new_names = list([name for name in row_dict.keys() if row_dict[name][0]<row_dict[name][1]])
    #print(new_names)

    with open ('filtered_strain_names.txt', 'w',newline='') as csvfile2:
        my_writer = csv.writer(csvfile2, delimiter='\t')
        for key in new_names:
            my_writer.writerow([key])

if __name__ == '__main__' :
    main()
