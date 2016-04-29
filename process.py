import os
from numpy import *

def get_headers(file_name):
	heads_file = open(file_name)
	heads = []
	with heads_file as f:
		for line in f:
			string = line.split()
			for x in string:
				heads.append(x)
	return heads

def check_col(headers, x):
	for h in headers:
		if x == h:
			return True
	return False

def clean(headers, rootdir):
	matrix = []
	for subdir, dirs, files in os.walk(rootdir):
	    for file in files:
	        filepath = subdir + os.sep + file
	        if filepath.endswith(".txt"):
	            file = open(filepath)
	            data = []
	            cols = []
	            head = True
	            with file as f:
					for line in f:
						string = line.split()
						j = 0
						for x in string:
							if head:
								if check_col(headers,x):
									cols.append(True)
								else:
									cols.append(False)
							elif cols[j]:

								data.append(float(x))
							j = j + 1
						head = False
					if len(data) > 0:
						matrix.append(data)
	return matrix

def write_data(matrix, file_name):
	with open(file_name, 'w') as file:
		file.writelines('\t'.join(str(i)) + '\n' for i in matrix)

tumor_rootdir = "data_new/protein/tumor/Expression-Protein/MDA__MDA_RPPA_Core/Level_2"
tumor_output = "data_new/protein/all_tumor_protein.txt"
normal_rootdir = "data_new/protein/normal/Expression-Protein/MDA__MDA_RPPA_Core/Level_2"
normal_output = "data_new/protein/all_normal_protein.txt"
header_file = "data_new/protein/headers.txt"

heads = get_headers(header_file)
print "%d protein samples" % len(heads)

tumor = clean(heads, tumor_rootdir)
print "%d tumor samples" % len(tumor)
normal = clean(heads, normal_rootdir)
print "%d normal samples" % len(normal)

# write_data(tumor, tumor_output)
# write_data(normal, normal_output)
savetxt(tumor_output, tumor, delimiter='\t')
savetxt(normal_output, normal, delimiter='\t')