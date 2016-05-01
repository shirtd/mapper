import os
from numpy import *
import re

def splt(delimiters, string, maxsplit=0):
    regexPattern = '|'.join(map(re.escape, delimiters))
    return re.split(regexPattern, string, maxsplit)

def get_protein(string):
	delimiters = ["R-V_", "R-C_", "R-E_", "G-V_", "G-C_", "G-E_", "M-V_", "M-C_", "M-E_" ]
	stng = splt(delimiters, string)
	# if len(stng) < 2:
	# 	print "error" + string
	return stng[0]

def get_headers(file_name):
	delimiters = ["R-V_", "R-C_", "R-E_", "G-V_", "G-C_", "G-E_", "M-V_", "M-C_", "M-E_" ]
	heads_file = open(file_name)
	heads = []
	with heads_file as f:
		for line in f:
			string = line.split()
			for x in string:
				heads.append(get_protein(x))
	return heads

def check_col(headers, x):
	for h in headers:
		if x == h:
			return True
	return False

def clean(headers, rootdir):
	matrix = []
	samples = []
	for subdir, dirs, files in os.walk(rootdir):
	    for file in files:
	        filepath = subdir + os.sep + file
	        if filepath.endswith(".txt"):
	            fle = open(filepath)
	            data = []
	            cols = []
	            head = True
	            with fle as f:
					for line in f:
						string = line.split()
						j = 0
						for x in string:
							if head:
								if check_col(headers,get_protein(x)):
									cols.append(True)
									samples.append(file)
								else:
									cols.append(False)
							elif cols[j]:
								data.append(float(x))
							j = j + 1
						head = False
					if len(data) > 0:
						matrix.append(data)
	return [matrix, samples]

# def get_barcodes(samples, sample_map):
# 	file = open(sample_map)
# 	barcodes = []
# 	with file as f:
# 		for line in f:
# 			string = line.split()
# 			if len(string) == 2:
# 				for s in samples:
# 					if string[0] == s:
# 						# print s + " -> " + string[1]
# 						barcodes.append(string[1])

# 	return barcodes

# def write_matrix(matrix, file_name):
# 	with open(file_name, 'w') as file:
# 		file.writelines('\t'.join(str(i)) + '\n' for i in matrix)

# def write_list(lst, file_name):
# 	file = open(file_name, 'w')
# 	for item in lst:
# 		file.write("%s\n" % item)



tumor_rootdir = "data_new/protein/tumor/Expression-Protein/MDA__MDA_RPPA_Core/Level_2"
# tumor_sample_map = "data_new/protein/tumor/FILE_SAMPLE_MAP.txt"
# tumor_barcode_output = "data_new/protein/tumor/tumor_barcodes.txt"
tumor_output = "data_new/protein/all_tumor_protein.txt"
normal_rootdir = "data_new/protein/normal/Expression-Protein/MDA__MDA_RPPA_Core/Level_2"
# normal_sample_map = "data_new/protein/normal/FILE_SAMPLE_MAP.txt"
normal_output = "data_new/protein/all_normal_protein.txt"
header_file = "data_new/protein/headers.txt"

heads = get_headers(header_file)
print "%d proteins" % len(heads)

tumor = clean(heads, tumor_rootdir)
tumor_mat = tumor[0]
# tumor_samples = tumor[1]
# tumor_barcodes = get_barcodes(tumor_samples, tumor_sample_map)
# print "%d tumor sample barcodes" % len(tumor_barcodes)
# write_list(tumor_barcodes, tumor_barcode_output)

print "%d tumor samples" % len(tumor_mat)
normal = clean(heads, normal_rootdir)
normal_mat = normal[0]
# normal_samples = normal[1]
print "%d normal samples" % len(normal_mat)

# print normal_mat
# print tumor_mat

# write_data(tumor, tumor_output)
# write_data(normal, normal_output)
savetxt(tumor_output, tumor_mat, delimiter='\t')
savetxt(normal_output, normal_mat, delimiter='\t')