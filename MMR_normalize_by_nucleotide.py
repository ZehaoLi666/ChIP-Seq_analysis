# Goes through txt file with 3 columns: 
# {chromosome} {nucleotide} {read_count}
# and for each nucleotide position, divides the read count by the numbers
# of million mapped reads (MMR).
# Applicable for normalizing ChIP-seq or RNA-seq samples for generating
# a normalized genome browser track.
# September 18, 2020.

from sys import argv

if len(argv) != 4:
	print("Run with 3 arguments: input_norm.txt output_norm.txt {number_of_MMR}")
	exit(0)
else:
	pass

script, input_filename, output_filename, MMR_number = argv

MMR_number = float(MMR_number)
write_page = open(output_filename, 'w')

# Open input file and go through line-by-line, dividing read counts by
# the provided number of million mapped reads.
with open(input_filename) as input_file:
	for line in input_file:
		line = line.rstrip()
		current_line = line.split('\t') 
		
		chromosome = current_line[0]
		nucleotide= current_line[1]
			
		read_count = float(current_line[2].rstrip())	
		read_count = round(read_count / MMR_number, 3)
		
		write_page.write(f"{chromosome}\t{nucleotide}\t{read_count}\n")
	
			
write_page.close()

