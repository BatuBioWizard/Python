import re
import sys
import os

if len(sys.argv) != 3:
    print("You should put an input like this: python myScript.py mySamFile.sam myInputTable.txt")
    sys.exit(1)  

gene_sam_file = sys.argv[1]
gene_txt_file = sys.argv[2]

print('You are using SAM file: {}'.format(gene_sam_file))
print('You are using TXT file: {}'.format(gene_txt_file))

def junction_finder(pos, cigar):
    new_pos = pos
    junction_start_point = []
    junction_end_point = []

    matches = re.finditer(r'(\d+)([MDN])', cigar)

    for match in matches:
        num = int(match.group(1)) 
        alph = match.group(2)

        if alph in ('M', 'D'):
            new_pos += num

        elif alph == 'N':
            junction_start_point.append(new_pos)
            new_pos += num
            junction_end_point.append(new_pos - 1)  
    return junction_start_point, junction_end_point

gene_location_dict = {}
try:
        with open(gene_txt_file) as fileA:
            for seqA in fileA:
                data = seqA.strip().split('\t')

                if len(data) == 3:

                    gene_id = data[0]
                    genomic_location = data[2]
                    genomic_location = genomic_location.split(':')

                    if len(genomic_location) == 2:
                        chromosome, gene_location = genomic_location[0], genomic_location[1]
                        coordinates, strand = gene_location.split('(')
                        start_point, end_point = coordinates.split("..")
                        strand = strand.rstrip(')')
                        start_point = int(start_point.replace(",", ""))
                        end_point = int(end_point.replace(",", ""))
                        gene_location_dict[gene_id] = {
                        'chromosome': chromosome,
                        'start_point': start_point,
                        'end_point': end_point }
except FileNotFoundError:
        infoLogger.error(f'File {gene_txt_file} not found. Please try again')
except ValueError:
        infoLogger.error(f'Invalid value type, please check {start_point} and {end_point} assert in the incorrect format Int/Str')
        raise SystemExit(1)

sam_data_dict = {}
reference_count_dict = {}
output_file = "output.txt"

try:
    with open(output_file, 'w') as output_file, open(gene_sam_file) as fileB:
        output_file.write("Gene id\tJunction Start\tJunction End\tNumber of reads supporting the junction\n")
        prev_gene_id = None
        for seqB in fileB:
            if seqB.startswith('HWI'):
                data = seqB.strip().split('\t')
                if len(data) >= 20:
                    rname = data[2]
                    pos = int(data[3])
                    cigar = data[5]
                    NH_position = data[-1]
                    NH_values = NH_position.split(':')
                    NH = int(NH_values[2])
                    reference_pos_start, reference_pos_end = junction_finder(pos, cigar)
                    # Find gene id based on reference positions
                    gene_id = None
                    for key, value in gene_location_dict.items():
                        if value['chromosome'] == rname:
                            for start, end in zip(reference_pos_start, reference_pos_end):
                                if value['start_point'] <= start <= value['end_point'] and value['start_point'] <= end <= value['end_point']:
                                    gene_id = key
                                    reference_count_dict[(gene_id, start, end)] = reference_count_dict.get((gene_id, start, end), 0) + 1
                                    break
                            if gene_id is not None:
                                break
        for (gene_id, start, end), count in reference_count_dict.items():
            output_file.write(f"{gene_id}\t{start}\t{end}\t{count}\n")
            if gene_id != prev_gene_id:
                output_file.write('\n')
            prev_gene_id = gene_id

except FileNotFoundError:
    infoLogger.error(f'File {gene_sam_file} not found. Please try again')  
except ValueError:
    infoLogger.error(f'Invalid value type, {NH} or {reference_pos_start} or {reference_pos_end} assert in the incorrect format Int/Str ')
    raise SystemExit(2)

