import argparse
import pysam

##
##
##

#inversion 1: 13-14Mb & 319~320Mb
#inversion 2: 26-27Mb & 56~57Mb
#inversion 3: 243~244Mb & 254~255Mb


transversion_locations = { 13e6:319e6, 26e6:56e6,  243e6:254e6}

def process_reads(bam_input):
    in_bam = pysam.AlignmentFile(bam_input,"rb")
    current_read = ""
    start_positions=[]
    reads=[]
    for read in in_bam:
        if read.mapping_quality >= 60:
            if read.reference_name == "smed_chr1":
                if current_read == read.query_name:
                    start_position = read.reference_start
                    start_positions.append(start_position)
                    reads.append(read)
                if current_read != read.query_name:
                    matched_one= False 
                    matched_two= False 
                    for key, items in transversion_locations.items():
                        for starts in start_positions:
                            if starts > key and starts < (key + 1e6):
                                matched_one = True
                            if starts > items and starts < (items+ 1e6):
                                matched_two= True
                        if matched_one and matched_two:
                            for read in reads:
                                print(key, items, read)
                        matched_one = False
                        matched_two = False
                    start_positions=[]
                    reads=[]
                    start_position = read.reference_start
                    start_positions.append(start_position)
                    reads.append(read)
                current_read = read.query_name



def main():
    parser = argparse.ArgumentParser(description="Identify planarian pacbio reads that span the junctions")
    parser.add_argument("bam_in")
    args = parser.parse_args()
    bam_input = args.bam_in
    process_reads(bam_input)
if __name__=="__main__":
    main()
