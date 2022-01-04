from kmers_dataset import *
from create_kmers_output import *

def main():
    print('Give inputs\n')
    k_mers, script_name, genome_file, genome_folder ,input_folder, output_folder = input().split()
    k_mers = int(k_mers)
    
    first_stage = create_kmers_output(k_mers, script_name, genome_file, genome_folder ,input_folder, output_folder)
    first_stage.fasta_generator()
    first_stage.folder_loop()
    

    second_stage = kmers_dataset(k_mers, output_folder, genome_folder)
    headers = second_stage.create_csv_headers()
    second_stage.create_dataset(headers)


if __name__ == "__main__":
    main()