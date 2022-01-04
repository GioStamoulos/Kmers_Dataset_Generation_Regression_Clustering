import os
class create_kmers_output:

    def __init__(self, k_mers, script_name, genome_file, genome_folder ,input_folder, output_folder):
        self.k_mers = k_mers
        self.script_name = script_name
        self.genome_file = genome_file
        self.genome_folder = genome_folder
        self.input_folder = input_folder
        self.output_folder = output_folder

    # folder_loop --> insert fasta files to count_kmers perl script
    def folder_loop(self):
        try:
            os.mkdir(self.output_folder)
        except:
            pass
        for entry in os.scandir(self.input_folder):
            filename = str(os.path.basename(entry))
            os.system('perl {0} -s {1} < {2} > {3}'.format(self.script_name,self.k_mers, '{0}/{1}'.format(self.input_folder, filename), '{0}/{1}.csv'.format(self.output_folder, filename[:-6])))
             
    #fasta_generator function --> split genome's blast result output to fasta files
    def fasta_generator(self):
        f1 = open(self.genome_file, 'r').read()
        f2 = f1.split('>')
        f2.pop(0)
        n_iter = len(f2)
        print(n_iter, 'fasta files')
        os.mkdir('{0}/{1}'.format(self.genome_folder, self.k_mers))
        os.mkdir(self.input_folder)
        for i in range(n_iter):
            try:
                name_index = f2[i].index(' ')
            except:
                name_index = f2[i].index('\n')
            with open('{0}\{1}.fasta'.format(self.input_folder, f2[i][:name_index]), "w") as nf:
                nf.write(">")
                nf.write(f2[i])
                nf.close