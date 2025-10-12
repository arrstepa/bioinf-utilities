import os

def convert_multiline_fasta_to_oneline(input_fasta):
    """
    A tool for converting multiline fasta files to oneline format

    Arguments:
    input_fastq: str

    Returns str
    """
    oneline_fasta = []
    current_seq = ''
    with open(input_fasta) as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if current_seq:
                    oneline_fasta.append(current_seq)
                    current_seq = ''
                oneline_fasta.append(line)
            else:
                current_seq += line
        if current_seq:
            oneline_fasta.append(current_seq)
    output_fasta = os.path.join(os.path.dirname(input_fasta), 'output_fasta.txt')
    with open(output_fasta, 'a') as file:
        for line in oneline_fasta:
            file.write(line + "\n")
