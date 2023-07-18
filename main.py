'''
Aya khaled Mohamed Abdel Raheem  (20198016)
Amira Hamdy Sayed                (20198013)
Fatma Mohamed Abdelgafour okasha (20198064)
Yousef Ahmed Mohamed Nasar       (20198103)
Mahmoud Khaled Helmy Abo-Elmagd  (20188045)
'''

##note:
##Used Python version should be at least 3.10 because of using Match/Case keyword


from Bio import SeqIO
# from Bio import pairwise2
from Bio.Align import PairwiseAligner
from Bio.Blast import NCBIWWW, NCBIXML
import sys
import getopt


def count_gc_nucleotides(seq):
    g_and_c_occurrences = 0
    for nucleotide in seq:
        if nucleotide in 'GCgc':
            g_and_c_occurrences += 1
    return g_and_c_occurrences


def gc(seq):
    corrected_seq = filter_nbases(seq)
    g_and_c_occurrences = count_gc_nucleotides(corrected_seq)

    return (g_and_c_occurrences / len(corrected_seq)) * 100


def transcribe(seq):
    return seq.translate(str.maketrans('tT', 'uU'))

# =========================================================================


def sequence_aligment(seq1, seq2, o=""):
    # pairwise2 is deprecated so we use pairwisealigner instead
    # alignments = pairwise2.align.globalxx(seq1, seq2)
    alignments = PairwiseAligner().align(seq1, seq2)
    if o == "":
        for alignment in alignments:
            print(alignment)
    else:
        outputfile = open(o, "w")
        for alignment in alignments:
            alignment = str(alignment) + '\n'
            outputfile.write(alignment)
        outputfile.close()
        print('Done.')


def seq_alignment_files(file1, file2, o=""):
    inputfile1 = list(SeqIO.parse(file1, "fasta"))
    inputfile2 = list(SeqIO.parse(file2, "fasta"))
    seq1 = str(inputfile1[0].seq)
    seq2 = str(inputfile2[0].seq)
    sequence_aligment(seq1, seq2, o)

# =========================================================================


def calc_nbases(seq):
    # This command takes a seq and calculates its nbases
    my_seq = list(seq)
    count = 0
    for i in my_seq:
        if i in ['n', 'N']:
            count += 1
    if count != 0:
        print("DNA sequence has %d undefined bases " % count)
    else:
        print("DNA sequence has no undefined")


def convert_to_fasta(file_path):
    # file path for genbank file
    # This command converts the input genbank file with multiple records onto a fasta formatted file.
    #  The output is to be written in a different output fasta file.
    with open(file_path) as input_handle, open("converted_file.fasta", "w") as output_handle:
        sequences = SeqIO.parse(input_handle, "genbank")
        count = SeqIO.write(sequences, output_handle, "fasta")
        print("Converted %i records" % count)


def reverse_complement(seq):
    # This command takes a seq as a string and returns its reverse complement
    my_seq = list(seq)
    my_seq = my_seq[::-1]
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N', 'a': 't', 'c': 'g',
                      'g': 'c', 't': 'a', 'n': 'n'}
    my_seq = [basecomplement[base] for base in my_seq]
    return ''.join(my_seq)

# =========================================================================


def seq_alignment_online(sequence, file_name=""):
    # sequence to try it with: ggtaagtcctctagtacaaacacccccaatattgtgatataattaaaattatattcatattctgttgccagaaaaaacacttttaggctatattagagccatcttctttgaagcgttgtc
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)
    blast_record = NCBIXML.read(result_handle)

    if file_name == "":
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('e value:', hsp.expect)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)
    else:
        f = open(file_name, "w")
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                f.write(str(alignment.title)+'\n')
                f.write(str(alignment.length)+'\n')
                f.write(str(hsp.expect)+'\n')
                f.write(str(hsp.query)+'\n')
                f.write(str(hsp.match)+'\n')
                f.write(str(hsp.sbjct)+'\n')
        f.close()


def is_valid(sequence,type):
    "This function takes sequence & it's type and returns boolean value explains if seq is valid or not "
    protein_seq="ARNDCQEGHILKMFPSTWYVarndcqeghilkmfpstwyv"
    dna_seq="acgtACGT"
    rna_seq="acguACGU"
    for i in range(len(sequence)):
        if type=="protein" and sequence[i] in protein_seq :
            check=True
        elif type=="dna" and sequence[i] in dna_seq:
            check=True
        elif  type=="rna" and sequence[i] in rna_seq:
            check=True
        else :
            check = False
            break

    return check




# =========================================================================


def merge_fasta(*files, file=''):
    """This command takes any number of fasta files (at least two) and merge their contents into one fasta output
       file"""
    for i in files:
        if ".fasta" not in i:
            print("YOU MUST ENTER FASTA FILES!!")
            return False
    if file == '':
        for i in files:
            input_lines = list(SeqIO.parse(i, "fasta"))
            for j in input_lines:
                print(str(j))

    else:
        with open(file, 'a') as file:
            for i in files:
                input_lines = list(SeqIO.parse(i, "fasta"))
                file.writelines(input_lines)
            file.close()


def filter_nbases(Seq):
    """This command takes a seq and returns the Seq after removing n bases"""
    Seq = Seq.upper()
    for i in Seq:
        if i not in 'ATCGN':
            return 'Invalid Seq'
        else:
            Seq = Seq.replace("N", "")
    return Seq

# =========================================================================


# Main
argv_len = len(sys.argv)
if argv_len < 2:
    sys.exit('You should enter a command.')

command = sys.argv[1]

one_arg_commands = ['gc', 'transcribe', 'reverse_complement',
                    'calc_nbases', 'filter_nbases', 'online_alignment', 'convert_to_fasta']
two_args_commands = ['is_valid', 'seq_alignment', 'seq_alignment_files']
at_least_two_args_commands = ['merge_fasta']

if argv_len < 3 and command in one_arg_commands:
    sys.exit('Command require 1 arg')
if argv_len < 4 and command in two_args_commands:
    sys.exit('Command require 2 args')
if argv_len < 4 and command in at_least_two_args_commands:
    sys.exit('Command require at least 2 args')

args = sys.argv[2:]
opts, args = getopt.gnu_getopt(sys.argv[1:], 'o:')
command_args = args[1:]
output_option = opts[0][1] if len(opts) > 0 else ''

match command:
    case "gc":
        seq = command_args[0]
        print('GC percentage:', gc(seq))

    case "transcribe":
        seq = command_args[0]
        print('Transcription:', transcribe(seq))

    case "reverse_complement":
        seq = command_args[0]
        print('Reverse complement:', reverse_complement(seq))

    case "calc_nbases":
        seq = command_args[0]
        calc_nbases(seq)

    case "is_valid":
        seq = command_args[0]
        seq_type = command_args[1]

        print('Valid' if is_valid(seq, seq_type) else 'Invalid')

    case "filter_nbases":
        seq = command_args[0]
        print('Nbases filtered seq:', filter_nbases(seq))

    case "seq_alignment":
        seq1 = command_args[0]
        seq2 = command_args[1]

        sequence_aligment(seq1, seq2, output_option)

    case "seq_alignment_files":
        try:
            file1 = open(command_args[0])
            file2 = open(command_args[1])

            seq_alignment_files(file1, file2, output_option)
        except IOError:
            sys.exit('A file does not exists.')

    case "online_alignment":
        seq = command_args[0]

        seq_alignment_online(seq, output_option)

    case "merge_fasta":
        files = command_args

        merge_fasta(*files, output_option)

    case "convert_to_fasta":
        file_name = command_args[0]

        convert_to_fasta(file_name)

    case _:
        sys.exit("Incorrect command.")
