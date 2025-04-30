# type: ignore #code of the poject
#Lynn Limbach & Max Henrotin
import sys

#retreive the filepath from the desired .fastq file
file_path = ""      # fastq/J29_B_CE_IonXpress_005.fastq  is one example of argument
if len(sys.argv) == 2:
    file_path = sys.argv[1]
else:
    print("Please provide the path of the FASTQ file.")
    sys.exit(1)

# this method extract in an array all the ACTG sequences from the .fastq file
def extract_sequences(path):
    sequences = []
    try:
        with open(path, "r") as f:
            lines = f.readlines()
            for i in range(1, len(lines), 4):  # extract line 2 from each 4 line bloc (fastq format)
                seq = lines[i].strip()
                if all(c in "ACTG" for c in seq):  # Verify the content of the file is in the right format
                    sequences.append(seq)
                else :
                    print("the provided .fastq file is not in the right format")
                    sys.exit(1)
        return sequences
    except FileNotFoundError:
        print(f"Error: The file {path} was not found.")
        sys.exit(1)
    except OSError as e:
        print(f"Error: Unable to open the file {path}. {e}")
        sys.exit(1)

#method that filter only the sequences within a given lenght region
def keep_seq_lenght(sequences, minlenght, maxlenght):
    if minlenght <= maxlenght :
        filtered_sequences = []
        for seq in sequences:
            sequence_lenght = len(seq)
            if sequence_lenght >= minlenght and sequence_lenght <= maxlenght:
                filtered_sequences.append(seq)
        return filtered_sequences
    else :
        print("minlenght > maxlenght")
        sys.exit(1)

#print a sequence list
def print_sequence_list(sequences):
    for seq in sequences:
        print("- " + seq + "\n")


#Execution
sequences = extract_sequences(file_path)
sequences296 = keep_seq_lenght(sequences, 296, 296)
print_sequence_list(sequences296)
print(len(sequences))
print(len(sequences296))