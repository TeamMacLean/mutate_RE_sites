Certainly! Here's a README that explains how to use the script:

## Restriction Enzyme Site Mutation Script

The Restriction Enzyme Site Mutation Script is a Python program that recognizes restriction enzyme sites in nucleotide sequences and performs mutations or replacements of nucleotides in those sites. The script ensures that the mutations do not change the amino acid and selects the new codon with the highest score of occurrence.

### Prerequisites

- Python 3.x
- Biopython package (`pip install biopython`)

### Usage

1. Clone the repository or download the script to your local machine.
2. Open a terminal or command prompt and navigate to the directory where the script is located.

```
$ cd /path/to/restriction-enzyme-script/
```

3. Run the script using the following command:

```
$ python script.py --input <input_file> --re <enzyme_file> --codon <codon_file> --output <output_file> --report <report_file>
```

Replace the `<input_file>`, `<enzyme_file>`, `<codon_file>`, `<output_file>`, and `<report_file>` placeholders with the appropriate file paths or names.

### Arguments

- `--input`: Specifies the path to the input FASTA file containing nucleotide sequences.
- `--re`: Specifies the path to the file containing the list of restriction enzymes and their corresponding restriction sites.
- `--codon`: Specifies the path to the file containing codons and their occurrence scores.
- `--output`: Specifies the name of the output file that will contain the mutated sequences. (Default: output.fasta)
- `--report`: Specifies the name of the report file that will contain information about the changes made to the sequences. (Default: report.txt)

### Output

- The script will generate an output file (`output.fasta`) containing the mutated sequences.
- It will also generate a report file (`report.txt`) that provides information about the changes made to the sequences, including the original restriction site, the mutated site, and the start and end positions of the changes.

### Examples

Here are a few examples of how to run the script:

```
$ python script.py --input sequences.fasta --re enzymes.txt --codon codons.txt --output mutated_sequences.fasta --report changes_report.txt
```

```
$ python script.py --input input.fasta --re restriction_sites.txt --codon codon_scores.txt
```

### Notes

- Make sure the input FASTA file, enzyme file, and codon file are in the correct format as specified in the script.
- The script uses the Biopython library, so ensure it is installed before running the script.
- The mutated sequences will be saved in the output file specified.
- The changes made to the sequences will be recorded in the report file specified.

That's it! You can now use the Restriction Enzyme Site Mutation Script to recognize restriction enzyme sites and perform nucleotide mutations or replacements.
