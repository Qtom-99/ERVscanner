import argparse
import re

def load_dict_file(dictionary_file):
    """
    Load the 6-column dictionary file into a dictionary mapping Accession to Name and Classification.
    """
    dict_name = {}
    dict_classification = {}
    
    with open(dictionary_file, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith("Accession"):  # Skip header
                continue
            columns = line.strip().split('\t')
            if len(columns) >= 3:
                accession = columns[0]
                name = columns[1]
                classification = columns[2]
                dict_name[accession] = name
                dict_classification[accession] = classification
                
    return dict_name, dict_classification

def process_input_file(input_file, output_file, dict_name, dict_classification):
    """
    Process the input file to:
    - Remove decimal part of DF* identifiers in column 2 and 3.
    - Replace identifiers in column 2 and 3 with corresponding values from the dictionary.
    """
    pattern = re.compile(r"(.*?)\.\d+$")  # Regex to match DF* identifiers with decimal part
    
    with open(input_file, 'r', encoding='utf-8') as f_in, open(output_file, 'w', encoding='utf-8') as f_out:
        for line in f_in:
            columns = line.strip().split('\t')
            if len(columns) < 3:
                continue  # Skip lines with fewer than 3 columns
            
            # Remove decimal part in column 2 and 3
            for i in [1, 2]:
                match = pattern.match(columns[i])
                if match:
                    columns[i] = match.group(1)
            
            # Replace column 2 and 3 using dictionary
            columns[1] = dict_name.get(columns[1], columns[1])
            columns[2] = dict_classification.get(columns[2], columns[2])
            
            # Write the updated line to output
            f_out.write('\t'.join(columns) + '\n')

def main():
    parser = argparse.ArgumentParser(description="Process a 3-column tab-separated file using a 6-column dictionary.")
    parser.add_argument("input_file", help="Path to the input 3-column file.")
    parser.add_argument("dictionary_file", help="Path to the 6-column dictionary file.")
    parser.add_argument("output_file", help="Path to the output file.")
    
    args = parser.parse_args()
    
    # Load dictionary
    dict_name, dict_classification = load_dict_file(args.dictionary_file)
    
    # Process input file
    process_input_file(args.input_file, args.output_file, dict_name, dict_classification)

if __name__ == "__main__":
    main()
