#%% This is still a work in progress

#%%
import random
import time

codon_usage_full = {
    'F': ['TTT', 'TTC'], 'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'I': ['ATT', 'ATC', 'ATA'], 'M': ['ATG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'], 'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'], 'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'A': ['GCT', 'GCC', 'GCA', 'GCG'], 'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA'], 'H': ['CAT', 'CAC'],
    'Q': ['CAA', 'CAG'], 'N': ['AAT', 'AAC'],
    'K': ['AAA', 'AAG'], 'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'], 'C': ['TGT', 'TGC'],
    'W': ['TGG'], 'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG']
}

# New protein sequence
protein_sequence_twice_100 = 'MIQCGSSTCAFLIREPIPVWWDIKWMHGNLLLGKKMRNVECMTIRQVITPAGLVMFKTWHWRNLFNSHSLATQVCIFMALAPAQDHMEPTHYVPTRFDGRMIQCGSSTCAFLIREPIPVWWDIKWMHGNLLLGKKMRNVECMTIRQVITPAGLVMFKTWHWRNLFNSHSLATQVCIFMALAPAQDHMEPTHYVPTRFDGR'
protein_sequence_twice_1000 = 'AWMIGVSEYNIRVSYINHPFIVVDVRNWLPWYSTQQMLLIDNLSARDKQTQSTSMCDFTNDHDADNTGGQPDIPHWVKGASQILVSHFWDTTPLFSGWRKYTCQGMFKEQESYQVSHQMWQHNGTCATMLAVQVKWKVCSQGTEPMIICMMECRRELNTEWLPNCKSGTWCGFFSFVRRKGQRTKWYKLKLNRYPMTPKHCDGSKSNMWSDIHIMDSQDHRAKKNYFMNQKMPEKNIYSFEGLVRKFVMGKFMNMFDASDQGTHIYQCDYDTDPQHIENYWNFAMNYWSPTITDQEGQLKKYMALIITKGWVGALFPMQWEAQFHHSWCLLMEQQVCHDKCAVIKQLFASTEGQQYAFVINYTNQFAKCRDEGFWGTQIPWDNELLKSIQCVSGLFEITKIDVVGLARMGCIRADGPQCTWYAINPVPVKVINLLSHNFNDWKIEMFNQSCLSCMRIRQPADMFNWWLDCGVPRMFILYAILIKLWARGLDKTIHIGGTCRVKKMTRYTMVKGRFVSFMGMMVYANNGMVDCCIDLERQGADSNMTGYCHPYGPPNLWHNNTARGKQVTKTNTAKWTSENNCAQNDGVPACNDQSMYSHGWNSEHRLWKLALTEVYFANVLFKMQHAYKDFCEQKHHDGYGFTEFDYSYNYRWHMAKYENDKLNEFVVVQAIQVRNPWRRNWRDGMENIVKYSTYYAYIELGRLLGDRHGVTWSHCLLATMDHFPTCGWITYQKNDLDSFGIDFRVDHMCHHKFRMKPIPVEWCSLQETRGKAMRSHLMLVEYPPWVLRRHNMPMWDMSAATKNPAWIVMQHWQNHIPPVHMKGHGAYMCHCYEEYWWYGALKQRANYIFDWYNLYFGCTICLSRPCLPHAAHSIDPLDWCCYMDEEWQMHMIMLTITSGDPGLCLVDGYLHGTHQKDRVYPNKTDLWWQKEFQLWIGFELTQEWLPMVYRDRFPMFSLWDDMSDAFGWVFFAWNYLTLVVRNTMACYLWPMQMKIANQFAWMIGVSEYNIRVSYINHPFIVVDVRNWLPWYSTQQMLLIDNLSARDKQTQSTSMCDFTNDHDADNTGGQPDIPHWVKGASQILVSHFWDTTPLFSGWRKYTCQGMFKEQESYQVSHQMWQHNGTCATMLAVQVKWKVCSQGTEPMIICMMECRRELNTEWLPNCKSGTWCGFFSFVRRKGQRTKWYKLKLNRYPMTPKHCDGSKSNMWSDIHIMDSQDHRAKKNYFMNQKMPEKNIYSFEGLVRKFVMGKFMNMFDASDQGTHIYQCDYDTDPQHIENYWNFAMNYWSPTITDQEGQLKKYMALIITKGWVGALFPMQWEAQFHHSWCLLMEQQVCHDKCAVIKQLFASTEGQQYAFVINYTNQFAKCRDEGFWGTQIPWDNELLKSIQCVSGLFEITKIDVVGLARMGCIRADGPQCTWYAINPVPVKVINLLSHNFNDWKIEMFNQSCLSCMRIRQPADMFNWWLDCGVPRMFILYAILIKLWARGLDKTIHIGGTCRVKKMTRYTMVKGRFVSFMGMMVYANNGMVDCCIDLERQGADSNMTGYCHPYGPPNLWHNNTARGKQVTKTNTAKWTSENNCAQNDGVPACNDQSMYSHGWNSEHRLWKLALTEVYFANVLFKMQHAYKDFCEQKHHDGYGFTEFDYSYNYRWHMAKYENDKLNEFVVVQAIQVRNPWRRNWRDGMENIVKYSTYYAYIELGRLLGDRHGVTWSHCLLATMDHFPTCGWITYQKNDLDSFGIDFRVDHMCHHKFRMKPIPVEWCSLQETRGKAMRSHLMLVEYPPWVLRRHNMPMWDMSAATKNPAWIVMQHWQNHIPPVHMKGHGAYMCHCYEEYWWYGALKQRANYIFDWYNLYFGCTICLSRPCLPHAAHSIDPLDWCCYMDEEWQMHMIMLTITSGDPGLCLVDGYLHGTHQKDRVYPNKTDLWWQKEFQLWIGFELTQEWLPMVYRDRFPMFSLWDDMSDAFGWVFFAWNYLTLVVRNTMACYLWPMQMKIANQF'

def translate_codon_to_protein(codon_sequence, codon_dict):
    """
    Translate a codon sequence into its corresponding amino acid sequence.

    Args:
        codon_sequence (str): The nucleotide sequence to be translated.
        codon_dict (dict): A dictionary mapping amino acids to their corresponding codons.

    Returns:
        str: The translated amino acid sequence.
    """

    # Create a reverse dictionary where each codon maps to its corresponding amino acid.
    # This is done by iterating over each amino acid and its list of codons in the provided dictionary.
    codons = [codon_sequence[i:i+3] for i in range(0, len(codon_sequence), 3)]
    
    # Initialized list for amino acids
    amino_acids = []

    # Translate each codon in the sequence to its corresponding amino acid using the reverse dictionary.
    reverse_codon_dict = {c: aa for aa, codons_list in codon_dict.items() for c in codons_list}
    for codon in codons:
        amino_acids.append(reverse_codon_dict[codon])

    # Join the list of amino acids into a single string and return.
    return ''.join(amino_acids)

def generate_codon_sequences(sequence, codon_usage, max_alternatives=1000, input_type='codon'):
    """
    Generate alternative codon sequences based on a given amino acid or codon sequence.

    Parameters:
        sequence (str): Input nucleotide or amino acid sequence.
        codon_usage (dict): Dictionary of codon usage. Keys are amnio acids and value are corresponding coding sequences
        max_alternatives (int, optional): Maximum number of alternative sequences to generate. Defaults to 1000.
        input_type (str, optional): Either 'codon' or 'amino_acid' to specify the type of the input sequence.

    Returns:
        list: List of alternative codon sequences.
    """

    # Validate input_type
    if input_type not in ['codon', 'amino_acid']:
        raise ValueError("input_type must be either 'codon' or 'amino_acid'.")

    # Reverse dictionary to map codons to amino acids
    reverse_codon_dict = {codon: aa for aa, codons in codon_usage.items() for codon in codons}

    if input_type == 'codon':
        # Split the input nucleotide sequence into codons
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
        
        # Check for invalid codons
        for codon in codons:
            if codon not in reverse_codon_dict:
                raise ValueError(f"Invalid codon: {codon}")
        
        # Check for stop codons in the middle of the sequence
        if any(reverse_codon_dict[codon] == '*' for codon in codons[:-1]):
            raise ValueError("Stop codon found in the middle of the sequence.")
        
        # Translate the nucleotide sequence into its corresponding amino acid sequence
        amino_acid_sequence = [reverse_codon_dict[codon] for codon in codons]

    else:  # input_type == 'amino_acid'
        amino_acid_sequence = list(sequence)

        # Check for invalid amino acids
        valid_amino_acids = set(codon_usage.keys())
        for aa in amino_acid_sequence:
            if aa not in valid_amino_acids:
                raise ValueError(f"Invalid amino acid: {aa}")

        # Check for stop codons in the middle of the sequence
        if '*' in amino_acid_sequence[:-1]:
            raise ValueError("Stop codon found in the middle of the sequence.")
    
    # Generate alternative codon sequences based on the amino acid sequence
    alternative_sequences = []
    for _ in range(max_alternatives):
        codon_sequence = ''.join(random.choice(codon_usage[aa]) for aa in amino_acid_sequence)
        alternative_sequences.append(codon_sequence)

    # Remove duplicates and sort by GC content
    alternative_sequences = list(set(alternative_sequences))
    alternative_sequences.sort(key=lambda seq: abs(calculate_gc_content(seq) - 0.5))

    # Returns the string 'alternative_sequences'
    return alternative_sequences

def reverse_complement(sequence):
    """
    Convert a DNA sequence into its reverse complement.
      
    Parameters:
    - sequence (str): The DNA sequence to be reverse complemented.

    Returns:
    - str: The reverse complemented sequence.
    """

    # Create a dictionary that maps each DNA base to its complementary base.
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    # For each base in the reversed sequence, find its complement.
    # The 'reversed()' function iterates over the sequence in reverse order.
    # The 'complement[base]' lookup gives the complementary base for each base in the sequence.
    # Finally, the ''.join() function combines the complementary bases to form the reverse complemented sequence.
    return ''.join(complement[base] for base in reversed(sequence))

# Define the function to calculate the GC content of a sequence
def calculate_gc_content(sequence):
    """
    This function calculates the GC content of a nucleotide sequence, defined as the proportion of G's and C's
    in the sequence.

    Parameters:
    - sequence (str): A string representing the nucleotide sequence.

    Returns:
    - float: The GC content of the sequence.
    
    Note:
    - Is faster then gc_content_statistics() 
    """

    # Count the number of 'G' bases in the sequence
    g_content = sequence.count('G')
    
    # Count the number of 'C' bases in the sequence
    c_content = sequence.count('C')

    # Calculate the sum of G's and C's
    total_gc_content = g_content + c_content

    # Calculate the proportion of G's and C's by dividing the sum by the total length of the sequence
    return total_gc_content / len(sequence)

def gc_content_statistics(sequence, window_size=30, target_range=(40, 60)):
    """
    Calculate the overall GC content and the number of windows outside the target GC content range.

    Parameters:
    - sequence (str): The DNA sequence to analyze.
    - window_size (int): The size of the window to consider when calculating local GC content. Default is 30.
    - target_range (tuple of (float, float)): The minimum and maximum acceptable GC content percentages. Default is (40, 60).

    Returns:
    - list of [float, int]: The overall GC content of the entire sequence and the number of windows that fall outside the target range.

     Note:
    - Is slower then calculate_gc_statistics() 
    """
    
    # Calculate overall GC content of the entire sequence
    gc_content = calculate_gc_content(sequence)
    
    # Initialize count for windows outside target range
    outliers = 0
    
    # Loop through the DNA sequence using a sliding window approach
    for i in range(len(sequence) - window_size + 1):
        # Extract the current window from the sequence
        window = sequence[i:i+window_size]
        
        # Calculate the GC content for the current window
        window_gc_content = calculate_gc_content(window) * 100
        
        # Check if the window's GC content falls outside the target range
        if window_gc_content < target_range[0] or window_gc_content > target_range[1]:
            outliers += 1

    return [gc_content, outliers]

def find_repeats(sequence, min_length=1, max_length=None):
    """
    Identify repeated substrings in a given sequence.

    Parameters:
    - sequence (str): The input sequence to search for repeated substrings.
    - min_length (int, optional): The minimum length of the repeated substrings. Defaults to 1.
    - max_length (int, optional): The maximum length of the repeated substrings. If not provided, it defaults to the length of the sequence.

    Returns:
    - list of tuples: Each tuple contains a repeated substring and its positions (1-based) in the sequence.
    - float: The percentage of positions in the sequence that are covered by repeated substrings.
    """
    
    # Identify the length of the input sequence
    length = len(sequence)

    # Dictionary to store repeated substrings and their positions
    repeats = {}

    # If max_length is not provided or exceeds the sequence length, set it to the sequence length
    if max_length is None:
        max_length = length

    # Loop over the range of possible substring lengths (from min_length to max_length)
    for l in range(min_length, max_length+1):
        # Loop over the starting positions of the substrings
        for i in range(length - l + 1):
            # Extract the substring of length l starting from position i
            repeat = sequence[i:i+l]
            
            # Check if the substring is already in the repeats dictionary
            if repeat in repeats:
                # If yes, append the new position (1-based) to its list of positions
                repeats[repeat].append(i + 1)
            else:
                # If not, add the substring to the dictionary with its position
                repeats[repeat] = [i + 1]

    # Filter out the substrings that appear only once in the sequence
    repeats = {repeat: positions for repeat, positions in repeats.items() if len(positions) > 1}

    # Create a list of tuples where each tuple is (substring, positions)
    repeats_with_positions = [(substring, positions) for substring, positions in repeats.items()]

    # Sort the list of tuples in descending order based on the length of the substring
    repeats_with_positions.sort(key=lambda x: len(x[0]))

    # Calculate the coverage of repeated substrings in the sequence
    # A set is used to ensure unique positions and avoid double counting
    covered_positions = set()
    for repeat, positions in repeats_with_positions:
        for pos in positions:
            # Update the set with all positions covered by the current repeated substring
            covered_positions.update(range(pos - 1, pos - 1 + len(repeat)))

    # Calculate the percentage of positions in the sequence that are covered by repeated substrings
    coverage = (len(covered_positions) / len(sequence)) * 100 

    # Return the list of repeated substrings with their positions and the coverage
    return repeats_with_positions, coverage

def find_hairpins(seq, min_stem_size=10, max_stem_size=20, min_loop_size=3):
    """
    Annotate potential DNA hairpins in the sequence and return a detailed description of each hairpin.

    Parameters:
    - seq (str): The DNA sequence to search for hairpins.
    - min_stem_size (int): The minimum size of the hairpin stem.
    - max_stem_size (int): The maximum size of the hairpin stem.
    - min_loop_size (int): The minimum size of the loop between hairpin stems.

    Returns:
    - list: A list of annotated hairpins.
    """

    # Initializes the list for the hairpins
    hairpins = []

    # Determine the range of stem sizes to search, from the minimum to the lesser of
    # the maximum or half the sequence length. The list is then reversed to prioritize
    # longer stems.
    stem_sizes = list(range(min_stem_size, min(max_stem_size + 1, len(seq) // 2 + 1)))[::-1]

    # Determine the range of loop sizes, from the minimum to the maximum possible 
    # considering the length of the sequence and the minimum stem size.
    loop_sizes = list(range(min_loop_size, len(seq) - 2 * min_stem_size + 1))

    # Iterate through possible stem sizes
    for stem_size in stem_sizes:
        # For each stem size, iterate through the sequence to find potential stems
        for i in range(len(seq) - 2 * stem_size):
            stem = seq[i:i + stem_size]
            rev_comp = reverse_complement(stem)

            # For each potential stem, iterate through possible loop sizes
            for loop_size in loop_sizes:
                # Calculate the starting index for searching the reverse complement
                start_idx = i + stem_size + loop_size
                reverse_start_idx = seq[start_idx:].find(rev_comp)

                # If a reverse complement is found, add the hairpin details to the annotations
                if reverse_start_idx != -1:
                    # Adjusting reverse_start_idx to be relative to the entire sequence
                    reverse_start_idx += start_idx
                    # Appending the stem, its reverse complement, and the 1-based indices of the forward and reverse stems
                    hairpins.append((stem, rev_comp, i + 1, reverse_start_idx + 1))
      
                    #Break to avoid redundant findings for the same stem and move on to the next stem
                    break
                    
    return hairpins

def get_codons_from_full_sequence(full_sequence, subsequence, position):
    """
    Extracts the codon sequence (or sequences) from the full DNA sequence that encompasses a given subsequence and its position.

    Parameters:
    - full_sequence (str): The full DNA sequence from which codons will be extracted.
    - subsequence (str): The subsequence whose encompassing codons are to be identified.
    - position (int): The 1-based starting position of the subsequence in the full_sequence.

    Returns:
    - str: The codon sequence that encompasses the given subsequence.
    - int: The 1-based starting position of the codon sequence in the full DNA sequence.
    """

    # Calculate the start position of the codon sequence.
    codon_start = (position - 1) // 3 * 3
    
    # Calculate the end position of the codon sequence.
    codon_end = codon_start + (len(subsequence) + 2) // 3 * 3

    # Extract the codon sequence from the full sequence.
    codon_sequence = full_sequence[codon_start:codon_end]

    # Convert from 0-based to 1-based index
    codon_start = codon_start + 1

    # Return the codon sequence, its start position, and its end position.
    return codon_sequence, codon_start

def validate_potential_replacement(sequence, repeat, alternative_codon_sequence, position):
    """
    Validates if replacing a repeated sequence with an alternative codon sequence in the main sequence introduces new repeats.
    
    Parameters:
    - sequence (str): The full DNA sequence.
    - repeat (str): The repeated sequence that is intended to be replaced.
    - alternative_codon_sequence (str): The sequence proposed to replace the repeat.
    - position (int): The 1-based starting position of the repeated sequence in the full DNA sequence.

    Returns:
    - bool: True if the replacement does not introduce new repeats, False otherwise.
    - str: The updated sequence after the replacement if no new repeats are introduced, else the original sequence.
    """

    # Adjust the position to be 0-based.
    position = position - 1

    # Locate the alternative_codon_sequence in the sequence by position
    start = position
    end = position + len(alternative_codon_sequence)

    # Get indexes of the basepairs left and right of the stretch in the sequence if possible
    flanking_bp = len(alternative_codon_sequence) - 1
    left_index = max(0, start - flanking_bp)
    right_index = min(len(sequence), end + flanking_bp)

    # Build the 'input_stretch plus' the basepairs left and right that will be named for instance 'extended_stretch'
    left_stretch = sequence[left_index:start]
    right_stretch = sequence[end:right_index]
    extended_stretch = left_stretch + alternative_codon_sequence + right_stretch

    # Update sequence with input_stretch
    updated_sequence = sequence[:start] + alternative_codon_sequence + sequence[end:]
    updated_sequence_rv = reverse_complement(updated_sequence)
     
    # Start with a window size equal to the length of the repeat and go up to the size of the check region
    for window_size in range(len(repeat), len(extended_stretch) + 1):

        # Sliding window approach to check for repeats in the check_region
        for i in range(len(extended_stretch) - window_size + 1):
            sub_stretch = extended_stretch[i:i+window_size]

            # Check the sub_stretch against the entire updated_sequence
            if updated_sequence.count(sub_stretch) > 1 or updated_sequence_rv.count(sub_stretch) > 1:
                return False, sequence 
            
    # Return True and the updated sequence if no sub_stretch is found more than once
    return True, updated_sequence  

def optimize_gc_content(sequence, codon_usage, gc_window_size=30, gc_step_size=15, gc_target_content_range=(30,70), gc_iterations=10, max_alternatives=500):
    """
    Iteratively optimizes the GC content of a sequence by replacing windows of the sequence with sequences
    that have a GC content within a specified range.

    Parameters:
        sequence (str): The original nucleotide sequence.
        codon_usage (dict): Dictionary mapping amino acids to their corresponding codons.
        gc_window_size (int): The size of the window to consider when optimizing GC content.
        gc_step_size (int): The number of bases to move the window at each step.
        gc_target_content_range (tuple): A tuple specifying the lower and upper bounds of the target GC content.
        gc_iterations (int): The number of iterations that of the algoritm

    Returns:
        str: The optimized nucleotide sequence.
    """
    
    # Generate protein sequence from original input codon sequence
    original_protein = translate_codon_to_protein(sequence, codon_usage)

    # Convert the sequence string into a list of nucleotide characters for ease of modification
    sequence_list = list(sequence)
    
    # Compute the total length of the sequence
    sequence_length = len(sequence_list)
    
    # Compute the indices for each window
    window_indices = [(i, min(i + gc_window_size, sequence_length)) for i in range(0, sequence_length, gc_step_size)]
   
    if len(sequence_list) - window_indices[-1][0] < gc_window_size:
        window_indices[-1] = (window_indices[-2][0], sequence_length)

    # Perform optimization over multiple iterations
    for loop in range(gc_iterations):

        # Generate protein sequence from 'sequence' and compare it to 'original_protein'
        current_protein = translate_codon_to_protein(sequence, codon_usage)
        if current_protein != original_protein:
            raise ValueError(f"Warning: The translated protein sequence after loop {loop} does not match the original protein sequence!")
            break

        # Loop over each window defined by the start and end indices
        for start, end in window_indices:
            
            # Extract the current window from the sequence
            window = "".join(sequence_list[start:end])
            
            # Calculate the GC content of the current window
            gc_content = calculate_gc_content(window)
            
            # If the GC content of the window is outside the target range, generate candidate sequences
            if not gc_target_content_range[0] <= gc_content <= gc_target_content_range[1]:
                
                # Generate alternative codon sequences for the current window
                candidate_sequences = generate_codon_sequences(window, codon_usage, max_alternatives=max_alternatives, input_type='codon')
                
                # Filter the candidate sequences to only include those within the target GC content range
                valid_candidate_sequences = [candidate_sequence for candidate_sequence in candidate_sequences if gc_target_content_range[0] <= (calculate_gc_content(candidate_sequence) * 100) <= gc_target_content_range[1]]
                
                # If there are no valid candidates, continue to the next window
                if not valid_candidate_sequences:
                    continue        

                # From the valid sequences, select the one closest to 50% GC content
                best_candidate = min(valid_candidate_sequences, key=lambda seq: abs(calculate_gc_content(seq) - 0.5))
                
                # Replace the current window in the sequence with the best candidate
                sequence_list[start:end] = list(best_candidate)
        
        # Convert the sequence list back into a string
        optimized_sequence = "".join(sequence_list)
        # Optimized sequence
        sequence = optimized_sequence if optimized_sequence else sequence
 
        # Get more statistics of the GC-content
        statistics = gc_content_statistics(sequence, gc_window_size, gc_target_content_range)

        # Make statistics a procentage
        overal_content = statistics[0] * 100

        # Print out progress
        print(f"- Iteration {loop+1} - Overal GC%: {overal_content:.2f}% - {gc_window_size}bp streches <{gc_target_content_range[0]}% and >{gc_target_content_range[1]}% :{statistics[1]}")
       
        # Check if GC-content is fully optimized
        fully_optimized = False
        if(round(statistics[0], 2)  == 0.50 and statistics[1] == 0):
            
            # Generate protein sequence from 'sequence' and compare it to 'original_protein'
            current_protein = translate_codon_to_protein(sequence, codon_usage)
            if current_protein != original_protein:
                raise ValueError(f"Warning: The protein sequence after loop {loop} does not match the original protein sequence!")
                break

            # Print out that the GC-content is maximally optimized
            print("- The GC-content is maximally optimized!")

            # Return the sequence
            return sequence

    # Return the optimized sequence, or the original sequence if no optimization was performed
    print("- Finished: The GC-content is maximally optimized!")

    # Generate protein sequence from 'sequence' and compare it to 'original_protein'
    current_protein = translate_codon_to_protein(sequence, codon_usage)
    if current_protein != original_protein:
        raise ValueError(f"Warning: The translated protein sequence after loop {loop} does not match the original protein sequence!")
    
    # Return the sequence
    return sequence

def reduce_repeats(sequence, codon_usage, iterations=1, min_length_repeats=9, max_length_repeats=20, max_number_of_alternatives_per_repeat=1000):
    """
    Attempts to reduce the number of repeated sequences in a given DNA sequence by replacing them with alternative codon sequences that encode the same amino acids.
    
    Parameters:
    - sequence (str): The initial DNA sequence to be optimized.
    - codon_usage (dict): A dictionary mapping amino acids to their corresponding codons.
    - iterations (int, optional): The number of times the repeat reduction process should be repeated. Defaults to 1.
    - min_length_repeats (int, optional): The minimum length of repeats to consider. Defaults to 9.
    - max_length_repeats (int, optional): The maximum length of repeats to consider. Defaults to 20.
    - max_number_of_alternatives_per_repeat (int, optional): The maximum number of alternative codon sequences to generate for each repeat. Defaults to 1000.
    
    Returns:
    - str: The optimized DNA sequence with reduced repeats.

    Note:
    - Increasing to number for max_length_repeats results in longer computation time, but only marginally improves the outcome. It is better to increase the iterations.
    - Increasing to number of max_number_of_alternatives results in longer computation time, but only marginally improves the outcome. It is better to increase the iterations.
    """

    # Generate protein sequence from original input codon sequence
    original_protein = translate_codon_to_protein(sequence, codon_usage)

    
    for loop in range(iterations):
        
        # Get the repeats
        repeats, coverage = find_repeats(sequence, min_length=min_length_repeats, max_length=max_length_repeats)

        # Print the progress
        print(f"- Iteration {loop+1} - Coverage >{min_length_repeats-1}bp repeats in sequence: {coverage:.2f}%")

        for current_repeat in range(len(repeats)):
         
            # Store the output of repeats in different variable for readibility
            repeat = repeats[current_repeat][0]
            positions = repeats[current_repeat][1]
            
            # If the repeat does not exist in the sequence go the next repeat
            if repeat not in sequence:
                continue

            # Loop trough the different occurences of the repeat and resolve them unit only once occurrence is left
            for i in range(len(positions)):
                
                # Check if the repeating sequence is not repeating anymore
                if sequence.count(repeat) <= 1:
                    break

                # Select the correct position of the repeat
                position = positions[i]
               
                # Generate codon sequence for the current repeat
                codons, codon_position = get_codons_from_full_sequence(sequence, repeat, position)
                
                # Generate alternative codon sequences
                alternative_codons_sequences = generate_codon_sequences(codons, codon_usage, max_alternatives=max_number_of_alternatives_per_repeat)
             
                # Loop through the alternative codon sequences
                found_valid_replacement = False
                for alternative_codon_sequence in alternative_codons_sequences:
                    
                    # Validate the curren alternative codon sequence
                    valid, updated_sequence = validate_potential_replacement(sequence, repeat, alternative_codon_sequence, codon_position)
                    
                    # If the proposed sequence resolves the repeat, accept the change
                    if valid:
                        
                        # Update the current sequence
                        sequence = updated_sequence
                        
                        # Exit the current loop as soon as a valid replacement is found
                        break  
           
        # Generate protein sequence from 'sequence' and compare it to 'original_protein'
        current_protein = translate_codon_to_protein(sequence, codon_usage)

        # In case the nucleotide sequence does not translate to the same protein anymore 
        if current_protein != original_protein:
            raise ValueError(f"Warning: The translated protein sequence after loop {loop+1} does not match the original protein sequence!")
        
        if len(repeats) == 0:     
            print(f"- Finished: All >{min_length_repeats-1}bp repeats are removed!")
            break

    return sequence

def reduce_hairpins(sequence, codon_usage, iterations=5, min_stem_size=10, max_stem_size=20, min_loop_size=3, max_number_of_alternatives=1000):
    """
    Attempts to reduce the number of hairpin structures in a given DNA sequence by replacing them with alternative codon sequences that encode the same amino acids.

    Parameters:
    - sequence (str): The initial DNA sequence to be optimized.
    - codon_usage (dict): A dictionary mapping amino acids to their corresponding codons.
    - iterations (int, optional): The number of times the hairpin reduction process should be repeated. Defaults to 5.
    - min_stem_size (int, optional): The minimum size of the hairpin stem to consider. Defaults to 10.
    - max_stem_size (int, optional): The maximum size of the hairpin stem to consider. Defaults to 20.
    - min_loop_size (int, optional): The minimum size of the hairpin loop. Defaults to 3.
    - max_number_of_alternatives (int, optional): The maximum number of alternative codon sequences to generate for each hairpin stem. Defaults to 1000.
    
    Returns:
    - str: The optimized DNA sequence with reduced hairpin structures.

    Notes:
    - Increasing to number of max_stem_size results in longer computation time, but only marginally improves the outcome. It is better to increase the iterations.
    - Increasing to number of max_number_of_alternatives results in longer computation time, but only marginally improves the outcome. It is better to increase the iterations.
    """

    # Generate protein sequence from original input codon sequence
    original_protein = translate_codon_to_protein(sequence, codon_usage)

    for loop in range(iterations):

        # Generate protein sequence from 'sequence' and compare it to 'original_protein'
        current_protein = translate_codon_to_protein(sequence, codon_usage)
        if current_protein != original_protein:
            raise ValueError(f"Warning: The translated protein sequence after loop {loop} does not match the original protein sequence!")
            break

        # Get the hairpins
        hairpins = find_hairpins(sequence, min_stem_size=min_stem_size, max_stem_size=max_stem_size, min_loop_size=min_loop_size)
        
        print(f"- Iteration {loop+1} - Number of hairpins: {len(hairpins)}")

        # If there are no hairpins anymore, break the loop
        if len(hairpins) == 0:
            break

        for hairpin in hairpins:

            # Reassign variable for readability 
            stem_fw, stem_rv, position_fw, position_rv = hairpin

            # Generate codon sequence for the current stem
            codons_fw, codon_position_fw = get_codons_from_full_sequence(sequence, stem_fw, position_fw)
                
            # Generate alternative codon sequences
            alternative_codons_sequences_fw = generate_codon_sequences(codons_fw, codon_usage, max_alternatives=max_number_of_alternatives)
            
            # Loop through the alternative codon sequences
            found_valid_replacement = False
            for alternative_codon_sequence in alternative_codons_sequences_fw:
                    
                # Validate the current alternative codon sequence
                valid, updated_sequence = validate_potential_replacement(sequence, stem_fw, alternative_codon_sequence, codon_position_fw)
                    
                # If the proposed sequence resolves the hairpin, accept the change
                if valid:
                    sequence = updated_sequence
                    found_valid_replacement = True

                    # Exit the current loop as soon as a valid replacement is found
                    break  

            # If no valid direct replacement was found, try the reverse complement of the stem
            if not found_valid_replacement:
                codons_rv, codon_position_rv = get_codons_from_full_sequence(sequence, stem_rv, position_rv)
                alternative_codons_sequences_rv = generate_codon_sequences(codons_rv, codon_usage, max_alternatives=max_number_of_alternatives)

                # Loop through the alternative codon sequences
                for alternative_codon_sequence in alternative_codons_sequences_rv:
                    
                    # Validate the current alternative codon sequence
                    valid, updated_sequence = validate_potential_replacement(sequence, stem_rv, alternative_codon_sequence, codon_position_rv)
                    
                    # If the proposed sequence resolves the hairpin, accept the change
                    if valid:
                        sequence = updated_sequence
                        break  # Exit the current loop as soon as a valid replacement is found
          
    # Generate protein sequence from 'sequence' and compare it to 'original_protein'
    current_protein = translate_codon_to_protein(sequence, codon_usage)
    if current_protein != original_protein:
        raise ValueError(f"Warning: The translated protein sequence after loop {loop+1} does not match the original protein sequence!")
    
    print(f"- Finished: All hairpins with stem size of >{min_stem_size - 1}bp are removed!")


    return sequence
  
print("STEP 1: Translating the protein sequence to a nucleotide sequence...")
codon_sequence = generate_codon_sequences(protein_sequence_twice_1000, codon_usage_full, max_alternatives=5000, input_type='amino_acid')
print("")

# Give a general overview of properties of the starting sequence

print("STEP 2: Changing the codon-sequence to optimize the GC-content...")
codon_sequence = optimize_gc_content(codon_sequence[0], codon_usage_full, gc_window_size=30, gc_step_size=6, gc_target_content_range=(30,70), gc_iterations=10, max_alternatives=1000)
print("")

print("STEP 3: Changing the codon-sequence to remove  repeats...")
codon_sequence = reduce_repeats(codon_sequence, codon_usage_full, iterations=10, min_length_repeats=12, max_length_repeats=15, max_number_of_alternatives_per_repeat=1000)
print(codon_sequence)
print("")

