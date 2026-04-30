# =========================
# DNA ANALYSIS FUNCTIONS
# =========================

# 1. Load FASTA file
def load_fasta(filename):
    """
    Reads a FASTA file and returns name + DNA sequence
    """
    with open(filename) as f:
        dna=""
        for line in f:
            if line[0]!=">":
                 dna+=line.rstrip()
            else:
                header=line.split()
                name=header[0][1:]  
    return name, dna 

def gc_content(seq):
    """
    Returns GC percentage of DNA sequence
    """
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

# 3. Reverse complement
def reverse_complement(seq):
    """
    Returns reverse complement of DNA
    """
    complement = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
    return "".join(complement[b] for b in reversed(seq))


# 4. Motif search
def find_motif(seq, motif):
    """
    Returns all positions of a motif in DNA
    """
    positions = []

    for i in range(len(seq) - len(motif) + 1):
        if seq[i:i+len(motif)] == motif:
            positions.append(i)

    return positions


# 5. ORF detection
def find_orfs(seq):
    """
    Finds Open Reading Frames (ATG -> stop codon)
    """
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]

    orfs = []

    for i in range(len(seq) - 2):
        if seq[i:i+3] == start_codon:

            for j in range(i+3, len(seq)-2, 3):
                codon = seq[j:j+3]

                if codon in stop_codons:
                    orfs.append((i, j+3))
                    break

    return orfs


# 6. DNA → RNA (transcription)
def transcribe_dna_to_rna(dna):
    """
    Converts DNA into mRNA
    """
    return dna.replace("T", "U")


# 7. RNA → Protein (translation)
codon_table = {
    "AUA":"I","AUC":"I","AUU":"I","AUG":"M",
    "ACA":"T","ACC":"T","ACG":"T","ACU":"T",
    "AAC":"N","AAU":"N","AAA":"K","AAG":"K",
    "AGC":"S","AGU":"S","AGA":"R","AGG":"R",
    "CUA":"L","CUC":"L","CUG":"L","CUU":"L",
    "CCA":"P","CCC":"P","CCG":"P","CCU":"P",
    "CAC":"H","CAU":"H","CAA":"Q","CAG":"Q",
    "CGA":"R","CGC":"R","CGG":"R","CGU":"R",
    "GUA":"V","GUC":"V","GUG":"V","GUU":"V",
    "GCA":"A","GCC":"A","GCG":"A","GCU":"A",
    "GAC":"D","GAU":"D","GAA":"E","GAG":"E",
    "GGA":"G","GGC":"G","GGG":"G","GGU":"G",
    "UCA":"S","UCC":"S","UCG":"S","UCU":"S",
    "UUC":"F","UUU":"F","UUA":"L","UUG":"L",
    "UAC":"Y","UAU":"Y","UAA":"(stop)","UAG":"(stop)",
    "UGC":"C","UGU":"C","UGA":"(stop)","UGG":"W"
}


def translate_rna(rna):
    """
    Converts mRNA into protein sequence
    """
    protein = ""

    for i in range(0, len(rna)-2, 3):
        codon = rna[i:i+3]

        if codon in codon_table:
            aa = codon_table[codon]

            if aa == "_":
                break

            protein += aa

    return protein

def naive_with_rc(p,t):
    appearences=[]
    for i in range(len(t)-len(p)+1):
        for j in range(len(p)):
            match=True
            if t[i+j]!=p[j]:
                match=False
        if match==True:
            break
        appearences.append(i)

    coappearences=[]
    cp=reverse_complement

    for i in range(len(t)-len(cp)+1):
        for j in range(len(cp)):
            match=True
            if t[i+j]!=cp[j]:
                match=False
        if match==True:
            break
        coappearences.append(i)
    return appearences,coappearences

            
        
    

 