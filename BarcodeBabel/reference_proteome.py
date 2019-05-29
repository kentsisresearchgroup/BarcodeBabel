from Bio import SeqIO
from pkg_resources import resource_filename

path_to_reference_proteome_fasta = resource_filename('BarcodeBabel', 'data/20181017_UPR_homo_cRAP_tatfly.fasta')


class ReferenceProteome():
    def __init__(self,
                 path_to_fasta):
        """Basic, will speed this up if needed"""
        self.load_from_fasta_file(path_to_fasta)
        self.path_to_fasta = path_to_fasta

    def load_from_fasta_file(self, path_to_fasta):
        self.reference_sequences = []
        for seq_record in SeqIO.parse(path_to_fasta, "fasta"):
            self.reference_sequences.append(str(seq_record.seq))

    def check_if_in_reference_proteome(self, peptide):
        for reference_seq in self.reference_sequences:
            if peptide in reference_seq:
                return True
        return False


reference_proteome = ReferenceProteome(path_to_fasta=path_to_reference_proteome_fasta)


def not_in_reference_proteome(peptide):
    """Return True only if this peptide sequence is not present in the reference proteome"""
    return not reference_proteome.check_if_in_reference_proteome(peptide)
