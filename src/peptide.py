# Importation of python modules
from Bio import SeqIO
import random

class AminoAcidSequence:
    """
    Represents a sequence of amino acids and provides methods for analyzing and modifying it.
    """

    def __init__(self):
        """
        Initializes an empty AminoAcidSequence object.
        """

    def read_fasta_to_list(self, fasta_file):
        """
        Reads amino acid sequences from a FASTA file and returns them as a list.

        Args:
            fasta_file (str): The path to the FASTA file containing amino acid sequences.

        Returns:
            list: A list of amino acids extracted from the FASTA file.

        Note:
            This method uses the Biopython library's SeqIO module to parse the FASTA file.
        """
        amino_acids = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            amino_acids.extend(list(record.seq))
        return amino_acids

    def replace_amino_acids(self, amino_acids):
        """
        Replaces amino acids in a given sequence with their corresponding polarity values.

        Args:
            amino_acids (list): A list of amino acids.

        Returns:
            list: A list of amino acids with hydrophobic (H) and polar (P) replacements.

        Note:
            Hydrophobic amino acids: "AILMFVPGWC"
            Polar amino acids: "RNDQHEKSTY"
            Other amino acids remain unchanged.
        """
        # Define sets of hydrophobic and polar amino acids
        hydrophobic = "AILMFVPGWC"
        polar = "RNDQHEKSTY"

        # Initialize an empty list to store the replaced amino acids
        replaced_sequence = []

        # Iterate through the input amino acids and replace them based on polarity
        for aa in amino_acids:
            if aa in hydrophobic:
                replaced_sequence.append('H')  # Replace hydrophobic amino acids with 'H'
            elif aa in polar:
                replaced_sequence.append('P')  # Replace polar amino acids with 'P'
            else:
                replaced_sequence.append(aa)  # Keep other amino acids unchanged

        # Return the list of amino acids with replacements
        return replaced_sequence

    def create_amino_acid_info_list(self, modified_sequence):
        """
        Creates a list of dictionaries representing amino acids' information and positions.

        Args:
            modified_sequence (list): A list of amino acids with hydrophobic (H) and polar (P) replacements.

        Returns:
            list: A list of dictionaries, each containing amino acid information.

        Note:
            This method randomly arranges amino acids in a grid while ensuring no two amino acids occupy the same position.
            The resulting list contains dictionaries with keys 'id', 'x', 'y', and 'polarity'.
        """
        while True:
            amino_acid_info_list = []
            x, y = 0, 0  # Initial coordinates

            for index, aa in enumerate(modified_sequence, start=1):
                # Randomly select a direction (up, down, left, or right)
                direction = random.choice([(0, 1), (0, -1), (1, 0), (-1, 0)])

                # Calculate new coordinates based on the selected direction
                new_x, new_y = x + direction[0], y + direction[1]

                # Check if at least one direction is available
                if (x, y + 1) in [(info["x"], info["y"]) for info in amino_acid_info_list] and \
                        (x, y - 1) in [(info["x"], info["y"]) for info in amino_acid_info_list] and \
                        (x + 1, y) in [(info["x"], info["y"]) for info in amino_acid_info_list] and \
                        (x - 1, y) in [(info["x"], info["y"]) for info in amino_acid_info_list]:
                    continue

                # Ensure no two amino acids have the same position
                while (new_x, new_y) in [(info["x"], info["y"]) for info in amino_acid_info_list]:
                    direction = random.choice([(0, 1), (0, -1), (1, 0), (-1, 0)])
                    new_x, new_y = x + direction[0], y + direction[1]

                # Update the coordinates
                x, y = new_x, new_y

                amino_acid_info = {
                    "id": index,
                    "x": x,
                    "y": y,
                    "polarity": aa
                }
                amino_acid_info_list.append(amino_acid_info)

            if len(amino_acid_info_list) != len(modified_sequence):
                continue
            else:
                return amino_acid_info_list

    def calculate_energy(self, amino_acid_info_list):
        """
        Calculates the energy of the amino acid sequence based on hydrophobic interactions.

        Args:
            amino_acid_info_list (list): A list of dictionaries representing amino acid information.

        Returns:
            float: The calculated energy of the amino acid sequence.

        Note:
            Energy is calculated based on the interaction of hydrophobic amino acids (H) with neighboring amino acids.
            If neighboring hydrophobic amino acids are not adjacent to each other, energy is reduced by 0.5.
        """
        energy = 0

        for aa in amino_acid_info_list:
            if aa["polarity"] == "H":
                x, y = aa["x"], aa["y"]

                # Check the neighborhood (up, down, left, right)
                neighbors = [
                    (x, y + 1),  # Up
                    (x, y - 1),  # Down
                    (x + 1, y),  # Right
                    (x - 1, y)  # Left
                ]

                # Implement energy if neighbor H is not adjacent to aa
                for aa_com in amino_acid_info_list:
                    if aa_com["polarity"] == 'H':
                        if (aa_com["x"], aa_com["y"]) in neighbors and aa_com["id"] != aa["id"] + 1 and \
                                aa_com["id"] != (aa["id"] - 1):
                            energy -= 0.5
        return energy
