# Importation of python modules
import random
import math
import copy
import src.peptide as peptide
from src.moves import ProteinMoveGenerator

class MCSimulation:
    """
    Monte Carlo simulation for protein conformational changes.

    Args:
        amino_acid_info_list (list): List of dictionaries containing amino acid information.

    Attributes:
        amino_acid_info_list (list): List of dictionaries containing amino acid information.

    Methods:
        select_move(all_possible_moves): Select moves with the same key as the first element of the chosen term.
        generate_new_conformation(selected_moves): Generate a new conformation based on selected moves.
        calculate_energy_difference(new_conformation): Calculate the energy difference between original and new conformations.
        calculate_probability_of_move(energy_difference, temperature): Calculate the probability of making a move.
        make_selected_move(new_conformation, probability_move): Make or reject a selected move.
        mc_simulation(temperature, move_mode, nb_iterations): Perform Monte Carlo simulation.
    """
    def __init__(self, amino_acid_info_list):
        """
        Initialize the MCSimulation instance with the provided amino acid information.

        Args:
            amino_acid_info_list (list): List of dictionaries containing amino acid information.
        """
        self.amino_acid_info_list = amino_acid_info_list

    def select_move(self, all_possible_vshd_moves):
        """
        Select moves with the same key as the first element of the chosen term.

        Args:
            all_possible_moves (list): List of possible moves to select from.

        Returns:
            list: List of selected moves.
        """
        # Check if the list of possible moves is empty
        if not all_possible_vshd_moves:
            return []

        # Shuffle the list of possible moves to introduce randomness
        random.shuffle(all_possible_vshd_moves)

        # Get the key of the first element in the chosen term
        first_element_key = next(iter(all_possible_vshd_moves[0]), None)

        # Select moves with the same key as the first element
        matching_moves = [move for move in all_possible_vshd_moves if next(iter(move), None) == first_element_key]

        return matching_moves

    def generate_new_conformation(self, selected_moves):
        """
        Generate a new conformation based on selected moves.

        Args:
            selected_moves (list): List of selected moves to apply to the conformation.

        Returns:
            list: New conformation after applying selected moves.
        """
        # Create a deep copy of the current amino acid info list to modify
        new_conformation = copy.deepcopy(self.amino_acid_info_list)

        # Iterate through the selected moves and update the new conformation
        for move in selected_moves:
            if "new_x" in move and "new_y" in move:
                amino_acid_id = move.get("id")

                # Find the index of the amino acid in the new conformation
                index_to_update = next((i for i, aa in enumerate(new_conformation) if aa["id"] == amino_acid_id), None)

                # Update the position of the amino acid if found
                if index_to_update is not None:
                    new_conformation[index_to_update]["x"] = move["new_x"]
                    new_conformation[index_to_update]["y"] = move["new_y"]

        return new_conformation

    def calculate_energy_difference(self, new_conformation):
        """
        Calculate the energy difference between original and new conformations.

        Args:
            new_conformation (list): New conformation to calculate the energy difference for.

        Returns:
            float: Energy difference between original and new conformations.
        """
        # Create an instance of the AminoAcidSequence class
        sequence = peptide.AminoAcidSequence()

        # Calculate the energy of the original and new conformations
        original_energy = sequence.calculate_energy(self.amino_acid_info_list)
        new_energy = sequence.calculate_energy(new_conformation)

        # Calculate the energy difference
        energy_difference = new_energy - original_energy
        return energy_difference

    def calculate_probability_of_move(self, energy_difference, temperature):
        """
        Calculate the probability of making a move based on energy difference and temperature.

        Args:
            energy_difference (float): Energy difference between original and new conformations.
            temperature (float): Temperature for the simulation.

        Returns:
            float: Probability of making the move.
        """
        # Boltzmann constant
        K_B = 0.0019872

        # Check if the energy difference is non-positive (favorable move)
        if energy_difference <= 0:
            probability = 1
        else:
            # Calculate the probability based on the Metropolis criterion
            probability = math.exp(-energy_difference / (temperature * K_B))
        return probability

    def make_selected_move(self, new_conformation, probability_move):
        """
        Make or reject a selected move based on the calculated probability.

        Args:
            new_conformation (list): New conformation to apply the move to.
            probability_move (float): Probability of making the move.

        Returns:
            list: Updated conformation (either accepted or rejected).
        """
        # Generate a random value between 0 and 1
        random_value = random.random()

        # Check if the random value is less than the calculated probability
        if random_value < probability_move:
            # Accept the move by updating the current conformation
            self.amino_acid_info_list = copy.deepcopy(new_conformation)

        # Return the updated conformation, whether accepted or rejected
        return self.amino_acid_info_list

    def mc_simulation(self, temperature, move_mode, nb_iterations):
        """
        Perform Monte Carlo simulation for protein conformational changes.

        Args:
            temperature (float): Temperature for the simulation.
            move_mode (str): Mode for selecting moves ("VSHD" or "Pull move").
            nb_iterations (int): Number of simulation iterations.

        Returns:
            tuple: A tuple containing the final conformation and a dictionary of energy history.
        """
        # Create an empty dictionary to store energy history
        energy_history = {}

        # Iterate through the specified number of simulation iterations
        for iteration in range(1, nb_iterations + 1):
            sequence = peptide.AminoAcidSequence()

            # Calculate and store the current energy of the conformation
            current_energy = sequence.calculate_energy(self.amino_acid_info_list)
            energy_history[iteration - 1] = current_energy

            # Determine the move mode and generate a new conformation accordingly
            if move_mode == "VSHD":
                all_possible_vshd_moves = ProteinMoveGenerator(
                    self.amino_acid_info_list).generate_all_possible_vshd_moves()
                selected_moves = self.select_move(all_possible_vshd_moves)
                new_conformation = self.generate_new_conformation(selected_moves)
            elif move_mode == "Pull move":
                new_conformation = ProteinMoveGenerator(self.amino_acid_info_list).generate_pull_move()

            # Calculate the energy difference between current and new conformations
            energy_difference = self.calculate_energy_difference(new_conformation)

            # Calculate the probability of making the move based on energy difference and temperature
            probability_move = float(self.calculate_probability_of_move(energy_difference, temperature))

            # Apply or reject the move based on the calculated probability
            self.amino_acid_info_list = self.make_selected_move(new_conformation, probability_move)

            # Calculate and store the current energy again at the final iteration
            if iteration == nb_iterations:
                sequence = peptide.AminoAcidSequence()
                current_energy = sequence.calculate_energy(self.amino_acid_info_list)
                energy_history[iteration] = current_energy

        # Return the final conformation and the energy history
        return self.amino_acid_info_list, energy_history

