from src import peptide
import copy
import random
class ProteinMoveGenerator:
    """
    A class for generating possible protein conformation moves.

    Attributes:
        amino_acid_info_list (list): List of dictionaries containing information about amino acid positions.
    """

    def __init__(self, amino_acid_info_list):
        """
        Initialize the ProteinMoveGenerator object.

        Args:
            amino_acid_info_list (list): List of dictionaries with amino acid information.
        """
        self.amino_acid_info_list = amino_acid_info_list

    def calculate_vector(self, point1, point2):
        """
        Calculate the vector between two points.

        Args:
            point1 (dict): Dictionary containing x and y coordinates of the first point.
            point2 (dict): Dictionary containing x and y coordinates of the second point.

        Returns:
            tuple: Tuple containing the x and y components of the vector.
        """
        return point2["x"] - point1["x"], point2["y"] - point1["y"]

    def get_neighbors(self, amino_acid):
        """
        Get the neighboring positions of an amino acid.

        Args:
            amino_acid (dict): Dictionary containing x and y coordinates of the amino acid.

        Returns:
            list: List of neighboring positions.
        """
        x, y = amino_acid["x"], amino_acid["y"]
        return [
            (x, y + 1),  # Up
            (x, y - 1),  # Down
            (x + 1, y),  # Right
            (x - 1, y)   # Left
        ]

    def find_free_positions(self, neighbors):
        """
        Find free positions among neighboring positions.

        Args:
            neighbors (list): List of neighboring positions.

        Returns:
            list: List of free positions.
        """
        return [(x, y) for x, y in neighbors if (x, y) not in [(info["x"], info["y"]) for info in self.amino_acid_info_list]]

    def generate_moves(self, free_positions, amino_acid_id):
        """
        Generate move dictionaries for free positions.

        Args:
            free_positions (list): List of free positions.
            amino_acid_id (int): ID of the amino acid for which moves are generated.

        Returns:
            list: List of move dictionaries.
        """
        moves = []

        for i, (x, y) in enumerate(free_positions, start=1):
            moves.append({
                f"end_move_{i}": i,
                "id": amino_acid_id,
                "new_x": x,
                "new_y": y
            })

        return moves

    def generate_end_moves(self):
        """
        Generate possible "end moves" for the protein conformation.

        Returns:
            list: List of dictionaries representing possible end moves.
        """
        # Initialize an empty list to store possible moves
        possible_moves = []

        # Find the amino acid with ID 2 in the amino acid info list
        amino_acid_id_2 = next((aa for aa in self.amino_acid_info_list if aa["id"] == 2), None)

        # Get the penultimate amino acid in the amino acid info list
        penultimate_amino_acid = self.amino_acid_info_list[-2]

        # Find neighbors of amino acid with ID 2 and the penultimate amino acid
        amino_acid_id_2_neighbors = self.get_neighbors(amino_acid_id_2)
        penultimate_amino_acid_neighbors = self.get_neighbors(penultimate_amino_acid)

        # Find free positions for amino acid with ID 2 and the penultimate amino acid
        free_positions_id_2 = self.find_free_positions(amino_acid_id_2_neighbors)
        free_positions_penultimate = self.find_free_positions(penultimate_amino_acid_neighbors)

        # Generate moves for amino acid with ID 2 and the last amino acid in the sequence
        possible_moves.extend(self.generate_moves(free_positions_id_2, 1))
        possible_moves.extend(self.generate_moves(free_positions_penultimate, self.amino_acid_info_list[-1]["id"]))

        # Return the list of possible end moves
        return possible_moves

    def generate_corner_moves(self):
        """
        Generate possible "corner moves" for the protein conformation.

        Returns:
            list: List of dictionaries representing possible corner moves.
        """
        # Initialize an empty list to store possible corner moves
        possible_moves = []

        # Initialize a counter for corner moves
        corner_move_counter = 1

        # Iterate through the amino acid info list, excluding the last two amino acids
        for i in range(0, len(self.amino_acid_info_list) - 2):
            # Get information about the current amino acid and its neighbors
            amino_acid_i = self.amino_acid_info_list[i]
            amino_acid_i_plus_1 = self.amino_acid_info_list[i + 1]
            amino_acid_i_plus_2 = self.amino_acid_info_list[i + 2]

            # Calculate vectors AB and BC
            AB_x, AB_y = self.calculate_vector(amino_acid_i, amino_acid_i_plus_1)
            BC_x, BC_y = self.calculate_vector(amino_acid_i_plus_1, amino_acid_i_plus_2)

            # Calculate the dot product of vectors AB and BC
            dot_product = AB_x * BC_x + AB_y * BC_y

            # Check if vectors AB and BC are perpendicular (dot product is 0)
            if dot_product == 0:
                # Calculate the coordinates of the corner
                x_corner = amino_acid_i["x"] + BC_x
                y_corner = amino_acid_i["y"] + BC_y

                # Check if the corner position is occupied by another amino acid
                is_occupied = any(info["x"] == x_corner and info["y"] == y_corner for info in self.amino_acid_info_list)

                # If the corner position is not occupied, add it to possible_moves
                if not is_occupied:
                    possible_moves.append({
                        f"corner_move_{corner_move_counter}": corner_move_counter,
                        "id": amino_acid_i_plus_1["id"],
                        "new_x": x_corner,
                        "new_y": y_corner
                    })

                    # Increment the corner move counter
                    corner_move_counter += 1

        # Return the list of possible corner moves
        return possible_moves

    def generate_crankshaft_moves(self):
        """
        Generate possible "crankshaft moves" for the protein conformation.

        Returns:
            list: List of dictionaries representing possible crankshaft moves.
        """
        possible_moves = []
        move_number = 1  # Initialize move number

        # Iterate from the first amino acid to the n-4 amino acid
        for i in range(len(self.amino_acid_info_list) - 4):
            amino_acid_i = self.amino_acid_info_list[i]
            amino_acid_i_plus_3 = self.amino_acid_info_list[i + 3]

            # Calculate the distance between amino acid i and i+3
            distance = abs(amino_acid_i["x"] - amino_acid_i_plus_3["x"]) + abs(amino_acid_i["y"] - amino_acid_i_plus_3["y"])

            if distance == 1:
                # Determine the coordinates of potential positions
                if amino_acid_i["x"] == amino_acid_i_plus_3["x"]:
                    positions = [
                        {"x": amino_acid_i["x"] + 1, "y": amino_acid_i["y"]},
                        {"x": amino_acid_i["x"] - 1, "y": amino_acid_i["y"]},
                        {"x": amino_acid_i_plus_3["x"] + 1, "y": amino_acid_i_plus_3["y"]},
                        {"x": amino_acid_i_plus_3["x"] - 1, "y": amino_acid_i_plus_3["y"]}
                    ]
                else:
                    positions = [
                        {"x": amino_acid_i["x"], "y": amino_acid_i["y"] + 1},
                        {"x": amino_acid_i["x"], "y": amino_acid_i["y"] - 1},
                        {"x": amino_acid_i_plus_3["x"], "y": amino_acid_i_plus_3["y"] + 1},
                        {"x": amino_acid_i_plus_3["x"], "y": amino_acid_i_plus_3["y"] - 1}
                    ]

                # Check if exactly two positions are free
                free_positions = [pos for pos in positions if pos not in [{"x": aa["x"], "y": aa["y"]} for aa in self.amino_acid_info_list]]

                if len(free_positions) == 2:
                    # Add crankshaft moves to the list with the "crankshaft_move" number
                    possible_moves.append({
                        f"crankshaft_move_{move_number}": move_number,
                        "id": self.amino_acid_info_list[i + 1]["id"],
                        "new_x": free_positions[0]["x"],
                        "new_y": free_positions[0]["y"]
                    })
                    possible_moves.append({
                        f"crankshaft_move_{move_number}": move_number,
                        "id": self.amino_acid_info_list[i + 2]["id"],
                        "new_x": free_positions[1]["x"],
                        "new_y": free_positions[1]["y"]
                    })

                    move_number += 1  # Increment move number

        return possible_moves

    def generate_all_possible_vshd_moves(self):
        """
        Generate all possible moves for the protein conformation, including end moves, corner moves, and crankshaft moves.

        Returns:
            list: List of dictionaries representing all possible moves.
        """
        all_possible_moves = []

        # Generate end moves
        end_moves = self.generate_end_moves()
        all_possible_moves.extend(end_moves)

        # Generate corner moves
        corner_moves = self.generate_corner_moves()
        all_possible_moves.extend(corner_moves)

        # Generate crankshaft moves
        crankshaft_moves = self.generate_crankshaft_moves()
        all_possible_moves.extend(crankshaft_moves)

        return all_possible_moves

    def generate_pull_move(self):
        """
        Generate a possible "pull move" for the protein conformation.

        Returns:
            list: A new conformation of amino acid info after applying the pull move.
        """
        # Create a deep copy of the current amino acid info list to modify
        new_conformation = copy.deepcopy(self.amino_acid_info_list)

        # Loop until a valid pull move is generated
        while True:
            # Select a random amino acid (i1) from the new conformation
            i1 = random.choice(new_conformation)

            # Determine the adjacent amino acid (i2) based on the selected amino acid (i1)
            if i1["id"] == 1:
                i2 = new_conformation[1]
            elif i1["id"] == len(new_conformation):
                i2 = new_conformation[-2]
            else:
                i2 = random.choice([new_conformation[i1["id"] - 2], new_conformation[i1["id"]]])

            # Calculate two possible positions (position_1 and position_2) for i1 based on i2's position
            if i1["x"] == i2["x"]:
                position_1 = {"x": i1["x"] + 1, "y": i1["y"]}
                position_2 = {"x": i1["x"] - 1, "y": i1["y"]}
            else:
                position_1 = {"x": i1["x"], "y": i1["y"] + 1}
                position_2 = {"x": i1["x"], "y": i1["y"] - 1}

            # Determine the adjacent amino acid (adjacent_amino_acid) for i1
            if i1["id"] == 1:
                adjacent_amino_acid = None
            elif i1["id"] == len(new_conformation):
                adjacent_amino_acid = None
            else:
                adjacent_amino_acid = new_conformation[i1["id"] - 2] if i1["id"] < i2["id"] else new_conformation[i1["id"]]

            x1, y1 = position_1["x"], position_1["y"]
            x2, y2 = position_2["x"], position_2["y"]

            # Check if positions 1 and 2 are occupied by other amino acids
            if any(aa["x"] == x1 and aa["y"] == y1 and aa != adjacent_amino_acid for aa in new_conformation) and \
                    any(aa["x"] == x2 and aa["y"] == y2 and aa != adjacent_amino_acid for aa in new_conformation):
                continue
            else:
                # Select one of the unoccupied positions as the move position (selected_position)
                selected_position = random.choice([position_1, position_2])
                if any(aa["x"] == selected_position["x"] and aa["y"] == selected_position["y"] and \
                       aa != adjacent_amino_acid for aa in new_conformation):
                    selected_position = position_1 if selected_position == position_2 else position_2

                # Determine the new move position (selected_move_position) based on i2's position
                if i1["x"] == i2["x"]:
                    position_3 = {"x": i2["x"] + 1, "y": i2["y"]}
                    position_4 = {"x": i2["x"] - 1, "y": i2["y"]}
                    selected_move_position = position_3 if selected_position["x"] == position_3["x"] else position_4
                else:
                    position_3 = {"x": i2["x"], "y": i2["y"] + 1}
                    position_4 = {"x": i2["x"], "y": i2["y"] - 1}
                    selected_move_position = position_3 if selected_position["y"] == position_3["y"] else position_4

                # Check if the selected move position is occupied by other amino acids
                if any(aa["x"] == selected_move_position["x"] and aa["y"] == selected_move_position["y"] for \
                       aa in new_conformation):
                    continue
                else:
                    # Update the position of i1 to the selected move position
                    new_conformation[i1["id"] - 1]["x"] = selected_move_position["x"]
                    new_conformation[i1["id"] - 1]["y"] = selected_move_position["y"]

                    # Determine the indices of i1 and i2 in the new conformation
                    i1_index = new_conformation[i1["id"] - 1]["id"]
                    i2_index = new_conformation[i2["id"] - 1]["id"]

                    # Check if i1 is at the beginning or end of the sequence
                    if i1_index == new_conformation[0]["id"] or i1_index == new_conformation[-1]["id"]:
                        return new_conformation
                    else:
                        # Adjust the positions of other amino acids in the chain
                        for i in range(i1_index - 2, -1, -1) if i1_index < i2_index else range(i1_index, len(new_conformation)):
                            if i1_index < i2_index:
                                j = i + 1
                            else:
                                j = i - 1
                            current_amino_acid = new_conformation[i]
                            moved_amino_acid = new_conformation[j]

                            # Calculate the Manhattan distance between the current and moved amino acids
                            distance = abs(current_amino_acid["x"] - moved_amino_acid["x"]) + abs(
                                current_amino_acid["y"] - moved_amino_acid["y"])

                            # Move the amino acid closer to the target position while maintaining distance
                            while distance > 1:
                                if current_amino_acid["x"] < moved_amino_acid["x"]:
                                    if any(aa["x"] == (current_amino_acid["x"] + 1) and aa["y"] ==
                                            current_amino_acid["y"] for aa in new_conformation):
                                        current_amino_acid["x"] = current_amino_acid["x"]
                                    else:
                                        if moved_amino_acid["y"] == (current_amino_acid["y"] + 2):
                                            if any(aa["x"] == (current_amino_acid["x"] + 1) and (aa["y"] ==
                                                   current_amino_acid["y"] + 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["x"] += 1
                                        elif moved_amino_acid["y"] == (current_amino_acid["y"] - 2):
                                            if any(aa["x"] == (current_amino_acid["x"] + 1) and (aa["y"] ==
                                                   current_amino_acid["y"] - 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["x"] += 1
                                        else:
                                            current_amino_acid["x"] += 1
                                if current_amino_acid["x"] > moved_amino_acid["x"]:
                                    if any(aa["x"] == (current_amino_acid["x"] - 1) and aa["y"] ==
                                            current_amino_acid["y"] for aa in new_conformation):
                                        current_amino_acid["x"] = current_amino_acid["x"]
                                    else:
                                        if moved_amino_acid["y"] == (current_amino_acid["y"] + 2):
                                            if any(aa["x"] == (current_amino_acid["x"] - 1) and (aa["y"] ==
                                                   current_amino_acid["y"] + 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["x"] -= 1
                                        elif moved_amino_acid["y"] == (current_amino_acid["y"] - 2):
                                            if any(aa["x"] == (current_amino_acid["x"] - 1) and (aa["y"] ==
                                                   current_amino_acid["y"] - 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["x"] -= 1
                                        else:
                                            current_amino_acid["x"] -= 1
                                if current_amino_acid["y"] < moved_amino_acid["y"]:
                                    if any(aa["x"] == current_amino_acid["x"] and aa["y"] == (
                                            current_amino_acid["y"] + 1) for aa in new_conformation):
                                        current_amino_acid["y"] = current_amino_acid["y"]
                                    else:
                                        if moved_amino_acid["x"] == (current_amino_acid["x"] + 2):
                                            if any(aa["x"] == (current_amino_acid["x"] + 1) and (aa["y"] ==
                                                   current_amino_acid["y"] + 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["y"] += 1
                                        elif moved_amino_acid["x"] == (current_amino_acid["x"] - 2):
                                            if any(aa["x"] == (current_amino_acid["x"] - 1) and (aa["y"] ==
                                                   current_amino_acid["y"] + 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["y"] += 1
                                        else:
                                            current_amino_acid["y"] += 1
                                if current_amino_acid["y"] > moved_amino_acid["y"]:
                                    if any(aa["x"] == current_amino_acid["x"] and aa["y"] == (
                                            current_amino_acid["y"] - 1) for aa in new_conformation):
                                        current_amino_acid["y"] = current_amino_acid["y"]
                                    else:
                                        if moved_amino_acid["x"] == (current_amino_acid["x"] + 2):
                                            if any(aa["x"] == (current_amino_acid["x"] + 1) and (aa["y"] ==
                                                   current_amino_acid["y"] - 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["y"] -= 1
                                        elif moved_amino_acid["x"] == (current_amino_acid["x"] - 2):
                                            if any(aa["x"] == (current_amino_acid["x"] - 1) and (aa["y"] ==
                                                   current_amino_acid["y"] - 1) for aa in new_conformation):
                                                current_amino_acid["x"] = current_amino_acid["x"]
                                            else:
                                                current_amino_acid["y"] -= 1
                                        else:
                                            current_amino_acid["y"] -= 1

                                # Update the position of the current amino acid
                                new_conformation[current_amino_acid["id"] - 1]["x"] = current_amino_acid["x"]
                                new_conformation[current_amino_acid["id"] - 1]["y"] = current_amino_acid["y"]

                                # Recalculate the distance for the next iteration
                                distance = abs(current_amino_acid["x"] - moved_amino_acid["x"]) + abs(
                                    current_amino_acid["y"] - moved_amino_acid["y"])

                        # Return the modified conformation after applying the pull move
                        return new_conformation