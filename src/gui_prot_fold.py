# Importation of python modules
import tkinter as tk
from tkinter import filedialog
from tkinter import ttk
import matplotlib.backends.backend_tkagg as tkagg
import matplotlib.pyplot as plt
import time
import src.peptide as peptide
import src.simulation as simulation

class ProteinFolderApp:
    """
    A class representing the Protein Folder application.

    Attributes:
        root (tk.Tk): The root tkinter window.
        temperature (float): The simulation temperature.
        nb_iterations (int): The number of simulation iterations.
        move_mode (str): The selected move mode.
        energy_history (list): A list to store energy history during simulation.

    Methods:
        __init__(self, root): Initializes the ProteinFolderApp instance.
        initialize_gui(self): Initializes the GUI elements of the application.
    """

    def __init__(self, root):
        """
        Initializes the ProteinFolderApp instance.

        Args:
            root (tk.Tk): The root tkinter window.
        """
        self.root = root
        root.title("Prot Fold")

        # Initialize default values
        self.temperature = 1.0
        self.nb_iterations = 10
        self.move_mode = "Pull move"
        self.energy_history = []

        # Initialize GUI elements
        self.initialize_gui()

    def initialize_gui(self):
        """
        Initializes the GUI elements of the application.
        """
        # Create canvas
        self.canvas_width = 600
        self.canvas_height = 600
        self.canvas = tk.Canvas(self.root, width=self.canvas_width, height=self.canvas_height)
        self.canvas.pack()

        # Draw a gray square
        square_points = [23, 23, 482, 23, 482, 577, 23, 577]
        gray_square = self.draw_square(square_points, "gray")

        # Create a "Browse" button
        browse_button = tk.Button(self.root, text="Browse", command=self.browse_file)
        browse_button.place(x=495, y=25)  # Adjust the position as needed

        # Create a Temperature label and input box
        temperature_label = tk.Label(self.root, text="Temperature:")
        temperature_label.place(x=495, y=65)
        self.temperature_entry = tk.Entry(self.root, width=16)
        self.temperature_entry.place(x=495, y=95)

        # Create an Nb of iterations label and input box
        iterations_label = tk.Label(self.root, text="Nb of MC steps:")
        iterations_label.place(x=495, y=135)
        self.iterations_entry = tk.Entry(self.root, width=16)
        self.iterations_entry.place(x=495, y=165)

        # Create a Move mode dropdown list
        move_mode_label = tk.Label(self.root, text="Move mode:")
        move_mode_label.place(x=495, y=205)
        self.move_mode_var = tk.StringVar()
        move_mode_combobox = ttk.Combobox(self.root, textvariable=self.move_mode_var, values=["Pull move", "VSHD"], width=13, state="readonly")
        move_mode_combobox.set("Pull move")  # Set default choice
        move_mode_combobox.place(x=495, y=235)

        # Create a "Start" button
        start_button = tk.Button(self.root, text="Start", command=self.start_simulation)
        start_button.place(x=495, y=275)  # Adjust the position as needed

        # Label for "Stats:"
        stats_label = tk.Label(self.root, text="Stats:")
        stats_label.place(x=495, y=355)

        # Label for "Current energy:"
        current_energy_label = tk.Label(self.root, text="Current energy:")
        current_energy_label.place(x=495, y=395)

        # Create a text variable for displaying the current energy
        self.current_energy_text = tk.StringVar()

        # Text widget to display current energy
        current_energy_display = tk.Label(self.root, textvariable=self.current_energy_text)
        current_energy_display.place(x=495, y=425)

        # Create an "Energy history" button (initially disabled)
        self.energy_history_button = tk.Button(self.root, text="Energy history", command=self.plot_energy_history, state="disabled")
        self.energy_history_button.place(x=495, y=465)  # Adjust the position as needed

        # Label for "Calculation time:"
        calculation_time_label = tk.Label(self.root, text="Calculation time:")
        calculation_time_label.place(x=495, y=505)

        # Create a text variable for displaying the calculation time
        self.calculation_time_text = tk.StringVar()

        # Text widget to display the calculation time
        calculation_time_display = tk.Label(self.root, textvariable=self.calculation_time_text)
        calculation_time_display.place(x=495, y=535)

    def draw_circle(self, x, y, color, radius):
        """
        Draw a filled circle on the canvas.

        Args:
            x (int): The x-coordinate of the circle's center.
            y (int): The y-coordinate of the circle's center.
            color (str): The fill color of the circle.
            radius (int): The radius of the circle.

        Returns:
            int: The ID of the drawn circle object on the canvas.
        """
        return self.canvas.create_oval(x - radius, y - radius, x + radius, y + radius, fill=color)

    def draw_square(self, points, color):
        """
        Draw a filled square on the canvas.

        Args:
            points (list): A list of coordinates defining the square's corners.
            color (str): The fill color of the square.

        Returns:
            int: The ID of the drawn square object on the canvas.
        """
        return self.canvas.create_polygon(points, fill=color)

    def draw_pink_line(self, x1, y1, x2, y2):
        """
        Draw a pink line on the canvas.

        Args:
            x1 (int): The x-coordinate of the starting point of the line.
            y1 (int): The y-coordinate of the starting point of the line.
            x2 (int): The x-coordinate of the ending point of the line.
            y2 (int): The y-coordinate of the ending point of the line.

        Returns:
            int: The ID of the drawn line object on the canvas.
        """
        return self.canvas.create_line(x1, y1, x2, y2, fill='pink', width=2)

    def draw_black_line(self, x1, y1, x2, y2):
        """
        Draw a black line on the canvas.

        Args:
            x1 (int): The x-coordinate of the starting point of the line.
            y1 (int): The y-coordinate of the starting point of the line.
            x2 (int): The x-coordinate of the ending point of the line.
            y2 (int): The y-coordinate of the ending point of the line.

        Returns:
            int: The ID of the drawn line object on the canvas.
        """
        return self.canvas.create_line(x1, y1, x2, y2, fill='black', width=2)

    def display_conformation(self, conformation_data):
        """
        Display the conformation of amino acids on the canvas.

        Args:
            conformation_data (list): A list of dictionaries representing amino acid conformation data.

        This method clears the canvas and then draws circles representing amino acids at their specified positions.
        'H' amino acids are drawn in red, while others are drawn in blue. Pink lines connect adjacent 'H' amino acids,
        and black lines connect all adjacent amino acids.

        """
        circle_radius = 10
        self.canvas.delete("all")  # Clear the canvas
        num_amino_acids = len(conformation_data)

        # Draw a gray square
        square_points = [23, 23, 482, 23, 482, 577, 23, 577]
        gray_square = self.draw_square(square_points, "gray")

        if num_amino_acids > 1:
            # Calculate the circle spacing based on the provided constraints
            circle_spacing_x = (self.canvas_width - 180 - 2 * circle_radius) / (
                    max(conformation_data, key=lambda x: x['x'])['x'] -
                    min(conformation_data, key=lambda x: x['x'])['x'])
            circle_spacing_y = (self.canvas_height - 85 - 2 * circle_radius) / (
                    max(conformation_data, key=lambda x: x['y'])['y'] -
                    min(conformation_data, key=lambda x: x['y'])['y'])

            # Use the smaller of the two spacings to ensure the peptide fits within the window
            circle_spacing = min(circle_spacing_x, circle_spacing_y)

            # Calculate the circle radius as 1/4 of the spacing
            circle_radius = circle_spacing / 4

            # Calculate the offsets to position the peptide within the window
            x_offset = 25 + circle_radius
            y_offset = 25 + circle_radius

            # Create a set to store the positions of 'H' amino acids
            h_positions = set()

            for i in range(num_amino_acids):
                data = conformation_data[i]
                x = x_offset + (data['x'] - min(conformation_data, key=lambda x: x['x'])['x']) * circle_spacing
                y = y_offset + (max(conformation_data, key=lambda x: x['y'])['y'] - data['y']) * circle_spacing

                if data['polarity'] == 'H':
                    h_positions.add((x, y))

                if data['polarity'] == 'H':
                    color = 'red'
                else:
                    color = 'blue'

                circle = self.draw_circle(x, y, color, circle_radius)

            # Draw pink lines between 'H' amino acids that are neighbors
            for i in range(num_amino_acids):
                data = conformation_data[i]
                x = x_offset + (data['x'] - min(conformation_data, key=lambda x: x['x'])['x']) * circle_spacing
                y = y_offset + (max(conformation_data, key=lambda x: x['y'])['y'] - data['y']) * circle_spacing

                if data['polarity'] == 'H':
                    for j in range(i + 1, num_amino_acids):
                        neighbor_data = conformation_data[j]
                        neighbor_x = x_offset + (neighbor_data['x'] - min(conformation_data, key=lambda x: x['x'])[
                            'x']) * circle_spacing
                        neighbor_y = y_offset + (max(conformation_data, key=lambda x: x['y'])['y'] - neighbor_data[
                            'y']) * circle_spacing

                        # Check if the neighbor is an 'H' amino acid and is a neighbor in positions (x+1,y), (x-1,y), (x,y+1), (x,y-1)
                        if neighbor_data['polarity'] == 'H' and (
                                (neighbor_x == x + circle_spacing and neighbor_y == y) or
                                (neighbor_x == x - circle_spacing and neighbor_y == y) or
                                (neighbor_x == x and neighbor_y == y + circle_spacing) or
                                (neighbor_x == x and neighbor_y == y - circle_spacing)):
                            # Draw a pink line between two adjacent 'H' amino acids
                            self.draw_pink_line(x, y, neighbor_x, neighbor_y)

            # Draw the black lines over the circles
            for i in range(1, num_amino_acids):
                data = conformation_data[i]
                prev_data = conformation_data[i - 1]
                x = x_offset + (data['x'] - min(conformation_data, key=lambda x: x['x'])['x']) * circle_spacing
                y = y_offset + (max(conformation_data, key=lambda x: x['y'])['y'] - data['y']) * circle_spacing
                x_prev = x_offset + (
                        prev_data['x'] - min(conformation_data, key=lambda x: x['x'])['x']) * circle_spacing
                y_prev = y_offset + (
                        max(conformation_data, key=lambda x: x['y'])['y'] - prev_data['y']) * circle_spacing
                # Draw a black line between all adjacent amino acids
                self.draw_black_line(x_prev, y_prev, x, y)

    def calculate_and_display_energy(self, conformation, text_variable):
        """
        Calculate and display the energy of the given conformation.

        Args:
            conformation (list): A list of dictionaries representing amino acid conformation data.
            text_variable (tk.StringVar): The text variable for displaying the energy.

        This method calculates the energy of the given conformation using an external function and sets the calculated energy
        as text in the specified text variable.

        """
        sequence = peptide.AminoAcidSequence()
        energy = sequence.calculate_energy(conformation)
        text_variable.set(f"{energy:.2f}")

    def browse_file(self):
        """
        Open a file dialog to browse and load a FASTA file for conformation analysis.

        This method opens a file dialog to allow the user to select a FASTA file. If a valid file is selected, it reads the
        amino acid sequence from the file, processes it, and updates the initial conformation, current energy, and displays
        the conformation on the canvas.

        """
        # Open a file dialog to select a FASTA file with a .fasta extension
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", "*.fasta")])

        # Check if a valid file was selected
        if file_path:
            # Create an instance of AminoAcidSequence to work with sequences
            sequence = peptide.AminoAcidSequence()

            # Read the contents of the selected FASTA file into a list
            amino_acids_list = sequence.read_fasta_to_list(file_path)

            # Replace any amino acids in the list if necessary (not shown in the code, but assumed)
            modified_sequence = sequence.replace_amino_acids(amino_acids_list)

            # Create an amino acid information list from the modified sequence
            self.initial_conformation = sequence.create_amino_acid_info_list(modified_sequence)

            # Calculate and display the energy of the initial conformation
            self.calculate_and_display_energy(self.initial_conformation, self.current_energy_text)

            # Display the initial conformation on the canvas
            self.display_conformation(self.initial_conformation)

    def start_simulation(self):
        """
        Start the simulation based on user-provided parameters.

        This method initiates a simulation based on the user-provided temperature, number of iterations, and move mode.
        It records the simulation's execution time and updates the current energy, calculation time, and displays the
        final conformation on the canvas.

        """
        # Retrieve user-provided temperature, number of iterations, and move mode
        temperature = float(self.temperature_entry.get())
        nb_iterations = int(self.iterations_entry.get())
        move_mode = self.move_mode_var.get()

        # Record the start time of the simulation
        start_time = time.time()

        # Run the Monte Carlo simulation using the provided parameters
        final_conformation, self.energy_history = simulation.MCSimulation(self.initial_conformation).mc_simulation( \
            temperature, move_mode, nb_iterations)

        # Record the end time of the simulation
        end_time = time.time()

        # Calculate the elapsed time for the simulation
        elapsed_time = end_time - start_time

        # Update the calculation time text to display the elapsed time in seconds
        self.calculation_time_text.set(f"{elapsed_time:.2f} seconds")

        # Display the final conformation on the canvas
        self.display_conformation(final_conformation)

        # Create an instance of AminoAcidSequence to calculate the current energy
        sequence = peptide.AminoAcidSequence()

        # Calculate the energy of the final conformation
        current_energy = sequence.calculate_energy(final_conformation)

        # Update the current energy text to display the calculated energy
        self.current_energy_text.set(f"{current_energy}")

        # Enable the energy history button to allow users to view energy history
        self.energy_history_button.config(state="normal")

    def plot_energy_history(self):
        """
        Plot the energy history of the simulation.

        This method creates a plot of the energy history using Matplotlib, displaying the energy changes over iterations.
        It opens a new Tkinter window to show the plot.

        """
        if hasattr(self, 'energy_history') and self.energy_history:
            iterations = list(self.energy_history.keys())
            energies = list(self.energy_history.values())

            # Create the energy history plot
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.plot(iterations, energies, marker='o', linestyle='-', color='b')
            ax.set_title("Energy History")
            ax.set_xlabel("Iterations")
            ax.set_ylabel("Energy")
            ax.grid(True)

            # Create a new Tkinter window for the plot
            plot_window = tk.Toplevel(self.root)

            # Embed the Matplotlib plot in a Tkinter canvas within the new window
            canvas_matplotlib = tkagg.FigureCanvasTkAgg(fig, master=plot_window)
            canvas_matplotlib.get_tk_widget().pack()

            # Show the plot
            canvas_matplotlib.draw()

    def run(self):
        """
        Start the Tkinter main loop.

        This method starts the Tkinter main loop, which is required to run the graphical user interface (GUI).

        """
        self.root.mainloop()

if __name__ == "__main__":
    root = tk.Tk()
    app = ProteinFolderApp(root)
    app.run()
