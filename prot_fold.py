# Importation of python modules
import tkinter as tk
import src.gui_prot_fold as gui

# Check if this script is the main module being executed
if __name__ == "__main__":
    # Create a Tkinter root window
    root = tk.Tk()

    # Create an instance of the ProteinFolderApp GUI application
    app = gui.ProteinFolderApp(root)

    # Run the GUI application
    app.run()