# Modified from the work of James Draper
# https://codereview.stackexchange.com/questions/162920/file-selection-button-for-jupyter-notebook
# This file is licensed under the Creative Commons Attribution-ShareAlike
# license (CC BY-SA 3.0)
# https://creativecommons.org/licenses/by-sa/3.0/

import traitlets
from ipywidgets import Button, Box, Label
from IPython.display import display
from tkinter import Tk, filedialog


class SaveButton(Button):
    save_function = None

    def __init__(self):
        super(SaveButton, self).__init__()

        # Add the selected_files trait
        self.add_traits(path=traitlets.traitlets.Unicode())

        # Create the button.
        self.description="Save graph"
        self.icon="save"
        # Set on click behavior.
        self.on_click(self.save_file)

    def select_file(self):
        """Generate instance of tkinter.filedialog.

        Parameters
        ----------
        button : obj:
            An instance of ipywidgets.widgets.Button
        """
        # Create Tk root
        root = Tk()
        # Hide the main window
        root.withdraw()
        # Raise the root to the top of all windows.
        root.call('wm', 'attributes', '.', '-topmost', True)
        # Path to use to save the file
        self.path = filedialog.asksaveasfilename()

    def save_file(self, button=None):
        if self.save_function is not None:
            self.select_file()
            self.save_function(self.path)
