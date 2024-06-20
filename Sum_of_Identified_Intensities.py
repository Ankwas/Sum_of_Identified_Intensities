# ----- Packages -----
from pyteomics import mass, tandem                                                      # Used to read xml file and calculate y- and b-ions
import customtkinter as ctk                                                             # Used for user interface
import tkinter as tk                                                                    # Used for user interface
from tkinter import Tk, ttk, Text, filedialog, messagebox, DISABLED, END, Scrollbar     # Used for user interface
import matplotlib.pyplot as plt                                                         # Used to make plots
from matplotlib.backends.backend_tkagg import NavigationToolbar2Tk, FigureCanvasTkAgg   # Used to make plots
from matplotlib.backends.backend_agg import FigureCanvasAgg                             # Used to make plots
import matplotlib as mpl                                                                # Used to make mouse hover label on plot
import mplcursors                                                                       # Used to make mouse hover label on plot
# ----- Functions -----
   
def peptide_info_extracter(xml_file_path, list_of_strings):
    """
    Funktion seachec an xml file for peptides containing one of multipel strings
    It returns info on the peptides containing thouse strings
    Args:
        xml_file_paht (string), String with the location of the file 
        List_of_strings (list), A list containing the strings the peptides must include
    Returs:
        list, A list containing the info about the peptides containing one of the strings
    """
    matching_spectra = []
    for string in list_of_strings:   
        with tandem.read(xml_file_path) as reader:
            for spectrum in reader:
                if string in spectrum['protein'][0]['peptide']['seq']:
                    id_ = spectrum['id']
                    max_i = spectrum['maxI']
                    expect = spectrum['expect']
                    start = spectrum['protein'][0]['peptide']['start']
                    hyperscore = spectrum['protein'][0]['peptide']['hyperscore']
                    y_ions = spectrum['protein'][0]['peptide']['y_ions']
                    b_ions = spectrum['protein'][0]['peptide']['b_ions']
                    seq = spectrum['protein'][0]['peptide']['seq']
                    try:
                        aa = spectrum['protein'][0]['peptide']['aa']
                    except KeyError:
                        aa = [{'type': 'None', 'at': 0, 'modified': 0}]
                    protein = spectrum['protein'][0]['peptide']['id']
                    xdata = spectrum['support']['fragment ion mass spectrum']['Xdata']['values']
                    ydata = spectrum['support']['fragment ion mass spectrum']['Ydata']['values']
                    charge = spectrum['support']['fragment ion mass spectrum']['charge']
                    spectrum_info = {
                        "ID": id_,
                        "Modification site": string,
                        "Max intensity": max_i,
                        "Expect": expect,
                        "Start" : start,
                        "Charge" : charge,
                        "Hyperscore": hyperscore,
                        "Number of y-ions": y_ions,
                        "Number of b-ions": b_ions,
                        "Sequence": seq,
                        "Modified aa": aa,
                        "m/z values": xdata,
                        "intensity:": ydata
                    }
                    matching_spectra.append(spectrum_info)
    return matching_spectra

def y_ion_calculator(peptide_sequence, start_position, modified_aa):
    """
    Takes a peptide sequen, the starting position of the peptide,
    a list of witch aa are modified and thire position and the charge of the spectrum
    and calculates a list of all the m/z values of y-ions in the peptide
    Args:
        peptide_sequence (string), the sequence of the peptide 
        start_position (int), the position of the first aa in the peptide in the protein 
        modified_aa (list), a list of modifications to the peptide, with position and weight 
    Returs:
        list, a list of all posibel mz values for y-ios and dubbelt chargeg y-ions
    """
    y_ion_mz = []
    for i in range(len(peptide_sequence)):
        y_ion_sequence = peptide_sequence[i:]
        y_ion_mass = mass.calculate_mass(sequence=y_ion_sequence, ion_type='y')
        for mod_aa in modified_aa:
            mod_position = mod_aa['at']
            mod_weight = mod_aa['modified']
            if mod_position >= start_position + i:
                y_ion_mass += mod_weight
        y_ion_mz.append(round(y_ion_mass+1.008,3))
        y_ion_mz.append(round(y_ion_mass+2.018,3)/2)
    return y_ion_mz

def b_ion_calculator(peptide_sequence, start_position, modified_aa):
    """
    Takes a peptide sequen, the starting position of the peptide,
    a list of witch aa are modified and thire position and the charge of the spectrum
    and calculates a list of all the m/z values of b-ions in the peptide
    Args:
        peptide_sequence (string), the sequence of the peptide 
        start_position (int), the position of the first aa in the peptide in the protein 
        modified_aa (list), a list of modifications to the peptide, with position and weight 
    Returs:
        list, a list of all posibelmz values for b-ios and dubbelt chargeg b-ions
    """
    b_ion_masses = []
    modification_rules = []
    for mod_aa in modified_aa:
        if 'at' in mod_aa: 
            start = mod_aa['at'] - start_position
            if start < 0:
                start_position += 1
            weight = mod_aa['modified']
            modification_rules.append({'start': start, 'weight': weight}) 
    for i in range(len(peptide_sequence)):
        b_ion_sequence = peptide_sequence[:i+1]
        b_ion_mass = mass.calculate_mass(sequence=b_ion_sequence, ion_type='b')
        for rule in modification_rules:
            if i+1 > rule['start']: 
                b_ion_mass += rule['weight']
        b_ion_masses.append(round(b_ion_mass+1.011,3))
        b_ion_masses.append(round(b_ion_mass+2.022,3)/2)
    return b_ion_masses

def calculate_soii(peptide_data, expect_value, min_hyperscore, y_ion_dict, b_ion_dict, precision):
    """
    Finds the y- and b-ions and calculates the sum of all the identifyid ions
    Args:
        peptide_data (dict), a dictionary containing all data on the peptides from the xml file 
        expect_value (float), the maximum expect value 
        min_hyperscore (float), the minimum hyperscor values 
        y_ion_dict (dict), A dictionary of all theoretical y-ions for each peptide
        b_ion_dict (dict), A dictionary of all theoretical b-ions for each peptide
        precision (float), the precision of the y- and b-ions
    Returs:
        dict, A dictionary containing the sum of all identifyed intensitys for each ID
    """
    soii_dict = {}
    for item in peptide_data:
        if item["Expect"] <= expect_value and item["Hyperscore"] >= min_hyperscore:
            y_ion_list = y_ion_dict.get(item['ID'], [])
            b_ion_list = b_ion_dict.get(item['ID'], [])
            sum_intensity = 0
            for mz_value, intensity in zip(item['m/z values'], item['intensity:']):
                for y_ion in y_ion_list:
                    if abs(mz_value - y_ion) <= precision:
                        sum_intensity += intensity
                        break
                else:
                    for b_ion in b_ion_list:
                        if abs(mz_value - b_ion) <= precision:
                            sum_intensity += intensity
                            break
            soii_dict[item['ID']] = sum_intensity
    return soii_dict

def summarize_soii(peptide_data, soii_dict, modification_mass):
    """
    A funktions that makes a summeraisation of the SoII based on:
    What modification sites they contain, if it is modified 
    Args:
        peptide_data (dict),
        soii_dict (dict),
        modification_mass (list),
    Returs:
        dict, a dictionary containing the sum of all identifyed intensetys sorted by:
                modification site and wether or not they contain the modification of intrest   
    """
    summary = {}
    for item in peptide_data:
        id_ = item['ID']
        if id_ in soii_dict:
            max_intensity = item['Max intensity']
            modification_site = item['Modification site']
            modified_aa = item['Modified aa']
            modification_mass_detected = False
            for aa in modified_aa:
                if 'modified' in aa and aa['modified'] == modification_mass:
                    modification_mass_detected = True
                    break
            soii_sum = soii_dict[id_] * max_intensity
            if modification_site not in summary:
                summary[modification_site] = {
                    'modification_mass_detected': [],
                    'sum_of_soii_multiplied_max_intensity': []
                }
            summary[modification_site]['modification_mass_detected'].append(modification_mass_detected)
            summary[modification_site]['sum_of_soii_multiplied_max_intensity'].append(soii_sum)
    sorted_summary = dict(sorted(summary.items()))
    for site in sorted_summary:
        modification_mass_detected_list = sorted_summary[site]['modification_mass_detected']
        sum_of_soii_multiplied_max_intensity_list = sorted_summary[site]['sum_of_soii_multiplied_max_intensity']
        sorted_indices = sorted(range(len(modification_mass_detected_list)), key=lambda k: modification_mass_detected_list[k], reverse=True)
        sorted_summary[site]['modification_mass_detected'] = [modification_mass_detected_list[i] for i in sorted_indices]
        sorted_summary[site]['sum_of_soii_multiplied_max_intensity'] = [sum_of_soii_multiplied_max_intensity_list[i] for i in sorted_indices]
    return sorted_summary

def check_file_input(file_path, query):
    """ 
    Checks that files are inputted correctly
    Args:
        file_path (string),
        query (string)
    Returs:
        bool 
    """
    if not file_path:
        messagebox.showerror("Error", f"Please select a {query}.")
        return False  
    return True

def check_value_input(value, query):
    """ 
    Checks that values are inputted correctly and are difrent from 0
    Args:
        value (string),
        query (string),
    Returs:
        value (float),
    """
    error_occurred = False 
    if not value:
        messagebox.showerror("Error", f"Please enter a value for {query}.")
        error_occurred = True 
    try:
        value = float(value)
        if value == 0:
            raise ValueError
    except ValueError:
        messagebox.showerror("Error", f"{query} must be a valid non-zero number.")
        error_occurred = True
    if error_occurred:
        return None
    return float(value)
    
def check_inputs(xml_file_path, substrings, min_hyperscore, modification_mass, expect_value, precision):
    """
    Gets all the inputs from the buttons
    Checks that all inputs are filled out and valid.
    Then calls the main function with all the inputs     
    """
    xml_file_valid = check_file_input(xml_file_path, "xml file")
    min_hyperscore_valid = check_value_input(min_hyperscore, "minimum hyperscore")
    modification_mass_valid = check_value_input(modification_mass, "modification mass")
    expect_value_valid = check_value_input(expect_value, "maximum expect value")/1000
    precision_valid = check_value_input(precision, "precision for y- and b- ions ")
    if not all([xml_file_valid, min_hyperscore_valid, modification_mass_valid, expect_value_valid, precision_valid]):
        return  
    main(xml_file_path, substrings, min_hyperscore_valid, modification_mass_valid, expect_value_valid, precision_valid)


def main(xml_file_path, list_of_substrings, min_hyperscore, modification_mass, expect_value, precision):
    """
    Main function to process the XML file extract relevant data and calculate SoII.
    Then feads the relevant inputs into the UI for to see the results 
    Args:
        xml_file_path (str): The path to the XML file.
        substrings (list): List of substrings to search for.
        min_hyperscore (float): The minimum hyperscore value.
        modification_mass (float): The value of modification mass.
    """
    peptide_data = peptide_info_extracter(xml_file_path, list_of_substrings)
    y_ion_dict = {}
    b_ion_dict = {}
    for item in peptide_data:
        y_ion_mz = y_ion_calculator(item['Sequence'], item['Start'], item['Modified aa'])
        b_ion_mz = b_ion_calculator(item['Sequence'], item['Start'], item['Modified aa'])
        y_ion_list = []
        b_ion_list = []
        for value in item['m/z values']:
            for y_ion in y_ion_mz:
                if abs(value - y_ion) <= precision:
                    y_ion_list.append(value)
                    break
            for b_ion in b_ion_mz:
                if abs(value - b_ion) <= precision:
                    b_ion_list.append(value)
                    break
        y_ion_dict[item['ID']] = y_ion_list
        b_ion_dict[item['ID']] = b_ion_list
    soii_dict = calculate_soii(peptide_data, expect_value, min_hyperscore, y_ion_dict, b_ion_dict, precision)
    soii_summary = summarize_soii(peptide_data, soii_dict, modification_mass)
    results_ui(soii_summary,soii_dict, peptide_data, y_ion_dict, b_ion_dict, precision, expect_value, min_hyperscore, modification_mass)

# ----- User interface -----
# ----- Resutls UI -----
def display_summary_tab(treeview, summary_data):
    """
    Display summary of soii data in the tab 
    """
    treeview.delete(*treeview.get_children())
    for col in treeview["columns"]:
        treeview.heading(col, anchor='center')
        treeview.column(col, anchor='center')
    for modification_site, data in summary_data.items():
        total_soii_site = sum(data['sum_of_soii_multiplied_max_intensity'])
        treeview.insert('', 'end', values=[modification_site, "", "", ""]) 
        if total_soii_site != 0:
            t_soii = sum(data['sum_of_soii_multiplied_max_intensity'][i] for i in range(len(data['modification_mass_detected'])) if data['modification_mass_detected'][i])
            t_percentage = (t_soii / total_soii_site) * 100
            f_percentage = ((total_soii_site - t_soii) / total_soii_site) * 100
        else:
            t_percentage = f_percentage = 0
        treeview.insert('', 'end', values=["", "T", t_soii, f"{t_percentage:.0f}%"])
        treeview.insert('', 'end', values=["", "F", total_soii_site - t_soii, f"{f_percentage:.0f}%"])

def plot_spectrum(selected_id, peptide_data, y_ion_dict, b_ion_dict, precision, plot_frame):
    """
    Plot spectra in tab 2.   
    """
    data = next(item for item in peptide_data if item["ID"] == selected_id)
    mz_values = data["m/z values"]
    intensities = data["intensity:"]
    for widget in plot_frame.winfo_children():
        widget.destroy()
    fig, ax = plt.subplots(figsize=(16,8))
    for mz, intensity in zip(mz_values, intensities):
        is_y_ion = any(abs(mz - y_ion) < precision for y_ion in y_ion_dict.get(selected_id, []))
        is_b_ion = any(abs(mz - b_ion) < precision for b_ion in b_ion_dict.get(selected_id, []))
        if is_y_ion and is_b_ion:
            ax.stem([mz], [intensity], linefmt='purple', markerfmt='None', basefmt='k-')
        elif is_y_ion:
            ax.stem([mz], [intensity], linefmt='r-', markerfmt='None', basefmt='k-')
        elif is_b_ion:
            ax.stem([mz], [intensity], linefmt='b-', markerfmt='None', basefmt='k-')
        else:
            ax.stem([mz], [intensity], linefmt='k-', markerfmt='None', basefmt='k-')
    
    ax.set_xlabel("m/z values")
    ax.set_ylabel("Intensity in % of max intensety")
    ax.set_title(f"Spectrum for ID: {selected_id}")
    ax.grid(True)

    cursor = mplcursors.cursor(hover=True)
    @cursor.connect("add")
    def on_hover(sel):
        if isinstance(sel.artist, mpl.container.StemContainer):
            x = sel.target[0]
            y = sel.target[1]
        elif isinstance(sel.artist, mpl.lines.Line2D):
            index = sel.index
            x = sel.artist.get_xdata()[index]
            y = sel.artist.get_ydata()[index]
        else:
            x, y, *_ = sel.artist.get_offsets()[sel.index]
        sel.annotation.set_text(f"x={x:.2f}, y={y:.2f}")

    canvas = FigureCanvasTkAgg(fig, master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    toolbar_frame = tk.Frame(plot_frame)
    toolbar = NavigationToolbar2Tk(canvas, toolbar_frame)
    toolbar.update()
    toolbar_frame.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

def display_plot_tab(peptide_data, y_ion_dict, b_ion_dict, precision, plot_frame, id_listbox, soii_dict):
    """
    Display plot tab.
    """
    for widget in plot_frame.winfo_children():
        widget.destroy()
    
    plot_container = tk.Frame(plot_frame)
    plot_container.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

    larger_font = ('Arial', 14)

    text_row_plot = tk.Label(plot_container, text="ID's with a ✓ are included in the SoII", justify="center", font=larger_font)
    text_row_plot.grid(column=0, row=0)
    
    text_row_plot = tk.Label(plot_container, text="Blue = b-ions", justify="center", fg='blue', font=larger_font)
    text_row_plot.grid(column=1, row=0)

    text_row_plot = tk.Label(plot_container, text="Red = y-ions", justify="center", fg='red', font=larger_font)
    text_row_plot.grid(column=2, row=0)

    text_row_plot = tk.Label(plot_container, text="Black = not identified", justify="center", font=larger_font)
    text_row_plot.grid(column=3, row=0)

    plot_canvas = tk.Canvas(plot_container)
    plot_canvas.grid(column=0, row=1, columnspan=4, sticky="nsew")

    def on_select(event):
        index = event.widget.curselection()[0]
        selected_item = event.widget.get(index)
        selected_id = selected_item.split(":")[1].split(" ")[1].strip()
        plot_spectrum(selected_id, peptide_data, y_ion_dict, b_ion_dict, precision, plot_canvas)

    id_listbox.bind("<<ListboxSelect>>", on_select)
    id_listbox.delete(0, tk.END)
    for item in peptide_data:
        id_ = item["ID"]
        if id_ in soii_dict:
            id_ += " [✓]"
        else:
            id_ += " [⠀]"
        id_listbox.insert(tk.END, f"ID: {id_}")


def display_data_tab(tab3, peptide_data, soii_dict, expect_value, min_hyperscore, modification_mass):
    """
    Display data in tab.
    """
    for widget in tab3.winfo_children():
        widget.destroy()
    text_above_left_tab_3 = ttk.Label(tab3, text="")
    text_above_left_tab_3.pack(side='top')
    text_left = tk.Text(tab3, wrap='word')
    text_left.pack(fill='both', expand=True)
    text_left_above = f"Expect value cut of: {expect_value*1000} e-3    Minimum hyperscore: {min_hyperscore}"
    text_above_left_tab_3.config(text=text_left_above)
    headers = ["ID", "Modification site", "Modified", "Expect value", "Hyperscore", "Included in SoII"]
    widths = [8, 20, 10, 15, 10, 15]
    header_format = "{:<{w[0]}} {:<{w[1]}} {:<{w[2]}} {:<{w[3]}} {:<{w[4]}} {:<{w[5]}}\n".format(*headers, w=widths)
    text_left_data = header_format
    for item in peptide_data:
        id_ = item["ID"]
        expect = item["Expect"]
        hyperscore = item["Hyperscore"]
        modification_site = item["Modification site"]
        modified_aa = item.get("Modified aa", [])
        contains_modification_mass = any(aa.get("modified") == modification_mass for aa in modified_aa)
        contains_modification_mass = "True" if contains_modification_mass else "False"
        included_in_soii = "✓" if id_ in soii_dict else "⠀"
        row = "{:<{w[0]}} {:<{w[1]}} {:<{w[2]}} {:<{w[3]}} {:<{w[4]}} {:<{w[5]}}\n".format(id_, modification_site, contains_modification_mass, expect, hyperscore, included_in_soii, w=widths)
        text_left_data += row
    text_left.insert('end', text_left_data)

def results_ui(soii_summary, soii_dict, peptide_data, y_ion_dict, b_ion_dict, precision, expect_value, min_hyperscore, modification_mass):
    """
    Creates a user interface with a tabed window, each tab containing difrent information.
    The user can select from three tabs in the top of the window:
        the "summary tab", the "plot tab" and the "data tab"
    The user is placed on the "summery tab" as default when this function is called
    Args:

    """
    root = tk.Tk()
    root.geometry("1000x600") 
    root.title("Results of SoII")

    notebook = ttk.Notebook(root)
    notebook.pack(fill='both', expand=True)

    tab1 = ttk.Frame(notebook)
    tab2 = ttk.Frame(notebook)
    tab3 = ttk.Frame(notebook)
    notebook.add(tab1, text="Sum of Identified Intensities")
    notebook.add(tab2, text="Plots of intensity")
    notebook.add(tab3, text="Data for SoII before and after filtering")

    text_above_left = ttk.Label(tab1, text="")
    text_above_left.pack(side='top')
    treeview = ttk.Treeview(tab1, columns=['Modification site', 'Modified', 'SoII', 'Percentage'], show='headings', style='Custom.Treeview')
    treeview.heading('Modification site', text='Modification site', anchor='center')
    treeview.heading('Modified', text='Modified', anchor='center')
    treeview.heading('SoII', text='SoII', anchor='center')
    treeview.heading('Percentage', text='Percentage', anchor='center')
    treeview.pack(expand=True, fill='both')

    id_listbox = tk.Listbox(tab2)
    id_listbox.pack(side=tk.LEFT, fill=tk.Y)
    plot_frame = tk.Frame(tab2)
    plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    text_above_left_tab_3 = ttk.Label(tab3, text="")
    text_above_left_tab_3.pack(side='top')
    text_left = tk.Text(tab3, wrap='word')
    text_left.pack(fill='both', expand=True)

    def on_tab_changed(event):
        """
        Calls the function to display the selected tab on the event that a tab is selected 
        """
        tab_index = notebook.index(notebook.select())
        if tab_index == 0:
            display_summary_tab(treeview, soii_summary)
        elif tab_index == 1:
            display_plot_tab(peptide_data, y_ion_dict, b_ion_dict, precision, plot_frame, id_listbox, soii_dict)
        elif tab_index == 2:
            display_data_tab(tab3, peptide_data, soii_dict, expect_value, min_hyperscore, modification_mass)

    notebook.bind("<<NotebookTabChanged>>", on_tab_changed)

    root.mainloop()  
# ----- Setup UI -----
def indput_ui():
    """
    Function to create the user interface when the program starts.
    This user interface contains input fields to alow the user to input values
    And buttons to set inputs to predetermint values
    The "Calculate Sum of Identified Intensities" button starts the program.
    """
    root = ctk.CTk()
    root.title("Sum of Identified Intensities calculator")

    # Functions for  each buttons 
    def browse_input_xml_file():
        """
        Makes a button that opens a file explorer that only shows xml files to select as input
        """
        filename = filedialog.askopenfilename(filetypes=[("XML files", "*.xml")])
        entry_input_xml_file.delete(0, ctk.END)
        entry_input_xml_file.insert(0, filename)

    def set_target_difference():
        """
        Makes a button that sets the input to kallistatins glycosylation sites
        """
        entry_target_difference.delete(0, ctk.END)
        entry_target_difference.insert(0, "NSS,NLT,NDT,NTT")

    def set_O18_mass():
        """
        Makes a button that sets the input to the mass of a 18O deglycosylation
        """
        entry_modification_mass.delete(0, ctk.END)
        entry_modification_mass.insert(0, "2.98900")

    def run_program():
        """
        Gets all the inuts made in input_ui and gives them to check_inputs
        """
        check_inputs(entry_input_xml_file.get(),
                     entry_target_difference.get().split(","),
                     entry_min_hyperscore.get(),
                     entry_modification_mass.get(),
                     entry_expect_value.get(),
                     entry_precision.get())
    
    # Create input fields and buttons
    label_input_xml_file = ctk.CTkLabel(root, text="Input XML File:", justify='right')
    label_input_xml_file.grid(row=2, column=0)
    entry_input_xml_file = ctk.CTkEntry(root)
    entry_input_xml_file.grid(row=2, column=1)
    button_browse_input_xml = ctk.CTkButton(root, text="Browse", command=browse_input_xml_file)
    button_browse_input_xml.grid(row=2, column=2, sticky="ew")

    label_target_difference = ctk.CTkLabel(root, text="   Amino acids sequences of interest (use ',' to separate):")
    label_target_difference.grid(row=3, column=0)
    entry_target_difference = ctk.CTkEntry(root)
    entry_target_difference.grid(row=3, column=1)
    button_set_target_difference = ctk.CTkButton(root, text="Set to kallistatins glycosylation sites:", command=set_target_difference)
    button_set_target_difference.grid(row=3, column=2, sticky="ew")

    label_modification_mass = ctk.CTkLabel(root, text="Mass of modification of interest:")
    label_modification_mass.grid(row=4, column=0)
    entry_modification_mass = ctk.CTkEntry(root)
    entry_modification_mass.grid(row=4, column=1)
    button_set_O18_mass = ctk.CTkButton(root, text="Set to \u00B9\u2078O modification mass", command=set_O18_mass)
    button_set_O18_mass.grid(row=4, column=2, sticky="ew")

    label_min_hyperscore = ctk.CTkLabel(root, text="Minimum hyper score (ex 20):")
    label_min_hyperscore.grid(row=5, column=0)
    entry_min_hyperscore = ctk.CTkEntry(root)
    entry_min_hyperscore.grid(row=5, column=1)

    label_expect_value = ctk.CTkLabel(root, text="Maximum expect score in *10\u207B \u00B3 (ex 1 for 1*10\u207B \u00B3 ):")
    label_expect_value.grid(row=6, column=0)
    entry_expect_value = ctk.CTkEntry(root)
    entry_expect_value.grid(row=6, column=1)

    label_precision = ctk.CTkLabel(root, text="Precision for y- and b- ions (ex 0.1 ):")
    label_precision.grid(row=7, column=0)
    entry_precision = ctk.CTkEntry(root)
    entry_precision.grid(row=7, column=1)

    button_run = ctk.CTkButton(root, text="Calculate Sum of Identified Intensities", command = run_program)
    button_run.grid(row=8, column=1)

    root.mainloop()

indput_ui()
