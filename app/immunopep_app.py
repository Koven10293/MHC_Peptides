import pandas as pd
import requests
import numpy as np
import tkinter as tk
from tkinter import filedialog as fd
from tkinter import ttk

from matplotlib import pyplot as plt
from matplotlib.figure import Figure 
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg 
import seaborn as sns

import time

class Window:
    
    def __init__(self, master):
        self.master = master
        master.title("IPD")
        master.geometry("500x500")
        self.peptide_csv = None
        self.csv_status = "No CSV selected"
        self.input_frame = None
        self.peptide_frame = None
        self.protein_frame = None
        self.filepath = None
        self.filename = None
        
        #Canvas for plot
        self.canvas = None
        
        self.frm_top = tk.Frame(self.master)
        self.frm_2 = tk.Frame(self.master)
        self.frm_bot = tk.Frame(self.master)
        
        #Progress bar for processing progress
        self.process_progress = ttk.Progressbar(master, orient = "horizontal", length = 100, mode = "determinate")
        
        #Button that allows user to upload a csv
        self.btn_upload = tk.Button(self.frm_top, text = 'Select CSV', width = 25, command = self.get_csv)
        #Button that triggers processing of csv into peptide and protein frames
        self.btn_process = tk.Button(self.frm_top, text = 'Process CSV', command = self.start_process)
        
        #input field for protein charts
        self.protein_selection = tk.Entry(self.master)
        
        #VIS BUTTONS
        #Button that visualises heat map of 
        self.btn_heat = tk.Button(self.frm_bot, text = "Heat Map", command = self.plot_fig_locations_heat)
        #Button that displays protein distribution
        self.btn_prot = tk.Button(self.frm_bot, text = "Protein", command = self.plot_fig_prot_dist_heat)
        #Button that displays length distribution
        self.btn_len = tk.Button(self.frm_bot, text = "Length Distribution", command = self.plot_fig_len_dist)
        
        #CSV status label 
        self.lbl_is_csv = tk.Label(self.frm_top, text = self.csv_status, width = 25)
        
        self.frm_top.pack()
        self.frm_2.pack()
        self.frm_bot.pack(side = "bottom")
        self.btn_upload.pack(side = 'left')
        self.btn_heat.pack()
        self.btn_prot.pack(side = "left")
        self.btn_len.pack(side = "right")
        self.lbl_is_csv.pack()
        self.protein_selection.pack()
        
    def start_process(self):
        """
        Input:
        -self: Class instance
        
        Output: None
        
        Description: Called by button to start processing selected CSV
        """
        self.process_frame(self.input_frame)
    
    def clean(self, a_frame):
        """
        Input: 
        -self: Class instance
        a_frame: A dataframe generated from the .csv output of *mass spec* 
        
        Output: Dataframe with rows containing NA values for Accession removed and charges removed from peptides
        
        Description: Cleans dataframe for processing
        """
        #remove all #""# accessions
        clean_frame = a_frame.copy()
        clean_frame["Accession"] = clean_frame["Accession"].str.replace(r'\#.+\#[^\:]+\:', '')
        clean_frame["Accession"] = clean_frame["Accession"].str.replace(r'\:\#.+\#[^\:]+$', '')
        clean_frame["Accession"] = clean_frame["Accession"].str.replace(r'\#.+\#[^\:]+$', 'NaN')
        clean_frame = clean_frame.replace("NaN", np.nan)
        #Remove charges from peptides
        clean_frame["Peptide"] = clean_frame["Peptide"].str.replace(r'[\(\+\d+\.\d+\)]', '')
        #remove NA accession entries
        clean_frame = clean_frame.dropna(subset = ["Accession"]).reset_index()
        return clean_frame
    
    def get_ids(self, a_frame):
        """
        Input:
        -self: Class instance
        -a_frame: a dataframe of peptides with accession
        
        Output: 
        - ids: List of protein IDs to query uniprot
        
        Description: 
        """
        ids = []
        for row in range(len(a_frame)):
            accession = a_frame["Accession"][row]
            if type(accession) != float:
                pipes = []
                pair = []
                check = 0
                for i in range(len(accession)):
                    if accession[i] == "|":
                        pair.append(i)
                        check += 1
                    else:
                        pass
                    if check == 2:
                        pipes.append(pair.copy())
                        pair = []
                        check = 0
                start = pipes[0][0] + 1
                end = pipes[0][1]
                code = accession[start:end]
                ids.append(code)
            else:
                pass
        return ids

    def get_seqs(self, id_list):
        """
        Input:
        -self: Class instance
        
        Output:
        -sequences: list of fasta sequences returned by uniprot sorted alphabetically by their corresponding codes
        -codes: list of codes of sequences returned by uniprot sorted alphabetically
        
        Description:
        """
        #Returns a list of sequences, sorted by accession code
        query = "+OR+id:".join(id_list)
        query = "id:" + query
        tries = 0
        uniprot = requests.get("https://www.uniprot.org/uniprot/", params = {'query':query, 'format':'FASTA'}, headers = {"contact":"jfeh0001@student.monash.edu"})
        while str(uniprot) != "<Response [200]>":
            try:
                tries <= 3
            except:
                print("Uniprot not responding")
            uniprot = requests.get("https://www.uniprot.org/uniprot/", params = {'query':query, 'format':'FASTA'}, headers = {"contact":"jfeh0001@student.monash.edu"})
        sequences = uniprot.text.split("\n>")
        if len(sequences[0]) > 0:
            sequences[0] = sequences[0].replace(">", "", 1)
            sequences.sort()
            codes = []
            for i in range(len(sequences)):
                sequences[i] = sequences[i].split("\n")
                codes.append(sequences[i][0].split("|")[1])
                sequences[i] = sequences[i][1:len(sequences[i])]
                sequences[i] = "".join(sequences[i])
        else:
            sequences = []
            codes = []
        return sequences, codes
    
    def get_csv(self):
        """
        Input:
        -self: Class instance
        
        Output:
        
        Description: Activated at button press. Allows user to choose a .csv file for analysis.
        """
        self.filepath = fd.askopenfilename(initialdir = "/",title = "Select file",filetypes = (("csv files","*.csv"),("all files","*.*")))
        print(self.filepath)
        name_start = 0
        name_end = 0
        for i in range(len(self.filepath)):
            if self.filepath[i] == "/":
                name_start = i
            elif self.filepath[i] == ".":
                name_end = i
            else:
                pass
        self.filename = self.filepath[name_start + 1:name_end]
        self.peptide_csv = self.filepath
        self.input_frame = pd.read_csv(self.peptide_csv)
        self.btn_process.pack(side = "top")
        self.lbl_is_csv['text'] = "File: " + self.filename
        
            
    def process_frame(self, a_frame):
        """
        Input:
        
        Output:
        
        Description:
        """
        self.process_progress.pack()
        self.master.update()
        self.clean_frame = self.clean(a_frame)#[0:100])
        self.peptide_frame, self.protein_frame = self.get_locations(self.clean_frame)
        pep_times = [0]*len(self.protein_frame)
        for i in range(len(self.protein_frame)):
            pep_times[i] = len(self.peptide_frame[self.peptide_frame["ID"] == self.protein_frame["Protein"][i]])
        self.protein_frame["Peptides"] = pep_times
        #Save summary of key information (pep_report)
        self.lbl_is_csv['text'] = self.filename + "Processed"
        self.summary_report()
        #Show initial overview figure
        figure = self.fig_locations_heat()
        self.canvas = FigureCanvasTkAgg(figure, self.master)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        #Remove loading bar        
        self.process_progress.destroy()
        #save frames
        #protein
        prot_save = open(self.filename + "_protein_frame.csv", "w")
        prot_save.write(self.protein_frame.to_csv())
        prot_save.close()
        #peptide
        pep_save = open(self.filename + "_peptide_frame.csv", "w")
        pep_save.write(self.peptide_frame.to_csv())
        pep_save.close()
        
    
    def summary_report(self):
        report_name = self.filename + "_report.txt"
        report_string = ""
        report_string += str("Number of peptides: " + str(len(self.peptide_frame)) + "\n")
        report_string += str("Number of proteins: " + str(len(self.protein_frame)) + "\n")
        report_string += str("Average peptides per protein: " + str(self.protein_frame["Peptides"].mean()) + "\n")
        report_string += str("Protein with most peptides: " + str(self.protein_frame["Protein"][self.protein_frame["Peptides"].idxmax()]) + "\n")
        report = open(report_name, "w")
        report.write(report_string)
        report.close()
    
    def get_locations(self, a_frame):
        """
        Input:
        
        Output:
        
        Description:
        """
        sort_frame = a_frame.sort_values(by = ["Accession"]).reset_index()
        seq_list = []
        all_codes = []
        subset = 0
        sleeper = 1
        #Function to get ids
        id_list = self.get_ids(sort_frame)
        sort_frame["ID"] = id_list
        #Protein retrieval##################################################
        #Fetch uniprot fasta files in batches of 50
        while subset < len(a_frame)-50:
            print(subset)
            sub_list = id_list[subset:subset+50]
            #Function to get sequences
            sequences, returned_codes = self.get_seqs(sub_list)
            index = 0
            #check all codes were returned from uniprot
            if len(sub_list) == len(returned_codes):
                pass
            else:
                missing = 0
                unique_ids = []
                last_val = None
                for value in sub_list:
                    if value != last_val:
                        unique_ids.append(value)
                        last_val = value
                    else:
                        pass
                while index < len(unique_ids):
                    if index == len(returned_codes):
                        returned_codes.append("No Sequence " + str(unique_ids[index]))
                        sequences.append("No Sequence " + str(unique_ids[index]))
                    elif unique_ids[index] == returned_codes[index]:
                        pass
                    else:
                        returned_codes = returned_codes[0:index] + ["No Sequence " + str(unique_ids[index])] + returned_codes[index:len(returned_codes)]
                        sequences = sequences[0:index] + ["No Sequence " + str(unique_ids[index])] + sequences[index:len(sequences)]
                    index += 1
            if subset == 0:
                pass
            else: 
                if all_codes[-1] == returned_codes[0]:
                    all_codes.pop()
                    seq_list.pop()
                else:
                    pass
            seq_list = seq_list + sequences
            all_codes = all_codes + returned_codes
            subset += 50
            self.process_progress['value'] = subset/len(a_frame)*100
            self.master.update()
            time.sleep(sleeper)
        sub_list = id_list[subset:len(id_list)]
        #Function to get sequences
        sequences, returned_codes = self.get_seqs(sub_list)
        if len(sub_list) == len(returned_codes):
                pass
        else:
            missing = 0
            unique_ids = []
            last_val = None
            index = 0
            for value in sub_list:
                if value != last_val:
                    unique_ids.append(value)
                    last_val = value
                else:
                    pass
            while index < len(unique_ids):
                if index == len(returned_codes):
                    returned_codes.append("No Sequence " + str(unique_ids[index]))
                    sequences.append("No Sequence " + str(unique_ids[index]))
                elif unique_ids[index] == returned_codes[index]:
                    pass
                else:
                    returned_codes = returned_codes[0:index] + ["No Sequence " + str(unique_ids[index])] + returned_codes[index:len(returned_codes)]
                    sequences = sequences[0:index] + ["No Sequence " + str(unique_ids[index])] + sequences[index:len(sequences)]
                index += 1
        if subset == 0:
                pass
        else: 
            if all_codes[-1] == returned_codes[0]:
                all_codes.pop()
                seq_list.pop()
            else:
                pass
        seq_list = seq_list + sequences
        all_codes = all_codes + returned_codes
        #Make proteins dataframe
        prot_frame = pd.DataFrame(all_codes, columns = ["Protein"])
        prot_frame["Sequence"] = seq_list
        #Peptide-protein matching###############################################################################
        loc_list = []
        rel_list = []
        rel_len = []
        row = 0
        sequence = 0
        offset = 0
        last_id = None
        errors = []
        for i in range(len(sort_frame)):
            pep = sort_frame["Peptide"][i]
            #check if current protein is a duplicate and increment offset by one if it is
            if id_list[i] == last_id:
                offset += 1
            else:
                pass
            last_id = id_list[i]
            if i-offset < len(seq_list):
                prot = seq_list[i - offset]
            else:
                pass
            pos = prot.find(pep)
            loc_list.append(pos)
            rel_list.append(pos/len(prot))
            rel_len.append(len(pep)/len(prot))
        sort_frame["Location"] = loc_list
        sort_frame["Relative Location"] = rel_list
        sort_frame["Relative Length"] = rel_len
        self.process_progress['value'] = 100 
        return (sort_frame, prot_frame)

    def protein_dist(self, protein):
        peptides = self.peptide_frame[self.peptide_frame["ID"] == protein]
        peptides.reset_index(drop = True, inplace = True)
        locations = [0]*100
        for pep in range(len(peptides)):
            start = int(peptides["Relative Location"][pep]//0.01)
            end = start + int(peptides["Relative Length"][pep]//0.01)
            for position in range(start, end + 1):
                locations[position] += 1
        return locations        
    
    #FIGURE PLOTTING####################################################################################################
    
    def fig_locations_heat(self):
        pos_freqs = plt.hist(self.peptide_frame["Relative Location"], bins = 100, range = [0,1])[0]
        pos_freqs_ox = plt.hist(self.peptide_frame[self.peptide_frame["PTM"]=='Oxidation (M)']["Relative Location"], bins = 100, range = [0,1])[0]
        pos_freqs_da = plt.hist(self.peptide_frame[self.peptide_frame["PTM"]=='Deamidation (NQ)']["Relative Location"], bins = 100, range = [0,1])[0]
        figure = Figure(figsize=(6, 6))
        ax = figure.subplots()
        sns.heatmap([pos_freqs, pos_freqs_ox, pos_freqs_da], yticklabels = ["All peptides", "Oxidation", "Deamidation"], cmap = "Reds", ax=ax)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_title("Peptide Distribution", fontsize = 20)
        ax.set_xlabel('Relative Peptide Position', fontsize=12)
        ax.set_ylabel('Protein', fontsize=15)
        return figure
    def plot_fig_locations_heat(self):
        self.canvas.get_tk_widget().destroy()
        figure = self.fig_locations_heat()
        self.canvas = FigureCanvasTkAgg(figure, self.master)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        
    def fig_prot_dist_heat(self):
        proteins = self.protein_selection.get()
        proteins = proteins.replace(" ", "")
        proteins = proteins.split(",")
        locations = []
        print(proteins)
        print(self.peptide_frame[self.peptide_frame["ID"] == proteins[1]])
        if len(proteins) > 0:
            #check all proteins are in list
            for i in proteins:
                if len(self.peptide_frame[self.peptide_frame["ID"] == i]) > 0:
                    pass
                else:
                    print("Error Message: Protein Not Found")
                locations.append(self.protein_dist(i))
            figure = Figure(figsize=(6, 6))
            ax = figure.subplots()
            chart = sns.heatmap(locations, cmap = "Reds", yticklabels = proteins, ax=ax)
            chart.set_yticklabels(chart.get_yticklabels(), rotation=0)
            ax.set_title("Peptide Distribution For Proteins", fontsize = 20)
            ax.set_xlabel('Relative Peptide Position', fontsize=12)
            ax.set_ylabel('Protein', fontsize=12)
            return figure
        else:
            #Print an error message
            pass
    def plot_fig_prot_dist_heat(self):
        self.canvas.get_tk_widget().destroy()
        figure = self.fig_prot_dist_heat()
        self.canvas = FigureCanvasTkAgg(figure, self.master)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        
    def fig_len_dist(self):
        len_dist = plt.hist(self.peptide_frame["Length"], bins = self.peptide_frame["Length"].max(), range = (1,self.peptide_frame["Length"].max()))[0]
        figure = Figure(figsize=(6, 6))
        ax = figure.subplots()
        #sns.lineplot(x = range(1,len(len_dist) + 1), y = len_dist, ax=ax)
        ax.hist(self.peptide_frame["Length"], bins = self.peptide_frame["Length"].max(), range = (1,self.peptide_frame["Length"].max()))
        ax.set_title("Peptide Length Distribution", fontsize = 20)
        ax.set_xlabel('Peptide Length', fontsize=12)
        ax.set_ylabel('Peptide Count', fontsize=12)
        return figure
    def plot_fig_len_dist(self):
        self.canvas.get_tk_widget().destroy()
        figure = self.fig_len_dist()
        self.canvas = FigureCanvasTkAgg(figure, self.master)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=1)
        
def main():
    root = tk.Tk()
    window = Window(root)
    root.mainloop()
  
if __name__ == "__main__":
    main()
