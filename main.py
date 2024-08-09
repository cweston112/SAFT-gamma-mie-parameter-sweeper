import numpy as np
import pandas as pd
from icecream import ic
import math
import tkinter as tk
import customtkinter

# Setting the colour theme for the GUI
customtkinter.set_default_color_theme("green")
# Setting the light appearance mode
customtkinter.set_appearance_mode("light")
# Defining a class for the single group parameters


class Single_Group:
    def __init__(self, group, number_of_spheres, shape_factor, sigma, sigma_born, lambda_r, lambda_a, epsilon,
                 number_of_association_group_types, number_of_hydrogen_sites, number_of_type_1_electron_sites,
                 number_of_type_2_electron_sites, number_of_self_bonding_sites, charge):
        self.group = group
        self.number_of_spheres = number_of_spheres
        self.shape_factor = shape_factor
        self.sigma = sigma
        self.sigma_born = sigma_born
        self.lambda_r = lambda_r
        self.lambda_a = lambda_a
        self.epsilon = epsilon
        self.number_of_association_group_types = number_of_association_group_types
        self.number_of_hydrogen_sites = number_of_hydrogen_sites
        self.number_of_type_1_electron_sites = number_of_type_1_electron_sites
        self.number_of_type_2_electron_sites = number_of_type_2_electron_sites
        self.number_of_self_bonding_sites = number_of_self_bonding_sites
        self.charge = charge
# Defining a class for the cross parameters


class Cross_Values:
    def __init__(self, group_1, group_2, cross_epsilon, cross_lambda_r):
        self.group_1 = group_1
        self.group_2 = group_2
        self.cross_epsilon = cross_epsilon
        self.cross_lambda_r = cross_lambda_r
        self.CR_sigma = None
        self.cross_epsilon_CR = None
        self.cross_lambda_r_CR = None

    def CR_lambda_R(self, single_group_1: Single_Group, single_group_2: Single_Group):
        if self.cross_lambda_r == "CR":
            self.cross_lambda_r = 3 + math.sqrt((single_group_1.lambda_r - 3)*(single_group_2.lambda_r - 3))
            self.cross_lambda_r_CR = True
        else:
            self.cross_lambda_r = self.cross_lambda_r
            self.cross_lambda_r_CR = False
        return self.cross_lambda_r

    def CR_epsilon(self, single_group_1: Single_Group, single_group_2: Single_Group):
        self.CR_sigma = (single_group_1.sigma * single_group_2.sigma) / 2
        if self.cross_epsilon == "CR":
            self.cross_epsilon = ((math.sqrt((single_group_1.sigma**3)*(single_group_2.sigma**3))/(self.CR_sigma ** 3)) *
                                  math.sqrt(single_group_1.epsilon*single_group_2.epsilon))
            self.cross_epsilon_CR = True
        else:
            self.cross_epsilon = self.cross_epsilon
            self.cross_epsilon_CR = False
        return self.cross_epsilon
# Defining a class for the association data


class Association_Values:
    def __init__(self, group_1, site_1, group_2, site_2, epsilon_association, bonding_volume):
        self.group_1 = group_1
        self.site_1 = site_1
        self.group_2 = group_2
        self.site_2 = site_2
        self.epsilon_association = epsilon_association
        self.bonding_volume = bonding_volume


# Extracting the single group data
single_group_data = pd.read_csv('Single component parameters (CSV).csv', index_col=0)
# Replacing the empty cells in the CSV with n/a
single_group_data = single_group_data.fillna("n/a")
# Extracting the group names
groups = single_group_data.index.values
N = groups.size
pure_groups_params = {}
# Creating a dictionary containing the single component groups parameters
for i in range(N):
    pure_groups_params[single_group_data.index[i]] = Single_Group(single_group_data.index[i],
                                                           single_group_data.iloc[i]["Number of spheres"],
                                                           single_group_data.iloc[i]["Shape factor"],
                                                           single_group_data.iloc[i]["Sigma (Angstroms)"],
                                                           single_group_data.iloc[i]["Sigma (Born)"],
                                                           single_group_data.iloc[i]["Lambda R"],
                                                           single_group_data.iloc[i]["Lambda A"],
                                                           single_group_data.iloc[i]["epsilon (like like)"],
                                                           single_group_data.iloc[i]
                                                           ["Number of association group types"],
                                                           single_group_data.iloc[i]["Number of hydrogen sites"],
                                                           single_group_data.iloc[i]
                                                           ["Number of electron sites (type 1)"],
                                                           single_group_data.iloc[i]
                                                           ["Number of electron sites (type 2)"],
                                                           single_group_data.iloc[i]["Number of self bonding sites"],
                                                           single_group_data.iloc[i]["Charge"])

# Extracting the cross component data
cross_data = pd.read_csv("Cross epsilon and lambda R (CSV).csv", index_col=[0,1])
# Extracting group names
cross_groups = cross_data.index.values
N_c = cross_groups.size
cross_groups_params = {}

for i in range(N_c):
    cross_groups_params[(cross_groups[i][0], cross_groups[i][1])] = Cross_Values(cross_groups[i][0], cross_groups[i][1],
                                                                                 cross_data.iloc[i]["Cross epsilon"],
                                                                                 cross_data.iloc[i]["Cross lambda R"])
    pure_group_0 = pure_groups_params[str(cross_groups[i][0])]
    pure_group_1 = pure_groups_params[str(cross_groups[i][1])]
    cross_groups_params[(cross_groups[i][0], cross_groups[i][1])].cross_epsilon = cross_groups_params[(cross_groups[i][0], cross_groups[i][1])].CR_epsilon(pure_group_0, pure_group_1)
    cross_groups_params[(cross_groups[i][0], cross_groups[i][1])].cross_lambda_r = cross_groups_params[(cross_groups[i][0], cross_groups[i][1])].CR_lambda_R(pure_group_0, pure_group_1)

    cross_groups_params[(cross_groups[i][1], cross_groups[i][0])] = Cross_Values(cross_groups[i][1], cross_groups[i][0],
                                                                                 cross_data.iloc[i]["Cross epsilon"],
                                                                                 cross_data.iloc[i]["Cross lambda R"])
    pure_group_0 = pure_groups_params[str(cross_groups[i][1])]
    pure_group_1 = pure_groups_params[str(cross_groups[i][0])]
    cross_groups_params[(cross_groups[i][1], cross_groups[i][0])].cross_epsilon = cross_groups_params[
        (cross_groups[i][1], cross_groups[i][0])].CR_epsilon(pure_group_0, pure_group_1)
    cross_groups_params[(cross_groups[i][1], cross_groups[i][0])].cross_lambda_r = cross_groups_params[
        (cross_groups[i][1], cross_groups[i][0])].CR_lambda_R(pure_group_0, pure_group_1)

# Extracting the association data
association_data = pd.read_csv("Association (CSV).csv", index_col=[0,1,2,3])
# Extracting group names
N_a = len(association_data.index)
association_groups = [["",""] for i in range(N_a)]
for i in range(N_a):
    association_groups[i][0] = association_data.index.values[i][0]
    association_groups[i][1] = association_data.index.values[i][2]
site_on_group_1 = ["" for i in range(N_a)]
site_on_group_2 = ["" for i in range(N_a)]

association_groups_params = {}
for j in range(N_a):
    site_on_group_1[j] = association_data.index.values[j][1]
    site_on_group_2[j] = association_data.index.values[j][3]

for i in range(N_a):
    association_groups_params[(association_groups[i][0], site_on_group_1[i], association_groups[i][1], site_on_group_2[i])] = Association_Values(
        association_groups[i][0], site_on_group_1[i], association_groups[i][1],
        site_on_group_2[i], association_data.iloc[i]["epsilon (assoc.)"],
        association_data.iloc[i]["Bonding volume (cubic angstroms)"])
    association_groups_params[(association_groups[i][1], site_on_group_2[i], association_groups[i][0], site_on_group_1[i])] = Association_Values(
        association_groups[i][1], site_on_group_2[i], association_groups[i][0],
        site_on_group_1[i], association_data.iloc[i]["epsilon (assoc.)"],
        association_data.iloc[i]["Bonding volume (cubic angstroms)"])
ic(len(association_groups_params))
# Now that I have all the data stored in structs - I can begin building the UI

class App():
    def __init__(self):
        # Initialising the app window and building a tab view
        self.root = customtkinter.CTk()
        self.root.geometry("1000x850")
        self.tabview = customtkinter.CTkTabview(master = self.root, width = 950, height = 800, fg_color="#ddfa9d")
        self.tabview.pack()
        # Creating an initial tab
        self.initial_tab = self.tabview.add("Initial")
        self.root.resizable(False, False)
        # Titling the tab window
        self.root.title("SAFT gamma MIE parameters")
        # Getting a list of all the pure groups we have data for
        self.group_options = groups
        # Creating lists contain the fields, tabs and groups created by loops in the program
        self.group_labels = []
        self.group_fields = []
        self.tabs = []
        self.pure_titles = []
        self.pure_bodies = []
        self.pure_association = []
        self.pure_association_out = ""
        self.cross_titles = []
        self.cross_bodies = []
        self.cross_association = []
        self.cross_association_out = ""

        def output_parameters():
            # Obtaining the selected number of groups
            self.num_groups_selected = int(self.num_groups_entry.get())
            # Obtaining the selected groups
            self.selected_groups = []
            for i in range(self.num_groups_selected):
                self.dummy = self.group_fields[i].get()
                self.selected_groups.append(self.dummy)
            # Creating a tab for each pure group
            for pure in self.selected_groups:
                self.dummy_1 = self.tabview.add(pure)
                self.tabs.append(self.dummy_1)
            # Creating a tab for each pair of groups
            self.group_pairs = [(a, b) for idx, a in enumerate(self.selected_groups) for b in self.selected_groups[idx + 1:]]
            for pair in self.group_pairs:
                self.dummy_2 = self.tabview.add(str(pair[0]) + " and " + str(pair[1]))
                self.tabs.append(self.dummy_2)
            # Outputting the parameter data for the pure components
            self.pass_1 = 0
            self.pass_2 = 0
            for i in range(len(self.selected_groups)):
                self.pure_out_title = "Group: "+ str(self.selected_groups[i])
                self.pure_out_text = ("Number of spheres: "+ str(pure_groups_params[self.selected_groups[i]].number_of_spheres) +
                                      "\n Shape factor: " + str(pure_groups_params[self.selected_groups[i]].shape_factor) +
                                      "\n Sigma (Angstroms): " + str(pure_groups_params[self.selected_groups[i]].sigma) +
                                      "\n Sigma (Born): " + str(pure_groups_params[self.selected_groups[i]].sigma_born) +
                                      "\n Lambda r: " + str(pure_groups_params[self.selected_groups[i]].lambda_r) +
                                      "\n Lambda a: " + str(pure_groups_params[self.selected_groups[i]].lambda_a) +
                                      "\n Epsilon: " + str(pure_groups_params[self.selected_groups[i]].epsilon) +
                                      "\n Number of association group types: " + str(pure_groups_params[self.selected_groups[i]].number_of_association_group_types) +
                                      "\n Number of hydrogen sites: " + str(pure_groups_params[self.selected_groups[i]].number_of_hydrogen_sites) +
                                      "\n Number of electron sites (type 1): " + str(pure_groups_params[self.selected_groups[i]].number_of_type_1_electron_sites) +
                                      "\n Number of electron sites (type 2): " + str(pure_groups_params[self.selected_groups[i]].number_of_type_2_electron_sites) +
                                      "\n Number of self association sites: " + str(pure_groups_params[self.selected_groups[i]].number_of_self_bonding_sites) +
                                      "\n Charge: " + str(pure_groups_params[self.selected_groups[i]].charge))
                try:
                    self.pure_association_out = ("\n\n SELF ASSOCIATION DATA:\n Site 1: " + str(association_groups_params[(self.selected_groups[i], "e1", self.selected_groups[i], "H")].site_1) +
                                                 "\n Site 2: " + str(association_groups_params[(self.selected_groups[i], "e1", self.selected_groups[i], "H")].site_2) +
                                                 "\n Epsilon (assoc.): " + str(association_groups_params[(self.selected_groups[i], "e1", self.selected_groups[i], "H")].epsilon_association) +
                                                 "\n Bonding volume (cubic angstroms): " + str(association_groups_params[(self.selected_groups[i], "e1", self.selected_groups[i], "H")].bonding_volume))
                except Exception:
                    self.pass_1 += self.pass_1
                    try:
                        self.pure_association_out = ("\n\n SELF ASSOCIATION DATA:\n Site 1: " + str(
                            association_groups_params[
                                (self.selected_groups[i], "e2", self.selected_groups[i], "H")].site_1) +
                                                     "\n Site 2: " + str(association_groups_params[(
                                self.selected_groups[i], "e2", self.selected_groups[i], "H")].site_2) +
                                                     "\n Epsilon (assoc.): " + str(association_groups_params[(
                                self.selected_groups[i], "e2", self.selected_groups[i], "H")].epsilon_association) +
                                                     "\n Bonding volume (cubic angstroms): " + str(
                                    association_groups_params[(
                                    self.selected_groups[i], "e2", self.selected_groups[i], "H")].bonding_volume))
                    except Exception:
                        self.pass_2 += self.pass_2
                        try:
                            self.pure_association_out = ("\n\n SELF ASSOCIATION DATA:\n Site 1: " + str(
                                association_groups_params[
                                    (self.selected_groups[i], "e*", self.selected_groups[i], "e*")].site_1) +
                                                         "\n Site 2: " + str(association_groups_params[(
                                    self.selected_groups[i], "e*", self.selected_groups[i], "e*")].site_2) +
                                                         "\n Epsilon (assoc.): " + str(association_groups_params[(
                                    self.selected_groups[i], "e*", self.selected_groups[i],
                                    "e*")].epsilon_association) +
                                                         "\n Bonding volume (cubic angstroms): " + str(
                                        association_groups_params[(
                                        self.selected_groups[i], "e*", self.selected_groups[i],
                                        "e*")].bonding_volume))
                        except Exception:
                            self.pure_association_out = "\n\n SELF ASSOCIATION DATA:\n No self association data "
                if self.pass_1 == 0:
                    try:
                        self.pure_association_out = self.pure_association_out + ("\n\n SELF ASSOCIATION DATA:\n Site 1: " + str(
                            association_groups_params[
                                (self.selected_groups[i], "e2", self.selected_groups[i], "H")].site_1) +
                                                     "\n Site 2: " + str(association_groups_params[(
                                    self.selected_groups[i], "e2", self.selected_groups[i], "H")].site_2) +
                                                     "\n Epsilon (assoc.): " + str(association_groups_params[(
                                    self.selected_groups[i], "e2", self.selected_groups[i], "H")].epsilon_association) +
                                                     "\n Bonding volume (cubic angstroms): " + str(
                                    association_groups_params[(
                                        self.selected_groups[i], "e2", self.selected_groups[i], "H")].bonding_volume))
                    except Exception:
                        self.pass_2 += self.pass_2

                self.pure_out_title_wid = customtkinter.CTkLabel(self.tabs[i], text=self.pure_out_title,
                                                                 font=("Times New Roman", 50))
                self.pure_out_title_wid.pack(side="top", fill="both", expand=True)
                self.pure_titles.append(self.pure_out_title_wid)
                self.pure_out_body_wid = customtkinter.CTkLabel(self.tabs[i], text=self.pure_out_text,
                                                                font=("Times New Roman", 20))
                self.pure_out_body_wid.pack(side="top", fill="both", expand=True)
                self.pure_bodies.append(self.pure_out_body_wid)
                self.pure_assoc_wid = customtkinter.CTkLabel(self.tabs[i], text=self.pure_association_out,
                                                             font=("Times New Roman", 20))
                self.pure_assoc_wid.pack(side="top", fill="both", expand=True)
                self.pure_association.append(self.pure_assoc_wid)

            self.pass_1 = 0
            self.pass_2 = 0

            for k, groups in enumerate(self.group_pairs):
                group1, group2 = groups
                ic(group1, group2, k)
                j = self.num_groups_selected + k
                self.group_out_title = "Groups: " + str(group1) + " - " + str(group2)
                self.group_out_text = ("Cross epsilon combining rules: " + str(cross_groups_params[(group1, group2)].cross_epsilon_CR) +
                                      "\n Cross epsilon: " + str(cross_groups_params[(group1, group2)].cross_epsilon) +
                                      "\n Cross lambda r combining rules: " + str(cross_groups_params[(group1, group2)].cross_lambda_r_CR) +
                                      "\n Cross lambda r: " + str(cross_groups_params[(group1, group2)].cross_lambda_r))

                try:
                    self.cross_association_out = ("\n\n CROSS ASSOCIATION DATA:\n Group 1: " + group1 + "\n Group 2: " + group2 +
                                                 "\n Site 1: " + str(association_groups_params[(group1, "e1", group2, "H")].site_1) +
                                                 "\n Site 2: " + str(association_groups_params[(group1, "e1", group2, "H")].site_2) +
                                                 "\n Epsilon (assoc.): " + str(association_groups_params[(group1, "e1", group2, "H")].epsilon_association) +
                                                 "\n Bonding volume (cubic angstroms): " + str(association_groups_params[(group1, "e1", group2, "H")].bonding_volume))
                except Exception:
                    self.pass_1 += self.pass_1
                    try:
                        self.cross_association_out = (
                                    "\n\n CROSS ASSOCIATION DATA:\n Group 1: " + group1 + "\n Group 2: " + group2 +
                                    "\n Site 1: " + str(association_groups_params[(group1, "e2", group2, "H")].site_1) +
                                    "\n Site 2: " + str(association_groups_params[(group1, "e2", group2, "H")].site_2) +
                                    "\n Epsilon (assoc.): " + str(
                                association_groups_params[(group1, "e2", group2, "H")].epsilon_association) +
                                    "\n Bonding volume (cubic angstroms): " + str(
                                association_groups_params[(group1, "e2", group2, "H")].bonding_volume))
                    except Exception:
                        self.pass_2 += self.pass_2
                        self.cross_association_out = "\n\n CROSS ASSOCIATION DATA:\n No cross association data "
                if self.pass_1 == 0:
                    try:
                        self.cross_association_out = self.cross_association_out + (
                                "\n\n CROSS ASSOCIATION DATA:\n Group 1: " + group1 + "\n Group 2: " + group2 +
                                "\n Site 1: " + str(association_groups_params[(group1, "e2", group2, "H")].site_1) +
                                "\n Site 2: " + str(association_groups_params[(group1, "e2", group2, "H")].site_2) +
                                "\n Epsilon (assoc.): " + str(
                            association_groups_params[(group1, "e2", group2, "H")].epsilon_association) +
                                "\n Bonding volume (cubic angstroms): " + str(
                            association_groups_params[(group1, "e2", group2, "H")].bonding_volume))
                    except Exception:
                        self.pass_2 += self.pass_2
                self.pass_1 = 0
                self.pass_2 = 0
                try:
                    self.cross_association_out = self.cross_association_out + ("\n\n CROSS ASSOCIATION DATA:\n Group 1: " + group1 + "\n Group 2: " + group2 +
                                                 "\n Site 1: " + str(association_groups_params[(group1, "H", group2, "e1")].site_1) +
                                                 "\n Site 2: " + str(association_groups_params[(group1, "H", group2, "e1")].site_2) +
                                                 "\n Epsilon (assoc.): " + str(association_groups_params[(group1, "H", group2, "e1")].epsilon_association) +
                                                 "\n Bonding volume (cubic angstroms): " + str(association_groups_params[(group1, "H", group2, "e1")].bonding_volume))
                except Exception:
                    self.pass_1 += self.pass_1
                    try:
                        self.cross_association_out = self.cross_association_out + (
                                    "\n\n CROSS ASSOCIATION DATA:\n Group 1: " + group1 + "\n Group 2: " + group2 +
                                    "\n Site 1: " + str(association_groups_params[(group1, "H", group2, "e2")].site_1) +
                                    "\n Site 2: " + str(association_groups_params[(group1, "H", group2, "e2")].site_2) +
                                    "\n Epsilon (assoc.): " + str(
                                association_groups_params[(group1, "H", group2, "e2")].epsilon_association) +
                                    "\n Bonding volume (cubic angstroms): " + str(
                                association_groups_params[(group1, "H", group2, "e2")].bonding_volume))
                    except Exception:
                        self.pass_2 += self.pass_2
                if self.pass_1 == 0:
                    try:
                        self.cross_association_out = self.cross_association_out + (
                                "\n\n CROSS ASSOCIATION DATA:\n Group 1: " + group1 + "\n Group 2: " + group2 +
                                "\n Site 1: " + str(association_groups_params[(group1, "H", group2, "e2")].site_1) +
                                "\n Site 2: " + str(association_groups_params[(group1, "H", group2, "e2")].site_2) +
                                "\n Epsilon (assoc.): " + str(
                            association_groups_params[(group1, "H", group2, "e2")].epsilon_association) +
                                "\n Bonding volume (cubic angstroms): " + str(
                            association_groups_params[(group1, "H", group2, "e2")].bonding_volume))
                    except Exception:
                        self.pass_2 += self.pass_2
                ic(self.cross_association_out)

                self.group_out_title_wid = customtkinter.CTkLabel(self.tabs[j], text=self.group_out_title,
                                                                 font=("Times New Roman", 50))

                self.group_out_title_wid.pack(side="top", fill="both", expand=True)
                self.cross_titles.append(self.group_out_title_wid)

                self.group_out_body_wid = customtkinter.CTkLabel(self.tabs[j], text=self.group_out_text,
                                                                 font=("Times New Roman", 20))
                self.group_out_body_wid.pack(side="top", fill="both", expand=True)
                self.cross_bodies.append(self.group_out_body_wid)
                self.cross_assoc_wid = customtkinter.CTkLabel(self.tabs[j], text=self.cross_association_out,
                                                             font=("Times New Roman", 20))
                self.cross_assoc_wid.pack(side="top", fill="both", expand=True)
                self.cross_association.append(self.cross_assoc_wid)

        def num_groups_submit():
            self.num_groups_selected = int(self.num_groups_entry.get())
            for i in range(self.num_groups_selected):
                self.dummy_label = customtkinter.CTkLabel(self.initial_tab, text = "Group " + str(i+1) + ":",
                                                          font = ("Times New Roman", 20))
                self.dummy_label.grid(row = i+1, column = 0, padx = 10, pady = 10)
                self.group_labels.append(self.dummy_label)
                self.dummy_field = customtkinter.CTkComboBox(self.initial_tab,values=self.group_options)
                self.dummy_field.grid(row = i+1, column = 1, padx = 10, pady = 10)
                self.group_fields.append(self.dummy_field)
            self.group_submit = customtkinter.CTkButton(self.initial_tab, text = "submit", command = output_parameters)
            self.group_submit.grid(row = self.num_groups_selected+1, column = 1, padx = 10, pady = 10)



        self.num_groups_label = customtkinter.CTkLabel(self.initial_tab,
                                                       text="Please enter the number of groups you wish to consider: ",
                                                       font=("Times New Roman", 20))
        self.num_groups_label.grid(row=0, column=0, padx=10, pady=10, sticky=customtkinter.W)
        self.num_groups_entry = customtkinter.CTkEntry(self.initial_tab, placeholder_text="(integer)")
        self.num_groups_entry.grid(row=0, column=1, padx=10, pady=10, sticky=customtkinter.W)
        self.num_groups_button = customtkinter.CTkButton(self.initial_tab, text = "submit", command = num_groups_submit)
        self.num_groups_button.grid(row=0, column=2, padx=10, pady=10, sticky=customtkinter.W)

        self.root.mainloop()
        return


if __name__ == "__main__":
    App()
