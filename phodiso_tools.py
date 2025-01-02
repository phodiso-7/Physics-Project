#!/usr/bin/env python
# coding: utf-8
#
import ROOT
from ROOT import TLorentzVector
import pylhe
import pandas as pd
from tabulate import tabulate
import numpy as np
from functools import reduce
try:
    import uproot3
except Exception as e:
    print(f"Error importing uproot3: {e}, skipping functionalities.")
try:
    import awkward as ak
except ImportError:
    print("awkward is not installed, skipping functionalities.")
import matplotlib.pyplot as plt
import warnings
import math
import ast

ROOT.gErrorIgnoreLevel = ROOT.kError  ##supress warnings

# list of functions : open_root, cos_cut , compute_delta_R_jets , Delta_R_cut_jets , PT_cut_jets

#def list_of_functions(a=None):
    #return ['open_root', 'cos_cut' , 'compute_delta_R_jets' , 'Delta_R_cut_jets' , 'PT_cut_jets']
    
# k  


def list_of_functions(ras=None):
    a = '-> open_root(directory)'
    b = '-> open_plot(masses = masses_list, signal_paths = signal_paths_list, background_paths = background_path_list, signal_dfs=signal_dfs_list , background_dfs= background_dfs_list)'
    c = '-> cos_cut(df, pos_threshold, neg_threshold)'
    d = '-> compute_delta_R_jets(df_dijets)'
    e = '-> Delta_R_cut_jets(df_dijet_delta_R_computed, delta_R_criteria)'
    f = '-> PT_cut_jets(df, PT_criteria,condition=greater/smaller _string)'
    g = '-> PT_cut_leptons(df, PT_criteria,condition=greater/smaller _string)'
    all_stuff = [a, b, c, d, e, f,g]
    
    # Print each item with a newline in between
    print("\n".join(all_stuff))
    

    
def save_to_csv(df, path):
    # Save the DataFrame to a CSV file
    df.to_csv(f'{path}.csv', index=False)  # Use index=False to exclude row indices if not needed
    print("File has been saved successfully.")

def open_root_csv(path):
    df = pd.read_csv(path)
    def convert_to_list(element):
        # Check if the element is a string and appears to be a list
        if isinstance(element, str) and element.startswith('[') and element.endswith(']'):
            try:
                # Convert string to list using ast.literal_eval
                converted = ast.literal_eval(element)
                # Ensure all elements are numbers (int or float)
                if isinstance(converted, list) and all(isinstance(x, (int, float)) for x in converted):
                    return converted
            except (ValueError, SyntaxError):
                pass
        return element  # Return the original element if conversion fails
    # Apply the conversion function to all elements in the DataFrame
    df_converted = df.applymap(convert_to_list)
    return df_converted



def open_root(file_path, tree_name="Delphes"):
    """Extracts lepton and jet data from the given ROOT file."""
    # Open the ROOT file
    warnings.filterwarnings("ignore")
    file = ROOT.TFile.Open(file_path)
    
    # Check if the file is open
    if not file or file.IsZombie():
        raise Exception(f"Failed to open file {file_path}.")
    
    # Access the tree
    tree = file.Get(tree_name)
    if not tree:
        raise Exception(f"Failed to get the tree {tree_name} from file {file_path}.")
        
    #Reader for Particles
    reader_particles = ROOT.TTreeReader(tree)
    particles_pid_reader = ROOT.TTreeReaderArray('int')(reader_particles, "Particle.PID") #maybe int
    particles_mass_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.Mass") # Adjust the type if necessary
    particles_energy_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.E")  # E is often used for total energy
    particles_eta_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.Eta")
    particles_pt_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.PT")
    particles_phi_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.Phi")

    # Reader for electrons
    reader_electron = ROOT.TTreeReader(tree)
    electron_pt_reader = ROOT.TTreeReaderArray('float')(reader_electron, "Electron.PT")
    electron_phi_reader = ROOT.TTreeReaderArray('float')(reader_electron, "Electron.Phi")
    electron_eta_reader = ROOT.TTreeReaderArray('float')(reader_electron, "Electron.Eta")
    electron_charge_reader = ROOT.TTreeReaderArray('int')(reader_electron, "Electron.Charge")

    # Reader for muons
    reader_muon = ROOT.TTreeReader(tree)
    muon_pt_reader = ROOT.TTreeReaderArray('float')(reader_muon, "Muon.PT")
    muon_phi_reader = ROOT.TTreeReaderArray('float')(reader_muon, "Muon.Phi")
    muon_eta_reader = ROOT.TTreeReaderArray('float')(reader_muon, "Muon.Eta")
    muon_charge_reader = ROOT.TTreeReaderArray('int')(reader_muon, "Muon.Charge")

    # Reader for jets
    reader_jet = ROOT.TTreeReader(tree)
    jet_pt_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.PT")
    jet_phi_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.Phi")
    jet_eta_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.Eta")
    jet_mass_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.Mass")
    jet_tau_tagged_reader = ROOT.TTreeReaderArray('int')(reader_jet, "Jet.TauTag")
    jet_b_tagged_reader = ROOT.TTreeReaderArray('int')(reader_jet, "Jet.BTag")
    
    # Reader for Missing ET
    missing_et_reader = ROOT.TTreeReader(tree)
    missing_et_met_reader = ROOT.TTreeReaderArray('float')(missing_et_reader, "MissingET.MET")
    missing_et_phi_reader = ROOT.TTreeReaderArray('float')(missing_et_reader, "MissingET.Phi")
    missing_et_eta_reader = ROOT.TTreeReaderArray('float')(missing_et_reader, "MissingET.Eta")


    # Initialize lists to store data
    particles_data =[]
    electron_data = []
    muon_data = []
    jet_data = []
    missing_et_data = []
    
# Loop over events for particles
    event_number = 0
    while reader_particles.Next():
        # Read data from the tree
        particles_pid = [pid for pid in particles_pid_reader]
        particles_mass = [mass for mass in particles_mass_reader]
        particles_energy = [energy for energy in particles_energy_reader]
        particles_eta = [eta for eta in particles_eta_reader]
        particles_pt = [pt for pt in particles_pt_reader]
        particles_phi = [phi for phi in particles_phi_reader]
        
        # Append the data for each event
        particles_data.append({
            'Event': event_number,
            'Particle_PID': particles_pid,
            'Particle_Mass': particles_mass,
            'Particle_E': particles_energy,
            'Particle_Eta': particles_eta,
            'Particle_PT': particles_pt,
            'Particle_Phi': particles_phi
        })
        # Increment the event counter
        event_number += 1
    
    # Loop over events for electrons
    event_number = 0
    while reader_electron.Next():
        electron_pts = [pt for pt in electron_pt_reader]
        electron_phis = [phi for phi in electron_phi_reader]
        electron_etas = [eta for eta in electron_eta_reader]
        electron_charges = [charge for charge in electron_charge_reader]
        electron_thetas = [2 * np.arctan(np.exp(-eta)) for eta in electron_etas]
        electron_data.append({
            'Event': event_number,
            'Electron_PT': electron_pts,
            'Electron_Phi': electron_phis,
            'Electron_Eta': electron_etas,
            'Electron_Charge': electron_charges,
            'Electron_Theta': electron_thetas
        })
        event_number += 1

    # Loop over events for muons
    event_number = 0
    while reader_muon.Next():
        muon_pts = [pt for pt in muon_pt_reader]
        muon_phis = [phi for phi in muon_phi_reader]
        muon_etas = [eta for eta in muon_eta_reader]
        muon_charges = [charge for charge in muon_charge_reader]
        muon_thetas = [2 * np.arctan(np.exp(-eta)) for eta in muon_etas]
        muon_data.append({
            'Event': event_number,
            'Muon_PT': muon_pts,
            'Muon_Phi': muon_phis,
            'Muon_Eta': muon_etas,
            'Muon_Charge': muon_charges,
            'Muon_Theta': muon_thetas
        })
        event_number += 1

    # Loop over events for jets
    event_number = 0
    while reader_jet.Next():
        jet_pts = [pt for pt in jet_pt_reader]
        jet_phis = [phi for phi in jet_phi_reader]
        jet_etas = [eta for eta in jet_eta_reader]
        jet_masses = [mass for mass in jet_mass_reader]
        jet_tau_tagged = [tag for tag in jet_tau_tagged_reader]
        jet_b_tagged = [b for b in jet_b_tagged_reader]
        jet_data.append({
            'Event': event_number,
            'Jet_PT': jet_pts,
            'Jet_Phi': jet_phis,
            'Jet_Eta': jet_etas,
            'Jet_Mass': jet_masses,
            'Jet_TauTagged': jet_tau_tagged,
            'Jet_BTagged': jet_b_tagged
        })
        event_number += 1
        
    # Loop over events for Missing ET
    event_number = 0
    while missing_et_reader.Next():
        missing_met_values = [met for met in missing_et_met_reader]
        missing_phi_values = [phi for phi in missing_et_phi_reader]
        missing_eta_values = [eta for eta in missing_et_eta_reader]
        missing_et_data.append({
            'Event': event_number,
            'MissingET_MET': missing_met_values,
            'MissingET_Phi': missing_phi_values,
            'MissingET_Eta': missing_eta_values
        })
        event_number += 1

    # Convert to DataFrames
    df_particles = pd.DataFrame(particles_data)
    df_electrons = pd.DataFrame(electron_data)
    df_muons = pd.DataFrame(muon_data)
    df_jets = pd.DataFrame(jet_data)
    df_missing_et = pd.DataFrame(missing_et_data)


    # Merge the DataFrames on the 'Event' column
    df_merged = reduce(lambda left, right: pd.merge(left, right, on='Event', how='outer'), [df_particles,df_electrons, df_muons, df_jets,df_missing_et])

    # Combine PTs of electrons and muons, handling cases where lists may be empty
    df_merged['Combined_PT'] = df_merged['Electron_PT'].apply(lambda x: x if x else []) + df_merged['Muon_PT'].apply(lambda x: x if x else [])
    df_merged['Combined_Charge'] = df_merged['Electron_Charge'].apply(lambda x: x if x else []) + df_merged['Muon_Charge'].apply(lambda x: x if x else [])

    # Calculate Lepton_Theta as the combined theta from electrons and muons
    df_merged['Electron_Theta'] = df_merged['Electron_Theta'].apply(lambda x: x if x else [])
    df_merged['Electron_cos_theta'] = df_merged['Electron_Theta'].apply(lambda theta: np.cos(theta) if theta is not None else None)
    df_merged['Muon_Theta'] = df_merged['Muon_Theta'].apply(lambda x: x if x else [])
    df_merged['Muon_cos_theta'] = df_merged['Muon_Theta'].apply(lambda theta: np.cos(theta) if theta is not None else None)
    df_merged['Combined_Theta'] = df_merged['Electron_Theta'].apply(lambda x: x if x else []) + df_merged['Muon_Theta'].apply(lambda x: x if x else [])
    
    #Fix MET, it was made of lists
    df_merged['MissingET_MET'] = df_merged['MissingET_MET'].apply(
    lambda x: x[0] if isinstance(x, list) and x else np.nan
)
    
    df_merged['Combined_cos_theta'] = df_merged['Electron_cos_theta'].apply(list) + df_merged['Muon_cos_theta'].apply(list)

    df_merged['Combined_PT_Theta'] = df_merged.apply(
        lambda row: list(zip(row['Combined_PT'], row['Combined_Theta'])) if row['Combined_PT'] and row['Combined_Theta'] else [], axis=1
    )
    
    df_merged['Combined_PT_Charge'] = df_merged.apply(
        lambda row: list(zip(row['Combined_PT'], row['Combined_Charge'])) if row['Combined_PT'] and row['Combined_Charge'] else [], axis=1
    )
    
    # Sort by PT in descending order, keeping the corresponding Theta and Charge values aligned
    df_merged['Sorted_PT_Theta'] = df_merged['Combined_PT_Theta'].apply(
        lambda pairs: sorted(pairs, key=lambda x: x[0], reverse=True)
    )
    
    df_merged['Sorted_PT_Charge'] = df_merged['Combined_PT_Charge'].apply(
        lambda pairs: sorted(pairs, key=lambda x: x[0], reverse=True)
    )
    
    # Extract leading and subleading Max_PT, Max_Theta, and Charge based on the sorted order
    df_merged['Leading_lep_PT'] = df_merged['Sorted_PT_Theta'].apply(lambda pairs: pairs[0][0] if pairs else None)
    df_merged['Leading_lep_Theta'] = df_merged['Sorted_PT_Theta'].apply(lambda pairs: pairs[0][1] if pairs else None)
    df_merged['Leading_lep_Charge'] = df_merged['Sorted_PT_Charge'].apply(lambda pairs: pairs[0][1] if pairs else None)
    
    df_merged['Sub_leading_lep_PT'] = df_merged['Sorted_PT_Theta'].apply(
        lambda pairs: pairs[1][0] if len(pairs) > 1 else None
    )
    df_merged['Sub_leading_lep_Theta'] = df_merged['Sorted_PT_Theta'].apply(
        lambda pairs: pairs[1][1] if len(pairs) > 1 else None
    )
    df_merged['Sub_leading_lep_Charge'] = df_merged['Sorted_PT_Charge'].apply(
        lambda pairs: pairs[1][1] if len(pairs) > 1 else None
    )
    
    ###Separating the data by taggings
    def separate_tagged_untagged_jets(df):
        
    # Extract the necessary columns
        jet_pts = df['Jet_PT']
        jet_eta = df['Jet_Eta']
        jet_phi = df['Jet_Phi']
        jet_tagged = df['Jet_TauTagged']
        
        # Initialize lists to store the separated jet data
        tagged_pts = []
        tagged_eta = []
        tagged_phi = []
        
        untagged_pts = []
        untagged_eta = []
        untagged_phi = []
        
        # Separate the data based on the tagging
        for pts, eta, phi, tagged in zip(jet_pts, jet_eta, jet_phi, jet_tagged):
            if tagged == 1:  # Check for tagged (1)
                tagged_pts.append(pts)
                tagged_eta.append(eta)
                tagged_phi.append(phi)
            else:  # Check for untagged (0)
                untagged_pts.append(pts)
                untagged_eta.append(eta)
                untagged_phi.append(phi)
        
        # Update DataFrame with separated data
        df['Jet_TauTagged_PT'] = [
            [pt for pt, tag in zip(pts, tags) if tag == 1]
            for pts, tags in zip(jet_pts, jet_tagged)
        ]
        
        df['Jet_TauTagged_Eta'] = [
            [eta for eta, tag in zip(etas, tags) if tag == 1]
            for etas, tags in zip(jet_eta, jet_tagged)
        ]
        
        df['Jet_TauTagged_Phi'] = [
            [phi for phi, tag in zip(phis, tags) if tag == 1]
            for phis, tags in zip(jet_phi, jet_tagged)
        ]
        
        df['Jet_UT_PT'] = [
            [pt for pt, tag in zip(pts, tags) if tag == 0]
            for pts, tags in zip(jet_pts, jet_tagged)
        ]
        
        df['Jet_UT_Eta'] = [
            [eta for eta, tag in zip(etas, tags) if tag == 0]
            for etas, tags in zip(jet_eta, jet_tagged)
        ]
        
        df['Jet_UT_Phi'] = [
            [phi for phi, tag in zip(phis, tags) if tag == 0]
            for phis, tags in zip(jet_phi, jet_tagged)
        ]
        
        return df

    # Apply the function to the DataFrame
    
    df_merged = separate_tagged_untagged_jets(df_merged)
    
    ## Leading Jet_UT_PT :
    
    #Is it hadronic or leptonic tau :
    
    df_merged['decay_mode'] = df_merged['Jet_TauTagged_PT'].apply(lambda x: 1 if len(x) > 0 else 0)
    
    # Calculate the leading (maximum) PT
    df_merged["Leading_Jet_UT_PT"] = df_merged["Jet_UT_PT"].apply(lambda x: max(x) if isinstance(x, list) and x else np.nan)

    # Calculate the subleading (second largest) PT directly
    df_merged["Subleading_Jet_UT_PT"] = df_merged["Jet_UT_PT"].apply(lambda x: sorted(x, reverse=True)[1] if isinstance(x, list) and len(x) > 1 else np.nan)
    
    ###Masses of the Tau Tagged
    # Zip Jet_Mass and Jet_TauTagged columns
    df_merged['Jet_TauTagged_Mass'] = df_merged.apply(lambda row: [mass for mass, tag in zip(row['Jet_Mass'], row['Jet_TauTagged']) if tag == 1], axis=1)
    df_merged['Jet_UT_Mass'] = df_merged.apply(lambda row: [mass for mass, tag in zip(row['Jet_Mass'], row['Jet_TauTagged']) if tag == 0], axis=1)
    
    # Calculate Jet_TauTagged_ET using the formula
    df_merged['Jet_TauTagged_ET'] = df_merged.apply(lambda row: [(mass**2 + pt**2)**0.5 for mass, pt in zip(row['Jet_TauTagged_Mass'], row['Jet_TauTagged_PT'])], axis=1) 


    
    ### Cos theta of the leading tau tagged jet
    
    df_merged['TauTaggedTheta'] = df_merged['Jet_TauTagged_Eta'].apply(lambda etas: [ 2 * np.arctan(np.exp(-eta)) for eta in etas])
    
    df_merged['TauTaggedCos_Theta'] = df_merged['TauTaggedTheta'].apply(lambda thetas: [np.cos(theta) for theta in thetas])
    
    # Combine TauTaggedPT with TauTaggedCos_Theta using zip to pair PT with cos(theta)
    df_merged['TauTagged_PT_CosTheta'] = df_merged.apply(lambda row: list(zip(row['Jet_TauTagged_PT'], row['TauTaggedCos_Theta'])), axis=1)
        
    # Find the Leading Tau Tagged PT
    df_merged['Leading_TauTaggedPT'] = df_merged['Jet_TauTagged_PT'].apply(lambda pts: max(pts) if pts else None)
    
    # Find the corresponding Leading Cos Theta of the Leading Tau Tagged Jet
    df_merged['LeadingTauTaggedCos_Theta'] = df_merged['TauTagged_PT_CosTheta'].apply(lambda pairs: max(pairs, key=lambda x: x[0])[1] if pairs else None) 
    
    #Combining the PT'S and ET'S
    df_merged['TauTagged_PT_ET'] = df_merged.apply(lambda row: list(zip(row['Jet_TauTagged_PT'], row['Jet_TauTagged_ET'])), axis=1)
    
    #Now finding the et of the leading tau jet
    df_merged['LeadingTauTagged_ET'] = df_merged['TauTagged_PT_ET'].apply(lambda pairs: max(pairs, key=lambda x: x[0])[1] if pairs else None) 
    
    df_merged['Sub_LeadingTauTagged_ET'] = df_merged['TauTagged_PT_ET'].apply(
    lambda pairs: min(pairs, key=lambda x: x[0])[1] if len(pairs) > 1 else None
)

    
####Jet Invariant mass 
    df = df_merged
#     df['Jet_UT_PX'] = df.apply(lambda row: [ a * np.cos(b) for a, b  in zip(row['Jet_UT_PT'], row['Jet_UT_Phi'] )], axis=1)
#     df['Jet_UT_PY'] = df.apply(lambda row: [ a * np.sin(b) for a, b  in zip(row['Jet_UT_PT'], row['Jet_UT_Phi'] )], axis=1)
#     df['Jet_UT_PZ'] = df.apply(lambda row: [ a * np.sinh(b) for a, b  in zip(row['Jet_UT_PT'], row['Jet_UT_Phi'] )], axis=1)

#     df['Jet_UT_E'] = df.apply(lambda row: [ np.sqrt(a**2 + b**2 + c**2) for a, b , c in zip(row['Jet_UT_PT'], row['Jet_UT_PZ'],row['Jet_UT_Mass'] )], axis=1)

   # This will only reduce it to four or more events, then only selecte the top for

    
    df = df[df['Jet_UT_PT'].apply(len)>=4]

    df['Jet_UT_PT'] = df['Jet_UT_PT'].apply(lambda x : x[:4])




#attempt two
    df['Jet_UT_PX'] = df.apply(lambda row: [ a * np.sin(b) for a, b  in zip(row['Jet_UT_PT'], row['Jet_UT_Phi'] )], axis=1)
    df['Jet_UT_PY'] = df.apply(lambda row: [ a * np.cos(b) for a, b  in zip(row['Jet_UT_PT'], row['Jet_UT_Phi'] )], axis=1)
    df['Jet_UT_PZ'] = df.apply(lambda row: [ a * np.sinh(b) for a, b   in zip(row['Jet_UT_PT'], row['Jet_UT_Eta'] )], axis=1)
#     df['Jet_UT_PZ'] = df.apply(lambda row: [ np.sqrt(a**2 + c**2) * np.sinh(b) for a, b ,c   in zip(row['Jet_UT_PT'], row['Jet_UT_Eta'] , row['Jet_UT_PT'] )], axis=1)

    df['Jet_UT_E'] = df.apply(lambda row: [ np.sqrt(a**2 + b**2 + c**2) for a, b , c in zip(row['Jet_UT_PT'], row['Jet_UT_PZ'],row['Jet_UT_Mass'] )], axis=1)

#     df['Jet_UT_E'] = df.apply(lambda row: [ np.cosh(b)* np.sqrt(a**2 + c**2) for a, b , c in zip(row['Jet_UT_PT'], row['Jet_UT_Eta'],row['Jet_UT_Mass'] )], axis=1)

    # Handle empty lists by returning 0 for sum if the list is empty
    df['energy_total'] = df['Jet_UT_E'].apply(lambda x: sum(x) if x else 0)
    df['px_total'] = df['Jet_UT_PX'].apply(lambda x: sum(x) if x else 0)
    df['py_total'] = df['Jet_UT_PY'].apply(lambda x: sum(x) if x else 0)
    df['pz_total'] = df['Jet_UT_PZ'].apply(lambda x: sum(x) if x else 0)

    # Calculate the invariant mass using the formula: M = (E^2 - (px^2 + py^2 + pz^2)) ^0.5
    df['jet_invariant_mass'] = df.apply(
        lambda row: np.sqrt(row['energy_total']**2 - (row['px_total']**2 + row['py_total']**2 + row['pz_total']**2)) 
        if row['energy_total'] > 0 else np.nan, axis=1)
    
    def calculate_kinematic_variables(row):
        total_vector = TLorentzVector()

        # Iterate over each jet in the list
        for pt, eta, phi, energy in zip(row['Jet_UT_PT'], row['Jet_UT_Eta'], row['Jet_UT_Phi'], row['Jet_UT_E']):
            jet_vector = TLorentzVector()
            # Set the Lorentz vector using (pt, eta, phi, energy)
            jet_vector.SetPtEtaPhiE(pt, eta, phi, energy)
            total_vector += jet_vector

        # Extract the kinematic variables
        invariant_mass = total_vector.M()
        theta = total_vector.Theta()  # Polar angle in radians
        eta = total_vector.Eta()      # Pseudorapidity
        pt = total_vector.Pt()        # Transverse momentum
        phi = total_vector.Phi()      # Azimuthal angle in radians
        energy = total_vector.E()     # Energy

        # Return the results as a series
        return pd.Series({
            'M_jets': invariant_mass,
            'theta_jets': theta,
            'eta_jets': eta,
            'pt_jets': pt,
            'phi_jets': phi,
            'energy_jets': energy
        })

    # Apply the function to each row and create new columns for M, theta, eta, pt, and phi
    df[['M_jets', 'theta_jets', 'eta_jets', 'pt_jets', 'phi_jets' , 'energy_jets']] = df.apply(calculate_kinematic_variables, axis=1)

    
    df_merged = df
    
    #### Counting leptons and jets
    
    count_lep = [len(k) for k in df_merged["Combined_PT"]]
    count_jet = [len(k) for k in df_merged["Jet_PT"]]

    df_merged["Lepton Count"] = count_lep
    df_merged["Jet Count"] = count_jet
    
    
    
    # Drop the columns that are no longer needed
    drop_cols = [
        "Sorted_PT_Theta", "Combined_PT_Theta", "Combined_Theta",
        "Electron_Eta", "Muon_Eta",
        "Sorted_PT_Charge",'TauTagged_PT_CosTheta','TauTagged_PT_ET'
    ]
    df_merged = df_merged.drop(columns=drop_cols)
    

    

    # Clean up
    file.Close()
    

    return df_merged 


#file_path = '/home/physics/MG5_aMC_v2_9_4/6septau_background_tau_tau/Events/run_01/tag_1_delphes_events.root'

#df = open_root(file_path)

#df [ df['Jet_TauTagged_PT'].apply(len)==2] [['Event','Jet_PT', 'MissingET_MET', 'Jet_TauTagged_PT', 'Jet_TauTagged_Eta', 'Jet_UT_PT','TauTaggedCos_Theta', 'Leading_TauTaggedPT','LeadingTauTaggedCos_Theta', 'Lepton Count', 'Jet Count']]
    


def open_plot(masses, signal_paths=None, background_paths=None, signal_dfs=None, background_dfs=None):
    """
    Plots histograms of PT, Theta, and cos(Theta) for signal and background data.

    Args:
    - masses: List of mass values used in titles.
    - signal_paths: List of paths to the signal ROOT files (optional if signal_dfs is provided).
    - background_paths: List of paths to the background ROOT files (optional if background_dfs is provided).
    - signal_dfs: List of DataFrames for signal data (optional if signal_paths is provided).
    - background_dfs: List of DataFrames for background data (optional if background_paths is provided).
    """
    
    if signal_dfs is not None:
        # Use provided DataFrames for signal data
        dfs = signal_dfs
    elif signal_paths is not None:
        # Extract data for signal files
        dfs = [open_root(signal_path) for signal_path in signal_paths]
    else:
        raise ValueError("Either signal_paths or signal_dfs must be provided.")

    if background_dfs is not None:
        # Use provided DataFrames for background data
        df_background = pd.concat(background_dfs, ignore_index=True)
        background_provided = True
    elif background_paths is not None:
        # Extract data for background files and concatenate
        background_dfs = [open_root(file_path) for file_path in background_paths]
        df_background = pd.concat(background_dfs, ignore_index=True)
        background_provided = True
    else:
        background_provided = False

    # Iterate over the extracted signal DataFrames and plot histograms
    for idx, df in enumerate(dfs):
        fig, ax3 = plt.subplots()
        fig, ax4 = plt.subplots()
        fig, ax5 = plt.subplots()
        fig, ax6 = plt.subplots()
        fig, ax7 = plt.subplots()
        fig, ax8 = plt.subplots()
        #fig, ax9 = plt.subplots()


        # Plot Max_PT for signal and background
        ax3.hist(df['Leading_lep_PT'], bins=60, edgecolor='red', histtype='step', label='Signal')  # signal
        if background_provided:
            ax3.hist(df_background['Leading_lep_PT'], bins=60, edgecolor='black', histtype='step', label='Background')  
        ax3.set_xlabel(r'$p_T$')
        ax3.set_ylabel('Frequency')
        ax3.set_title(f'PTs of the Leading Leptons for MHP = {masses[idx]}')
        ax3.legend()

        # Plot Second_Max_PT for signal and background
        ax4.hist(df['Sub_leading_lep_PT'], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax4.hist(df_background['Sub_leading_lep_PT'], bins=60, edgecolor='black', histtype='step', label='Background')
        ax4.set_xlabel(r'$p_T$')
        ax4.set_ylabel('Frequency')
        ax4.set_title(f'PTs of the Sub-Leading Leptons for MHP = {masses[idx]}')
        ax4.legend()

        # Plot Max_Theta for signal and background
        ax5.hist(df['Leading_lep_Theta'], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax5.hist(df_background['Leading_lep_Theta'], bins=60, edgecolor='black', histtype='step', label='Background')
        ax5.set_xlabel(r'$\theta$')
        ax5.set_ylabel('Frequency')
        ax5.set_title(f'Thetas of the Leading Leptons for MHP = {masses[idx]}')
        ax5.legend()

        # Plot Second_Max_Theta for signal and background
        ax6.hist(df['Sub_leading_lep_Theta'] , bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax6.hist(df_background['Sub_leading_lep_Theta'], bins=60, edgecolor='black', histtype='step', label='Background')
        ax6.set_xlabel(r'$\theta$')
        ax6.set_ylabel('Frequency')
        ax6.set_title(f'Thetas of the Sub-Leading Leptons for MHP = {masses[idx]}')
        ax6.legend()
        
      
        
        ax7.hist(df["Leading_Jet_UT_PT"], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax7.hist(df_background["Leading_Jet_UT_PT"], bins=60, edgecolor='black', histtype='step', label='Background')
        ax7.set_xlabel(r'$PT_j$')
        ax7.set_ylabel('Frequency')
        ax7.set_title(f'leading Jet PT for MHP = {masses[idx]}')
        ax7.legend()
        
        ax8.hist(df["Subleading_Jet_UT_PT"], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax8.hist(df_background["Subleading_Jet_UT_PT"], bins=60, edgecolor='black', histtype='step', label='Background')
        ax8.set_xlabel(r'$PT_j$')
        ax8.set_ylabel('Frequency')
        ax8.set_title(f'Sub leading jet PT for MHP = {masses[idx]}')
        ax8.legend()
      
        '''# Plot abs cos_Theta for the leading tau tagged jet
        ax7.hist(df['LeadingTauTaggedCos_Theta'].abs(), bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax7.hist(df_background['LeadingTauTaggedCos_Theta'].abs(), bins=60, edgecolor='black', histtype='step', label='Background')
        ax7.set_xlabel(r'$|\cos(\theta_{\tau})|$')
        ax7.set_ylabel('Frequency')
        ax7.set_title(f'|Cos theta| of the Leading tau tagged jet\'  for MHP = {masses[idx]}')
        ax7.legend()
        
         # Plot ET for the leading tau tagged jet
        ax8.hist(df['LeadingTauTagged_ET'], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax8.hist(df_background['LeadingTauTagged_ET'], bins=60, edgecolor='black', histtype='step', label='Background')
        ax8.set_xlabel(r'$E_{\tau}$')
        ax8.set_ylabel('Frequency')
        ax8.set_title(f'ET of the Leading tau tagged jet\'  for MHP = {masses[idx]}')
        ax8.legend()
        
        # Plot ET for the SUBleading tau tagged jet
        ax9.hist(df['Sub_LeadingTauTagged_ET'], bins=60, edgecolor='crimson', histtype='step', label='Signal')
        if background_provided:
            ax9.hist(df_background['Sub_LeadingTauTagged_ET'], bins=60, edgecolor='black', histtype='step', label='Background')
        ax9.set_xlabel(r'$E_{\tau}$')
        ax9.set_ylabel('Frequency')
        ax9.set_title(f'ET of the Sub_Leading tau tagged jet\'  for MHP = {masses[idx]}')
        ax9.legend()'''


        # Show plots
        #plt.tight_layout()
        plt.show()


def open_plot_jets(masses, signal_paths=None, background_paths=None, signal_dfs=None, background_dfs=None):
    """
    Plots histograms of PT, Theta, and cos(Theta) for signal and background data.

    Args:
    - masses: List of mass values used in titles.
    - signal_paths: List of paths to the signal ROOT files (optional if signal_dfs is provided).
    - background_paths: List of paths to the background ROOT files (optional if background_dfs is provided).
    - signal_dfs: List of DataFrames for signal data (optional if signal_paths is provided).
    - background_dfs: List of DataFrames for background data (optional if background_paths is provided).
    """
    
    if signal_dfs is not None:
        # Use provided DataFrames for signal data
        dfs = signal_dfs
    elif signal_paths is not None:
        # Extract data for signal files
        dfs = [open_root(signal_path) for signal_path in signal_paths]
    else:
        raise ValueError("Either signal_paths or signal_dfs must be provided.")

    if background_dfs is not None:
        # Use provided DataFrames for background data
        df_background = pd.concat(background_dfs, ignore_index=True)
        background_provided = True
    elif background_paths is not None:
        # Extract data for background files and concatenate
        background_dfs = [open_root(file_path) for file_path in background_paths]
        df_background = pd.concat(background_dfs, ignore_index=True)
        background_provided = True
    else:
        background_provided = False

    # Iterate over the extracted signal DataFrames and plot histograms
    for idx, df in enumerate(dfs):
        
        fig, ax7 = plt.subplots()
        fig, ax8 = plt.subplots()
        fig, ax9 = plt.subplots()
        fig, ax11 = plt.subplots()
        fig, ax12 = plt.subplots()
        fig, ax13 = plt.subplots()
        fig, ax10 = plt.subplots()
        fig, ax2 = plt.subplots()
        fig, ax1 = plt.subplots()




        
        
        ax7.hist(df["Leading_Jet_UT_PT"], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax7.hist(df_background["Leading_Jet_UT_PT"], bins=60, edgecolor='black', histtype='step', label='Background')
        ax7.set_xlabel(r'$PT_j$')
        ax7.set_ylabel('Frequency')
        ax7.set_title(f'leading Jet PT for MHP = {masses[idx]}')
        ax7.legend()
        
        ax8.hist(df["Subleading_Jet_UT_PT"], bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax8.hist(df_background["Subleading_Jet_UT_PT"], bins=60, edgecolor='black', histtype='step', label='Background')
        ax8.set_xlabel(r'$PT_j$')
        ax8.set_ylabel('Frequency')
        ax8.set_title(f'Sub leading jet PT for MHP = {masses[idx]}')
        ax8.legend()
        
        #Plot of phi
        
#         ax12.hist(df['Leading_Jet_UT_Phi'], bins=60, edgecolor='red', histtype='step', label='Signal')
#         if background_provided:
#             ax12.hist(df_background['Leading_Jet_UT_Phi'], bins=60, edgecolor='black', histtype='step', label='Background')
#         ax12.set_xlabel(r'$\phi$')
#         ax12.set_ylabel('Frequency')
#         ax12.set_title(f'Azimuthal angle of the leading jets for MHP = {masses[idx]}')
#         ax12.legend()
        
        
        #Plot of theta
        
        ax13.hist(df['Leading_Jet_UT_Theta'], bins=60, edgecolor='crimson', histtype='step', label='Signal')
        if background_provided:
            ax13.hist(df_background['Leading_Jet_UT_Theta'], bins=60, edgecolor='black', histtype='step', label='Background')
        ax13.set_xlabel(r'$\theta$')
        ax13.set_ylabel('Frequency')
        ax13.set_title(f'Polar angle of the leading jet : MHP = {masses[idx]}')
        ax13.legend()
        
        
        
        #misssing energy
#         ax11.hist(df['MissingET_MET'], bins=60, edgecolor='red', histtype='step', label='Signal')
#         if background_provided:
#             ax11.hist(df_background['MissingET_MET'], bins=60, edgecolor='black', histtype='step', label='Background')
#         ax11.set_xlabel(r'$E_t$')
#         ax11.set_ylabel('Frequency')
#         ax11.set_title(f'Missing energy of the hadronic channel for MHP = {masses[idx]}')
#         ax11.legend()
        
        
      
        # Plot abs cos_Theta for the leading tau tagged jet
#         ax9.hist(df['LeadingTauTaggedCos_Theta'].abs(), bins=60, edgecolor='red', histtype='step', label='Signal')
#         if background_provided:
#             ax9.hist(df_background['LeadingTauTaggedCos_Theta'].abs(), bins=60, edgecolor='black', histtype='step', label='Background')
#         ax9.set_xlabel(r'$|\cos(\theta_{\tau})|$')
#         ax9.set_ylabel('Frequency')
#         ax9.set_title(f'|Cos theta| of the Leading tau tagged jet\'  for MHP = {masses[idx]}')
#         ax9.legend()
        
#          # Plot ET for the leading tau tagged jet
#         ax10.hist(df['LeadingTauTagged_ET'], bins=60, edgecolor='red', histtype='step', label='Signal')
#         if background_provided:
#             ax10.hist(df_background['LeadingTauTagged_ET'], bins=60, edgecolor='black', histtype='step', label='Background')
#         ax10.set_xlabel(r'$E_{\tau}$')
#         ax10.set_ylabel('Frequency')
#         ax10.set_title(f'ET of the Leading tau tagged jet\'  for MHP = {masses[idx]}')
#         ax10.legend()
        
#         # Plot ET for the SUBleading tau tagged jet
#         ax1.hist(df['Sub_LeadingTauTagged_ET'], bins=60, edgecolor='crimson', histtype='step', label='Signal')
#         if background_provided:
#             ax1.hist(df_background['Sub_LeadingTauTagged_ET'], bins=60, edgecolor='black', histtype='step', label='Background')
#         ax1.set_xlabel(r'$E_{\tau}$')
#         ax1.set_ylabel('Frequency')
#         ax1.set_title(f'ET of the Sub_Leading tau tagged jet\'  for MHP = {masses[idx]}')
#         ax1.legend()


        # Show plots
        #plt.tight_layout()
        plt.show()


    
def open_plot_tau(masses, signal_paths=None, background_paths=None, signal_dfs=None, background_dfs=None ,offset=[0,0,0],density=True):
    """
    Plots histograms of PT, Theta, and cos(Theta) for signal and background data.

    Args:
    - masses: List of mass values used in titles.
    - signal_paths: List of paths to the signal ROOT files (optional if signal_dfs is provided).
    - background_paths: List of paths to the background ROOT files (optional if background_dfs is provided).
    - signal_dfs: List of DataFrames for signal data (optional if signal_paths is provided).
    - background_dfs: List of DataFrames for background data (optional if background_paths is provided).
    """
    
    
    if signal_dfs is not None:
        # Use provided DataFrames for signal data
        dfs = signal_dfs
    elif signal_paths is not None:
        # Extract data for signal files
        dfs = [open_root(signal_path) for signal_path in signal_paths]
    else:
        raise ValueError("Either signal_paths or signal_dfs must be provided.")

    if background_dfs is not None:
        # Use provided DataFrames for background data
        df_background = pd.concat(background_dfs, ignore_index=True)
        background_provided = True
    elif background_paths is not None:
        # Extract data for background files and concatenate
        background_dfs = [open_root(file_path) for file_path in background_paths]
        df_background = pd.concat(background_dfs, ignore_index=True)
        background_provided = True
    else:
        background_provided = False

    # Iterate over the extracted signal DataFrames and plot histograms
    for idx, df in enumerate(dfs):
        fig, ax7 = plt.subplots()
        fig, ax8 = plt.subplots()
        fig, ax9 = plt.subplots()


        
        
        # Plot abs cos_Theta for the leading tau tagged jet
        ax7.hist(df['LeadingTauTaggedCos_Theta'].abs() ,density=density, bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax7.hist(df_background['LeadingTauTaggedCos_Theta'].abs() + offset[0] ,density=density, bins=60, edgecolor='black', histtype='step', label='Background')
        ax7.set_xlabel(r'$|\cos(\theta_{\tau})|$')
        ax7.set_ylabel('Frequency')
        ax7.set_title(f'|Cos theta| of the Leading tau tagged jet\'  for MHP = {masses[idx]}')
        ax7.legend()
        
         # Plot ET for the leading tau tagged jet
        ax8.hist(df['LeadingTauTagged_ET'] ,density=True, bins=60, edgecolor='red', histtype='step', label='Signal')
        if background_provided:
            ax8.hist(df_background['LeadingTauTagged_ET'] + offset[1] ,density=density, bins=60, edgecolor='black', histtype='step', label='Background')
        ax8.set_xlabel(r'$E_{\tau}$')
        ax8.set_ylabel('Frequency')
        ax8.set_title(f'ET of the Leading tau tagged jet\'  for MHP = {masses[idx]}')
        ax8.legend()
        
        # Plot ET for the SUBleading tau tagged jet
        ax9.hist(df['Sub_LeadingTauTagged_ET'] , density=density ,bins=60, edgecolor='crimson', histtype='step', label='Signal')
        if background_provided:
            ax9.hist(df_background['Sub_LeadingTauTagged_ET']+ offset[2],density=density, bins=60, edgecolor='black', histtype='step', label='Background')
        ax9.set_xlabel(r'$E_{\tau}$')
        ax9.set_ylabel('Frequency')
        ax9.set_title(f'ET of the Sub_Leading tau tagged jet\'  for MHP = {masses[idx]}')
        ax9.legend()

        # Show plots
        #plt.tight_layout()
        plt.show()

    


def cos_cut(df, pos_threshold, neg_threshold):
    """
    Filters the DataFrame based on the conditions for positive and negative lepton charges,
    and outputs a summary table showing the number of events before and after the cut,
    along with signal efficiency and rejection.

    Args:
        df (pd.DataFrame): The input DataFrame containing 'Combined_Charge' and 'Combined_cos_theta'.
        pos_threshold (float): The maximum value for the positive lepton's cos_theta.
        neg_threshold (float): The minimum value for the negative lepton's cos_theta.

    Returns:
        pd.DataFrame: The filtered DataFrame meeting the specified conditions.
    """
    # Initial count of input events
    input_events = len(df)
    
    # Apply the cos_cut to filter the DataFrame
    df_filtered = df[
        df.apply(
            lambda row: all(
                (charge < 0 and cos_theta > neg_threshold) or (charge > 0 and cos_theta < pos_threshold)
                for charge, cos_theta in zip(row['Combined_Charge'], row['Combined_cos_theta'])
            ),
            axis=1,
        )
    ]
    
    # Count of events after the cut
    output_events = len(df_filtered)
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"cos(l+)<{pos_threshold}&cos(l-)>{neg_threshold}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_filtered
    
    
    
    
def compute_delta_R_jets(df):
    """
    Computes the Delta R between jets based on their Phi and Eta values in the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing 'Jet_Eta' and 'Jet_Phi' columns.

    Returns:
        pd.DataFrame: The DataFrame with a new column 'Delta_R_jets' containing the computed Delta R values.
    """
    # Ensure columns are lists
    df['Jet_Eta'] = df['Jet_Eta'].apply(lambda x: x if isinstance(x, list) else [])
    df['Jet_Phi'] = df['Jet_Phi'].apply(lambda x: x if isinstance(x, list) else [])
    
    # Calculate Delta R
    df['Delta_R_jets'] = df.apply(
        lambda row: np.sqrt(
            (row['Jet_Eta'][0] - row['Jet_Eta'][1])**2 + (row['Jet_Phi'][0] - row['Jet_Phi'][1])**2
        ) if len(row['Jet_Eta']) == 2 and len(row['Jet_Phi']) == 2 else np.nan,
        axis=1
    )
    
    return df

def compute_delta_R_jets_untagged(df):
    """
    Computes the Delta R between jets based on their Phi and Eta values in the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing 'Jet_Eta' and 'Jet_Phi' columns.

    Returns:
        pd.DataFrame: The DataFrame with a new column 'Delta_R_jets' containing the computed Delta R values.
    """
    # Ensure columns are lists
    df['Jet_UT_Eta'] = df['Jet_UT_Eta'].apply(lambda x: x if isinstance(x, list) else [])
    df['Jet_UT_Phi'] = df['Jet_UT_Phi'].apply(lambda x: x if isinstance(x, list) else [])
    
    # Calculate Delta R
    df['Delta_R_jets_UT'] = df.apply(
        lambda row: np.sqrt(
            (row['Jet_UT_Eta'][0] - row['Jet_UT_Eta'][1])**2 + (row['Jet_UT_Phi'][0] - row['Jet_UT_Phi'][1])**2
        ) if len(row['Jet_UT_Eta']) == 2 and len(row['Jet_UT_Phi']) == 2 else np.nan,
        axis=1
    )
    
    return df
    
    
def Delta_R_cut_jets(df,del_R_cut):
    #count input
    input_events = len(df)
    
    df_final = df [df["Delta_R_jets"] < del_R_cut]
    
    # Count of events after the cut
    output_events = len(df_final)
    
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"Delta_R<{del_R_cut}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final
    
    
def PT_cut_jets(df, PT_cut, condition='greater'):
    """
    Filters the DataFrame based on the PT values and specified condition, and prints a summary table.

    Args:
        df (pd.DataFrame): The input DataFrame containing the 'Jet_PT' column.
        PT_cut (float): The threshold value for filtering.
        condition (str): The condition for filtering, either 'greater' or 'less'. Defaults to 'greater'.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    # Count input events
    input_events = len(df)
    
    # Apply filter based on the condition
    if condition == 'greater':
        df_final = df[df['Jet_PT'].apply(lambda pts: all(pt > PT_cut for pt in pts))]
    elif condition == 'smaller':
        df_final = df[df['Jet_PT'].apply(lambda pts: all(pt < PT_cut for pt in pts))]
    else:
        raise ValueError("Condition must be either 'greater' or 'less'")
    
    # Count of events after the cut
    output_events = len(df_final)
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"Jet_PT {'>' if condition == 'greater' else '<'} {PT_cut}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final


def PT_cut_leptons(df, PT_cut, condition='greater'):
    """
    Filters the DataFrame based on the PT values and specified condition, and prints a summary table.

    Args:
        df (pd.DataFrame): The input DataFrame containing the 'Jet_PT' column.
        PT_cut (float): The threshold value for filtering.
        condition (str): The condition for filtering, either 'greater' or 'less'. Defaults to 'greater'.

    Returns:
        pd.DataFrame: The filtered DataFrame.
    """
    # Count input events
    input_events = len(df)
    
    # Apply filter based on the condition
    if condition == 'greater':
        df_final = df[df['Combined_PT'].apply(lambda pts: all(pt > PT_cut for pt in pts))]
    elif condition == 'smaller':
        df_final = df[df['Combined_PT'].apply(lambda pts: all(pt < PT_cut for pt in pts))]
    else:
        raise ValueError("Condition must be either 'greater' or 'less'")
    
    # Count of events after the cut
    output_events = len(df_final)
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"Lepton_PT {'>' if condition == 'greater' else '<'} {PT_cut}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final
    
    
def jet_primary_cut(df):
    #count input
    input_events = len(df)
    
    
    
    df_final = df [( df['Jet Count'] >=4 ) ]
    
    df_final = df_final [(df_final["Jet_TauTagged_PT"].apply(len) >= 1) ]
    
    #df_final = df_final [( df_final['Jet_PT'].apply(len) >=5 ) ]
    
            
        
    
    # Count of events after the cut
    output_events = len(df_final)
    
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        ["Entirely Hadronic (tau inc.)", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final
    
    
def tau_jet_cut(df,count):
    #count input
    input_events = len(df)
    
    df_final = df [(df["Jet_TauTagged_PT"].apply(len) >= count) ]
    
    if count == 1 :
        cut_crit = "Hadronic Tau decay"
    elif count == 2 :
        cut_crit = "Tau di_jet"
            
        
    
    # Count of events after the cut
    output_events = len(df_final)
    
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [cut_crit, output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final



def missing_et_cut(df,criteria):
    #count input
    input_events = len(df)
    
    df_final = df[df["MissingET_MET"].apply(lambda x: any(item > criteria for item in x))]
    
    
    # Count of events after the cut
    output_events = len(df_final)
    
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"Missing_ET>{criteria}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final
    
    
    
def qcd_jet_lepton_cut(df,criteria):
    #count input
    input_events = len(df)
    
#     df_final = df_merged[df_merged["MissingET_MET"].apply(lambda x: any(item > criteria for item in x))]
    
    ##PT cut for External Jets
    #drop the rows that meet the criteria hence the ~ in front of the condition
    df_final = df[~df['Jet_UT_PT'].apply(lambda x: any(abs(item) > criteria for item in x))]
    
    
    # Count of events after the cut
    output_events = len(df_final)
    
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"|Lepton/Jet PT|>{criteria}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_final
    
    
    
def cos_cut_tau(df, criteria):
    # Initial count of input events
    input_events = len(df)
    
    # Apply the cos_cut to filter the DataFrame
    df_filtered = df[ df['LeadingTauTaggedCos_Theta'].abs() < criteria ]
    # Count of events after the cut
    output_events = len(df_filtered)
    
    # Calculate Signal Efficiency and Rejection
    signal_efficiency = output_events / input_events if input_events > 0 else 0
    rejection = 1 - signal_efficiency

    # Display the summary table
    table = [
        ["Cut", "Events", "Efficiency", "Rejection"],
        ["input events", input_events, "-", "-"],
        [f"|cos(theta_tau<{criteria}", output_events, f"{signal_efficiency:.4f}", f"{rejection:.4f}"]
    ]
    print(tabulate(table, headers="firstrow", tablefmt="grid"))

    return df_filtered
    
    


def all_tau_cuts(df,missing_et_cut_value, qcd_jet_lepton_cut_value, cos_cut_tau_value):
    # Apply the tau_jet_cut
    #Hadronic
    df_2 = tau_jet_cut(df, 1)
    #Di_tau
    df_3 = tau_jet_cut(df_2, 2)
    
    # Apply the missing_et_cut
    df_4 = missing_et_cut(df_3, missing_et_cut_value)
    
    # Apply the qcd_jet_lepton_cut
    df_5 = qcd_jet_lepton_cut(df_4, qcd_jet_lepton_cut_value)
    
    # Apply the cos_cut_tau
    df_6 = cos_cut_tau(df_5, cos_cut_tau_value)
    
    return df_2, df_3, df_4, df_5, df_6
# Usage example
#df_2, df_3, df_4, df_5, df_6 = all_tau_cuts(df,5, 5, 0.5)


# Extract the leading jet's UT_phi and UT_eta
# Extract the leading jet's UT_phi and UT_eta
def extract_leading_jet_info(row):
    pt_list = row['Jet_UT_PT']
    phi_list = row['Jet_UT_Phi']
    eta_list = row['Jet_UT_Eta']
    
    if isinstance(pt_list, list) and pt_list:
        leading_pt = max(pt_list)
        leading_index = pt_list.index(leading_pt)
        
        leading_phi = phi_list[leading_index] if leading_index < len(phi_list) else np.nan
        leading_eta = eta_list[leading_index] if leading_index < len(eta_list) else np.nan
        
        return pd.Series([leading_phi, leading_eta])
    else:
        return pd.Series([np.nan, np.nan])



def open_root_collab(file_path, tree_name="Delphes"):
    """Extracts lepton and jet data from the given ROOT file."""
    # Open the ROOT file

    
    
    file = uproot3.open(file_path)

    tree = file[tree_name]

    # To inspect the structure of 'Electron', we list its contents
    electron_data = tree.arrays(['Electron*'], namedecode="utf-8")
    # Optional: Convert the data to a DataFrame if desired
    df_electrons = pd.DataFrame({key: value for key, value in electron_data.items()}) [['Electron.PT', 'Electron.Phi', 'Electron.Eta', 'Electron.Charge']]

    # # Extract and convert Muon data
    # muon_data = tree.arrays(['Muon*'], namedecode="utf-8")
    # df_muons = pd.DataFrame({key: value for key, value in muon_data.items()}) [['Muon.PT', 'Muon.Phi', 'Muon.Eta', 'Muon.Charge']]
    
    
    # Extract and convert Particle data individually then make data frame
    
    particle_pid = tree['Particle.PID'].array()
    particle_pt = tree['Particle.PT'].array()
    particle_phi = tree['Particle.Phi'].array()
    particle_eta = tree['Particle.Eta'].array()
    particle_mass = tree['Particle.Mass'].array()
    particle_energy = tree['Particle.E'].array()
    
    #Combine into DataFrame
    df_particles = pd.DataFrame({ 'particle_pid': particle_pid ,'particle.PT': particle_pt, 'particle.Phi': particle_phi, 'particle.Eta': particle_eta, 'particle.mass': particle_mass , 'particle_energy' : particle_energy})
    
    muon_pt = tree['Muon.PT'].array()
    muon_phi = tree['Muon.Phi'].array()
    muon_eta = tree['Muon.Eta'].array()
    muon_charge = tree['Muon.Charge'].array()

    # Extract and convert Muon data individually
    muon_pt = tree['Muon.PT'].array()
    muon_phi = tree['Muon.Phi'].array()
    muon_eta = tree['Muon.Eta'].array()
    muon_charge = tree['Muon.Charge'].array()

    #Combine into DataFrame
    df_muons = pd.DataFrame({'Muon.PT': muon_pt, 'Muon.Phi': muon_phi, 'Muon.Eta': muon_eta, 'Muon.Charge': muon_charge})



    # # Extract and convert Missing ET data
    missing_et_data = tree.arrays(['MissingET*'], namedecode="utf-8")
    df_missing_et = pd.DataFrame({key: value for key, value in missing_et_data.items()})[['MissingET.MET', 'MissingET.Phi', 'MissingET.Eta']]

    # # Extract and convert Jet data
    jet_pt_data = tree['Jet.PT'].array()
    jet_eta_data = tree['Jet.Eta'].array()
    jet_tau_tagged_data = tree['Jet.TauTag'].array()
    jet_mass_data = tree['Jet.Mass'].array()
    jet_phi_data = tree['Jet.Phi'].array()

    # Convert to DataFrame
    df_jets = pd.DataFrame({
    'Jet_PT': jet_pt_data,
    'Jet_Eta': jet_eta_data,
    'Jet_TauTagged': jet_tau_tagged_data,
    'Jet_Mass': jet_mass_data,
    'Jet_Phi': jet_phi_data})
  

    # Rename columns for the Electron DataFrame
    df_electrons.columns = ['Electron_PT', 'Electron_Phi', 'Electron_Eta', 'Electron_Charge']
    df_muons.columns = ['Muon_PT', 'Muon_Phi', 'Muon_Eta', 'Muon_Charge']
    df_missing_et.columns = ['MissingET_MET', 'MissingET_Phi', 'MissingET_Eta']



    # Merge the DataFrames on the 'Event' column
    df_merged = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True, how='outer'), [df_particles,df_electrons, df_muons, df_jets, df_missing_et])

    # Convert array-like elements to lists using map for each column
    df_merged = df_merged.apply(lambda col: col.map(lambda x: list(x) if isinstance(x, (np.ndarray, list)) else x))

    # Convert string representations of lists back to actual lists
    for col in df_merged.columns:
        df_merged[col] = df_merged[col].apply(
            lambda x: ast.literal_eval(x) if isinstance(x, str) and x.startswith('[') and x.endswith(']') else x
        )
    
    # Now, convert the elements within the lists to floats
    for col in df_merged.columns:
        df_merged[col] = df_merged[col].apply(
            lambda x: list(map(float, x)) if isinstance(x, list) else x
        )
        


    #calculate the theta of the charged scalar "
    df_merged['theta_Higgs'] = df_merged['particle.Eta'].apply(lambda eta_list: 2 * np.arctan(np.exp(-eta_list[4])) if len(eta_list) > 4 else np.nan)
    
    # Remove the 'theta_Higgs' column temporarily
    theta_higgs = df_merged.pop('theta_Higgs')

    # Insert 'theta_Higgs' as the first column
    df_merged.insert(4, 'theta_Higgs', theta_higgs)
    
    #Fix MET, it was made of lists
    df_merged['MissingET_MET'] = df_merged['MissingET_MET'].apply(
    lambda x: x[0] if isinstance(x, list) and x else np.nan
)


    
    # Combine PTs of electrons and muons, handling cases where lists may be empty
    df_merged['Combined_PT'] = df_merged['Electron_PT'].apply(lambda x: x if x else []) + df_merged['Muon_PT'].apply(lambda x: x if x else [])
    df_merged['Combined_Charge'] = df_merged['Electron_Charge'].apply(lambda x: x if x else []) + df_merged['Muon_Charge'].apply(lambda x: x if x else [])

    # Calculate Lepton_Theta as the combined theta from electrons and muons
    # df_merged['Electron_Theta'] = df_merged['Electron_Eta'].apply(lambda x: 2 * np.arctan(np.exp(-x)) if x else [])
    df_merged['Electron_Theta'] = df_merged['Electron_Eta'].apply(
    lambda etas: [2 * np.arctan(np.exp(-eta)) for eta in etas] if etas else [])
    df_merged['Electron_cos_theta'] = df_merged['Electron_Theta'].apply(lambda theta: np.cos(theta) if theta is not None else None)
    df_merged['Muon_Theta']= df_merged['Muon_Eta'].apply(lambda etas: [ 2 * np.arctan(np.exp(-eta)) for eta in etas])
    # df_merged['Muon_Theta'] = df_merged['Muon_Eta'].apply(lambda x: 2 * np.arctan(np.exp(-x)) if x else [])
    df_merged['Muon_cos_theta'] = df_merged['Muon_Theta'].apply(lambda theta: np.cos(theta) if theta is not None else None)
    df_merged['Combined_Theta'] = df_merged['Electron_Theta'].apply(lambda x: x if x else []) + df_merged['Muon_Theta'].apply(lambda x: x if x else [])
    
    df_merged['Combined_cos_theta'] = df_merged['Electron_cos_theta'].apply(list) + df_merged['Muon_cos_theta'].apply(list)

    df_merged['Combined_PT_Theta'] = df_merged.apply(
        lambda row: list(zip(row['Combined_PT'], row['Combined_Theta'])) if row['Combined_PT'] and row['Combined_Theta'] else [], axis=1
    )
    
    df_merged['Combined_PT_Charge'] = df_merged.apply(
        lambda row: list(zip(row['Combined_PT'], row['Combined_Charge'])) if row['Combined_PT'] and row['Combined_Charge'] else [], axis=1
    )
    
    # Sort by PT in descending order, keeping the corresponding Theta and Charge values aligned
    df_merged['Sorted_PT_Theta'] = df_merged['Combined_PT_Theta'].apply(
        lambda pairs: sorted(pairs, key=lambda x: x[0], reverse=True)
    )
    
    df_merged['Sorted_PT_Charge'] = df_merged['Combined_PT_Charge'].apply(
        lambda pairs: sorted(pairs, key=lambda x: x[0], reverse=True)
    )
    
    # Extract leading and subleading Max_PT, Max_Theta, and Charge based on the sorted order
    df_merged['Leading_lep_PT'] = df_merged['Sorted_PT_Theta'].apply(lambda pairs: pairs[0][0] if pairs else None)
    df_merged['Leading_lep_Theta'] = df_merged['Sorted_PT_Theta'].apply(lambda pairs: pairs[0][1] if pairs else None)
    df_merged['Leading_lep_Charge'] = df_merged['Sorted_PT_Charge'].apply(lambda pairs: pairs[0][1] if pairs else None)
    
    df_merged['Sub_leading_lep_PT'] = df_merged['Sorted_PT_Theta'].apply(
        lambda pairs: pairs[1][0] if len(pairs) > 1 else None
    )
    df_merged['Sub_leading_lep_Theta'] = df_merged['Sorted_PT_Theta'].apply(
        lambda pairs: pairs[1][1] if len(pairs) > 1 else None
    )
    df_merged['Sub_leading_lep_Charge'] = df_merged['Sorted_PT_Charge'].apply(
        lambda pairs: pairs[1][1] if len(pairs) > 1 else None
    )
    
    ###Separating the data by taggings
    
    
    def separate_tagged_untagged_jets(df):
      # Extract the necessary columns
        jet_pts = df['Jet_PT']
        jet_eta = df['Jet_Eta']
        jet_phi = df['Jet_Phi']
        jet_tagged = df['Jet_TauTagged']
        # Initialize lists to store the separated jet data
        tagged_pts = []
        tagged_eta = []
        tagged_phi = []
        untagged_pts = []
        untagged_eta = []
        untagged_phi = []

        # Separate the data based on the tagging
        for pts, eta, phi, tagged in zip(jet_pts, jet_eta, jet_phi, jet_tagged):
            if tagged == 1:  # Check for tagged (1)
                tagged_pts.append(pts)
                tagged_eta.append(eta)
                tagged_phi.append(phi)
            else:  # Check for untagged (0)
                untagged_pts.append(pts)
                untagged_eta.append(eta)
                untagged_phi.append(phi)
        # Update DataFrame with separated data
        df['Jet_TauTagged_PT'] = [
            [pt for pt, tag in zip(pts, tags) if tag == 1]
            for pts, tags in zip(jet_pts, jet_tagged)
        ]
        df['Jet_TauTagged_Eta'] = [
            [eta for eta, tag in zip(etas, tags) if tag == 1]
            for etas, tags in zip(jet_eta, jet_tagged)
        ]
        df['Jet_TauTagged_Phi'] = [
            [phi for phi, tag in zip(phis, tags) if tag == 1]
            for phis, tags in zip(jet_phi, jet_tagged)
        ]
        df['Jet_UT_PT'] = [
            [pt for pt, tag in zip(pts, tags) if tag == 0]
            for pts, tags in zip(jet_pts, jet_tagged)
        ]
        df['Jet_UT_Eta'] = [
            [eta for eta, tag in zip(etas, tags) if tag == 0]
            for etas, tags in zip(jet_eta, jet_tagged)
        ]
        df['Jet_UT_Phi'] = [
            [phi for phi, tag in zip(phis, tags) if tag == 0]
            for phis, tags in zip(jet_phi, jet_tagged)
        ]
        return df
    

      

    # Apply the function to the DataFrame
    
    df_merged = separate_tagged_untagged_jets(df_merged)
    
    # Calculate the leading (maximum) PT
    df_merged["Leading_Jet_UT_PT"] = df_merged["Jet_UT_PT"].apply(lambda x: max(x) if isinstance(x, list) and x else np.nan)

    # Calculate the subleading (second largest) PT directly
    df_merged["Subleading_Jet_UT_PT"] = df_merged["Jet_UT_PT"].apply(lambda x: sorted(x, reverse=True)[1] if isinstance(x, list) and len(x) > 1 else np.nan)
    
    
    ##apply function
    df_merged[['Leading_Jet_UT_Phi', 'Leading_Jet_UT_Eta']] = df_merged.apply(extract_leading_jet_info, axis=1)
    
    df_merged['Leading_Jet_UT_Theta'] = df_merged['Leading_Jet_UT_Eta'].apply(
    lambda eta: 2 * np.arctan(np.exp(-eta)) if not pd.isna(eta) else np.nan
)

    
     

    
    
    
    
    ###Masses of the Tau Tagged
    # Zip Jet_Mass and Jet_TauTagged columns
    df_merged['Jet_TauTagged_Mass'] = df_merged.apply(lambda row: [mass for mass, tag in zip(row['Jet_Mass'], row['Jet_TauTagged']) if tag == 1], axis=1)
    df_merged['Jet_UT_Mass'] = df_merged.apply(lambda row: [mass for mass, tag in zip(row['Jet_Mass'], row['Jet_TauTagged']) if tag == 0], axis=1)
    
    # Calculate Jet_TauTagged_ET using the formula
    df_merged['Jet_TauTagged_ET'] = df_merged.apply(lambda row: [(mass**2 + pt**2)**0.5 for mass, pt in zip(row['Jet_TauTagged_Mass'], row['Jet_TauTagged_PT'])], axis=1) 


    
    ### Cos theta of the leading tau tagged jet
    
    df_merged['TauTaggedTheta'] = df_merged['Jet_TauTagged_Eta'].apply(lambda etas: [ 2 * np.arctan(np.exp(-eta)) for eta in etas])
    
    df_merged['TauTaggedCos_Theta'] = df_merged['TauTaggedTheta'].apply(lambda thetas: [np.cos(theta) for theta in thetas])
    
    # Combine TauTaggedPT with TauTaggedCos_Theta using zip to pair PT with cos(theta)
    df_merged['TauTagged_PT_CosTheta'] = df_merged.apply(lambda row: list(zip(row['Jet_TauTagged_PT'], row['TauTaggedCos_Theta'])), axis=1)
        
    # Find the Leading Tau Tagged PT
    df_merged['Leading_TauTaggedPT'] = df_merged['Jet_TauTagged_PT'].apply(lambda pts: max(pts) if pts else None)
    
    # Find the corresponding Leading Cos Theta of the Leading Tau Tagged Jet
    df_merged['LeadingTauTaggedCos_Theta'] = df_merged['TauTagged_PT_CosTheta'].apply(lambda pairs: max(pairs, key=lambda x: x[0])[1] if pairs else None) 
    
    #Combining the PT'S and ET'S
    df_merged['TauTagged_PT_ET'] = df_merged.apply(lambda row: list(zip(row['Jet_TauTagged_PT'], row['Jet_TauTagged_ET'])), axis=1)
    
    #Now finding the et of the leading tau jet
    df_merged['LeadingTauTagged_ET'] = df_merged['TauTagged_PT_ET'].apply(lambda pairs: max(pairs, key=lambda x: x[0])[1] if pairs else None) 
    
    df_merged['Sub_LeadingTauTagged_ET'] = df_merged['TauTagged_PT_ET'].apply(
    lambda pairs: min(pairs, key=lambda x: x[0])[1] if len(pairs) > 1 else None
)

    
    
    

    
    #### Counting leptons and jets
    
    count_lep = [len(k) for k in df_merged["Combined_PT"]]
    count_jet = [len(k) for k in df_merged["Jet_PT"]]

    df_merged["Lepton Count"] = count_lep
    df_merged["Jet Count"] = count_jet
    
    
    
    #Drop the columns that are no longer needed
    drop_cols = [
        "Sorted_PT_Theta", "Combined_PT_Theta", "Combined_Theta",
        "Electron_Eta", "Muon_Eta",
        "Sorted_PT_Charge",'TauTagged_PT_CosTheta','TauTagged_PT_ET'
    ]
    df_merged = df_merged.drop(columns=drop_cols)

 
    

    # return df_merged 
    return df_merged
    
    
    
def open_lhe(file_path):

    # Initialize an empty list to store particle data for each event
    event_data = []

    # Open and parse the LHE file
    events = pylhe.read_lhe(file_path)

    # Initialize event counter
    event_counter = 0

    # Loop over the events and extract particle information
    for event in events:
        # Initialize lists to store kinematic data for the particles in this event
        particle_ids = []
        px_list = []
        py_list = []
        pz_list = []
        energy_list = []
        mass_list = []
        status_list = []
        mother1_list = []
        mother2_list = []

        # Loop through all particles in the event and collect their properties
        for particle in event.particles:
            particle_ids.append(particle.id)
            px_list.append(particle.px)
            py_list.append(particle.py)
            pz_list.append(particle.pz)
            energy_list.append(particle.e)
            mass_list.append(particle.m)
            status_list.append(particle.status)
            mother1_list.append(particle.mother1)
            mother2_list.append(particle.mother2)
            
        #the full momenta
        p_list = [np.sqrt(px**2 + py**2 + pz**2) for px, py, pz in zip(px_list, py_list, pz_list)]
        
        # Append the lists to the event_data list, where each row corresponds to an event
        event_data.append({
            'event_id': event_counter,  # Track the event number
            'particle_ids': particle_ids[2:],  # List of particle IDs for this event
            'px': px_list[2:],  # List of px for each particle in the event
            'py': py_list[2:],  # List of py for each particle in the event
            'pz': pz_list[2:],  # List of pz for each particle in the event
            'p':p_list[2:],
            'energy': energy_list[2:],  # List of energy for each particle in the event
            'mass': mass_list[2:],  # List of mass for each particle in the event
            'status': status_list[2:],  # List of statuses for each particle
            'mother1': mother1_list[2:],  # List of mother1 indices for each particle
            'mother2': mother2_list[2:]   # List of mother2 indices for each particle
        })

        # Increment the event counter after each event
        event_counter += 1

    # Convert the list of dictionaries into a Pandas DataFrame
    df = pd.DataFrame(event_data)
    df['pt'] = df.apply(lambda row: [np.sqrt(row['px'][i]**2 + row['py'][i]**2) for i in range(len(row['px']))], axis=1)
    df['theta'] = df.apply(lambda row: [np.arctan2(np.sqrt(px**2 + py**2), pz) for px, py, pz in zip(row['px'], row['py'], row['pz'])], axis=1)

    
    # Define the particle IDs for which we want to extract kinematics
    particles_of_interest = {
        'w': [24.0,-24],
        'z': [23.0,-23],
        'h+': 990002.0,
        'h-': -990002.0,
        'tau': [15.0,-15],
        'v_tau': [16.0,-16],
        'electron' : [-11,11],
        'muon' : [13,-13],
        'v_mu_e' : [12,-12,14,-14],
        'leptons': [-13,13,-11,11],
        'jets': [1, -1, 2, -2, 3.0, -3.0, 4.0, -4.0]
        }
    
    

    # Helper function to extract the first match's kinematics or return NaN if not found
    def get_first_kinematics(pids, kinematics, target_pids):
        # Check if target_pids is a list of possible particle IDs
        if isinstance(target_pids, list):
            # Iterate over all target_pids and find the first match in particle_ids
            for target_pid in target_pids:
                matching_kinematics = [kinematics[i] for i in range(len(pids)) if pids[i] == target_pid]
                if matching_kinematics:  # If we found a match
                    return matching_kinematics[0]  # Return the first match
        else:
            # If only one target particle ID, just use it directly
            matching_kinematics = [kinematics[i] for i in range(len(pids)) if pids[i] == target_pids]
            if matching_kinematics:
                return matching_kinematics[0]

        # If no match is found, return NaN
        return np.nan
    
    particles = ['w', 'z', 'h+', 'h-', 'tau', 'v_tau', 'electron', 'muon','v_mu_e']

    # Loop over particles and create new columns for each particle's kinematic properties
    for particle in particles:
        df[f'p_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['p'], particles_of_interest[particle]), axis=1)
        df[f'energy_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['energy'], particles_of_interest[particle]), axis=1)
        df[f'theta_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['theta'], particles_of_interest[particle]), axis=1)
        df[f'mass_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['mass'], particles_of_interest[particle]), axis=1)
        df[f'pt_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['pt'], particles_of_interest[particle]), axis=1)    
        
    #Loop over the tau and tau neutrino to get px,py,pz
    
    for particle in ['tau', 'v_tau'] :
        df[f'px_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['px'], particles_of_interest[particle]), axis=1)
        df[f'py_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['py'], particles_of_interest[particle]), axis=1)
        df[f'pz_{particle}'] = df.apply(lambda row: get_first_kinematics(row['particle_ids'], row['pz'], particles_of_interest[particle]), axis=1)
        

        #Loop over the jets and leptons
        
    particles2 = ['leptons','jets'] #jets/leptons
    
    for particle in particles2:
        df[f'energy_{particle}'] = df.apply(lambda row: [row['energy'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'mass_{particle}'] = df.apply(lambda row: [row['mass'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'theta_{particle}'] = df.apply(lambda row: [row['theta'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'p_{particle}'] = df.apply(lambda row: [row['p'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'px_{particle}'] = df.apply(lambda row: [row['px'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'py_{particle}'] = df.apply(lambda row: [row['py'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'pz_{particle}'] = df.apply(lambda row: [row['pz'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        df[f'pt_{particle}'] = df.apply(lambda row: [row['pt'][i] for i in range(len(row['particle_ids'])) if row['particle_ids'][i] in particles_of_interest[particle]], axis=1)
        
    
    #the muon and electron then its neutrinos (the column 'leptons exclude the tau)
    
    

    # For jets, we create a list of pz and energy because there can be multiple jets per event
    df["Leading_Jet_UT_PT"] = df['pt_jets'].apply(lambda x: max(x) if x else np.nan)


    # Function to find the corresponding theta for the lepton/jet with the highest pt
    def get_leading_jet_theta(pt_jets, theta_jets):
        if pt_jets and theta_jets:  # Ensure both lists are not empty
            # Zip the pt and theta values together
            pt_theta_pairs = list(zip(pt_jets, theta_jets))

            # Find the jet with the highest pt and its corresponding theta
            leading_jet = max(pt_theta_pairs, key=lambda x: x[0])  # x[0] is the pt, x[1] is theta

            # Return the corresponding theta
            return leading_jet[1]  # Return the theta corresponding to the highest pt
        return np.nan  # Return NaN if pt_jets or theta_jets is empty

    # Apply the function to each row to get the leading jet theta (based on highest pt)
    df["Leading_Jet_UT_Theta"] = df.apply(lambda row: get_leading_jet_theta(row['pt_jets'], row['theta_jets']), axis=1)
    df["Leading_Jet_UT_Cos_Theta"] = np.abs(np.cos(np.array(list(df["Leading_Jet_UT_Theta"]))))
   

    df["Subleading_Jet_UT_PT"] = df['pt_jets'].apply(lambda x: sorted(x, reverse=True)[1] if len(x) > 1 else np.nan)
    
    ###leading leptons :
    
    df["Leading_Lep_PT"] = df['pt_leptons'].apply(lambda x: max(x) if x else np.nan)
    df["Leading_Lep_energy"] = df.apply(lambda row: get_leading_jet_theta(row['pt_leptons'], row['energy_leptons']), axis=1)


    #inavariant masses for the jets
    
    # Handle empty lists by returning 0 for sum if the list is empty
    df['energy_total'] = df['energy_jets'].apply(lambda x: sum(x) if x else 0)
    df['px_total'] = df['px_jets'].apply(lambda x: sum(x) if x else 0)
    df['py_total'] = df['py_jets'].apply(lambda x: sum(x) if x else 0)
    df['pz_total'] = df['pz_jets'].apply(lambda x: sum(x) if x else 0)

    # Calculate the invariant mass using the formula: M = (E^2 - (px^2 + py^2 + pz^2)) ^0.5
    df['jet_invariant_mass'] = df.apply(
        lambda row: np.sqrt(row['energy_total']**2 - (row['px_total']**2 + row['py_total']**2 + row['pz_total']**2)) 
        if row['energy_total'] > 0 else np.nan, axis=1)
    
    def calculate_invariant_mass(row):
        total_vector = TLorentzVector()

        # Iterate over each jet in the list
        for px, py, pz, energy in zip(row['px_jets'], row['py_jets'], row['pz_jets'], row['energy_jets']):
            jet_vector = TLorentzVector()
            jet_vector.SetPxPyPzE(px, py, pz, energy)
            total_vector += jet_vector

        # Return the invariant mass of the system
        return total_vector.M()

    # Apply the function to each row and create a new column 'M_jets' with the result
    df['M_jets'] = df.apply(calculate_invariant_mass, axis=1)



    # Sum of the total energy of the tau and neutrino
    df['energy_total_tau_nu'] = df['energy_tau'] + df['energy_v_tau'] 

    # Sum of the momentum components
    df['px_total_tau_nu'] = df['px_tau'] + df['px_v_tau']
    df['py_total_tau_nu'] = df['py_tau'] + df['py_v_tau']
    df['pz_total_tau_nu'] = df['pz_tau'] + df['pz_v_tau']

    # Calculate the invariant mass using the formula: M = ( E^2 - (px^2 + py^2 + pz^2) ^0.5
    df['tau_nu_invariant_mass'] = np.sqrt(df['energy_total_tau_nu']**2 - (df['px_total_tau_nu']**2 + df['py_total_tau_nu']**2 + df['pz_total_tau_nu']**2))
    
    ###### Unsure stuff ##tests
    
    df['tau_invariant_mass'] = np.sqrt(df['energy_tau']**2 - (df['px_tau']**2 + df['py_tau']**2 + df['pz_tau']**2))
    
#     df['nu_invariant_mass'] = np.sqrt(df['energy_v_tau']**2 - (df['px_total_nu']**2 + df['py_total_nu']**2 + df['pz_total_nu']**2))
    
#px_jets','py_jets','pz_jets',

    df["Jet Count"] = df["p_jets"].apply(len)
    df["Lepton Count"] = df["p_leptons"].apply(len)
    
    drop_cols = ['event_id','px', 'py','status', 'mother1', 'mother2','px_total','py_total',
                 'pz_total','energy_total',
                 'px_leptons','py_leptons','pz_leptons'
                ]
    df = df.drop(columns=drop_cols)



    return df



def open_root_original(file_path, tree_name="Delphes"):
    """Extracts lepton and jet data from the given ROOT file."""
    # Open the ROOT file
    warnings.filterwarnings("ignore")
    file = ROOT.TFile.Open(file_path)
    
    # Check if the file is open
    if not file or file.IsZombie():
        raise Exception(f"Failed to open file {file_path}.")
    
    # Access the tree
    tree = file.Get(tree_name)
    if not tree:
        raise Exception(f"Failed to get the tree {tree_name} from file {file_path}.")
        
    #Reader for Particles
    reader_particles = ROOT.TTreeReader(tree)
    particles_pid_reader = ROOT.TTreeReaderArray('int')(reader_particles, "Particle.PID") #maybe int
    particles_mass_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.Mass") # Adjust the type if necessary
    particles_energy_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.E")  # E is often used for total energy
    particles_eta_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.Eta")
    particles_pt_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.PT")
    particles_phi_reader = ROOT.TTreeReaderArray('float')(reader_particles, "Particle.Phi")

    # Reader for electrons
    reader_electron = ROOT.TTreeReader(tree)
    electron_pt_reader = ROOT.TTreeReaderArray('float')(reader_electron, "Electron.PT")
    electron_phi_reader = ROOT.TTreeReaderArray('float')(reader_electron, "Electron.Phi")
    electron_eta_reader = ROOT.TTreeReaderArray('float')(reader_electron, "Electron.Eta")
    electron_charge_reader = ROOT.TTreeReaderArray('int')(reader_electron, "Electron.Charge")

    # Reader for muons
    reader_muon = ROOT.TTreeReader(tree)
    muon_pt_reader = ROOT.TTreeReaderArray('float')(reader_muon, "Muon.PT")
    muon_phi_reader = ROOT.TTreeReaderArray('float')(reader_muon, "Muon.Phi")
    muon_eta_reader = ROOT.TTreeReaderArray('float')(reader_muon, "Muon.Eta")
    muon_charge_reader = ROOT.TTreeReaderArray('int')(reader_muon, "Muon.Charge")

    # Reader for jets
    reader_jet = ROOT.TTreeReader(tree)
    jet_pt_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.PT")
    jet_phi_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.Phi")
    jet_eta_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.Eta")
    jet_mass_reader = ROOT.TTreeReaderArray('float')(reader_jet, "Jet.Mass")
    jet_tau_tagged_reader = ROOT.TTreeReaderArray('int')(reader_jet, "Jet.TauTag")
    jet_b_tagged_reader = ROOT.TTreeReaderArray('int')(reader_jet, "Jet.BTag")
    
    # Reader for Missing ET
    missing_et_reader = ROOT.TTreeReader(tree)
    missing_et_met_reader = ROOT.TTreeReaderArray('float')(missing_et_reader, "MissingET.MET")
    missing_et_phi_reader = ROOT.TTreeReaderArray('float')(missing_et_reader, "MissingET.Phi")
    missing_et_eta_reader = ROOT.TTreeReaderArray('float')(missing_et_reader, "MissingET.Eta")


    # Initialize lists to store data
    particles_data =[]
    electron_data = []
    muon_data = []
    jet_data = []
    missing_et_data = []
    
# Loop over events for particles
    event_number = 0
    while reader_particles.Next():
        # Read data from the tree
        particles_pid = [pid for pid in particles_pid_reader]
        particles_mass = [mass for mass in particles_mass_reader]
        particles_energy = [energy for energy in particles_energy_reader]
        particles_eta = [eta for eta in particles_eta_reader]
        particles_pt = [pt for pt in particles_pt_reader]
        particles_phi = [phi for phi in particles_phi_reader]
        
        # Append the data for each event
        particles_data.append({
            'Event': event_number,
            'Particle_PID': particles_pid,
            'Particle_Mass': particles_mass,
            'Particle_E': particles_energy,
            'Particle_Eta': particles_eta,
            'Particle_PT': particles_pt,
            'Particle_Phi': particles_phi
        })
        # Increment the event counter
        event_number += 1
    
    # Loop over events for electrons
    event_number = 0
    while reader_electron.Next():
        electron_pts = [pt for pt in electron_pt_reader]
        electron_phis = [phi for phi in electron_phi_reader]
        electron_etas = [eta for eta in electron_eta_reader]
        electron_charges = [charge for charge in electron_charge_reader]
        electron_thetas = [2 * np.arctan(np.exp(-eta)) for eta in electron_etas]
        electron_data.append({
            'Event': event_number,
            'Electron_PT': electron_pts,
            'Electron_Phi': electron_phis,
            'Electron_Eta': electron_etas,
            'Electron_Charge': electron_charges,
            'Electron_Theta': electron_thetas
        })
        event_number += 1

    # Loop over events for muons
    event_number = 0
    while reader_muon.Next():
        muon_pts = [pt for pt in muon_pt_reader]
        muon_phis = [phi for phi in muon_phi_reader]
        muon_etas = [eta for eta in muon_eta_reader]
        muon_charges = [charge for charge in muon_charge_reader]
        muon_thetas = [2 * np.arctan(np.exp(-eta)) for eta in muon_etas]
        muon_data.append({
            'Event': event_number,
            'Muon_PT': muon_pts,
            'Muon_Phi': muon_phis,
            'Muon_Eta': muon_etas,
            'Muon_Charge': muon_charges,
            'Muon_Theta': muon_thetas
        })
        event_number += 1

    # Loop over events for jets
    event_number = 0
    while reader_jet.Next():
        jet_pts = [pt for pt in jet_pt_reader]
        jet_phis = [phi for phi in jet_phi_reader]
        jet_etas = [eta for eta in jet_eta_reader]
        jet_masses = [mass for mass in jet_mass_reader]
        jet_tau_tagged = [tag for tag in jet_tau_tagged_reader]
        jet_b_tagged = [b for b in jet_b_tagged_reader]
        jet_data.append({
            'Event': event_number,
            'Jet_PT': jet_pts,
            'Jet_Phi': jet_phis,
            'Jet_Eta': jet_etas,
            'Jet_Mass': jet_masses,
            'Jet_TauTagged': jet_tau_tagged,
            'Jet_BTagged': jet_b_tagged
        })
        event_number += 1
        
    # Loop over events for Missing ET
    event_number = 0
    while missing_et_reader.Next():
        missing_met_values = [met for met in missing_et_met_reader]
        missing_phi_values = [phi for phi in missing_et_phi_reader]
        missing_eta_values = [eta for eta in missing_et_eta_reader]
        missing_et_data.append({
            'Event': event_number,
            'MissingET_MET': missing_met_values,
            'MissingET_Phi': missing_phi_values,
            'MissingET_Eta': missing_eta_values
        })
        event_number += 1

    # Convert to DataFrames
    df_particles = pd.DataFrame(particles_data)
    df_electrons = pd.DataFrame(electron_data)
    df_muons = pd.DataFrame(muon_data)
    df_jets = pd.DataFrame(jet_data)
    df_missing_et = pd.DataFrame(missing_et_data)


    # Merge the DataFrames on the 'Event' column
    df_merged = reduce(lambda left, right: pd.merge(left, right, on='Event', how='outer'), [df_particles,df_electrons, df_muons, df_jets,df_missing_et])

    # Combine PTs of electrons and muons, handling cases where lists may be empty
    df_merged['Combined_PT'] = df_merged['Electron_PT'].apply(lambda x: x if x else []) + df_merged['Muon_PT'].apply(lambda x: x if x else [])
    df_merged['Combined_Charge'] = df_merged['Electron_Charge'].apply(lambda x: x if x else []) + df_merged['Muon_Charge'].apply(lambda x: x if x else [])

    # Calculate Lepton_Theta as the combined theta from electrons and muons
    df_merged['Electron_Theta'] = df_merged['Electron_Theta'].apply(lambda x: x if x else [])
    df_merged['Electron_cos_theta'] = df_merged['Electron_Theta'].apply(lambda theta: np.cos(theta) if theta is not None else None)
    df_merged['Muon_Theta'] = df_merged['Muon_Theta'].apply(lambda x: x if x else [])
    df_merged['Muon_cos_theta'] = df_merged['Muon_Theta'].apply(lambda theta: np.cos(theta) if theta is not None else None)
    df_merged['Combined_Theta'] = df_merged['Electron_Theta'].apply(lambda x: x if x else []) + df_merged['Muon_Theta'].apply(lambda x: x if x else [])
    
    #Fix MET, it was made of lists
    df_merged['MissingET_MET'] = df_merged['MissingET_MET'].apply(
    lambda x: x[0] if isinstance(x, list) and x else np.nan
)
    
    df_merged['Combined_cos_theta'] = df_merged['Electron_cos_theta'].apply(list) + df_merged['Muon_cos_theta'].apply(list)

    df_merged['Combined_PT_Theta'] = df_merged.apply(
        lambda row: list(zip(row['Combined_PT'], row['Combined_Theta'])) if row['Combined_PT'] and row['Combined_Theta'] else [], axis=1
    )
    
    df_merged['Combined_PT_Charge'] = df_merged.apply(
        lambda row: list(zip(row['Combined_PT'], row['Combined_Charge'])) if row['Combined_PT'] and row['Combined_Charge'] else [], axis=1
    )
    
    # Sort by PT in descending order, keeping the corresponding Theta and Charge values aligned
    df_merged['Sorted_PT_Theta'] = df_merged['Combined_PT_Theta'].apply(
        lambda pairs: sorted(pairs, key=lambda x: x[0], reverse=True)
    )
    
    df_merged['Sorted_PT_Charge'] = df_merged['Combined_PT_Charge'].apply(
        lambda pairs: sorted(pairs, key=lambda x: x[0], reverse=True)
    )
    
    # Extract leading and subleading Max_PT, Max_Theta, and Charge based on the sorted order
    df_merged['Leading_lep_PT'] = df_merged['Sorted_PT_Theta'].apply(lambda pairs: pairs[0][0] if pairs else None)
    df_merged['Leading_lep_Theta'] = df_merged['Sorted_PT_Theta'].apply(lambda pairs: pairs[0][1] if pairs else None)
    df_merged['Leading_lep_Charge'] = df_merged['Sorted_PT_Charge'].apply(lambda pairs: pairs[0][1] if pairs else None)
    
    df_merged['Sub_leading_lep_PT'] = df_merged['Sorted_PT_Theta'].apply(
        lambda pairs: pairs[1][0] if len(pairs) > 1 else None
    )
    df_merged['Sub_leading_lep_Theta'] = df_merged['Sorted_PT_Theta'].apply(
        lambda pairs: pairs[1][1] if len(pairs) > 1 else None
    )
    df_merged['Sub_leading_lep_Charge'] = df_merged['Sorted_PT_Charge'].apply(
        lambda pairs: pairs[1][1] if len(pairs) > 1 else None
    )
    
    ###Separating the data by taggings
    def separate_tagged_untagged_jets(df):
        
    # Extract the necessary columns
        jet_pts = df['Jet_PT']
        jet_eta = df['Jet_Eta']
        jet_phi = df['Jet_Phi']
        jet_tagged = df['Jet_TauTagged']
        
        # Initialize lists to store the separated jet data
        tagged_pts = []
        tagged_eta = []
        tagged_phi = []
        
        untagged_pts = []
        untagged_eta = []
        untagged_phi = []
        
        # Separate the data based on the tagging
        for pts, eta, phi, tagged in zip(jet_pts, jet_eta, jet_phi, jet_tagged):
            if tagged == 1:  # Check for tagged (1)
                tagged_pts.append(pts)
                tagged_eta.append(eta)
                tagged_phi.append(phi)
            else:  # Check for untagged (0)
                untagged_pts.append(pts)
                untagged_eta.append(eta)
                untagged_phi.append(phi)
        
        # Update DataFrame with separated data
        df['Jet_TauTagged_PT'] = [
            [pt for pt, tag in zip(pts, tags) if tag == 1]
            for pts, tags in zip(jet_pts, jet_tagged)
        ]
        
        df['Jet_TauTagged_Eta'] = [
            [eta for eta, tag in zip(etas, tags) if tag == 1]
            for etas, tags in zip(jet_eta, jet_tagged)
        ]
        
        df['Jet_TauTagged_Phi'] = [
            [phi for phi, tag in zip(phis, tags) if tag == 1]
            for phis, tags in zip(jet_phi, jet_tagged)
        ]
        
        df['Jet_UT_PT'] = [
            [pt for pt, tag in zip(pts, tags) if tag == 0]
            for pts, tags in zip(jet_pts, jet_tagged)
        ]
        
        df['Jet_UT_Eta'] = [
            [eta for eta, tag in zip(etas, tags) if tag == 0]
            for etas, tags in zip(jet_eta, jet_tagged)
        ]
        
        df['Jet_UT_Phi'] = [
            [phi for phi, tag in zip(phis, tags) if tag == 0]
            for phis, tags in zip(jet_phi, jet_tagged)
        ]
        
        return df

    # Apply the function to the DataFrame
    
    df_merged = separate_tagged_untagged_jets(df_merged)
    
    ## Leading Jet_UT_PT :
    
    #Is it hadronic or leptonic tau :
    
    df_merged['decay_mode'] = df_merged['Jet_TauTagged_PT'].apply(lambda x: 1 if len(x) > 0 else 0)
    
    # Calculate the leading (maximum) PT
    df_merged["Leading_Jet_UT_PT"] = df_merged["Jet_UT_PT"].apply(lambda x: max(x) if isinstance(x, list) and x else np.nan)

    # Calculate the subleading (second largest) PT directly
    df_merged["Subleading_Jet_UT_PT"] = df_merged["Jet_UT_PT"].apply(lambda x: sorted(x, reverse=True)[1] if isinstance(x, list) and len(x) > 1 else np.nan)
    
    ###Masses of the Tau Tagged
    # Zip Jet_Mass and Jet_TauTagged columns
    df_merged['Jet_TauTagged_Mass'] = df_merged.apply(lambda row: [mass for mass, tag in zip(row['Jet_Mass'], row['Jet_TauTagged']) if tag == 1], axis=1)
    df_merged['Jet_UT_Mass'] = df_merged.apply(lambda row: [mass for mass, tag in zip(row['Jet_Mass'], row['Jet_TauTagged']) if tag == 0], axis=1)
    
    # Calculate Jet_TauTagged_ET using the formula
    df_merged['Jet_TauTagged_ET'] = df_merged.apply(lambda row: [(mass**2 + pt**2)**0.5 for mass, pt in zip(row['Jet_TauTagged_Mass'], row['Jet_TauTagged_PT'])], axis=1) 


    
    ### Cos theta of the leading tau tagged jet
    
    df_merged['TauTaggedTheta'] = df_merged['Jet_TauTagged_Eta'].apply(lambda etas: [ 2 * np.arctan(np.exp(-eta)) for eta in etas])
    
    df_merged['TauTaggedCos_Theta'] = df_merged['TauTaggedTheta'].apply(lambda thetas: [np.cos(theta) for theta in thetas])
    
    # Combine TauTaggedPT with TauTaggedCos_Theta using zip to pair PT with cos(theta)
    df_merged['TauTagged_PT_CosTheta'] = df_merged.apply(lambda row: list(zip(row['Jet_TauTagged_PT'], row['TauTaggedCos_Theta'])), axis=1)
        
    # Find the Leading Tau Tagged PT
    df_merged['Leading_TauTaggedPT'] = df_merged['Jet_TauTagged_PT'].apply(lambda pts: max(pts) if pts else None)
    
    # Find the corresponding Leading Cos Theta of the Leading Tau Tagged Jet
    df_merged['LeadingTauTaggedCos_Theta'] = df_merged['TauTagged_PT_CosTheta'].apply(lambda pairs: max(pairs, key=lambda x: x[0])[1] if pairs else None) 
    
    #Combining the PT'S and ET'S
    df_merged['TauTagged_PT_ET'] = df_merged.apply(lambda row: list(zip(row['Jet_TauTagged_PT'], row['Jet_TauTagged_ET'])), axis=1)
    
    #Now finding the et of the leading tau jet
    df_merged['LeadingTauTagged_ET'] = df_merged['TauTagged_PT_ET'].apply(lambda pairs: max(pairs, key=lambda x: x[0])[1] if pairs else None) 
    
    df_merged['Sub_LeadingTauTagged_ET'] = df_merged['TauTagged_PT_ET'].apply(
    lambda pairs: min(pairs, key=lambda x: x[0])[1] if len(pairs) > 1 else None
)

    

    
    #### Counting leptons and jets
    
    count_lep = [len(k) for k in df_merged["Combined_PT"]]
    count_jet = [len(k) for k in df_merged["Jet_PT"]]

    df_merged["Lepton Count"] = count_lep
    df_merged["Jet Count"] = count_jet
    
    
    
    # Drop the columns that are no longer needed
    drop_cols = [
        "Sorted_PT_Theta", "Combined_PT_Theta", "Combined_Theta",
        "Electron_Eta", "Muon_Eta",
        "Sorted_PT_Charge",'TauTagged_PT_CosTheta','TauTagged_PT_ET'
    ]
    df_merged = df_merged.drop(columns=drop_cols)
    

    

    # Clean up
    file.Close()
    

    return df_merged 



