import os
import json
import csv


output_csv_file = r'C:\output\file\destination'

# Path to the root directory containing subfolders with JSON files
json_root_dir = r'C:path\to\TCGA-COAD\Clinical\Clinical_Supplement'

# Path to the has_all_keys.txt file
has_all_keys_file = "has_all_keys.txt"

# Output CSV file path with full path
output_csv_file = r'C:output\file.csv'

# Function to extract relevant values from a JSON file
def extract_values_from_json(json_file):
    with open(json_file, 'r', encoding='utf-8') as file:
        json_data = json.load(file)
    
    # Extract patient data
    patient_data = json_data.get("coad:tcga_bcr", {}).get("coad:patient", {})
    
    # Extract relevant fields
    vital_status = patient_data.get("clin_shared:vital_status", {}).get("#text")
    days_to_death = patient_data.get("clin_shared:days_to_death", {}).get("#text")
    days_to_followup = patient_data.get("clin_shared:days_to_last_followup", {}).get("#text")
    
    # Return the extracted values
    return vital_status, days_to_death, days_to_followup

# Function to find JSON files recursively
def find_json_file(filename, root_directory):
    # Traverse the root directory and subdirectories
    for subdir, _, files in os.walk(root_directory):
        for file in files:
            if file == filename:
                print(f"Found file: {file}")  # Debug print
                return os.path.join(subdir, file)
    print(f"File not found: {filename}")  # Debug print if file is not found
    return None

# Function to process JSON files and write data to CSV
def extract_and_save_to_csv(json_root_directory, has_all_keys_file, output_csv_file):
    # Open the CSV file for writing
    with open(output_csv_file, mode='w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        
        # Write the header row
        writer.writerow(["filename", "vital_status", "days_to_death", "days_to_followup"])
        print(f"Writing CSV to: {output_csv_file}")  # Debug print
        
        # Read the has_all_keys.txt file and process each JSON file
        with open(has_all_keys_file, 'r') as keys_file:
            for json_filename in keys_file:
                json_filename = json_filename.strip()  # Remove any leading/trailing whitespace
                print(f"Processing file: {json_filename}")  # Debug print
                
                # Find the JSON file in the subdirectories
                json_file_path = find_json_file(json_filename, json_root_directory)
                
                if json_file_path:
                    # Extract values from the JSON file
                    vital_status, days_to_death, days_to_followup = extract_values_from_json(json_file_path)
                    
                    # Write the extracted values to the CSV file
                    writer.writerow([json_filename, vital_status, days_to_death, days_to_followup])
                    print(f"Wrote data for: {json_filename}")  # Debug print
                else:
                    print(f"File not found for: {json_filename}")

# Run the extraction and save to CSV process
extract_and_save_to_csv(json_root_dir, has_all_keys_file, output_csv_file)
