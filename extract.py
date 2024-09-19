# import os
# import xmltodict
# import json

# # Define the root directory containing XML files
# xml_root_dir = r'C:\Users\bizor\Documents\GDCdata\TCGA-COAD\Clinical\Clinical_Supplement'

# # Function to convert XML to JSON
# def xml_to_json(xml_file):
#     with open(xml_file, 'r', encoding='utf-8') as file:
#         xml_content = file.read()
#         json_content = xmltodict.parse(xml_content)
#         return json_content

# # Function to check keys in a JSON file
# def check_keys(json_data):
#     patient_data = json_data.get("coad:tcga_bcr", {}).get("coad:patient", {})
#     vital_status = patient_data.get("clin_shared:vital_status", {}).get("#text")
#     days_to_death = patient_data.get("clin_shared:days_to_death", {}).get("#text")
#     days_to_followup = patient_data.get("clin_shared:days_to_last_followup", {}).get("#text")

#     if vital_status and days_to_death and days_to_followup:
#         return "Has all keys"
#     else:
#         return "Missing keys"

# # Process all XML files in the directory and subdirectories, convert to JSON, and categorize based on keys
# def process_and_categorize_xml_to_json(root_directory):
#     has_all_keys = []
#     missing_keys = []
    
#     for subdir, _, files in os.walk(root_directory):
#         for filename in files:
#             if filename.endswith('.xml'):
#                 # Convert XML to JSON
#                 xml_file = os.path.join(subdir, filename)
#                 json_data = xml_to_json(xml_file)
                
#                 # Save JSON file
#                 json_file = os.path.splitext(xml_file)[0] + '.json'
#                 with open(json_file, 'w', encoding='utf-8') as file:
#                     json.dump(json_data, file, indent=4)
                
#                 # Categorize JSON based on the presence of vital keys
#                 result = check_keys(json_data)
#                 if result == "Has all keys":
#                     has_all_keys.append(filename)
#                 else:
#                     missing_keys.append(filename)
    
#     # Save the categorization results to two separate files
#     with open("has_all_keys.txt", "w") as has_all_file:
#         for filename in has_all_keys:
#             has_all_file.write(f"{filename}\n")
    
#     with open("missing_keys.txt", "w") as missing_keys_file:
#         for filename in missing_keys:
#             missing_keys_file.write(f"{filename}\n")

# # Run the conversion and categorization process
# process_and_categorize_xml_to_json(xml_root_dir)

# import os
# import xmltodict
# import json

# # Define the root directory containing XML files
# xml_root_dir = r'C:\Users\bizor\Documents\GDCdata\TCGA-COAD\Clinical\Clinical_Supplement'

# # Function to convert XML to JSON
# def xml_to_json(xml_file):
#     with open(xml_file, 'r', encoding='utf-8') as file:
#         xml_content = file.read()
#         json_content = xmltodict.parse(xml_content)
#         return json_content

# # Enhanced function to check keys in a JSON file
# def check_keys(json_data, filename):
#     patient_data = json_data.get("coad:tcga_bcr", {}).get("coad:patient", {})
    
#     vital_status = patient_data.get("clin_shared:vital_status", {}).get("#text")
#     days_to_death = patient_data.get("clin_shared:days_to_death", {}).get("#text")
#     days_to_followup = patient_data.get("clin_shared:days_to_last_followup", {}).get("#text")

#     # Log if any key is missing or contains problematic values
#     if not vital_status:
#         print(f"{filename} is missing vital_status")
#     if not days_to_death:
#         print(f"{filename} is missing days_to_death")
#     if not days_to_followup:
#         print(f"{filename} is missing days_to_followup")
    
#     # Handle "Not Applicable" for days_to_death if the patient is alive
#     if vital_status == "Alive" and (days_to_death == "Not Applicable" or not days_to_death) and days_to_followup:
#         return "Has all keys"
    
#     # General check for all required keys
#     if vital_status and days_to_death and days_to_followup:
#         return "Has all keys"
#     else:
#         return "Missing keys"

# # Process all XML files in the directory and subdirectories, convert to JSON, and categorize based on keys
# def process_and_categorize_xml_to_json(root_directory):
#     has_all_keys = []
#     missing_keys = []
    
#     for subdir, _, files in os.walk(root_directory):
#         for filename in files:
#             if filename.endswith('.xml'):
#                 # Convert XML to JSON
#                 xml_file = os.path.join(subdir, filename)
#                 json_data = xml_to_json(xml_file)
                
#                 # Save JSON file
#                 json_file = os.path.splitext(xml_file)[0] + '.json'
#                 with open(json_file, 'w', encoding='utf-8') as file:
#                     json.dump(json_data, file, indent=4)
                
#                 # Categorize JSON based on the presence of vital keys
#                 result = check_keys(json_data, filename)
#                 if result == "Has all keys":
#                     has_all_keys.append(filename)
#                 else:
#                     missing_keys.append(filename)
    
#     # Save the categorization results to two separate files
#     with open("has_all_keys.txt", "w") as has_all_file:
#         for filename in has_all_keys:
#             has_all_file.write(f"{filename}\n")
    
#     with open("missing_keys.txt", "w") as missing_keys_file:
#         for filename in missing_keys:
#             missing_keys_file.write(f"{filename}\n")

# # Run the conversion and categorization process
# process_and_categorize_xml_to_json(xml_root_dir)


import os
import xmltodict
import json

# Define the root directory containing XML files
xml_root_dir = r'C:\path\to\Clinical_Supplement'

# Function to convert XML to JSON
def xml_to_json(xml_file):
    with open(xml_file, 'r', encoding='utf-8') as file:
        xml_content = file.read()
        json_content = xmltodict.parse(xml_content)
        return json_content

# Enhanced function to check keys in a JSON file
def check_keys(json_data, filename):
    patient_data = json_data.get("coad:tcga_bcr", {}).get("coad:patient", {})
    
    vital_status = patient_data.get("clin_shared:vital_status", {}).get("#text")
    days_to_death = patient_data.get("clin_shared:days_to_death", {}).get("#text")
    days_to_followup = patient_data.get("clin_shared:days_to_last_followup", {}).get("#text")

    # Log if any key is missing or contains problematic values
    if not vital_status:
        print(f"{filename} is missing vital_status")
    if not days_to_death:
        print(f"{filename} is missing days_to_death")
    if not days_to_followup:
        print(f"{filename} is missing days_to_followup")
    
    # Handle "Not Applicable" for days_to_death if the patient is alive
    if vital_status == "Alive" and (days_to_death == "Not Applicable" or not days_to_death) and days_to_followup:
        return "Has all keys"
    
    # General check for all required keys
    if vital_status and days_to_death and days_to_followup:
        return "Has all keys"
    else:
        return "Missing keys"

# Process all XML files in the directory and subdirectories, convert to JSON, and categorize based on keys
def process_and_categorize_xml_to_json(root_directory):
    has_all_keys = []
    missing_keys = []
    
    for subdir, _, files in os.walk(root_directory):
        for filename in files:
            if filename.endswith('.xml'):
                # Convert XML to JSON
                xml_file = os.path.join(subdir, filename)
                json_data = xml_to_json(xml_file)
                
                # Save JSON file
                json_file = os.path.splitext(xml_file)[0] + '.json'
                with open(json_file, 'w', encoding='utf-8') as file:
                    json.dump(json_data, file, indent=4)
                
                # Categorize JSON based on the presence of vital keys
                result = check_keys(json_data, filename)
                if result == "Has all keys":
                    # Store the JSON filename instead of the XML filename
                    has_all_keys.append(os.path.basename(json_file))
                else:
                    missing_keys.append(os.path.basename(json_file))
    
    # Save the categorization results to two separate files
    with open("has_all_keys.txt", "w") as has_all_file:
        for json_filename in has_all_keys:
            has_all_file.write(f"{json_filename}\n")
    
    with open("missing_keys.txt", "w") as missing_keys_file:
        for json_filename in missing_keys:
            missing_keys_file.write(f"{json_filename}\n")

# Run the conversion and categorization process
process_and_categorize_xml_to_json(xml_root_dir)
