import os

# Set the folder path
folder_path = "/home/markus/Malabio/DiffSBDD/results/LRRK2/inhibitors"

# Get a list of all files in the folder
files = os.listdir(folder_path)

# Iterate over the files and change mol in the last part of the file name to sdf
for filename in files:
    # rename the file to 8TZC instead of mol
    new_file = filename.replace('G2019S', 'LRRK2')

    # Form the full file paths
    new_file = os.path.join(folder_path, new_file)
    old_file = os.path.join(folder_path, filename)

    # Rename the file
    os.rename(old_file, new_file)
    print(f"Renamed: {filename} -> {new_file}")

# # Iterate over the files and remove the extra extension
# for filename in files:
#     # Split the file name and extension
#     name_parts = filename.split('.')
    
#     # Check if there are two or more extensions (e.g., 'file.txt.txt')
#     if len(name_parts) > 2:
#         # Remove the last extension
#         new_filename = '.'.join(name_parts[:-1])
        
#         # Ensure there is no trailing period
#         if new_filename.endswith('.'):
#             new_filename = new_filename[:-1]
        
#         # Form the full file paths
#         old_file = os.path.join(folder_path, filename)
#         new_file = os.path.join(folder_path, new_filename)
        
#         # Rename the file
#         os.rename(old_file, new_file)
#         print(f"Renamed: {filename} -> {new_filename}")
