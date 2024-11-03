import re

output = """
1885 LEU A
1886 GLY A
1887 ASP A
1888 GLY A
1889 SER A
1890 PHE A
1891 GLY A
1892 SER A
1893 VAL A
1895 ARG A
1902 GLU A
1903 VAL A
1904 ALA A
1905 VAL A
1906 LYS A
1933 ILE A
1947 MET A
1948 GLU A
1949 LEU A
1950 ALA A
1951 SER A
1952 LYS A
1953 GLY A
1955 LEU A
1956 ASP A
1957 ARG A
1994 ASP A
1996 LYS A
1998 HIS A
1999 ASN A
2000 VAL A
2001 LEU A
2002 LEU A
2003 PHE A
2016 ALA A
2017 ASP A
2018 TYR A
2019 SER A
2020 ILE A
2033 GLU A
2034 GLY A
2035 THR A
"""

# Find all unique elements in the output (all unique lines)

# Split the output into lines
lines = output.splitlines()

# Create a set of unique lines
unique_lines = list(set(lines))

# Remove the empty lines
unique_lines = [line for line in unique_lines if line]

# Find the integers in the lines
integers = [int(re.findall(r'\d+', line)[0]) for line in unique_lines]

# Sort the lines by the integers
sorted_lines = [line for _, line in sorted(zip(integers, unique_lines))]

# Print all lines on separate lines
for line in sorted_lines:
    print(line)

# Make one line pocket resi list e.g. A:1 A:2 A:3 
resi_list = [f"{line[-1]}:{int(line.split()[0])}" for line in sorted_lines]
resi_list_str = "\", \"".join(resi_list)
resi_list_str = "\"" + resi_list_str + "\""

resi_list_str2 = " ".join(resi_list)

# Print the resi list
print(resi_list_str)
print(resi_list_str2)