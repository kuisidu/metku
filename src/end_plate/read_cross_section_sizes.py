# Imported libraries
from re import findall
from re import match

with open('ArcelorMittal_I_cross-sections.txt', 'r') as f:
    lines = f.read().splitlines()

    for line in lines:
        name = []
        if match(r"^IPE", line):
            name = findall(r"^IPE [0-9][0-9]*", line)       # Find profile name
        elif match(r"^HE", line):
            name = findall(r"^HE*\s[0-9]*\s[ABCM]*", line)       # Find profile name

        if name != []:
            name = name[0].replace("\t", "")           # Format profile name
            line = line.replace(name, "")               # Remove profile name from dimension list
            dimensions = findall(r"[\d,]*", line)       # Find profile dimensions
            i = 0
            while i < len(dimensions):
                if dimensions[i] == "":
                    dimensions.remove(dimensions[i])        # Remove spaces from list
                    continue
                i = i + 1

            # Print profile information in tuple table form
            print("\"" + name + "\": {" +
                  "\"h\": " + str(float(dimensions[1].replace(",", "."))) + ", " +
                  "\"b\": " + str(float(dimensions[2].replace(",", "."))) + ", " +
                  "\"t_w\": " + str(float(dimensions[3].replace(",", "."))) + ", " +
                  "\"t_f\": " + str(float(dimensions[4].replace(",", "."))) + ", " +
                  "\"r\": " + str(float(dimensions[5].replace(",", "."))) + "},")
    f.close()