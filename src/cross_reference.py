import os
from pybedtools import BedTool

# This script finds overlaps between bed files in two directories and writes the results to a new file.
for bed_file in os.listdir(r"C:\Users\yiann\Desktop\projects\TE_annotation\Repeats_IGV\Repeats_IGV\wntA"):
    if not bed_file.endswith(".bed"):
        continue
    bed_file_path = os.path.join(r"C:\Users\yiann\Desktop\projects\TE_annotation\Repeats_IGV\Repeats_IGV\wntA", bed_file)
    name1 = bed_file.replace(".bed", "")
    for prediction_file in os.listdir(r"C:\Users\yiann\Desktop\projects\TE_annotation\TEs_predictions"):
        if not prediction_file.endswith(".bed"):
            continue
        prediction_file_path = os.path.join(r"C:\Users\yiann\Desktop\projects\TE_annotation\TEs_predictions", prediction_file)
        name2 = prediction_file.replace("_TEs_predictions.bed", "")
        if name1 == name2:
            model = BedTool(prediction_file_path)
            rm = BedTool(bed_file_path)
            overlap = model.intersect(rm, wa=True, wb=True, f=0.8, e=True)
        for entry in overlap:
            new_file = os.path.join(r"C:\Users\yiann\Desktop\projects\TE_annotation\TEs_predictions", name1 + "_overlap.bed")
            with open(new_file, 'a') as output_file:
                output_file.write(str(entry) + '\n')
            
        
