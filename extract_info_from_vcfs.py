import glob
import pandas as pd
import vcf
import sys
import os
import re

"""
loop though vcf files in the folder "sys.argv[1]"
extract key formation, "AD" "DP" "ALT" "RED" "PS"...
sum up the DP on the same position, sum up AD with the same position and same alternative allels
"""

##sys.argv[1] is the fold path to all vcf files
folder_path = sys.argv[1] + '/vcfs'

##create a list of all paths of vcf files
file_list = glob.glob(folder_path + "/*")

#a dictioanry { "ref pos alt" :  [ad, ref, alt, pos] }
dict_count = {}

# a dictioanry includes depth count for each position
dict_pos = {}

#loop through vcf files
for file_path in file_list:
    #read file
    vcf_reader = vcf.Reader(open(file_path))
    #loop record in vcf the file
    for record in vcf_reader:
        #extract key information
        dp = record.INFO['DP']
        ads = record.INFO['AD']
        alts = record.ALT
        ref = record.REF
        pos = record.POS

        #only keey s protein gene info
        if pos < 21563 or pos > 25384:
            continue
      
        #check if this position is detected in previous record
        #if yes, add the dp values to existing values
        #if no, add the key(nucleutide position) and value(sequence depth) to dictionary
        if pos in dict_pos:
            dict_pos[pos] += dp
        else:     
            dict_pos[pos] = dp
        
        # alts is a list with alternative allels
        # ads is a list with allele depth, ads[0] is the depth of reference allele depth
        # count is used for record the postion of alternative allels in the list
        count = 1
        # loop alternative allels for each record
        for alt in alts:
            key = str(ref) + str(pos) +str(alt)
            #if this allele mutation has been detected, add the allele depth to existing value in dict
            if key in dict_count:
                #ads[count] represent the depth to this specific allele           
                dict_count[key][0] += ads[count]
            #if not, create a new one
            else:
                dict_count[key] = [ads[count], ref, alt, pos]           
            count += 1

#convert the dictionary to dataframe
col_names = ["AD", 'REF', 'ALT', 'POS', 'DP']
result = pd.DataFrame(columns=col_names)
#loop though the dictionary
for key in dict_count:
    pos_num = re.findall(r'\d+', key)
    pos_dp = dict_pos[int(pos_num[0])]  
    new_row = dict_count[key]
    # the sequencing depth is saved seperately, append it 
    new_row.append(int(pos_dp))
    result.loc[key] = new_row

#specify the output name
output_name = os.path.basename(sys.argv[1])
result.to_csv(output_name+'.csv')
