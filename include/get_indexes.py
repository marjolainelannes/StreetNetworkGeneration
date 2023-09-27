##################################################################################
# Study: Disaggregated car fleet model
# Purpose: Function to get indexes of an key outdata of a list
# Author: Marjolaine Lannes
# Creation date: October 20, 2021
# Note: ...
##################################################################################

def get_indexes(my_list, e):
    output=[]
    for idx,element in enumerate(my_list) :
        if element == e :
            output.append(idx)
    return(output)

# output = [idx for idx,key in enumerate(my_list) if key == e]

