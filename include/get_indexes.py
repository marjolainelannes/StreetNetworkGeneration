def get_indexes(my_list, e):
    output=[]
    for idx,element in enumerate(my_list) :
        if element == e :
            output.append(idx)
    return(output)
