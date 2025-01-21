import copy
import numpy as np

def convert_to_basis_string(arr):

    arr = copy.deepcopy(arr)
    arr = np.array(arr, dtype=object)
    arr[arr == 0] = 'Z'
    arr[arr == 1] = 'X'

    return arr