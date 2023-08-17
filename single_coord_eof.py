import numpy as np
from argo_nn.TS_Parser import *
import pickle

def single_coord():
    N = 215
    L = 58
    Y_t = np.zeros((N, L))
    Y_s = np.zeros((N, L))

    sum = 0
    for year in range(2004, 2024):
        for month in range(1, 12):
            if (year == 2023 and month > 6):
                break
            
            data = TS_Grabber(year, month)
            t, s = data.get_profiles(data.lat[0], data.lon[0])
            # -64.5, 20.5

            Y_t[sum] = t
            Y_s[sum] = s
            sum += 1

    print('sum:', sum)
    print('temperature matrix has dimensions', Y_t.shape)
    U_t, amp_t, var_t = get_decomp(Y_t)
    print(U_t)
    
    print('\nsalinity matrix has dimensions', Y_s.shape)
    U_s, amp_s, var_s = get_decomp(Y_s)
    print(U_s)

    with open('testing-notebooks/ut.pkl', 'wb') as f:
        pickle.dump(U_t, f)

    with open('testing-notebooks/us.pkl', 'wb') as f:
        pickle.dump(U_s, f)

single_coord()