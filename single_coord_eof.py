import numpy as np
from argo_nn.TS_Parser import *
import pickle

def single_coord():
    N = 234
    L = 33
    Y_t = np.zeros((N, L))
    Y_s = np.zeros((N, L))

    sum = 0
    for year in range(2004, 2024):
        for month in range(1, 13):
            if (year == 2023 and month > 6):
                break
            
            data = TS_Grabber(year, month)
            t, s = data.get_profiles(lon=135, lat=30)

            Y_t[sum] = t
            Y_s[sum] = s
            sum += 1

    Y_t = Y_t.T
    Y_s = Y_s.T
    print('sum:', sum)
    print('temperature matrix has dimensions', Y_t.shape)
    U_t, amp_t, var_t = get_decomp(Y_t)
    
    print('\nsalinity matrix has dimensions', Y_s.shape)
    U_s, amp_s, var_s = get_decomp(Y_s)

    with open('testing-notebooks/ut.pkl', 'wb') as f:
        pickle.dump(U_t, f)

    with open('testing-notebooks/ampt.pkl', 'wb') as f:
        pickle.dump(amp_t, f)

    with open('testing-notebooks/us.pkl', 'wb') as f:
        pickle.dump(U_s, f)

    with open('testing-notebooks/amps.pkl', 'wb') as f:
        pickle.dump(amp_s, f)

single_coord()