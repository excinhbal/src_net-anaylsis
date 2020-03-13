
import numpy as np

def convert_to_adjacency_matrix(syn_array, N):

    assert np.shape(syn_array) == tuple([N*(N-1)]), \
        "Synpase array doesn't have the expected shape"

    # construct i,j arrays
    i = np.repeat(range(N), N)
    j = np.tile(range(N), N)

    idx = (i!=j)
    i,j = i[idx], j[idx]

    adj_matrix = np.zeros((N,N))

    for syn_el, row_id, col_id in zip(syn_array,i,j):
        adj_matrix[row_id, col_id] = syn_el

    return adj_matrix
        
    

# import pickle
# with open('synee_e.p', 'rb') as pfile:
#     sa= picke.load(pfile)
# with open('synee_a.p', 'rb') as pfile:
#     sa= picke.load(pfile)
# with open('synee_a.p', 'rb') as pfile:
#     sa= picklee.load(pfile)
# with open('synee_a.p', 'rb') as pfile:
#     sa= pickle.load(pfile)
# sa
# sa.keys()
# sa['i']
# np.shape(sa['i'])
# sa.keys()
# np.shape(sa['syn_active'])
# np.shape(sa['syn_active'][-1])
# A = sa['syn_active'][-1]
# np.reshape(A, (400,400))
# 400*399
# np.reshape(A, (400,399))
# np.reshape(sa['i'], (400,399))
# np.reshape(sa['j'], (400,399))
# np.reshape(sa['j'], (399,400))
# np.reshape(sa['j'], (400,399))
# for i,j in zip(sa['i'], sa['j']):
#     print(i,j)
# sx = np.zeros((400,400))
# for i,j in zip(sa['i'], sa['j']):
#     sx[i,j] = i*10**3+j
# xs
# sx
# print(sx)
# np.repeat([0,1,2], 3)
# np.tile([0,1,2], 3)
# x = np.tile([0,1,2], 3)
# x =  np.repeat([0,1,2], 3)
# y =  np.tile([0,1,2], 3)
# x
# y
# idx = (x == y)
# idx
# np.invert(idx)
# x[np.invert(idx)]
# y[np.invert(idx)]
# xl = np.repeat(range(400), 400)
# yl = np.tile(range(400), 400)
# idxl = (xl==yl)
# xx = xl[idxl]
# xx
# cl
# xl
# yl
# len(xl)
# yl
# len(yl)
# idx
# idxl
# xx = xl[np.invert(idxl)]
# yy = yl[np.invert(idxl)]
# xx
# yy
# np.testing.assert_array_equal(xx, sa['i'])
# np.testing.assert_array_equal(yy, sa['i'])
# np.testing.assert_array_equal(yy, sa['j'])
# %hist



def topright_axis_off(ax):

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')


if __name__ == "__main__":
    
    convert_to_adjacency_matrix([0,1], 2)
