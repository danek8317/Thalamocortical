import numpy as np
import h5py as h5

cell_range = [0,1000,1050,1140,1230,1320,1560,2360,2560,3060,3160,3260,3360]
pop_names = ['pyrRS23','pyrFRB23','bask23','axax23','LTS23', 
             'spinstel4', 'tuftIB5', 'tuftRS5', 'nontuftRS6', 
             'bask56', 'axax56', 'LTS56']
num_cmpts = [74, 74, 59, 59, 59, 59, 61, 61, 50, 59, 59, 59]
num_cells = np.diff(cell_range) / 10 #10% MODEL
total_cmpts = list(num_cmpts*num_cells)

def place_electrodes_1D(n):
    '''places n number of electrodes in ID along the column'''
    ele_x = np.ones((n,1))*25.
    ele_y = np.linspace(-2050., 450.0, num=n).reshape(n,1)
    ele_z = np.ones((n,1))*25.
    return np.hstack((ele_x, ele_y, ele_z))

def place_electrodes_2D(nx, ny):
    '''place nx*ny electrodes next to the column - like an MEA'''
    tot_ele = nx*ny
    zz = np.ones((tot_ele,1))*-25.
    xx, yy = np.mgrid[-375:375:np.complex(0,nx), -2250:450:np.complex(0,ny)]
    xx = xx.reshape(tot_ele, 1)
    yy = yy.reshape(tot_ele, 1)
    return np.hstack((xx, yy, zz))

def place_electrodes_3D(nx, ny, nz):
    '''place nx*ny*nz electrodes in the column - like utah array'''
    tot_ele = nx*ny*nz
    #zz = np.ones((tot_ele,1))*-25.
    xx, yy, zz = np.mgrid[-375:375:np.complex(0,nx),
                          -2250:450:np.complex(0,ny), 
                          -25:25:np.complex(0,nz)]
    xx = xx.reshape(tot_ele, 1)
    yy = yy.reshape(tot_ele, 1)
    zz = zz.reshape(tot_ele, 1)
    return np.hstack((xx, yy, zz))

def fetch_mid_pts(h, pop_name):
    all_pts = h['/data/static/morphology/'+pop_name]
    x = (all_pts['x0']+all_pts['x1']) / 2.
    y = (all_pts['y0']+all_pts['y1']) / 2.
    z = (all_pts['z0']+all_pts['z1']) / 2.
    x = x.reshape(x.size,1)
    y = y.reshape(y.size,1)
    z = z.reshape(z.size,1)
    return np.hstack((x, y, z))

# print place_electrodes_1D(10)
# print place_electrodes_2D(3, 3)
#h = h5.File('traub.h5', 'r')
h = h5.File('/media/cchintaluri/PersonalBackup/data/hela_data/small_awake05/1/traub_syn.h5', 'r')

ele_config = '1D'
ANIMATE = False

if ele_config == '1D':
    num_ele = 28
    ele_pos = place_electrodes_1D(num_ele)
    pot_filename = 'pot_sum_1D_'+str(num_ele)+'.npy'
elif ele_config == '2D':
    num_x, num_y = 16, 20
    num_ele = num_x*num_y
    ele_pos = place_electrodes_2D(num_x, num_y)
    pot_filename = 'pot_sum_2D_'+str(num_ele)+'.npy'
elif ele_config == '3D':
    num_x, num_y, num_z = 16, 20, 5
    num_ele = num_x*num_y*num_z
    ele_pos = place_electrodes_3D(num_x, num_y, num_z)
    pot_filename = 'pot_sum_3D_'+str(num_ele)+'.npy'
