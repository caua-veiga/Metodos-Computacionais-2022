import matplotlib.pyplot as plt 
import numpy as np 

def Jacobi(grid, target=1e-6):
    M = grid.shape[0] # grid side
    #target = 1e-6 # Target accuracy
        
    #Create arrays to hold potential values
    phi = grid.copy() # 100x100 zeros grid
        
    # Main loop
    delta = 1.0
    while delta>target:
        phiprime = phi.copy()
        
        phi[1:-1,1:-1] = (phi[1:-1,2:] + phi[1:-1,:-2] + phi[2:,1:-1] + phi[:-2,1:-1])/4

        # Calculate maximum difference from old values 
        delta = np.amax(abs(phi-phiprime))     

    return phi 

def gsidel(grid, target=1e-6):
    M = grid.shape[0] # grid side
    #target = 1e-6 # Target accuracy
    
    #Create arrays to hold potential values
    phi = grid.copy() # 100x100 zeros grid
    
    # Main loop
    delta = 1.0
    while delta>target:
        phiprime = phi.copy()
        # Calculate new values of the potential 
        for i in range(M-1):
            for j in range(M-1):
                if i==0 or i==M or j==0 or j==M:
                    pass
                else:
                    phi[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])/4

        # Calculate maximum difference from old values 
        delta = np.amax(abs(phi-phiprime))     

    return phi 


def SOR(grid, w, target=1e-6):
    M = grid.shape[0] # grid side
    #target = 1e-6 # Target accuracy
    
    #Create arrays to hold potential values
    phi = grid.copy() # 100x100 zeros grid
    
    
    # Main loop
    delta = 1.0
    while delta>target:
        phiprime = phi.copy()
        # Calculate new values of the potential 
        for i in range(M-1):
            for j in range(M-1):
                if i==0 or i==M or j==0 or j==M:
                    pass
                else:
                    phi[i,j] = (phi[i+1,j] + phi[i-1,j] + phi[i,j+1] + phi[i,j-1])*(1+w)/4 - w*phi[i,j]

        # Calculate maximum difference from old values 
        delta = np.amax(abs(phi-phiprime))   
    
    return phi


def plot(image, fsize=(15,8),gray=True):

    plt.figure(figsize=fsize)
    plt.imshow(image)
    if gray:
        plt.gray()
    plt.show()