import numpy as np 
import matplotlib.pyplot as plt 

def Jacobi(grid, dx , dy, target=1e-6,count=False):
    #Create arrays to hold potential values
    phi = grid.copy() 
    
    if count:
        ct = 0
        # Main loop
        delta = 1.0
        while delta>target:
            phiprime = phi.copy()
            
            phi[1:-1,1:-1] = ((dy**2 * (phi[1:-1, 2:] + phi[1:-1, 0:-2]) + dx**2 * (phi[2:, 1:-1] + phi[0:-2, 1:-1])
                                -grid[1:-1,1:-1] * dx**2 * dy**2) / (2* (dx**2 + dy**2)))

            ct+=1 # Conta iterações
            # Calculate maximum difference from old values 
            delta = np.amax(abs(phi-phiprime))     

        return phi,ct


    # Main loop
    delta = 1.0
    while delta>target:
        phiprime = phi.copy()
        
        phi[1:-1,1:-1] = ((dy**2 * (phi[1:-1, 2:] + phi[1:-1, 0:-2]) + dx**2 * (phi[2:, 1:-1] + phi[0:-2, 1:-1])
                            -grid[1:-1,1:-1] * dx**2 * dy**2) / (2* (dx**2 + dy**2)))

        # Calculate maximum difference from old values 
        delta = np.amax(abs(phi-phiprime))     

    return phi 


def gsidel(grid, dx, dy, target=1e-6,count=False,REDBLACK=False):
    ny = grid.shape[0]
    nx = grid.shape[1]

    #Create arrays to hold potential values
    phi = grid.copy()


    if count:
        ct = 0
        if REDBLACK: # CASO QUEIRA FAZER COM ORDENAÇÃO RED-BLACK
            while delta>target:
                phiprime = phi.copy()
                
                # Update Red-Dots
                # Lines where the first dot is red (excluding the boundaries)
                phi[1:-1:2,1:-1:2]=(dy**2*(phi[1:-1:2,2::2]+phi[1:-1:2,0:-2:2]) + dx**2*(phi[0:-2:2,1:-1:2]+phi[2::2,1:-1:2])
                                    - grid[1:-1:2,1:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))
                    
                # Lines where the first dot is black (excluding the boundaries)
                phi[2:-1:2,2:-1:2]=(dy**2*(phi[2:-1:2,3::2]+phi[2:-1:2,1:-2:2]) + dx**2*(phi[1:-2:2,2:-1:2]+phi[3::2,2:-1:2])
                                    - grid[2:-1:2,2:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))
                
                
                # Update Black-Dots (using the new red-dot's values)
                
                # Lines where the first dot is black (excluding the boundaries)
                phi[2:-1:2,1:-1:2]=(dy**2*(phi[2:-1:2,2::2]+phi[2:-1:2,0:-2:2]) + dx**2*(phi[1:-2:2,1:-1:2]+phi[3::2,1:-1:2])
                                    - grid[2:-1:2,1:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))
                    
                
                # Lines where the first dot is red (excluding the boundaries)
                phi[1:-1:2,2:-1:2]=(dy**2*(phi[1:-1:2,3::2]+phi[1:-1:2,1:-2:2]) + dx**2*(phi[0:-2:2,2:-1:2]+phi[2::2,2:-1:2])
                                    - grid[1:-1:2,2:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))  
                        
                ct +=1 # Conta iterações 

                # Calculate maximum difference from old values 
                delta = np.amax(abs(phi-phiprime)) 
                
                
            return phi,ct


        # Main loop
        delta = 1.0
        while delta>target:
            phiprime = phi.copy()
            # Calculate new values of the potential 
            for i in range(1,ny-1):
                for j in range(1,nx-1):
                    phi[i,j] = ( (phi[i+1,j] + phi[i-1,j])/dy**2 + (phi[i,j+1] + phi[i,j-1])/dx**2 - grid[i,j]* dx**2 * dy**2)/(2/dx**2 + 2/dy**2)

            ct +=1 # Conta iterações

            # Calculate maximum difference from old values 
            delta = np.amax(abs(phi-phiprime))     

        return phi,ct

    # Main loop
    delta = 1.0

    if REDBLACK: # CASO QUEIRA FAZER COM ORDENAÇÃO RED-BLACK
        while delta>target:
            phiprime = phi.copy()
            
            # Update Red-Dots
            # Lines where the first dot is red (excluding the boundaries)
            phi[1:-1:2,1:-1:2]=(dy**2*(phi[1:-1:2,2::2]+phi[1:-1:2,0:-2:2]) + dx**2*(phi[0:-2:2,1:-1:2]+phi[2::2,1:-1:2])
                                - grid[1:-1:2,1:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))
                
            # Lines where the first dot is black (excluding the boundaries)
            phi[2:-1:2,2:-1:2]=(dy**2*(phi[2:-1:2,3::2]+phi[2:-1:2,1:-2:2]) + dx**2*(phi[1:-2:2,2:-1:2]+phi[3::2,2:-1:2])
                                - grid[2:-1:2,2:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))
            
            
            # Update Black-Dots (using the new red-dot's values)
            
            # Lines where the first dot is black (excluding the boundaries)
            phi[2:-1:2,1:-1:2]=(dy**2*(phi[2:-1:2,2::2]+phi[2:-1:2,0:-2:2]) + dx**2*(phi[1:-2:2,1:-1:2]+phi[3::2,1:-1:2])
                                - grid[2:-1:2,1:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))
                
            
            # Lines where the first dot is red (excluding the boundaries)
            phi[1:-1:2,2:-1:2]=(dy**2*(phi[1:-1:2,3::2]+phi[1:-1:2,1:-2:2]) + dx**2*(phi[0:-2:2,2:-1:2]+phi[2::2,2:-1:2])
                                - grid[1:-1:2,2:-1:2]*dx**2*dy**2) / (2*(dx**2+dy**2))  
                    

            # Calculate maximum difference from old values 
            delta = np.amax(abs(phi-phiprime))     
            
        return phi



    while delta>target:
        phiprime = phi.copy()
        # Calculate new values of the potential 
        for i in range(1,ny-1):
            for j in range(1,nx-1):
                phi[i,j] = ( (phi[i+1,j] + phi[i-1,j])/dy**2 + (phi[i,j+1] + phi[i,j-1])/dx**2 - grid[i,j]* dx**2 * dy**2 )/(2/dx**2 + 2/dy**2)

        # Calculate maximum difference from old values 
        delta = np.amax(abs(phi-phiprime))     

    return phi 


def SOR(grid, dx, dy, w, target=1e-6,count=False,REDBLACK=False):
    
    ny = grid.shape[0]
    nx = grid.shape[1]

    #Create arrays to hold potential values
    phi = grid.copy()


    if count:
        ct = 0

        if REDBLACK:  # CASO QUEIRA FAZER COM ORDENAÇÃO RED-BLACK
            while delta>target:
                phiprime = phi.copy()
                
                # Update Red-Dots
                # Lines where the first dot is red (excluding the boundaries)
                phi[1:-1:2,1:-1:2]=(dy**2*(phi[1:-1:2,2::2]+phi[1:-1:2,0:-2:2]) + dx**2*(phi[0:-2:2,1:-1:2]+phi[2::2,1:-1:2])
                                    - grid[1:-1:2,1:-1:2]*dx**2*dy**2)*(1+w)/(2*(dx**2+dy**2)) - w*phi[1:-1:2,1:-1:2]
                    
                # Lines where the first dot is black (excluding the boundaries)
                phi[2:-1:2,2:-1:2]=(dy**2*(phi[2:-1:2,3::2]+phi[2:-1:2,1:-2:2]) + dx**2*(phi[1:-2:2,2:-1:2]+phi[3::2,2:-1:2])
                                    - grid[2:-1:2,2:-1:2]*dx**2*dy**2)*(1+w)/(2*(dx**2+dy**2)) - w*phi[2:-1:2,2:-1:2]
                
                
                # Update Black-Dots (using the new red-dot's values)
                
                # Lines where the first dot is black (excluding the boundaries)
                phi[2:-1:2,1:-1:2]=(dy**2*(phi[2:-1:2,2::2]+phi[2:-1:2,0:-2:2]) + dx**2*(phi[1:-2:2,1:-1:2]+phi[3::2,1:-1:2])
                                    - grid[2:-1:2,1:-1:2]*dx**2*dy**2)*(1+w)/ (2*(dx**2+dy**2)) - w*phi[2:-1:2,1:-1:2]
                    
                
                # Lines where the first dot is red (excluding the boundaries)
                phi[1:-1:2,2:-1:2]=(dy**2*(phi[1:-1:2,3::2]+phi[1:-1:2,1:-2:2]) + dx**2*(phi[0:-2:2,2:-1:2]+phi[2::2,2:-1:2])
                                    - grid[1:-1:2,2:-1:2]*dx**2*dy**2)*(1+w)/ (2*(dx**2+dy**2)) - w*phi[1:-1:2,2:-1:2]

                ct += 1 # conta iterações
                        
                # Calculate maximum difference from old values 
                delta = np.amax(abs(phi-phiprime))
                
            return phi,ct

        # Main loop
        delta = 1.0
        while delta>target:
            phiprime = phi.copy()
            # Calculate new values of the potential 
            for i in range(1,ny-1):
                for j in range(1,nx-1):
                    phi[i,j] = ( (phi[i+1,j] + phi[i-1,j])/dy**2 + (phi[i,j+1] + phi[i,j-1])/dx**2 - grid[i,j]* dx**2 * dy**2 )*(1+w)/(2/dx**2 + 2/dy**2) - w*phi[i,j]
            
            ct += 1 # conta iterações

            # Calculate maximum difference from old values 
            delta = np.amax(abs(phi-phiprime))     

        return phi,ct 

    # Main loop
    delta = 1.0

    if REDBLACK:  # CASO QUEIRA FAZER COM ORDENAÇÃO RED-BLACK
        while delta>target:
            phiprime = phi.copy()
            
            # Update Red-Dots
            # Lines where the first dot is red (excluding the boundaries)
            phi[1:-1:2,1:-1:2]=(dy**2*(phi[1:-1:2,2::2]+phi[1:-1:2,0:-2:2]) + dx**2*(phi[0:-2:2,1:-1:2]+phi[2::2,1:-1:2])
                                - grid[1:-1:2,1:-1:2]*dx**2*dy**2)*(1+w)/(2*(dx**2+dy**2)) - w*phi[1:-1:2,1:-1:2]
                
            # Lines where the first dot is black (excluding the boundaries)
            phi[2:-1:2,2:-1:2]=(dy**2*(phi[2:-1:2,3::2]+phi[2:-1:2,1:-2:2]) + dx**2*(phi[1:-2:2,2:-1:2]+phi[3::2,2:-1:2])
                                - grid[2:-1:2,2:-1:2]*dx**2*dy**2)*(1+w)/(2*(dx**2+dy**2)) - w*phi[2:-1:2,2:-1:2]
            
            
            # Update Black-Dots (using the new red-dot's values)
            
            # Lines where the first dot is black (excluding the boundaries)
            phi[2:-1:2,1:-1:2]=(dy**2*(phi[2:-1:2,2::2]+phi[2:-1:2,0:-2:2]) + dx**2*(phi[1:-2:2,1:-1:2]+phi[3::2,1:-1:2])
                                - grid[2:-1:2,1:-1:2]*dx**2*dy**2)*(1+w)/ (2*(dx**2+dy**2)) - w*phi[2:-1:2,1:-1:2]
                
            
            # Lines where the first dot is red (excluding the boundaries)
            phi[1:-1:2,2:-1:2]=(dy**2*(phi[1:-1:2,3::2]+phi[1:-1:2,1:-2:2]) + dx**2*(phi[0:-2:2,2:-1:2]+phi[2::2,2:-1:2])
                                - grid[1:-1:2,2:-1:2]*dx**2*dy**2)*(1+w)/ (2*(dx**2+dy**2)) - w*phi[1:-1:2,2:-1:2]
                    

            # Calculate maximum difference from old values 
            delta = np.amax(abs(phi-phiprime))
            
        return phi

    while delta>target:
        phiprime = phi.copy()
        # Calculate new values of the potential 
        for i in range(1,ny-1):
            for j in range(1,nx-1):
                phi[i,j] = ( (phi[i+1,j] + phi[i-1,j])/dy**2 + (phi[i,j+1] + phi[i,j-1])/dx**2 - grid[i,j]* dx**2 * dy**2 )*(1+w)/(2/dx**2 + 2/dy**2) - w*phi[i,j]

        # Calculate maximum difference from old values 
        delta = np.amax(abs(phi-phiprime))     

    return phi 


def plot(image, fsize=(15,8),gray=True):

    plt.figure(figsize=fsize)
    plt.imshow(image)
    if gray:
        plt.gray()
    plt.show()