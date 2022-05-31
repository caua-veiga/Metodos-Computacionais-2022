import numpy as np
from scipy import interpolate


class adveccao_solver():
    
    def __init__(self, a, tmin, tmax, xmin, xmax, Nx, c ):
        self.a = a
        self.tmin = tmin 
        self.tmax = tmax 
        self.xmin = xmin 
        self.xmax = xmax 
        #self.Nt = Nt
        self.Nx = Nx 
        self.c = c 

    def cond_inicial(self, x, n_func):
        '''
        n_func: 0 == função degrau ;
                1 == sin^2 ;
                2 == cartola ; (activation function)
        '''

        # função que define condições iniciais
        def um(x):
            """Atribui valor de 1.0 a todas as abcissas menores que 0.1"""
            f = np.zeros_like(x)
            f[np.where(x <= 0.1)] = 1.0
            return f

        def dois(x):
            """Uma função suave sin^2 entre x_left e x_right"""
            f       = np.zeros_like(x)
            x_left  = 0.25
            x_right = 0.75
            xm      = (x_right-x_left)/2.0
            f       = np.where((x>x_left) & (x<x_right),\
                            np.sin(np.pi*(x-x_left)/(x_right-x_left))**4,f) 
            return f

        def tres(x):
            """Uma função cartola entre x_left and x_right"""
            f       = np.zeros_like(x)
            x_left  = 0.25
            x_right = 0.75
            xm      = (x_right+x_left)/2.0
            width   = (x_right-x_left)/2.0
            f[np.where(abs(x-xm)< width)] = 1.0
            return f 
        
        func = [um,dois,tres][n_func]
        
        return func(x)

    def ftbs(self,u):
        '''
        Forward time backward space
        '''
        u[1:-1] = (1-self.c)*u[1:-1] + self.c * u[:-2]
        return u[1:-1]

    def ftcs(self,u):
        '''
        Forward time centered space
        '''
        u[1:-1] = u[1:-1] - self.c/2*(u[2:] - u[:-2])
        return u[1:-1]

    def lax_wendroff(self,u): 
        '''
        Lax-Wendroff
        '''
        u[1:-1] = self.c/2.0*(1+self.c)*u[:-2] + (1-self.c**2)*u[1:-1] - self.c/2.0*(1-self.c)*u[2:]
        return u[1:-1]

    def lax_friedrich(self,u):
        '''
        Lax-Friedrich Advection
        '''
        u[1:-1] = (u[:-2] +u[2:])/2.0 -  self.c*(u[2:] - u[:-2])/2.0
        return u[1:-1] 

    def select_solver(self, n):
        '''
        n: 0 == Forward time backward space;
           1 ==  Forward time centered space;
           2 == Lax-Wendroff;
           3 == Lax-Friedrich Advection
        '''

        return [self.ftbs, self.ftcs, self.lax_wendroff, self.lax_friedrich][n]

    def solver(self, initial_cond, method):
        '''
        Return: Solução calculada, solução analítica
        '''

        solver = self.select_solver(method)

        # Discretização
        x    = np.linspace(self.xmin, self.xmax, self.Nx+1)  # discretização do espaço
        dx   = float((self.xmax-self.xmin)/self.Nx)          # passo espacial
        dt   = self.c/self.a*dx                         # passo temporal calculado do critério de estabilidade 
        Nt   = int((self.tmax-self.tmin)/dt)            # número de passos no tempo
        time = np.linspace(self.tmin, self.tmax, Nt)    # vector de instantes de tempo


        uanalytical = np.zeros((len(time), len(x))) # armazena a solução analítica

        u = self.cond_inicial(x, initial_cond)
        import matplotlib.pyplot as plt
        plt.plot(u)
        un = np.zeros((len(time), len(x)))   # armazena a solução numérica

        for i, t in enumerate(time[1:]):
            uanalytical[i,:] = self.cond_inicial(x-self.a*t, initial_cond)  # calcula a solução analítica para este passo temporal
                
            u_bc = interpolate.interp1d(x[-2:], u[-2:]) # interpolar na fronteira direita
            
            u[1:-1] = solver(u[:])       # calcula a solução numérica nos pontos interiores
            u[-1] = u_bc(x[-1] - self.a*dt)   # interpola ao longo da característica para obter o valor na fronteira
            
            un[i,:] = u[:]               # armazena a solução para fazer gráfico
        
        return un, uanalytical

