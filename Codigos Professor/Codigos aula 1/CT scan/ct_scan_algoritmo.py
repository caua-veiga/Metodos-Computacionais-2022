#-*- coding: iso-8859-1 -*-
#
# 
# usa modulo imutils e requer ficheiro sinogram.png com ... sinograma a inverter!
#
"Importar pacotes/fun��es"
import numpy as np 
#import imutils
import matplotlib.pyplot as plt
from skimage.transform import rotate ## Rotina para rodar imagem
import scipy.fftpack as fft          ## Fast Fourier Transform
import scipy.misc                    ## Contem pacote para gravar arrays numpy como imagem .png

from PIL import Image, ImageChops


## Fun��es    


"M�todo da transformada de Radon - transforma uma imagem num sinograma (N�o � usada para a reconstru��o - foi"
"assim que o sinograma de que partimos foi gerado"
def radon(image, steps):        
    #Constroi a Transformada de Radon usando 'steps' projec��es da 'imagem'. 
    height, width = np.shape(image)

    dTheta = 180/steps                                # Incremento do �ngulo das rota��es.
    thetas = np.arange(0,180,dTheta)


    R = np.sqrt(min(image.shape)**2 + max(image.shape)**2)//2
    padded = np.zeros([int(2*R),int(2*R)],dtype=float)
    padded[int(R-height/2):int(R+height/2),int(R-width/2):int(R+width/2)]= image
    
    xx,yy = np.arange(0,padded.shape[0]),np.arange(0,padded.shape[1])
    
    projections = np.zeros((steps,len(xx)))    # Array numpy para acumular projec��es.
    
    
    for j, ang in enumerate(thetas):
        # rotate the image by theta
        rotated = rotate(padded, ang)
        
            
        projections[j,:] = rotated.sum(0)
        
        # c�digo a escrever
    return projections 
    


"Transforma o sinograma para o dom�nio das frequ�ncias usando transformada de Fourier"
def fft_translate(projs):
    # Cria FFTs 1-D de um sinograma
	# ????? c�digo a escrever
    return fft.rfft(projs, axis=1)



"Filtra as  projec��es usando um filtro em rampa"
def ramp_filter(ffts):
    # Filtro Rampa de um array 2-D de FFTs 1-D (FFTs 1D ao longo das linhas)
    ramp =  np.floor(np.arange(0.5, ffts.shape[1]//2 + 0.1, 0.5))  # filtro
    return ffts * ramp

"Regressa ao dom�nio espacial usando uma transformada de Fourier inversa"
def inverse_fft_translate(operator):
    return fft.irfft(operator, axis=1) # fft filtrada



"Reconstr�i a imagem retroprojectando as projec��es filtradas (ou n�o)"
def back_project(sinograma):
    laminograma = np.zeros((sinograma.shape[1],sinograma.shape[1])) ##
    dTheta =   -180.0 / sinograma.shape[0] # dTheta negativo pois estamos a tomar os passos contrarios que fizemos no
                                        # sinograma           ##
    for i in range(sinograma.shape[0]):

        temp = rotate(np.tile(sinograma[i],(sinograma.shape[1],1)), dTheta*i)  # pode usar rotina rotate do scipy!!
        laminograma += temp

        if i in np.arange(0,sinograma.shape[0],sinograma.shape[0]//10):
            plt.figure(figsize=(4,4))
            plt.title(f'Número de Projeções = {i}')
            plt.imshow(laminograma,cmap='gray')
            plt.show()
    return laminograma



## Programa principal, que usa rotinas anteriores

def tomo(impath,oq):
    # To make sure that the input directory exist
    key = False
    while not key:
        try:
            # ler imagem  do paciente (aqui ainda queremos gerar sinograma) e criar array numpy
            imagem =   np.array(Image.open(impath).convert("L"))          
            key = True
        except Exception as e:
            print(e)
            impath = input("Enter a valid image path: ")


    if oq==0:
        sinograma = radon(imagem, np.amax(imagem.shape))
    if oq==1:
        sinograma = imagem
    #imutils.imshow(sinograma)
    ###################
    # Daqui em diante, reconstru�ao propriamente dita. Na pr�ctica � aqui que come�amos, dado um sinograma.
    #
    "Importar a imagem como array numpy e mostrar a imagem do sinograma original"
    print("Sinograma Original")
    #sinogram =  ?????  # ler sinograma, plot
    plt.figure(figsize=(8,8))
    plt.title('Sinograma Original',fontweight='bold')
    plt.imshow(sinograma,cmap='gray')
    plt.show()


    "Tentar reconstruir a imagem directamente do sinograma, sem usar qualquer tipo de filtragem (retroprojec��o n�o filtrada)"
    print("Reconstruc��o sem filtragem")
    unfiltered_reconstruction =  back_project(sinograma)# retroprojec��o que d� esborratado. Plot, gravar, imagem
    plt.figure(figsize=(8,8))
    plt.title('Backprojeção sem filtragem',fontweight='bold')
    plt.imshow(unfiltered_reconstruction,cmap='gray')
    plt.show()


    "Usar a FFT para passar o sinograma para o dom�nio das frequ�ncias e imprimir o resultado"
    frequency_domain_sinogram =  fft_translate(sinograma)                      # Plot, gravar, imagem      #


    "Filtragem das projec��es no dom�nio das frequ�ncias, multiplicando cada uma pelo filtro em rampa (nas frequ�ncias)"
    filtered_frequency_domain_sinogram = ramp_filter(frequency_domain_sinogram)   # filtrar espectro. Plot, gravar, imagem 


    "Usar a FFT inversa para regressar ao dom�nio espacial"
    filtered_spatial_domain_sinogram = inverse_fft_translate(filtered_frequency_domain_sinogram)     # de volta ao dom�nio espacial. Plot, gravar, imagem 


    "Reconstru��o  da imagem original 2D por retroprojec��o filtrada"
    reconstructed_image =  back_project(filtered_spatial_domain_sinogram) # retroprojectar a filtrada. Plot, gravar, imagem 
    plt.figure(figsize=(8,8))
    plt.title('Backprojeção com filtro rampa no sinograma',fontweight='bold')
    plt.imshow(reconstructed_image,cmap='gray')
    plt.show() 


if __name__ == "__main__":
    print(tomo(input('Path para a imagem: '),int(input('A imagem é um sinograma? (0 -Não / 1 -Sim): '))))

