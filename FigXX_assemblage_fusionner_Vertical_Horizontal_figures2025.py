#### read first the code of the figure ############
##### Abundance , see below for biovolume #########
###################################################
import cv2
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from PIL import Image
# Chemin du dossier contenant les images
path='C:/Users\syawo\Downloads\Figures_pdf_Angola_Biogeosciences_submission2025/Interp_Cluster/'
pathos= os.listdir(path)
images_files = glob.glob(path + "*2025.png")  # Get all PNG files in the current directory
# Charger les images
c1=0
for filename in images_files:
    c1=c1+1
    exec('img'+ str(c1) + ' = cv2.imread(filename)')
    #exec('img' + str(c1)+ '= cv2.cvtColor(' + 'img'+ str(c1) + ', cv2.COLOR_BGR2RGB)')  # Convert to correct color format

'''# Resize both images to the same width
width = min(img1.shape[1], img2.shape[1])
img1 = cv2.resize(img1, (width, img1.shape[0]))
img2 = cv2.resize(img2, (width, img2.shape[0]))'''
# crope image
#cropped_img = img[:, :-20]  # Adjust value as needed
# Crop the image (left, upper, right, lower)
# Assembler verticalement
h1, w1 = img1.shape[:2]
#v_concat1 = np.vstack((img1[:-380, 80:w1]  , img2[:, 80:w1])) # Haut, bas, droite, gauche
#v_concat2 = np.vstack((img3[:-380, 80:w1]  , img4[:, 80:w1]))
a=15
v_concat1 = np.vstack((img1[:-385, 50-a:w1-90-a], img3[:, 50-a:w1-90-a])) # Haut:bass, droite:gauche
v_concat2 = np.vstack((img2[:-385, 50-a:w1-90-a], img4[:, 50-a:w1-90-a]))

# Assembler horizontalement
h_concat = np.hstack((v_concat1, v_concat2))
# save image without showing it
cv2.imwrite(path+'Fig6abundBis.png', h_concat)  # Save image

#### read first the code of the figure #
###### Biovolume #######################
########################################
import cv2
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from PIL import Image
# Chemin du dossier contenant les images
path='C:/Users\syawo\Downloads\Figures_pdf_Angola_Biogeosciences_submission2025/Interp_ClusterBiovolume/'
pathos= os.listdir(path)
images_files = glob.glob(path + "*2025.png")  # Get all PNG files in the current directory
# Charger les images
c1=0
for filename in images_files:
    c1=c1+1
    exec('img'+ str(c1) + ' = cv2.imread(filename)')
    #exec('img' + str(c1)+ '= cv2.cvtColor(' + 'img'+ str(c1) + ', cv2.COLOR_BGR2RGB)')  # Convert to correct color format

'''# Resize both images to the same width
width = min(img1.shape[1], img2.shape[1])
img1 = cv2.resize(img1, (width, img1.shape[0]))
img2 = cv2.resize(img2, (width, img2.shape[0]))'''
# crope image
#cropped_img = img[:, :-20]  # Adjust value as needed
# Crop the image (left, upper, right, lower)
# Assembler verticalement
h1, w1 = img1.shape[:2]
#v_concat1 = np.vstack((img1[:-380, 80:w1]  , img2[:, 80:w1])) # Haut, bas, droite, gauche
#v_concat2 = np.vstack((img3[:-380, 80:w1]  , img4[:, 80:w1]))
a=15
v_concat1 = np.vstack((img1[:-385, 50-a:w1-90-a], img3[:, 50-a:w1-90-a])) # Haut:bass, droite:gauche
v_concat2 = np.vstack((img2[:-385, 50-a:w1-90-a], img4[:, 50-a:w1-90-a]))
# Assembler horizontalement
h_concat = np.hstack((v_concat1, v_concat2))
# save image without showing it
cv2.imwrite(path+'Fig6biovBis.png', h_concat)  # Save image

###############################################################
############## Physical parameters 0-100m #############################
import cv2
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from PIL import Image
# Chemin du dossier contenant les images
path='C:/Users\syawo\Downloads\Figures_pdf_Angola_Biogeosciences_submission2025\Interp_PhysqiueEnv/'
pathos= os.listdir(path)
images_files = glob.glob(path + "*2025lARS100m*.png")  # Get all PNG files in the current directory
# Charger les images
c1=0
for filename in images_files:
    c1=c1+1
    exec('img'+ str(c1) + ' = cv2.imread(filename)')
    #exec('img' + str(c1)+ '= cv2.cvtColor(' + 'img'+ str(c1) + ', cv2.COLOR_BGR2RGB)')  # Convert to correct color format

'''# Resize both images to the same width
width = min(img1.shape[1], img2.shape[1])
img1 = cv2.resize(img1, (width, img1.shape[0]))
img2 = cv2.resize(img2, (width, img2.shape[0]))'''
# crope image
#cropped_img = img[:, :-20]  # Adjust value as needed
# Crop the image (left, upper, right, lower)
# Assembler verticalement
h1, w1 = img1.shape[:2]
#v_concat1 = np.vstack((img1[:-380, 80:w1]  , img2[:, 80:w1])) # Haut, bas, droite, gauche
#v_concat2 = np.vstack((img3[:-380, 80:w1]  , img4[:, 80:w1]))
a=15
#v_concat1 = np.vstack((img1[:-300, 50-a:w1-90-a], img3[:, 50-a:w1-90-a])) # Haut:bass, droite:gauche
#v_concat2 = np.vstack((img2[:-300, 50-a:w1-90-a], img4[:, 50-a:w1-90-a]))

v_concat1 = np.vstack((img6[:-385, 50-a:w1-90-a], img2[:-385, 50-a:w1-90-a], img1[:, 50-a:w1-90-a])) # Haut:bass, droite:gauche
v_concat2 = np.vstack((img5[:-385, 50-a:w1-90-a], img3[:-385, 50-a:w1-90-a], img4[:, 50-a:w1-90-a]))

# Assembler horizontalement
h_concat = np.hstack((v_concat1, v_concat2))
# save image without showing it
cv2.imwrite(path+'Fig3BisLARS.png', h_concat)  # Save image

###############################################################
############## Physical parameters 0-1000m #############################
import cv2
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from PIL import Image
# Chemin du dossier contenant les images
path='C:/Users\syawo\Downloads\Figures_pdf_Angola_Biogeosciences_submission2025\Interp_PhysqiueEnv/'
pathos= os.listdir(path)
images_files = glob.glob(path + "*lARS1000*.png")  # Get all PNG files in the current directory
# Charger les images
c1=0
for filename in images_files:
    c1=c1+1
    exec('img'+ str(c1) + ' = cv2.imread(filename)')
    #exec('img' + str(c1)+ '= cv2.cvtColor(' + 'img'+ str(c1) + ', cv2.COLOR_BGR2RGB)')  # Convert to correct color format

'''# Resize both images to the same width
width = min(img1.shape[1], img2.shape[1])
img1 = cv2.resize(img1, (width, img1.shape[0]))
img2 = cv2.resize(img2, (width, img2.shape[0]))'''
# crope image
#cropped_img = img[:, :-20]  # Adjust value as needed
# Crop the image (left, upper, right, lower)
# Assembler verticalement
h1, w1 = img1.shape[:2]
#v_concat1 = np.vstack((img1[:-380, 80:w1]  , img2[:, 80:w1])) # Haut, bas, droite, gauche
#v_concat2 = np.vstack((img3[:-380, 80:w1]  , img4[:, 80:w1]))
a=50
#v_concat1 = np.vstack((img1[:-300, 50-a:w1-90-a], img3[:, 50-a:w1-90-a])) # Haut:bass, droite:gauche
#v_concat2 = np.vstack((img2[:-300, 50-a:w1-90-a], img4[:, 50-a:w1-90-a]))

v_concat1 = np.vstack((img3[:-345, 50-a:w1-90-a], img2[:-345, 50-a:w1-90-a], img1[:, 50-a:w1-90-a])) # Haut:bass, droite:gauche
#v_concat2 = np.vstack((img5[:-345, 50-a:w1-90-a], img3[:-345, 50-a:w1-90-a], img4[:, 50-a:w1-90-a]))

# Assembler horizontalement
#h_concat = np.hstack((v_concat1, v_concat2))
# save image without showing it
cv2.imwrite(path+'Fig4LARS.png', v_concat1)  # Save image

###################################################################
#################### Changer le fond ##############################
###################################################################
import cv2
import numpy as np

path='C:/Users\syawo\Downloads\Figures_pdf_Angola_Biogeosciences_submission2025/'
# Charger l'image
image = cv2.imread(path+"FIG5.png")
h1, w1 = image.shape[:2]
v1 = image[:, :w1-2800]  # Haut:bass, droite:gauche
v2 = image[:, 1500:w1]
v_concat1 = np.hstack((v1,v2))
cv2.imwrite(path+'Fig5checkCUT.png', v_concat1)  # Save image
cv2.imwrite(path+'Fig5checkCUT1.png', v1)  # Save image
cv2.imwrite(path+'Fig5checkCUT2.png', v2)  # Save image
# after using image windows to change contraste
v1f = cv2.imread(path+"Fig5checkCUT1fond3.png")
v2f = cv2.imread(path+"Fig5checkCUT2fond.png")

v1f = cv2.imread(path+"Fig5checkCUT1_keep.png")
v2f = cv2.imread(path+"Fig5checkCUT2.png")
v_concat1f = np.hstack((v1f,v2f))
cv2.imwrite(path+'Fig5fond_used5_keep.png', v_concat1f)  # Save image

''''
# Convertir en HSV pour segmenter l'arrière-plan
hsv = cv2.cvtColor(v1, cv2.COLOR_BGR2HSV)
# Définir la plage de couleur pour détecter le fond (ex: blanc)
lower_bound = np.array([0, 0, 200])  # Blanc minimum
upper_bound = np.array([180, 50, 255])  # Blanc maximum
# Créer un masque
mask = cv2.inRange(hsv, lower_bound, upper_bound)
# Inverser le masque pour garder l'objet
mask_inv = cv2.bitwise_not(mask)
# Créer une image de fond (ex: bleu)
#background = np.full(image.shape, (255, 0, 0), dtype=np.uint8)  # Bleu
background = np.full(v1.shape,(255, 255, 255) , dtype=np.uint8)  # gris
#(0, 157, 102)
# Appliquer le masque sur l'image originale et le fond
fg = cv2.bitwise_and(v1, v1, mask=mask_inv)
bg = cv2.bitwise_and(background, background, mask=mask)
# Fusionner les deux
final = cv2.add(fg, bg)
# Sauvegarder et afficher l'image
cv2.imwrite(path+"output.jpg", final) '''
