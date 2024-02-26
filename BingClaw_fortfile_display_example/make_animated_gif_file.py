import PIL
from PIL import Image
import numpy as np

scenario = "png_files"
image_frames = []

steps = np.arange( 1, 180, 1 )

for k in steps:
    stepstring = "{:04}".format( k )
    filename   = "./" + scenario + "/" + "q" + stepstring + ".png"
    print ( filename )
    newimage   = PIL.Image.open(r'' + filename)
    image_frames.append( newimage )
    ##DONOTUSE! newimage.show()

outfile = "./" + "BingCLAW" + ".gif"
image_frames[0].save( outfile, format = 'GIF',
                      append_images = image_frames[1:],
                      save_all = True, duration = 90, loop = 0 )

