import numpy as np
import scipy.misc.pilutil as smp
import stepminer

# Create a 1024x1024x3 array of 8 bit unsigned integers
data = np.zeros( (1024,1024,3), dtype=np.uint8 )

data[512,512] = [254,0,0]       # Makes the middle pixel red
data[512,513] = [0,0,255]       # Makes the next pixel blue

img = smp.toimage( data )       # Create a PIL image
img.save("test.png")                      # View in default viewer


# 00 (0): lowY, lowX
# 01 (1): lowY, highX
# 10 (2): highY, lowX
# 11 (3): highY, highX
CLS = {
  CLASSES_I['UNL']: {'name': 'UNL', 'color': (0,0,0), 'draw': (0,1,2,3)},
  CLASSES_I['HIH']: {'name': 'HIH', 'color': (0,0,0), 'draw': (0,1,2,3)},
  CLASSES_I['PC']:  {'name': 'UNL', 'color': (0,0,0), 'draw': (0,1,2,3)},
  CLASSES_I['LIL']: {'name': 'UNL', 'color': (0,0,0), 'draw': (0,1,2,3)},
  CLASSES_I['HIL']: {'name': 'UNL', 'color': (0,0,0), 'draw': (0,1,2,3)},
  CLASSES_I['NC'],  {'name': 'UNL', 'color': (0,0,0), 'draw': (0,1,2,3)},
  CLASSES_I['LIH'], {'name': 'UNL', 'color': (0,0,0), 'draw': (0,1,2,3)},
}
def draw_glyph(cls, d=1):
  
