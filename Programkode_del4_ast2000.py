import numpy as np
from PIL import Image

img = Image.open('sample0000.png')
pixels = np.array(img)
width = len(pixels[0, :])
redpixs = [(255, 0, 0) for i in range(width)]
# img2 = Image.fromarray(pixels)
print(img)
# img.save('sample0000_withredline.png')
