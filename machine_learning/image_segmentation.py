


from __future__ import print_function
import tensorflow as tf
import skimage.filter as filter
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import sys


image_file = '3.png'

image_string = tf.gfile.FastGFile(image_file).read()

ksize_rows = 299
ksize_cols = 299

# strides_rows and strides_cols determine the distance between
#+ the centers of two consecutive patches.
strides_rows = 299 # 128
strides_cols = 299 # 128

sess = tf.InteractiveSession()

image = tf.image.decode_image(image_string, channels=3)

# The size of sliding window
ksizes = [1, ksize_rows, ksize_cols, 1] 

# How far the centers of 2 consecutive patches are in the image
strides = [1, strides_rows, strides_cols, 1]

rates = [1, 1, 1, 1] # sample pixel consecutively

padding='VALID' # or 'SAME'

image = tf.expand_dims(image, 0)
image_patches = tf.extract_image_patches(image, ksizes, strides, rates, padding)

shapes = tf.shape(image_patches).eval()
nr, nc = shapes[1], shapes[2]

fig = plt.figure()
gs = gridspec.GridSpec(nr, nc)
gs.update(wspace=0.01, hspace=0.01)

for i in range(nr):
  for j in range(nc):
  	x = tf.reshape(image_patches[0,i,j,], [ksize_rows, ksize_cols, 3])
  	x = tf.squeeze(tf.image.rgb_to_grayscale(x))
  	x = sess.run(x)
  	ax = plt.subplot(gs[i*nc+j])
  	plt.axis('off')
  	ax.set_xticklabels([])
  	ax.set_yticklabels([])
  	ax.set_aspect('auto')
  	cutoff = filter.threshold_yen(x)
  	#cutoff = filter.threshold_isodata(x)
  	#cutoff = filter.threshold_otsu(x)
  	plt.imshow(x<cutoff, cmap='gray')
  	cancer_cell_area = (x < cutoff).sum() / float(ksize_rows * ksize_cols)
  	print('processed {},{} patch, {}, cancer cell area: {}.'.format(i,j, i*nc+j, cancer_cell_area),file=sys.stderr)
plt.savefig('image_patches_gray.png', bbox_inches='tight',dpi=120)

#patch1 = tf.reshape(image_patches[0,0,0,], [ksize_rows, ksize_cols, 3])
#patch1_gray = sess.run(tf.squeeze(tf.image.rgb_to_grayscale(patch1)))
#yen_cutoff = filter.threshold_yen(patch1_gray)


