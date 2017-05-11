from __future__ import print_function
import tensorflow as tf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import sys

def plot_image_patches(x, ksize_rows=299, ksize_cols=299):
  nr = x.shape[1]
  nc = x.shape[2]
  # figsize: width and height in inches. can be changed to make
  #+output figure fit well.
  #fig = plt.figure(figsize=(nr, nc)) 
  fig = plt.figure()
  gs = gridspec.GridSpec(nr, nc)
  gs.update(wspace=0.01, hspace=0.01)

  for i in range(nr):
    for j in range(nc):
      ax = plt.subplot(gs[i*nc+j])
      plt.axis('off')
      ax.set_xticklabels([])
      ax.set_yticklabels([])
      ax.set_aspect('auto')
      plt.imshow(x[0,i,j,].reshape(ksize_rows, ksize_cols, 3))
  return fig

def plot_image_patches2(image_patches, sess, ksize_rows=299, ksize_cols=299):
  #x = sess.run(image_patches)
  #nr = x.shape[1]
  #nc = x.shape[2]
  #del x
  a = sess.run(tf.shape(image_patches))
  nr, nc = a[1], a[2]
  print('width: {}; height: {}'.format(nr, nc), file=sys.stderr)
  # figsize: width and height in inches. can be changed to make
  #+output figure fit well. The default often works well.
  #fig = plt.figure(figsize=(nr, nc)) 
  fig = plt.figure()
  gs = gridspec.GridSpec(nr, nc)
  gs.update(wspace=0.01, hspace=0.01)

  for i in range(nr):
    for j in range(nc):
      ax = plt.subplot(gs[i*nc+j])
      plt.axis('off')
      ax.set_xticklabels([])
      ax.set_yticklabels([])
      ax.set_aspect('auto')
      patch = tf.reshape(image_patches[0,i,j,], [ksize_rows, ksize_cols, 3])
      #patch = tf.image.random_brightness(patch, 0.3)
      #patch = tf.image.random_contrast(patch, 0.1, 0.9)
      #patch = tf.image.random_saturation(patch, 0.1, 0.9)
      #patch = tf.image.random_hue(patch, 0.4)
      #patch = tf.image.random_flip_up_down(patch, 0.4)
      plt.imshow(sess.run(patch))
      print('processed {},{} patch, {}.'.format(i,j, i*nc+j),file=sys.stderr)
  return fig



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

# Method 1:
#x=sess.run(image_patches)
#print(x.shape, file=sys.stderr)
#fig = plot_image_patches(x)

# Method 2:
fig = plot_image_patches2(image_patches, sess)

plt.savefig('image_patches.png', bbox_inches='tight',dpi=300) # use dpi to control image size, e.g. 800

plt.close(fig)

sess.close()




