# Extract image patches in `tensorflow`

When applying deep learning to classify gigapixel pathology images, we often apply trained models (e.g. `Resnet` and `inception` etc.) to every image patches. The prediction for all patches of each image forms a heatmap.

Here in this example, I try to briefly present an example to extract all image patches from a pathology image.

```python
import tensorflow as tf

image_file = '996.jpg'
image_string = tf.gfile.FastGFile(image_file).read()

ksize_rows = 299
ksize_cols = 299
strides_rows = 128
strides_cols = 128

sess = tf.InteractiveSession()

image = tf.image.decode_image(image_string, channels=3)

# The size of sliding window
ksizes = [1, ksize_rows, ksize_cols, 1] 

# How far the centers of 2 consecutive patches are in the image
strides = [1, strides_rows, strides_cols, 1]

# The document is unclear. However, an intuitive example posted on StackOverflow illustrate its behaviour clearly. 
# http://stackoverflow.com/questions/40731433/understanding-tf-extract-image-patches-for-extracting-patches-from-an-image
rates = [1, 1, 1, 1] # sample pixel consecutively

# padding algorithm to used
padding='SAME'

image = tf.expand_dims(image, 0)
image_patches = tf.extract_image_patches(image, ksizes, strides, rates, padding)

# print image shape of image patches
print sess.run(image_patches).shape

# retrieve the 1st patches
patch1 = image_patches[0,0,0,]

# reshape
patch1 = tf.reshape(patch1, [ksize_rows, ksize_cols, 3])

# visualize image
import matplotlib.pyplot as plt
plt.imshow(sess.run(patch1))
plt.show()

# close session
sess.close()
```

The extracted image patches can be fed to fine-tuned model such as `inception_v3`. To feed image patches to fine-tuned `ResNet` model, change `ksize_rows` and `ksize_cols` to 224, which is the image size used by `ResNet`, 

## An awesome [example](http://stackoverflow.com/questions/40731433/understanding-tf-extract-image-patches-for-extracting-patches-from-an-image) from stack overflow to understand `tf.extract_image_patches`. 

Here is how the method works:

- `ksize` is used to decide the dimensions of each patch, or in other words, how many pixels each patch should contain.

- `strides` denotes the length of the gap between the start of one patch and the start of the next consecutive patch within the original image.

- `rates` is a number that essentially means our patch should jump by rates pixels in the original image for each consecutive pixel that ends up in our patch. (The example below helps illustrate this.)

- `padding` is either "VALID", which means every patch must be fully contained in the image, or "SAME", which means patches are allowed to be incomplete (the remaining pixels will be filled in with zeroes).

Here is some sample code with output to help demonstrate how it works:

```python
import tensorflow as tf

n = 10
# images is a 1 x 10 x 10 x 1 array that contains the numbers 1 through 100 in order
images = [[[[x * n + y + 1] for y in range(n)] for x in range(n)]]

# We generate four outputs as follows:
# 1. 3x3 patches with stride length 5
# 2. Same as above, but the rate is increased to 2
# 3. 4x4 patches with stride length 7; only one patch should be generated
# 4. Same as above, but with padding set to 'SAME'
with tf.Session() as sess:
  print tf.extract_image_patches(images=images, ksizes=[1, 3, 3, 1], strides=[1, 5, 5, 1], rates=[1, 1, 1, 1], padding='VALID').eval(), '\n\n'
  print tf.extract_image_patches(images=images, ksizes=[1, 3, 3, 1], strides=[1, 5, 5, 1], rates=[1, 2, 2, 1], padding='VALID').eval(), '\n\n'
  print tf.extract_image_patches(images=images, ksizes=[1, 4, 4, 1], strides=[1, 7, 7, 1], rates=[1, 1, 1, 1], padding='VALID').eval(), '\n\n'
  print tf.extract_image_patches(images=images, ksizes=[1, 4, 4, 1], strides=[1, 7, 7, 1], rates=[1, 1, 1, 1], padding='SAME').eval()
```

The content of `images`
```python
In [183]: images
Out[183]: 
[[[[1], [2], [3], [4], [5], [6], [7], [8], [9], [10]],
  [[11], [12], [13], [14], [15], [16], [17], [18], [19], [20]],
  [[21], [22], [23], [24], [25], [26], [27], [28], [29], [30]],
  [[31], [32], [33], [34], [35], [36], [37], [38], [39], [40]],
  [[41], [42], [43], [44], [45], [46], [47], [48], [49], [50]],
  [[51], [52], [53], [54], [55], [56], [57], [58], [59], [60]],
  [[61], [62], [63], [64], [65], [66], [67], [68], [69], [70]],
  [[71], [72], [73], [74], [75], [76], [77], [78], [79], [80]],
  [[81], [82], [83], [84], [85], [86], [87], [88], [89], [90]],
  [[91], [92], [93], [94], [95], [96], [97], [98], [99], [100]]]]
```

Output from the 1st `print`:
```python
print tf.extract_image_patches(images=images, ksizes=[1, 3, 3, 1], strides=[1, 5, 5, 1], rates=[1, 1, 1, 1], padding='VALID').eval()

[[[[ 1  2  3 11 12 13 21 22 23]
   [ 6  7  8 16 17 18 26 27 28]]

  [[51 52 53 61 62 63 71 72 73]
   [56 57 58 66 67 68 76 77 78]]]]
```

Here `[ 1  2  3 11 12 13 21 22 23]` is the 1st patch from `images`, e.g.:
```python
[1], [2], [3]
[11], [12], [13]
[21], [22], [23]
```

Output from the 2nd `print`:
```python
  print tf.extract_image_patches(images=images, ksizes=[1, 3, 3, 1], strides=[1, 5, 5, 1], rates=[1, 2, 2, 1], padding='VALID').eval()

[[[[  1   3   5  21  23  25  41  43  45]
   [  6   8  10  26  28  30  46  48  50]]

  [[ 51  53  55  71  73  75  91  93  95]
   [ 56  58  60  76  78  80  96  98 100]]]] 
```

Here `[  1   3   5  21  23  25  41  43  45]` is reshaped as the 1st patch. By comparing with the content of `images`, you  will understand what `rates` does.
```python
[1], [3], [5]
[21], [23], [25]
[41], [43], [45]
```



