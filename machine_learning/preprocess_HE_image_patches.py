#!/storage1/lixc/Software/Install/anaconda2-2.4.0/bin/python
# Copyright 2016 The TensorFlow Authors. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# ==============================================================================
r"""Downloads and converts Flowers data to TFRecords of TF-Example protos.

This module downloads the Flowers data, uncompresses it, reads the files
that make up the Flowers data and creates two TFRecord datasets: one for train
and one for test. Each TFRecord dataset is comprised of a set of TF-Example
protocol buffers, each of which contain a single image and label.

The script should take about a minute to run.

"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import argparse
parser = argparse.ArgumentParser(description="Extract image patches of fixed size from H&E whole slide image (WSI), \
  perform augmentation for extracted patches and convert them to tfrecord.")
parser.add_argument('--models_slim_path', dest='models_slim_path', help='The absolute path to `models/slim` \
  where `dataset_utils` module will be loaded. Downloaded from https://github.com/tensorflow/models/tree/master/slim')
parser.add_argument('--image_dir', dest='image_dir', help="The root directory of images")
parser.add_argument('--out_dir', dest='out_dir', help="The output directory")
parser.add_argument('--dataset_name', dest='dataset_name', help='The name of the dataset [flowers]', type=str, default='flowers')
parser.add_argument('--data_label', dest='data_label', help="The data label [train]", type=str, choices=['train','validation','test'], default='train')
parser.add_argument('--num_shards', dest='num_shards', help="The number of shards per dataset split [1]", type=int, default=1)
parser.add_argument('--splitting', dest='splitting', help="If enabled, split data into train and validation sets", action="store_true")
parser.add_argument('--num_validation', dest='num_validation', help="The number of samples in validation set if the `splitting` option is enabled [300]", type=int, default=300)
parser.add_argument('--height', dest='height', help='the height of image patches [299]', type=int, default=299)
parser.add_argument('--width', dest='width', help='the width of image patches [299]', type=int, default=299)
parser.add_argument('--strides_rows', dest='strides_rows', help='vertical distance between the centers of two image patches [299]', type=int, default=299)
parser.add_argument('--strides_cols', dest='strides_cols', help='horizontal distance between the centers of two image patches [299]', type=int, default=299)

args = parser.parse_args()
print(args)

if args.height != args.width:
  raise ValueError('The `height` and `width` must be equal. Input height: {}, width: {}.'.format(args.height, args.width))
if args.strides_rows != args.strides_cols:
  raise ValueError('The `strides_rows` and `strides_cols` must be equal. Input height: {}, width: {}.'.format(args.strides_rows, args.strides_cols))

import math
import os
import random
import sys

if args.models_slim_path is None: 
  sys.path.append('/Users/lixiangchun/GitHub/tensorflow/models/slim/')
else:
  sys.path.append(args.models_slim_path)

import tensorflow as tf

from datasets import dataset_utils

# The URL where the Flowers data can be downloaded.
_DATA_URL = 'http://download.tensorflow.org/example_images/flower_photos.tgz'

# The number of images in the validation set.
_NUM_VALIDATION = args.num_validation

# Seed for repeatability.
_RANDOM_SEED = 0

# The number of shards per dataset split.
_NUM_SHARDS = args.num_shards

dataset_name = args.dataset_name


class ImageReader(object):
  """Helper class that provides TensorFlow image coding utilities."""

  def __init__(self):
    # Initializes function that decodes RGB JPEG data.
    self._decode_jpeg_data = tf.placeholder(dtype=tf.string)
    self._decode_jpeg = tf.image.decode_jpeg(self._decode_jpeg_data, channels=3)

  def read_image_dims(self, sess, image_data):
    image = self.decode_jpeg(sess, image_data)
    return image.shape[0], image.shape[1]

  def decode_jpeg(self, sess, image_data):
    image = sess.run(self._decode_jpeg,
                     feed_dict={self._decode_jpeg_data: image_data})
    assert len(image.shape) == 3
    assert image.shape[2] == 3
    return image


def _get_filenames_and_classes(image_dir):
  """Returns a list of filenames and inferred class names.

  Args:
    dataset_dir: A directory containing a set of subdirectories representing
      class names. Each subdirectory should contain PNG or JPG encoded images.

  Returns:
    A list of image file paths, relative to `dataset_dir` and the list of
    subdirectories, representing class names.
  """
  flower_root = image_dir
  directories = []
  class_names = []
  for filename in os.listdir(flower_root):
    path = os.path.join(flower_root, filename)
    if os.path.isdir(path):
      directories.append(path)
      class_names.append(filename)

  photo_filenames = []
  for directory in directories:
    for filename in os.listdir(directory):
      path = os.path.join(directory, filename)
      photo_filenames.append(path)

  return photo_filenames, sorted(class_names)


def _get_dataset_filename(dataset_dir, split_name, shard_id):
  output_filename = '%s_%s_%05d-of-%05d.tfrecord' % (
      dataset_name, split_name, shard_id, _NUM_SHARDS)
  return os.path.join(dataset_dir, output_filename)

def _encode_image(image):
  return tf.image.encode_jpg(image)



def _rotate_flip_core(sess, patch_image, rot90_k, tfrecord_writer, height, width, class_id):
  k = 0
  rotated_image = tf.image.rot90(patch_image, k=rot90_k)
  rotated_image_data = tf.image.encode_jpeg(rotated_image)
  example = dataset_utils.image_to_tfexample(sess.run(rotated_image_data), 'jpg', height, width, class_id)
  tfrecord_writer.write(example.SerializeToString())
  k += 1

  flipped_image = tf.image.flip_left_right(patch_image)
  flipped_image_data = tf.image.encode_jpeg(flipped_image)
  example = dataset_utils.image_to_tfexample(sess.run(flipped_image_data), 'jpg', height, width, class_id)
  tfrecord_writer.write(example.SerializeToString())
  k += 1

  return k

def _rotate_flip(sess, patch_image, tfrecord_writer, height, width, class_id):
  k1 = _rotate_flip_core(sess, patch_image, 1, tfrecord_writer, height, width, class_id)
  k2 = _rotate_flip_core(sess, patch_image, 2, tfrecord_writer, height, width, class_id)
  k3 = _rotate_flip_core(sess, patch_image, 3, tfrecord_writer, height, width, class_id)
  k = k1 + k2 + k3
  return k

def _image_random_X(sess, patch_image, tfrecord_writer, height, width, class_id):
  k = 0
  ## random_brightness
  patch_image_data = tf.image.encode_jpeg(tf.image.random_brightness(patch_image, max_delta=64/255))
  example = dataset_utils.image_to_tfexample(sess.run(patch_image_data), 'jpg', height, width, class_id)
  tfrecord_writer.write(example.SerializeToString())
  k += 1

  ## random_saturation
  patch_image_data = tf.image.encode_jpeg(tf.image.random_saturation(patch_image, lower=0, upper=0.25))
  example = dataset_utils.image_to_tfexample(sess.run(patch_image_data), 'jpg', height, width, class_id)
  tfrecord_writer.write(example.SerializeToString())
  k += 1

  ## random_hue
  atch_image_data = tf.image.encode_jpeg(tf.image.random_hue(patch_image, max_delta=0.04))
  example = dataset_utils.image_to_tfexample(sess.run(patch_image_data), 'jpg', height, width, class_id)
  tfrecord_writer.write(example.SerializeToString())
  k += 1

  ## random_constrast
  patch_image_data = tf.image.encode_jpeg(tf.image.random_contrast(patch_image, lower=0, upper=0.75))
  example = dataset_utils.image_to_tfexample(sess.run(patch_image_data), 'jpg', height, width, class_id)
  tfrecord_writer.write(example.SerializeToString())
  k += 1
  
  return k

def _convert_dataset(split_name, filenames, class_names_to_ids, dataset_dir, ksize_rows=299, ksize_cols=299, strides_rows=128, strides_cols=128):
  """Converts the given filenames to a TFRecord dataset.

  Args:
    split_name: The name of the dataset, either 'train' or 'validation'.
    filenames: A list of absolute paths to png or jpg images.
    class_names_to_ids: A dictionary from class names (strings) to ids
      (integers).
    dataset_dir: The directory where the converted datasets are stored.
    ksize_rows, ksize_cols: the height and width of extracted image patches
    strides_rows, strides_cols: the distance between the centers of two consecutive patches
  """
  assert split_name in ['train', 'validation']

  num_per_shard = int(math.ceil(len(filenames) / float(_NUM_SHARDS)))
  number_of_images = 0

  # The size of sliding window
  ksizes = [1, ksize_rows, ksize_cols, 1]

  # How far the centers of 2 consecutive patches are in the image
  strides = [1, strides_rows, strides_cols, 1]

  rates = [1, 1, 1, 1] # sample pixel consecutively

  padding = 'VALID' # or 'SAME'

  total_image_number = len(filenames)
  d = dict.fromkeys(class_names_to_ids.values())

  with tf.Graph().as_default():

    with tf.Session('') as sess:

      for shard_id in range(_NUM_SHARDS):
        output_filename = _get_dataset_filename(
            dataset_dir, split_name, shard_id)

        with tf.python_io.TFRecordWriter(output_filename) as tfrecord_writer:
          start_ndx = shard_id * num_per_shard
          end_ndx = min((shard_id+1) * num_per_shard, len(filenames))
          for i in range(start_ndx, end_ndx):

            try:
              # Read the filename:
              image_data = tf.gfile.FastGFile(filenames[i], 'r').read()
            except Exception as e:
              sys.stderr.write("Error in decoding image {} into tensor - {}.".format(filenames[i], str(e)))
              continue

            image = tf.image.decode_image(image_data, channels=3)
            image = tf.expand_dims(image, 0)

            image_patches = tf.extract_image_patches(image, ksizes, strides, rates, padding)
            image_patch_shape = sess.run(tf.shape(image_patches))
            nrows, ncols = image_patch_shape[1], image_patch_shape[2]
            #print('{},{}'.format(nrows,ncols), file=sys.stderr)

            class_name = os.path.basename(os.path.dirname(filenames[i]))
            class_id = class_names_to_ids[class_name]

            for nr in range(nrows):
              for nc in range(ncols):
                patch_image = tf.reshape(image_patches[0,nr,nc,], [ksize_rows, ksize_cols, 3])
                height, width = ksize_rows, ksize_cols

                k = 0
                
				        # original image patch
                patch_image_data = tf.image.encode_jpeg(patch_image)
                example = dataset_utils.image_to_tfexample(sess.run(patch_image_data), 'jpg', height, width, class_id)
                tfrecord_writer.write(example.SerializeToString())
                k += 1

                flipped_image = tf.image.flip_left_right(patch_image)
                flipped_image_data = tf.image.encode_jpeg(flipped_image)
                example = dataset_utils.image_to_tfexample(sess.run(flipped_image_data), 'jpg', height, width, class_id)
                tfrecord_writer.write(example.SerializeToString())
                k += 1

                k += _rotate_flip(sess, patch_image, tfrecord_writer, height, width, class_id)
                k += _image_random_X(sess, patch_image, tfrecord_writer, height, width, class_id)

                number_of_images += k
                d[class_id] += k

                sys.stdout.write('\r>> Converting image {}/{} shard {}, patch {}/{},{}/{}, total patch: {}.'.format(i+1, total_image_number, shard_id, nr, nrows, nc, ncols, number_of_images))
                sys.stdout.flush()

            #sys.stdout.write('\r>> Converting image %d/%d shard %d' % (i+1, len(filenames), shard_id))
            #sys.stdout.flush()

  sys.stdout.write('\n')
  sys.stdout.flush()

  with open('{}/class_id_number.txt'.format(dataset_dir), 'w') as f:
    for k, v in d.iteritems():
      print('{}:{}'.format(k, v), file=f)

  with open('{}/{}_number_of_images.txt'.format(dataset_dir, split_name), 'w') as f:
    f.write(str(number_of_images))

def _clean_up_temporary_files(dataset_dir):
  """Removes temporary files used to create the dataset.

  Args:
    dataset_dir: The directory where the temporary files are stored.
  """
  filename = _DATA_URL.split('/')[-1]
  filepath = os.path.join(dataset_dir, filename)
  tf.gfile.Remove(filepath)

  tmp_dir = os.path.join(dataset_dir, 'flower_photos')
  tf.gfile.DeleteRecursively(tmp_dir)


def _dataset_exists(dataset_dir):
  for split_name in ['train', 'validation']:
    for shard_id in range(_NUM_SHARDS):
      output_filename = _get_dataset_filename(
          dataset_dir, split_name, shard_id)
      if not tf.gfile.Exists(output_filename):
        return False
  return True


def run(image_dir, out_dir = 'tfrecord_out', data_label = 'train', ksize_rows=299, ksize_cols=299, strides_rows=128, strides_cols=128):
  """Runs the conversion operation.

  Args:
    image_dir: The dataset directory where the dataset is stored.
    out_dir: The output directory
    data_label: data label
    ksize_rows, ksize_cols: the height and width of image patches to extract
    strides_rows, strides_cols: the distance between two consecutive image patches
  """
  if not tf.gfile.Exists(out_dir):
    tf.gfile.MakeDirs(out_dir)
  photo_filenames, class_names = _get_filenames_and_classes(image_dir)
  class_names_to_ids = dict(zip(class_names, range(len(class_names))))

  # Divide into train and test:
  random.seed(_RANDOM_SEED)
  random.shuffle(photo_filenames)

  if args.splitting:
    training_filenames = photo_filenames[_NUM_VALIDATION:]
    validation_filenames = photo_filenames[:_NUM_VALIDATION]

    # First, convert the training and validation sets.
    _convert_dataset('train', training_filenames, class_names_to_ids, out_dir, ksize_rows, ksize_cols, strides_rows, strides_cols)
    _convert_dataset('validation', validation_filenames, class_names_to_ids, out_dir, ksize_rows, ksize_cols, strides_rows, strides_cols)
  else:
    _convert_dataset(data_label, photo_filenames, class_names_to_ids, out_dir, ksize_rows, ksize_cols, strides_rows, strides_cols)

  # Finally, write the labels file:
  labels_to_class_names = dict(zip(range(len(class_names)), class_names))
  dataset_utils.write_label_file(labels_to_class_names, out_dir)

  #_clean_up_temporary_files(dataset_dir)
  print('\nFinished converting the dataset!')

## modified from models/slim/datasets/download_and_convert_flowers.py
run(args.image_dir, args.out_dir, args.data_label, args.height, args.width, args.strides_rows, args.strides_cols)
