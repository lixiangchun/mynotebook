
## classify_image.py
Make prediction for multiple jpeg images or images stored in TFRecord file using pretrained tensorflow models, for example, `inception_v3` and `resnet_v1` etc.

```python
mynotebook/machine_learning$ python classify_image.py -h
usage: classify_image.py [-h] [--num_classes NUM_CLASSES] [--infile INFILE]
                         [--tfrecord [TFRECORD]] [--notfrecord]
                         [--outfile OUTFILE] [--model_name MODEL_NAME]
                         [--preprocessing_name PREPROCESSING_NAME]
                         [--checkpoint_path CHECKPOINT_PATH]
                         [--eval_image_size EVAL_IMAGE_SIZE]

optional arguments:
  -h, --help            show this help message and exit
  --num_classes NUM_CLASSES
                        The number of classes.
  --infile INFILE       Image file, one image per line.
  --tfrecord [TFRECORD]
                        Input file is formatted as TFRecord.
  --notfrecord
  --outfile OUTFILE     Output file for prediction probabilities.
  --model_name MODEL_NAME
                        The name of the architecture to evaluate.
  --preprocessing_name PREPROCESSING_NAME
                        The name of the preprocessing to use. If left as
                        `None`, then the model_name flag is used.
  --checkpoint_path CHECKPOINT_PATH
                        The directory where the model was written to or an
                        absolute path to a checkpoint file.
  --eval_image_size EVAL_IMAGE_SIZE
                        Eval image size.

# Predict jpeg images, input file `a.txt` includes path to images.
mynotebook/machine_learning$ python classify_image.py --num_classes 1001 --infile a.txt --model_name inception_v3 --checkpoint_path ../../tensorflow/models/slim/pretrained_models/inception_v3.ckpt --outfile out.txt

# Predict for TFRecord input
mynotebook/machine_learning$ python classify_image.py --num_classes 1001 --infile ../../tensorflow/models/slim/flowers/flowers_train_00000-of-00005.tfrecord --model_name inception_v3 --checkpoint_path ../../tensorflow/models/slim/pretrained_models/inception_v3.ckpt --outfile out.txt --tfrecord
```

## Extract learned features from pretrained model
The pretrained deep learning model can be used to extract learned feature for an input image. Take `inception_v3` as an example, the following code is able to extract bottleneck features from pretrained`inception_v3` model.

What we need to do is to replace a  few lines of code in `classify_image.py`.


```python
#logits, _ = network_fn(processed_images)
logits, end_points = network_fn(processed_images)

# The keys of end_points, i.e. features that can be extracted.
print(end_points.keys())

probabilities = tf.nn.softmax(logits)

# 'PreLogits' is the bottleneck feature
PreLogits = tf.squeeze(end_points['PreLogits'])
```

```python
#probs = sess.run(probabilities, feed_dict={image_string:x})
probs, prelogits = sess.run([probabilities, PreLogits], feed_dict={image_string:x})

print(prelogits.shape) # the dimension of bottleneck features
print(prelogits) # print bottleneck feature
```

In summary, extract bottleneck features within `slim` is very convenient. Caution to be taken is that different models have different bottleneck features. If you are using `ResNet` as a feature extractor, you need to modify `end_points['PreLogits']` accordingly.