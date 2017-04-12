
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


