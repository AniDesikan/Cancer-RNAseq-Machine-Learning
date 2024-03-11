#################################################
# IMPORTS
#################################################
PARAM_GRID = {
        'angle': [0, 90, 180, 270],
        'flip_code': [-1, 0, 1],
        'crop_size': [50, 100, 150, 200, 250, 300, 350, 400, 450, 500],
        'shear_factor': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        'rotate_prop': [0.1, 0.2, 0.3],
        'flip_prop': [0.1, 0.2, 0.3],
        'crop_scale_prop': [0.1, 0.2, 0.3],
        'color_number': [0, 1, 2, 3, 4, 5],
        'kernel_size': [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]
    }
BATCH_SIZE = 200
import gc
import cv2
import itertools
import random
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score, f1_score, roc_auc_score, roc_curve, auc
import pandas as pd
import numpy as np
import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import glob
from PIL import Image
from sklearn.metrics import accuracy_score
from machines.Xception import Xception, BasicConv, ResidualUnit, MiddleFlowUnit, ExitFlow
import tensorflow as tf
from tensorflow.keras import models,layers,losses,metrics,optimizers,Input
# gpus = tf.config.experimental.list_physical_devices('GPU')
# print(gpus)
# if gpus:
#     for gpu in gpus:
#         tf.config.experimental.set_memory_growth(gpu, True)
from machines.Xception import Xception
from machines.ResNet import ResNet_152
from machines.DenseNet import DenseNet_264
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.metrics import accuracy_score
from data_augmentation.image_processor_all import ImageProcessor,ImageProcessor2,augment_images_with_proportions_2,augment_images_with_proportions,combined_augment_function
from sklearn.metrics import roc_auc_score, precision_score, recall_score, f1_score
from features.image_features import *
# from machines.Machines import Machine

#################################################
# START OF FUNCTIONS
#################################################

import os
print(os.getenv("CUDA_VISIBLE_DEVICES"))

# opt = tf.keras.mixed_precision.LossScaleOptimizer(opt)

# policy = mixed_precision.Policy('mixed_float16')
# mixed_precision.set_policy(policy)

#################################################
# ADAPTED IMAGE AUGMENTATION FUNCTIONS
#################################################
 

def remove_images_from_folder(folder_path):
    # List all files in the folder
    files = os.listdir(folder_path)

    # Iterate over each file in the folder
    for file in files:
        file_path = os.path.join(folder_path, file)
        
        # Check if the file is a regular file and ends with an image extension
        if os.path.isfile(file_path) and file.endswith(('.jpg', '.jpeg', '.png', '.gif')):
            os.remove(file_path)  # Remove the file

    print(f"All images removed from the folder: {folder_path}")

#################################################
# AUGMENTATION FUNCTIONS
#################################################

#This function gives a dictionary of random parameters from the PARAM_GRID above
# "Random"
def get_augmentation_parameters(random_params, user_params=None):
    param_combinations = [dict(zip(PARAM_GRID.keys(), values)) for values in itertools.product(*PARAM_GRID.values())]
    print("param_combinations starting...")
    for param_dict in param_combinations:
        param_dict['segment_prop'] = 1 - (param_dict['rotate_prop'] + param_dict['flip_prop'] + param_dict['crop_scale_prop'])
    print("param_combinations finished!")
    if random_params:
        param_dict = random.choice(param_combinations)
        print("Selected random param_dict:")
        print(param_dict)
    else:
        if user_params == None:
            param_dict = param_combinations[0]
            print("Selected first param_dict:")
            print(param_dict)
        else:
            param_dict = {}
            for key, value in zip(PARAM_GRID.keys(), user_params):
                param_dict[key] = value
            print("Selected user-defined param_dict:")
            print(param_dict)
    return param_dict

# immediately make a param dict for use in augmentation so that it's usable for both normal and tumor datasets

# random_param_dict_selected = get_augmentation_parameters(True)

#################################################
# IMAGE DATAFRAME CONSTRUCTION FUNCTIONS
#################################################

# Function to load images and create a DataFrame, used when you need to load images into a dataframe without augmenting them

def put_images_in_df(images, tumor_label):
    image_arrays = []
    for image_path in images:
        with Image.open(image_path) as img:
            img_array = np.array(img)
            # Expand dimensions to make it compatible with conv2d
            img_array = np.expand_dims(img_array, axis=0)
            image_arrays.append(img_array)
    data = {'Image_Array': image_arrays, 'Tumor': [tumor_label] * len(images)}
    df = pd.DataFrame(data)
    
    return df

# Function that you use when you need to load augmented image arrays into a pandas dataframe

def create_image_dataframe(image_arrays, tumor_label):
    # Convert image arrays to a list of dictionaries
    data = [{'Image_Array': np.expand_dims(img_array, axis=0), 'Tumor': tumor_label} for img_array in image_arrays]
    # Create a DataFrame from the list of dictionaries
    df = pd.DataFrame(data)
    return df

# This function uses put images in df if you aren't using augmentation, and uses create_image_dataframe if you are

#This applies the same augmentation strategies to both normal and cancer datasets, old version did not
def make_image_df(folder, tumor_label, number_images):
    images = glob.glob(folder + "/*.jpg")
    df = put_images_in_df(images[:number_images], tumor_label)
    return df

#################################################
# MACHINE MODELLING FUNCTIONS
#################################################

#lr & early stop
def my_lr(epoch,init_lr=0.001,attenuation_rate=0.94,attenuation_step=2):
    lr = init_lr
    lr = lr * attenuation_rate ** (epoch // attenuation_step)
    lr = max(2e-06, lr)
    return lr

def run_model(model, train_ds, test_ds, y_val, epochs):
    train_count = len(train_ds)
    test_count = len(test_ds)
    val_steps = test_count//BATCH_SIZE
    steps_per_epoch = train_count//BATCH_SIZE

    if model == "Xception":
        model = Xception()
        model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
                loss=tf.keras.losses.binary_crossentropy,
                metrics=['acc'])
        # model = Machine(model_name='xception', learning_rate=0.001)
        
    elif model == "Densenet":
        input_shape = (512, 512, 3)  
        inputs = tf.keras.layers.Input(shape=input_shape)
        outputs = DenseNet_264(inputs)
        model = models.Model(inputs=inputs, outputs=outputs)
        model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
                loss=tf.keras.losses.binary_crossentropy,
                metrics=['acc'])
        # model = Machine(model_name='densenet', learning_rate=0.001)
    elif model == "Resnet":
        model = ResNet_152()
        model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.001),
              loss=tf.keras.losses.binary_crossentropy,
              metrics=['acc'])
        # model = Machine(model_name='resnet', learning_rate=0.001)

    lr_callback = tf.keras.callbacks.LearningRateScheduler(my_lr, verbose=False)
    es_callback = tf.keras.callbacks.EarlyStopping(
        monitor='val_acc', patience=15, mode='max', restore_best_weights=True)
    

    history = model.fit(train_ds, epochs=epochs,
                    steps_per_epoch=steps_per_epoch,
                    validation_data=test_ds,
                    validation_steps=val_steps,
                    callbacks=[lr_callback,es_callback])

    loss, accuracy = model.evaluate(test_ds, steps=val_steps)
    # print(f"Loss: {loss}, Accuracy: {accuracy}")
    y_pred = model.predict(test_ds)
    y_pred_binary = (y_pred > 0.5).astype(int)
    precision = precision_score(y_val, y_pred_binary)
    recall = recall_score(y_val, y_pred_binary)
    f1 = f1_score(y_val, y_pred_binary)
    roc_auc = roc_auc_score(y_val, y_pred)
    print({'loss': loss, 'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f1': f1, 'roc_auc': roc_auc})
    # print(f"Precision: {precision}, Recall: {recall}, F1 Score: {f1}, ROC-AUC: {roc_auc}")
    return {'loss': loss, 'accuracy': accuracy, 'precision': precision, 'recall': recall, 'f1': f1, 'roc_auc': roc_auc}

#################################################
# END OF FUNCTIONS
#################################################

def count_images_in_folder(folder_path):
    # List files in the folder
    files = os.listdir(folder_path)

    # Filter only image files (assuming extensions like .jpg, .png, etc.)
    image_files = [f for f in files if f.endswith(('.jpg', '.jpeg', '.png', '.gif'))]

    # Count the number of image files
    num_images = len(image_files)
    return num_images

def create_dataset(X, y, batch_size):
    num_batches = len(X) // batch_size
    initial_ds = tf.data.Dataset.from_tensor_slices((X[:batch_size], y[:batch_size]))

    for i in range(1, num_batches):
        next_ds = tf.data.Dataset.from_tensor_slices((X[i*batch_size:(i+1)*batch_size], 
                                                       y[i*batch_size:(i+1)*batch_size]))
        initial_ds = initial_ds.concatenate(next_ds)

    if len(X) % batch_size != 0:
        remaining_ds = tf.data.Dataset.from_tensor_slices((X[num_batches*batch_size:], 
                                                            y[num_batches*batch_size:]))
        initial_ds = initial_ds.concatenate(remaining_ds)

    return initial_ds

def create_multiple_datasets(X, y, batch_size, num_slices):
    datasets = []
    for i in range(num_slices):
        start_idx = i * batch_size
        end_idx = (i + 1) * batch_size
        dataset = create_dataset(X[start_idx:end_idx], y[start_idx:end_idx], batch_size)
        datasets.append(dataset)
    
    # Concatenate all datasets
    final_dataset = datasets[0]
    for dataset in datasets[1:]:
        final_dataset = final_dataset.concatenate(dataset)
    
    return final_dataset

# Here folder 1 has to be the normal samples, and folder 2 has to be the cancer samples
def pipeline_image_df(folder1, folder2, augmentation, num_images):
    if num_images == None:
        num_images_normal = count_images_in_folder(folder1)
        num_images_scc = count_images_in_folder(folder2)
        num_images = min(num_images_normal,num_images_scc)
    if augmentation == True:
        params = get_augmentation_parameters(True)
        folder1_augmented = folder1 + "_augmented"
        folder2_augmented = folder2 + "_augmented"

        if not os.path.exists(folder1_augmented):
            os.makedirs(folder1_augmented)
        else:
            remove_images_from_folder(folder1_augmented)
        if not os.path.exists(folder2_augmented):
            os.makedirs(folder2_augmented)
        else:  
            remove_images_from_folder(folder2_augmented)
        
        
        print("normal_df Augmentation")
        print("augment function")
        combined_augment_function(folder1, folder1_augmented, **params)
        num_images = count_images_in_folder(folder1_augmented)
        print(f"Number of images in normal augmented folder: {num_images}")

        print("scc_df Augmentation")
        print("augment function")
        combined_augment_function(folder2, folder2_augmented, **params)
        num_images = count_images_in_folder(folder2_augmented)
        print(f"Number of images in scc augmented folder: {num_images}")

        normal_df = make_image_df(folder1_augmented, number_images=num_images * 4, tumor_label=0)
        scc_df = make_image_df(folder2_augmented, number_images=num_images * 4, tumor_label=1)
        return normal_df, scc_df, params
    else:
        normal_df = make_image_df(folder1, number_images=num_images, tumor_label=0)
        scc_df = make_image_df(folder2, number_images=num_images, tumor_label=1)
        params = {'angle': False, 'flip_code': False, 'crop_size': False, 'shear_factor': False, 'rotate_prop':False, 'flip_prop':False, 'crop_scale_prop': False, 'color_number': False, 'kernel_size': False, 'segment_prop': False}
        return normal_df, scc_df, params

def pipeline_features(normal_df, scc_df, participation_factors, params = None):
    normal_images = normal_df['Image_Array']
    print("Starting image features")
    normal_features_list = Image_Features.generate_image_features(normal_images, participation_factors, None)
    print("Image features done!")
    normal_features_df = pd.DataFrame(normal_features_list)
    normal_df = pd.concat([normal_df, normal_features_df], axis=1)
    scc_images = scc_df['Image_Array']
    scc_features_list = Image_Features.generate_image_features(scc_images, participation_factors, None)
    scc_features_df = pd.DataFrame(scc_features_list)
    scc_df = pd.concat([scc_df, scc_features_df], axis=1)
    params.update(participation_factors)
    return normal_df, scc_df, params

def pipeline_machine(normal_df, scc_df, model, params = None, epochs = 20):
    results = []
    oral_image_data_augmented = pd.concat([normal_df, scc_df], ignore_index=True)
    X = np.asarray(oral_image_data_augmented['Image_Array'].tolist()).astype(np.float32)
    y = np.asarray(oral_image_data_augmented['Tumor']).astype(int)
    print("X Length: " + str(len(X)))
    print("y Length: " + str(len(y)))
    print("split")
    X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.2, random_state=42)
    y_train = np.asarray(y_train).astype('float32').reshape((-1, 1))
    y_val = np.asarray(y_val).astype('float32').reshape((-1, 1))
    print("X_train Length: " + str(len(X_train)))
    print("X_val Length: " + str(len(X_val)))
    print("tensor_slices")
    with tf.device("CPU"):
        # train_ds = (X_train, y_train, 32)
        # test_ds = create_dataset(X_val, y_val, 32)
        train_data = (X_train, y_train)
        train_ds = tf.data.Dataset.from_tensor_slices(train_data)
        test_data = (X_val, y_val)
        test_ds = tf.data.Dataset.from_tensor_slices(test_data)
        # train_ds = tf.data.Dataset.from_tensor_slices(X_train, y_train)
        # test_ds = tf.data.Dataset.from_tensor_slices(X_val, y_val)
    print("dataset lengths")
    print(len(train_ds))
    print(len(test_ds))

    metrics = run_model(model, train_ds, test_ds, y_val, epochs=epochs)
    results.append({'params': params, 'metrics': metrics})
    return metrics
