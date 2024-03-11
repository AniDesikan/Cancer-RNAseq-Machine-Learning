# ML_Testing

The aim of this project is to create a backend for an app that can allow researchers to easily classify their data using modern machine learning techniques.

This is a repository showing machine learning testing on various cancer RNA sequencing datasets taken from the Xena UCSC portal, as well as whole slide cancer images. 

The tm_function.py intakes a differential expression csv of enhanced genes, and allows the user to select how many genes they want to input into the machine.

After that, the csv from the tm_function is input into the binary_augmentation.py functions, where the data is bootstrapped and more samples are created using linear interpolation and noise injection, allowing machine learning to be more accurate.

Finally, the binary_machines.py intakes the csv obtained previously (currently just the csv from the tm function) and trains  machine of the user's choice on the rna sequencing data.

This is all tied together in the binary_pipeline functions which are then called in the main.py function, which calls image_pipeline_functions, where a similar process is done for the images.


(Work in progress)
