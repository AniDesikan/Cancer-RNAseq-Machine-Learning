# ML_Testing

This is a repository showing machine learning testing on various cancer RNA sequencing datasets taken from the Xena UCSC portal. The tm_function.py intakes a differential expression csv of enhanced genes, and allows the user to select for only those genes
  * In the future, the user will be able to select for genes of pvalue less than a number they desire

After that, the csv from the tm_function is input into the binary_augmentation.py functions, where the data is bootstrapped and more samples are created, allowing machine learning to be more accurate

Finally, the binary_machines.py intakes the csv obtained previously (currently just the csv from the tm function) and trains an SVM machine on the rna sequencing data. Currently we are able to obtain 99% accuracy in determine whether a sample is a tumor or not from the machine.

(Work in progress)
