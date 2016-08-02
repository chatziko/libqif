# Semantic Example

This example is useful for studying how priors constructed using the knowledge of semantic labels from 
different location based datasets (eg. foursquare) can affect the adversary's expected loss and in turn the user's privacy.

## Using this example

Installation of libqif is necessary. Details regarding this can be found in the README of libqif.

* `Get datasets`: Download the relevent location based datasets. Please refer readme in [datasets](https://github.com/susheels/libqif/tree/master/samples/semantic/preproc_foursquare/datasets) folder.
* `Preprocessing`: To be able to use libqif for the study,one has to constuct global and user priors from the dataset. Use the python code in [preproc_foursquare](https://github.com/susheels/libqif/tree/master/samples/semantic/preproc_foursquare). Details regarding the usage is in it's readme. Redo this procedure for different cases of labels to get it's specific global and user prior.
* `Study`: Use the c++ file test_foursquare.cpp in [semantic](https://github.com/susheels/libqif/tree/master/samples/semantic) which has the complete workflow of generating laplace mechanisms, strategies and finnaly the privacy metric. 
* `Visualize`: Use the graph_foursquare.py in [semantic](https://github.com/susheels/libqif/tree/master/samples/semantic) for generating boxplots.



