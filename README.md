# ML_model_for_TE_annotation

In this project, I have created a deep learning model that can detect DNA transposons and LTR Retrotransposons in any genome file. 

## Motivation

The motivation of this project is the nearly impossible task of mapping transposable elements in different taxa. 
RepeatMasker needs an already created library of transposable elements to map any transposable elements in a particular genome file. 
This is pretty rare when the research is being done on organisms such as butterflies, where nobody has really mapped any transposable elements in those organisms before.

## Tech stack

All the scripts are written in Python. A virtual environment (Anaconda) was used in order to download and install all the libraries needed, since they weren't available in VS Code as extensions. 
I propose that anyone who wants to use the code in this repository do the same or use Linux instead of Windows.


## Getting started

1. Clone the repo:
```bash
git clone https://github.com/yiannisadreas/ML_model_for_TE_annotation.git
cd ML_model_for_TE_annotation
