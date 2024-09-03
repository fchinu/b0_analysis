# Macros for BDT training + application

## Requirements
In order to execute the training macro in this folder, the following python libraries are needed:
- [hipe4ml](https://github.com/hipe4ml/hipe4ml)
- [hipe4ml_converter](https://github.com/hipe4ml/hipe4ml_converter)

## Training

### Samples
Run the preliminary workflow on **MC** data on Hyperloop and download the derived AOD. Then, run the reduced workflow using this derived AOD as input. This will produce new derived AOD tables that we will use as input for our BDT training macro.

*Note 1: we run the reduced workflow with ML activated for D mesons only, not for B0.*

*Note 2: the sample preparation is done inside the (reduced) task as the signal is tagged with `fOriginMcRec=1` and the background with `fOriginMcRec=0`.*

### Perform training
In order to perform the training and produce the BDT models to be used on real data, the following script shall be used:
```
python3 ml_training_xor_application.py --config config_ml_B0ToDPi.yml --train
```
Given the output directory set in `config_ml_B0ToDPi.yml`, a directory is created for each pT bin (i.e. each model trained) and filled with:
- plots at data preparation level: variables distributions and correlations
- plots at training-testing level: BDT output scores, ROC curves, precision recall, features importance
- trained models in files containing `ModelHandler` prefix (in `pickle` and `onnx` formats). *The `.onnx` can be used for ML inference in O2Physics selection task.*
- model applied to test set in file containing `ModelApplied` suffix

## Application

### Samples
Run the preliminary workflow on **real** data on Hyperloop and download the derived AOD. Then, run the reduced workflow using this derived AOD as input. This will produce new derived AOD tables that we will use as input for our BDT application macro.

*Note: we run the reduced workflow with ML activated for D mesons only, not for B0.*

### Perform application
In order to apply BDT models to real data and produce a `.root` file, the following script shall be used:
```
python3 ml_training_xor_application.py --config config_ml_B0ToDPi.yml --apply
```

*Note: in order to perform KDE fit, one also needs to apply BDT models on MC data with `checkDecayTypeMc` activated.*
