# Landmark prediction model

Deep-learning pipeline that automatically places anatomical landmarks (keypoints)
on full-body DXA (dual-energy X-ray absorptiometry) images. These landmarks are
the raw input from which all downstream pelvic phenotypes are derived.

The model is an **HRNet** (High-Resolution Network) trained to predict **23
keypoints** across the hip and head regions.

## Files

| File | Type | Description |
| --- | --- | --- |
| `0_train_data_prepare.ipynb` | Notebook | Verify manual-annotation accuracy and build the training set. Converts COCO-format annotations into the 23-keypoint landmark representation for the hip and head regions. |
| `1_model_landmarks_pred.ipynb` | Notebook | Run the trained HRNet model to predict landmarks. Centrally crops each image, predicts on the validation set, and then generates predictions for all images. |
| `hrnet.py` | Module | HRNet model architecture (`hrnet18s`, `hrnet18`, `hrnet32`), adapted from [HRNet-Image-Classification](https://github.com/HRNet/HRNet-Image-Classification). |

## Order

Run `0_train_data_prepare.ipynb` first, then `1_model_landmarks_pred.ipynb`.
The predicted landmarks feed into `../data_analysis` (step 2 onward).

> **Note:** The notebooks were run against UK Biobank imaging data using absolute
> local paths. Paths must be adapted to your own environment, and the underlying
> images/annotations are access-controlled (see the top-level README).
