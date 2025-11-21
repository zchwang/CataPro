# CataPro
<img src="models/logo.png">

Predicting enzyme kinetic parameters is a crucial task in enzyme discovery and enzyme engineering. Here, we propose a new enzyme kinetic parameter prediction algorithm called CataPro, based on protein language models, small molecule language models, and molecular fingerprints. We collected the latest turnover number (kcat), Michaelis constant (Km), and catalytic efficiency (kcat/Km) data from the BRENDA and SABIO-RK databases. By clustering these data based on 0.4 protein sequence similarity, we obtained the corresponding 10-fold cross-validation datasets. CataPro was trained on these unbiased 10-fold cross-validation datasets, demonstrating superior performance compared to previous predictors in predicting kcat, Km, and kcat/Km.

<img src="models/catapro.png">

---

## Installation 

## Setup a Python environment 

To ensure a clean and isolated setup, we recommend to use [uv](https://docs.astral.sh/uv/), a lightweight tool that simplifies Python environment and package management. If you don’t have it yet:

```p
# macOS / Linux
curl -LsSf https://astral.sh/uv/install.sh | sh  
```

```powershell
# Windows
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
$env:Path += ";$env:USERPROFILE\.local\bin"
```

Create and activate a virtual environment with uv:

```bash
# macOS / Linux
uv venv
source .venv/bin/activate
```

```powershell
# Windows
uv venv
.venv\Scripts\activate
```

## Install dependencies 

```bash
uv pip install torch transformers numpy pandas RDKit sentencepiece
```

### 2. Clone the CataPro repository

```bash
git clone https://github.com/zchwang/CataPro
```

### 3. Set up Git LFS 

CataPro uses [Git Large File Storage (LFS)](https://git-lfs.github.com/) to handle large model files. 
If you don't have Git LFS installed, you can install it using the following command:

```bash
git lfs install
```

### 4. Download the models

In addition, CataPro also relies on additional pre-trained models, including [prot_t5_xl_uniref50](https://huggingface.co/Rostlab/prot_t5_xl_uniref50) and [molt5-base-smiles2caption](https://huggingface.co/laituan245/molt5-base-smiles2caption). These two models are used for extracting features from enzymes and substrates, respectively.

> [!WARNING] 
> The models prot_t5_xl_uniref50 and molt5-base-smiles2caption required for CataPro are 64 and 1.9 GB, 
> respectively. 

```bash
# macOS / Linux
cd CataPro/models/

LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/Rostlab/prot_t5_xl_uniref50
cd prot_t5_xl_uniref50
git lfs pull

cd ..

LFS_SKIP_SMUDGE=1 git clone https://huggingface.co/laituan245/molt5-base-smiles2caption
cd molt5-base-smiles2caption
git lfs pull

cd ../..
```

```powershell
# Windows 
git -c filter.lfs.smudge= -c filter.lfs.required=false clone https://huggingface.co/Rostlab/prot_t5_xl_uniref50
cd prot_t5_xl_uniref50
git lfs pull

cd ..

git -c filter.lfs.smudge= -c filter.lfs.required=false clone https://huggingface.co/laituan245/molt5-base-smiles2caption
cd molt5-base-smiles2caption
git lfs pull

cd ../..
```

---

## Contact
Zechen Wang, PhD, Shandong University, wangzch97@gmail.com</p>

---

## Usage
### 1. Prepare the input files for inference
Enzyme and substrate information should be organized in a DataFrame created with pandas (in CSV format). Each enzyme-substrate pair must include the Enzyme_id, type (wild-type or mutant), the enzyme sequence, and the substrate's SMILES. The format is as follows:

|   Enzyme_id   |   type    |   sequence    |   smiles  |
|---------------|-----------|---------------|-----------|
|   Q6WZB0  |   wild    |   MTESPTTHHGAAPPDSV...    |   C(CC(C(=O)O)N)CN=C(N)N  |
|   B2MWN0  |   wild    |   MSSCQWSSFTRVSLSPF...    |   C(C(C(=O)O)N)S  |
|   C5AQP6  |   wild    |   MVLTRFPRVALTDGPTP...    |   C(C(C(=O)O)N)S  |

You can also refer to a sample file samples/sample_inp.csv

### 2. Next, you can use the following command to run CataPro to infer the kinetic parameters of the enzymatic reaction:

```bash
# In CataPro folder 
python inference/predict.py \
        -inp_fpath samples/sample_inp.csv \
        -model_dpath models \
        -batch_size 64 \
        -device cuda:0 \
        -out_fpath catapro_prediction.csv
```

Finally, the prediction results from CataPro are stored in the "catapro_prediction.csv" file. You can also run "bash run_catapro.sh" directly in the inference directory to achieve the above process.

---

## Question and Answer
To be updated ...
