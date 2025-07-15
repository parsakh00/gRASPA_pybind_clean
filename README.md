# gRASPA pybind: pybind extension for gRASPA (test)
* enjoy (?) gRASPA and GPU backend while on python and Jupyter Notebook!

## Advantages
* **Access**:  Get access to gRASPA's internal variables
* **Control**: See what happens during an MC move
* **Create**:  Create your own MC moves with python code and packaged MC move parts

## How to use
0. Install pybind11 via `pip install pybind11`
    * check [pybind installation](https://pybind11.readthedocs.io/en/stable/installing.html) for more information!
1. download gRASPA_pybind and [Zhaoli2042/gRASPA_fork](https://github.com/Zhaoli2042/gRASPA_fork) or [snurr-group/gRASPA](https://github.com/snurr-group/gRASPA)
2. copy the files in gRASPA's [src_clean](https://github.com/snurr-group/gRASPA/tree/main/src_clean) to gRASPA_pybind's [pybind_src](https://github.com/Zhaoli2042/gRASPA_pybind/tree/main/pybind_src)
3. Compile using `./BIND_NVC_COMPILE`
    * If successful, you can see a **shared library file (gRASPA.so)** being generated
4. Copy `gRASPA.so` file to the `JUPYTER_NOTEBOOKS` folder
5. Run `jupyter notebook` and open one of the example notebooks!
    * :memo: NOTE: you need to install jupyter notebook to run these examples!
    * check [jupyter notebook installation](https://jupyter.org/install#jupyter-notebook)

## :memo: NOTE
* the pybind extension is under rapid-development. New features and bug-fixes will be available soon!
* Please submit questions/bug of the extension via Issues
* Or contact me (Zhao Li) via zhaoli2023@u.northwestern.edu or zhaoli2042@gmail.com
