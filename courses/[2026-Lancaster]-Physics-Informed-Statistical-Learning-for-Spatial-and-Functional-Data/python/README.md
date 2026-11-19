<p align="center">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_dark.png">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_light.png">
  <img alt="Fallback image description" src="https://raw.githubusercontent.com/fdaPDE/.github/refs/heads/main/profile/fdaPDE_light.png">
</picture>
</p>

Material for the short course "Physics-Informed Statistical Learning for Spatial and Functional Data" held during the lecture series in "Advanced Topics in Statistical Machine Learning", Rome (IT), 18 February 2026. 

__Lectures: Laura Maria Sangalli <sup>1</sup>, Ilenia Di Battista <sup>1</sup>__.

&nbsp;&nbsp;&nbsp;&nbsp;_<sup>1</sup> MOX, Department of Mathematics, Politecnico di Milano_

# Getting started

`fdaPDE` is a C++ library interfacing with Python, one of the most widely used languages for data analysis. The version of `fdaPDE` used in this course is not yet available as an official Python package for all platforms. A testing version can be found on [TestPyPI](https://test.pypi.org/project/fdaPDE/2.0.0/). 

All notebooks can be executed on [Google Colab](https://colab.google/). More informations on the "Colab instructions" section. If you instead prefer to run the code locally, follow the "Installation instructions" section to set up your system.

## Colab instructions

To run the notebooks on [Google Colab](https://colab.google/), just open one of the link below. If necessary, grant permissions to Google to trust the notebook.

**Be sure to run the setup cell at the beginning of each notebook to download the data and install all the required dependencies**.

* [Spatial regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Rome]-Advanced-Topics-in-Statistical-Machine-Learning/notebook/SRPDE_2D.ipynb)
* [Generalized regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Rome]-Advanced-Topics-in-Statistical-Machine-Learning/notebook/GSRPDE_2D.ipynb)
* [Quantile regression](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Rome]-Advanced-Topics-in-Statistical-Machine-Learning/notebook/QSRPDE_2D.ipynb)
* [Density estimation](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Rome]-Advanced-Topics-in-Statistical-Machine-Learning/notebook/DEPDE_2D.ipynb)
* [functional Principal Component Analysis](https://colab.research.google.com/github/fdaPDE/course-materials/blob/main/courses/[2026-Rome]-Advanced-Topics-in-Statistical-Machine-Learning/notebook/fPCA_2D.ipynb)

## Installation instructions

`fdaPDE` requires Python 3.11 or higher. If you do not have Python installed, or if your version is too old, install or upgrade Python according to your operating system.

* **Linux**
  
  The commands reported here assume an Ubuntu (or Debian-based) distribution with the `apt` package manager. If you use a different distribution, install all system-level dependencies according to your distribution.
  
  First, check your Python version:
  
  ```
  python --version
  ```
  
  If needed, install or upgrade Python. In addition, make sure that `venv` is available on your system. You can install it with:
  
  ```
  sudo apt install -y python3-venv python3-dev
  ```

  From a folder of your choice, execute the following commands:
  
  ```
  # create a folder for the course
  mkdir rome26_fdapde_course
  cd rome26_fdapde_course
  
  # clone the course material
  git clone https://github.com/fdaPDE/course-materials
  cd course-materials/courses/\[2026-Rome\]-Advanced-Topics-in-Statistical-Machine-Learning/
  
  # create the Python virtual environment for this course
  python -m venv .venv
  source .venv/bin/activate
     
    # for non-Ubuntu users, the activation script may be suffixed by the shell extension
    source .venv/bin/activate.sh    # if your shell is bash or zsh
    source .venv/bin/activate.fish  # if your shell is fish
     
  # install the required Python dependencies
  pip install --upgrade pip
  pip install -r requirements.txt
  pip install notebook jupyterlab ipykernel
  
  # register the virtual environment 
  python -m ipykernel install --user --name=.venv --display-name "fdaPDE-py, Rome26"
  
  # install the fdaPDE package from TestPyPI
  pip install -i https://test.pypi.org/simple/ fdaPDE
  ```
  
  Some graphical dependencies require GDAL. If the installation of the requirements fails, make sure that GDAL is installed. If not, execute:
  
  ```
  sudo apt install -y libgdal-dev gdal-bin
  gdal-config --version     # this returns a version number in the format x.y.z
  pip install GDAL==x.y.z
  ```
  
  Finally, from the folder `courses/\[2026-Rome\]-Advanced-Topics-in-Statistical-Machine-Learning/`, launch jupyter by executing
  
  ```
  jupyter notebook
  ```
  
  If your browser does not open Jupyter automatically, inspect the output produced by the command above and copy-paste the link
`http://localhost:8888/tree?token=...` into your browser.
  
  Select a notebook from the `notebook/` folder, then select the "fdaPDE-py, Rome26" kernel by clicking on "Kernel" → "Change kernel..." and then selecting "fdaPDE-py, Rome26".


  You are now up and running. Enjoy the course!

* **MacOS**

  Check your python version 
  
  ```
  python3 --version
  ```
  If needed, install or upgrade Python. You can install Python either from the [official site](https://www.python.org/downloads/mac-osx/) or via [Homebrew](https://docs.brew.sh/Homebrew-and-Python).
  
  From a folder of your choice, execute the following commands
  
  ```
  # create a folder for the course
  mkdir rome26_fdapde_course
  cd rome26_fdapde_course
  
  # clone the course material
  git clone https://github.com/fdaPDE/course-materials
  cd course-materials/courses/\[2026-Rome\]-Advanced-Topics-in-Statistical-Machine-Learning/
  
  # create the Python virtual environment for this course
  python3 -m venv .venv
  source .venv/bin/activate
     
  # install the required python dependencies
  pip install --upgrade pip
  pip install -r requirements.txt
  pip install notebook jupyterlab ipykernel
  
  # register the virtual environment 
  python3 -m ipykernel install --user --name=.venv --display-name "fdaPDE-py, Rome26"
  
  # install the fdaPDE package from TestPyPI
  pip install -i https://test.pypi.org/simple/ fdaPDE
  ```
  
  Some graphical dependencies require GDAL. If the installation of the requirements fails, make sure that GDAL is installed. If not, execute:

  ```
  brew install gdal
  ```

  Finally, from the folder `courses/\[2026-Rome\]-Advanced-Topics-in-Statistical-Machine-Learning/`, launch jupyter by executing
  
  ```
  jupyter notebook
  ```
  
  If your browser does not open Jupyter automatically, inspect the output produced by the command above and copy-paste the link
`http://localhost:8888/tree?token=...` into your browser.
  
  Select a notebook from the `notebook/` folder, then select the "fdaPDE-py, Rome26" kernel by clicking on "Kernel" → "Change kernel..." and then selecting "fdaPDE-py, Rome26".


  You are now up and running. Enjoy the course!

* **Windows**

  Due to current incompatibilities between `fdaPDE` and MSVC, Windows users must rely on Windows Subsystem for Linux (WSL) to run the notebooks. We strongly recommend running the course material on [Google Colab](https://colab.google/) if the procedure below fails or if you prefer not to set up WSL on your machine.
  
  If you have never configured WSL before, open Command Prompt or PowerShell as Administrator and execute:
  
  ```
  wsl --install
  ```
  
  By default, the Ubuntu distribution will be installed. Follow the instructions to set up your account. Keep track of the chosen `username` and `password`.
  
  From WSL, update your system:
  
  ```
  sudo apt update
  sudo apt upgrade -y
  ```
  
  If prompted, enter the password of your Ubuntu account (the one chosen during setup).
  
  Install Python (the instruction below install `python3.12` from the deadsnakes PPA)
  
  ```
  sudo apt install -y build-essential
  sudo apt install -y software-properties-common
  sudo add-apt-repository ppa:deadsnakes/ppa
  sudo apt update
  sudo apt install -y python3.12 python3.12-venv python3.12-dev
  
  # check python version
  python3.12 --version
  
  # install pip
  sudo apt install -y python3-pip
  ```
  
  From a folder of your choice, execute the following commands
  
  ```
  # create a folder for the course
  mkdir rome26_fdapde
  cd rome26_fdapde
  
  # clone the course material
  git clone https://github.com/fdaPDE/course-materials
  cd course-materials/courses/\[2026-Rome\]-Advanced-Topics-in-Statistical-Machine-Learning/
  
  # create the Python virtual environment for this course
  python3.12 -m venv .venv
  source .venv/bin/activate
     
  # install the required python dependencies
  pip install --upgrade pip
  pip install -r requirements.txt
  pip install notebook jupyterlab ipykernel
  
  # register the virtual environment 
  python3.12 -m ipykernel install --user --name=.venv --display-name "fdaPDE-py, Rome26"
  
  # install the fdaPDE package from TestPyPI
  pip install -i https://test.pypi.org/simple/ fdaPDE
  ```
  
  Some graphical dependencies require GDAL. If the installation of the requirements fails, make sure that GDAL is installed. If not, execute:
  
  ```
  sudo apt install -y libgdal-dev gdal-bin
  gdal-config --version     # this returns a version number in the format x.y.z
  pip install GDAL==x.y.z
  ```
  
  Finally, from the folder `courses/\[2026-Rome\]-Advanced-Topics-in-Statistical-Machine-Learning/`, launch jupyter by executing
  
  ```
  jupyter notebook
  ```
  
  Inspect the output produced by the command above and copy-paste the link `http://localhost:8888/tree?token=...` into your browser.

   Select a notebook from the `notebook/` folder, then select the "fdaPDE-py, Rome26" kernel by clicking on "Kernel" → "Change kernel..." and then selecting "fdaPDE-py, Rome26".


  You are now up and running. Enjoy the course!
  
