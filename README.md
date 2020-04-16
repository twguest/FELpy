## SPB MODEL MAXWELL SETUP
The SPB model code, and accompanying WavefrontPropaGator (WPG) branch is located on github - contact Trey for access (trey.guest@desy.de). A workflow for installation to MAXWELL@DESY is outlined below.

Recommended storage: 

> /gfps/exfel/data/user/USERNAME/

In terminal, move to storage directory and activate conda

    cd /gfps/exfel/data/user/USERNAME
    load module anaconda/
   Clone the model repo and WPG branch via git
   

    git clone https://github.com/twguest/spb_model.git
    git clone https://github.com/twguest/WPG.git
   Create and activate the conda environment required to interpret the model:
   

    conda env create -f WPG/conda-env.yml
    conda activate optics
Copy the build folder containing the chosen SRW distribution (via Radiasoft-Sirepo) to WPG/ and unzip it.

    cp /gfps/exfel/data/user/guestt/build.zip /gfps/exfel/data/user/USERNAME/WPG/
    unzip WPG/build.zip
   Delete previous build data (to be ammended in future):
   

    rm WPG/build/tmp
    rm WPG/build/libs
  Compile the SRW library via WPG:
  

    cd WPG/
    make all



